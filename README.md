# visor_freqk

Benchmarking pipeline for [`freqk`](https://github.com/milesroberts-123/freqk/tree/main/src), a k-mer-based tool for estimating structural variant (SV) allele frequencies from pool-seq short-read data. Uses [VISOR](https://github.com/davidebolo1993/VISOR) to simulate haplotypes with known deletions at controlled frequencies against realistic GrENE-Net ecotype backgrounds drawn from the population VCF, then evaluates how accurately freqk recovers those frequencies across a range of SV sizes, sequencing depths, error rates, and k-mer sizes.

A parallel step benchmarks freqk against **vg giraffe + vg call**, a pangenome-graph SV genotyper, using the same simulated reads.

Compute jobs are submitted as chained SLURM dependency chains on an HPC cluster.

---

## Quick start

```bash
# One-time VCF chromosome rename (run once per cluster)
sbatch scripts/00_prep_vcf.sh scripts/config_sv_var_deletions.sh

# Single position, 10 ecotypes
python scripts/launch_experiment_var.py \
    --sv-type DEL --coverage 50 --sv-freq 0.50 --n-samples 10 \
    --sizes 1kb --positions 10000000

# Multiple SV sizes and replicates
python scripts/launch_experiment_var.py \
    --sv-type DEL --coverage 50 --sv-freq 0.50 --n-samples 10 \
    --sizes 1kb 5kb 10kb --positions 10000000 20000000 30000000

# Preview
python scripts/launch_experiment_var.py \
    --sv-type DEL --coverage 50 --sv-freq 0.50 --n-samples 10 \
    --sizes 1kb --positions 10000000 --dry-run
```

---

## Launcher arguments

### `launch_experiment_var.py`

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--sv-type` | yes | — | `DEL` (INS not yet supported) |
| `--coverage` | yes | — | Sequencing depth |
| `--sv-freq` | yes | — | Fraction of ecotypes carrying the SV |
| `--n-samples` | yes | — | Number of ecotypes drawn from GrENE-Net VCF |
| `--sizes` | no | `1kb` | SV sizes: `100bp` `500bp` `1kb` `5kb` `10kb` |
| `--positions` | no | `10000000` | One or more 0-based positions |
| `--error-rate` | no | `0.001` | Sequencing error rate |
| `--k` | no | `31` | k-mer size |
| `--dry-run` | no | off | Print config; do not submit |

---

## Pipeline steps

| Step | Script | Description |
|------|--------|-------------|
| 00a | `00_prep_vcf.sh` | Rename VCF chromosomes ("1" → "Chr1") |
| 00b | `00_apply_vcf.sh` | Build per-ecotype consensus FASTAs from the GrENE-Net population VCF |
| 01 | `01_make_beds.sh` | Write SV BED files for each size |
| 02 | `02_run_hack_var.sh` | VISOR HACk — inject SVs into N_SV ecotype haplotypes |
| 03 | `03_run_shorts_var.sh` | VISOR SHORtS — simulate pool-seq paired-end reads |
| 04 | `04_make_vcf_var.sh` | Build sequence-resolved truth VCFs |
| 05 | `05_freqk_var.sh` | Run freqk: index → var-dedup → ref-dedup → count → call |
| 05b | `05b_vg_giraffe.sh` | Benchmark: vg giraffe + vg call on the same reads |
| 06 | `06_compute_repeat_score.sh` | Annotate each rep position with its k-mer repeat score |

---

## Output structure

All parameters are encoded in output paths so runs never overwrite each other:

```
results/del/{pos_label}/var/cov{COV}_err{ERR}/{SIZE}/n{N}/f{FREQ}/k{K}/
```

Each terminal directory contains both freqk and vg giraffe allele-frequency TSVs:

```
var_del_{SIZE}_n{N}_f{FREQ}_err{ERR}.allele_frequencies.k{K}.tsv              # freqk
var_del_{SIZE}_n{N}_f{FREQ}_err{ERR}.vg_giraffe.allele_frequencies.k{K}.tsv   # vg giraffe
```

Both use the same schema: a single line `af_ref|af_alt`.

---

## vg giraffe benchmark (step 05b)

`05b_vg_giraffe.sh` runs the standard graph-SV-genotyping workflow on each simulated pool:

1. `vg autoindex --workflow giraffe` builds a pangenome graph from the reference FASTA + truth VCF (cached per `pos_label × SIZE` under `data/vg_indexes/`).
2. `vg giraffe` maps the simulated reads to the graph.
3. `vg pack` aggregates read support onto graph edges.
4. `vg call -v <truth.vcf> -a` genotypes the known deletion site, emitting per-allele depth (AD).
5. `af_alt = AD_alt / (AD_ref + AD_alt)`; written to the same directory as freqk's output.

`pang` is the only conda/mamba env required — it provides `vg`, `samtools`, and `bcftools`.

To run the benchmark on already-simulated reads (no VISOR re-run):

```bash
python scripts/launch_benchmarks.py --rep-ids rep1 rep2 \
    --sizes 1kb --coverage 50 --sv-freq 0.50 --n-samples 231
```

### Results vs freqk

`summaries/freqk_vs_vg_giraffe.ipynb` scans both methods' AF TSVs side-by-side, filters to (rep, coverage) cells where both methods have the full 5 sizes × 5 freqs grid, and partitions the comparison by coverage, SV size, and per-combo paired error.

Current matched set (`rep1`–`rep10`, `k=31`):

| Coverage | Reps with full grid | Combos |
|----------|---------------------|--------|
| 10× | rep1 – rep10 (10) | 250 |
| 20× | rep1 – rep10 (10) | 250 |
| 50× | rep1, rep2, rep3, rep4, rep10 (5) | 125 |

The cov=50× expansion is still in flight for `rep5`–`rep9`.

**Accuracy.** On combos where both methods detect the variant, vg giraffe's |estimate − truth| is closer to zero on **357/529 (~67 %)** of paired combos; freqk on 171.

**Detection.** freqk detects 100 % of matched combos. vg giraffe drops 8–20 % of combos at 10×–20× coverage (returns NA when read support on the graph is too thin), and 14–23 % at 50×.

**Runtime (median seconds per combo, all 5 SV sizes)**

| Method | 10× | 20× | 50× |
|--------|----:|----:|----:|
| freqk (1 kb only — older combos lack timing TSVs) | 96 | 165 | 371 |
| vg giraffe | 180–202 | 290–315 | 520–820 |

For vg, the breakdown is dominated by `vg pack` and `vg giraffe` mapping (each ~30–45 % of wall time once the graph index is cached), with `vg call` ~10–20 % and `autoindex` once-per-rep.

Plots emitted by the notebook (under `plots/`):

- `freqk_vs_vg_pooled.png` — pooled true vs estimated AF
- `freqk_vs_vg_by_coverage.png` — 2 methods × 3 coverages
- `freqk_vs_vg_by_size.png` — 2 methods × 5 SV sizes
- `freqk_vs_vg_r2_heatmap.png` — per-method R² + (vg − freqk) diff
- `freqk_vs_vg_detection_rate.png` — detection bars by cov × size
- `freqk_vs_vg_paired_error.png` — |vg err| vs |freqk err| same combo
- `freqk_vs_vg_runtime_total.png` — total wall time per combo (boxplots)
- `freqk_vs_vg_runtime_vg_steps.png` — vg per-step breakdown

---

## Checking pipeline timing

Job IDs are saved to `logs/pipeline_var_<TIMESTAMP>.jobids`. After jobs complete:

```bash
bash scripts/check_timing.sh
# or pass explicitly:
bash scripts/check_timing.sh logs/pipeline_var_20260316_143200.jobids
```

---

## Repeat score

Each genomic position in the analysis is annotated with a **repeat score** — the fraction of 31-mers in a 10 kb window starting at the SV position that appear **≥2 times anywhere in Chr1**:

```
repeat_score(pos) = |{ k-mers in ref[pos : pos+10000] that appear ≥2× in Chr1 }|
                    ─────────────────────────────────────────────────────────────
                              total k-mers in ref[pos : pos+10000]
```

The genome-wide multi-copy k-mer set is built once from the full Chr1 FASTA (`data/reference/Chr1.fa`), so the score reflects global uniqueness: a k-mer is only counted as unique if it occurs exactly once across the entire chromosome.

Positions are then classified as:

| Score | Class |
|-------|-------|
| < 0.10 | low repeat |
| 0.10 – 0.50 | moderate repeat |
| ≥ 0.50 | high repeat |

**Where it is computed:**
- **`scripts/sample_positions.py`** — at position-sampling time; stored in `data/positions_registry.tsv` under `global_repeat_score`
- **`scripts/compute_repeat_score.py`** / **`scripts/06_compute_repeat_score.sh`** — as a pipeline job for each `rep*` position, updating the registry in-place
- **`summaries/results_var.ipynb`** / **`summaries/results_rep.ipynb`** — re-computed from the FASTA at analysis time using the same 10 kb window, so notebook scores match the registry

---

## Tools & dependencies

- [VISOR](https://github.com/davidebolo1993/VISOR) — SV haplotype simulation and pool-seq read generation
- [freqk](https://github.com/milesroberts-123/freqk/tree/main/src) — k-mer-based SV allele frequency estimation
- [vg](https://github.com/vgteam/vg) — pangenome graph construction, mapping (giraffe), and genotyping (call)
- [bcftools](https://samtools.github.io/bcftools/) — VCF manipulation and consensus sequence generation
- SLURM — job scheduling
- Python 3, Bash
