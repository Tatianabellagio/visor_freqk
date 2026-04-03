# visor_freqk

Benchmarking pipeline for [`freqk`](https://github.com/grenenet/freqk), a k-mer-based tool for estimating structural variant (SV) allele frequencies from pool-seq short-read data. The pipeline uses [VISOR](https://davidlaewen.github.io/VISOR/) to simulate haplotypes with known deletions or insertions at controlled frequencies, then evaluates how accurately freqk recovers those frequencies across a range of SV sizes, sequencing depths, error rates, and k-mer sizes.

Two modes are supported:
- **Standard pipeline** — clean synthetic haplotypes (SV carrier + WT reference clone)
- **Variation-aware pipeline** — realistic GrENE-Net ecotype backgrounds drawn from the population VCF, modelling natural SNP diversity

Compute jobs are submitted as chained SLURM dependency chains on an HPC cluster.

---

## Quick start

### Standard pipeline

```bash
# Single position (default: 10 Mb), 50x coverage, 50% SV frequency
python scripts/launch_experiment.py --sv-type DEL --coverage 50 --freq 0.50

# With sequencing errors and custom k
python scripts/launch_experiment.py --sv-type INS --coverage 100 --freq 0.25 \
    --error-rate 0.01 --k 21

# Three independent replicates at different genomic positions
python scripts/launch_experiment.py --sv-type DEL --coverage 50 --freq 0.50 \
    --positions 10000000 20000000 30000000

# Preview without submitting
python scripts/launch_experiment.py --sv-type DEL --coverage 50 --freq 0.50 \
    --positions 10000000 20000000 --dry-run
```

### Variation-aware pipeline

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

### Standard (`launch_experiment.py`)

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--sv-type` | yes | — | `DEL` or `INS` |
| `--coverage` | yes | — | Sequencing depth |
| `--freq` | yes | — | SV allele frequency (0–1) |
| `--positions` | no | `10000000` | One or more 0-based positions (multiples of 1 Mb) |
| `--error-rate` | no | `0.0` | Sequencing error rate |
| `--k` | no | `31` | k-mer size for freqk |
| `--dry-run` | no | off | Print config; do not submit |

### Variation-aware (`launch_experiment_var.py`)

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
| 01 | `01_make_beds.sh` | Write SV BED files for each size |
| 02 | `02_run_hack.sh` | VISOR HACk — inject SVs into haplotype FASTAs |
| 03 | `03_run_shorts.sh` | VISOR SHORtS — simulate pool-seq paired-end reads |
| 04 | `04_make_vcf.sh` | Build sequence-resolved VCFs |
| 05 | `05_freqk.sh` | Run freqk: index → var-dedup → ref-dedup → count → call |

In the variation-aware pipeline, steps 00a (`prep_vcf`) and 00b (`apply_vcf`) precede these to build per-ecotype FASTAs from the GrENE-Net population VCF.

---

## Output structure

All parameters are encoded in output paths so runs never overwrite each other:

```
results/{sv}/{pos_label}/cov{COV}_err{ERR}/{SV_TYPE}/{SIZE}/f{FREQ}/k{K}/
results/{sv}/{pos_label}/var/cov{COV}_err{ERR}/{SV_TYPE}/{SIZE}/n{N}/f{FREQ}/k{K}/  # var pipeline
```

---

## Checking pipeline timing

Job IDs are saved to `logs/pipeline_var_<TIMESTAMP>.jobids`. After jobs complete:

```bash
bash scripts/check_timing.sh
# or pass explicitly:
bash scripts/check_timing.sh logs/pipeline_var_20260316_143200.jobids
```

---

## Tools & dependencies

- [VISOR](https://davidlaewen.github.io/VISOR/) — SV haplotype simulation and pool-seq read generation
- [freqk](https://github.com/grenenet/freqk) — k-mer-based SV allele frequency estimation
- [bcftools](https://samtools.github.io/bcftools/) — VCF manipulation and consensus sequence generation
- SLURM — job scheduling
- Python 3, Bash
