## VISOR + freqk pipeline (01–05)

This README describes how to run the parameterised benchmark pipeline (deletions **or** insertions) using VISOR + freqk.

---

### Quick start — using the launcher

The recommended way to run experiments is via `launch_experiment.py`, which generates a config and submits the Slurm pipeline automatically.

```bash
cd /home/tbellagio/scratch/visor_freqk

# Single position (default: 10 Mb)
python scripts/launch_experiment.py --sv-type DEL --coverage 50 --freq 0.50

# With sequencing errors and a non-default k
python scripts/launch_experiment.py --sv-type INS --coverage 100 --freq 0.25 \
    --error-rate 0.01 --k 21

# Replicates: three independent positions
python scripts/launch_experiment.py --sv-type DEL --coverage 50 --freq 0.50 \
    --positions 10000000 20000000 30000000

# Preview without submitting
python scripts/launch_experiment.py --sv-type DEL --coverage 50 --freq 0.50 \
    --positions 10000000 20000000 --dry-run
```

The launcher:
1. Derives a `POS_LABEL` for each position (e.g. `10000000` → `pos10mb`).
2. Generates one config file per position under `scripts/generated_configs/`.
3. Calls `bash scripts/run_pipeline.sh <config>` for each position, which submits the five Slurm jobs as a chained dependency chain (`01` → `02` → `03` → `04` → `05`).

**Launcher arguments:**

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--sv-type` | yes | — | `DEL` or `INS` |
| `--coverage` | yes | — | Sequencing depth (e.g. `50`) |
| `--freq` | yes | — | SV allele frequency (e.g. `0.50`) |
| `--positions` | no | `10000000` | One or more 0-based positions in bp (see constraint below) |
| `--error-rate` | no | `0.0` | Sequencing error rate (e.g. `0.01`) |
| `--k` | no | `31` | k-mer size for freqk |
| `--dry-run` | no | off | Print config and command, do not submit |

---

### Positions and replicates

Each genomic position is an independent **replicate** — the same SV is placed at a different location in the genome, giving a distribution of freqk estimates across different local sequence contexts rather than a single point.

#### Position → folder label mapping

Positions are converted to a human-readable `POS_LABEL` used in all output paths:

```
 10000000  →  pos10mb   →  data/reads/del/pos10mb/...
 20000000  →  pos20mb   →  data/reads/del/pos20mb/...
 30000000  →  pos30mb   →  data/reads/del/pos30mb/...
```

**Positions must be multiples of 1 Mb** (e.g. `5000000`, `10000000`, `25000000`).
Non-Mb positions (e.g. `12500000`) will produce a truncated label (`pos12mb`) that does not exactly represent the coordinate.

#### Safety rules for choosing positions

- Keep positions at least `10 000 bp` from chromosome ends (to fit the largest SV size, 10 kb).
- Keep positions at least `10 000 bp` apart from each other (so SVs do not overlap).
- Avoid known centromeres and highly repetitive regions if you want unambiguous k-mers.

#### Generated config name

For each position a config is written as:

```
scripts/generated_configs/config_<sv>_cov<COV>_f<FREQ%>_err<ERR>_k<K>_pos<N>mb.sh
```

Example for `--positions 10000000 20000000 30000000`:

```
config_del_cov50_f50_err0_k31_pos10mb.sh
config_del_cov50_f50_err0_k31_pos20mb.sh
config_del_cov50_f50_err0_k31_pos30mb.sh
```

#### WT clone is shared

The WT clone (reference with no variants) does not depend on position and is shared across all replicates for the same SV type:

```
data/haplotypes/del/_clone_WT/   ← shared across all DEL positions
data/haplotypes/ins/_clone_WT/   ← shared across all INS positions
```

---

### Manual config editing

If you prefer not to use the launcher, edit one of the static configs directly:

- Deletions: `scripts/config_sv_deletions.sh`
- Insertions: `scripts/config_sv_insertions.sh`

Key variables:

| Variable | Description |
|----------|-------------|
| `SV_START_0` | 0-based SV start position on `CHROM` |
| `POS_LABEL` | Human-readable label — **must match `SV_START_0`** (e.g. `pos10mb` for `10000000`) |
| `ANCHOR_POS` | DEL only: `SV_START_0 - 1` |
| `FREQ` | SV allele frequency |
| `COVERAGE` | Sequencing depth |
| `ERROR_RATE` | Sequencing error rate |
| `K` | k-mer size for freqk |

Then submit with:

```bash
bash scripts/run_pipeline.sh scripts/config_sv_deletions.sh
```

---

### Pipeline steps

#### Step 1 — Make BEDs (01)

Writes one BED file per SV size into `${BEDS}` (= `data/beds/{sv}/{pos_label}/`):

- `SV_TYPE="DEL"` → `hack_del_<SIZE>.bed`
- `SV_TYPE="INS"` → `hack_ins_<SIZE>.bed`

#### Step 2 — Run VISOR HACk (02)

- Creates the shared WT clone at `data/haplotypes/{sv}/_clone_WT/` (skipped if it already exists).
- Writes variant haplotypes under `${HAPS}` (= `data/haplotypes/{sv}/{pos_label}/`):
  - DEL: `del_<SIZE>/HAP1/` and `HAP2/`
  - INS: `ins_<SIZE>/HAP1/` and `HAP2/`

#### Step 3 — Simulate pool-seq reads (03)

Simulates paired-end reads using VISOR SHORtS. Output under `${READS}` (= `data/reads/{sv}/{pos_label}/`):

```
cov<COV>/freq_<SIZE>_f<FREQ%>_err<ERR>/
  r1.fq
  r2.fq
```

Example: `data/reads/del/pos10mb/cov50/freq_1kb_f50_err0/`

#### Step 4 — Build VCFs (04)

Creates sequence-resolved VCFs under `${VCF_DIR}` (= `data/vcf/{sv}/{pos_label}/`):

- DEL: `del_<SIZE>.vcf.gz` + `.tbi`
- INS: `ins_<SIZE>.vcf.gz` + `.tbi`

#### Step 5 — Run freqk (05)

Runs `index` → `var-dedup` → `ref-dedup` → `count` → `call` for each SV size.
Results go under `${RESULTS}` (= `results/{sv}/{pos_label}/`):

```
cov<COV>_err<ERR>/<SV_TYPE>/<SIZE>/f<FREQ%>/k<K>/
```

Example: `results/del/pos10mb/cov50_err0/DEL/1kb/f50/k31/`

---

### Pipeline output paths (no overlap)

Every variable that affects the data is encoded in the path so different runs never overwrite each other.

| Step | Output location | Parameters in path | Reused? |
|------|----------------|--------------------|---------|
| 01 | `${BEDS}/hack_{del\|ins}_<SIZE>.bed` | SV_TYPE, **POS_LABEL**, size | yes — skip if exists |
| 02 | `${HAPS}/{del\|ins}_<SIZE>/HAP1,HAP2` | SV_TYPE, **POS_LABEL**, size | yes — skip if exists |
| 02 | `data/haplotypes/{sv}/_clone_WT` | SV_TYPE only (shared across positions) | yes — skip if exists |
| 03 | `${READS}/cov<COV>/freq_<SIZE>_f<FREQ%>_err<ERR>/` | SV_TYPE, **POS_LABEL**, COVERAGE, size, FREQ, ERROR_RATE | no — own dir per run |
| 04 | `${VCF_DIR}/{del\|ins}_<SIZE>.vcf.gz` | SV_TYPE, **POS_LABEL**, size | no — overwrites |
| 05 | `${RESULTS}/cov<COV>_err<ERR>/<SV_TYPE>/<SIZE>/f<FREQ%>/k<K>/` | SV_TYPE, **POS_LABEL**, COVERAGE, ERROR_RATE, size, FREQ, K | no — own dir per run |

**Design-only data** (BEDs, haplotypes, VCFs) depends on SV geometry (`SV_TYPE`, `SV_START_0`, sizes) but not on pool-seq parameters. Steps 01 and 02 skip if outputs already exist, so different coverage/freq/error runs at the same position reuse the same haplotypes.

**Experiment-specific data** (reads, results) always gets its own directory because COVERAGE, FREQ, ERROR_RATE, and K are all in the path.

---

### Full path example

Running:

```bash
python scripts/launch_experiment.py \
    --sv-type DEL --coverage 50 --freq 0.50 --error-rate 0.01 \
    --positions 10000000 20000000
```

produces:

```
data/
  beds/del/
    pos10mb/hack_del_100bp.bed  ...  hack_del_10kb.bed
    pos20mb/hack_del_100bp.bed  ...  hack_del_10kb.bed
  haplotypes/del/
    _clone_WT/                        ← shared across all positions
    pos10mb/del_1kb/HAP1/  HAP2/
    pos20mb/del_1kb/HAP1/  HAP2/
  reads/del/
    pos10mb/cov50/freq_1kb_f50_err01/r1.fq  r2.fq
    pos20mb/cov50/freq_1kb_f50_err01/r1.fq  r2.fq
  vcf/del/
    pos10mb/del_1kb.vcf.gz
    pos20mb/del_1kb.vcf.gz

results/del/
  pos10mb/cov50_err01/DEL/1kb/f50/k31/
  pos20mb/cov50_err01/DEL/1kb/f50/k31/
```

---

## Variation-aware pipeline (var)

The variation-aware pipeline replaces the two synthetic haplotypes (SV + WT) with **N real ecotype haplotypes** drawn from the GrENET population VCF. This models realistic SNP background variation and allows freqk to be benchmarked under natural sequence diversity.

### How it differs from the standard pipeline

| Aspect | Standard pipeline | Variation-aware pipeline |
|--------|------------------|--------------------------|
| Haplotypes | 2 (one SV + one WT reference clone) | N ecotypes from GrENET VCF |
| SNP background | None (clean reference) | Per-ecotype SNPs from `bcftools consensus` |
| SV carriers | 1 clone (SV haplotype) | N_SV = round(N × SV_FREQ) randomly assigned |
| SHORtS input | 2 clone dirs | N clone dirs at equal fractions (100/N each) |
| VCF used by freqk | Deletion only | Deletion only (SNPs are baked into FASTAs, not in the freqk VCF) |

### Prerequisites — one-time VCF setup

The GrENET VCF uses chromosome names `"1"`, while the reference FASTA uses `"Chr1"`. Run `00_prep_vcf.sh` once to produce a renamed VCF used by all downstream steps:

```bash
sbatch scripts/00_prep_vcf.sh scripts/config_sv_var_deletions.sh
```

Produces: `data/reference/greneNet_final_v1.1.recode.Chr1.vcf.gz` (+ `.tbi`).
This step is idempotent — it exits immediately if the renamed VCF already exists.

---

### Quick start — variation-aware launcher

```bash
cd /home/tbellagio/scratch/visor_freqk

# Single position, 10 ecotypes, 50 % SV frequency, 1 kb deletion
python scripts/launch_experiment_var.py \
    --sv-type DEL --coverage 50 --sv-freq 0.50 --n-samples 10 \
    --sizes 1kb --positions 10000000

# Multiple SV sizes
python scripts/launch_experiment_var.py \
    --sv-type DEL --coverage 50 --sv-freq 0.50 --n-samples 10 \
    --sizes 1kb 5kb 10kb --positions 10000000

# Multiple positions (replicates at different genomic contexts)
python scripts/launch_experiment_var.py \
    --sv-type DEL --coverage 50 --sv-freq 0.50 --n-samples 10 \
    --sizes 1kb --positions 10000000 20000000 30000000

# Preview without submitting
python scripts/launch_experiment_var.py \
    --sv-type DEL --coverage 50 --sv-freq 0.50 --n-samples 10 \
    --sizes 1kb --positions 10000000 --dry-run
```

The launcher:
1. Generates one config file per position under `scripts/generated_configs/config_var_*.sh`.
2. Calls `bash scripts/run_pipeline_var.sh <config>` for each position.
3. Saves all Slurm job IDs to `logs/pipeline_var_<TIMESTAMP>.jobids` for timing queries.

**Launcher arguments:**

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--sv-type` | yes | — | `DEL` (INS not yet supported in var pipeline) |
| `--coverage` | yes | — | Sequencing depth (e.g. `50`) |
| `--sv-freq` | yes | — | Fraction of N_SAMPLES haplotypes carrying the SV (e.g. `0.50`) |
| `--n-samples` | yes | — | Number of ecotypes drawn from the GrENET VCF (e.g. `10`) |
| `--sizes` | no | `1kb` | One or more SV sizes: `100bp` `500bp` `1kb` `5kb` `10kb` |
| `--positions` | no | `10000000` | One or more 0-based positions in bp |
| `--error-rate` | no | `0.001` | Sequencing error rate |
| `--k` | no | `31` | k-mer size for freqk |
| `--dry-run` | no | off | Print config and command, do not submit |

---

### Variation-aware pipeline steps

#### Step 00a — Rename VCF chromosomes (`00_prep_vcf.sh`)

One-time step. Renames chromosomes in the GrENET VCF from `"1"` to `"Chr1"` to match the reference FASTA, and subsets to Chr1 only. Output: `RENAMED_VCF`. Idempotent.

#### Step 00b — Apply VCF to produce per-ecotype FASTAs (`00_apply_vcf.sh`)

Runs `bcftools consensus` on each of the N ecotype samples from `RENAMED_VCF`, producing one FASTA per ecotype:

```
data/haplotypes_var/del/<pos_label>/s_<sample>/h1.fa
```

Each FASTA is Chr1 with that ecotype's SNP background applied (`--haplotype 1`, appropriate for inbred Arabidopsis). Deletion coordinates remain valid because the GrENET VCF contains SNPs only (no indels → no coordinate shift).

#### Step 01 — Make BEDs (`01_make_beds.sh`)

Shared with the standard pipeline. Writes SV BED files to `data/beds/del/<pos_label>/`. Skipped if outputs already exist.

#### Step 02 — VISOR HACk on SV carriers (`02_run_hack_var.sh`)

Injects the deletion into N_SV = round(N × SV_FREQ) of the ecotype FASTAs:

```
data/haplotypes_var/del/<pos_label>/s_<sample>_sv_del_<SIZE>/h1.fa  ← SV clone
data/haplotypes_var/del/<pos_label>/s_<sample>/h1.fa                ← WT clone (untouched)
```

SV carriers are always the first N_SV samples in VCF header order.

#### Step 03 — Simulate pool-seq reads (`03_run_shorts_var.sh`)

Pools all N clones via VISOR SHORtS at equal fractions (100/N each). SV clones come first in the clone list, then WT clones. Output:

```
data/reads_var/del/<pos_label>/cov<COV>/var_del_<SIZE>_n<N>_f<FREQ%>_err<ERR>/
  r1.fq
  r2.fq
```

Example: `data/reads_var/del/pos10mb/cov50/var_del_1kb_n10_f50_err001/`

#### Step 04 — Build SV VCFs (`04_make_vcf_var.sh`)

Creates sequence-resolved deletion VCFs in the `/var` subdirectory (separate from the standard pipeline, for future-proofing):

```
data/vcf/del/<pos_label>/var/del_<SIZE>.vcf.gz  (+.tbi)
```

Runs independently of steps 00–03 (no dependency on haplotypes or reads).

#### Step 05 — Run freqk (`05_freqk_var.sh`)

Same freqk steps as the standard pipeline (`index` → `var-dedup` → `ref-dedup` → `count` → `call`). Results go to the `/var` subdirectory:

```
results/del/<pos_label>/var/cov<COV>_err<ERR>/DEL/<SIZE>/n<N>/f<FREQ%>/k<K>/
  var_del_<SIZE>_n<N>_f<FREQ%>_err<ERR>.allele_frequencies.k<K>.tsv  ← estimated AF
  var_del_<SIZE>_n<N>_f<FREQ%>_err<ERR>.counts_by_allele.k<K>.tsv
  var_del_<SIZE>_n<N>_f<FREQ%>_err<ERR>.raw_kmer_counts.k<K>.tsv
  del_<SIZE>.k<K>.freqk.{index,var_index,ref_index}
```

The `n<N>` subdirectory ensures runs with different numbers of ecotypes (e.g. `--n-samples 10` vs `--n-samples 100`) are fully isolated and never overwrite each other.

---

### Variation-aware job chain

```
00_prep_vcf ──┬── 00_apply_vcf ──┐
              └── 01_make_beds  ──┴── 02_run_hack_var ── 03_run_shorts_var ──┐
04_make_vcf_var (no deps, starts immediately) ────────────────────────────────┴── 05_freqk_var
```

`04_make_vcf_var` only needs the reference FASTA so it runs in parallel with the haplotype chain, saving wall time.

---

### Checking pipeline timing

After submitting via the launcher, job IDs are saved to `logs/pipeline_var_<TIMESTAMP>.jobids`. Once all jobs finish:

```bash
bash scripts/check_timing.sh
# auto-selects the latest .jobids file, or pass explicitly:
bash scripts/check_timing.sh logs/pipeline_var_20260316_143200.jobids
```

Reports elapsed time per step and total wall time (first job start → last job end).

---

### Full path example — variation-aware

Running:

```bash
python scripts/launch_experiment_var.py \
    --sv-type DEL --coverage 50 --sv-freq 0.50 --n-samples 10 \
    --sizes 1kb --positions 10000000
```

produces:

```
data/
  reference/
    greneNet_final_v1.1.recode.Chr1.vcf.gz   ← renamed VCF (one-time)
  haplotypes_var/del/pos10mb/
    s_<ecotype_0>/h1.fa                       ← WT consensus (SNP background)
    s_<ecotype_0>_sv_del_1kb/h1.fa            ← SV clone (deletion injected)
    ...                                        ← (10 WT + 5 SV dirs)
  reads_var/del/pos10mb/
    cov50/var_del_1kb_n10_f50_err001/r1.fq r2.fq
  vcf/del/pos10mb/var/
    del_1kb.vcf.gz  (+.tbi)

results/del/pos10mb/var/
  cov50_err001/DEL/1kb/n10/f50/k31/
    var_del_1kb_n10_f50_err001.allele_frequencies.k31.tsv
```

### Output path isolation

The var pipeline uses entirely separate directory roots from the standard pipeline:

| Data type | Standard pipeline | Variation-aware pipeline |
|-----------|------------------|--------------------------|
| Haplotypes | `data/haplotypes/` | `data/haplotypes_var/` |
| Reads | `data/reads/` | `data/reads_var/` |
| VCFs | `data/vcf/del/<pos>/` | `data/vcf/del/<pos>/var/` |
| Results | `results/del/<pos>/` | `results/del/<pos>/var/` |
| BEDs | `data/beds/del/<pos>/` | `data/beds/del/<pos>/` ← shared |
