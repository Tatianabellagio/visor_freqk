## VISOR + freqk pipeline (01–05)

This README describes how to run the parameterised benchmark pipeline (deletions **or** insertions) using VISOR + freqk.

---

### Quick start — using the launcher

The recommended way to run experiments is via `launch_experiment.py`, which generates a config and submits the Slurm pipeline automatically.

```bash
cd /home/tbellagio/scratch/pang/visor_freqk

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
