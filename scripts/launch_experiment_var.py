#!/usr/bin/env python3
"""
launch_experiment_var.py
========================
Generate variation-aware config files and submit the var pipeline to Slurm.
One config + one pipeline submission is created per (coverage × sv_freq × position) combination.

Usage
-----
  # Sample 5 random positions and run immediately
  python scripts/launch_experiment_var.py --sv-type DEL \\
      --coverage 10 20 50 --sv-freq 0.10 0.30 0.50 0.70 0.90 \\
      --n-samples 231 --sizes 100bp 500bp 1kb 5kb 10kb \\
      --random-n 5

  # Re-run previously registered replicates
  python scripts/launch_experiment_var.py --sv-type DEL \\
      --coverage 10 20 50 --sv-freq 0.10 0.30 0.50 \\
      --n-samples 231 --sizes 100bp 500bp 1kb 5kb 10kb \\
      --rep-ids rep1 rep2 rep3

  # Legacy: explicit bp positions (old results unaffected)
  python scripts/launch_experiment_var.py --sv-type DEL \\
      --coverage 10 20 --sv-freq 0.10 0.50 \\
      --n-samples 231 --positions 10000000 20000000

  # Dry run — print configs without submitting
  python scripts/launch_experiment_var.py --sv-type DEL --coverage 20 --sv-freq 0.50 \\
      --n-samples 231 --random-n 1 --dry-run

Position source (mutually exclusive — pick one)
------------------------------------------------
  --random-n N   Sample N new random positions, register them, submit immediately.
                 Uses hard constraints: ≥100kb from chromosome ends, no N-bases,
                 ≥100kb between positions.  Seed is auto-generated and printed.
                 Optional: --seed INT to make sampling reproducible.

  --rep-ids      Re-use positions already in data/positions_registry.tsv.

  --positions    Legacy mode: explicit bp coordinates → posXXX folder labels.
                 Old results are completely unaffected.

Required
--------
  --sv-type     DEL
  --coverage    one or more sequencing depths (e.g. --coverage 10 20 50)
  --sv-freq     one or more SV frequencies (e.g. --sv-freq 0.10 0.30 0.50)
  --n-samples   number of ecotype haplotypes from GrENET VCF

Optional
--------
  --sizes         SV sizes: 100bp 500bp 1kb 5kb 10kb  (default: 1kb)
  --seed          random seed for --random-n (default: random)
  --end-buffer    min distance from chromosome ends in bp (default: 100000)
  --min-sep       min distance between positions in bp (default: 100000)
  --error-rate    sequencing error rate (default: 0.001)
  --k             k-mer size for freqk (default: 31)
  --stagger-secs  seconds between pipeline submissions (default: 30)
  --dry-run       print configs and commands, do not submit
"""

import argparse
import csv
import datetime
import itertools
import random
import subprocess
import sys
import time
from pathlib import Path

# ── Fixed paths ───────────────────────────────────────────────────────────────
WORK     = "/home/tbellagio/scratch/visor_freqk"
FREQK    = "/home/tbellagio/scratch/bin/freqk"
FASTA    = Path(WORK) / "data/reference/Chr1.fa"
CHROM    = "Chr1"
REGISTRY = Path(WORK) / "data/positions_registry.tsv"

SV_START_0_DEFAULT = 10_000_000
MAX_DEL_SIZE       = 10_000   # largest deletion size — used for end-buffer constraint

# ── Supported SV sizes ────────────────────────────────────────────────────────
ALL_DEL_SIZES = {"100bp": 100, "500bp": 500, "1kb": 1000, "5kb": 5000, "10kb": 10000}
ALL_INS_SIZES = {"100bp": 100, "500bp": 500, "1kb": 1000, "5kb": 5000, "10kb": 10000}

REGISTRY_FIELDS = ["rep_id", "pos_bp", "chrom", "global_repeat_score", "seed", "date_added"]


# ── Registry helpers ──────────────────────────────────────────────────────────

def _load_registry() -> dict:
    """Return {rep_id: int(pos_bp)}."""
    if not REGISTRY.exists() or REGISTRY.stat().st_size == 0:
        return {}
    with open(REGISTRY) as f:
        return {row["rep_id"]: int(row["pos_bp"])
                for row in csv.DictReader(f, delimiter="\t")}


def _load_registry_rows() -> list:
    if not REGISTRY.exists() or REGISTRY.stat().st_size == 0:
        return []
    with open(REGISTRY) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _next_rep_id(rows: list) -> str:
    used = {int(r["rep_id"].replace("rep", "")) for r in rows}
    n = 1
    while n in used:
        n += 1
    return f"rep{n}"


def _append_registry(new_rows: list) -> None:
    write_header = not REGISTRY.exists() or REGISTRY.stat().st_size == 0
    REGISTRY.parent.mkdir(parents=True, exist_ok=True)
    with open(REGISTRY, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=REGISTRY_FIELDS, delimiter="\t")
        if write_header:
            writer.writeheader()
        for row in new_rows:
            writer.writerow(row)


# ── Random position sampling ──────────────────────────────────────────────────

def sample_positions(n: int, seed: int, end_buffer: int, min_sep: int) -> list:
    """
    Sample n random positions from Chr1 satisfying hard constraints.
    Runs entirely on the login node — no heavy computation:
      - FASTA is read once to check N-bases in candidate windows
      - repeat score is NOT computed here (stored as 'NA'; notebooks compute it)
    Returns list of (rep_id, pos_bp) tuples and writes to registry.
    """
    if not FASTA.exists():
        sys.exit(f"ERROR: reference FASTA not found: {FASTA}")

    # Read only the sequence — no k-mer counting, no repeat set
    print(f"Loading reference …", end=" ", flush=True)
    lines     = FASTA.read_text().splitlines()
    ref_seq   = "".join(l for l in lines if not l.startswith(">")).upper()
    chrom_len = len(ref_seq)
    print(f"done  ({chrom_len:,} bp)")

    existing_rows = _load_registry_rows()
    existing_pos  = {int(r["pos_bp"]) for r in existing_rows}
    print(f"Existing registry entries: {len(existing_rows)}")

    lo   = end_buffer
    hi   = chrom_len - MAX_DEL_SIZE - end_buffer
    STEP = 1_000
    rng  = random.Random(seed)
    candidates = list(range(lo, hi, STEP))
    rng.shuffle(candidates)

    accepted: list = []
    for pos in candidates:
        if len(accepted) == n:
            break
        # Hard constraint: no N-bases in the deletion window
        if "N" in ref_seq[pos : pos + MAX_DEL_SIZE]:
            continue
        # Hard constraint: minimum separation from existing + already accepted
        all_taken = existing_pos | set(accepted)
        if any(abs(pos - p) < min_sep for p in all_taken):
            continue
        accepted.append(pos)

    if len(accepted) < n:
        print(f"WARNING: only {len(accepted)}/{n} positions found — "
              f"try --min-sep {min_sep // 2}", file=sys.stderr)

    # Build registry rows (repeat score computed later by analysis notebooks)
    new_rows = []
    all_rows = existing_rows[:]
    for pos in accepted:
        rep_id = _next_rep_id(all_rows + new_rows)
        row    = {
            "rep_id":              rep_id,
            "pos_bp":              str(pos),
            "chrom":               CHROM,
            "global_repeat_score": "NA",   # computed in analysis notebooks
            "seed":                str(seed),
            "date_added":          datetime.date.today().isoformat(),
        }
        new_rows.append(row)
        print(f"  {rep_id:8s}  pos={pos:>12,}  ({pos/1e6:6.3f} Mb)")

    _append_registry(new_rows)
    print(f"Registered {len(new_rows)} new position(s) → {REGISTRY}\n")

    return [(r["rep_id"], int(r["pos_bp"])) for r in new_rows]


# ── Helpers ───────────────────────────────────────────────────────────────────

def format_error_label(error_rate: float) -> str:
    if error_rate == 0.0:
        return "0"
    s = f"{error_rate}"
    return s.split(".")[1] if "." in s else s


def pos_label_from_pos(pos: int) -> str:
    if pos % 1_000_000 == 0:
        return f"pos{pos // 1_000_000}mb"
    return f"pos{pos}"


# ── Config generation ─────────────────────────────────────────────────────────

def generate_config(sv_type: str, coverage: int, sv_freq: float, n_samples: int,
                    error_rate: float, k: int, pos: int, pos_label: str,
                    sizes: dict) -> str:
    sv         = sv_type.lower()
    size_var   = "DEL_SIZES" if sv_type == "DEL" else "INS_SIZES"
    size_lines = "\n".join(f'  ["{name}"]={val}' for name, val in sizes.items())
    anchor_pos = pos - 1

    return f"""\
#!/bin/bash
# Auto-generated by launch_experiment_var.py
# SV_TYPE={sv_type}  COVERAGE={coverage}  SV_FREQ={sv_freq}  N_SAMPLES={n_samples}
# ERROR_RATE={error_rate}  K={k}  POS={pos}  POS_LABEL={pos_label}

WORK={WORK}
REF=${{WORK}}/data/reference/Chr1.fa
VCF_FILE=${{WORK}}/data/reference/greneNet_final_v1.1.recode.vcf.gz
CHROM_MAP=${{WORK}}/data/reference/chrom_map.txt
RENAMED_VCF=${{WORK}}/data/reference/greneNet_final_v1.1.recode.Chr1.vcf.gz

SV_TYPE="{sv_type}"
CHROM="{CHROM}"
SV_START_0={pos}
ANCHOR_POS={anchor_pos}
POS_LABEL="{pos_label}"

N_SAMPLES={n_samples}
SV_FREQ={sv_freq}

declare -A {size_var}=(
{size_lines}
)

COVERAGE={coverage}
ERROR_RATE={error_rate}
FREQ=${{SV_FREQ}}

K={k}
FREQK={FREQK}

HAPS_WT=${{WORK}}/data/haplotypes_var/wt/n{n_samples}
HAPS_SV=${{WORK}}/data/haplotypes_var/{sv}/${{POS_LABEL}}
BEDS=${{WORK}}/data/beds/{sv}/${{POS_LABEL}}
READS_VAR=${{WORK}}/data/reads_var/{sv}/${{POS_LABEL}}
READS=${{READS_VAR}}
VCF_DIR=${{WORK}}/data/vcf/{sv}/${{POS_LABEL}}/var
RESULTS=${{WORK}}/results/{sv}/${{POS_LABEL}}/var
"""


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate variation-aware configs and submit the Slurm pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--sv-type",      required=True, choices=["DEL"])
    p.add_argument("--coverage",     required=True, type=int,   nargs="+")
    p.add_argument("--sv-freq",      required=True, type=float, nargs="+")
    p.add_argument("--n-samples",    required=True, type=int)
    p.add_argument("--sizes",        nargs="+", default=["1kb"],
                   choices=list(ALL_DEL_SIZES.keys()))
    p.add_argument("--error-rate",   default=0.001, type=float)
    p.add_argument("--k",            default=31,    type=int)
    p.add_argument("--stagger-secs", default=30,    type=int)
    p.add_argument("--dry-run",      action="store_true")
    p.add_argument("--seed",         default=None,  type=int,
                   help="Random seed for --random-n (default: auto-generated)")
    p.add_argument("--end-buffer",   default=100_000, type=int,
                   help="Min distance from chromosome ends in bp (default: 100000)")
    p.add_argument("--min-sep",      default=100_000, type=int,
                   help="Min distance between positions in bp (default: 100000)")

    pos_group = p.add_mutually_exclusive_group()
    pos_group.add_argument("--random-n",  type=int,
                           help="Sample N new random positions, register, and submit")
    pos_group.add_argument("--rep-ids",   nargs="+",
                           help="Re-use registered rep IDs (e.g. rep1 rep2)")
    pos_group.add_argument("--positions", nargs="+", type=int,
                           help="(legacy) explicit bp positions on Chr1")
    return p.parse_args()


def main():
    args = parse_args()

    # ── Resolve positions and labels ──────────────────────────────────────────
    if args.random_n:
        seed = args.seed if args.seed is not None else random.randint(0, 2**32 - 1)
        print(f"\nSampling {args.random_n} random position(s)  [seed={seed}]")
        sampled    = sample_positions(args.random_n, seed,
                                      args.end_buffer, args.min_sep)
        pos_entries = [(pos, rep_id) for rep_id, pos in sampled]
        mode = "random"

    elif args.rep_ids:
        registry = _load_registry()
        missing  = [r for r in args.rep_ids if r not in registry]
        if missing:
            sys.exit(f"ERROR: rep_id(s) not in registry: {', '.join(missing)}\n"
                     f"  Registry: {REGISTRY}")
        pos_entries = [(registry[r], r) for r in args.rep_ids]
        mode = "rep"

    else:
        positions   = args.positions if args.positions else [SV_START_0_DEFAULT]
        pos_entries = [(pos, pos_label_from_pos(pos)) for pos in positions]
        mode = "pos"

    coverages = args.coverage
    sv_freqs  = args.sv_freq
    all_sizes = ALL_DEL_SIZES if args.sv_type == "DEL" else ALL_INS_SIZES
    sizes     = {name: all_sizes[name] for name in args.sizes}

    scripts_dir = Path(__file__).parent
    config_dir  = scripts_dir / "generated_configs"
    config_dir.mkdir(exist_ok=True)
    err_label = format_error_label(args.error_rate)
    pipeline  = scripts_dir / "run_pipeline_var.sh"

    combinations = sorted(
        itertools.product(coverages, sv_freqs, pos_entries),
        key=lambda x: (x[2][1], x[0], x[1]),
    )
    total = len(combinations)

    print(f"\nLaunching {total} pipeline run(s)  [mode: {mode}]:")
    print(f"  Positions  : {[label for _, label in pos_entries]}")
    print(f"  Coverages  : {coverages}")
    print(f"  SV freqs   : {sv_freqs}")
    print(f"  Sizes      : {list(sizes.keys())}")
    print(f"  N samples  : {args.n_samples}")
    print(f"  Error rate : {args.error_rate}  K={args.k}")
    if not args.dry_run:
        print(f"  Stagger    : {args.stagger_secs}s between submissions")
    print()

    for i, (cov, freq, (pos, pos_label)) in enumerate(combinations, 1):
        freq_label  = int(round(freq * 100))
        config_text = generate_config(
            sv_type=args.sv_type, coverage=cov, sv_freq=freq,
            n_samples=args.n_samples, error_rate=args.error_rate,
            k=args.k, pos=pos, pos_label=pos_label, sizes=sizes,
        )
        config_name = (f"config_var_{args.sv_type.lower()}"
                       f"_n{args.n_samples}_f{freq_label}_cov{cov}"
                       f"_err{err_label}_k{args.k}_{pos_label}.sh")
        config_path = config_dir / config_name
        config_path.write_text(config_text)
        config_path.chmod(0o644)

        print(f"[{i}/{total}] {'─'*54}")
        print(f"  Position   : {pos:,} bp  →  {pos_label}")
        print(f"  Coverage   : {cov}x   SV freq: {freq}"
              f"  ({int(round(freq * args.n_samples))}/{args.n_samples} haplotypes)")
        print(f"  Sizes      : {', '.join(sizes.keys())}")
        print(f"  Config     : {config_path}")

        cmd = ["bash", str(pipeline), str(config_path)]
        print(f"  Command    : {' '.join(cmd)}")

        if args.dry_run:
            print("  [dry-run] Not submitting.\n")
        else:
            result = subprocess.run(cmd)
            if result.returncode != 0:
                print(f"  ERROR: pipeline returned {result.returncode}, aborting.")
                sys.exit(result.returncode)
            print()
            if i < total and args.stagger_secs > 0:
                print(f"  Waiting {args.stagger_secs}s …")
                time.sleep(args.stagger_secs)

    if not args.dry_run:
        print(f"\nAll {total} pipeline(s) submitted.")


if __name__ == "__main__":
    main()
