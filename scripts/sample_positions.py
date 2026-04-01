#!/usr/bin/env python3
"""
sample_positions.py
===================
Randomly sample N genomic positions for SV deletion experiments.

Hard constraints enforced on every sampled position:
  1. pos >= END_BUFFER from chromosome start
  2. pos + MAX_DEL_SIZE + END_BUFFER <= chrom_len
     (the entire largest deletion + one buffer must fit before the chromosome end)
  3. No N-bases anywhere in [pos, pos + MAX_DEL_SIZE)
  4. |pos - existing_pos| >= MIN_SEP for all already-registered positions
     (prevents SV windows from overlapping and keeps positions well-spaced)

Each accepted position receives the next available rep_id (rep1, rep2, …)
and is appended to data/positions_registry.tsv.  The repeat score stored in
the registry is informational only — it is NOT used to filter positions.

Usage
-----
  # Sample 20 new positions and write them to the registry
  python scripts/sample_positions.py --n 20 --seed 42

  # Preview without modifying the registry
  python scripts/sample_positions.py --n 20 --seed 42 --dry-run

  # Pass the printed --rep-ids line directly to launch_experiment_var.py
  python scripts/sample_positions.py --n 5 --seed 42 | tail -1
"""

import argparse
import csv
import datetime
import random
import sys
from pathlib import Path

# ── Paths (relative to this script) ──────────────────────────────────────────
WORK     = Path(__file__).resolve().parent.parent
FASTA    = WORK / "data/reference/Chr1.fa"
REGISTRY = WORK / "data/positions_registry.tsv"

# ── Default constraint parameters ────────────────────────────────────────────
MAX_DEL_SIZE = 10_000    # largest deletion size in the pipeline (bp)
END_BUFFER   = 100_000   # minimum distance from chromosome ends (bp)
MIN_SEP      = 100_000   # minimum distance between any two registered positions (bp)
K            = 31        # k-mer size (matches the rest of the pipeline)
RS_WINDOW    = 10_000    # window used for global repeat score (informational)

REGISTRY_FIELDS = ["rep_id", "pos_bp", "chrom", "global_repeat_score", "date_added"]


# ── Repeat score helpers (informational — not used for filtering) ─────────────

def build_repeat_set(ref_seq: str, k: int = K) -> frozenset:
    """Return frozenset of k-mers that appear ≥2× anywhere in ref_seq."""
    seen_once: set = set()
    repeated:  set = set()
    seq_clean = ref_seq.replace("N", "")
    for i in range(len(seq_clean) - k + 1):
        km = seq_clean[i : i + k]
        if km in seen_once:
            repeated.add(km)
        else:
            seen_once.add(km)
    del seen_once
    return frozenset(repeated)


def global_repeat_score(seq: str, repeat_set: frozenset, k: int = K) -> float:
    """Fraction of k-mers in seq that appear ≥2× genome-wide."""
    seq   = seq.replace("N", "")
    total = len(seq) - k + 1
    if total <= 0:
        return float("nan")
    return sum(1 for i in range(total) if seq[i : i + k] in repeat_set) / total


# ── Registry I/O ──────────────────────────────────────────────────────────────

def load_registry() -> list:
    if not REGISTRY.exists() or REGISTRY.stat().st_size == 0:
        return []
    with open(REGISTRY) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def next_rep_id(rows: list) -> str:
    """Return the next unused rep_id (rep1, rep2, …)."""
    used = {int(r["rep_id"].replace("rep", "")) for r in rows}
    n = 1
    while n in used:
        n += 1
    return f"rep{n}"


def append_to_registry(rows: list) -> None:
    write_header = not REGISTRY.exists() or REGISTRY.stat().st_size == 0
    REGISTRY.parent.mkdir(parents=True, exist_ok=True)
    with open(REGISTRY, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=REGISTRY_FIELDS, delimiter="\t")
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--n",            required=True, type=int,
                   help="Number of new positions to sample")
    p.add_argument("--seed",         required=True, type=int,
                   help="Random seed (for reproducibility)")
    p.add_argument("--end-buffer",   default=END_BUFFER,   type=int,
                   help=f"Min distance from chromosome ends in bp (default: {END_BUFFER})")
    p.add_argument("--min-sep",      default=MIN_SEP,      type=int,
                   help=f"Min distance between any two positions in bp (default: {MIN_SEP})")
    p.add_argument("--max-del-size", default=MAX_DEL_SIZE, type=int,
                   help=f"Largest deletion size in the pipeline in bp (default: {MAX_DEL_SIZE})")
    p.add_argument("--dry-run",      action="store_true",
                   help="Print sampled positions but do not write to the registry")
    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if not FASTA.exists():
        sys.exit(f"ERROR: reference FASTA not found: {FASTA}")

    # ── Load reference ────────────────────────────────────────────────────────
    print(f"Loading reference … ", end="", flush=True)
    lines     = FASTA.read_text().splitlines()
    ref_seq   = "".join(l for l in lines if not l.startswith(">")).upper()
    chrom_len = len(ref_seq)
    print(f"done  ({chrom_len:,} bp)")

    # ── Build genome-wide repeat set (informational) ──────────────────────────
    print("Building genome-wide k-mer repeat set … ", end="", flush=True)
    repeat_set = build_repeat_set(ref_seq)
    print(f"done  ({len(repeat_set):,} multi-copy k-mers)")

    # ── Load existing registry ────────────────────────────────────────────────
    existing     = load_registry()
    existing_pos = {int(r["pos_bp"]) for r in existing}
    print(f"\nExisting registry entries : {len(existing)}")

    # ── Validate constraint range ─────────────────────────────────────────────
    lo = args.end_buffer
    hi = chrom_len - args.max_del_size - args.end_buffer
    if lo >= hi:
        sys.exit(f"ERROR: valid range [{lo}, {hi}) is empty — "
                 f"reduce --end-buffer or --max-del-size")

    print(f"Valid sampling range      : [{lo:,}, {hi:,}) bp  "
          f"({(hi-lo)/1e6:.1f} Mb available)\n")

    # ── Sample ────────────────────────────────────────────────────────────────
    # Step through candidates at 1 kb resolution (finer = slower but more choices)
    STEP = 1_000
    rng  = random.Random(args.seed)
    candidates = list(range(lo, hi, STEP))
    rng.shuffle(candidates)

    accepted:  list = []
    n_n_fail:  int  = 0
    n_sep_fail: int = 0

    for pos in candidates:
        if len(accepted) == args.n:
            break

        # Hard constraint 1: no N-bases in the deletion window
        if "N" in ref_seq[pos : pos + args.max_del_size]:
            n_n_fail += 1
            continue

        # Hard constraint 2: minimum separation from existing + already accepted
        all_taken = existing_pos | {p for p in accepted}
        if any(abs(pos - p) < args.min_sep for p in all_taken):
            n_sep_fail += 1
            continue

        accepted.append(pos)

    print(f"Sampling result: {len(accepted)} / {args.n} requested")
    print(f"  Rejected — N-bases   : {n_n_fail}")
    print(f"  Rejected — too close : {n_sep_fail}")

    if len(accepted) < args.n:
        print(f"\nWARNING: only {len(accepted)} positions satisfied all constraints.",
              file=sys.stderr)
        print(f"  Try: --min-sep {args.min_sep // 2}  or  --end-buffer {args.end_buffer // 2}",
              file=sys.stderr)

    # ── Build registry rows ───────────────────────────────────────────────────
    print()
    new_rows = []
    all_rows = existing[:]   # used for next_rep_id bookkeeping
    for pos in accepted:
        rep_id = next_rep_id(all_rows + new_rows)
        rs     = global_repeat_score(ref_seq[pos : pos + RS_WINDOW], repeat_set)
        row    = {
            "rep_id":              rep_id,
            "pos_bp":              str(pos),
            "chrom":               "Chr1",
            "global_repeat_score": f"{rs:.4f}",
            "date_added":          datetime.date.today().isoformat(),
        }
        new_rows.append(row)
        print(f"  {rep_id:8s}  pos={pos:>12,}  ({pos/1e6:6.3f} Mb)  "
              f"repeat_score={rs:.3f}")

    # ── Print launch command ──────────────────────────────────────────────────
    print()
    print("--rep-ids " + " ".join(r["rep_id"] for r in new_rows))

    # ── Write registry ────────────────────────────────────────────────────────
    if args.dry_run:
        print("\n[dry-run] Registry not modified.")
    else:
        append_to_registry(new_rows)
        print(f"\nAppended {len(new_rows)} entries to: {REGISTRY}")


if __name__ == "__main__":
    main()
