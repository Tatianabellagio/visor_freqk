#!/usr/bin/env python3
"""
compute_repeat_score.py
=======================
Compute the genome-wide repeat score for a registered position and update
the registry in-place.  Intended to run as a Slurm job, not on the login node.

Usage
-----
  python scripts/compute_repeat_score.py --rep-id rep1 --pos-bp 14725000
"""

import argparse
import csv
import fcntl
import sys
import tempfile
import os
from pathlib import Path

WORK     = Path("/home/tbellagio/scratch/visor_freqk")
FASTA    = WORK / "data/reference/Chr1.fa"
REGISTRY = WORK / "data/positions_registry.tsv"
K        = 31
WINDOW   = 10_000   # must match the window used in the analysis notebooks

# Canonical column order — must match launch_experiment_var.py
REGISTRY_FIELDS = ["rep_id", "pos_bp", "chrom", "global_repeat_score", "seed", "date_added"]


def build_repeat_set(ref_seq: str, k: int = K) -> frozenset:
    """frozenset of k-mers appearing ≥2× in ref_seq."""
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
    seq   = seq.replace("N", "")
    total = len(seq) - k + 1
    if total <= 0:
        return float("nan")
    return sum(1 for i in range(total) if seq[i : i + k] in repeat_set) / total


def update_registry(rep_id: str, score: float) -> None:
    """Update global_repeat_score for rep_id in the registry (file-locked)."""
    lock_path = REGISTRY.with_suffix(".lock")
    with open(lock_path, "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        try:
            if not REGISTRY.exists():
                sys.exit(f"ERROR: registry not found: {REGISTRY}")

            with open(REGISTRY, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                # Use canonical fields regardless of what the header says —
                # guards against column-count mismatches (e.g. missing 'seed' column)
                rows = [{k: v for k, v in row.items() if k is not None}
                        for row in reader]

            # Ensure all rows have every canonical field (fill missing with '')
            for row in rows:
                for field in REGISTRY_FIELDS:
                    row.setdefault(field, "")

            updated = False
            for row in rows:
                if row["rep_id"] == rep_id:
                    row["global_repeat_score"] = f"{score:.4f}"
                    updated = True
                    break

            if not updated:
                sys.exit(f"ERROR: rep_id '{rep_id}' not found in registry")

            # Write to a temp file in the same directory, then rename atomically.
            # This prevents data loss if writing fails mid-way.
            tmp_fd, tmp_path = tempfile.mkstemp(dir=REGISTRY.parent, suffix=".tmp")
            try:
                with os.fdopen(tmp_fd, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=REGISTRY_FIELDS,
                                            delimiter="\t", extrasaction="ignore")
                    writer.writeheader()
                    writer.writerows(rows)
                os.replace(tmp_path, REGISTRY)  # atomic on POSIX
            except Exception:
                os.unlink(tmp_path)
                raise

        finally:
            fcntl.flock(lf, fcntl.LOCK_UN)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rep-id", required=True, help="Registry rep_id (e.g. rep1)")
    p.add_argument("--pos-bp", required=True, type=int, help="0-based position in bp")
    args = p.parse_args()

    if not FASTA.exists():
        sys.exit(f"ERROR: reference FASTA not found: {FASTA}")

    print(f"Loading reference …", end=" ", flush=True)
    lines   = FASTA.read_text().splitlines()
    ref_seq = "".join(l for l in lines if not l.startswith(">")).upper()
    print(f"done  ({len(ref_seq):,} bp)")

    print(f"Building genome-wide k-mer repeat set …", end=" ", flush=True)
    repeat_set = build_repeat_set(ref_seq)
    print(f"done  ({len(repeat_set):,} multi-copy k-mers)")

    rs = global_repeat_score(ref_seq[args.pos_bp : args.pos_bp + WINDOW], repeat_set)
    print(f"global_repeat_score({args.rep_id}, pos={args.pos_bp:,}) = {rs:.4f}")

    update_registry(args.rep_id, rs)
    print(f"Registry updated: {REGISTRY}")


if __name__ == "__main__":
    main()
