"""
06_eval_cactus_em.py
Compare cactus_em predicted alt_freq at the deletion locus to truth (= SV_FREQ).

Loads the per-record TSV produced by per_sample_driver and pulls the records
within ±50 bp of the simulated deletion start. Reports per-mode (global,
window) the predicted alt_freq, the truth, and the absolute error.

Usage:
    python 06_eval_cactus_em.py --config config_cactus_em.sh --size 1kb
"""
from __future__ import annotations
import argparse, os, re, json
from pathlib import Path
import pandas as pd
import numpy as np


def parse_config(path):
    """Parse a bash-style config file into a dict (var=value lines)."""
    cfg = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r'^([A-Za-z_][A-Za-z_0-9]*)="?([^"#]*?)"?(?:\s*#.*)?$', line)
            if m:
                k, v = m.group(1), m.group(2).strip()
                v = re.sub(r'\$\{([A-Za-z_][A-Za-z_0-9]*)\}',
                           lambda m: cfg.get(m.group(1), m.group(0)), v)
                cfg[k] = v
    return cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config_cactus_em.sh")
    ap.add_argument("--size", default="1kb",
                    help="deletion size key from DEL_SIZES (default 1kb)")
    ap.add_argument("--window-bp", type=int, default=200,
                    help="report alt_freq for records within ±this many bp "
                         "of the deletion start (default 200)")
    args = ap.parse_args()

    cfg = parse_config(os.path.join(os.path.dirname(__file__), args.config))
    chrom = cfg["CHROM"]
    pos = int(cfg["SV_START_0"])
    sv_freq_truth = float(cfg["SV_FREQ"])
    n_samples = int(cfg["N_SAMPLES"])
    n_sv = round(n_samples * sv_freq_truth)
    truth = n_sv / n_samples
    coverage = cfg.get("COVERAGE", "?")
    err_label = cfg.get("ERROR_RATE", "0.001").split(".")[-1] or "0"
    freq_label = round(sv_freq_truth * 100)

    run_tag = f"n{n_samples}_f{freq_label}_err{err_label}"
    results_root = Path(cfg["RESULTS"])
    results_dir = results_root / f"cov{coverage}" / f"var_del_{args.size}_{run_tag}"

    print(f"=== {cfg['POS_LABEL']} / {args.size} / cov{coverage} ===")
    print(f"  truth: alt_freq = {n_sv}/{n_samples} = {truth:.4f}")
    print(f"  results dir: {results_dir}")

    rows = []
    for mode in ["global", "window"]:
        tsv = results_dir / f"cactus_em_{mode}.tsv"
        if not tsv.exists():
            print(f"\n  [{mode}] not found: {tsv}")
            continue
        df = pd.read_csv(tsv, sep="\t")
        # Records near the deletion
        m = (df.chrom == chrom) & (df.pos >= pos - args.window_bp) & \
            (df.pos <= pos + args.window_bp)
        near = df[m].copy()
        if near.empty:
            print(f"\n  [{mode}] no cactus VCF records within ±{args.window_bp} bp "
                  f"of {chrom}:{pos}")
            continue
        # Mean alt_freq over near-deletion records (excluding any duplicates)
        af_mean = near.alt_freq.mean()
        af_max = near.alt_freq.max()
        af_min = near.alt_freq.min()
        n_records = len(near)
        err = af_mean - truth
        print(f"\n  [{mode}]  n_records={n_records}  "
              f"af_mean={af_mean:.4f} (truth={truth:.4f}, err={err:+.4f}), "
              f"af range=[{af_min:.4f}, {af_max:.4f}]")
        rows.append({
            "mode": mode, "n_records": n_records,
            "af_mean": af_mean, "af_min": af_min, "af_max": af_max,
            "truth": truth, "abs_err": abs(err),
        })

    if rows:
        out = pd.DataFrame(rows)
        out_path = results_dir / "cactus_em_eval.tsv"
        out.to_csv(out_path, sep="\t", index=False)
        print(f"\nWrote summary → {out_path}")


if __name__ == "__main__":
    main()
