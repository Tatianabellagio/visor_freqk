"""Slide plot — R² of true vs estimated AF as a function of SV size.

Uses results_rep.csv (k=31, ~31 randomly picked rep positions, mostly
low-repeat, ~100% detection) so the only thing varying is SV size.

One dot per (size × coverage) bin with 95% bootstrap CIs (1000 resamples).
Viridis palette for coverage, clean style.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

REP_CSV  = '../summaries/results_rep.csv'
PLOT_DIR = '../plots'
os.makedirs(PLOT_DIR, exist_ok=True)

SIZE_ORDER = ['100bp', '500bp', '1kb', '5kb', '10kb']
COV_ORDER  = [10, 20, 50]
_viridis = cm.get_cmap('viridis')
COV_COLORS = {10: _viridis(0.0), 20: _viridis(0.5), 50: _viridis(1.0)}

N_BOOT = 1000
RNG    = np.random.default_rng(0)

df = pd.read_csv(REP_CSV)
df = df[df['k_label'] == 'k31'].copy()
df['detected'] = df['af_alt'].notna()
df['af_alt_imputed'] = df['af_alt'].fillna(0.0)


def r2(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    if len(y_true) < 2 or np.var(y_true) == 0:
        return np.nan
    ss_res = ((y_true - y_pred) ** 2).sum()
    ss_tot = ((y_true - y_true.mean()) ** 2).sum()
    return max(0.0, 1.0 - ss_res / ss_tot)


def bootstrap_r2_ci(y_true, y_pred, n_boot=N_BOOT, rng=RNG):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    n = len(y_true)
    if n < 2:
        return (np.nan, np.nan)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boots[i] = r2(y_true[idx], y_pred[idx])
    return (np.nanpercentile(boots, 2.5), np.nanpercentile(boots, 97.5))


rows = []
for cov in COV_ORDER:
    for s in SIZE_ORDER:
        sub = df[(df['coverage'] == cov) & (df['size'] == s)]
        point = r2(sub['freq_nominal'], sub['af_alt_imputed'])
        lo, hi = bootstrap_r2_ci(sub['freq_nominal'].values,
                                 sub['af_alt_imputed'].values)
        rows.append({
            'coverage': cov,
            'size': s,
            'n_total': len(sub),
            'n_detected': int(sub['detected'].sum()),
            'detection_rate': sub['detected'].mean() if len(sub) else np.nan,
            'r2': point,
            'r2_lo': lo,
            'r2_hi': hi,
        })
summary = pd.DataFrame(rows)
print(summary.to_string(index=False))

x = np.arange(len(SIZE_ORDER))
X_JITTER = {10: -0.12, 20: 0.0, 50: 0.12}

fig, ax = plt.subplots(figsize=(8, 4.6), constrained_layout=True)
ax.set_facecolor('white')

for cov in COV_ORDER:
    s = summary[summary['coverage'] == cov].set_index('size').reindex(SIZE_ORDER)
    y     = s['r2'].values
    y_lo  = s['r2_lo'].values
    y_hi  = s['r2_hi'].values
    xx    = x + X_JITTER[cov]
    yerr  = np.vstack([y - y_lo, y_hi - y])

    ax.errorbar(xx, y, yerr=yerr, fmt='o',
                color=COV_COLORS[cov],
                ms=10, mec='white', mew=1.0,
                elinewidth=1.2, capsize=3.5, capthick=1.0,
                label=f'{cov}×', zorder=3)

ax.set_xticks(x)
ax.set_xticklabels(SIZE_ORDER, fontsize=10)
ax.set_xlabel('Deletion size', fontsize=11)
ax.set_ylabel('R²  (true AF  vs  estimated AF)', fontsize=11)
ax.set_ylim(0, 1.05)
ax.grid(axis='y', alpha=0.35, linewidth=0.6)
ax.grid(axis='x', visible=False)
ax.set_axisbelow(True)
for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)
ax.tick_params(labelsize=9)

ax.legend(title='Coverage', fontsize=9, title_fontsize=9,
          loc='lower right', frameon=False)

n_pos = df['rep_id'].nunique()
ax.set_title(f'Effect of SV size on AF estimation  —  '
             f'{n_pos} random rep positions (k=31)',
             fontsize=11)

plt.savefig(f'{PLOT_DIR}/slide_sv_size_r2.pdf', bbox_inches='tight')
plt.savefig(f'{PLOT_DIR}/slide_sv_size_r2.png', bbox_inches='tight', dpi=300)
print(f'\nWrote {PLOT_DIR}/slide_sv_size_r2.pdf and .png')
