"""Slide plot — R² heatmap by repeat level × coverage.

Replaces the dot-plot attempt (too few positions per bin for a strip plot).
Same viridis color scheme as the SV-size plot.

Data: results_var_selectedpos_w_repeatscore.csv (k=31). End-to-end R² with
non-detected imputed as AF=0. Uses Pearson² (np.corrcoef) to match the
original heatmap in presentation_plots.ipynb.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SELPOS_CSV = '../summaries/results_var_selectedpos_w_repeatscore.csv'
PLOT_DIR   = '../plots'
os.makedirs(PLOT_DIR, exist_ok=True)

REGION_ORDER = ['low repeat', 'moderate repeat', 'high repeat']
COV_ORDER    = [10, 20, 50]

df = pd.read_csv(SELPOS_CSV)
df = df[df['k_label'] == 'k31'].copy()
df['af_alt_imputed'] = df['af_alt'].fillna(0.0)

r2_mat = np.full((len(REGION_ORDER), len(COV_ORDER)), np.nan)
for i, rg in enumerate(REGION_ORDER):
    for j, cov in enumerate(COV_ORDER):
        sub = df[(df['region_type'] == rg) & (df['coverage'] == cov)]
        if len(sub) >= 2 and sub['freq_nominal'].nunique() >= 2:
            corr = np.corrcoef(sub['freq_nominal'], sub['af_alt_imputed'])[0, 1]
            if not np.isnan(corr):
                r2_mat[i, j] = corr ** 2

print('R² matrix (rows=region, cols=coverage):')
print(pd.DataFrame(r2_mat, index=REGION_ORDER,
                   columns=[f'{c}×' for c in COV_ORDER]).round(3))

fig, ax = plt.subplots(figsize=(5.6, 4.2), constrained_layout=True)

im = ax.imshow(r2_mat, cmap='viridis', vmin=0, vmax=1, aspect='auto')

# Cell annotations (white text on dark cells, black on light)
for i in range(len(REGION_ORDER)):
    for j in range(len(COV_ORDER)):
        val = r2_mat[i, j]
        if np.isnan(val):
            continue
        color = 'white' if val < 0.55 else 'black'
        ax.text(j, i, f'{val:.3f}',
                ha='center', va='center',
                color=color, fontsize=14)

ax.set_xticks(range(len(COV_ORDER)))
ax.set_xticklabels([f'{c}×' for c in COV_ORDER], fontsize=11)
ax.set_yticks(range(len(REGION_ORDER)))
ax.set_yticklabels(REGION_ORDER, fontsize=11)
ax.set_xlabel('Pool coverage', fontsize=12)
ax.tick_params(length=0)

for spine in ax.spines.values():
    spine.set_visible(False)

cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.03)
cbar.set_label('R²  (true AF  vs  estimated AF)', fontsize=10)
cbar.ax.tick_params(labelsize=9)
cbar.outline.set_visible(False)

plt.savefig(f'{PLOT_DIR}/slide_repeat_r2.pdf', bbox_inches='tight')
plt.savefig(f'{PLOT_DIR}/slide_repeat_r2.png', bbox_inches='tight', dpi=300)
print(f'\nWrote {PLOT_DIR}/slide_repeat_r2.pdf and .png')
