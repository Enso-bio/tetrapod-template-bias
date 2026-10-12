"""Figure for the post hoc calibration of the per-gene tests (amendment E2b)."""
import importlib
import json

import numpy as np

from common import RESULTS

fig_mod = importlib.import_module('09_figures')
plt, save = fig_mod.plt, fig_mod.save


def main():
    z = np.load(RESULTS / 'posthoc_E2b_null_pvalues.npz', allow_pickle=True)
    info = json.loads((RESULTS / 'posthoc_E2b_pergene_calibration.json').read_text())
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for ax, key, p, title in [(axes[0], 'lambda_ML', z['p_ml'], 'λ estimated per gene (as in the tests)'),
                              (axes[1], 'lambda_fixed_null', z['p_fixed'], 'λ fixed at the null-model value')]:
        counts = (p < 0.05).sum(axis=1)
        r = info[key]
        ax.hist(counts, bins=40, color=fig_mod.AXIS, edgecolor='white', linewidth=0.8)
        ax.axvline(r['observed_nominal'], color=fig_mod.ORANGE, lw=2)
        ax.axvline(0.05 * info['n_tests'], color=fig_mod.MUTED, lw=1, ls='--')
        ax.annotate(f"observed = {r['observed_nominal']}\ncalibrated P = {r['p_calibrated_count']:.2f}",
                    (r['observed_nominal'], ax.get_ylim()[1] * 0.92), xytext=(6, 0), textcoords='offset points',
                    fontsize=7.5, color=fig_mod.INK2, va='top')
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(f"Genes with P < 0.05 (of {info['n_tests']:,})")
    axes[0].set_ylabel(f"Simulated null data sets (of {info['B']:,})")
    fig.tight_layout()
    save(fig, 'figS3_pergene_calibration')


if __name__ == '__main__':
    main()
