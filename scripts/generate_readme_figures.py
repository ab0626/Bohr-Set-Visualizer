"""
Generate PNG figures for README.md (run from repo root):
  python scripts/generate_readme_figures.py
"""
from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bohr_set import bohr_mask, bohr_size


def main() -> None:
    out_dir = os.path.join(ROOT, "docs", "images")
    os.makedirs(out_dir, exist_ok=True)

    N = 180
    thetas = [7, 31]
    eps_grid = np.linspace(0.02, 0.45, 120)
    sizes = [bohr_size(N, thetas, float(e)) for e in eps_grid]
    eps_cur = 0.08
    mask = bohr_mask(N, thetas, eps_cur)

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 12,
            "figure.facecolor": "white",
            "axes.facecolor": "#fafafa",
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=150)
    ax.plot(eps_grid, sizes, color="#1f77b4", lw=2, label=r"$|\Lambda_{\Theta,\varepsilon}|$")
    ax.axvline(eps_cur, color="#c44e52", ls="--", lw=1.2, label=r"example $\varepsilon$")
    ax.set_xlabel(r"radius $\varepsilon$ on $\mathbb{R}/\mathbb{Z}$")
    ax.set_ylabel(r"cardinality $|\Lambda_{\Theta,\varepsilon}|$")
    ax.set_title(rf"Bohr set size vs $\varepsilon$ in $\mathbb{{Z}}_{{{N}}}$, $\Theta={tuple(thetas)}$")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.35)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "bohr_cardinality_vs_epsilon.png"))
    plt.close(fig)

    fig2, ax2 = plt.subplots(figsize=(9.5, 2.8), dpi=150)
    green = "#2ca02c"
    gray = "#e8e8e8"
    colors = [green if mask[i] else gray for i in range(N)]
    ax2.bar(range(N), mask.astype(int), width=1.0, color=colors, edgecolor="none")
    ax2.set_xlim(0, N)
    ax2.set_ylim(0, 1.15)
    ax2.set_xlabel(r"$x \in \mathbb{Z}_N$")
    ax2.set_title(rf"Indicator $\mathbf{{1}}_\Lambda(x)$ for $\varepsilon={eps_cur}$, same $\Theta$")
    ax2.set_yticks([0, 1])
    fig2.tight_layout()
    fig2.savefig(os.path.join(out_dir, "bohr_indicator_zn.png"))
    plt.close(fig2)

    fig3, ax3 = plt.subplots(figsize=(7, 4), dpi=150)
    ratio_up = []
    ratio_mid = []
    sigmas = np.linspace(0.02, 0.15, 40)
    for sig in sigmas:
        e_in = eps_cur * (1 - sig)
        e_out = eps_cur * (1 + sig)
        if e_in <= 0 or e_out >= 0.5:
            ratio_up.append(np.nan)
            ratio_mid.append(np.nan)
            continue
        s0 = bohr_size(N, thetas, eps_cur)
        si = bohr_size(N, thetas, e_in)
        so = bohr_size(N, thetas, e_out)
        ratio_up.append(so / s0 if s0 else np.nan)
        ratio_mid.append(si / s0 if s0 else np.nan)
    ax3.plot(sigmas, ratio_up, label=r"$|\Lambda_{\varepsilon(1+\sigma)}|/|\Lambda_\varepsilon|$", color="#1f77b4")
    ax3.plot(sigmas, ratio_mid, label=r"$|\Lambda_{\varepsilon(1-\sigma)}|/|\Lambda_\varepsilon|$", color="#ff7f0e")
    ax3.axhline(1.0, color="#333", lw=0.8)
    ax3.set_xlabel(r"$\sigma$ (relative radius perturbation)")
    ax3.set_ylabel("size ratio")
    ax3.set_title(r"Regularity-style stability: sandwich $\Lambda_{\varepsilon(1\pm\sigma)}$ vs $\Lambda_\varepsilon$")
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.35)
    fig3.tight_layout()
    fig3.savefig(os.path.join(out_dir, "regularity_size_ratios.png"))
    plt.close(fig3)

    print("Wrote PNGs to", out_dir)


if __name__ == "__main__":
    main()
