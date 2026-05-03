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

from bohr_set import bohr_mask, bohr_size, dual_group_fourier_magnitude


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

    # Illustrative quantitative arc (not literal theorem constants)
    fig4, ax4 = plt.subplots(figsize=(7.2, 4.2), dpi=150)
    n = np.logspace(1.5, 5, 80)
    doubly_log = 1.0 / np.maximum(np.log(np.maximum(np.log(n), np.e)), 1e-6) ** 0.4
    quasi = np.exp(-(np.log(n) ** 0.35))
    ax4.loglog(n, doubly_log / doubly_log.max(), label=r"illustrative $(\log\log |G|)^{-c}$ regime", color="#7f7f7f", lw=2)
    ax4.loglog(n, quasi / quasi.max(), label=r"illustrative $\exp(-(\log|G|)^{\Omega(1)})$ saving", color="#1f77b4", lw=2)
    ax4.set_xlabel(r"group size $|G|$ (illustrative)")
    ax4.set_ylabel(r"normalized density savings (schematic)")
    ax4.set_title("Quantitative arc (schematic — see paper for definitions)")
    ax4.legend(loc="lower left", fontsize=9)
    ax4.grid(True, which="both", alpha=0.3)
    fig4.tight_layout()
    fig4.savefig(os.path.join(out_dir, "quantitative_arc_schematic.png"))
    plt.close(fig4)

    # Fourier magnitudes of Bohr indicator (dual group of Z_N)
    f_ind = mask.astype(float)
    mag = dual_group_fourier_magnitude(f_ind)
    fig5, ax5 = plt.subplots(figsize=(8, 3.8), dpi=150)
    ax5.plot(np.arange(N), mag, color="#2ca02c", lw=1.0)
    for th in thetas:
        ax5.axvline(th, color="#d62728", ls=":", lw=1.0, alpha=0.8)
    ax5.set_xlabel(r"dual frequency index $r \in \mathbb{Z}_N$")
    ax5.set_ylabel(r"$|\widehat{\mathbf{1}_\Lambda}(r)|$")
    ax5.set_title(rf"Fourier magnitudes for same $\Theta={tuple(thetas)}$, $\varepsilon={eps_cur}$")
    fig5.tight_layout()
    fig5.savefig(os.path.join(out_dir, "fourier_bohr_indicator.png"))
    plt.close(fig5)

    # Non-abelian extension pipeline (conceptual)
    fig6, ax6 = plt.subplots(figsize=(8.5, 2.4), dpi=150)
    ax6.set_xlim(0, 10)
    ax6.set_ylim(0, 1)
    ax6.axis("off")
    boxes = [
        (0.2, 0.35, 1.6, 0.55, "finite group G"),
        (2.3, 0.35, 1.8, 0.55, "large abelian\nsubgroup H ⊆ G"),
        (4.7, 0.35, 1.9, 0.55, "corners theorem\non H"),
        (7.2, 0.35, 1.8, 0.55, "average / lift\nto G"),
    ]
    for x, y, w, h, txt in boxes:
        ax6.add_patch(plt.Rectangle((x, y), w, h, fill=False, lw=1.8, ec="#333"))
        ax6.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=11)
    for i in range(3):
        ax6.annotate(
            "",
            xy=(boxes[i + 1][0], 0.625),
            xytext=(boxes[i][0] + boxes[i][2], 0.625),
            arrowprops=dict(arrowstyle="->", lw=1.5, color="#555"),
        )
    ax6.set_title(r"Non-abelian extension (Sec. 1.1 — conceptual)", fontsize=12, pad=12)
    fig6.tight_layout()
    fig6.savefig(os.path.join(out_dir, "nonabelian_pipeline.png"))
    plt.close(fig6)

    print("Wrote PNGs to", out_dir)


if __name__ == "__main__":
    main()
