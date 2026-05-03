# Bohr Set Visualizer

## Description

The **Bohr Set Visualizer** is a **Streamlit** research app plus Python numerics for **Bohr sets** $\Lambda_{\Theta,\varepsilon}$ inside $\mathbb{Z}_N$, motivated by the corners theorem proof for general abelian groups in Jaber–Liu–Lovett–Ostuni–Sawhney ([arXiv:2504.07006](https://arxiv.org/abs/2504.07006)). It turns dense definitions (regular Bohr sets, relative sifting, $(d,\eta)$ chains, $\ell_1$-spreadness, dual Fourier spectrum, and qualitative links to Sec. 7 and Sec. 8) into **plots and sliders** you can change live.

**Scope.** The implementation fixes $G=\mathbb{Z}_N$ with explicit characters $\chi_\theta(x)=e^{2\pi i \theta x/N}$; the README also summarizes how such abelian inputs relate to **extensions beyond abelian groups** and **coloring** results stated in the paper.

**GitHub “About” (one line, copy-paste):**  
*Interactive Streamlit tool: Bohr sets in ℤ_N, regularity, sub-Bohr sifting, increment chains, ℓ₁-spread and Fourier views, pair-health proxy, and coloring pedagogy—aligned with arXiv:2504.07006 (quasipolynomial corners).*

**Math in this file** uses GitHub’s `$…$` / `$$…$$` syntax so formulas render in the GitHub web UI.

---

## Paper reference (PDF)

| Copy | Path |
|------|------|
| **Current (Cursor workspace)** | `c:\Users\adith\AppData\Roaming\Cursor\User\workspaceStorage\e0f943467b0a7cf1dc5dcf6b6f4194bc\pdfs\b0e9945c-06fe-4c99-89b3-c9e09200ecec\2504.07006v2 (1).pdf` |
| Other workspace copies | `…\0b984680-ef1a-48bb-bf61-1afb12277129\…`, `…\ba6baa42-86cb-4ab0-8c4b-3d1c72ba8949\…`, `…\b6b3e01f-9e8f-4de2-bc0a-3587d1e9a389\…`, `…\32e81826-9a93-4e51-9a32-2e7cd4f2312a\…` |

All refer to the same arXiv revision: **2504.07006v2** ([abstract](https://arxiv.org/abs/2504.07006)).

---

## Section 6 structure (paper)

The corners proof for **general abelian groups** replaces subspaces with Bohr sets; Section 6 develops Bohr sets, spreadness, and pseudorandomization. Dependency outline:

```mermaid
flowchart TB
  S6["Sec. 6 — Bohr sets, algebraic spreadness, pseudorandomization"]
  S6 --> A["6.1 Bohr sets"]
  A --> D61["Def 6.1 Λ_{Θ,ε}, dilation cΛ"]
  A --> L62["Lemma 6.2 |Λ| ≥ ε^d|G|"]
  A --> D63["Def 6.3 regular Bohr sets"]
  A --> L64["Lemma 6.4 ∃α Λ regular"]
  A --> L65["Lemma 6.5 shift invariance"]
  A --> D66["Def 6.6 (d,η)-small / exact sequences"]
  S6 --> B["Later in Sec. 6: grid norms, spreadness, …"]
  B --> D613["Def 6.13 ℓ₁-spreadness"]
  B --> D615["Def 6.15 algebraic spreadness"]
```

---

## Non-abelian groups (Section 1.1): where $\mathbb{Z}_N$ fits

The main theorem is stated for **finite abelian** $G$, but the introduction explains a striking extension: **quasipolynomial** corner bounds for abelian groups imply analogous bounds for **every finite group** $G$, including non-abelian ones. The pipeline (Fox-type argument, as discussed in Sec. 1.1) is:

1. **Find a large abelian subgroup** $H \le G$ (every large finite group contains a proportional abelian piece).
2. **Run the abelian theory** on $H$ — Bohr sets, relative sifting, grid norms, and the diagonal “$\ell_1$-spread” phenomena you simulate here are native to that step.
3. **Average / lift** from $H$ back to $G$ so that dense corner-free structure in $G$ cannot evade the abelian obstruction.

**Implementation link.** This repository only draws **$G = \mathbb{Z}_N$**, but that model is exactly the kind of **abelian building block** that feeds the general proof: Bohr sets are the right **containers** inside $H$; the lift to arbitrary $G$ is conceptual, not a different codebase.

---

## Coloring and the “tower” warning (Section 8)

The paper’s **coloring** corollary (Corollary 1.2, Sec. 8) addresses **$G \times G \times G$** with about **$L \asymp \log\log\log|G|$** colors, forcing a **monochromatic corner**. The authors stress that with **only** doubly-logarithmic density savings (as in Shkredov’s classical regime), one would face **tower-type** losses in how many colors can be handled — i.e. the quantitative strength of the new bound is what makes a **polylogarithmic-in-log** color count possible.

**Visualizer add-on.** Under **Tab 1 → “Coloring strip”**, the app assigns an **$L$-coloring** of $\mathbb{Z}_N$ (either $x \bmod L$ or i.i.d. random) and reports how many colors appear on the current Bohr set $\Lambda$, whether $\Lambda$ is monochromatic, and the exact probability $L^{1-|\Lambda|}$ that a **fully random** coloring makes $\Lambda$ monochromatic. That toy is **1-dimensional**; it is meant only to build intuition for why **better-than–$\log\log$** density decay matters for multicolor statements.

---

## Quantitative arc: why 2025 is a breakthrough

| Aspect | Shkredov (2006) corner bound | Jaber et al. (2025) |
|--------|------------------------------|----------------------|
| **Strength** | Doubly logarithmic: density $\ll (\log\log |G|)^{-c}$ | **Quasipolynomial**: $|A| \le |G|^2 \exp\bigl(-(\log|G|)^{\Omega(1)}\bigr)$ |
| **Core mechanism** | $L^2$ / energy-style increment | **Relative sifting**, **Bohr containers**, **Gowers grid norms**, **asymmetric** treatment of the diagonal $D$ ($\ell_1$-spread) |
| **Diagonal** | (Symmetric finite-field intuition) | **“Diagonal drop”**: $D$ is not assumed algebraically spread like $X,Y$ |

---

## How this repo sits in a larger “suite” (conceptual map)

| Focus | Role in the story |
|-------|-------------------|
| **Behrend-type constructions** | **Lower bounds** — how dense a corner-free set can still be. |
| **Grid norms / Von Neumann lemmas** | **Triggers** — when lack of uniformity forces structured increments. |
| **Bohr sets (this repo)** | **Containers** — soft structured neighborhoods inside general abelian groups. |
| **Exactly-$N$ / communication** | **Applications** — complexity consequences stated in the paper’s abstract. |

---

## Definitions from the paper (Section 6.1)

**Definition 6.1 (Bohr set).** Let $\varepsilon \in \mathbb{R}^+$, let $G$ be a finite abelian group, and $\Theta = (\Theta_1,\ldots,\Theta_d)$ with $\Theta_i \in \widehat{G}$ homomorphisms $G \to \mathbb{R}/\mathbb{Z}$. The Bohr set is

$$
\Lambda = \Lambda_{\Theta,\varepsilon} = \bigcap_{i=1}^{d} \Big\{ x \in G : \|\Theta_i(x)\|_{\mathbb{R}/\mathbb{Z}} \le \varepsilon \Big\},
$$

where $\|x\|_{\mathbb{R}/\mathbb{Z}} = \min_{z \in \mathbb{Z}} |x - z|$. For $c > 0$, **dilation** is $c\Lambda_{\Theta,\varepsilon} = \Lambda_{\Theta,\,c\varepsilon}$. The paper calls $d$ the **dimension** and $\varepsilon$ the **radius** $\nu(\Lambda)$.

**Lemma 6.2.** If $\Lambda = \Lambda_{\Theta,\varepsilon}$ has dimension $d$ and radius $\varepsilon$, then $|\Lambda| \ge \varepsilon^d |G|$.

**Definition 6.3 (Regular).** A Bohr set $\Lambda = \Lambda_{\Theta,\varepsilon}$ of dimension $d$ is **regular** if for all $|c| \le 1/(100d)$,

$$
1 - 100d|c| \;\le\; \frac{|(1+c)\Lambda|}{|\Lambda|} \;\le\; 1 + 100d|c|.
$$

**Lemma 6.5.** Let $f$ be $1$-bounded and $\Lambda$ a **regular** Bohr set of dimension $d$. If $|c| \le 1/(100d)$ and $n' \in c\Lambda$, then

$$
\mathbb{E}_{n \in \Lambda} f(n) \;=\; \mathbb{E}_{n \in \Lambda} f(n + n') \;+\; O(cd).
$$

This captures **approximate shift-invariance** when the shift lies in a smaller dilate $c\Lambda$.

---

## Cyclic group model ($G = \mathbb{Z}_N$)

Characters are $\chi_\theta(x) = e^{2\pi i \theta x/N}$ with $\theta \in \{0,\ldots,N-1\}$. This repo implements

$$
x \in \Lambda_{\Theta,\varepsilon}
\quad\Longleftrightarrow\quad
\big\| (\theta x)/N \big\|_{\mathbb{R}/\mathbb{Z}} < \varepsilon \quad \forall\, \theta \in \Theta,
$$

with either an $\ell_\infty$ or $\ell_1$ combination over coordinates (see `bohr_set.py`). That matches Definition 6.1 once $\Theta_i(x)$ is read mod $1$.

---

## Research simulation: mapping to the paper

The Streamlit app is not only drawing Bohr sets; it **simulates mechanisms** from the proof: relative sifting inside a sparse container, nested Bohr “zoom” chains, and **Definition 6.13** $\ell_1$-spreadness.

### Sub-Bohr sifting and relative density (Sections 2.2, 3.5)

**Problem (informal).** Classical sifting bounds degrade with the density of a majorant in all of $G$. If $A$ is tiny in $G$, crude estimates lose useful factors.

**What the tool reports.** For an outer Bohr container $B_1$ and sparse $A \subset B_1$, compare **global** density $|A|/|G|$ to **relative** density $|A|/|B_1|$, and track how much of $A$ falls into a smaller inner Bohr neighborhood $B_2$. That is the right mental model for **relative sifting**: locating a sub-instance where a set that is sparse globally is still **dense enough relative to $B_1$** (and refinable inside $B_2$) to feed a density-increment argument.

### Increment chains and $(d,\eta)$-small sequences (Sections 6.1, 6.6)

**Logic.** Definition 6.6 fixes a **$(d,\eta)$-small sequence**: same frequency set (rank $d$), nested Bohr sets with radii shrinking so $\nu(B_{i+1})/\nu(B_i) \le \eta$.

**What the tool plots.** Radii $\varepsilon_i = \varepsilon_0 \eta^i$ at fixed $\Theta$, with cardinalities $|\Lambda_{\Theta,\varepsilon_i}|$. Dimension stays fixed while the radius collapses — the same “zoom” iteration that supports the density-increment loop (as in the proof structure toward results such as Theorem 7.9 in the paper).

### $\ell_1$-spread heatmap (Sections 6.1, 6.3; Definition 6.13)

**Breakthrough (informal).** The diagonal-type object $D$ is not treated with the same **algebraic spreadness** used for the other sides; the authors use **$\ell_1$-spreadness** (Definition 6.13) as a workable substitute.

**What the tool computes.** For $f = \mathbf{1}_A$ (extended by zero off $A$), it evaluates the deviation averaged over $x \sim B_1$:

$$
\mathbb{E}_{x \sim B_1} \Big| \mathbb{E}_{y \sim B_2}[f(x+y)] - \mathbb{E}[f] \Big|,
$$

with $\mathbb{E}[f] = \mathbb{E}_{x \sim B_1}[f(x)]$, and compares it to $\varepsilon_{\mathrm{tolerance}} \cdot \mathbb{E}[f]$. The bar plot over $x \in \mathbb{Z}_N$ shows **where** $f$ fails to look translation-averaged along $B_2$ — empirical feedback for the **log-potential** style analysis that rules out persistent “dips below the mean” on bad pieces of $D$.

### Lab components vs paper concepts

| Lab component | Paper concept | Purpose |
|---------------|---------------|---------|
| **Sub-Bohr sifting** | Relative sifting (see Sec. 2.2, Sec. 3.5) | Find structure inside a sparse majorant $B_1$ by measuring density relative to $B_1$ and refinement inside $B_2$. |
| **Increment chain** | $(d,\eta)$-small / exact sequences (Def. 6.6) | Model the iterative Bohr “zoom” where rank is fixed and radius shrinks. |
| **$\ell_1$-spread analysis** | Definition 6.13 | Test stability of the diagonal-type direction when full algebraic spreadness is unavailable. |
| **Fourier spectrum (dual group)** | Sec. 6.2, Appendix A (almost periodicity) | $|\widehat{\mathbb{1}_\Lambda}(r)|$ over $r \in \mathbb{Z}_N$; compare energy at low modes and at the fixed $\Theta$ to see sparse frequency support. |
| **Balanced strip** | Balanced $f - \mathbb{E}f$ (Sec. 2, Sec. 4) | Visualize $\mathbb{1}_\Lambda - \alpha$ (global or on $B_1$) to see mean-zero “clumpiness” used in Gowers-type arguments. |
| **Pair health map (proxy)** | Sec. 7.1 (well-conditioned pairs), Lemma 7.5 (not poor) | Coarse 0–3 score from slice-stability of $X,Y$ and a $D$-tube check; **not** a full $(B_i,B_8,B_9,K,K)$ grid norm. |

### Presentation note ($\ell_1$-spread and the heatmap)

> **Why check $\ell_1$-spreadness on a heatmap?** In Section 6.3 the authors use a **log-potential** argument to show that if a set is not $\ell_1$-spread, it can be **partitioned** into pieces that are (cf. Lemma 6.21 and the surrounding recursion). The heatmap marks residues $x$ where the pointwise term $\big|\mathbb{E}_{y \sim B_2}[f(x+y)] - \mathbb{E}[f]\big|$ is large — i.e. where the obstructions to spreadness live and where “bad pieces” would be carved out in that partitioning step.

---

## Figures (generated)

These PNGs live under `docs/images/` and are referenced with **relative** paths so they render on GitHub.

**Cardinality vs radius** $|\Lambda_{\Theta,\varepsilon}|$ for a fixed $\Theta \subset \mathbb{Z}_{180}$:

![|Λ| vs ε](docs/images/bohr_cardinality_vs_epsilon.png)

**Membership** in $\mathbb{Z}_N$ (same parameters):

![indicator on Z_N](docs/images/bohr_indicator_zn.png)

**Regularity-style ratios** (sandwich $\Lambda_{\varepsilon(1\pm\sigma)}$ vs $\Lambda_\varepsilon$): compare qualitatively to Definition 6.3 — bounded relative change in radius should imply bounded relative change in size for **regular** Bohr sets.

![regularity ratios](docs/images/regularity_size_ratios.png)

Regenerate figures after changing defaults:

```powershell
python scripts/generate_readme_figures.py
```

---

## Code map

| File | Role |
|------|------|
| `bohr_set.py` | Torus norm, $\Lambda_{\Theta,\varepsilon}$, regularity sandwich, random walk / TV, sub-Bohr sifting, Bohr chains, $\ell_1$-spread, **dual Fourier magnitudes**, **balanced function**, **pair-health proxy matrix** |
| `app.py` | Streamlit UI: nested $B_2 \subset B_1$, sparse $A$, Fourier spectrum, balanced strip, **coloring strip (Sec. 8 pedagogy)**, pair-health heatmap, earlier tabs |
| `scripts/generate_readme_figures.py` | Builds `docs/images/*.png` |

---

## Setup

```powershell
cd "c:\Users\adith\Bohr Set Visualizer"
pip install -r requirements.txt
streamlit run app.py
```

Python 3.10+ recommended.

---

## License

No license file is included; add one if you redistribute.
