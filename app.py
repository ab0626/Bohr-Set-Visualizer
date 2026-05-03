"""
Bohr Set Visualizer: Streamlit UI for Λ_{Θ,ε} ⊂ Z_N, regularity, and random walks.
Run: streamlit run app.py
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from bohr_set import (
    NormType,
    balanced_function,
    bohr_chain_sizes,
    bohr_mask,
    bohr_set,
    bohr_size,
    cyclic_group_coloring,
    dual_group_fourier_magnitude,
    l1_spread_analysis,
    pair_health_proxy_matrix,
    random_walk_empirical,
    regularity_check,
    sample_sparse_subset_of_mask,
    shift_invariance_tv,
    subsample_indices,
    sub_bohr_sifting_report,
    uniform_on_bohr_tv,
)

st.set_page_config(page_title="Bohr Set Visualizer", layout="wide")
st.title("Bohr set Λ_{Θ,ε} in ℤ_N")

col_a, col_b = st.columns(2)
with col_a:
    N = st.number_input("N (group order)", min_value=2, max_value=5000, value=180, step=1)
with col_b:
    norm: NormType = st.selectbox("Ball in T^d", ("linf", "l1"), format_func=lambda x: "max (ℓ∞)" if x == "linf" else "mean (ℓ1)")

st.subheader("Frequency set Θ")
c1, c2, c3 = st.columns(3)
with c1:
    dim_d = st.slider("Dimension d = |Θ|", min_value=1, max_value=8, value=2)
with c2:
    theta_mode = st.radio("Θ generation", ("Arithmetic progression", "Random", "Manual top modes"))
with c3:
    seed = st.number_input("RNG seed (random Θ)", min_value=0, value=42, step=1)

rng = np.random.default_rng(int(seed))
if theta_mode == "Arithmetic progression":
    start = st.number_input("First θ", min_value=1, max_value=N - 1, value=min(7, N - 1), step=1)
    step = st.number_input("Step", min_value=1, max_value=max(1, N // 2), value=min(11, N - 1), step=1)
    thetas = [(start + j * step) % N for j in range(dim_d)]
elif theta_mode == "Random":
    thetas = list(rng.choice(np.arange(1, N), size=dim_d, replace=False))
else:
    # Low-frequency modes: 1,2,...,d — good for large symmetric Bohr sets
    thetas = [((j + 1) * max(1, N // (dim_d + 2))) % N for j in range(dim_d)]
    thetas = [t if t != 0 else 1 for t in thetas]

st.caption(f"Θ = {thetas}")

with st.expander(
    "Nested Bohr **B₂ ⊂ B₁** + sparse **A ⊂ B₁** (tabs: sub-Bohr sifting · chain · ℓ₁-spread)",
    expanded=False,
):
    st.markdown(
        "Same Θ throughout; **ε_inner < ε_outer** so **ν(B₂) ≤ ν(B₁)**. "
        "Sparse **A** is Bernoulli subsampling inside **B₁** (relative density in **B₁**)."
    )
    eps_outer = st.slider("ε_outer = ν(B₁)", min_value=0.02, max_value=0.49, value=0.12, step=0.002, format="%.3f", key="eo")
    _ei_max = min(float(eps_outer) - 0.004, 0.485)
    eps_inner = st.slider(
        "ε_inner = ν(B₂)",
        min_value=0.002,
        max_value=max(0.01, _ei_max),
        value=min(0.035, max(0.01, _ei_max)),
        step=0.002,
        format="%.3f",
        key="ei",
    )
    st.divider()
    rho_a = st.slider(
        "Bernoulli density |A|/|B₁| (each point of B₁ independently)",
        min_value=0.01,
        max_value=1.0,
        value=0.12,
        step=0.01,
        key="rho_a",
    )
    sift_seed = st.number_input("RNG seed for A", min_value=0, value=12345, step=1, key="sift_seed")

m_b1 = bohr_mask(int(N), thetas, float(eps_outer), norm=norm)
m_b2 = bohr_mask(int(N), thetas, float(eps_inner), norm=norm)
mask_a = sample_sparse_subset_of_mask(m_b1, float(rho_a), np.random.default_rng(int(sift_seed)))

eps = st.slider("Radius ε (in ℝ/ℤ)", min_value=0.002, max_value=0.49, value=0.08, step=0.002, format="%.3f")

mask = bohr_mask(int(N), thetas, float(eps), norm=norm)
lam = np.nonzero(mask)[0]
size = int(lam.size)

st.metric("|Λ_{Θ,ε}|", f"{size}  ({100.0 * size / N:.2f}% of ℤ_N)")

tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs(
    [
        "Size vs ε · membership",
        "Regularity (Def. 6.3 style)",
        "Random walk · shift invariance",
        "Sub-Bohr sifting (relative density)",
        "Increment chain (Def. 6.6 style)",
        "ℓ₁-spread heatmap (Def. 6.13)",
        "Fourier spectrum (dual group)",
        "Pair health (Sec. 7.1 proxy)",
    ]
)

with tab1:
    cleft, cright = st.columns((1, 1))
    with cleft:
        n_eps = st.slider("Curve resolution", 30, 200, 80)
    eps_grid = np.linspace(0.01, 0.49, n_eps)
    sizes = [bohr_size(int(N), thetas, float(e), norm=norm) for e in eps_grid]
    fig1 = go.Figure()
    fig1.add_trace(go.Scatter(x=eps_grid, y=sizes, mode="lines", name="|Λ|"))
    fig1.add_vline(x=eps, line_dash="dash", annotation_text="current ε")
    fig1.update_layout(
        title="|Λ_{Θ,ε}| vs ε",
        xaxis_title="ε",
        yaxis_title="Cardinality",
        height=420,
        margin=dict(l=40, r=20, t=50, b=40),
    )
    st.plotly_chart(fig1, use_container_width=True)

    with cright:
        bal_mode = st.radio(
            "Membership strip display",
            ("Indicator 1_Λ", "Balanced 1_Λ − α"),
            horizontal=True,
            key="mem_strip_mode",
        )
        alpha_mean = st.radio(
            "α for balanced strip",
            ("E on ℤ_N (|Λ|/N)", "E on B₁ (|Λ∩B₁|/|B₁|)"),
            horizontal=True,
            key="alpha_mode",
        )
        f_lam = mask.astype(float)
        if bal_mode.startswith("Balanced"):
            if alpha_mean.startswith("E on ℤ"):
                y_bar = balanced_function(f_lam, None)
            else:
                y_bar = balanced_function(f_lam, m_b1)
            col_bar = ["#2ecc71" if y_bar[i] >= 0 else "#e74c3c" for i in range(int(N))]
            title_bar = "Balanced 1_Λ − α (green ≥ 0, red < 0)"
        else:
            y_bar = f_lam
            col_bar = ["#2ecc71" if mask[i] else "#ecf0f1" for i in range(int(N))]
            title_bar = "Indicator 1_Λ (green = in Λ)"
        fig_bar = go.Figure()
        fig_bar.add_trace(
            go.Bar(
                x=list(range(int(N))),
                y=y_bar,
                name="1_{Λ}",
                marker_color=col_bar,
            )
        )
        fig_bar.update_layout(
            title=title_bar,
            xaxis_title="x ∈ ℤ_N",
            height=420,
            showlegend=False,
            bargap=0,
        )
        st.plotly_chart(fig_bar, use_container_width=True)

    with st.expander("Coloring strip (Sec. 8 — pedagogical)", expanded=False):
        st.markdown(
            r"""
The paper’s **coloring** corollary (Sec. 8, Corollary 1.2 style) lives on **$G \times G \times G$** with **$L \asymp \log\log\log \lvert G\rvert$** colors.
This toy colors **one copy** of $\mathbb{Z}_N$ and asks how many colors hit **$\Lambda$** and whether **$\Lambda$** is **monochromatic** under a fixed rule.
It is **not** a 3D corner finder; it illustrates why **quasipolynomial** density savings beat **doubly logarithmic** ones for tower-type color-size tradeoffs.
            """
        )
        Lcol = st.slider("Number of colors L", 2, min(64, max(3, int(N))), 6, key="L_colors")
        cmode = st.radio("Color rule", ("x mod L (deterministic)", "Uniform random"), horizontal=True, key="c_mode")
        cseed = st.number_input("Seed (random mode)", 0, 10**9, 101, key="c_seed")
        _rng_c = np.random.default_rng(int(cseed))
        mode_c = "random" if cmode.startswith("Uniform") else "mod"
        cstat = cyclic_group_coloring(mask, int(Lcol), mode_c, _rng_c)
        _pal = (
            "#636efa,#ef553b,#00cc96,#ab63fa,#ffa15a,#19d3f3,#ff6692,#b6e880,#ff97ff,#fecb52"
        ).split(",")
        bar_c = [_pal[int(cstat.colors[i]) % len(_pal)] for i in range(int(N))]
        opa = [1.0 if mask[i] else 0.35 for i in range(int(N))]
        fig_c = go.Figure(
            go.Bar(
                x=list(range(int(N))),
                y=np.ones(int(N)),
                marker=dict(color=bar_c, opacity=opa),
                showlegend=False,
            )
        )
        fig_c.update_layout(
            title="L-coloring of ℤ_N (full opacity = in Λ)",
            xaxis_title="x",
            yaxis=dict(showticklabels=False, range=[0, 1.05]),
            height=220,
            bargap=0,
        )
        st.plotly_chart(fig_c, use_container_width=True)
        c1, c2, c3 = st.columns(3)
        c1.metric("Distinct colors on Λ", cstat.num_colors_on_lambda)
        c2.metric("Λ monochromatic?", "yes" if cstat.lambda_monochromatic else "no")
        c3.metric("P(random i.i.d. coloring makes Λ mono)", f"{cstat.p_random_monochromatic:.2e}")
        st.caption(
            r"If the abelian theorem only gave doubly-log savings, color-size guarantees would degrade to **tower-type** dependencies; "
            r"quasipolynomial density decay is what makes an **$L \approx \log\log\log \lvert G\rvert$** scale feasible in the paper."
        )

    st.subheader("|Λ| vs dimension d")
    d_max = st.slider("Scan d from 1 to …", min_value=2, max_value=min(8, max(2, N // 2)), value=min(5, dim_d + 2))
    if theta_mode == "Arithmetic progression":
        base_t = [(start + j * step) % N for j in range(d_max)]
    elif theta_mode == "Random":
        base_t = list(rng.choice(np.arange(1, N), size=d_max, replace=False))
    else:
        base_t = [((j + 1) * max(1, N // (d_max + 2))) % N for j in range(d_max)]
        base_t = [t if t != 0 else 1 for t in base_t]
    ds = list(range(1, d_max + 1))
    sizes_d = [bohr_size(int(N), base_t[:d], float(eps), norm=norm) for d in ds]
    fig_d = go.Figure(go.Scatter(x=ds, y=sizes_d, mode="lines+markers"))
    fig_d.update_layout(
        title=f"|Λ| vs d at fixed ε = {eps:.3f}",
        xaxis_title="d",
        yaxis_title="|Λ|",
        height=360,
    )
    st.plotly_chart(fig_d, use_container_width=True)

with tab2:
    st.markdown(
        r"""
**Regularity check:** For small **$\sigma$**, compare **$\Lambda_{\varepsilon(1-\sigma)}$**, **$\Lambda_{\varepsilon}$**, **$\Lambda_{\varepsilon(1+\sigma)}$**.
If the relative symmetric difference between the inner and outer sandwich is small compared to **$\lvert \Lambda_{\varepsilon}\rvert$**,
then perturbing **$\varepsilon$** slightly changes membership in a predictable way (no sudden mass loss/gain).
        """
    )
    sigma = st.slider("σ (relative width)", min_value=0.001, max_value=0.2, value=0.05, step=0.001)
    try:
        rep = regularity_check(int(N), thetas, float(eps), float(sigma), norm=norm)
        m1, m2, m3 = st.columns(3)
        m1.metric("|Λ_{ε(1−σ)}|", rep.size_inner)
        m2.metric("|Λ_ε|", rep.size_eps)
        m3.metric("|Λ_{ε(1+σ)}|", rep.size_outer)
        st.write(
            f"**Relative symmetric difference** (inner vs outer): **{rep.relative_symmetric_difference:.4f}** "
            f"(smaller ⇒ more stable under ε perturbations)"
        )
        st.write(
            f"log ratios: log|Λ_{{ε(1+σ)}}| − log|Λ_ε| = **{rep.log_size_ratio_up:.4f}**, "
            f"log|Λ_ε| − log|Λ_{{ε(1−σ)}}| = **{rep.log_size_ratio_down:.4f}**"
        )
        if rep.is_plausible_regular:
            st.success("Passes loose regularity heuristic.")
        else:
            st.warning("Large boundary ratio — ε may sit near a jump in |Λ|.")
    except ValueError as e:
        st.error(str(e))

with tab3:
    st.markdown(
        r"""
**Random walk:** **$X_0 = 0$**, **$X_{t+1} = X_t + S_t \pmod{N}$** with **$S_t$** uniform on **$\Lambda_{\Theta,\varepsilon}$**.
Long-run visit frequencies test mixing and **shift invariance**: compare empirical **$P$** to **$P$** shifted by **$a \in \Lambda$**.
        """
    )
    steps = st.number_input("Walk length (steps)", min_value=100, max_value=5_000_000, value=80_000, step=1000)
    wseed = st.number_input("Walk seed", min_value=0, value=7, step=1)
    shift_a = st.number_input("Shift a (mod N) for TV(P, P+a)", min_value=0, max_value=int(N) - 1, value=0, step=1)

    if st.button("Run walk"):
        try:
            tidx, positions, visits = random_walk_empirical(
                int(N), thetas, float(eps), int(steps), seed=int(wseed), norm=norm
            )
            tv_uni = uniform_on_bohr_tv(visits, thetas, float(eps), norm=norm)
            tv_shift = shift_invariance_tv(visits, thetas, float(eps), int(shift_a), norm=norm)
            st.write(f"**TV to uniform on Λ:** {tv_uni:.4f}  (→ 0 if visits mix on Λ)")
            st.write(f"**TV to shift by a={shift_a}:** {tv_shift:.4f}  (→ 0 if empirical law is shift-invariant)")

            # Also report a few random a in Λ
            lam_arr = bohr_set(int(N), thetas, float(eps), norm=norm)
            if len(lam_arr) > 1:
                aa = rng.choice(lam_arr, size=min(5, len(lam_arr)), replace=False)
                tvs = [shift_invariance_tv(visits, thetas, float(eps), int(a), norm=norm) for a in aa]
                st.caption(f"Sample shifts in Λ: {list(map(int, aa))} → TV = {[round(v, 4) for v in tvs]}")

            fig_t = go.Figure()
            subs = slice(None, None, max(1, len(tidx) // 5000))
            fig_t.add_trace(
                go.Scatter(x=tidx[subs], y=positions[subs], mode="lines", line=dict(width=0.6), name="X_t")
            )
            fig_t.update_layout(
                title="Trajectory (subsampled if long)",
                xaxis_title="t",
                yaxis_title="X_t mod N",
                height=400,
            )
            st.plotly_chart(fig_t, use_container_width=True)

            top = np.argsort(visits)[-min(30, N) :][::-1]
            fig_v = go.Figure(
                go.Bar(
                    x=[str(int(x)) for x in top],
                    y=visits[top],
                    name="visits",
                )
            )
            fig_v.update_layout(title="Top visited residues", height=360)
            st.plotly_chart(fig_v, use_container_width=True)
        except ValueError as e:
            st.error(str(e))

with tab4:
    st.markdown(
        r"""
**Sub-Bohr (relative) sifting.** Majorant **$B_1$** is sparse in **$G$**, but one searches for structured subsets **$A \subset B_1$**
and refines by a smaller Bohr neighborhood **$B_2$**. Compare global density **$\lvert A\rvert/\lvert G\rvert$** to **$\lvert A \cap B_2\rvert/\lvert B_2\rvert$** (“density inside the inner Bohr tube”).
        """
    )
    if int(m_b1.sum()) == 0:
        st.warning("B₁ is empty — increase ε_outer.")
    elif np.any(m_b2 & ~m_b1):
        st.error("B₂ is not contained in B₁ — lower ε_inner below ε_outer.")
    else:
        try:
            rep = sub_bohr_sifting_report(m_b1, m_b2, mask_a)
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("|B₁|", rep.size_b1)
            c2.metric("|B₂|", rep.size_b2)
            c3.metric("|A|", rep.size_a)
            c4.metric("|A ∩ B₂|", rep.size_a_in_b2)
            st.write(
                f"**Global** |A|/|G| = **{rep.density_a_in_g:.4f}** · "
                f"**In B₁** |A|/|B₁| = **{rep.density_a_in_b1:.4f}** · "
                f"**Fraction of A landing in B₂** |A∩B₂|/|A| = **{rep.fraction_a_in_b2:.4f}** · "
                f"**Relative to B₂** |A∩B₂|/|B₂| = **{rep.density_a_cap_b2_in_b2:.4f}**"
            )
            colors = []
            for i in range(int(N)):
                if mask_a[i]:
                    colors.append("#e74c3c" if m_b2[i] else "#c0392b")
                elif m_b2[i]:
                    colors.append("#3498db")
                elif m_b1[i]:
                    colors.append("#aed6f1")
                else:
                    colors.append("#ecf0f1")
            fig_s = go.Figure(
                go.Bar(
                    x=list(range(int(N))),
                    y=np.ones(int(N)),
                    marker_color=colors,
                    showlegend=False,
                )
            )
            fig_s.update_layout(
                title="Membership (pale=B₁∖A, blue=B₂∖A, dark red=A∖B₂, bright red=A∩B₂)",
                xaxis_title="x",
                yaxis=dict(showticklabels=False, range=[0, 1.05]),
                height=240,
                bargap=0,
                margin=dict(t=40, b=40),
            )
            st.plotly_chart(fig_s, use_container_width=True)
        except ValueError as e:
            st.error(str(e))

with tab5:
    st.markdown(
        r"""
**Increment chain (Def. 6.6 style).** A **$(d,\eta)$-small sequence** uses the same frequencies, nested radii with **$\nu(B_{i+1})/\nu(B_i) \le \eta$**.
Here **$\varepsilon_i = \varepsilon_0 \eta^i$** with fixed **$\Theta$** (same rank). Compare masses **$\lvert \Lambda_{\Theta,\varepsilon_i}\rvert$** as the Bohr window zooms in — the density-increment proof repeatedly passes to smaller regular Bohr sets.
        """
    )
    eps0_chain = st.slider("Starting ε₀ (chain)", min_value=0.02, max_value=0.49, value=float(eps_outer), step=0.002, format="%.3f", key="eps0_chain")
    eta = st.slider("η (radius multiplier per step)", min_value=0.5, max_value=0.98, value=0.85, step=0.01, key="eta_chain")
    n_steps = st.slider("Chain length (steps)", min_value=1, max_value=20, value=8, key="n_chain")
    try:
        e_chain, s_chain = bohr_chain_sizes(int(N), thetas, float(eps0_chain), float(eta), int(n_steps), norm=norm)
        if len(e_chain) == 0:
            st.error("Chain produced no valid radii (check ε₀ and η).")
        else:
            st.caption(f"Radii ε_i = ε₀·η^i for i = 0…{len(e_chain)-1}")
            fig_ch = go.Figure()
            fig_ch.add_trace(
                go.Scatter(
                    x=list(range(len(s_chain))),
                    y=s_chain,
                    mode="lines+markers",
                    name="|Λ_i|",
                )
            )
            fig_ch.update_layout(
                title="|Λ_{Θ, ε_i}| along nested Bohr chain",
                xaxis_title="i",
                yaxis_title="Cardinality",
                height=420,
            )
            st.plotly_chart(fig_ch, use_container_width=True)
            tbl = {"i": list(range(len(e_chain))), "ε_i": e_chain.tolist(), "|Λ_i|": s_chain.tolist()}
            st.dataframe(tbl, use_container_width=True)
    except ValueError as e:
        st.error(str(e))

with tab6:
    st.markdown(
        r"""
**Definition 6.13** ($(B_1,B_2,\varepsilon)$ **ℓ₁-spread**): for **f : B₁ → [0,1]**,

$$
\mathbb{E}_{x \sim B_1}\,\Bigl\lvert\,\mathbb{E}_{y \sim B_2}[f(x+y)] - \mathbb{E}[f]\,\Bigr\rvert \;\le\; \varepsilon \cdot \mathbb{E}[f],
$$

with **E[f] = E_{x∼B₁}[f(x)]**. Below, **f = 1_A** (extended by **0** off **A**). The heatmap shows **pointwise**
$\bigl\lvert \mathbb{E}_{y \sim B_2}[f(x+y)] - \mathbb{E}[f]\bigr\rvert$ over **x ∈ ℤ_N**; the tracker compares the **left-hand average over x ~ B₁** to **ε · E[f]** (your ε_tolerance).
        """
    )
    eps_l1 = st.slider(
        "ε_tolerance (compare LHS to ε · E[f])",
        min_value=0.001,
        max_value=0.5,
        value=0.08,
        step=0.005,
        format="%.3f",
        key="eps_l1",
    )
    if int(m_b1.sum()) == 0 or int(m_b2.sum()) == 0:
        st.warning("Need nonempty B₁ and B₂.")
    elif np.any(m_b2 & ~m_b1):
        st.error("B₂ not ⊆ B₁ — adjust ε_inner < ε_outer.")
    else:
        f_vec = mask_a.astype(np.float64)
        rep_l1 = l1_spread_analysis(int(N), m_b1, m_b2, f_vec, float(eps_l1))
        st.metric(
            "E_{x~B1} |local(x) − E[f]|  (LHS)",
            f"{rep_l1.l1_left_side:.6f}",
            delta=f"ε·E[f] = {rep_l1.rhs_threshold:.6f}",
        )
        if rep_l1.mean_f_b1 > 0:
            if rep_l1.passes_def_if:
                st.success("LHS ≤ ε · E[f] at chosen tolerance (finite-sample analogue of Def. 6.13).")
            else:
                st.info("LHS exceeds ε · E[f] — indicator of A may not be ℓ₁-spread at this (B₁,B₂,ε).")
        fig_h = go.Figure()
        fig_h.add_trace(
            go.Bar(
                x=list(range(int(N))),
                y=rep_l1.pointwise,
                marker_color=["#9b59b6" if m_b1[i] else "#dfe6e9" for i in range(int(N))],
                name="|local−E[f]|",
            )
        )
        fig_h.update_layout(
            title="Pointwise deviation (purple highlights x ∈ B₁)",
            xaxis_title="x ∈ ℤ_N",
            yaxis_title=r"|E_y f(x+y) − E[f]|",
            height=380,
            bargap=0,
            showlegend=False,
        )
        st.plotly_chart(fig_h, use_container_width=True)
        st.caption(
            f"max over B₁: **{rep_l1.max_pointwise:.4f}** · E_B₁[f]=**{rep_l1.mean_f_b1:.4f}**"
        )

with tab7:
    st.markdown(
        r"""
**Dual characters of** $\mathbb{Z}_N$: $\chi_r(x)=e^{-2\pi i r x/N}$. For $f:\mathbb{Z}_N\to\mathbb{R}$ the transform is
$\widehat{f}(r)=\sum_x f(x)\,\chi_r(x)$. This tab plots **$\bigl\lvert\widehat{\mathbb{1}_{\Lambda}}(r)\bigr\rvert$** (and optional other indicators).
Peaks at small $r$ and at the defining modes $\Theta$ illustrate that Bohr sets are **concentrated in low / structured frequency** (cf. almost periodicity, Sec. 6.2 and Appendix A in the paper).
        """
    )
    which_f = st.selectbox(
        "Function f to transform",
        ("1_Λ (current ε)", "1_A (sparse in B₁)", "1_{B₁}", "1_{B₂}"),
        key="fourier_which",
    )
    if which_f.startswith("1_Λ"):
        f_sp = mask.astype(float)
    elif which_f.startswith("1_A"):
        f_sp = mask_a.astype(float)
    elif "B₁" in which_f or "B1" in which_f:
        f_sp = m_b1.astype(float)
    else:
        f_sp = m_b2.astype(float)
    mag = dual_group_fourier_magnitude(f_sp)
    r_axis = np.arange(int(N))
    fig_ft = go.Figure()
    fig_ft.add_trace(
        go.Scatter(x=r_axis, y=mag, mode="lines", name=r"$|\widehat{f}(r)|$", line=dict(width=1.2))
    )
    for th in thetas:
        fig_ft.add_vline(
            x=th,
            line_dash="dot",
            line_color="rgba(200,100,0,0.5)",
            annotation_text=f"θ={th}" if int(N) < 80 else None,
        )
    fig_ft.update_layout(
        title="Fourier magnitudes (orthonormal character basis)",
        xaxis_title="r (dual index in ℤ_N)",
        yaxis_title=r"|ˆf(r)|",
        height=440,
    )
    st.plotly_chart(fig_ft, use_container_width=True)
    st.caption(
        "Vertical dotted lines mark frequencies used in Θ. Large mass near low modes is expected for ‘almost periodic’ / Bohr-type functions."
    )

with tab8:
    st.markdown(
        r"""
**Pedagogical proxy (not a full Gowers grid-norm check).** Lemma 7.1’s **well-conditioned** pairs require small **Gowers grid norms**
$\lVert \mathbb{1}_X(x+\cdot)\rVert_{(B_i,\ldots)}$ simultaneously for several $i$ — too expensive to compute here.

Instead we score each $(x,y)\in B_1\times B_2$ by **three slice checks** inspired by the surrounding discussion:
**X**-slice and **Y**-slice stability under averaging along **B₂**, and **D**-stability along sums **x + y** (compare Sec. 7.1 with Lemma 7.5’s “not poor” regime).

Scores $0,1,2,3$ count how many checks pass at your thresholds $(\varepsilon_s,\varepsilon_L)$.
        """
    )
    ph_eps_s = st.slider("ε_s (slice tolerance vs δ)", min_value=0.01, max_value=0.5, value=0.12, step=0.01, key="ph_eps_s")
    ph_eps_l = st.slider("ε_L (D-tube tolerance)", min_value=0.01, max_value=0.5, value=0.15, step=0.01, key="ph_eps_l")
    ph_delta = st.radio("δ_X, δ_Y, δ_D relative to", ("ℤ_N (global mean)", "B₁ (conditional mean)"), horizontal=True, key="ph_delta")
    _opts = {
        "A (sparse)": mask_a,
        "B₁": m_b1,
        "B₂": m_b2,
        "Λ (current ε)": mask,
    }
    _names = list(_opts.keys())
    cph1, cph2, cph3 = st.columns(3)
    with cph1:
        mX = _opts[st.selectbox("Mask for X", _names, index=0, key="ph_x")]
    with cph2:
        mY = _opts[st.selectbox("Mask for Y", _names, index=1, key="ph_y")]
    with cph3:
        mD = _opts[st.selectbox("Mask for D", _names, index=2, key="ph_d")]

    max_side = st.slider("Max heatmap side (subsample B₁ / B₂ if large)", 20, 200, 72, 4, key="ph_max")

    if int(m_b1.sum()) == 0 or int(m_b2.sum()) == 0:
        st.warning("Need nonempty B₁ and B₂ (see expander).")
    elif np.any(m_b2 & ~m_b1):
        st.error("B₂ not ⊆ B₁.")
    else:
        xi0 = np.nonzero(m_b1)[0]
        yi0 = np.nonzero(m_b2)[0]
        xi = subsample_indices(xi0, int(max_side))
        yi = subsample_indices(yi0, int(max_side))
        dm = "G" if ph_delta.startswith("ℤ") else "B1"
        ph = pair_health_proxy_matrix(
            int(N),
            m_b1,
            m_b2,
            mX,
            mY,
            mD,
            float(ph_eps_s),
            float(ph_eps_l),
            delta_mode=dm,
            row_indices=xi,
            col_indices=yi,
        )
        st.metric("Fraction of pairs with score = 3", f"{ph.fraction_perfect:.3f}")
        st.metric("Mean score / 3", f"{ph.mean_normalized:.3f}")
        st.caption(
            f"δ_X={ph.delta_X:.4f}, δ_Y={ph.delta_Y:.4f}, δ_D={ph.delta_D:.4f} "
            f"(rows × cols = {len(ph.x_indices)} × {len(ph.y_indices)})"
        )
        fig_ph = go.Figure(
            data=go.Heatmap(
                z=ph.score,
                x=[str(int(y)) for y in ph.y_indices],
                y=[str(int(x)) for x in ph.x_indices],
                colorscale="Viridis",
                colorbar=dict(title="score"),
                zmin=0,
                zmax=3,
            )
        )
        fig_ph.update_layout(
            title="Pair health proxy (0–3 components satisfied)",
            xaxis_title="y ∈ B₂ (subsample)",
            yaxis_title="x ∈ B₁ (subsample)",
            height=min(520, 20 * max(8, len(ph.x_indices) // 2)),
        )
        st.plotly_chart(fig_ph, use_container_width=True)
