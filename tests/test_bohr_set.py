"""Tests for bohr_set.py (core numerics; UI not exercised here)."""

from __future__ import annotations

import numpy as np

from bohr_set import (
    balanced_function,
    bohr_chain_sizes,
    bohr_mask,
    bohr_size,
    cyclic_group_coloring,
    dual_group_fourier_magnitude,
    l1_spread_analysis,
    pair_health_proxy_matrix,
    regularity_check,
    sample_sparse_subset_of_mask,
    subsample_indices,
    torus_norm,
)


def test_torus_norm_half():
    assert torus_norm(0.0) == 0.0
    assert abs(float(torus_norm(0.25)) - 0.25) < 1e-9
    assert abs(float(torus_norm(-0.25)) - 0.25) < 1e-9


def test_bohr_mask_nested():
    N, th, e_big, e_small = 120, [5, 17], 0.15, 0.05
    m_big = bohr_mask(N, th, e_big)
    m_small = bohr_mask(N, th, e_small)
    assert np.all(~m_small | m_big)


def test_bohr_chain_monotone_sizes():
    N, th = 100, [3, 11]
    eps, eta, steps = 0.2, 0.88, 6
    e_arr, s_arr = bohr_chain_sizes(N, th, eps, eta, steps)
    assert len(e_arr) >= 2
    assert np.all(np.diff(e_arr) < 0)
    assert np.all(np.diff(s_arr) <= 0)


def test_regularity_check_finite():
    rep = regularity_check(90, [7, 19], 0.12, 0.04)
    assert rep.size_eps > 0
    assert not np.isnan(rep.log_size_ratio_up)


def test_l1_spread_analysis_basic():
    N = 48
    m1 = bohr_mask(N, [5], 0.2)
    m2 = bohr_mask(N, [5], 0.08)
    rng = np.random.default_rng(0)
    a = sample_sparse_subset_of_mask(m1, 0.5, rng)
    f = a.astype(float)
    rep = l1_spread_analysis(N, m1, m2, f, 0.5)
    assert rep.pointwise.shape == (N,)
    assert rep.mean_f_b1 >= 0


def test_fourier_dc_equals_sum():
    N = 36
    m = bohr_mask(N, [4], 0.18).astype(float)
    mag = dual_group_fourier_magnitude(m)
    assert mag.shape == (N,)
    assert abs(mag[0] - float(m.sum())) < 1e-9


def test_balanced_mean_zero_global():
    m = np.array([1, 0, 1, 0, 0], dtype=float)
    b = balanced_function(m, None)
    assert abs(float(b.mean())) < 1e-12


def test_coloring_mod_mono_probability():
    N = 24
    lam = np.zeros(N, dtype=bool)
    lam[[2, 5, 11]] = True
    rng = np.random.default_rng(42)
    st = cyclic_group_coloring(lam, 4, "mod", rng)
    assert st.colors.shape == (N,)
    k = int(lam.sum())
    assert abs(st.p_random_monochromatic - float(4 ** (1 - k))) < 1e-9


def test_pair_health_matrix_shape():
    N = 40
    m1 = bohr_mask(N, [3], 0.22)
    m2 = bohr_mask(N, [3], 0.09)
    rng = np.random.default_rng(1)
    mx = m1 & (rng.random(N) < 0.5)
    xi = subsample_indices(np.nonzero(m1)[0], 8)
    yi = subsample_indices(np.nonzero(m2)[0], 6)
    ph = pair_health_proxy_matrix(
        N, m1, m2, mx.astype(float), m1.astype(float), m2.astype(float),
        0.2, 0.2, "G", row_indices=xi, col_indices=yi,
    )
    assert ph.score.shape == (len(xi), len(yi))


def test_app_imports_without_running_streamlit():
    """Smoke import: Streamlit app module loads (may warn if streamlit missing)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("app", "app.py")
    assert spec and spec.loader
    # Do not execute full app (starts session context); parse only
    import ast
    with open("app.py", encoding="utf-8") as f:
        ast.parse(f.read())
