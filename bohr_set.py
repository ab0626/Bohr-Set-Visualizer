"""
Bohr sets in the cyclic group Z_N.

A Bohr set Λ_{Θ,ε} ⊂ Z_N (with Θ = {θ_1,…,θ_d} ⊂ Z_N) consists of all x such that
each character χ_θ(x) = exp(2πi θx/N) is "close to 1" in T = R/Z, i.e.
||(θx)/N||_T < ε in the max (ℓ_∞) sense on the d-torus, where ||t||_T is the
distance from t mod 1 to 0 in [0, 1/2].
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Literal, Sequence

NormType = Literal["linf", "l1"]
ColoringMode = Literal["mod", "random"]


def torus_norm(t: np.ndarray | float) -> np.ndarray | float:
    """Distance in R/Z to 0, in [0, 1/2]."""
    t = np.asarray(t, dtype=float)
    u = t - np.floor(t)
    return np.minimum(u, 1.0 - u)


def bohr_mask(
    N: int,
    thetas: Sequence[int],
    eps: float,
    norm: NormType = "linf",
) -> np.ndarray:
    """
    Boolean mask of length N: x is in Λ_{Θ,ε} iff mask[x] is True.
    """
    if N < 1:
        raise ValueError("N must be positive")
    if eps <= 0 or eps > 0.5:
        raise ValueError("ε should lie in (0, 1/2] for a nontrivial Bohr ball")
    thetas = list(thetas)
    if not thetas:
        return np.ones(N, dtype=bool)

    x = np.arange(N, dtype=np.int64)
    # Phase distances in T for each θ
    dists = np.stack([torus_norm(th * x / N) for th in thetas], axis=0)
    if norm == "linf":
        m = np.max(dists, axis=0)
    else:
        m = np.sum(dists, axis=0) / max(len(thetas), 1)
    return m < eps


def bohr_set(N: int, thetas: Sequence[int], eps: float, norm: NormType = "linf") -> np.ndarray:
    """Return sorted array of x in Z_N in Λ_{Θ,ε}."""
    m = bohr_mask(N, thetas, eps, norm=norm)
    return np.nonzero(m)[0].astype(np.int64)


def bohr_size(N: int, thetas: Sequence[int], eps: float, norm: NormType = "linf") -> int:
    return int(bohr_mask(N, thetas, eps, norm=norm).sum())


@dataclass
class RegularityReport:
    """
    Heuristic regularity (Def. 6.3 style): compare sizes at ε and at ε scaled by
    (1±σ). If the relative boundary |Λ_{ε(1+σ)} \\ Λ_{ε(1-σ)}| / |Λ_ε| is small,
    small radius perturbations predictably control mass (no sharp jumps).
    """

    size_eps: int
    size_inner: int
    size_outer: int
    relative_symmetric_difference: float
    log_size_ratio_up: float
    log_size_ratio_down: float

    @property
    def is_plausible_regular(self) -> bool:
        """Loose flag: symmetric-difference mass / size at ε is below 0.5."""
        return self.relative_symmetric_difference < 0.5 and self.size_eps > 0


def regularity_check(
    N: int,
    thetas: Sequence[int],
    eps: float,
    sigma: float,
    norm: NormType = "linf",
) -> RegularityReport:
    """
    For small σ, compare Λ_{Θ, ε(1-σ)}, Λ_{Θ, ε}, Λ_{Θ, ε(1+σ)}.
    """
    e0 = eps
    e_in = max(e0 * (1.0 - sigma), 1e-9)
    e_out = min(e0 * (1.0 + sigma), 0.5 - 1e-9)
    if e_in >= e0 or e_out <= e0:
        raise ValueError("σ too large for given ε; need ε(1-σ) < ε < ε(1+σ) within (0,1/2)")

    m0 = bohr_mask(N, thetas, e0, norm=norm)
    mi = bohr_mask(N, thetas, e_in, norm=norm)
    mo = bohr_mask(N, thetas, e_out, norm=norm)

    s0, si, so = int(m0.sum()), int(mi.sum()), int(mo.sum())
    # Symmetric difference between inner and outer sandwiching ε
    sym = int((mi ^ mo).sum())
    rel = sym / s0 if s0 else float("inf")

    def _lr(a: int, b: int) -> float:
        if a <= 0 or b <= 0:
            return float("nan")
        return float(np.log(b) - np.log(a))

    return RegularityReport(
        size_eps=s0,
        size_inner=si,
        size_outer=so,
        relative_symmetric_difference=rel,
        log_size_ratio_up=_lr(s0, so),
        log_size_ratio_down=_lr(si, s0),
    )


def random_walk_empirical(
    N: int,
    thetas: Sequence[int],
    eps: float,
    steps: int,
    seed: int | None = None,
    norm: NormType = "linf",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Random walk: X_0 = 0, X_{t+1} = X_t + S_t (mod N) with S_t uniform on Λ_{Θ,ε}.
    Returns (times, positions, visit_counts on Z_N).
    """
    rng = np.random.default_rng(seed)
    lam = bohr_set(N, thetas, eps, norm=norm)
    if len(lam) == 0:
        raise ValueError("Bohr set is empty; increase ε or change Θ")
    pos = 0
    positions = np.zeros(steps + 1, dtype=np.int64)
    visits = np.zeros(N, dtype=np.int64)
    visits[0] = 1
    for t in range(steps):
        s = int(rng.choice(lam))
        pos = (pos + s) % N
        positions[t + 1] = pos
        visits[pos] += 1
    tidx = np.arange(steps + 1)
    return tidx, positions, visits


def shift_invariance_tv(
    visit_counts: np.ndarray,
    thetas: Sequence[int],
    eps: float,
    a: int,
    norm: NormType = "linf",
) -> float:
    """
    Total variation distance between normalized visit law P and P shifted by a
    on Z_N: (1/2) sum_x |P(x) - P(x-a)|. Small values support approximate shift
    invariance of the empirical measure on long walks (Lemma 6.5 style).
    """
    N = len(visit_counts)
    lam_mask = bohr_mask(N, thetas, eps, norm=norm)
    p = visit_counts.astype(float)
    tot = p.sum()
    if tot <= 0:
        return float("nan")
    p = p / tot
    q = np.roll(p, a % N)
    return 0.5 * float(np.abs(p - q).sum())


def uniform_on_bohr_tv(visit_counts: np.ndarray, thetas: Sequence[int], eps: float, norm: NormType = "linf") -> float:
    """TV distance from visits to uniform on Λ_{Θ,ε}."""
    N = len(visit_counts)
    m = bohr_mask(N, thetas, eps, norm=norm)
    k = int(m.sum())
    if k == 0:
        return float("nan")
    u = np.zeros(N)
    u[m] = 1.0 / k
    p = visit_counts.astype(float)
    tot = p.sum()
    if tot <= 0:
        return float("nan")
    p = p / tot
    return 0.5 * float(np.abs(p - u).sum())


# --- Sub-Bohr sifting, Bohr chains, ℓ₁-spreadness (paper §6) ---------------------------------


@dataclass
class SubBohrSiftReport:
    """B1 = outer Bohr, B2 = inner (smaller radius), A ⊂ B1 sparse."""

    size_g: int
    size_b1: int
    size_b2: int
    size_a: int
    size_a_in_b2: int
    density_a_in_g: float
    density_a_in_b1: float
    fraction_a_in_b2: float  # |A ∩ B2| / |A|
    density_a_cap_b2_in_b2: float  # |A ∩ B2| / |B2|


def sample_sparse_subset_of_mask(
    mask_b1: np.ndarray,
    relative_density: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Random A ⊂ {x : mask_b1[x]} by Bernoulli(relative_density) on each point of B1.
    Returns boolean length-N mask for A.
    """
    m = mask_b1.astype(bool)
    if not m.any():
        return np.zeros(len(mask_b1), dtype=bool)
    p = float(np.clip(relative_density, 0.0, 1.0))
    draw = rng.random(m.sum()) < p
    out = np.zeros(len(mask_b1), dtype=bool)
    out[np.nonzero(m)[0]] = draw
    return out


def sub_bohr_sifting_report(
    mask_b1: np.ndarray,
    mask_b2: np.ndarray,
    mask_a: np.ndarray,
) -> SubBohrSiftReport:
    """Statistics for relative sifting: how much of sparse A ⊂ B1 lands in inner B2 ⊂ B1."""
    N = len(mask_b1)
    b1 = mask_b1.astype(bool)
    b2 = mask_b2.astype(bool)
    a = mask_a.astype(bool)
    if np.any(a & ~b1):
        raise ValueError("A must be contained in B1")
    if not np.all((~b2) | b1):
        raise ValueError("B2 must be contained in B1 (use smaller ε for inner Bohr)")
    size_g = N
    sb1, sb2 = int(b1.sum()), int(b2.sum())
    sa = int(a.sum())
    a_in_b2 = int((a & b2).sum())
    return SubBohrSiftReport(
        size_g=size_g,
        size_b1=sb1,
        size_b2=sb2,
        size_a=sa,
        size_a_in_b2=a_in_b2,
        density_a_in_g=sa / size_g if size_g else 0.0,
        density_a_in_b1=sa / sb1 if sb1 else 0.0,
        fraction_a_in_b2=a_in_b2 / sa if sa else 0.0,
        density_a_cap_b2_in_b2=a_in_b2 / sb2 if sb2 else 0.0,
    )


def bohr_chain_sizes(
    N: int,
    thetas: Sequence[int],
    eps0: float,
    eta: float,
    steps: int,
    norm: NormType = "linf",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Radii ε_i = ε0 * η^i (same Θ), nested Bohr sets B_i = Λ_{Θ, ε_i}.
    Returns (epsilons, sizes) arrays of length steps+1.
    """
    if eta <= 0 or eta >= 1:
        raise ValueError("η should lie in (0, 1) for shrinking radii")
    out_eps = []
    out_sz = []
    for i in range(steps + 1):
        e = eps0 * (eta**i)
        if e <= 0 or e > 0.5:
            break
        out_eps.append(e)
        out_sz.append(bohr_size(N, thetas, float(e), norm=norm))
    return np.array(out_eps, dtype=float), np.array(out_sz, dtype=np.int64)


@dataclass
class L1SpreadReport:
    """Definition 6.13 style: E_{x~B1} |E_{y~B2} f(x+y) - E[f]| vs ε·E[f]."""

    mean_f_b1: float
    l1_left_side: float  # E_{x~B1} |...|
    rhs_threshold: float  # ε_user * E[f] for comparison
    passes_def_if: bool  # l1_left_side <= rhs_threshold
    max_pointwise: float
    pointwise: np.ndarray  # length N, |E_y f(x+y) - E[f]| for each x


def l1_spread_analysis(
    N: int,
    mask_b1: np.ndarray,
    mask_b2: np.ndarray,
    f: np.ndarray,
    epsilon_ratio: float,
) -> L1SpreadReport:
    """
    f: values in [0,1] on Z_N (e.g. indicator of A). Expectations over uniform B1, B2 as in Def. 6.13.
    E[f] = E_{x~B1}[f(x)]. For each x, local(x) = E_{y~B2}[f(x+y)] with y uniform on B2 as a subset of G.
    """
    b1 = mask_b1.astype(bool)
    b2 = mask_b2.astype(bool)
    f = np.asarray(f, dtype=float).reshape(N)
    if f.shape[0] != N:
        raise ValueError("f length must equal N")
    idx2 = np.nonzero(b2)[0].astype(np.int64)
    k2 = int(idx2.size)
    if k2 == 0:
        z = np.zeros(N)
        return L1SpreadReport(0.0, float("nan"), float("nan"), False, float("nan"), z)
    sb1 = int(b1.sum())
    if sb1 == 0:
        z = np.zeros(N)
        return L1SpreadReport(0.0, float("nan"), float("nan"), False, float("nan"), z)

    mean_b1 = float((f * b1).sum() / sb1)
    acc = np.zeros(N, dtype=float)
    for y in idx2:
        acc += np.roll(f, -int(y))
    local = acc / k2
    pointwise = np.abs(local - mean_b1)
    l1_mean = float((pointwise * b1).sum() / sb1)
    rhs = float(epsilon_ratio) * mean_b1 if mean_b1 > 0 else 0.0
    passes = l1_mean <= rhs if mean_b1 > 0 else l1_mean == 0
    return L1SpreadReport(
        mean_f_b1=mean_b1,
        l1_left_side=l1_mean,
        rhs_threshold=rhs,
        passes_def_if=passes,
        max_pointwise=float(np.max(pointwise[b1]) if b1.any() else 0.0),
        pointwise=pointwise,
    )


# --- Fourier spectrum, balanced functions, pair-health proxies (presentation layer) ----------


def dual_group_fourier_magnitude(f: np.ndarray) -> np.ndarray:
    """
    For f : Z_N → ℂ (real-valued input), return |\\hat f(r)| with respect to the orthonormal
    characters χ_r(x) = e^{-2π i r x / N} (numpy.fft convention).
    """
    v = np.asarray(f, dtype=np.complex128).reshape(-1)
    return np.abs(np.fft.fft(v))


def balanced_function(f: np.ndarray, domain_mask: np.ndarray | None = None) -> np.ndarray:
    """
    Return f − α with α = E[f] on ℤ_N, or α = E[f | domain] if domain_mask is given.
    """
    f = np.asarray(f, dtype=float).reshape(-1)
    N = f.size
    if domain_mask is None:
        alpha = float(f.mean()) if N else 0.0
    else:
        m = domain_mask.astype(bool).reshape(-1)
        w = int(m.sum())
        alpha = float((f * m).sum() / w) if w else 0.0
    return f - alpha


def bohr_slice_expect_vector(f: np.ndarray, b2_idx: np.ndarray, N: int) -> np.ndarray:
    """For each t ∈ Z_N, return E_{y~B2}[f(t+y)] with y uniform on the set B2 (indices in b2_idx)."""
    f = np.asarray(f, dtype=float).reshape(N)
    idx = np.asarray(b2_idx, dtype=np.int64)
    if idx.size == 0:
        return np.zeros(N, dtype=float)
    acc = np.zeros(N, dtype=float)
    for y in idx:
        acc += np.roll(f, -int(y))
    return acc / float(idx.size)


@dataclass
class PairHealthMatrix:
    """Finite-dimensional proxy for Sec. 7.1 style slice stability + D stability (not full grid norms)."""

    x_indices: np.ndarray  # rows: points in B1
    y_indices: np.ndarray  # cols: points in B2
    score: np.ndarray  # (len_x, len_y) in {0,1,2,3}
    fraction_perfect: float  # P[all 3 checks pass]
    mean_normalized: float  # E[score] / 3
    delta_X: float
    delta_Y: float
    delta_D: float


@dataclass
class BohrColoringStats:
    """Pedagogical: L-coloring of ℤ_N, restrict to Bohr mask Λ."""

    colors: np.ndarray  # length N, values in 0 .. L-1
    num_colors_on_lambda: int
    lambda_monochromatic: bool
    p_random_monochromatic: float  # if each vertex i.i.d. uniform L, P(Λ mono | |Λ|>0)


def cyclic_group_coloring(
    mask_lambda: np.ndarray,
    n_colors: int,
    mode: ColoringMode,
    rng: np.random.Generator,
) -> BohrColoringStats:
    """
    Assign each x ∈ ℤ_N a color in {0,…,L-1}: either c(x)=x mod L or uniform random.
    Report how many distinct colors appear on Λ and whether Λ is monochromatic.
    For random mode, P(all |Λ| points same color) = L^{1-|Λ|} (empty Λ → 0).
    """
    N = int(len(mask_lambda))
    L = int(n_colors)
    if L < 2:
        raise ValueError("need at least 2 colors")
    if mode == "mod":
        cols = (np.arange(N, dtype=np.int64) % L).astype(np.int64)
    else:
        cols = rng.integers(0, L, size=N, endpoint=False).astype(np.int64)
    lam = np.nonzero(mask_lambda)[0]
    if lam.size == 0:
        return BohrColoringStats(cols, 0, False, 0.0)
    sub = cols[lam]
    uniq = int(len(np.unique(sub)))
    mono = uniq == 1
    p_mono = float(L ** (1 - lam.size)) if lam.size > 0 else 0.0
    return BohrColoringStats(cols, uniq, mono, p_mono)


def subsample_indices(idx: np.ndarray, max_n: int) -> np.ndarray:
    """Deterministic stride subsample for large Bohr sets (heatmap caps)."""
    idx = np.asarray(idx, dtype=np.int64).reshape(-1)
    if idx.size <= max_n:
        return idx
    step = max(1, idx.size // max_n)
    return idx[::step][:max_n]


def pair_health_proxy_matrix(
    N: int,
    mask_b1: np.ndarray,
    mask_b2: np.ndarray,
    mask_x: np.ndarray,
    mask_y: np.ndarray,
    mask_d: np.ndarray,
    eps_s: float,
    eps_l: float,
    delta_mode: Literal["G", "B1"] = "G",
    row_indices: np.ndarray | None = None,
    col_indices: np.ndarray | None = None,
) -> PairHealthMatrix:
    """
    Pedagogical proxy only: Lemma 7.1 uses Gowers **grid norms** on Bohr slices — far too heavy to compute here.

    Instead:
      • **Slice stability**: |E_{z~B2} 1_X(x+z) − δ_X|, |E_{z~B2} 1_Y(y+z) − δ_Y| small (ε_s-scale).
      • **Diagonal stability**: |E_{z~B2} 1_D(x+y+z) − δ_D| small (ε_L-scale).

    Each contributes one point toward a 3-point score per pair (x,y) ∈ B1×B2.
    """
    N = int(N)
    b1 = mask_b1.astype(bool)
    b2 = mask_b2.astype(bool)
    ix = np.asarray(mask_x, dtype=float).reshape(N)
    iy = np.asarray(mask_y, dtype=float).reshape(N)
    id_ = np.asarray(mask_d, dtype=float).reshape(N)
    idx2 = np.nonzero(b2)[0].astype(np.int64)

    if delta_mode == "G":
        delta_X = float(ix.mean())
        delta_Y = float(iy.mean())
        delta_D = float(id_.mean())
    else:
        sb1 = max(int(b1.sum()), 1)
        delta_X = float((ix * b1).sum() / sb1)
        delta_Y = float((iy * b1).sum() / sb1)
        delta_D = float((id_ * b1).sum() / sb1)

    sx = bohr_slice_expect_vector(ix, idx2, N)
    sy = bohr_slice_expect_vector(iy, idx2, N)
    sd = bohr_slice_expect_vector(id_, idx2, N)

    thr_x = max(2.0 * eps_s * max(delta_X, 1e-12), 1e-9)
    thr_y = max(2.0 * eps_s * max(delta_Y, 1e-12), 1e-9)
    thr_d = max(2.0 * eps_l * max(delta_D, 1e-12), 1e-9)

    xi = row_indices if row_indices is not None else np.nonzero(b1)[0].astype(np.int64)
    yi = col_indices if col_indices is not None else np.nonzero(b2)[0].astype(np.int64)
    xi = np.asarray(xi, dtype=np.int64).reshape(-1)
    yi = np.asarray(yi, dtype=np.int64).reshape(-1)
    nx, ny = int(xi.size), int(yi.size)
    score = np.zeros((nx, ny), dtype=np.float64)
    for a in range(nx):
        x = int(xi[a])
        ok_x = abs(sx[x] - delta_X) <= thr_x
        for b in range(ny):
            y = int(yi[b])
            ok_y = abs(sy[y] - delta_Y) <= thr_y
            t = (x + y) % N
            ok_d = abs(sd[t] - delta_D) <= thr_d
            score[a, b] = float(ok_x + ok_y + ok_d)

    total = float(score.size) if score.size else 1.0
    fraction_perfect = float((score >= 2.99).sum() / total)
    mean_normalized = float(0.0 if total == 0 else score.sum() / (3.0 * total))

    return PairHealthMatrix(
        x_indices=xi,
        y_indices=yi,
        score=score,
        fraction_perfect=fraction_perfect,
        mean_normalized=mean_normalized,
        delta_X=delta_X,
        delta_Y=delta_Y,
        delta_D=delta_D,
    )
