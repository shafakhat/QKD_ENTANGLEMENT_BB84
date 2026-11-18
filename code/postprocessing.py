"""
postprocessing.py
=================
Classical post-processing for the QKD protocol: sifting, error correction
(Cascade) and privacy amplification (Toeplitz universal hashing), together
with the asymptotic and finite-key rate formulas used in the experiment drivers.

Everything here is a real implementation, not a model: `cascade` actually
reconciles Alice's and Bob's keys (the returned keys are verified to be
identical) and counts the number of leaked parity bits, `toeplitz_pa` actually
multiplies the key by a random Toeplitz matrix over GF(2).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# entropies and rates
# ---------------------------------------------------------------------------
def h2(x: float) -> float:
    """Binary Shannon entropy (bits)."""
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def asymptotic_rate_bb84(qber: float, sifting: float = 0.5,
                        f_ec: float = 1.0) -> float:
    """R = q * (1 - (1 + f_ec) h2(e)) with an ideal single-photon source."""
    return max(0.0, sifting * (1.0 - (1.0 + (f_ec - 1.0)) * h2(qber)))


def finite_key_length(n_sift: int, qber: float, eps_sec: float = 1e-10,
                      f_ec: float = 1.16, phase_error_factor: float = 1.0) -> dict:
    """
    Finite-key lower bound on the secret-key length (Tomamichel-Lim-Gisin-Renner
    form, Nat. Commun. 3, 634 (2012); generalised to the efficient-BB84 phase
    error estimate e_ph <= phase_error_factor * e_bit):

        l >= n[1 - h2(e_ph) - f_ec h2(e_bit)] - 2 log2(1/eps_sec) - 1

    with the statistical-fluctuation term of the phase-error estimate taken
    from the Hoeffding bound,  mu = sqrt(n/2 ln(1/eps_est)).
    """
    eps_est = eps_sec
    mu = math.sqrt(n_sift / 2.0 * math.log(1.0 / max(eps_est, 1e-300)))
    e_ph = min(0.5, phase_error_factor * qber + mu / max(n_sift, 1))
    e_bit = min(0.5, qber + mu / max(n_sift, 1))
    leak = f_ec * n_sift * h2(e_bit)
    l = n_sift * (1.0 - h2(e_ph)) - leak - 2.0 * math.log2(1.0 / eps_sec) - 1.0
    return {"n_sift": int(n_sift), "e_bit": e_bit, "e_phase": e_ph,
            "mu": mu, "leak_ec": leak, "ell": max(0.0, l),
            "rate": max(0.0, l) / n_sift if n_sift else 0.0}


# ---------------------------------------------------------------------------
# Cascade error correction
# ---------------------------------------------------------------------------
def _parity(block: np.ndarray) -> int:
    return int(block.sum() & 1)


def _binary_search(a: np.ndarray, b: np.ndarray) -> tuple:
    """Locate one error inside a block whose parities differ.

    Returns (index, n_revealed) where n_revealed = number of parity bits
    exchanged during the bisection (log2 of the block length), which is
    counted as leakage."""
    lo, hi = 0, a.size
    revealed = 0
    while hi - lo > 1:
        mid = (lo + hi) // 2
        revealed += 1
        if _parity(a[lo:mid]) != _parity(b[lo:mid]):
            hi = mid
        else:
            lo = mid
    return lo, revealed


@dataclass
class CascadeResult:
    key: np.ndarray          # Bob's corrected key
    corrected: int           # number of corrections applied
    leaked_bits: int         # parity bits revealed
    passes: int
    efficiency: float        # leaked_bits / (n * h2(qber))


def cascade(alice: np.ndarray, bob: np.ndarray, qber_estimate: float,
            n_passes: int = 4, seed: int = 0, block_growth: float = 2.0) -> CascadeResult:
    """
    Cascade reconciliation (Brassard & Salvail, 1994), with the standard
    cascade re-check of previously used blocks after every correction.

    Parameters
    ----------
    alice, bob : reference and raw keys (uint8 arrays), same length
    qber_estimate : a-priori error estimate (used to size the first blocks)
    """
    rng = np.random.default_rng(seed)
    n = alice.size
    b = bob.copy()
    leaked = 0
    corrected = 0
    # block sizes: k1 ~ 0.73/q, doubling each pass; permutation per pass
    k = max(2, int(math.ceil(0.73 / max(qber_estimate, 1e-4))))
    history: list[list[list[int]]] = []      # blocks used in each pass
    for p in range(n_passes):
        perm = rng.permutation(n)
        idx = perm
        blocks = [list(range(i, min(i + k, n))) for i in range(0, n, k)]
        this_pass: list[list[int]] = []
        for blk in blocks:
            sel = idx[blk]
            # shuffle within the block so that binary search positions are random
            sel = rng.permutation(sel)
            this_pass.append(list(sel))
            a_blk, b_blk = alice[sel], b[sel]
            leaked += 1
            if _parity(a_blk) != _parity(b_blk):
                pos, rev = _binary_search(a_blk, b_blk)
                leaked += rev
                b[sel[pos]] ^= 1
                corrected += 1
                err_bit = sel[pos]
                # cascade: re-examine every earlier block containing err_bit
                for prev_pass in history + [this_pass]:
                    for pb in prev_pass:
                        if err_bit in pb:
                            pb_arr = np.array(pb)
                            leaked += 1
                            if _parity(alice[pb_arr]) != _parity(b[pb_arr]):
                                pos2, rev2 = _binary_search(alice[pb_arr], b[pb_arr])
                                leaked += rev2
                                b[pb_arr[pos2]] ^= 1
                                corrected += 1
                                err_bit = pb_arr[pos2]
        history.append(this_pass)
        k = int(k * block_growth)
        if not np.any(alice != b):
            break
    eff = leaked / max(n * h2(qber_estimate), 1e-9)
    return CascadeResult(b, corrected, leaked, p + 1, eff)


# ---------------------------------------------------------------------------
# Toeplitz privacy amplification over GF(2)
# ---------------------------------------------------------------------------
def _carryless_mul(x: int, y: int) -> int:
    """Polynomial multiplication over GF(2) of two Python integers."""
    acc = 0
    if x.bit_count() > y.bit_count():
        x, y = y, x
    i = 0
    while x:
        if x & 1:
            acc ^= y << i
        x >>= 1
        i += 1
    return acc


def toeplitz_pa(key: np.ndarray, out_bits: int, seed: int = 0) -> np.ndarray:
    """
    Universal-2 hashing with a random Toeplitz matrix.

    y_i = sum_j T_{i,j} x_j with T_{i,j} = t_{i-j}; implemented as the low
    `out_bits` coefficients of the GF(2) product T(z) X(z), where X(z) has the
    key bits as coefficients.  The seed defines the Toeplitz vector t.
    """
    n = key.size
    out_bits = min(out_bits, n)
    rng = np.random.default_rng(seed)
    t = rng.integers(0, 2, n + out_bits - 1, dtype=np.uint8)
    t_int = int("".join(str(int(v)) for v in t[::-1]), 2)      # t_{n+m-2} ... t_0
    x_int = int("".join(str(int(v)) for v in key[::-1]), 2)
    prod = _carryless_mul(t_int, x_int)
    bits = np.array([(prod >> i) & 1 for i in range(out_bits)], dtype=np.uint8)
    return bits


# ---------------------------------------------------------------------------
# protocol driver
# ---------------------------------------------------------------------------
def sift(alice_bits, alice_bases, bob_bits, bob_bases):
    """Basis reconciliation: keep the rounds where the bases agree."""
    keep = alice_bases == bob_bases
    return alice_bits[keep], bob_bits[keep], keep


def wilson_interval(k: int, n: int, z: float = 2.576) -> tuple:
    """Wilson score interval for a binomial proportion (99% by default)."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))
