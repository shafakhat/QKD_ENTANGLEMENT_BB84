"""
run_qkd.py
==========
All QKD experiments.  Each experiment writes the numbers it measures into
`results/qkd_results.json` and the companion CSV files.

    E1  CHSH / E91 Bell-inequality scan vs the quantum-channel error p_ch
    E2  Entanglement-based BB84: full protocol (random bases, sifting, Cascade
        error correction, Toeplitz privacy amplification, statistical battery
        on the final key)
    E3  Prepare-and-measure BB84: measured QBER vs p_ch
    E4  Intercept-resend eavesdropping: QBER and key rate vs attack fraction
    E5  Key rate vs fibre length (link budget with measured QBER)
    E6  Finite-key analysis

Reproducing on Classiq cloud
----------------------------
Every circuit below is the semantic equivalent of a QMOD model from
`qmod_models.py`.  In Classiq Studio the same numbers are obtained with

    from classiq import synthesize, execute, ExecutionPreferences
    qprog = synthesize(qmod_models.model_chsh(a, b))
    res = execute(qprog, execution_preferences=ExecutionPreferences(
                      num_shots=200_000)).result()
    counts = res[0].value.parsed_counts

(the sampler `circuits.py` is used here only so that the whole study runs
offline and without an account token).
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import circuits as C
import postprocessing as P
import nist_sp800_22 as nist

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(REPO, "results")
os.makedirs(RES, exist_ok=True)

LAM = lambda p: 1.0 - 4.0 * p / 3.0          # Bloch shrink per qubit


# ---------------------------------------------------------------------------
# E1  CHSH
# ---------------------------------------------------------------------------
def chsh_scan(p_values, shots=200_000, p_ro=0.0, seed0=1000):
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    rows = []
    for i, p in enumerate(p_values):
        sim = AerSimulator(method="density_matrix", seed_simulator=seed0 + i,
                           max_parallel_threads=2)
        Es, sig = [], []
        for name, a, b in C.CHSH_SETTINGS:
            qc = C.chsh_circuit(a, b, p)
            counts = sim.run(transpile(qc, sim), shots=shots).result().get_counts()
            E, s, _ = C.chsh_correlator(counts, shots)
            Es.append(E)
            sig.append(s)
        S = C.chsh_S(Es)
        sS = float(np.sqrt(sum(s ** 2 for s in sig)))
        lam = LAM(p)
        rows.append({
            "p_ch": p, "readout_err": p_ro, "shots_per_setting": shots,
            "E_A0B0": Es[0], "E_A0B1": Es[1], "E_A1B0": Es[2], "E_A1B1": Es[3],
            "sigma_E": sig,
            "S": S, "sigma_S": sS,
            "S_theory": 2 * np.sqrt(2) * lam ** 2,
            "violation_sigma": (S - 2.0) / sS,
            "ideal_gap_sigma": (2 * np.sqrt(2) * lam ** 2 - S) / sS,
            "lambda_from_S": float(np.sqrt(max(S, 0) / (2 * np.sqrt(2)))),
            "qber_from_S": (1 - max(S, 0) / (2 * np.sqrt(2))) / 2,
            "qber_ent_theory": (1 - lam ** 2) / 2,
        })
        print(f"  p_ch={p:.4f}  S={S:.4f} +- {sS:.4f}  theory={rows[-1]['S_theory']:.4f}"
              f"  (+{(S-2)/sS:.1f} sigma over 2)")
    return rows


# ---------------------------------------------------------------------------
# E2  entanglement-based BB84
# ---------------------------------------------------------------------------
def _ent_bb84_batch_circuit(batch, alice_basis, bob_basis, p_ch):
    from qiskit import QuantumCircuit
    nq = 2 * batch
    qc = QuantumCircuit(nq, nq)
    for k in range(batch):
        q0, q1 = 2 * k, 2 * k + 1
        qc.h(q0)
        qc.cx(q0, q1)
    C.add_channel(qc, range(nq), p_ch)
    if alice_basis == 1:
        for k in range(batch):
            qc.h(2 * k)
    if bob_basis == 1:
        for k in range(batch):
            qc.h(2 * k + 1)
    qc.measure(range(nq), range(nq))
    return qc


def sample_ent_bb84(n_rounds, p_ch=0.0, p1=0.0, p_ro=0.0, batch=12, seed=None):
    """
    Sample `n_rounds` rounds of entanglement-based BB84 with random and
    independent bases on both sides.  Returns (alice_bits, bob_bits,
    alice_bases, bob_bases).
    """
    from qiskit import transpile
    rng = np.random.default_rng(seed)
    a_bases = rng.integers(0, 2, n_rounds, dtype=np.uint8)
    b_bases = rng.integers(0, 2, n_rounds, dtype=np.uint8)
    a_bits = np.zeros(n_rounds, dtype=np.uint8)
    b_bits = np.zeros(n_rounds, dtype=np.uint8)
    sim = C.simulator(p1=p1, p_ro=p_ro, method="stabilizer", seed=seed)
    for ba in (0, 1):
        for bb in (0, 1):
            sel = np.flatnonzero((a_bases == ba) & (b_bases == bb))
            if sel.size == 0:
                continue
            qc = _ent_bb84_batch_circuit(batch, ba, bb, p_ch)
            shots = int(np.ceil(sel.size / batch))
            mem = sim.run(transpile(qc, sim), shots=shots,
                          memory=True).result().get_memory()
            arr = np.frombuffer("".join(mem).encode(), dtype=np.uint8) - ord("0")
            arr = arr.reshape(shots, 2 * batch)[:, ::-1]     # column q = qubit q
            a = arr[:, 0::2].ravel()
            b = arr[:, 1::2].ravel()
            take = sel[: a.size]
            a_bits[take] = a[: take.size]
            b_bits[take] = b[: take.size]
    return a_bits, b_bits, a_bases, b_bases


def run_ent_bb84(p_ch, n_rounds=1 << 20, p_ro=0.0, seed=7, do_battery=True):
    t0 = time.time()
    a, b, ab, bb = sample_ent_bb84(n_rounds, p_ch=p_ch, p_ro=p_ro, seed=seed)
    ka, kb, keep = P.sift(a, ab, b, bb)
    n_sift = ka.size
    errors = int(np.count_nonzero(ka != kb))
    qber = errors / n_sift
    ci = P.wilson_interval(errors, n_sift)
    out = {"p_ch": p_ch, "p_ro": p_ro, "n_rounds": n_rounds,
           "n_sift": int(n_sift), "sifting_fraction": n_sift / n_rounds,
           "qber": qber, "qber_ci99": ci, "errors": errors,
           "qber_theory": (1 - LAM(p_ch) ** 2) / 2,
           "t_sampling": round(time.time() - t0, 1)}

    # --- Cascade error correction
    t1 = time.time()
    casc = P.cascade(ka, kb, qber, n_passes=4, seed=seed)
    out["cascade"] = {"corrected": casc.corrected, "leaked_bits": casc.leaked_bits,
                      "efficiency_f": casc.efficiency, "passes": casc.passes,
                      "residual_errors": int(np.count_nonzero(ka != casc.key)),
                      "keys_verified_equal": bool(np.array_equal(ka, casc.key)),
                      "t_seconds": round(time.time() - t1, 1)}

    # --- privacy amplification
    t2 = time.time()
    # guard the reconciliation-efficiency factor: with zero measured errors the
    # ratio leaked/(n h2(e)) is undefined, and f_ec < 1 is not physically
    # meaningful, so f_ec is clipped to [1, 2]
    f_ec = float(min(max(casc.efficiency, 1.0), 2.0)) if qber > 0 else 1.0
    out["f_ec_used"] = f_ec
    fk = P.finite_key_length(n_sift, qber, eps_sec=1e-10, f_ec=f_ec)
    out["finite_key"] = fk
    n_out = int(min(fk["ell"], n_sift - 1))
    final = P.toeplitz_pa(casc.key, n_out, seed=seed + 1)
    out["final_key_bits"] = int(final.size)
    out["final_key_rate"] = final.size / n_rounds
    out["t_pa"] = round(time.time() - t2, 1)

    # --- statistics of the final key
    if do_battery and final.size >= 20000:
        t3 = time.time()
        n = final.size
        m_serial = int(max(2, np.floor(np.log2(n)) - 4))
        res = nist.run_battery(final)
        out["battery"] = {k: {"p_value": v["p_value"], "pass": v["pass"]}
                          for k, v in res.items()}
        out["battery_pass"] = sum(1 for v in out["battery"].values() if v["pass"])
        out["battery_applicable"] = sum(1 for v in out["battery"].values()
                                        if v["p_value"] is not None)
        out["battery_total"] = len(out["battery"])
        out["battery_seconds"] = round(time.time() - t3, 1)
        out["final_key_ones"] = float(final.mean())
        ent = nist.approximate_entropy(final, min(10, m_serial - 1))[1]
        out["final_key_apen"] = ent.get("ApEn")
    out["t_total"] = round(time.time() - t0, 1)
    os.makedirs(os.path.join(RES, "keys"), exist_ok=True)
    np.save(os.path.join(RES, "keys", f"qkd_key_pch{p_ch}.npy"), final)
    np.save(os.path.join(RES, "keys", f"qkd_raw_pch{p_ch}.npy"), np.stack([ka, casc.key]))
    return out


# ---------------------------------------------------------------------------
# E3  prepare-and-measure BB84
# ---------------------------------------------------------------------------
def pm_qber(p_values, n_rounds=200_000, seed=11):
    rows = []
    for i, p in enumerate(p_values):
        eZ = C.sample_error_mask(n_rounds, 0, 0, p_ch=p, seed=seed + i)
        eX = C.sample_error_mask(n_rounds, 1, 1, p_ch=p, seed=seed + 100 + i)
        rows.append({"p_ch": p,
                     "qber_Z": float(eZ.mean()), "qber_X": float(eX.mean()),
                     "qber_theory": (1 - LAM(p)) / 2,
                     "qber_ent_theory": (1 - LAM(p) ** 2) / 2,
                     "n_rounds": n_rounds})
        print(f"  p_ch={p:.4f}  QBER_Z={rows[-1]['qber_Z']:.4f} "
              f"QBER_X={rows[-1]['qber_X']:.4f}  theory={(1-LAM(p))/2:.4f}")
    return rows


# ---------------------------------------------------------------------------
# E4  intercept-resend eavesdropping
# ---------------------------------------------------------------------------
def eve_attack(f_values, p_ch=0.0, model="intercept-resend"):
    """
    Intercept-resend: Eve guesses the basis correctly with probability 1/2, in
    which case she forwards the correct state; otherwise Bob's outcome is
    random, so the QBER contributed by the attack is 1/4 per attacked round.
    """
    rows = []
    for f in f_values:
        qber = f * 0.25 + (1 - f) * (1 - LAM(p_ch) ** 2) / 2
        rows.append({"attack_fraction": f, "qber": qber,
                     "asymptotic_rate": P.asymptotic_rate_bb84(qber),
                     "rate_ent": max(0.0, 1 - 2 * P.h2(qber))})
    f_star = (0.11 - (1 - LAM(p_ch) ** 2) / 2) / 0.25     # rate -> 0 at QBER=11%
    return rows, {"f_star": float(max(0.0, min(1.0, f_star))),
                  "qber_at_f_star": 0.11}


# ---------------------------------------------------------------------------
# E5  key rate vs fibre length
# ---------------------------------------------------------------------------
def distance_model(qber0, alpha_db=0.2, mu=0.5, eta_d=0.3, p_dc=1e-6,
                   f_ec=1.16, sifting=0.5, lengths_km=None):
    """
    BB84 decoy-free link budget:
        eta    = eta_d * 10^(-alpha L / 10)
        p_det  = 1 - exp(-eta mu)             (detection probability per pulse)
        e_det  = (qber0 * eta mu + 0.5 * p_dc) / (eta mu + p_dc)
        R      = sifting * p_det * [1 - h2(e) - f_ec * h2(e)]

    The secret fraction is the Devetak-Winter expression for the
    entanglement-based protocol: the h2(e_det) term is the privacy-amplification
    cost and f_ec * h2(e) is the error-correction leakage with f_ec the
    measured/assumed reconciliation inefficiency (EC only).  An earlier version
    of this function omitted the privacy-amplification term, which overstated
    the rate by ~13 % at zero distance; it is included here so that the
    asymptotic model and the finite-key bound of Eq. (finitekey) use the same
    convention.
    The misalignment term qber0 is the value measured on the simulator.
    """
    if lengths_km is None:
        lengths_km = np.arange(0, 301, 10)
    rows = []
    for L in lengths_km:
        eta = eta_d * 10 ** (-alpha_db * L / 10)
        p_det = 1 - np.exp(-eta * mu)
        e = (qber0 * eta * mu + 0.5 * p_dc) / (eta * mu + p_dc)
        r_asym = (sifting * p_det * (1.0 - P.h2(e) - f_ec * P.h2(e))
                  if e < 0.11 else 0.0)
        rows.append({"length_km": float(L), "eta": eta, "p_detect": p_det,
                     "qber": e, "rate_asymptotic": max(0.0, r_asym)})
    return rows


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(quick=False):
    results = {"meta": {"date": time.strftime("%Y-%m-%d %H:%M:%S"),
                        "simulator": "qiskit-aer (density_matrix / stabilizer)",
                        "note": "semantic equivalents of the QMOD models"}}
    t0 = time.time()

    print("\n== E1  CHSH scan ==")
    p_vals = [0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
    if quick:
        p_vals = [0.0, 0.01, 0.02]
    results["E1_chsh"] = chsh_scan(p_vals, shots=200_000 if not quick else 20_000)

    print("\n== E1b CHSH with 1% readout error ==")
    results["E1b_chsh_readout"] = chsh_scan([0.0, 0.01, 0.02], shots=100_000,
                                            p_ro=0.01, seed0=2000)

    print("\n== E2  entanglement-based BB84 (full protocol) ==")
    results["E2_protocol"] = []
    for p in ([0.0, 0.01, 0.02] if not quick else [0.01]):
        r = run_ent_bb84(p, n_rounds=1 << (20 if not quick else 15), seed=7)
        print(f"  p_ch={p}: sifted={r['n_sift']} QBER={r['qber']:.4f} "
              f"(theory {r['qber_theory']:.4f}) "
              f"f={r['f_ec_used']:.3f}"
              + ("*" if not r.get("errors") else "") + " "
              f"final={r['final_key_bits']} bits, battery "
              f"{r.get('battery_pass','-')}/{r.get('battery_applicable','-')}")
        results["E2_protocol"].append(r)

    print("\n== E3  prepare-and-measure BB84 ==")
    results["E3_pm"] = pm_qber([0.0, 0.005, 0.01, 0.02, 0.03],
                              n_rounds=200_000 if not quick else 20_000)

    print("\n== E4  intercept-resend attack ==")
    fs = np.round(np.arange(0, 1.001, 0.05), 3).tolist()
    # the attack is superimposed on the measured device error, so that f* is the
    # interception fraction a real link would tolerate, not the value for an
    # error-free device
    _p_dev = results["E2_protocol"][1]["p_ch"] if len(results["E2_protocol"]) > 1 else 0.0
    results["E4_eve"], results["E4_threshold"] = eve_attack(fs, p_ch=_p_dev)
    print("  threshold f*:", results["E4_threshold"])

    print("\n== E5  rate vs distance ==")
    results["E5_distance"] = distance_model(
        qber0=results["E2_protocol"][0]["qber"] if results["E2_protocol"] else 0.0)
    results["E5_distance_noisy"] = distance_model(
        qber0=results["E2_protocol"][1]["qber"] if len(results["E2_protocol"]) > 1 else 0.01)
    print("  0 km rate:", results["E5_distance"][0]["rate_asymptotic"])

    print("\n== E6  finite-key analysis ==")
    fk = []
    for ns in [1e3, 1e4, 1e5, 1e6, 1e7, 1e8]:
        for q in [0.01, 0.03, 0.05]:
            d = P.finite_key_length(int(ns), q, eps_sec=1e-10, f_ec=1.16)
            fk.append({"n_sift": int(ns), "qber": q, **{k: d[k] for k in
                       ("e_bit", "e_phase", "leak_ec", "ell", "rate")}})
    results["E6_finite_key"] = fk

    # Crossover lengths: the smallest sifted-key length at which the
    # finite-key bound of Eq. (finitekey) is positive, and the asymptotic rate
    # it approaches.  Both are needed to describe the finite-size penalty
    # honestly (the bound is vacuous below the crossover).
    def _nmin(e, f_ec=1.16, eps=1e-10):
        lo, hi = 1.0, 1e8
        if P.finite_key_length(int(hi), e, eps_sec=eps, f_ec=f_ec)["ell"] <= 0:
            return None
        for _ in range(60):
            mid = math.sqrt(lo * hi)
            if P.finite_key_length(int(mid), e, eps_sec=eps, f_ec=f_ec)["ell"] > 0:
                hi = mid
            else:
                lo = mid
        return int(round(hi))

    results["E6_nmin"] = {f"{e:.2f}": _nmin(e) for e in (0.01, 0.03, 0.05)}
    results["E6_asymptotic"] = {f"{e:.2f}": 1.0 - P.h2(e) - 1.16 * P.h2(e)
                                for e in (0.01, 0.03, 0.05)}
    print("  crossover n:", results["E6_nmin"],
          " asymptotic:", {k: round(v, 4) for k, v in results["E6_asymptotic"].items()})
    results["meta"]["wall_seconds"] = round(time.time() - t0, 1)

    out_json = "qkd_results_quick.json" if quick else "qkd_results.json"
    with open(os.path.join(RES, out_json), "w") as f:
        json.dump(results, f, indent=2, default=float)

    # companion csvs
    import csv
    suffix = "_quick" if quick else ""
    with open(os.path.join(RES, f"qkd_chsh_scan{suffix}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results["E1_chsh"][0].keys()))
        w.writeheader()
        for r in results["E1_chsh"]:
            w.writerow({k: (v if not isinstance(v, list) else str(v))
                        for k, v in r.items()})
    with open(os.path.join(RES, f"qkd_distance{suffix}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results["E5_distance"][0].keys()))
        w.writeheader()
        w.writerows(results["E5_distance"])
    with open(os.path.join(RES, f"qkd_finite_key{suffix}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results["E6_finite_key"][0].keys()))
        w.writeheader()
        w.writerows(results["E6_finite_key"])
    print(f"\nall QKD experiments finished in {results['meta']['wall_seconds']} s")
    return results


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
