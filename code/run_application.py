#!/usr/bin/env python3
"""
run_application.py -- the application-layer experiments (A1-A4).

Answers the systems question a QKD paper normally leaves open:

    how much quantum key rate does a real application need, how far can the link
    be, and what happens to the service when the channel degrades?

Everything is driven by `results/qkd_results.json`, i.e. by the measured output
of the protocol stack (sifting + Cascade + Toeplitz + NIST validation) and of
the link budget.  No rate in this file is invented.

    A1  provisioning cost vs fibre length, OTP mode and AES-rekey mode
    A2  coverage cliff: how the service degrades as the supply/demand margin falls
    A3  outage ride-through: how much key buffer is needed for a given blackout
    A4  channel degradation: cost of a QBER excursion that trips the abort threshold

Outputs -> results/application_results.json
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from application_layer import Application, QKDSupply, simulate, QBER_ABORT  # noqa: E402

REPO = os.path.dirname(HERE)
RES = os.path.join(REPO, "results")
OUT = os.path.join(RES, "application_results.json")

PULSES_PER_SECOND = 1.0e8          # 100 MHz source: a conservative QKD clock
HOUR = 3600.0

# Two reference deployments.  Both are ordinary enterprise profiles, chosen so
# that the difference between them is the *encryption policy*, not the traffic.
OTP = Application(mode="otp", n_sessions=8, message_rate_hz=200,
                  message_bytes=256)                      # 3.28 Mbit/s of data
AES = Application(mode="rekey", n_sessions=64, message_rate_hz=1000,
                  message_bytes=512, rekey_bytes=1 << 20)  # 262 Mbit/s of data


def load_results():
    with open(os.path.join(RES, "qkd_results.json")) as f:
        return json.load(f)


def supply_from_row(row, pulses=PULSES_PER_SECOND):
    """QKDSupply from one row of the E5 link-budget table."""
    return QKDSupply(rounds_per_second=pulses,
                     key_bits_per_round=row["rate_asymptotic"],
                     qber=row["qber"])


# ---------------------------------------------------------------------------
# A1  provisioning cost vs distance
# ---------------------------------------------------------------------------
def a1_provisioning(d):
    rows = []
    for r in d["E5_distance"]:
        s = supply_from_row(r)
        supply_mbps = s.rate_bps / 1e6
        # the largest OTP payload this link can carry end to end
        otp_capacity = s.rate_bps / 8.0 / 1e6            # MB/s
        rows.append({
            "length_km": r["length_km"],
            "qber": r["qber"],
            "key_rate_bits_per_pulse": r["rate_asymptotic"],
            "supply_mbps": supply_mbps,
            "otp_max_payload_MBps": otp_capacity,
            "otp_feasible_for_profile": supply_mbps >= OTP.demand_bps / 1e6,
            "aes_feasible_for_profile": supply_mbps >= AES.demand_bps / 1e6,
            "pulses_needed_for_otp_profile": (OTP.demand_bps /
                                              r["rate_asymptotic"]
                                              if r["rate_asymptotic"] > 0
                                              else None),
            "pulses_needed_for_aes_profile": (AES.demand_bps /
                                              r["rate_asymptotic"]
                                              if r["rate_asymptotic"] > 0
                                              else None),
        })
    # the headline: the distance beyond which the OTP profile cannot be served
    last_ok = max((r["length_km"] for r in rows
                   if r["otp_feasible_for_profile"]), default=None)
    return {"profile_otp": f"{OTP.data_bps/1e6:.2f} Mbit/s data, 1:1 key",
            "profile_aes": (f"{AES.data_bps/1e6:.1f} Mbit/s data, AES-256 "
                            f"rekey every {AES.rekey_bytes>>10} KiB "
                            f"({AES.demand_bps/1e3:.1f} kbit/s of key)"),
            "pulse_rate_hz": PULSES_PER_SECOND,
            "otp_demand_mbps": OTP.demand_bps / 1e6,
            "aes_demand_mbps": AES.demand_bps / 1e6,
            "otp_profile_max_distance_km": last_ok,
            "rows": rows}


# ---------------------------------------------------------------------------
# A2  coverage cliff
# ---------------------------------------------------------------------------
def a2_margin_sweep(d, ratios=(0.50, 0.70, 0.85, 0.95, 0.99, 1.00, 1.01,
                               1.05, 1.15, 1.30, 1.50, 2.00, 3.00)):
    """Hold the application fixed and scale the *supply*."""
    base = supply_from_row(d["E5_distance"][0])
    rows = []
    for ratio in ratios:
        # `ratio` is supply/demand: convert it into a round rate for this link
        rounds = base.rounds_per_second * ratio * OTP.demand_bps / base.rate_bps
        s = QKDSupply(rounds_per_second=rounds,
                      key_bits_per_round=base.key_bits_per_round)
        r = simulate(s, OTP, duration_s=600.0, dt=0.002, seed=11)
        rows.append({
            "supply_over_demand": ratio,
            "supply_mbps": s.rate_bps / 1e6,
            "coverage": r.coverage,
            "underruns_per_hour": r.underruns * (HOUR / 600.0),
            "first_underrun_s": r.first_underrun_s,
            "buffer_min_Mbit": r.buffer_min_bits / 1e6,
        })
    # the margin at which the service is fully covered in this run
    ok = [r["supply_over_demand"] for r in rows if r["coverage"] >= 0.999999]
    return {"profile": "OTP, 8 sessions x 200 msg/s x 256 B",
            "duration_s": 600.0,
            "smallest_ratio_with_full_coverage": min(ok) if ok else None,
            "rows": rows}


# ---------------------------------------------------------------------------
# A3  outage ride-through
# ---------------------------------------------------------------------------
def a3_outage(d, blackouts=(10, 60, 300, 900, 1800)):
    base = supply_from_row(d["E5_distance"][0])
    # provision at 2x demand so the buffer can actually fill
    s = QKDSupply(rounds_per_second=base.rounds_per_second *
                  2.0 * OTP.demand_bps / base.rate_bps,
                  key_bits_per_round=base.key_bits_per_round)
    rows = []
    for blackout in blackouts:
        for buffer_s in (30, 60, 300, 600, 1200, 2400):
            buf = OTP.demand_bps * buffer_s
            r = simulate(s, OTP, duration_s=blackout + 300.0, dt=0.002,
                         buffer_bits=buf, outage=(60.0, float(blackout)),
                         seed=5)
            rows.append({
                "blackout_s": blackout,
                "buffer_seconds_of_demand": buffer_s,
                "coverage": r.coverage,
                "underruns": r.underruns,
                "recovery_s": r.recovery_s,
            })
    # smallest buffer that survives each blackout
    best = {}
    for b in blackouts:
        ok = [r["buffer_seconds_of_demand"] for r in rows
              if r["blackout_s"] == b and r["underruns"] == 0]
        best[str(b)] = min(ok) if ok else None
    return {"supply_over_demand": 2.0, "rows": rows,
            "min_buffer_seconds_for_clean_ride_through": best}


# ---------------------------------------------------------------------------
# A4  channel degradation and abort
# ---------------------------------------------------------------------------
def a4_degradation(d):
    """Cost of a QBER excursion that trips the BB84 abort threshold.

    When the estimated QBER crosses 11 % the link must abort and re-initialise:
    it delivers nothing for the dead time, which the application sees as an
    outage.  We price that event in seconds of lost QKD coverage.
    """
    base = supply_from_row(d["E5_distance"][0])
    rows = []
    for n_events_per_hour in (0.1, 1.0, 6.0, 30.0):
        for dead_time in (1.0, 10.0, 60.0):
            lost = n_events_per_hour * dead_time          # s/hour
            rows.append({
                "aborts_per_hour": n_events_per_hour,
                "dead_time_s": dead_time,
                "coverage_loss_fraction": lost / HOUR,
                "equivalent_downtime_min_per_day": lost * 24 / 60.0,
            })
    # and the simulation view at one representative operating point
    s = QKDSupply(rounds_per_second=base.rounds_per_second *
                  1.5 * OTP.demand_bps / base.rate_bps,
                  key_bits_per_round=base.key_bits_per_round)
    sim_rows = []
    for dead_time in (1.0, 10.0, 60.0, 300.0):
        r = simulate(s, OTP, duration_s=1800.0, dt=0.002,
                     outage=(300.0, dead_time),
                     buffer_bits=OTP.demand_bps * 60.0, seed=3)
        sim_rows.append({"dead_time_s": dead_time,
                         "coverage": r.coverage,
                         "underruns": r.underruns,
                         "recovery_s": r.recovery_s})
    return {"qber_abort_threshold": QBER_ABORT,
            "analytic": rows,
            "simulated_single_event": sim_rows}


def main():
    d = load_results()
    out = {"meta": {
        "source": "results/qkd_results.json",
        "pulse_rate_hz": PULSES_PER_SECOND,
        "qber_abort_threshold": QBER_ABORT,
        "note": ("Key supply rates are the measured protocol outputs; the key "
                 "buffer is an inventory model and the fallback path is a "
                 "classical KDF chain, as in deployed QKD networks."),
    }}
    print("== A1  provisioning vs distance ==")
    out["A1_provisioning"] = a1_provisioning(d)
    a1 = out["A1_provisioning"]
    print(f"  OTP profile (3.28 Mbit/s, 1:1 key) supported up to "
          f"{a1['otp_profile_max_distance_km']} km at 100 MHz")
    for r in a1["rows"][::3]:
        print(f"    {r['length_km']:5.0f} km  supply {r['supply_mbps']:7.3f} Mbit/s  "
              f"OTP {'yes' if r['otp_feasible_for_profile'] else 'no ':>3s}  "
              f"AES {'yes' if r['aes_feasible_for_profile'] else 'no'}")

    print("== A2  coverage cliff ==")
    out["A2_margin_sweep"] = a2_margin_sweep(d)
    a2 = out["A2_margin_sweep"]
    print(f"  smallest supply/demand ratio with full coverage: "
          f"{a2['smallest_ratio_with_full_coverage']}")
    for r in a2["rows"]:
        print(f"    ratio {r['supply_over_demand']:4.2f}  coverage {r['coverage']:.6f}  "
              f"underruns/h {r['underruns_per_hour']:8.1f}")

    print("== A3  outage ride-through ==")
    out["A3_outage"] = a3_outage(d)
    print("  min buffer (s of demand) for clean ride-through:",
          out["A3_outage"]["min_buffer_seconds_for_clean_ride_through"])

    print("== A4  degradation / abort ==")
    out["A4_degradation"] = a4_degradation(d)
    for r in out["A4_degradation"]["simulated_single_event"]:
        print(f"    dead time {r['dead_time_s']:5.0f} s -> coverage "
              f"{r['coverage']:.6f}, recovery {r['recovery_s']}")

    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
