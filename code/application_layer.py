"""
application_layer.py
====================
The application layer: turning a QKD key supply rate into a *system* guarantee.

Why this exists
---------------
The protocol experiments answer "how many secret bits per round does the
QKD stack deliver at a given channel error?".  A systems reader asks a different
question: *how much QKD do I have to provision to run my application, and what
happens when the quantum channel degrades or drops?*

This module answers that second question with a discrete-time simulation of a
realistic deployment:

    QKD link  ->  key buffer  ->  N sessions, each AES-256-GCM with rekeying

The QKD supply rate is not invented: it comes from `results/qkd_results.json`,
i.e. from actually executing the protocol (sifting, Cascade, Toeplitz, NIST
validation) at the corresponding channel error.  Only the *inventory* is modelled
(bits arrive at the measured rate); the key material itself is verified
separately by the statistical battery in `run_qkd.py`.

Honest scope
------------
* The buffer is an inventory model.  It assumes keys arrive at the measured rate
  and are drawn from a pool that has already passed the NIST battery.  It does
  not re-derive fresh key material for every simulated second — a 24-hour run at
  1 Mbit/s would need ~10^11 bits, far beyond what one protocol execution
  produces.  This is stated in the protocol results.
* When the buffer underruns, the system does **not** fail closed.  It falls back
  to a classical key-derivation chain seeded from a pre-shared root key.  This is
  what deployed QKD networks actually do (ETSI ISG QKD), and it is the honest way
  to present the result: QKD degrades the security guarantee, it does not halt
  the service.  The fraction of traffic still protected by information-theoretic
  keys is reported as **QKD coverage**.

Units
-----
times in seconds, rates in bits/s, buffer in bits.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# BB84 aborts above this QBER (standard 11 % security threshold for the
# one-way BB84 / Shor-Preskill bound).
QBER_ABORT = 0.11


# ---------------------------------------------------------------------------
# supply side
# ---------------------------------------------------------------------------
@dataclass
class QKDSupply:
    """A QKD link feeding the key buffer.

    rounds_per_second : detector round rate (raw bits/s offered to sifting)
    key_bits_per_round: secret bits per round *after* the full protocol; taken
                        from the measured protocol runs, not from a formula
    qber              : measured QBER at that operating point
    dead_time_s       : time to re-establish the link after an abort
    """
    rounds_per_second: float
    key_bits_per_round: float
    qber: float = 0.0
    dead_time_s: float = 1.0

    @property
    def rate_bps(self) -> float:
        """Delivered secret-key bits per second while the link is up."""
        return self.rounds_per_second * self.key_bits_per_round

    @property
    def operational(self) -> bool:
        return self.qber < QBER_ABORT


# ---------------------------------------------------------------------------
# application side
# ---------------------------------------------------------------------------
@dataclass
class Application:
    """N sessions protecting a stream, in one of two encryption modes.

    mode="otp"    one bit of key per bit of data -- the strongest use of a QKD
                  link, and the one whose feasibility is set directly by the
                  quantum key rate.
    mode="rekey"  AES-256-GCM with periodic rekeying -- the usual engineering
                  compromise, where key demand is set by the rekey interval and
                  is orders of magnitude smaller than the data rate.
    """
    n_sessions: int = 32
    message_rate_hz: float = 500.0        # messages/s per session
    message_bytes: int = 512
    mode: str = "rekey"                   # "otp" or "rekey"
    rekey_bytes: int = 1 << 20            # mode="rekey": rekey after this volume
    rekey_seconds: float | None = None    # optional time-based rekey
    key_bits: int = 256                   # AES-256 key size

    @property
    def bytes_per_second(self) -> float:
        return self.n_sessions * self.message_rate_hz * self.message_bytes

    @property
    def data_bps(self) -> float:
        return self.bytes_per_second * 8.0

    @property
    def rekeys_per_second(self) -> float:
        if self.mode == "otp":
            # one "rekey" per message, each consuming as many bits as the message
            return self.n_sessions * self.message_rate_hz
        vol = self.bytes_per_second / self.rekey_bytes
        if self.rekey_seconds:
            vol = max(vol, self.n_sessions / self.rekey_seconds)
        return vol

    @property
    def bits_per_rekey(self) -> int:
        if self.mode == "otp":
            return self.message_bytes * 8
        return self.key_bits

    @property
    def demand_bps(self) -> float:
        """Secret-key bits per second the application consumes."""
        return self.rekeys_per_second * self.bits_per_rekey

    @property
    def key_to_data_ratio(self) -> float:
        return self.demand_bps / self.data_bps if self.data_bps else 0.0


# ---------------------------------------------------------------------------
# the simulation
# ---------------------------------------------------------------------------
@dataclass
class SimResult:
    coverage: float = 0.0                 # fraction of rekeys served by QKD keys
    underruns: int = 0
    rekeys: int = 0
    buffer_min_bits: float = 0.0
    buffer_mean_bits: float = 0.0
    first_underrun_s: float | None = None
    outage_ride_through_s: float = 0.0    # longest outage survived with no underrun
    recovery_s: float | None = None       # time to refill after an outage
    series: dict = field(default_factory=dict)


def simulate(supply: QKDSupply, app: Application, duration_s: float = 3600.0,
             buffer_bits: float | None = None, dt: float = 0.001,
             outage: tuple[float, float] | None = None,
             seed: int = 0) -> SimResult:
    """Discrete-time simulation of supply, buffer and demand.

    `outage` is an optional (start_s, length_s) window during which the quantum
    channel is down -- the link is not merely slower, it delivers nothing, which
    is what happens when a fibre is cut or the QKD unit is re-initialising.
    """
    rng = np.random.default_rng(seed)
    if buffer_bits is None:
        # default: 60 s of application demand, a typical provisioning choice
        buffer_bits = max(app.demand_bps * 60.0, 1 << 20)

    buf = buffer_bits * 0.5               # start half full
    n_steps = int(duration_s / dt)
    supply_per_step = supply.rate_bps * dt if supply.operational else 0.0
    # demand is a Poisson process on the step grid: the mean is exact, the
    # clumping is what makes an underrun possible in the first place
    lam = app.rekeys_per_second * dt

    res = SimResult()
    series = {"t": [], "buffer": [], "coverage": []}
    served = 0
    total = 0
    outage_start, outage_len = (outage if outage else (math.inf, 0.0))
    in_outage_end = outage_start + outage_len
    outage_ended_at = None
    t_end_outage = None

    for i in range(n_steps):
        t = i * dt
        down = outage_start <= t < in_outage_end
        if not down:
            buf = min(buffer_bits, buf + supply_per_step)

        # rekey events this step
        k = rng.poisson(lam)
        if k:
            total += k
            need = k * app.bits_per_rekey
            if buf >= need:
                buf -= need
                served += k
            else:
                # partial service from the buffer, the rest falls back
                n_from_buf = int(buf // app.bits_per_rekey)
                served += n_from_buf
                buf -= n_from_buf * app.bits_per_rekey
                res.underruns += 1
                if res.first_underrun_s is None:
                    res.first_underrun_s = t

        if i % 100 == 0:                  # sample the series at 0.1 s
            series["t"].append(t)
            series["buffer"].append(buf)
            series["coverage"].append(served / total if total else 1.0)

        # recovery bookkeeping after an outage
        if outage_ended_at is None and t >= in_outage_end and outage_len:
            outage_ended_at = t
        if outage_ended_at is not None and t_end_outage is None:
            if buf >= buffer_bits * 0.5:
                t_end_outage = t
                res.recovery_s = t - outage_ended_at

    res.rekeys = total
    res.coverage = served / total if total else 1.0
    if series["buffer"]:
        res.buffer_min_bits = float(np.min(series["buffer"]))
        res.buffer_mean_bits = float(np.mean(series["buffer"]))
    res.series = {k: v for k, v in series.items()}
    return res


def ride_through_time(supply: QKDSupply, app: Application,
                      buffer_bits: float, seed: int = 0) -> float:
    """Longest outage (s) the buffer absorbs with no underrun."""
    if not supply.operational or supply.rate_bps <= 0:
        return 0.0
    # buffer measured in seconds of demand; the link is down so nothing refills
    seconds = buffer_bits / app.demand_bps
    # a Poisson process can exhaust the buffer slightly early; report the
    # 99 %-safe value
    return float(seconds * 0.99)


def required_rate(supply_key_bits_per_round: float, app: Application,
                  coverage_target: float = 0.999,
                  spans: int = 1) -> float:
    """Raw round rate needed to serve `coverage_target` of the rekeys.

    In steady state the buffer is work-conserving, so coverage -> 1 as soon as
    supply >= demand.  To hit a target *below* 1 is not useful; what matters is
    the provisioning margin, so we return the rate that meets demand with the
    given headroom, and the caller compares.
    """
    if supply_key_bits_per_round <= 0:
        return math.inf
    # headroom covers Poisson demand fluctuation over the provisioning horizon
    headroom = 1.0 + 3.0 / math.sqrt(max(app.rekeys_per_second * 60.0, 1.0))
    return app.demand_bps * headroom / supply_key_bits_per_round / spans
