# QKD_ENTANGLEMENT_BB84
Experimental code and measured results for an **entanglement-based BB84 quantum key
distribution** study: a Bell-pair source is simulated, its CHSH violation is
measured, the same pairs are run through a full QKD protocol (sifting, Cascade
error correction, Toeplitz privacy amplification), and the resulting keys are
validated with the NIST SP 800-22 statistical battery. The study also measures
how the scheme degrades as a function of channel noise, intercept–resend
eavesdropping, fibre length and finite key size.

Everything here is CPU-only classical simulation. No quantum hardware is required
and no hardware results are claimed.

---

## 1. What the experiment does

`code/run_qkd.py` runs six experiments and writes every number it measures to
JSON/CSV:

| experiment | question it answers |
|---|---|
| **E1** CHSH scan | How far does the measured Bell parameter `S` fall as the channel depolarises (`p_ch` = 0 → 0.04, nine points)? |
| **E2** protocol run | With 1 048 576 rounds, how many bits survive sifting, error correction and privacy amplification, and at what QBER? |
| **E3** prepare-and-measure | Does a prepare-and-measure BB84 link show the same measured QBER as the entanglement-based one? |
| **E4** intercept–resend | At what fraction `f` of intercepted rounds does the secret-key rate collapse? |
| **E5** fibre link budget | What secret-key rate survives 0 → 300 km of fibre at the measured detector efficiency? |
| **E6** finite-key analysis | How large must the sifted key be before a positive key length can be certified at all? |

`code/run_application.py` adds the application layer (A1–A4): a key-buffer service
model that asks how much key must be stored to ride through outages at a given
link length and demand rate.

---

## 2. What you need

| software | version | why it is needed |
|---|---|---|
| **Python** | 3.10 or newer | everything runs on Python |
| **NumPy** | ≥ 1.26 | array maths, the key material itself |
| **SciPy** | ≥ 1.11 | incomplete-gamma / special functions in the statistical tests |
| **Qiskit** | ≥ 1.0 | circuit construction |
| **Qiskit Aer** | ≥ 0.14 | the local simulator that actually executes the circuits |
| **Matplotlib** | ≥ 3.8 | only for `visualize_in_studio.py` |
| **Classiq SDK** | ≥ 1.30, *optional* | only for `qmod_models.py`, `run_classiq_synthesis.py` and the cloud route (needs a free Classiq account) |

Install everything that is required:

```bash
git clone https://github.com/shafakhat/qkd-entanglement-bb84.git
cd qkd-entanglement-bb84
python3 -m venv .venv && source .venv/bin/activate     # optional but recommended
pip install -r requirements.txt
```

`requirements.txt` also installs the optional Classiq SDK. If you do not want it,
`pip install numpy scipy matplotlib qiskit qiskit-aer` is enough to run every
experiment in this repository.

---

## 3. How to run it

### Smoke test (about 20 seconds)

```bash
python3 code/run_qkd.py --quick
```

Uses small shot counts. It proves the whole pipeline works and writes
`results/qkd_results_quick.json` plus `*_quick.csv`. Its numbers are noisier than
the shipped ones by design — use it to check the installation, not to compare
values.

### Full run (about 4 minutes)

```bash
python3 code/run_qkd.py
```

This reproduces the shipped `results/qkd_results.json` (the run that produced it
took 238 s on a laptop-class CPU). Expected output ends with:

```
all QKD experiments finished in <N> s
```

and the E6 line prints the finite-key crossover, which should read
`{'0.01': 1695, '0.03': 2759, '0.05': 5372}`.

### Application layer (about 10 seconds)

```bash
python3 code/run_application.py
```

writes `results/application_results.json`.

### Optional: Classiq cloud route

The QMOD models in `qmod/` are the same circuits written in Classiq's QMOD
language. To build and execute them on the hosted simulator instead of locally:

```bash
pip install classiq
classiq authenticate            # paste your free account token
python3 code/run_classiq_synthesis.py     # synthesise, save circuits, report depth/width
python3 code/visualize_in_studio.py       # printable circuit diagrams
