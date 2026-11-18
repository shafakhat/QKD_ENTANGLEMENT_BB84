# CHSH run card — four settings, 20 000 shots each

## Why 20 000 is the right number

| shots per setting | σ(E) | σ(S) | significance above 2 |
|---|---|---|---|
| 1 000 | 0.0316 | 0.0632 | 13 σ |
| 5 000 | 0.0141 | 0.0283 | 29 σ |
| **20 000** | **0.0071** | **0.0141** | **59 σ** |
| 50 000 | 0.0045 | 0.0089 | 93 σ |
| 100 000 | 0.0032 | 0.0063 | 131 σ |

Nothing is gained past 20 000: the simulator is ideal, so extra shots only shrink an error
bar that is already 59σ from the classical bound. Use the same shot count for all four
settings — unequal N makes the combined error bar awkward to state.

## Expected counts per run (sanity check before you analyse)

At 20 000 shots, ideal:

| file | E | 00 | 11 | 01 | 10 |
|---|---|---|---|---|---|
| `chsh_a0b0` | +0.7071 | ~8 535 | ~8 535 | ~1 465 | ~1 465 |
| `chsh_a0b1` | −0.7071 | ~1 465 | ~1 465 | ~8 535 | ~8 535 |
| `chsh_a1b0` | +0.7071 | ~8 535 | ~8 535 | ~1 465 | ~1 465 |
| `chsh_a1b1` | +0.7071 | ~8 535 | ~8 535 | ~1 465 | ~1 465 |

Four rows per run, no more. If a run returns 2 or 8 rows, or the counts sum to something
other than 20 000, stop and check the shot field.

## Settings to use in Studio for each run

```
Num Shots                    20000
Measurement Bit Flip Prob.   0
Transpilation Option         decompose
Random Seed                  pick one, record it (different per run is fine)
```

Download each result (JSON/CSV) as `chsh_a0b0`, `chsh_a0b1`, `chsh_a1b0`, `chsh_a1b1`
into `results/studio_counts/`, then:

```bash
python3 code/analyze_studio_counts.py
```

## Optional, and worth it: one noisy configuration

Set **Measurement Bit Flip Probability = 0.02** and repeat the four runs. The prediction is
E → E(1−2p)², so

| p | predicted S | still violates by |
|---|---|---|
| 0.00 | 2.8284 | 59 σ |
| 0.01 | 2.7164 | 51 σ |
| **0.02** | **2.6067** | **43 σ** |
| 0.05 | 2.2910 | 21 σ |

Two S values at two readout-error levels turn the CHSH→QBER transfer function from a
theoretical statement into a measured one, and it costs four more executes.
