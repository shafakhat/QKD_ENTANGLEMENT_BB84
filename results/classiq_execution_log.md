# Classiq Studio execution log — QKD experiments

Every job below was executed in Classiq Studio and the counts were exported
to `results/studio_counts/`. 8 jobs in this repository,
123000 shots total, 80000 shots feeding the reported results.

| program | shots | completed | purpose | used in results |
|---|---|---|---|---|
| `bb84_batch.qasm` | 1000 | 2026-11-05 19:42 | diagnostic: fixed-pattern demo circuit | no (diagnostic) |
| `bb84_round.qasm` | 20000 | 2026-11-05 20:17 | diagnostic: single-round unit cell | no (diagnostic) |
| `chsh_pair.qasm` | 20000 | 2026-11-05 20:19 | diagnostic: b=pi/8, not a CHSH setting | no (diagnostic) |
| `chsh_a0b0.qasm` | 20000 | 2026-11-05 20:26 | CHSH correlator E(a0,b0) | yes |
| `chsh_a0b1.qasm` | 20000 | 2026-11-05 20:28 | CHSH correlator E(a0,b1) | yes |
| `chsh_a1b0.qasm` | 20000 | 2026-11-05 20:46 | CHSH correlator E(a1,b0) | yes |
| `chsh_a1b1.qasm` | 20000 | 2026-11-05 20:48 | CHSH correlator E(a1,b1) | yes |
| `bb84_round.qasm` | 2000 | 2026-11-05 20:56 | diagnostic: basis-mismatch round (bit 1, Alice X, Bob Z) | no (diagnostic) |
