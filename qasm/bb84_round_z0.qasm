// bb84_round_z0  --  QKD Classiq model, OpenQASM 2.0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[1];
creg c[1];
measure q[0] -> c[0];
