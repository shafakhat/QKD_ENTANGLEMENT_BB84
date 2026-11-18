// chsh_pair  --  QKD Classiq model, OpenQASM 2.0
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
ry(-0.7853981633974483) q[0];
ry(-0.39269908169872414) q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
