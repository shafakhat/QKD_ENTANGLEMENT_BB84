// chsh_pair  --  QKD Classiq model, OpenQASM 3.0
OPENQASM 3.0;
include "stdgates.inc";
qubit[2] q;
bit[2] c;
h q[0];
cx q[0], q[1];
ry(-0.7853981633974483) q[0];
ry(-0.39269908169872414) q[1];
c[0] = measure q[0];
c[1] = measure q[1];
