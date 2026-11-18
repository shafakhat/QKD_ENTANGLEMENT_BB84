// bb84_round_z0  --  QKD Classiq model, OpenQASM 3.0
OPENQASM 3.0;
include "stdgates.inc";
qubit[1] q;
bit[1] c;
c[0] = measure q[0];
