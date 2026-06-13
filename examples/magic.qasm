OPENQASM 2.0;
include "qelib1.inc";

// q[0] = data qubit
// q[1] = magic state ancilla

qreg q[2];
creg c[1];
creg out[1]

// data state |+>
h q[0];

// magic state |A> = T|+>
h q[1];
t q[1];

cx q[0], q[1];

measure q[1] -> c[0];

if (c==1) s q[0];

// data state should be in T|+>
// measure in Y basis
sdg q[0];
h q[0];
measure q[0] -> out[0];

// result (measure basis: counts)
// with classical control
// Z: {"counts": { "1": 5021, "0": 4979 } }
// X: {"counts": { "0": 8546, "1": 1454 } }
// Y: {"counts": { "1": 1505, "0": 8495 } }
// without classical control
// Y: {"counts": { "1": 4970, "0": 5030 } }
// X: {"counts": { "1": 1498, "0": 8502 } }
// Z: {"counts": { "1": 4940, "0": 5060 } }