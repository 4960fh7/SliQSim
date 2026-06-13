OPENQASM 2.0;
include "qelib1.inc";

// q[0], q[1]: Controls
// q[4]: Target
// q[2], q[3]: Ancilla
qreg q[5];
creg c[2]; 
creg out[3];

h q[0];
h q[1];

h q[2]; 
cx q[0], q[3];
cx q[2], q[0];
cx q[2], q[1];
cx q[1], q[3];
tdg q[0];
tdg q[1];
t q[2];
t q[3];
cx q[1], q[3];
cx q[2], q[0];
cx q[2], q[1];
cx q[0], q[3];
h q[2];

s q[2];
cx q[2], q[4];
h q[2];

measure q[2] -> c[0];

if(c==1) cz q[0], q[1];

// sdg q[0];
// sdg q[1];
// sdg q[4];
// h q[0];
// h q[1];
// h q[4];

measure q[0] -> out[0];
measure q[1] -> out[1];
measure q[4] -> out[2];

// result
// with classical control
// circuit: ccx
// Z: {"counts": { "111": 2499, "010": 2482, "000": 2533, "001": 2486 } }
// X: {"counts": { "111": 1231, "000": 5047, "101": 1228, "100": 1252, "110": 1242 } }
// Y: {"counts": { "011": 3063, "110": 625, "100": 3186, "000": 653, "111": 617, "010": 616, "101": 628, "001": 612 } }
// without classical control
// circuit: ccx,csdg
// Z: {"counts": { "001": 2502, "010": 2465, "111": 2531, "000": 2502 } }
// X: {"counts": { "110": 629, "001": 593, "011": 590, "101": 635, "100": 3204, "010": 629, "111": 623, "000": 3097 } }
// Y: {"counts": { "110": 1271, "101": 1251, "100": 1233, "011": 2457, "000": 2517, "111": 1271 } }