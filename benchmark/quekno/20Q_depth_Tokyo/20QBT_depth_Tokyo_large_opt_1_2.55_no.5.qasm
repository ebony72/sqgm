OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
cx q[17],q[5];
cx q[0],q[1];
h q[18];
h q[5];
cx q[15],q[6];
h q[2];
h q[8];
h q[6];
h q[8];
cx q[15],q[6];
h q[6];
h q[6];
h q[15];
cx q[13],q[8];
h q[2];
h q[7];
h q[8];
h q[5];
h q[7];
h q[14];
h q[1];
h q[15];
h q[18];
cx q[15],q[3];
cx q[6],q[18];
cx q[2],q[13];
h q[5];
h q[14];
h q[6];
h q[5];
h q[1];
h q[1];
h q[0];
h q[8];
cx q[14],q[7];
h q[2];
cx q[18],q[10];
h q[14];
h q[1];
cx q[6],q[18];
cx q[14],q[7];
h q[2];
cx q[2],q[13];
h q[14];
h q[17];
h q[5];
h q[0];
cx q[14],q[7];
h q[13];
cx q[17],q[5];
h q[15];
h q[2];
cx q[0],q[1];
cx q[14],q[7];
cx q[0],q[1];
h q[5];
cx q[0],q[1];
cx q[15],q[3];
h q[10];
h q[8];
h q[17];
h q[8];
h q[8];
h q[0];
h q[8];
cx q[8],q[14];
cx q[0],q[13];
h q[12];
h q[14];
h q[3];
h q[16];
h q[8];
h q[6];
h q[6];
h q[0];
cx q[0],q[6];
h q[3];
cx q[13],q[19];
h q[5];
cx q[1],q[14];
h q[19];
h q[1];
h q[10];
cx q[6],q[19];
h q[7];
h q[14];
h q[3];
h q[13];
h q[14];
h q[16];
h q[6];
h q[12];
cx q[16],q[5];
cx q[19],q[3];
h q[6];
cx q[9],q[17];
h q[7];
cx q[1],q[10];
cx q[4],q[7];
h q[13];
h q[17];
h q[3];
h q[6];
h q[7];
cx q[13],q[12];
h q[13];
cx q[16],q[5];
h q[7];
h q[5];
h q[1];
h q[19];
cx q[19],q[3];
h q[9];
cx q[9],q[17];
cx q[1],q[14];
cx q[0],q[13];
h q[7];
h q[17];
h q[8];
cx q[9],q[17];
h q[7];
h q[10];
cx q[9],q[17];
h q[8];
cx q[16],q[5];
h q[16];
cx q[9],q[17];
cx q[4],q[7];
cx q[16],q[5];
cx q[9],q[17];
h q[12];
h q[8];
h q[12];
h q[8];
h q[10];
h q[8];