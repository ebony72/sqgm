OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
h q[14];
h q[3];
h q[2];
cx q[3],q[11];
h q[7];
cx q[2],q[15];
h q[4];
h q[1];
h q[3];
cx q[11],q[4];
cx q[3],q[11];
cx q[14],q[5];
h q[5];
cx q[9],q[19];
h q[19];
h q[12];
h q[4];
h q[4];
cx q[9],q[7];
h q[1];
h q[3];
cx q[1],q[5];
cx q[4],q[12];
h q[3];
cx q[15],q[0];
h q[9];
h q[2];
h q[2];
cx q[15],q[0];
cx q[9],q[7];
cx q[15],q[0];
cx q[15],q[0];
h q[11];
cx q[9],q[19];
h q[7];
cx q[9],q[7];
h q[19];
h q[1];
cx q[14],q[5];
h q[11];
h q[5];
h q[2];
h q[14];
h q[2];
h q[1];
h q[6];
h q[11];
h q[3];
cx q[3],q[12];
cx q[11],q[6];
h q[14];
h q[17];
h q[1];
cx q[7],q[5];
cx q[15],q[3];
cx q[17],q[15];
h q[3];
cx q[16],q[4];
cx q[14],q[19];
h q[6];
cx q[5],q[13];
h q[19];
h q[7];
h q[11];
cx q[1],q[12];
h q[14];
h q[19];
h q[12];
h q[12];
h q[19];
cx q[0],q[17];
cx q[5],q[3];
h q[13];
h q[11];
h q[0];
cx q[5],q[19];
h q[14];
cx q[5],q[19];
cx q[1],q[12];
cx q[0],q[17];
h q[14];
cx q[16],q[4];
h q[16];
cx q[3],q[12];
h q[13];
h q[15];
h q[14];
cx q[16],q[4];
h q[16];
h q[16];
h q[1];
h q[14];
h q[7];
h q[16];
cx q[11],q[6];
h q[4];
cx q[11],q[6];
h q[6];
h q[13];
h q[11];
h q[7];
h q[4];
h q[7];
h q[4];
h q[7];
h q[13];
h q[15];
h q[7];
h q[15];
h q[13];
h q[0];
h q[11];
h q[19];
h q[1];
h q[19];
cx q[1],q[10];
h q[1];
cx q[1],q[19];
h q[3];
h q[0];
h q[19];
h q[14];
h q[11];
h q[10];
cx q[14],q[3];
cx q[15],q[5];
cx q[5],q[14];
cx q[4],q[1];
h q[6];
cx q[9],q[5];
cx q[10],q[6];
h q[6];
cx q[19],q[10];
h q[19];
cx q[0],q[9];
h q[4];
h q[19];
cx q[11],q[6];
cx q[11],q[6];
h q[4];
h q[15];
h q[11];
h q[6];
cx q[11],q[6];
cx q[1],q[10];
cx q[1],q[19];
h q[10];
cx q[5],q[14];
h q[10];
cx q[14],q[3];
h q[14];
h q[3];
h q[4];
h q[5];
cx q[14],q[3];
cx q[0],q[9];
cx q[0],q[9];
h q[0];
cx q[15],q[5];
h q[0];
cx q[15],q[5];
cx q[15],q[5];
h q[14];
h q[4];
h q[3];
h q[9];
h q[9];
h q[19];
h q[6];
cx q[9],q[3];
cx q[4],q[19];
h q[9];
h q[16];
h q[16];
h q[2];
h q[13];
h q[2];
cx q[0],q[16];
cx q[15],q[6];
h q[12];
cx q[2],q[9];
h q[13];
cx q[5],q[7];
cx q[9],q[8];
cx q[13],q[10];
h q[4];
cx q[1],q[12];
h q[7];
cx q[19],q[6];
h q[7];
h q[19];
h q[12];
cx q[9],q[3];
h q[5];
h q[8];
h q[5];
h q[5];
cx q[2],q[9];
h q[4];
cx q[6],q[10];
h q[8];
h q[12];
h q[19];
cx q[0],q[16];
cx q[4],q[19];
h q[13];
cx q[6],q[10];
cx q[1],q[12];
h q[3];
h q[6];
h q[15];
cx q[5],q[7];
cx q[5],q[7];
h q[13];
cx q[13],q[10];
cx q[0],q[16];
h q[1];
cx q[0],q[16];
h q[15];
h q[12];
h q[15];
h q[15];
h q[18];
h q[12];
h q[1];
h q[1];
cx q[14],q[0];
h q[18];
cx q[5],q[12];
cx q[0],q[9];
h q[6];
cx q[6],q[13];
h q[5];
h q[13];
h q[4];
cx q[4],q[6];
h q[4];
cx q[18],q[10];
h q[6];
h q[1];
cx q[3],q[7];
cx q[1],q[7];
h q[5];
cx q[6],q[18];
h q[4];
h q[6];
h q[13];
cx q[1],q[3];
cx q[1],q[7];
h q[13];
h q[10];
cx q[0],q[9];
h q[0];
h q[0];
cx q[0],q[9];
cx q[18],q[10];
cx q[5],q[12];
h q[3];
h q[14];
cx q[5],q[12];
h q[13];
h q[14];
h q[4];
h q[14];
h q[14];
h q[14];
cx q[14],q[7];
h q[19];
cx q[4],q[11];
h q[6];
h q[15];
cx q[0],q[15];
h q[10];
h q[11];
h q[1];
cx q[7],q[6];
cx q[2],q[16];
h q[17];
cx q[7],q[12];
h q[17];
cx q[6],q[10];
cx q[1],q[0];
h q[16];
h q[16];
h q[11];
cx q[19],q[7];
cx q[17],q[10];
h q[15];
h q[6];
h q[6];
h q[4];
h q[17];
h q[15];
cx q[0],q[15];
h q[16];
cx q[7],q[12];
cx q[2],q[16];
cx q[4],q[11];
cx q[4],q[11];
h q[1];
h q[19];
h q[14];
h q[14];
h q[14];
h q[10];
h q[1];
h q[14];
