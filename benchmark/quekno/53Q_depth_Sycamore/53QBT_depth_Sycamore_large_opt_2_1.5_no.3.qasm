OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[52];
cx q[30],q[12];
h q[25];
h q[5];
h q[51];
h q[11];
cx q[24],q[51];
h q[29];
h q[6];
cx q[5],q[34];
cx q[6],q[37];
cx q[17],q[29];
cx q[15],q[35];
cx q[0],q[13];
h q[31];
cx q[25],q[38];
h q[5];
cx q[10],q[28];
h q[6];
cx q[26],q[10];
h q[30];
cx q[43],q[50];
h q[0];
cx q[7],q[48];
cx q[46],q[17];
cx q[21],q[27];
h q[34];
cx q[40],q[11];
cx q[21],q[27];
h q[29];
h q[26];
h q[25];
cx q[37],q[52];
h q[15];
h q[21];
h q[29];
h q[40];
h q[24];
cx q[37],q[31];
h q[27];
h q[17];
cx q[25],q[38];
h q[34];
h q[7];
h q[50];
h q[27];
h q[6];
h q[50];
h q[6];
cx q[0],q[13];
cx q[37],q[31];
h q[0];
cx q[26],q[10];
cx q[25],q[38];
h q[7];
h q[0];
h q[34];
h q[43];
h q[21];
cx q[24],q[51];
h q[35];
h q[15];
cx q[17],q[29];
cx q[7],q[48];
cx q[43],q[50];
h q[35];
cx q[7],q[48];
h q[40];
cx q[10],q[28];
h q[46];
h q[24];
h q[30];
h q[27];
h q[30];
cx q[43],q[50];
h q[51];
cx q[40],q[11];
cx q[15],q[35];
h q[5];
cx q[30],q[12];
h q[13];
h q[5];
cx q[15],q[35];
h q[13];
h q[52];
h q[52];
h q[26];
h q[21];
h q[46];
h q[3];
h q[28];
cx q[13],q[42];
h q[30];
cx q[32],q[33];
h q[13];
h q[52];
h q[24];
cx q[21],q[50];
cx q[6],q[52];
h q[50];
cx q[14],q[21];
h q[24];
cx q[15],q[50];
h q[30];
cx q[20],q[24];
h q[31];
cx q[18],q[35];
h q[28];
h q[51];
h q[6];
h q[33];
h q[32];
cx q[17],q[3];
cx q[49],q[31];
h q[21];
cx q[51],q[14];
cx q[28],q[20];
h q[17];
cx q[21],q[50];
h q[18];
h q[20];
cx q[11],q[46];
h q[15];
h q[52];
cx q[49],q[31];
h q[17];
h q[17];
cx q[30],q[5];
h q[33];
h q[3];
h q[49];
cx q[46],q[24];
h q[52];
h q[14];
h q[52];
cx q[11],q[46];
h q[24];
h q[15];
cx q[32],q[33];
h q[51];
h q[30];
cx q[14],q[21];
h q[18];
h q[33];
h q[51];
h q[32];
cx q[49],q[31];
cx q[30],q[5];
h q[50];
cx q[18],q[35];
h q[3];
cx q[18],q[35];
h q[6];
h q[13];
h q[13];
h q[28];
cx q[13],q[42];
h q[3];
h q[6];
cx q[14],q[22];
h q[9];
h q[9];
cx q[46],q[44];
cx q[15],q[43];
cx q[28],q[17];
cx q[20],q[14];
h q[3];
h q[43];
h q[46];
h q[29];
cx q[5],q[34];
h q[14];
cx q[7],q[26];
h q[46];
h q[29];
cx q[49],q[21];
cx q[19],q[23];
h q[28];
h q[16];
h q[46];
h q[26];
h q[31];
cx q[29],q[13];
cx q[9],q[16];
h q[20];
h q[47];
h q[49];
h q[25];
cx q[31],q[47];
h q[3];
cx q[25],q[10];
h q[6];
h q[19];
h q[15];
cx q[46],q[28];
h q[10];
h q[20];
h q[4];
cx q[28],q[17];
cx q[4],q[5];
cx q[3],q[24];
h q[43];
cx q[9],q[16];
h q[9];
cx q[6],q[34];
h q[19];
h q[3];
h q[31];
h q[9];
cx q[38],q[26];
cx q[3],q[24];
h q[9];
cx q[3],q[24];
cx q[6],q[34];
cx q[6],q[34];
cx q[28],q[17];
h q[15];
cx q[20],q[14];
cx q[3],q[24];
h q[44];
cx q[14],q[22];
h q[15];
h q[31];
cx q[49],q[21];
cx q[4],q[5];
h q[26];
h q[38];
cx q[46],q[44];
cx q[38],q[26];
h q[29];
cx q[20],q[14];
h q[21];
h q[44];
cx q[6],q[34];
h q[25];
h q[26];
cx q[4],q[5];
h q[29];
cx q[5],q[34];
h q[15];
h q[15];
cx q[6],q[34];
h q[49];
h q[19];
cx q[25],q[10];
h q[43];
h q[16];
h q[46];
h q[13];
cx q[31],q[47];
cx q[4],q[5];
cx q[15],q[43];
cx q[38],q[26];
cx q[19],q[23];
cx q[29],q[13];
cx q[49],q[21];
h q[23];
cx q[19],q[23];
cx q[25],q[10];
h q[7];
h q[49];
cx q[25],q[10];
h q[16];
h q[31];
h q[7];
cx q[31],q[47];
h q[25];
cx q[49],q[21];
h q[16];
cx q[29],q[13];
h q[22];
h q[7];
h q[10];
h q[7];
h q[7];
h q[7];
