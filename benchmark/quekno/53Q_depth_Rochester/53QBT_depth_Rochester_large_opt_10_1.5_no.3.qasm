OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[30];
h q[47];
h q[31];
h q[30];
h q[12];
h q[35];
cx q[31],q[30];
h q[19];
cx q[7],q[19];
cx q[46],q[44];
cx q[9],q[26];
h q[31];
h q[1];
h q[26];
h q[7];
h q[46];
h q[11];
cx q[51],q[40];
h q[0];
h q[46];
h q[51];
h q[46];
cx q[22],q[51];
cx q[6],q[28];
h q[7];
h q[26];
h q[26];
h q[26];
cx q[1],q[12];
cx q[35],q[21];
cx q[19],q[8];
cx q[38],q[47];
cx q[31],q[30];
cx q[0],q[37];
h q[22];
h q[31];
h q[32];
h q[6];
h q[8];
cx q[29],q[7];
cx q[41],q[11];
cx q[20],q[19];
cx q[32],q[36];
h q[26];
cx q[20],q[19];
h q[30];
h q[8];
cx q[41],q[11];
h q[0];
cx q[41],q[11];
cx q[7],q[19];
cx q[22],q[51];
h q[12];
h q[0];
h q[51];
h q[11];
h q[44];
cx q[38],q[47];
h q[0];
h q[36];
h q[35];
h q[47];
cx q[41],q[11];
cx q[32],q[36];
cx q[38],q[47];
h q[47];
h q[8];
h q[46];
h q[46];
cx q[6],q[28];
h q[38];
h q[1];
h q[35];
h q[0];
cx q[32],q[36];
cx q[35],q[21];
cx q[6],q[28];
cx q[32],q[36];
cx q[35],q[21];
cx q[1],q[12];
h q[22];
h q[12];
h q[1];
h q[9];
h q[29];
cx q[6],q[28];
h q[40];
cx q[6],q[28];
cx q[1],q[12];
h q[40];
h q[37];
h q[40];
h q[44];
h q[9];
h q[37];
h q[9];
h q[20];
h q[40];
h q[9];
h q[9];
h q[37];
h q[44];
h q[40];
h q[44];
h q[44];
h q[37];
h q[21];
h q[8];
h q[24];
h q[16];
h q[2];
cx q[20],q[8];
h q[6];
h q[19];
cx q[24],q[7];
h q[14];
h q[0];
cx q[21],q[2];
h q[38];
cx q[29],q[8];
h q[8];
cx q[14],q[20];
cx q[38],q[44];
h q[24];
cx q[46],q[49];
cx q[0],q[22];
h q[6];
h q[8];
h q[46];
h q[14];
h q[26];
cx q[8],q[19];
cx q[46],q[49];
cx q[16],q[24];
cx q[22],q[27];
h q[0];
h q[20];
h q[1];
cx q[34],q[42];
cx q[1],q[52];
cx q[14],q[20];
h q[42];
cx q[4],q[5];
h q[7];
h q[44];
cx q[6],q[26];
h q[8];
h q[34];
cx q[20],q[8];
h q[49];
h q[19];
h q[1];
h q[2];
cx q[7],q[29];
h q[1];
cx q[29],q[8];
h q[14];
h q[27];
cx q[6],q[26];
cx q[21],q[2];
h q[19];
h q[38];
h q[49];
h q[49];
h q[34];
cx q[0],q[22];
cx q[34],q[42];
h q[1];
cx q[34],q[42];
cx q[4],q[5];
cx q[22],q[27];
cx q[22],q[27];
h q[49];
cx q[1],q[52];
cx q[34],q[42];
h q[6];
cx q[22],q[27];
h q[49];
h q[44];
h q[42];
cx q[46],q[49];
cx q[6],q[26];
h q[26];
h q[38];
cx q[6],q[26];
cx q[24],q[7];
h q[42];
h q[42];
cx q[14],q[20];
cx q[22],q[27];
h q[6];
cx q[16],q[24];
cx q[0],q[22];
cx q[16],q[24];
h q[5];
h q[44];
h q[7];
h q[26];
h q[52];
cx q[38],q[44];
h q[19];
h q[16];
cx q[1],q[52];
h q[44];
cx q[24],q[7];
h q[2];
cx q[1],q[52];
h q[4];
h q[38];
cx q[16],q[24];
h q[5];
cx q[38],q[44];
h q[27];
h q[7];
cx q[21],q[2];
h q[2];
cx q[21],q[2];
cx q[21],q[2];
cx q[4],q[5];
h q[5];
h q[44];
h q[5];
h q[34];
h q[38];
cx q[4],q[5];
cx q[4],q[5];
h q[34];
h q[34];
h q[34];
cx q[7],q[19];
h q[40];
cx q[48],q[44];
cx q[34],q[23];
h q[19];
h q[40];
h q[19];
cx q[15],q[21];
h q[20];
h q[40];
cx q[45],q[20];
cx q[42],q[11];
cx q[50],q[43];
h q[6];
cx q[27],q[51];
h q[24];
cx q[37],q[24];
cx q[39],q[32];
h q[6];
h q[21];
h q[19];
h q[19];
h q[43];
h q[7];
h q[24];
h q[44];
h q[21];
h q[20];
h q[20];
cx q[43],q[12];
cx q[45],q[14];
h q[19];
cx q[3],q[37];
h q[32];
cx q[3],q[37];
cx q[24],q[16];
cx q[6],q[13];
cx q[45],q[14];
h q[11];
h q[50];
cx q[3],q[37];
cx q[41],q[11];
h q[7];
h q[15];
cx q[51],q[40];
h q[44];
cx q[48],q[44];
h q[34];
h q[43];
h q[20];
h q[12];
h q[24];
h q[15];
cx q[43],q[12];
h q[13];
h q[44];
cx q[50],q[43];
h q[14];
h q[23];
h q[32];
h q[41];
cx q[24],q[16];
h q[20];
cx q[45],q[14];
h q[6];
cx q[48],q[44];
cx q[15],q[21];
cx q[39],q[32];
h q[37];
h q[6];
cx q[41],q[11];
h q[13];
cx q[42],q[11];
cx q[34],q[23];
h q[41];
h q[51];
cx q[6],q[13];
h q[40];
cx q[15],q[21];
cx q[27],q[51];
cx q[34],q[23];
h q[32];
cx q[34],q[23];
cx q[39],q[32];
cx q[15],q[21];
h q[40];
h q[7];
h q[7];
h q[7];
h q[12];
h q[3];
cx q[50],q[12];
cx q[20],q[8];
cx q[51],q[35];
h q[12];
h q[25];
h q[3];
cx q[29],q[20];
h q[36];
cx q[39],q[33];
h q[33];
h q[36];
h q[3];
cx q[27],q[40];
h q[19];
h q[50];
h q[22];
cx q[3],q[16];
cx q[52],q[50];
h q[50];
h q[50];
h q[29];
h q[18];
cx q[25],q[49];
h q[19];
h q[36];
h q[27];
h q[36];
h q[40];
h q[51];
h q[14];
h q[22];
cx q[0],q[22];
cx q[18],q[36];
h q[25];
cx q[4],q[5];
h q[6];
h q[50];
h q[3];
cx q[32],q[39];
cx q[39],q[33];
cx q[19],q[20];
h q[50];
cx q[27],q[40];
h q[25];
cx q[40],q[51];
cx q[31],q[14];
cx q[6],q[30];
h q[30];
h q[36];
cx q[3],q[16];
cx q[18],q[36];
cx q[4],q[5];
cx q[29],q[20];
cx q[39],q[33];
h q[30];
h q[40];
h q[32];
h q[12];
h q[32];
h q[22];
h q[51];
h q[29];
h q[5];
h q[22];
h q[20];
cx q[31],q[14];
h q[31];
h q[20];
h q[30];
h q[25];
h q[52];
h q[19];
cx q[32],q[39];
h q[14];
h q[30];
h q[16];
h q[29];
cx q[3],q[16];
h q[32];
cx q[25],q[49];
cx q[39],q[33];
cx q[4],q[5];
cx q[32],q[39];
cx q[6],q[30];
cx q[40],q[51];
h q[19];
cx q[27],q[40];
h q[51];
h q[29];
h q[0];
cx q[20],q[8];
cx q[31],q[14];
h q[35];
cx q[31],q[14];
cx q[0],q[22];
cx q[0],q[22];
h q[52];
cx q[25],q[49];
cx q[4],q[5];
h q[33];
cx q[31],q[14];
h q[52];
cx q[4],q[5];
cx q[4],q[5];
h q[12];
h q[12];
h q[19];
h q[35];
h q[12];
h q[52];
h q[19];
h q[35];
h q[35];
h q[12];
h q[35];
h q[35];
cx q[40],q[17];
h q[27];
cx q[36],q[32];
h q[32];
h q[5];
h q[0];
cx q[13],q[50];
h q[0];
cx q[0],q[27];
cx q[32],q[33];
h q[17];
h q[13];
h q[27];
cx q[52],q[43];
cx q[34],q[42];
h q[40];
h q[44];
cx q[16],q[7];
h q[31];
cx q[12],q[41];
h q[5];
cx q[37],q[16];
h q[31];
cx q[12],q[5];
h q[36];
h q[36];
h q[7];
h q[16];
h q[16];
h q[0];
h q[42];
cx q[44],q[1];
h q[0];
h q[32];
cx q[14],q[31];
cx q[51],q[40];
h q[34];
h q[0];
h q[43];
cx q[37],q[16];
cx q[51],q[40];
h q[32];
cx q[4],q[12];
h q[52];
cx q[16],q[7];
cx q[51],q[40];
cx q[52],q[43];
cx q[51],q[40];
h q[14];
cx q[13],q[50];
cx q[12],q[5];
cx q[32],q[33];
h q[50];
cx q[4],q[12];
h q[42];
h q[42];
h q[27];
cx q[34],q[42];
h q[27];
cx q[34],q[42];
cx q[13],q[50];
h q[52];
h q[14];
cx q[14],q[31];
cx q[44],q[1];
h q[37];
h q[36];
h q[5];
h q[1];
h q[52];
cx q[44],q[1];
h q[17];
cx q[13],q[50];
h q[44];
h q[36];
h q[41];
h q[17];
h q[41];
h q[52];
h q[41];
h q[1];
h q[17];
h q[17];
h q[43];
h q[36];
h q[41];
h q[43];
h q[41];
h q[43];
cx q[14],q[1];
cx q[19],q[47];
cx q[28],q[50];
h q[19];
cx q[18],q[7];
h q[13];
h q[20];
cx q[5],q[41];
cx q[20],q[2];
h q[20];
h q[47];
h q[19];
h q[17];
h q[18];
cx q[24],q[40];
cx q[22],q[17];
cx q[40],q[25];
h q[22];
h q[36];
h q[18];
h q[42];
h q[47];
h q[50];
cx q[32],q[18];
h q[17];
cx q[52],q[4];
h q[24];
cx q[36],q[33];
h q[7];
h q[41];
cx q[41],q[42];
h q[41];
h q[22];
h q[5];
h q[52];
h q[47];
h q[33];
h q[18];
cx q[31],q[14];
h q[19];
h q[19];
cx q[41],q[42];
cx q[28],q[13];
cx q[31],q[14];
h q[52];
cx q[38],q[46];
h q[52];
cx q[19],q[47];
cx q[52],q[4];
h q[18];
cx q[40],q[25];
cx q[40],q[25];
h q[25];
h q[42];
cx q[20],q[2];
h q[5];
cx q[20],q[2];
h q[22];
h q[31];
h q[28];
h q[13];
cx q[31],q[14];
h q[33];
cx q[31],q[14];
h q[7];
h q[20];
cx q[22],q[17];
cx q[5],q[41];
cx q[24],q[40];
cx q[38],q[46];
h q[38];
h q[40];
cx q[52],q[4];
h q[25];
h q[1];
h q[1];
h q[7];
cx q[28],q[13];
cx q[28],q[13];
cx q[28],q[13];
cx q[38],q[46];
cx q[38],q[46];
h q[36];
cx q[38],q[46];
h q[2];
h q[7];
cx q[36],q[33];
cx q[36],q[33];
h q[32];
h q[32];
h q[24];
h q[7];
h q[1];
h q[1];
h q[1];
h q[50];
h q[50];
h q[50];
h q[50];
cx q[36],q[45];
h q[45];
h q[32];
h q[24];
h q[15];
cx q[32],q[7];
cx q[29],q[20];
h q[24];
h q[20];
cx q[11],q[41];
h q[44];
h q[45];
cx q[48],q[23];
h q[24];
cx q[24],q[15];
h q[1];
h q[13];
cx q[43],q[46];
h q[45];
h q[24];
cx q[44],q[50];
h q[15];
h q[9];
cx q[9],q[1];
cx q[13],q[52];
h q[29];
cx q[44],q[50];
h q[9];
h q[25];
h q[41];
cx q[2],q[20];
cx q[52],q[50];
h q[7];
cx q[9],q[1];
h q[41];
cx q[46],q[48];
cx q[15],q[25];
h q[11];
h q[20];
cx q[43],q[44];
h q[11];
cx q[43],q[44];
h q[20];
h q[41];
h q[50];
h q[9];
h q[13];
h q[13];
cx q[46],q[48];
h q[41];
h q[11];
cx q[44],q[50];
h q[1];
cx q[46],q[48];
h q[7];
h q[11];
cx q[2],q[20];
h q[1];
h q[23];
cx q[46],q[48];
cx q[46],q[48];
cx q[13],q[52];
cx q[13],q[52];
cx q[32],q[7];
cx q[11],q[41];
h q[29];
cx q[32],q[7];
cx q[36],q[45];
cx q[36],q[45];
h q[23];
h q[29];
h q[43];
h q[24];
h q[29];
h q[9];
h q[23];
h q[23];
h q[29];
h q[23];
h q[47];
h q[19];
cx q[29],q[19];
cx q[45],q[2];
h q[39];
h q[47];
cx q[26],q[28];
h q[1];
h q[50];
cx q[45],q[39];
cx q[15],q[21];
cx q[8],q[48];
h q[50];
cx q[46],q[44];
h q[28];
h q[26];
h q[15];
h q[34];
cx q[23],q[11];
cx q[28],q[50];
h q[29];
h q[19];
h q[50];
cx q[33],q[31];
h q[45];
h q[26];
h q[21];
cx q[1],q[9];
h q[8];
h q[33];
cx q[34],q[23];
cx q[49],q[47];
cx q[24],q[22];
h q[31];
h q[49];
cx q[33],q[31];
h q[47];
h q[28];
h q[50];
cx q[20],q[29];
h q[22];
cx q[15],q[21];
cx q[46],q[44];
cx q[49],q[47];
h q[45];
h q[15];
cx q[29],q[19];
cx q[1],q[9];
h q[39];
cx q[46],q[44];
cx q[45],q[39];
cx q[8],q[48];
cx q[24],q[22];
cx q[23],q[11];
cx q[8],q[48];
h q[31];
h q[46];
h q[24];
h q[48];
cx q[46],q[44];
h q[20];
h q[1];
cx q[24],q[22];
h q[31];
h q[28];
h q[2];
cx q[1],q[9];
cx q[23],q[11];
cx q[15],q[21];
cx q[34],q[23];
h q[11];
h q[8];
h q[26];
h q[2];
h q[33];
h q[2];
h q[33];
h q[26];
h q[2];
h q[6];
h q[45];
cx q[26],q[50];
h q[43];
h q[7];
cx q[49],q[48];
h q[4];
h q[51];
cx q[23],q[42];
cx q[44],q[43];
h q[38];
cx q[41],q[4];
h q[41];
cx q[0],q[51];
cx q[36],q[45];
cx q[24],q[25];
h q[6];
h q[4];
cx q[6],q[30];
h q[24];
h q[11];
h q[31];
cx q[38],q[23];
h q[7];
h q[11];
h q[36];
h q[11];
h q[44];
cx q[30],q[10];
cx q[41],q[4];
h q[45];
cx q[45],q[31];
cx q[52],q[5];
h q[11];
h q[49];
h q[31];
h q[36];
h q[25];
h q[43];
h q[51];
cx q[7],q[35];
cx q[38],q[23];
cx q[42],q[11];
cx q[52],q[5];
cx q[52],q[5];
cx q[30],q[10];
h q[38];
h q[43];
cx q[52],q[5];
cx q[26],q[50];
cx q[7],q[35];
h q[0];
cx q[7],q[35];
cx q[24],q[25];
h q[49];
h q[6];
cx q[52],q[5];
cx q[44],q[43];
cx q[26],q[50];
h q[0];
h q[6];
cx q[41],q[4];
h q[51];
h q[51];
cx q[38],q[23];
cx q[24],q[25];
h q[0];
cx q[49],q[48];
h q[36];
h q[25];
h q[26];
h q[49];
h q[50];
cx q[26],q[50];
h q[45];
h q[48];
h q[24];
cx q[44],q[50];
h q[24];
h q[40];
cx q[8],q[34];
h q[39];
cx q[28],q[13];
h q[38];
cx q[19],q[8];
cx q[32],q[40];
h q[6];
h q[50];
h q[31];
h q[39];
cx q[6],q[10];
h q[38];
cx q[34],q[48];
cx q[19],q[8];
h q[8];
h q[15];
h q[33];
h q[19];
h q[31];
cx q[21],q[15];
h q[21];
h q[44];
h q[19];
cx q[38],q[23];
cx q[25],q[24];
h q[6];
cx q[34],q[48];
cx q[23],q[42];
h q[8];
cx q[33],q[18];
cx q[39],q[31];
cx q[0],q[3];
h q[24];
h q[10];
h q[40];
h q[48];
h q[42];
h q[0];
cx q[19],q[8];
h q[0];
cx q[7],q[32];
h q[0];
h q[38];
cx q[33],q[18];
h q[3];
cx q[0],q[3];
cx q[25],q[24];
cx q[38],q[23];
cx q[7],q[32];
cx q[28],q[13];
cx q[32],q[40];
h q[15];
cx q[7],q[32];
h q[6];
cx q[6],q[10];
cx q[21],q[15];
cx q[23],q[42];
cx q[6],q[10];
h q[18];
cx q[39],q[31];
h q[31];
cx q[0],q[3];
h q[48];
h q[50];
cx q[21],q[15];
cx q[25],q[24];
cx q[21],q[15];
cx q[44],q[50];
h q[44];
cx q[33],q[18];
cx q[39],q[31];
h q[34];
cx q[33],q[18];
cx q[25],q[24];
h q[48];
cx q[44],q[50];
cx q[28],q[13];
h q[38];
cx q[28],q[13];
cx q[28],q[13];
h q[34];
h q[13];
h q[40];
h q[28];
h q[34];
h q[2];
cx q[22],q[24];
cx q[2],q[31];
h q[34];
h q[39];
cx q[35],q[19];
h q[22];
cx q[10],q[14];
h q[11];
h q[25];
cx q[41],q[4];
h q[35];
h q[2];
cx q[36],q[2];
h q[17];
h q[34];
h q[20];
cx q[29],q[34];
cx q[23],q[11];
cx q[40],q[20];
h q[2];
h q[10];
h q[25];
cx q[50],q[41];
h q[39];
cx q[41],q[4];
h q[50];
cx q[24],q[25];
h q[24];
h q[33];
h q[20];
h q[14];
h q[35];
cx q[7],q[33];
cx q[15],q[21];
cx q[41],q[12];
cx q[7],q[33];
cx q[17],q[22];
h q[19];
h q[36];
cx q[22],q[24];
cx q[2],q[39];
h q[14];
h q[41];
h q[34];
h q[40];
h q[21];
h q[31];
cx q[10],q[14];
h q[34];
cx q[15],q[21];
cx q[10],q[14];
cx q[23],q[11];
h q[14];
cx q[24],q[25];
cx q[7],q[33];
cx q[29],q[34];
cx q[35],q[19];
cx q[17],q[22];
cx q[40],q[20];
cx q[50],q[41];
cx q[40],q[20];
cx q[7],q[33];
cx q[40],q[20];
cx q[7],q[33];
h q[31];
cx q[15],q[21];
h q[31];
h q[36];
h q[35];
h q[35];
h q[19];
h q[12];
cx q[15],q[21];
cx q[23],q[11];
cx q[15],q[21];
h q[31];
cx q[23],q[11];
h q[11];
h q[12];
h q[19];
h q[10];
h q[4];
h q[4];
h q[23];
h q[4];
