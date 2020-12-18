#include "IL.h"
using namespace std;

extern PARAMETER para;


int IL::combine(IL &il, int k,int p,int b, int k1,int p1,int b1) {
	ion[0].combine(il.ion[0],k,p,b);
	ion[1].combine(il.ion[1],k1,p1,b1);
	return 1;
}


int IL::replace(IL &il) {
	ion[0].replace(il.ion[0]);
	ion[1].replace(il.ion[1]);

    smiles=il.smiles;
    molesmi=il.molesmi;
    nsubcomp=2;

	return 1;
}


int IL::mds2smi() {
	ion[0].mds2smi();
	ion[1].mds2smi();
	return 1;
}


int IL::pair() {
	int i,chgg[2]={0,0};
	for (i=0;i<ion[0].Cindex.size();i++) chgg[0]+=ion[0].data->a[ion[0].Mindex.at(i)].chg;
	for (i=0;i<ion[1].Cindex.size();i++) chgg[1]+=ion[1].data->a[ion[1].Mindex.at(i)].chg;
	chgg[0]=fabs(chgg[0]);
	chgg[1]=fabs(chgg[1]);
	
	numIL[0]=chgg[1];
	numIL[1]=chgg[0];
	return 1;
}


int IL::crossover(IL &il,int pp, int jj,int pp1, int jj1) {
	ion[0].crossover(il.ion[0],pp,jj);
	ion[1].crossover(il.ion[1],pp1,jj1);
	return 1;
}

int IL::exchange(int n,int id,int id2,int bond,int n_1,int id_1,int id2_1,int bond_1) {
	ion[0].exchange(n,id,id2,bond);
	ion[1].exchange(n_1,id_1,id2_1,bond_1);
	return 1;
}

int IL::exchange1(int n,int id,int bnd2par,int bnd2des,int n_1,int id_1,int bnd2par_1,int bnd2des_1) {
    ion[0].exchange(n,id,bnd2par,bnd2des);
    ion[1].exchange(n_1,id_1,bnd2par_1,bnd2des_1);
    return 1;
}

int IL::add(int pt, int id,int b,int pt1, int id1, int b1) {
	ion[0].add(pt,id,b);
	ion[1].add(pt1,id1,b1);
	return 1;
}


int IL::subtract(int n,int mode,int n1,int mode1) {
	ion[0].subtract(n,mode);
	ion[1].subtract(n1,mode1);
	return 1;
}


int IL::ring(int pt1,int pt2,int pt1_1,int pt2_1) {
	ion[0].ring(pt1,pt2);
	ion[1].ring(pt1_1,pt2_1);
	return 1;
}

int IL::prct() {
    ion[0].prct();
    ion[1].prct();
    return 1;
}

int IL::reset() {
    ion[0].reset();
    ion[1].reset();
    return 1;
}

int IL::input() {
    ion[0].input();
    ion[1].input();
    return 1;
}

int IL::print(ofstream &outs) {
    ion[0].print(outs);
    ion[1].print(outs);
    return 1;
}


int IL::canonicalize_SMILES() {
	ion[0].canonicalize_SMILES();
	ion[1].canonicalize_SMILES();
	return 1;
}


