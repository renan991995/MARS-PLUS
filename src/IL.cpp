#include "IL.h"
using namespace std;

extern PARAMETER para;

unsigned int IL::insertion(unsigned int n,unsigned int id,unsigned int bnd2par,unsigned int bnd2des,bool cistrans,unsigned int n1,unsigned int id1,unsigned int bnd2par1,unsigned int bnd2des1,bool cistrans1) {
    ion[0].insertion(n,id,bnd2par,bnd2des,cistrans);
    ion[1].insertion(n1,id1,bnd2par1,bnd2des1,cistrans1);
    return 1;
}

unsigned int IL::combination(IL &il,unsigned int k,unsigned int p,unsigned int b,unsigned int k1,unsigned int p1,unsigned int b1) {
	ion[0].combination(il.ion[0],k,p,b);
	ion[1].combination(il.ion[1],k1,p1,b1);
	return 1;
}

unsigned int IL::replace(IL &il) {
	ion[0].replace(il.ion[0]);
	ion[1].replace(il.ion[1]);

    smiles=il.smiles;
    molesmi=il.molesmi;
    nsubcomp=il.nsubcomp;
    frac=il.frac;


	return 1;
}


unsigned int IL::mds2smi(unsigned int chirality1,bool cistrans1,bool ct_on1,unsigned int chirality2,bool cistrans2,bool ct_on2) {
	ion[0].mds2smi(chirality1,cistrans1,ct_on1);
	ion[1].mds2smi(chirality2,cistrans2,ct_on2);
	return 1;
}


unsigned int IL::pair() {
	int i,chgg[2]={0,0};
	for (i=0;i<(int)ion[0].Cindex.size();i++) chgg[0]+=ion[0].data->a[ion[0].Mindex.at(i)].chg;
	for (i=0;i<(int)ion[1].Cindex.size();i++) chgg[1]+=ion[1].data->a[ion[1].Mindex.at(i)].chg;
	chgg[0]=fabs(chgg[0]);
	chgg[1]=fabs(chgg[1]);
	
	numIL[0]=(unsigned int)chgg[1];
	numIL[1]=(unsigned int)chgg[0];
	return 1;
}

unsigned int IL::crossover(IL &il,unsigned int pp, unsigned int jj,unsigned int pp1,unsigned int jj1) {
	ion[0].crossover(il.ion[0],pp,jj);
	ion[1].crossover(il.ion[1],pp1,jj1);
	return 1;
}


unsigned int IL::change_bnd(unsigned int n,unsigned int id,unsigned int id2,unsigned int bond,unsigned int n_1,unsigned int id_1,unsigned int id2_1,unsigned int bond_1) {
	ion[0].change_bnd(n,id,id2,bond);
	ion[1].change_bnd(n_1,id_1,id2_1,bond_1);
	return 1;
}

unsigned int IL::change_ele(unsigned int n,unsigned int id,unsigned int bnd2par,unsigned int bnd2des,unsigned int n_1,unsigned int id_1,unsigned int bnd2par_1,unsigned int bnd2des_1) {
    ion[0].change_ele(n,id,bnd2par,bnd2des);
    ion[1].change_ele(n_1,id_1,bnd2par_1,bnd2des_1);
    return 1;
}

unsigned int IL::addition(unsigned int pt,unsigned int id,unsigned int b,unsigned int pt1,unsigned int id1,unsigned int b1) {
	ion[0].addition(pt,id,b);
	ion[1].addition(pt1,id1,b1);
	return 1;
}

unsigned int IL::subtraction(unsigned int n,unsigned int bndfrm,unsigned int n1,unsigned int bndfrm1) {
	ion[0].subtraction(n,bndfrm);
	ion[1].subtraction(n1,bndfrm1);
	return 1;
}

unsigned int IL::cyclization(unsigned int pt1,unsigned int pt2,unsigned int bnd,unsigned int pt1_1,unsigned int pt2_1,unsigned int bnd_1) {
	ion[0].cyclization(pt1,pt2,bnd);
	ion[1].cyclization(pt1_1,pt2_1,bnd_1);
	return 1;
}

unsigned int IL::decyclization(unsigned int ring1,unsigned int ring2) {
    ion[0].decyclization(ring1);
    ion[1].decyclization(ring2);
    return 1;
}

unsigned int IL::prct() {
    ion[0].prct();
    ion[1].prct();
    return 1;
}

unsigned int IL::reset() {
    ion[0].reset();
    ion[1].reset();
    return 1;
}

unsigned int IL::input() {
    ion[0].input();
    ion[1].input();
    return 1;
}

unsigned int IL::printmds(ofstream &outs) {
    ion[0].printmds(outs);
    ion[1].printmds(outs);
    return 1;
}


unsigned int IL::canonicalize_SMILES() {
	ion[0].canonicalize_SMILES();
	ion[1].canonicalize_SMILES();
	return 1;
}

unsigned int IL::change_cistrans(unsigned int pos1,unsigned int w1,unsigned int pos2,unsigned int w2) {
    ion[0].change_cistrans(pos1,w1);
    ion[1].change_cistrans(pos2,w2);
    return 1;
}

unsigned int IL::change_chirality(unsigned int pos1,unsigned int pos2) {
    ion[0].change_chirality(pos1);
    ion[1].change_chirality(pos2);
    return 1;
}

