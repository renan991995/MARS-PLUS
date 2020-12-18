#ifndef IL_H
#define IL_H
#include <sys/types.h>
#include <dirent.h>
#include <cerrno>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>
#include <ctime>
#include <cstdio> 
#include <cctype>
#include <sstream>
#include "ATOM.h"
#include "MOLECULE.h"
#include "MATRIX.h"
#include "OPT.h"
#include "PARAMETER.h"


using namespace std;

class IL {
	public :
		IL() {
			smiles="";
			molesmi="";
			numIL[0]=numIL[1]=0;
			nsubcomp=2;

			for (int i=0;i<nsubcomp;i++) {
            	//ion[i].data.set_up();
				ion[i].data=NULL;
            	ion[i].chg=0;
            	ion[i].molesmi=="";
            	ion[i].parent=-1;
            	ion[i].if_circle=0;
				ion[i].natom=0;
            	
				ion[i].filenum=-10;
            	ion[i].clean();
            	ion[i].empty();
			}
			
		}
		~IL() {
			for (int i=0;i<nsubcomp;i++) {
				/*
                if (ion[i].data!=NULL) {
					delete ion[i].data;
					ion[i].data=NULL;
				}
				*/
                ion[i].clean();
                ion[i].empty();
            }
		}
		MOLECULE ion[2];
		string smiles;
		string molesmi;
		int nsubcomp;
		int numIL[2];
		
		int pair();
		int crossover(IL &il,int pp, int jj,int pp1, int jj1);
		int add(int pt, int id, int b,int pt1, int id1, int b1);
		int subtract(int n,int mode,int n1,int mode1);
		int combine(IL &il, int k,int p,int b,int k1,int p1,int b1);
		int exchange(int n,int id,int id2,int bond,int n_1,int id_1,int id2_1,int bond_1);
		int exchange1(int n,int id,int bnd2par,int bnd2des,int n_1,int id_1,int bnd2par_1,int bnd2des_1);
		int ring(int pt1,int pt2,int pt1_1,int pt2_1);
		int mds2smi();
		int replace(IL &);
		int reset();
		int prct();
		int input();
		int print(ofstream &outs);
		int canonicalize_SMILES();
};


#endif
