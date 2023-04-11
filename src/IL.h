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
#include "ELEMENTS.h"
#include "MOLECULE.h"
#include "PARAMETER.h"
using namespace std;

class IL {
	public :
		IL() {
            smiles="";
            molesmi="";
            numIL[0]=numIL[1]=0;

			comp_id=1;
            nsubcomp=2;
            frac=1.000000;

			for (unsigned int i=0;i<nsubcomp;i++) {

				ion[i].data=NULL;
            	ion[i].chg=0;
            	ion[i].smiles=ion[i].molesmi=="";
            	ion[i].parent=-1;
            	ion[i].if_circle=0;
				ion[i].natom=0;

                ion[i].frac=1.0/nsubcomp;
				ion[i].nsubcomp=i+1;
            	
            	ion[i].clean();
            	ion[i].empty();
			}
			
		}
		~IL() {
			for (unsigned int i=0;i<nsubcomp;i++) {
                ion[i].clean();
                ion[i].empty();
            }
		}
		MOLECULE ion[2];
		string smiles;
		string molesmi;
		double frac;
		unsigned int nsubcomp;
		unsigned int comp_id;
		unsigned int numIL[2];
		
		unsigned int pair();
		unsigned int crossover(IL &il,unsigned int pp,unsigned int jj,unsigned int pp1,unsigned int jj1);
		unsigned int addition(unsigned int pt,unsigned int id,unsigned int b,unsigned int pt1,unsigned int id1,unsigned int b1);
		unsigned int insertion(unsigned int n,unsigned int id,unsigned int bnd2par,unsigned int bnd2des,bool cistrans,unsigned int n1,unsigned int id1,unsigned int bnd2par1,unsigned int bnd2des1,bool cistrans1);
		unsigned int subtraction(unsigned int n,unsigned int bndfrm,unsigned int n1,unsigned int bndfrm1);
		unsigned int combination(IL &il,unsigned int k,unsigned int p,unsigned int b,unsigned int k1,unsigned int p1,unsigned int b1);
		unsigned int change_bnd(unsigned int n,unsigned int id,unsigned int id2,unsigned int bond,unsigned int n_1,unsigned int id_1,unsigned int id2_1,unsigned int bond_1);
		unsigned int change_ele(unsigned int n,unsigned int id,unsigned int bnd2par,unsigned int bnd2des,unsigned int n_1,unsigned int id_1,unsigned int bnd2par_1,unsigned int bnd2des_1);
		unsigned int cyclization(unsigned int pt1,unsigned int pt2,unsigned int bnd,unsigned int pt1_1,unsigned int pt2_1,unsigned int bnd_1);
		unsigned int decyclization(unsigned int ring1,unsigned int ring2);
		unsigned int mds2smi(unsigned int chirality1=2,bool cistrans1=0,bool ct_on1=0,unsigned int chirality2=2,bool cistrans2=0,bool ct_on2=0);
		unsigned int replace(IL &);
		unsigned int reset();
		unsigned int prct();
		unsigned int input();
		unsigned int printmds(ofstream &outs);
		unsigned int canonicalize_SMILES();
        unsigned int change_chirality(unsigned int pos1,unsigned int pos2);
        unsigned int change_cistrans(unsigned int pos1,unsigned int w1,unsigned int pos2,unsigned int w2);
};


#endif
