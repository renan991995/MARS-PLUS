/***************************************************************************************
                          MARS for Computer-aided Molecular Design
                                      Hsuan Hao Hsu
                                    Chen Hsuan Huang
                            Shiang Tai Lin (stlin@ntu.edu.tw)
      Department of Chemical Engineering, National Taiwan University, Taipei, Taiwan
                                   Copyright (c) 2019
                                   All rights Reserved
***************************************************************************************/
#ifndef MOLECULE_H
#define MOLECULE_H
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/data.h>
#include <openbabel/data_utilities.h>
#include <openbabel/obiter.h>
#include <openbabel/chiral.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/stereo/tetranonplanar.h>
#include <openbabel/stereo/tetraplanar.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/builder.h>
#include <openbabel/typer.h>
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
#include "MATRIX.h"
#include "OPT.h"
#include "PARAMETER.h"
//#include "IL.h"
using namespace std;
using namespace OpenBabel;

class PROPERTIES {
	public :
		PROPERTIES () {
            front=1;
            crowdist=0;
			Xco2=-999999;
			Pvap=-999999;
            lnKow=-999999;
			Henry=-999999;
            HOMO=-999999;
            LUMO=-999999;
            IP=-999999;
            EA=-999999;
            Egap=-999999;
            hardness=-999999;
            ephilic=-999999;
			lnKow=-999999;
            fitness=0;
            p_sel=0;
			P=100000;
			T=298.15;
        }
		~PROPERTIES() {}
		double crowdist;
		int front;
		double Xco2;
		double Pvap;
		double lnKow;
		double Henry;
		double HOMO;
		double LUMO;
		double IP;
		double EA;
		double hardness;
		double Egap;
		double ephilic;
		double fitness;
		double p_sel;
		double T;  //0502
		double P;  //0502
};	

class CALCflag : public PROPERTIES {
	public:
		CALCflag () {
            front=1;
            crowdist=0;
			Xco2=-999999;
			Pvap=-999999;
            lnKow=-999999;
			Henry=-999999;
            HOMO=-999999;
            LUMO=-999999;
            IP=-999999;
            EA=-999999;
            Egap=-999999;
            hardness=-999999;
            ephilic=-999999;
			lnKow=-999999;
            fitness=0;
            p_sel=0;			
            P=100000;
            T=298.15;
			
			filenum=-10;
			nosubmitscribd=-10;
			sameasancestor=0;
			tosubmit=1;
			frac=1.000000;
			tcal=0;
		}
		int nsubcomp;
		int comp_id;
		int filenum;
		int nosubmitscribd;
		double frac;
		bool sameasancestor;
		bool tosubmit;
		int tcal;
};
		
class MOLECULE : public CALCflag
{
	public:
		MOLECULE() {
            front=1;
            crowdist=0;
            lnKow=-999999;
			Xco2=-999999;
			Pvap=-999999;
			Henry=-999999;
            HOMO=-999999;
            LUMO=-999999;
            IP=-999999;
            EA=-999999;
            Egap=-999999;
            hardness=-999999;
            ephilic=-999999;
            lnKow=-999999;
            fitness=0;
            p_sel=0;
            P=100000;
            T=298.15;
            
            filenum=-10;
			nosubmitscribd=-10;
			sameasancestor=0;
			tosubmit=1;
			frac=1.000000;
			tcal=0;

			data=NULL;
			//data.set_up();
			chg=0;
			multiplicity=1;
			molesmi=="";
			ionic="";
			parent=-1;
			if_circle=0;
			dist=NULL;
			connect=NULL;
			order=NULL;
			atm=NULL;

			natom=0;
			Pindex.resize(0);
			Cindex.resize(0);
			Mindex.resize(0);
			Cyindex.resize(0);
			Rindex.resize(0);
			protect.resize(0);
			Bindex.resize(0);
			atomsmi.resize(0);
			ctsisomer.resize(0);
			//for (int i=0;i<1024;i++) {
			//	Bindex[i].resize(0);
			//	atomsmi[i].resize(0);
			//}

		}

		~MOLECULE() {
			clean();
			empty();
		}

		int print(); 
		int print(ofstream &ouf);
		int replace(MOLECULE &);
		int ring(int,int);
		void mds2smi(bool ct_on=1);
		int subtract(int,int);
		int add(int,int,int);
		int add();
		void reset();
		void clear();
		void empty();
		void clean();
		void wipe();
		void outmds(ostream &outs=cout);
		int exchange(int,int,int,int);
		void read(string);
		void smi2gjf(); //OBMol &mol
		void chir_and_stereo(OBMol &mol);
		void init();
		void input(bool ct_on=1);
		void check_bnd(); //OBMol &mol
		void smi2cod();
		void report();
		int rd_nps(); //OBMol &mol
		int crossover(MOLECULE &, int, int);
		int combine(MOLECULE &, int, int, int);												// give a exact point to combine
		//void recode(int);
        int prct(int,int);
        int prct();
		double prob();
		int del_unpaired_ring_no();
		int decyc_small_ring(int);
		int canonicalize_SMILES();
		int rechg();
		void chk_cistrans(int sposi,int lposi);
		//int change_stereo();
		//int change_cistrans();
		int change_cistrans(int pos,int he);
		//int change_stereo(int pos);
		vector<int> Pindex; 															// parent index
		vector<int> Cindex; 															// child index
		vector<int> Mindex; 															// molecule index
		vector<int> Rindex;
		vector<int> Cyindex;
		vector<int> protect;
		//vector<int> Bindex[1024]; 														// to measure bond to 0, do not output.
		//vector<int> atomsmi[1024];
		//vector< vector<string> > ctsisomer;
		vector< vector<string> > ctsisomer;
		vector< vector<int> > Bindex;
		vector< vector<int> > atomsmi;
		int chg;
		int multiplicity;
		//POOL data;		
		POOL *data;
		string smiles;
        string molesmi;                                 								// smiles representation
		string ionic;
		int natom;
		double **dist;
		int **connect;
		int **order;
		DEATOM *atm;
		int parent;
		int if_circle;
			
	private:
};



#endif
