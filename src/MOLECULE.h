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
#include <openbabel/bond.h>
#include <openbabel/bondtyper.h>
#include <openbabel/ring.h>
#include <openbabel/atom.h>
#include <openbabel/data.h>
#include <openbabel/data_utilities.h>
#include <openbabel/obiter.h>
//#include <openbabel/chiral.h>
#include <openbabel/elements.h>
#include <openbabel/distgeom.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/stereo/tetranonplanar.h>
#include <openbabel/stereo/tetraplanar.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/builder.h>
#include <openbabel/typer.h>
#include <openbabel/tokenst.h>
#include <openbabel/oberror.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mcdlutil.h>
#include <openbabel/forcefield.h>
#include <openbabel/groupcontrib.h>
//#include <openbabel/stereo/bindings.h>

#include <sys/types.h>
#include <dirent.h>
#include <cerrno>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>
#include <ctime>
#include <cstdio> 
#include <cctype>
#include <sstream>
#include <iterator>
//#include <random>
#include <memory>
#include <type_traits>
#include <new>
#include "ELEMENTS.h"
#include "PARAMETER.h"
using namespace std;
using namespace OpenBabel;

class MOLECULE
{
	public:
		MOLECULE() {
			data=NULL;
			chg=0;
			multiplicity=1;
			molesmi=="";
			parent=-1;
			if_circle=0;
			connect=NULL;
			atm=NULL;

			nsubcomp=1;
			comp_id=1;
			frac=1;

			natom=0;
			Pindex.resize(0);
			Cindex.resize(0);
			Mindex.resize(0);
			Cyindex.resize(0,vector<unsigned int>(0));
			Cybnd.resize(0);
			Rindex.resize(0);
			protect.resize(0);
			Bindex.resize(0,vector<unsigned int>(0));
			atomsmi.resize(0,vector<int>(0));
			ctsisomer.resize(0);
			chi.resize(0);
		}

		~MOLECULE() {
			clean();
			empty();
		}

		unsigned int cyclization(unsigned int,unsigned int,unsigned int);
        unsigned int subtraction(unsigned int,unsigned int bndfrm=0,bool cistrans=0);
        unsigned int insertion(unsigned int,unsigned int,unsigned int,unsigned int,bool cistrans=0);
        unsigned int addition(unsigned int,unsigned int,unsigned int,bool cistrans=0);
        unsigned int change_ele(unsigned int,unsigned int,unsigned int,unsigned int,bool cistrans=0);
        unsigned int change_bnd(unsigned int,unsigned int,unsigned int,unsigned int,bool cistrans=0);
        unsigned int crossover(MOLECULE &,unsigned int,unsigned int,bool cistrans=0);
        unsigned int combination(MOLECULE &,unsigned int,unsigned int,unsigned int);
		unsigned int decyclization(unsigned int);
        unsigned int change_cistrans(unsigned int,unsigned int);
        unsigned int change_chirality(unsigned int,unsigned int spec=2);

        unsigned int del_unpaired_ring_no();
        unsigned int decyc_small_ring(unsigned int);
		unsigned int chk_imine_ct(unsigned int,unsigned int,bool cistrans=0);
		unsigned int chk_cistrans(unsigned int,unsigned int,bool cistrans=0);
		unsigned int chk_chirality(unsigned int sposi,unsigned int lposi,unsigned int chirality=2);
        unsigned int prct(unsigned int,unsigned int);
        unsigned int prct();

        unsigned int printmds();
        unsigned int printmds(ofstream &ouf);
        unsigned int replace(MOLECULE &);
		unsigned int mds2smi(unsigned int chirality=2,bool cistrans=0,bool ct_on=0); //ct_on=1
		unsigned int reset();
		unsigned int clear();
		unsigned int empty();
		unsigned int clean();
		unsigned int wipe();
		unsigned int read(string);
		unsigned int input(bool ct_on=1);
		unsigned int smi2mds_OBabel();
		unsigned int canonicalize_SMILES();
		unsigned int rechg();
		unsigned int readmds(string);
		//unsigned int fcngroup2mds(unsigned int chirality=2,bool cistrans=0,bool ct_on=0);

		vector<unsigned int> Pindex; 															// parent index
		vector<unsigned int> Cindex; 															// child index
		vector<unsigned int> Mindex; 															// molecule index
		vector<unsigned int> Rindex;
		vector< vector<unsigned int> > Cyindex;
		vector<unsigned int> protect;
		vector<unsigned int> Cybnd;
		vector< vector<string> > ctsisomer;
		vector< vector<unsigned int> > Bindex;
		vector< vector<int> > atomsmi;
		vector<unsigned int> chi;

		int chg;
		int multiplicity;

		POOL *data;
		string smiles;
        string molesmi;                                 								// smiles representation
		unsigned int natom;
		unsigned int **connect;
		DEATOM *atm;
		int parent;
		unsigned int if_circle;
		unsigned int nsubcomp;
		unsigned int comp_id;
		double frac;
			
	private:
};



#endif
