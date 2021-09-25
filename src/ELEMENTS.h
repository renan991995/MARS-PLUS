#ifndef ELEMENTS_h
#define ELEMENTS_h
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/data.h>
#include <openbabel/data_utilities.h>
#include <openbabel/obiter.h>
//#include <openbabel/chiral.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/stereo/tetranonplanar.h>
#include <openbabel/stereo/tetraplanar.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/builder.h>
#include <openbabel/typer.h>
#include <openbabel/tokenst.h>
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
#include "PARAMETER.h"
using namespace std;

class ATOM
{
	public:
		ATOM() {
			order.resize(0);
			name=atm="";
			index=2;
			chg=0;
			bd[0]=bd[1]=bd[2]=0;
			norder=nbond=4;
			id=0;
			chigenic=ctgenic=0;
			//availability=0;
			probability=0;
		}
		~ATOM() {
			order.clear();
		}
		unsigned int id;
		string name;
		string atm;
		vector<unsigned int> order;
		unsigned int bd[3];
		unsigned int index;
		//bool availability;
		double probability;
		unsigned int norder;
		bool chigenic;
		bool ctgenic;
		int chg;
		unsigned int nbond; // how many char for name
};

class DEATOM {
	public :
		DEATOM() {
			name="";
			chg=0;
			nH=0;
			chirality=0;
			cistrans[0]=cistrans[1]="";
		}
		~DEATOM() {}
		string name;
		unsigned int nH;
		unsigned int chirality;
		string cistrans[2];
		int chg; 
}; 

class POOL
{
	public:
		POOL(){
			a.resize(0); 
			num=0;
		}
		~POOL(){ 
            if(a.size()) {
                for (unsigned int r1=0;r1<num;r1++) {
                    if(a.at(r1).order.size()) {
                        a.at(r1).order.resize(0);
                    }
                }
                a.resize(0);
            }

		}
		unsigned int num;
		vector<ATOM> a;		
		unsigned int set_up();
		unsigned int read_in();

	private:
};
#endif
