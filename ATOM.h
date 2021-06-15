#ifndef ATOM_h
#define ATOM_h
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
#include "MATRIX.h"
using namespace std;

class ATOM
{
	public:
		ATOM() {
			order.resize(0);
			//cybnd.resize(0);
			//cypos.resize(0);
			//order=NULL;
			name=atm="";
			index=2;
			nh=0;
			chg=0;
			probability=0;
			bd[0]=bd[1]=bd[2]=0;
			norder=nbond=4;
			rb=0;
			type=1;
			id=0;
			ang0=0;
			probability=0;
		}
		~ATOM() {
			order.clear();
			//cybnd.clear();
			//cypos.clear();
			//vector<int>().swap(order);
		}
		int id;
		string name;
		string atm;
		vector<int> order;
		//int *order;
		int nh;
		int bd[3];
		int type;
		//vector<int> cybnd;
		//vector< vector<int> > cypos;
		int index;
		double rb;
		double probability;
		double ang0;
		int norder;
		int chg;
		int nbond; // how many char for name
};

class DEATOM {
	public :
		DEATOM() {
			//x[0]=x[1]=x[2]=0;
			name="";
			r_bnd=0;
			nbnd=0;
			chg=0;
			nH=0;
			chirality=0;
			cistrans[0]=cistrans[1]="";
			//cistrans="";
		}
		~DEATOM() {}
		string name;
		double r_bnd;
		int nbnd;
		int nH;
		int chirality;
		//string cistrans;
		string cistrans[2];
		//double x[3];
		void find_r();
		int chg; //20200211
}; 

class POOL
{
	public:
		POOL(){
			a.resize(0); 
			//a=NULL;
			num=0;
		}
		~POOL(){ 
            //a.clear();
            if(a.size()) {
                for (int r1=0;r1<num;r1++) {
                    if(a.at(r1).order.size()) {
                        a.at(r1).order.resize(0);
                        //vector<int>().swap(a.at(r1).order);
                    }
                }
                a.resize(0);
                //vector<ATOM>().swap(a);
            }

            //vector<ATOM>().swap(a);
		}
		int num;
		//ATOM *a; //  types of atoms
		vector<ATOM> a;		
		void set_up();

	private:
};
#endif
