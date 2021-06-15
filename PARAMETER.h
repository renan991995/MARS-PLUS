#ifndef PARAMETER_H
#define PARAMETER_H
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
using namespace std;


class PARAMETER {
	public :
		PARAMETER () {
		    system("pwd > pwddir.txt");
		    ifstream PWD("pwddir.txt");
		    pwd="";
		    PWD >> pwd >> ws;
		    PWD.close();
			
			smidir = pwd+"/";
			out_route = pwd+"/Summary.txt";
			
			gssize=0;
			ion=0;		
			protect=0;
			enumeration=0;
			operation="All";
			obabel="";
			round=1;
			epoch=1;
			outdir=smidir+"mds/";
			ifwritemds=0;
			glbouf.close();
			glbouf1.close();
			stat.close();
		}
		~PARAMETER() {
			glbouf.close();
			glbouf1.close();
			stat.close();
		}
		int gssize;
		string guess;
		string out_route;
		string operation;
		string obabel;
		bool ion;
		bool protect;
		bool enumeration;
		bool ifwritemds;
		int round;
		int epoch;
		string smidir;
		string pwd;	
		string outdir;
		ofstream glbouf;
		ofstream glbouf1;
		ofstream stat;
};


#endif

