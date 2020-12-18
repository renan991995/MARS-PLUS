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
			operation="";
			obabel="";
		}
		~PARAMETER() {}
		int gssize;
		string guess;
		string out_route;
		string operation;
		string obabel;
		int ion;
		int protect;
		string smidir;
		string pwd;	
};


#endif

