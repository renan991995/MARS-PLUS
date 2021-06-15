#ifndef UTILITY_H
#define UTILITY_H
//#include <openbabel/obconversion.h>
//#include <openbabel/mol.h>
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
#include <cstring>
#include <sstream>
#include <bits/stdc++.h>
#include <algorithm>
#include <random>
#include <iterator>

#include "ELEMENTS.h"
#include "MATRIX.h"
#include "OPT.h"
#include "PARAMETER.h"
#include "MOLECULE.h"
#include "IL.h"
using namespace std;
using namespace OpenBabel;

void mark(ostream &);
void rd_para(char **&argv);
void rd_ILs(IL *&gsion,POOL *pol);
void rd_mols(MOLECULE *&gs,POOL *pol);
void canonicalize_mol(MOLECULE &A);
void canonicalize_IL(IL &A);
string getExt(string pathname);
int getdir(string dir, vector<string> &filenames, string fileextention);
int mk_datlist();
int smi2mol(string SMI,OBMol &mol);
string mol2smi(OBMol &mol);
int SMI_Enumerator(string SMI);
int cal_nmol(string r_inputlist);
int cal_nIL(string r_inputlist);
string rd_1molsmi(string r_inputlist,int num);
void rd_1ILsmi(string r_inputlist,int numcat,int numan,string &catsmi,string &ansmi);
void clearlog();

#endif




