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
#include "PARAMETER.h"
#include "MOLECULE.h"
#include "IL.h"
using namespace std;
using namespace OpenBabel;

unsigned int mark(ostream &);
unsigned int rd_para(char **&argv);
unsigned int output_para();
unsigned int rd_ILs(IL *&gsion,POOL *pol);
unsigned int rd_mols(MOLECULE *&gs,POOL *pol);
unsigned int canonicalize_mol(MOLECULE &A);
unsigned int canonicalize_IL(IL &A);
string getExt(string pathname);
unsigned int getdir(string dir, vector<string> &filenames, string fileextention);
unsigned int mk_datlist();
unsigned int smi2mol(string SMI,OBMol &mol);
string mol2smi(OBMol &mol);
unsigned int SMI_Enumerator(string SMI);
unsigned int cal_nmol(string r_inputlist);
unsigned int cal_nIL(string r_inputlist);
string rd_1molsmi(string r_inputlist,unsigned int num);
unsigned int rd_1ILsmi(string r_inputlist,unsigned int numcat,unsigned int numan,string &catsmi,string &ansmi);
unsigned int clearlog();

#endif




