#ifndef UTILITY_H
#define UTILITY_H
// #include <openbabel/obconversion.h>
// #include <openbabel/mol.h>
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
#include "ATOM.h"
#include "MATRIX.h"
#include "OPT.h"
#include "PARAMETER.h"
#include "MOLECULE.h"
#include "IL.h"
using namespace std;

void mark(ostream &);
void rd_para(char **&argv);
void rd_ILs(IL *&gsion,POOL *pol);
void rd_mols(MOLECULE *&gs,POOL *pol);
void canonicalize_mol(MOLECULE &A);
void canonicalize_IL(IL &A);
string getExt(string pathname);
int getdir(string dir, vector<string> &filenames, string fileextention);
int mk_datlist();

#endif




