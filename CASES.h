#ifndef CASES_H
#define CASES_H
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
#include "MOLECULE.h"
#include "IL.h"
#include "UTILITY.h"
using namespace std;


void exhaustive_combination(MOLECULE &A, MOLECULE &B);
void exhaustive_crossover(MOLECULE &A, MOLECULE &B);
void exhaustive_subtraction(MOLECULE &A);
void exhaustive_addition(MOLECULE &A);
void exhaustive_exchange(MOLECULE &A);
void exhaustive_cyclization(MOLECULE &A);
void CH4_to_Bz(MOLECULE &A);

#endif




