#ifndef CASES_IL_INDEPENDENT_H
#define CASES_IL_INDEPENDENT_H
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
#include "CASES_NEU.h"
using namespace std;


void exhaustive_combination1(IL &A, IL &B);
void exhaustive_crossover1(IL &A, IL &B);
void exhaustive_subtraction1(IL &A);
void exhaustive_addition1(IL &A);
void exhaustive_exchange1(IL &A);
void exhaustive_cyclization1(IL &A);
void exhaustive_changect1(IL &A);

#endif
