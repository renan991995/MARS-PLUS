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
#include "ELEMENTS.h"
#include "PARAMETER.h"
#include "MOLECULE.h"
#include "IL.h"
#include "UTILITY.h"
#include "CASES_NEU.h"
using namespace std;


unsigned int exhaustive_combination1(IL &A, IL &B);
unsigned int exhaustive_crossover1(IL &A, IL &B);
unsigned int exhaustive_subtraction1(IL &A);
unsigned int exhaustive_addition1(IL &A);
unsigned int exhaustive_change_ele1(IL &A);
unsigned int exhaustive_change_bnd1(IL &A);
unsigned int exhaustive_insertion1(IL &A);
unsigned int exhaustive_cyclization1(IL &A);
unsigned int exhaustive_change_cistrans1(IL &A);
unsigned int exhaustive_change_chirality1(IL &A);
unsigned int exhaustive_decyclization1(IL &A);

#endif
