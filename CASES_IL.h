#ifndef CASES_IL_H
#define CASES_IL_H
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


void exhaustive_combination(IL &A, IL &B);
void exhaustive_crossover(IL &A, IL &B);
void exhaustive_subtraction(IL &A);
void exhaustive_addition(IL &A);
void exhaustive_exchange(IL &A);
void exhaustive_cyclization(IL &A);
void swit(IL &A,IL &B);

#endif




