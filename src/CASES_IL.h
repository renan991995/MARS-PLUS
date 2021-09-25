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
#include "ELEMENTS.h"
#include "PARAMETER.h"
#include "MOLECULE.h"
#include "IL.h"
#include "UTILITY.h"
#include "CASES_NEU.h"
using namespace std;


unsigned int exhaustive_combination(IL &A, IL &B);
unsigned int exhaustive_crossover(IL &A, IL &B);
unsigned int exhaustive_subtraction(IL &A);
unsigned int exhaustive_addition(IL &A);
unsigned int exhaustive_change_ele(IL &A);
unsigned int exhaustive_cyclization(IL &A);
unsigned int comp_swap(IL &A,IL &B);

#endif




