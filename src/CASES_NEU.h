#ifndef CASES_NEU_H
#define CASES_NEU_H
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
using namespace std;


unsigned int exhaustive_combination(MOLECULE &A, MOLECULE &B);
unsigned int exhaustive_crossover(MOLECULE &A, MOLECULE &B);
unsigned int exhaustive_subtraction(MOLECULE &A);
unsigned int exhaustive_addition(MOLECULE &A);
unsigned int exhaustive_change_bnd(MOLECULE &A);
unsigned int exhaustive_change_ele(MOLECULE &A);
unsigned int exhaustive_insertion(MOLECULE &A);
unsigned int exhaustive_cyclization(MOLECULE &A);
unsigned int exhaustive_decyclization(MOLECULE &A);
unsigned int exhaustive_change_cistrans(MOLECULE &A);
unsigned int exhaustive_change_chirality(MOLECULE &A);

unsigned int CH4_to_Bz_path1(MOLECULE &A); //MOLECULE &A
unsigned int CH4_to_Bz_path2(MOLECULE &A); //MOLECULE &A
unsigned int CH4_to_aspirin(MOLECULE &A);
unsigned int CH4_to_alpha_carotene_path1(MOLECULE &A);
unsigned int CH4_to_alpha_carotene_path2(MOLECULE &A);
unsigned int CH4_to_alpha_carotene_path3(MOLECULE &A);

unsigned int Tamiflu_Corey(MOLECULE &A,MOLECULE &B);

#endif




