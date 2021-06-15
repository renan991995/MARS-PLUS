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
void exhaustive_change_bnd(MOLECULE &A);
void exhaustive_change_ele(MOLECULE &A);
void exhaustive_insertion(MOLECULE &A);
void exhaustive_cyclization(MOLECULE &A);
void exhaustive_decyclization(MOLECULE &A);
void exhaustive_change_cistrans(MOLECULE &A);
void exhaustive_change_chirality(MOLECULE &A);

void CH4_to_Bz_path1(MOLECULE &A); //MOLECULE &A
void CH4_to_Bz_path2(MOLECULE &A); //MOLECULE &A
void CH4_to_aspirin(MOLECULE &A);
void CH4_to_alpha_carotene_path1(MOLECULE &A);
void CH4_to_alpha_carotene_path2(MOLECULE &A);
void CH4_to_alpha_carotene_path3(MOLECULE &A);

#endif




