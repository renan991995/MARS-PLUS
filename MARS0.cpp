#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "ATOM.h"
#include "MATRIX.h"
#include "MOLECULE.h"
#include "OPT.h"

using namespace std;
void mark(ostream &);
int main(int argc, char **argv) {
  mark(cout);

  int i, j, k, dd3 = 0, enc = 0;
  MOLECULE *gs;
  gs = new MOLECULE[2];
  string smiles;
  ofstream out;

  cout << "Creating caffeine, from methane, C" << endl;
  out.close();
  gs[0].smiles = "methane";
  gs[0].input();
  out.open("caffeine.txt");
  mark(out);
  out << " 1 ";
  out << gs[0].molesmi << endl;

  out << " 2 ";
  gs[0].add(0, 7);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 3 ";
  gs[0].add(1, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 4 ";
  gs[0].add(2, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 5 ";
  gs[0].add(3, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 6 ";
  gs[0].add(4, 7);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 7 ";
  gs[0].ring(5, 0);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 8 ";
  gs[0].add(3, 7);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 9 ";
  gs[0].add(6, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "10 ";
  gs[0].add(7, 8);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "11 ";
  gs[0].ring(8, 4);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "12 ";
  gs[0].add(0, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "13 ";
  gs[0].exchange(9, 6, 2, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "14 ";
  gs[0].add(2, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "15 ";
  gs[0].exchange(10, 6, 2, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "16 ";
  gs[0].add(1, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "17 ";
  gs[0].add(6, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "18 ";
  gs[0].add(5, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out.close();
  cout << "Job completed!" << endl;

  cout << "Creating TRIS, from methane, C" << endl;
  out.close();
  gs[0].smiles = "methane";
  gs[0].input();
  out.open("TRIS.txt");
  mark(out);
  out << " 1 ";
  out << gs[0].molesmi << endl;

  out << " 2 ";
  gs[0].add(0, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 3 ";
  gs[0].add(0, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 4 ";
  gs[0].add(0, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 5 ";
  gs[0].add(1, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 6 ";
  gs[0].add(2, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 7 ";
  gs[0].add(3, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << " 8 ";
  gs[0].add(0, 7);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out.close();
  cout << "Job completed!" << endl;

  cout << "Creating ATP, from methane, C" << endl;
  out.close();
  gs[0].smiles = "methane";
  gs[0].input();
  out.open("ATP.txt");
  mark(out);
  out << gs[0].molesmi << endl;

  gs[0].add(0, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(1, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(2, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(3, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].ring(4, 0);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(0, 7);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(5, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(6, 8);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(7, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(8, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].ring(9, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out << "@@" << endl;
  gs[0].add(8, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(10, 8);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(11, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(12, 8);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].ring(13, 9);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(10, 7);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(1, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(2, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(3, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(17, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(18, 32);  // 18?
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(19, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(20, 32);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(21, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(22, 32);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(23, 5);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(19, 10);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(21, 10);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(23, 10);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(19, 6);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(21, 6);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(23, 6);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out.close();
  cout << "Job completed!" << endl;

  cout << "Creating phenanthrene, from methane, C" << endl;
  out.close();
  gs[0].smiles = "methane";
  gs[0].input();
  out.open("phenanthrene.txt");
  mark(out);
  out << gs[0].molesmi << endl;

  gs[0].add(0, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].exchange(1, 2, 2, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(1, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(2, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(3, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(4, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].ring(5, 0);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(0, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(6, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(7, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(8, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].ring(9, 1);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(6, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(10, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(11, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].add(12, 2);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  gs[0].ring(13, 7);
  gs[0].mds2smi();
  out << gs[0].molesmi << endl;

  out.close();
  cout << "Job completed!" << endl;
  // cout << "Results can be found in corresponding .txt files.
  // (combination.txt, "
  //         "crossover.txt, addition.txt, subtraction.txt, exchange, ring.txt,
  //         " "and benzene.txt)"
  //      << endl;
  delete[] gs;
  return 1;
}

void mark(ostream &out) {
  out << "*********************************************************************"
         "******************"
      << endl
      << "                   MARS: Molecular Assembling and Representation "
         "Suite                 "
      << endl
      << "                                      Hsuan Hao Hsu                  "
         "                  "
      << endl
      << "                                    Chen Hsuan Huang                 "
         "                  "
      << endl
      << "                            Shiang Tai Lin (stlin@ntu.edu.tw)        "
         "                  "
      << endl
      << "                      Computational Molecular Engineering Laboratory "
         "                  "
      << endl
      << "      Department of Chemical Engineering, National Taiwan "
         "University, Taipei, Taiwan   "
      << endl
      << "                                   Copyright (c) 2019                "
         "                  "
      << endl
      << "                                   All rights Reserved               "
         "                  "
      << endl
      << "*********************************************************************"
         "******************"
      << endl;
  return;
}
