#include <string.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <typeinfo>
#include <vector>
#include "ATOM.h"
#include "MATRIX.h"
#include "MOLECULE.h"
#include "OPT.h"

using namespace std;
void mark(ostream &);

void TRIS() {
	mark(cout);

	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;
	cout << "Creating TRIS from methane, C";

	gs[0].smiles = "methane";
	gs[0].input();
	out.open("_TRIS.txt");
	mark(out);
	out << " 1 ";
	out << gs[0].molesmi << endl;



	out << " 2 ";
	gs[0].add(0, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;
	gs[0].outmds(out);


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

	// out << " 9 ";
	// gs[0].add(7, 36);
	// gs[0].mds2smi();
	// out << gs[0].molesmi << endl;

	out.close();
	cout << "Job completed!" << endl;

	delete [] gs;
	gs=NULL;
}
void caffeine() {
	mark(cout);

	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;

	cout << "Creating caffeine, from methane, C" << endl;
	gs[0].smiles = "methane";
	gs[0].input();
	out.open("_caffeine.txt");
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
	gs[0].add(0, 2);  // 1 or 2 is ok
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

	delete [] gs;
	gs=NULL;
}

void ATP() {
	mark(cout);

	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;

	cout << "Creating ATP, from methane, C" << endl;
	gs[0].smiles = "methane";
	gs[0].input();
	out.open("_ATP.txt");
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

	gs[0].add(9, 8);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(10, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(11, 8);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(12, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].ring(13, 8);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;
	//  gs[0].outmds(out);

	gs[0].add(13, 7);
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

	delete[] gs;
	gs=NULL;
}
void phenanthrene() {
	mark(cout);

	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;

	cout << "Creating phenanthrene, from methane, C" << endl;
	gs[0].smiles = "methane";
	gs[0].input();
	out.open("_phenanthrene.txt");
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

	delete [] gs;
	gs=NULL;
}
void morphine() {
	mark(cout);

	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;

	cout << "Creating morphine, from methane, C" << endl;
	gs[0].smiles = "methane";
	gs[0].input();
	out.open("_morphine.txt");
	mark(out);
	out << gs[0].molesmi << endl;

	gs[0].add(0, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;
	gs[0].outmds(out);

	gs[0].add(1, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(2, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(3, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(4, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	out << "r1" << endl;
	gs[0].ring(5, 0);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(5, 5);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(0, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(7, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(7, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	out << "r2" << endl;
	gs[0].ring(9, 6);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(9, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(10, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(11, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	out << "r3" << endl;
	gs[0].ring(12, 8);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(8, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(13, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	out << "r4" << endl;

	gs[0].ring(14, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(14, 7);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;
	gs[0].outmds(out);
	// out.close();

	gs[0].add(0, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;
	gs[0].outmds(out);

	out << "r5" << endl;
	gs[0].ring(16, 15);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(15, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(4, 5);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;
	gs[0].outmds(out);

	cout << "Job completed!" << endl;

	delete [] gs;
	gs=NULL;
}
void aspirin() {
	mark(cout);

	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;
	cout << "Creating aspirin from methane, C" << endl;

	gs[0].smiles = "methane";
	gs[0].input();
	out.open("_aspirin.txt");
	mark(out);
	out << gs[0].molesmi << endl;

	gs[0].add(0, 2);
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

	gs[0].add(0, 5);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(6, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(7, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(7, 6);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(1, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(10, 6);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(10, 5);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	cout << "Job completed!" << endl;

	delete [] gs;
	gs=NULL;
}
void propenenitrile() {
	mark(cout);

	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;
	cout << "Creating propenenitrile from methane, C" << endl;

	gs[0].smiles = "methane";
	gs[0].input();
	out.open("_propenenitrile.txt");
	mark(out);
	out << gs[0].molesmi << endl;

	gs[0].add(0, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(1, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].exchange(2, 2, 2, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(0, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].exchange(3, 9, 3, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	out.close();
	cout << "Job completed!" << endl;

	delete[] gs;
	gs=NULL;
}
void b_ben() {
	mark(cout);
	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;
	cout << "Creating ben from benzene" << endl;
	out.close();

	gs[0].smiles = "benzene";
	gs[0].input();
	out.open("_ben.txt");
	mark(out);
	out << gs[0].molesmi << endl;

	gs[0].add(0, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;
	gs[0].outmds(out);

	gs[0].add(6, 5);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(7, 2);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].ring(8, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(8, 6);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	gs[0].add(6, 6);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;

	out.close();
	cout << "Job completed!" << endl;

	delete [] gs;
	gs=NULL;
}

void com(int num) {
	mark(cout);

	//  int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;
	cout << "Creating " << num << "C from methane, C" << endl;

	gs[0].smiles = "methane";
	gs[0].input();

	string filename;
	stringstream numstr;
	numstr << num;
	filename = "_com" + numstr.str() + ".txt";
	out.open(filename.c_str());
	mark(out);
	// out << gs[0].molesmi << endl;
	vector<MOLECULE> chlist;
	chlist.push_back(gs[0]);
	chlist.push_back(gs[0]);
	//  gs[0].outmds(out);

	/*
	   int i;
	   vector<string> chlist;
	// vector<MOLECULE> chlist;
	for (i = 0; i < num; i++) {
	if (gs[0].Cindex.size() < num) {
	out << "~ " << gs[0].Cindex.size() << endl;
	gs[0].add(i, 1);
	gs[0].mds2smi();
	// cout << typeid(gs[0].molesmi).name() << endl;
	chlist.push_back(gs[0].molesmi);
	//   out << chlist[i] << endl;

	// out << gs[0].molesmi << endl;
	}
	}
	out << gs[0].molesmi << endl;

	for (int i = 0; i < chlist.size(); i++) {
	cout << chlist.at(i) << endl;
	// cout << chlist.size() << endl;


	}
	*/
	int i;

	for (i = 0; i < num; i++) {
		if (gs[0].Cindex.size() < num) {
			// out << "~ " << gs[0].Cindex.size() << endl;
			gs[0].add(i, 1);
			gs[0].mds2smi();
			// cout << typeid(gs[0].molesmi).name() << endl;
			chlist.push_back(gs[0]);
			// out << chlist.at(i).molesmi << endl;
			// out << gs[0].molesmi << endl;
		}
	}
	out << "main chain" << endl;
	out << gs[0].molesmi << endl;

	for (int i = 0; i < chlist.size(); i++) {
		cout << chlist.at(i).molesmi << endl;
	}
	vector<MOLECULE> copy_chlist;
	copy_chlist.assign(chlist.begin(), chlist.end());
	// make_ring(chlist.at(num));
	// cout<<typeid(chlist.at(num)).name()<<endl;

	out << "ring" << endl;
	int m, n;
	for (m = 0; m < num / 2; m++) {
		for (n = m; n < num - 2; n++) {
			chlist[num].ring(m, n + 2);
			chlist[num].mds2smi();
			out << chlist[num].molesmi << endl;
			chlist.at(num) = copy_chlist.at(num);
		}
	}

	/*
	   switch (num) {
	   case 3:
	   chlist[1].ring(0, 2);
	   chlist[1].mds2smi();
	   out << chlist[1].molesmi << endl;
	   break;
	   case 4:
	   chlist[2].ring(0, 2);
	   chlist[2].mds2smi();
	   out << chlist[2].molesmi << endl;

	   chlist[2] = copy_chlist[2];
	   out << "//" << chlist[2].molesmi << endl;
	   chlist[2].ring(0, 3);
	   chlist[2].mds2smi();
	   out << chlist[2].molesmi << endl;
	   }
	   */
	/*
	   void make_ring(Molecule chain) {
	   MOLECULE copy_chain = chain;
	   int i, j;
	   for (i = 0; i < n / 2; i++) {
	   for (j = 0; j < i; j++) {
	   chain.ring(j, j + 2);
		// chain.mds2smi();
		// out << chain.molesmi << endl;

		chain = copy_chain;
	}
	}
	}
	*/

	out << "exchange, 1 pi" << endl;
	for (m = 1; m < num / 2 + 1; m++) {  // can't begin at 0
		chlist[num].exchange(m, 2, 2, 1);
		chlist[num].mds2smi();
		out << chlist[num].molesmi << endl;
		chlist.at(num) = copy_chlist.at(num);
	}

	vector<MOLECULE> conj_chlist;
	out << "conjugated pi" << endl;
	for (m = 1; m < num; m += 2) {
		chlist[num].exchange(m, 2, 2, 1);
		chlist[num].mds2smi();

		out << chlist[num].molesmi << endl;
	}
	conj_chlist.push_back(chlist.at(num));  // when have the most pi
	for (m = 1; m < num; m += 2) {
		chlist[num].exchange(m, 1, 1, -1);
		chlist[num].mds2smi();
		out << chlist[num].molesmi << endl;
	}

	if (num > 5) {
		out << "conjugated to benzoid" << endl;
		// cout << "/ " << conj_chlist[0].molesmi << endl;
		conj_chlist[0].ring(0, 5);
		conj_chlist[0].mds2smi();
		out << conj_chlist[0].molesmi << endl;
	}

	out << "exchange, three bond" << endl;
	for (m = 1; m < num / 2 + 1; m++) {  // can't begin at 0
		chlist[num].exchange(m, 3, 3, 2);
		chlist[num].mds2smi();
		out << chlist[num].molesmi << endl;
		chlist.at(num) = copy_chlist.at(num);
	}

/*
   out << "1 side" << endl;
   for (m = 1; m < num / 2 + 1; m++) {
   chlist[num - 1].add(m, 1);
   chlist[num - 1].mds2smi();
   out << chlist[num - 1].molesmi << endl;
   chlist.at(num - 1) = copy_chlist.at(num - 1);
   }
   */
/*
   out << "1 side" << endl;
   for (int p = 1; p < num / 2; p++) {  // position
   for (int k = 1; k <= p; k++) {  // times of doing (the length of side
   group) for (int j = 0; j < p; j++) {  // counter
// chlist[num - k].add((p + j), 1);   // wrong
if (j == 0)
chlist[num - k].add(p, 1);
else
chlist[num - k].add((num + j), 1);

chlist[num - k].mds2smi();

out << chlist[num - k].molesmi << endl;
}
chlist.at(num - k) = copy_chlist.at(num - k);
}
}
}
*/
// out << "1 side" << endl;
// for (int p = 1; p < num / 2; p++) {  // position
//   for (int j = 0; j < p; j++) {      // counter
//     // chlist[num - k].add((p + j), 1);   // wrong
//     if (j == 0)
//       chlist[num - p].add(p, 1);
//     else
//       chlist[num - p].add((num + j), 1);

//     chlist[num - p].mds2smi();

//     out << chlist[num - p].molesmi << endl;
//   }
//   chlist.at(num - p) = copy_chlist.at(num - p);
// }
}

void CH3COOH() {
	mark(cout);

	int i, j, k, dd3 = 0, enc = 0;
	MOLECULE *gs=NULL;
	gs = new MOLECULE[2];
	string smiles;
	ofstream out;
	cout << "Creating CH3COOH from methane, C";

	gs[0].smiles = "methane";
	gs[0].input();
	out.open("_CH3COOH.txt");
	//mark(out);

	out << gs[0].molesmi << endl;


	gs[0].add(0, 36);
	gs[0].mds2smi();
    out << gs[0].molesmi << endl;
	gs[0].add(0, 1);
    gs[0].mds2smi();
    out << gs[0].molesmi << endl;
    gs[0].add(0, 1);
    gs[0].mds2smi();
    out << gs[0].molesmi << endl;
    gs[0].add(0, 1);
    gs[0].mds2smi();
    out << gs[0].molesmi << endl;
    gs[0].add(2, 1);
    gs[0].mds2smi();
    out << gs[0].molesmi << endl;

	/*
    gs[0].add(2, 36);
    gs[0].mds2smi();
    out << gs[0].molesmi << endl;
	*/

	/*
	gs[0].add(0, 1);
	gs[0].mds2smi();
	out << gs[0].molesmi << endl;
    gs[0].add(1, 2);
    gs[0].mds2smi();
	out << gs[0].molesmi << endl;
    gs[0].add(1, 6);
    gs[0].mds2smi();
	out << gs[0].molesmi << endl;
    gs[0].add(1, 10);
    gs[0].mds2smi();
	out << gs[0].molesmi << endl;
	*/
	//  gs[0].outmds(out);


	out.close();
	cout << "Job completed!" << endl;

	delete[] gs;
	gs=NULL;
}

int main(int argc, char **argv) {
	int num=9;
	MOLECULE *mole=NULL;
	mole=new MOLECULE [num];
	caffeine();
	TRIS();
	phenanthrene();
	ATP();
	morphine();
	aspirin();
	propenenitrile();
	b_ben();
	CH3COOH();


	com(3);
	com(4);
	com(5);
	com(6);
	com(7);
	com(10);

	// cout << "Results can be found in corresponding .txt files.
	// (combination.txt, "
	//         "crossover.txt, addition.txt, subtraction.txt, exchange,
	//         ring.txt, " "and benzene.txt)"
	//      << endl;
	delete [] mole;
	mole=NULL;
	return 0;
}

void mark(ostream &out) {
	out << "*******************************************************************"
		"**"
		"******************"
		<< endl
		<< "                   MARS: Molecular Assembling and Representation "
		"Suite                 "
		<< endl
		<< "                                      Hsuan Hao Hsu                "
		"  "
		"                  "
		<< endl
		<< "                                    Chen Hsuan Huang               "
		"  "
		"                  "
		<< endl
		<< "                            Shiang Tai Lin (stlin@ntu.edu.tw)      "
		"  "
		"                  "
		<< endl
		<< "                      Computational Molecular Engineering "
		"Laboratory "
		"                  "
		<< endl
		<< "      Department of Chemical Engineering, National Taiwan "
		"University, Taipei, Taiwan   "
		<< endl
		<< "                                   Copyright (c) 2019              "
		"  "
		"                  "
		<< endl
		<< "                                   All rights Reserved             "
		"  "
		"                  "
		<< endl
		<< "*******************************************************************"
		"**"
		"******************"
		<< endl;
	return;
}
