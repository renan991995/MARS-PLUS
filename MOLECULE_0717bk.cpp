#include "MOLECULE.h"
#define p_circle 0.2
using namespace std;
using namespace OpenBabel;

extern PARAMETER para;

double MOLECULE::prob() {
	double p=double(rand())/RAND_MAX;
	return p;
}


int MOLECULE::replace(MOLECULE &mol) {
	if (mol.Cindex.size()==0) {
		cout<<"Error: It cannot be replaced by an empty molecule. smiles = " << mol.smiles << " | molesmi = " << mol.molesmi << endl;
		return 0;
		//exit(1);
	}
	else {
		int i=0,j=0;
		empty();
		clear();
		clean(); //20200708

		for (i=0;i<mol.Cindex.size();i++) {
			Cindex.push_back(mol.Cindex.at(i));
			Pindex.push_back(mol.Pindex.at(i));
			Rindex.push_back(mol.Rindex.at(i));
			Cyindex.push_back(mol.Cyindex.at(i));
			Mindex.push_back(mol.Mindex.at(i));
			if (para.protect) protect.push_back(mol.protect.at(i));
		}
		if_circle=mol.if_circle;

		smiles=mol.smiles;
		molesmi=mol.molesmi;
		nsubcomp=1;
		comp_id=mol.comp_id;
		chg=mol.chg;
		multiplicity=mol.multiplicity;

		reset();
		mds2smi();
        smiles=mol.smiles;
        molesmi=mol.molesmi;
	}
	return 1;
}

int MOLECULE::crossover(MOLECULE &aaa,int pp, int jj) {
	if (para.protect) {
		if (protect.at(pp)) return 0; //20200625
		if (aaa.protect.at(jj)) return 0; //20200625
	}
	int n;
	int m;
	int k=0;
	int a;
	int b;
	int y;
	int x;
	int c;
	int ref;
	int delta;
	int pos_i;
	int pos_j;
	vector<int> t,f;
	vector<int> ci_ref, mi_ref,pi_ref,ri_ref,cyi,pri;
	vector<int> cj_ref, mj_ref,pj_ref,rj_ref,cyj,prj;

	vector<bool> ifchg(4,0);
	if (para.ion) {
	for (int k1=0;k1<Cindex.size();k1++) {
		if (Cindex.at(k1)<=Cindex.at(pp) && data.a[Mindex.at(k1)].chg) {
			ifchg.at(0)=1;
			break;
		}
	}
	for (int k1=0;k1<Cindex.size();k1++) {
		if (Cindex.at(k1)>=Cindex.at(pp) && data.a[Mindex.at(k1)].chg) {
			ifchg.at(1)=1;
			break;
		}
	}
	for (int k1=0;k1<aaa.Cindex.size();k1++) {
		if (aaa.Cindex.at(k1)<=aaa.Cindex.at(jj) && data.a[aaa.Mindex.at(k1)].chg) {
			ifchg.at(2)=1;
			break;
		}
	}
	for (int k1=0;k1<aaa.Cindex.size();k1++) {
		if (aaa.Cindex.at(k1)>=aaa.Cindex.at(jj) && data.a[aaa.Mindex.at(k1)].chg) {
			ifchg.at(3)=1;
			break;
		}
	}
	bool no1=0;
	//cout << "cros_chg " << ifchg.at(0) << " " << ifchg.at(1) << " " << ifchg.at(2) << " " << ifchg.at(3) << endl;
	if (0) {
		if (ifchg.at(0)==1 && ifchg.at(2)==1 && ifchg.at(1)==0 && ifchg.at(3)==0) no1=1;
		if (ifchg.at(1)==1 && ifchg.at(3)==1 && ifchg.at(0)==0 && ifchg.at(2)==0) no1=1;
	}
	if (0) {
		if (ifchg.at(0)==1 && ifchg.at(3)==1 && ifchg.at(1)==0 && ifchg.at(2)==0) no1=1;
		if (ifchg.at(1)==1 && ifchg.at(2)==1 && ifchg.at(0)==0 && ifchg.at(3)==0) no1=1;
	}
	if (1) {
		if (ifchg.at(0)==1 && ifchg.at(3)==1) no1=1;
		if (ifchg.at(1)==0 && ifchg.at(2)==0) no1=1;
		if (ifchg.at(1)==1 && ifchg.at(2)==1) no1=1;
		if (ifchg.at(0)==0 && ifchg.at(3)==0) no1=1;
	}
	if (no1) {
		/*
		cout << "pp: " << pp << " | " 
			<< "jj: " << jj << " | "
			<< "cros_chg: fail cros " << ifchg.at(0) << " " << ifchg.at(1) << " " << ifchg.at(2) << " " << ifchg.at(3) << endl;
		*/
		return 0;
	}
	}

    vector<bool> ifcyc(4,0);
    if (para.protect && if_circle) {
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Cindex.at(k1)<=Cindex.at(pp) && Cyindex.at(k1)) {
                ifcyc.at(0)=1;
                break;
            }
        }
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Cindex.at(k1)>=Cindex.at(pp) && Cyindex.at(k1)) {
                ifcyc.at(1)=1;
                break;
            }
        }
        for (int k1=0;k1<aaa.Cindex.size();k1++) {
            if (aaa.Cindex.at(k1)<=aaa.Cindex.at(jj) && aaa.Cyindex.at(k1)) {
                ifcyc.at(2)=1;
                break;
            }
        }
        for (int k1=0;k1<aaa.Cindex.size();k1++) {
            if (aaa.Cindex.at(k1)>=aaa.Cindex.at(jj) && aaa.Cyindex.at(k1)) {
                ifcyc.at(3)=1;
                break;
            }
        }
    }
	bool no2=0;
    if (1) {
        if (ifcyc.at(0)==1 && ifcyc.at(1)==1) no2=1;
        if (ifcyc.at(2)==1 && ifcyc.at(3)==1) no2=1;
    }
    if (no2) {
		/*
        cout << "pp: " << pp << " | "
            << "jj: " << jj << " | "
            << "cros_cyc: fail cros " << ifcyc.at(0) << " " << ifcyc.at(1) << " " << ifcyc.at(2) << " " << ifcyc.at(3) << endl;
		*/
        return 0;
    }

	//double time1=time(NULL);

	k=0;
	if (aaa.Rindex.at(jj)==Rindex.at(pp)) {
		k=1;
		m=jj;
		n=pp;
		b=aaa.Pindex.at(m);
		a=Pindex.at(n);
		y=aaa.Cindex.at(m);
		x=Cindex.at(n);
		ri_ref.push_back(Rindex.at(n));
		rj_ref.push_back(aaa.Rindex.at(m));
		ci_ref.push_back(x);
		cj_ref.push_back(y);
		pi_ref.push_back(a);
		pj_ref.push_back(b);
		mi_ref.push_back(Mindex.at(n));
		cyi.push_back(Cyindex.at(n));
		if (para.protect) pri.push_back(protect.at(n));
		cyj.push_back(aaa.Cyindex.at(m));
		if (para.protect) prj.push_back(aaa.protect.at(m));
		mj_ref.push_back(aaa.Mindex.at(m));
	}

	if (k==0) return 0;
	if (k==1) {
		for (n=1;n<Cindex.size();n++) {
			if (Cindex.at(n) > ci_ref.at(0)) {
				ref = ci_ref.size();
				for (m=0;m<ref;m++) {
					if (Pindex.at(n) == ci_ref.at(m)) {
						ci_ref.push_back(Cindex.at(n));
						mi_ref.push_back(Mindex.at(n));
						pi_ref.push_back(Pindex.at(n));
						ri_ref.push_back(Rindex.at(n));
						cyi.push_back(Cyindex.at(n));
						if (para.protect) pri.push_back(protect.at(n));
					}
				}
			}
		}
		for (m=1;m<aaa.Cindex.size();m++) {
			if (aaa.Cindex.at(m) > cj_ref.at(0) ) {
				ref = cj_ref.size();
				for (n=0;n<ref;n++) {
					if (aaa.Pindex.at(m) == cj_ref.at(n)) {
						cj_ref.push_back(aaa.Cindex.at(m));
						mj_ref.push_back(aaa.Mindex.at(m));
						pj_ref.push_back(aaa.Pindex.at(m));
						rj_ref.push_back(aaa.Rindex.at(m));
						cyj.push_back(aaa.Cyindex.at(m));
						if (para.protect) prj.push_back(aaa.protect.at(m));
					}
				}
			}
		}
		for (n=0;n<pri.size();n++) if (pri.at(n)) return 0;
		for (n=0;n<prj.size();n++) if (prj.at(n)) return 0;
		a=0;
		for (n=0;n<cyi.size();n++) {
			if (cyi.at(n)) {
				cyi.at(n)+=if_circle;
				a++;
			}
		}
		if (a%2 && if_circle) {
			for (n=0;n<cyi.size();n++) if (cyi.at(n)) cyi.at(n)=if_circle;
		}
		else if (a%2) {
			for (n=0;n<cyi.size();n++) cyi.at(n)=0;
		}
		a=0;
		for (n=0;n<cyj.size();n++) {
			if (cyj.at(n)) {
				cyj.at(n)+=aaa.if_circle;
				a++;
			}
		}
		if (a%2 && aaa.if_circle) {
			for (n=0;n<cyj.size();n++) if (cyj.at(n)) cyj.at(n)=aaa.if_circle;
		}
		if (a%2) {
			for (n=0;n<cyj.size();n++) cyj.at(n)=0;
		}
		x=ci_ref[0];
		y=cj_ref[0];
		a=pi_ref[0];
		b=pj_ref[0];
		for (n=0;n<ci_ref.size();n++) {
			for (m=0;m<Cindex.size();m++) {
				if (Cindex.at(m) == ci_ref.at(n)) {
					Mindex.erase(Mindex.begin()+m);
					Cindex.erase(Cindex.begin()+m);
					Pindex.erase(Pindex.begin()+m);
					Rindex.erase(Rindex.begin()+m);
					Cyindex.erase(Cyindex.begin()+m);
					if (para.protect) protect.erase(protect.begin()+m);
				}
			}
		}
		for (m=0;m<cj_ref.size();m++) {
			for (n=0;n<aaa.Cindex.size();n++) {
				if (aaa.Cindex.at(n) == cj_ref.at(m)) {
					aaa.Mindex.erase(aaa.Mindex.begin()+n);
					aaa.Cindex.erase(aaa.Cindex.begin()+n);
					aaa.Pindex.erase(aaa.Pindex.begin()+n);
					aaa.Rindex.erase(aaa.Rindex.begin()+n);
					aaa.Cyindex.erase(aaa.Cyindex.begin()+n);
					if (para.protect) aaa.protect.erase(aaa.protect.begin()+n);
				}
			}
		}
		t.clear();
		for (n=1;n<ci_ref.size();n++) {
			for (m=0;m<ci_ref.size();m++) {
				if (pi_ref.at(n) == ci_ref.at(m)) {
					t.push_back(m);
				}
			}
		}
		f.clear();
		for (n=1;n<cj_ref.size();n++) {
			for (m=0;m<cj_ref.size();m++) {
				if (pj_ref.at(n) == cj_ref.at(m)) {
					f.push_back(m);
				}
			}
		}
		pi_ref.at(0) = b;
		pj_ref.at(0) = a;
		ci_ref.at(0) = y;
		cj_ref.at(0) = x;

		for (n=1;n<ci_ref.size();n++) {
			ci_ref.at(n) = ci_ref.at(n-1) + 1;
			pi_ref.at(n) = ci_ref.at(t.at(n-1));
		}
		for (n=1;n<cj_ref.size();n++) {
			cj_ref.at(n) = cj_ref.at(n-1) + 1;
			pj_ref.at(n) = cj_ref.at(f.at(n-1));
		}

		t.clear();
		f.clear();
		for (n=1;n<Cindex.size();n++) {
			for (m=0;m<Cindex.size();m++) {
				if (Pindex.at(n) == Cindex.at(m)) {
					t.push_back(m);
				}
			}
		}
		ref=-1;
		for (n=1;n<Cindex.size();n++) {
			if (Cindex.at(n) != Cindex.at(n-1)+1) {
				ref=n;
				Cindex.at(ref) = Cindex.at(ref-1) + cj_ref.size() + 1;
				break;
			}
		}
		if (ref!=-1) {
			for (n=ref+1;n<Cindex.size();n++) {
				Cindex.at(n) = Cindex.at(n-1) + 1;
			}

			for (n=1;n<Cindex.size();n++) {
				Pindex.at(n) = Cindex.at(t.at(n-1));
			}
		}
		t.clear();
		f.clear();
		for (n=1;n<aaa.Cindex.size();n++) {
			for (m=0;m<aaa.Cindex.size();m++) {
				if (aaa.Pindex.at(n) == aaa.Cindex.at(m)) {
					t.push_back(m);
				}
			}
		}
		ref=-1;
		for (n=1;n<aaa.Cindex.size();n++) {
			if (aaa.Cindex.at(n) != aaa.Cindex.at(n-1)+1) {
				ref = n;
				aaa.Cindex.at(ref) = aaa.Cindex.at(ref-1) + ci_ref.size() + 1;
				break;
			}
		}
		if (ref!=-1) {
			for (n=ref+1;n<aaa.Cindex.size();n++) {
				aaa.Cindex.at(n) = aaa.Cindex.at(n-1) + 1;
			}

			for (n=1;n<aaa.Cindex.size();n++) {
				aaa.Pindex.at(n) = aaa.Cindex.at(t.at(n-1));
			}
		}
	}
	n=0;
	ref=1;
	while(n<cj_ref.size()) {
		for (m=0;m<Cindex.size();m++) {
			if (Cindex.at(m) > cj_ref.at(n)) { // Cindex[m] > cj_ref[n] && Cindex[m-1] < cj_ref[n]  orig // 20191130
				if (Cindex.at(m-1) < cj_ref.at(n) && m >= 1) {
					Cindex.insert(Cindex.begin() + m, cj_ref.at(n));
					Pindex.insert(Pindex.begin() + m, pj_ref.at(n));
					Mindex.insert(Mindex.begin() + m, mj_ref.at(n));
					Rindex.insert(Rindex.begin() + m, rj_ref.at(n));
					Cyindex.insert(Cyindex.begin() + m, cyj.at(n));
					if (para.protect) protect.insert(protect.begin() + m, prj.at(n));
					n++;
					break;
				}
			} 
			else if (Cindex.at(Cindex.size()-1) < cj_ref.at(n)) {
				Cindex.push_back(cj_ref.at(n));
				Pindex.push_back(pj_ref.at(n));
				Mindex.push_back(mj_ref.at(n));
				Rindex.push_back(rj_ref.at(n));
				Cyindex.push_back(cyj.at(n));
				if (para.protect) protect.push_back(prj.at(n));
				n++;
				break;
			}
		}
	}
	m=0;
	ref=1;
	while (m<ci_ref.size()) {
		for (n=0;n<aaa.Cindex.size();n++) {
			if (aaa.Cindex.at(n) > ci_ref.at(m)) { // aaa.Cindex[n] > ci_ref[m] && aaa.Cindex[n-1] < ci_ref[m]  orig  //20191130
				if (aaa.Cindex.at(n-1) < ci_ref.at(m) && n >= 1) {
					ref=n;
					aaa.Cindex.insert(aaa.Cindex.begin() + ref, ci_ref.at(m));
					aaa.Pindex.insert(aaa.Pindex.begin() + ref, pi_ref.at(m));
					aaa.Mindex.insert(aaa.Mindex.begin() + ref, mi_ref.at(m));
					aaa.Rindex.insert(aaa.Rindex.begin() + ref, ri_ref.at(m));
					aaa.Cyindex.insert(aaa.Cyindex.begin() + ref, cyi.at(m));
					if (para.protect) aaa.protect.insert(aaa.protect.begin() + ref, pri.at(m));
					m++;
					break;
				}
			} 
			else if (aaa.Cindex.at(aaa.Cindex.size()-1) < ci_ref.at(m)) {
				ref=n;
				aaa.Cindex.push_back(ci_ref.at(m));
				aaa.Pindex.push_back(pi_ref.at(m));
				aaa.Mindex.push_back(mi_ref.at(m));
				aaa.Rindex.push_back(ri_ref.at(m));
				aaa.Cyindex.push_back(cyi.at(m));
				if (para.protect) aaa.protect.push_back(pri.at(m));
				m++;
				break;
			}
		}
	}

    vector<int>().swap(t);
    vector<int>().swap(f);
    vector<int>().swap(ci_ref);
    vector<int>().swap(mi_ref);
    vector<int>().swap(ri_ref);
    vector<int>().swap(cj_ref);
    vector<int>().swap(mj_ref);
    vector<int>().swap(rj_ref);
    vector<int>().swap(prj);
    vector<int>().swap(pri);


	//if (1) ring_no_chk(-1);
	
	del_unpaired_ring_no();
	aaa.del_unpaired_ring_no();
	decyc_small_ring(5);
	aaa.decyc_small_ring(5);

	reset();
	aaa.reset();

	//time1=time(NULL)-time1;
	//cout << "CROSSOVER FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}


void MOLECULE::read(string a) {
	ifstream inf((a+".mds").c_str());
	string b;
	int i,j,k;
	inf >> ws;
	while (!inf.eof()) {
		b="";
		inf>>b>>ws;
		if (b=="natom") inf>>j>>ws;
		else if (b=="Pindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Pindex.push_back(k);
			}
		}
		else if (b=="Cindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Cindex.push_back(k);
			}
		}
		else if (b=="Cyindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Cyindex.push_back(k);
			}
		}
		else if (b=="Rindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Rindex.push_back(k);
			}
		}
		else if (b=="if_circle") {
			inf>>k>>ws;
			if_circle=k;
		}
		else if (b=="Mindex") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				Mindex.push_back(k);
			}
		}
		else if (para.protect && b=="protection") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				protect.push_back(k);
			}
		}
	}
	inf.close();
	return;
}


int MOLECULE::ring(int pt1,int pt2) {
	if (if_circle>=9) return 0; //20200510
	if (pt1==pt2) return 0;
	if (Cyindex.at(pt1)) return 0;
	if (Cyindex.at(pt2)) return 0;
	if (para.protect && protect.at(pt1)) return 0;
	if (para.protect && protect.at(pt2)) return 0;

	int i,j,k,n=0,m=0,x,y;
	int b_pos[2][2];

	//long double time1=time(NULL);

	reset();

	b_pos[0][0]=0;
	b_pos[0][1]=1;
	b_pos[1][0]=2;
	b_pos[1][1]=3;  //20191129
	int b_circle;
	k=1;

	for (i=0;i<data.a[Mindex.at(pt1)].norder;i++) {
		if (Bindex[pt1].at(i)==1) k=0; //at
	}
	if (k) return 0;
	k=1;
	for (i=0;i<data.a[Mindex.at(pt2)].norder;i++) {
		if (Bindex[pt2].at(i)==1) k=0; //at
	}
	if (k) return 0;
	i=pt1;
	j=pt2;
	if ((Pindex.at(i)!=Pindex.at(j)) && (i!=j) && Pindex.at(i)!=Cindex.at(j) && Pindex.at(j)!=Cindex.at(i) && Cyindex.at(i)==0 && Cyindex.at(j)==0) { //at
		for (x=0;x<data.a[Mindex.at(i)].norder;x++) {
			if (Bindex[i].at(x)==1) {
				for (y=0;y<data.a[Mindex.at(j)].norder;y++) {
					if (Bindex[i].at(x) == Bindex[j].at(y)) {
						b_pos[0][0]=x;
						b_pos[0][1]=y;
						b_pos[1][0]=i;
						b_pos[1][1]=j;
						b_circle=Bindex[i].at(x); //at
						Bindex[i].at(x)=Bindex[j].at(y)=0; //at
						if (0) Cyindex.at(i)=Cyindex.at(j)=if_circle+1; //at
						if (1) { //20200509
							Cyindex.at(i)=Cyindex.at(i)*10+if_circle+1;
							Cyindex.at(j)=Cyindex.at(j)*10+if_circle+1;
						}
						m=1;
						if_circle++;
						break;
					}
				}
				if (m==1) break;
			}
		}
	}
	else {
		//ofstream log("molecule.log",ios::app);
		//log<<"Warning : this molecule "<<smiles<<" cannot be cyclized at point "<<pt1<<" and "<<pt2<<"."<<endl;
		//log.close();
		return 0;
	}

    del_unpaired_ring_no();
    decyc_small_ring(5);

	reset();

	//time1=time(NULL)-time1;
    //cout << "CYCLIZATION FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}

void MOLECULE::mds2smi()
{
	int j,q,n,w,i=Pindex.size(),t,tt,M;
	//print();
	t=0;
	for (q=0;q<Cindex.size();q++) {
		j=Mindex.at(q);
		j=data.a[j].name.length(); //20191011
		t+=j;
	}
	if (if_circle) t+=2*if_circle;
	t+=500;
	char *x=NULL;
	x=new char [t];
	for (q=0;q<t;q++) x[q]=' ';
	int xsize=t;

	if (if_circle) {
		n=0;
		vector<int> pos,m,num;
		char tmp[10];
		//for (int j1=0;j1<10;j1++) tmp[j1]=' ';
		for (q=0;q<i;q++) { 
			if (Cyindex.at(q) != 0) {
				m.push_back(q);
				pos.push_back(atomsmi[q].at(0));
				num.push_back(Cyindex.at(q));
			}
		}
		for (q=0;q<pos.size();q++) {
			if (q==(pos.size()-1)) break;
			for (j=0;j<(pos.size()-q);j++) {
				if (j==(pos.size()-1)) break;
				if (pos.at(j)>pos.at(j+1)) {  //at
					swap(pos.at(j),pos.at(j+1)); //at
					swap(m.at(j),m.at(j+1)); //at
					swap(num.at(j),num.at(j+1)); //at
				}
			}
		}
		stringstream nn;
		string nm="";
		for (w=0;w<m.size();w++) {
			if (data.a[Mindex.at(m.at(w))].chg) { //at
				if (1) { //20200509
					tt=(int)log10(num.at(w))+1;
					t=atomsmi[m.at(w)].at(5)+1; //at
					nn.str("");
					nn << num.at(w);
					nm=nn.str();
					for (int g=0;g<tt;g++) x[atomsmi[m.at(w)].at(5)+g+1] = nm[g];
				}
				if (0) {
					if (num.at(w)>10) {
						t=atomsmi[m.at(w)].at(5)+1; //at
						tt=2;
						sprintf(tmp,"%d%d",num.at(w));
						x[atomsmi[m.at(w)].at(5)+1]=tmp[0];
						x[atomsmi[m.at(w)].at(5)+2]=tmp[1];
					}
					else {
						tt=1;
						t=atomsmi[m.at(w)].at(5)+1;
						sprintf(tmp,"%d",num.at(w));
						x[atomsmi[m.at(w)].at(5)+1]=tmp[0];
					}
				}
			}
			else {
				if (1) {//20200509
                	tt=(int)log10(num.at(w))+1;
                	t=atomsmi[m.at(w)].at(0)+1; //at
                	nn.str("");
                	nn << num.at(w);
                	nm=nn.str();
					for (int g=0;g<tt;g++) x[atomsmi[m.at(w)].at(0)+g+1] = nm[g];
				}
				if (0) {
					if (num.at(w)>10) {
						t=atomsmi[m.at(w)].at(0)+1;
						tt=2;
						sprintf(tmp,"%d%d",num.at(w));
						x[atomsmi[m.at(w)].at(0)+1]=tmp[0];
						x[atomsmi[m.at(w)].at(0)+2]=tmp[1];
					}
					else {
						tt=1;
						t=atomsmi[m.at(w)].at(0)+1;
						sprintf(tmp,"%d",num.at(w));
						x[atomsmi[m.at(w)].at(0)+1] = tmp[0];
					}
				}
			}
			for (q=0;q<i;q++) {
				j=Mindex.at(q);
				j=data.a[j].name.length();
				for (n=0;n<j;n++) {
					if (atomsmi[q].at(n)>=t) atomsmi[q].at(n) = atomsmi[q].at(n)+tt;
				}
			}
		}
		pos.clear();
		m.clear();
		num.clear();
		vector<int>().swap(pos);
		vector<int>().swap(m);
		vector<int>().swap(num);
	}
	n=0;
	chg=0;
	for (q=0;q<Cindex.size();q++) chg+=data.a[Mindex.at(q)].chg;
	t=0;
	int tmp=0;
	string kk="";
	int slashnum=0;
	for (q=0;q<i;q++) {
		j = Mindex.at(q);
		M = Mindex.at(q);

		kk=data.a[M].name; //20191011
		j=kk.length(); //20191011

		t=t+j;
		tmp=0;

		if (data.a[M].chg && kk.length()<=data.a[M].nbond && data.a[M].index!=data.a[M].nbond) { // data.a[M].chg && kk.length()<=data.a[M].nbond && data.a[M].name[0]=='['
			for (tt=0;tt<data.a[M].norder;tt++) {
				tmp+=Bindex[q].at(tt);
			}
			if (0) {
				if (tmp==1) data.a[M].name[3]='1';
				else if (tmp==2) data.a[M].name[3]='2';
				else if (tmp==3) data.a[M].name[3]='3';
				else if (tmp==4) data.a[M].name[3]='4';
				else if (tmp==5) data.a[M].name[3]='5';
			}
			if (1) {
				if (tmp>0 && tmp<10) data.a[M].name[3]=char(tmp+48);
			}
		}
		tmp=data.a[Mindex.at(q)].index+1;  //at
		for (w=0;w<data.a[M].name.length();w++) { // data.a[M].nbond  //w<j
			n=atomsmi[q].at(w);
			if (w<data.a[M].nbond) {
				if (((w-tmp)%3)==0 && w>data.a[Mindex.at(q)].index && data.a[Mindex.at(q)].index<data.a[Mindex.at(q)].nbond) { // && data.a[Mindex.at(q)].index
					if (atomsmi[q].at(w-2)==(n-2) && atomsmi[q].at(w-1)==(n-1)) {
						if (0) {
							x[n+slashnum]=' ';
							x[n-1+slashnum]=' ';
							x[n-2+slashnum]=' ';
						}
						if (1) {
							x[n]=' ';
							x[n-1]=' ';
							x[n-2]=' ';
						}
					}
					else {
						if (0) x[n+slashnum]=data.a[Mindex.at(q)].name[w];
						if (1) x[n]=data.a[Mindex.at(q)].name[w];
					}
				}
				else {
					if (0) x[n+slashnum]=data.a[Mindex.at(q)].name[w];
					if (1) x[n]=data.a[Mindex.at(q)].name[w];
				}
			}
			else if (w>=data.a[M].nbond && w<data.a[M].name.length()) {
				if (0) x[n+slashnum]=data.a[Mindex.at(q)].name[w];
				if (1) x[n]=data.a[Mindex.at(q)].name[w];
			}
		}

		if (0) {
			if (0) x[n+1+slashnum]=x[n+2+slashnum]=x[n+3+slashnum]=x[n+4+slashnum]=' ';
			if (0) x[n+1]=x[n+2]=x[n+3]=x[n+4]=' ';
			int v=kk.length()-data.a[Mindex.at(q)].nbond,count=0;
			if (v>0) {
				while (count<v) {
					if (1) x[n+count+1+slashnum]=kk[data.a[Mindex.at(q)].nbond+count];  // 4 blanks
					if (0) x[n+count+1]=kk[data.a[Mindex.at(q)].nbond+count];
					count++;
				}
			}
			if (0) if(data.a[Mindex.at(q)].nbond<kk.length() && q>0) x[n+count+1+slashnum]=')';
			slashnum+=count;
		}

		data.a[M].name=kk;  //20191201
	}
	w=0;
	molesmi="";

	double a=prob();

	for (j=0;j<xsize;j++) {  //xsize
		if (0) {
			if (x[j]!=' ') {
				if (x[j]=='(') {
					if (a>=0.5) molesmi.push_back('/');
					else molesmi.push_back('\\');
					molesmi.push_back(x[j]);
				}
				else {
					molesmi.push_back(x[j]);
				}
			}
		}
		if (1) {
			if (x[j]!=' ') molesmi.push_back(x[j]);
		}
	}

	delete [] x;
	x=NULL;
	return;
}


void MOLECULE::clear() {
	int n;
	molesmi.clear();
	for (n=0;n<1024;n++) {
		Bindex[n].clear();
		atomsmi[n].clear();
	}
	return;
}


void MOLECULE::reset() {
	int i;
	int j;
	int n,m;
	int a,b;
	int c;
	int x;
	int u;
	int M;
	int smindex; // 20200509 comment ,t[2]={0,1},s[2]={2,3};
	clear();
	for (i=0;i<Cindex.size();i++) {
		b = Mindex.at(i);
		c = Pindex.at(i);
		x = Rindex.at(i);
		for (j=0;j<data.a[b].norder;j++) {
			Bindex[i].push_back(data.a[b].order[j]);
		}
		u=0;
		if (c==0) {
			for (j=0;j<data.a[b].name.length();j++) { // data.a[b].nbond //20191230
				atomsmi[i].push_back(j);
			}
		}
		else if (c>0) {
			M = Mindex.at(c-1);
			for (j=0;j<data.a[M].norder;j++) {
				if (Bindex[c-1].at(j) == x) {
					for (n=0;n<data.a[b].norder;n++) {
						if (Bindex[i].at(n) == x) {
							Bindex[i].at(n) = 0;
							Bindex[c-1].at(j) = 0;
							u=1;
							break;
						}
					}
					if (u==1) break;
				}
			}
			smindex=0;
			if (u) smindex=data.a[M].index+j*3;
			for (j=0;j<data.a[b].name.length();j++) {  // data.a[b].nbond //20191230
				a=j+atomsmi[c-1].at(smindex)+1;		
				atomsmi[i].push_back(a);
			}
			for (j=0;j<i;j++) {
				for (n=0;n<data.a[Mindex.at(j)].name.length();n++) { // data.a[Mindex.at(j)].nbond
					if (atomsmi[j].at(n) > atomsmi[c-1].at(smindex)) {
						atomsmi[j].at(n) = atomsmi[j].at(n)+data.a[b].name.length();  // + data.a[b].nbond
					}  
				}
			}
		}
	}


	if (1) {
		m=0;
		/*
		for (i=0;i<Cyindex.size();i++) {
		   if (Cyindex.at(i)>m) m=Cyindex.at(i);
		}
		*/
		//int t[2]={0,1},s[2]={2,3};
		vector<int> t(0);
		vector<int> s(0);
		for (i=0;i<Cyindex.size();i++) {
			if (Cyindex.at(i)>m) {
				if (1) { //20200509
                    int digits=(int)log10(Cyindex.at(i))+1;
                    for (int ap=1;ap<=digits;ap++) {
                        int nn=Cyindex.at(i)/(int)pow(10,ap-1);
                        nn=nn%10;
                        if (nn>m) m=nn;
                    }
				}
				if (0) {
					t.resize(2);
					t.at(0)=0;
					t.at(1)=1;
					s.resize(2);
					s.at(0)=2;
					s.at(1)=3;
					if (Cyindex.at(i)>10) {
						t[0]=Cyindex.at(i)/10;
						t[1]=Cyindex.at(i)%10;
						if (t[0]>m && t[0]>t[1]) m=t[0];
						else if (t[1]>t[0] && t[1]>m) m=t[1];
					}
					else m=Cyindex.at(i);
				}
			}
		}
		if (m==0) {
			if_circle=0;
			for (i=0;i<Cyindex.size();i++) Cyindex.at(i)=0;
			return;
		}
		if (if_circle) {
			vector<int> cyclic;
			cyclic.resize(m);
			for (i=0;i<cyclic.size();i++) cyclic.at(i)=0;
			for (i=0;i<Mindex.size();i++) {
				if (Cyindex.at(i))  {
					if (1) { //20200509
						int digits=(int)log10(Cyindex.at(i))+1;
						for (int ap=1;ap<=digits;ap++) {
							int nn=Cyindex.at(i)/(int)pow(10,ap-1);
							nn=nn%10;
							cyclic.at(nn-1)+=1;
						}
					}
					if (0) {
						if (Cyindex.at(i)>10) {
							j=Cyindex.at(i)/10;
							n=Cyindex.at(i)%10;
							cyclic.at(j-1)+=1;
							cyclic.at(n-1)+=1;
						}
						else cyclic.at(Cyindex.at(i)-1)+=1;
					}
				}
			}
			m=-1;
			for (i=0;i<cyclic.size();i++) {
				if (cyclic.at(i)==2) {
					u=0;
					for (j=0;j<Mindex.size();j++) {
						for (n=j+1;n<Mindex.size();n++) {
							if (1) { //20200509
								t.resize(0);
								s.resize(0);

                        		int digits=(int)log10(Cyindex.at(j))+1;
                        		for (int ap=1;ap<=digits;ap++) {
                            		int nn=Cyindex.at(j)/(int)pow(10,ap-1);
                            		nn=nn%10;
									t.push_back(nn);
                        		}
                                int digits1=(int)log10(Cyindex.at(n))+1;
                                for (int ap=1;ap<=digits1;ap++) {
                                    int nn=Cyindex.at(n)/(int)pow(10,ap-1);
                                    nn=nn%10;
                                    s.push_back(nn);
                                }

								
							}
							if (0) {
                    			t.resize(2);
                    			t.at(0)=0;
                    			t.at(1)=1;
                    			s.resize(2);
                    			s.at(0)=2;
                    			s.at(1)=3;

								if (Cyindex.at(j)>10) {
									t[0]=Cyindex.at(j)/10;
									t[1]=Cyindex.at(j)%10;
								}
								else {
									t[0]=t[1]=Cyindex.at(j);
								}
								if (Cyindex.at(n)>10) {
									s[0]=Cyindex.at(n)/10;
									s[1]=Cyindex.at(n)%10;
								}
								else {
									s[0]=s[1]=Cyindex.at(n);
								}
							}


							bool go1=0,go2=0;
							if (1) {
								for (int j1=0;j1<t.size();j1++) {
									for (int j2=0;j2<s.size();j2++) {
										if (t.at(j1)==s.at(j2)) {
											go1=1;
											break;
										}
									}
									if (go1) break;
								}
							}
							if (0) if (s[0]==t[0] || s[0]==t[1] || s[1]==t[0] || s[1]==t[1]) go1=1;


                            if (1) {
                                for (int j2=0;j2<s.size();j2++) {
                                    if (s.at(j2)==(i+1)) {
                                        go2=1;
                                        break;
                                    }
                                }
                            }
							if (0) if (s[0]==(i+1) || s[1]==(i+1)) go2=1;



							if (go1 && go2) { // ( (s[0]==t[0] || s[0]==t[1] || s[1]==t[0] || s[1]==t[1]) && (s[0]==(i+1) || s[1]==(i+1)) )
								for (m=0;m<data.a[Mindex.at(j)].norder;m++) {
									for (x=0;x<data.a[Mindex.at(n)].norder;x++) {
										if (Bindex[j].at(m)==1 && Bindex[n].at(x)==1) {
											Bindex[j].at(m)=0;
											Bindex[n].at(x)=0;
											u=1;
											break;
										}
									}
									if (u) break;
								}
								if (u==0) {
									u=1;
									Cyindex.at(j)=0;
									Cyindex.at(n)=0;
								}
								if (u) break;
							}
							if (u) break;
						}
						if (u) break;
					}
				} 
				else {
					for (j=0;j<Mindex.size();j++) {
						if (1) { //20200509
							if (Cyindex.at(j)) {
                                int digits=(int)log10(Cyindex.at(j))+1;
                                for (int ap=1;ap<=digits;ap++) {
                                    int nn=Cyindex.at(j)/(int)pow(10,ap-1);
                                    nn=nn%10;
									if (nn==(i+1)) {
										Cyindex.at(j)=Cyindex.at(j)%(int)pow(10,ap-1)+(Cyindex.at(j)/(int)pow(10,ap))*(int)pow(10,ap-1);									
									}
                                }
							}
						}
						if (0) {
							if (Cyindex.at(j)==(i+1)) {
								Cyindex.at(j)=0;
								break;
							}
						}
					}
				}
			}
		}

		n=0;

		/*
		   for (i=0;i<Cyindex.size();i++) {
		   if (Cyindex.at(i)) {
		   for (j=0;j<Cyindex.size();j++) {
		   if (Cyindex.at(i)==Cyindex.at(j) && i!=j) {
		   n=0;
		   break;
		   }
		   else n=1;
		   }
		   }
		   if (n) Cyindex.at(i)=0;
		   }
		   */

		if (if_circle) {
			m=0;
			for (i=0;i<Cyindex.size();i++) {
				//if (Cyindex.at(i)>m) m=Cyindex.at(i);
				if (Cyindex.at(i)>m) {
					if (1) {
                        int digits=(int)log10(Cyindex.at(i))+1;
                        for (int ap=1;ap<=digits;ap++) {
                            int nn=Cyindex.at(i)/(int)pow(10,ap-1);
                            nn=nn%10;
							if (nn>m) m=nn;
                        }
					}
					if (0) {
						if (Cyindex.at(i)>10) {
							j=Cyindex.at(i)/10;
							n=Cyindex.at(i)%10;
							if (j>n && j>m) m=j;
							else if (n>j&&n>m) m=n;
						}
						else m=Cyindex.at(i);
					}
				}
			}
			if_circle=m;
		}
	}

	if (para.protect) prct(); //20200101

	return;
}


void MOLECULE::init() {
	int i=0,j=0,k=0;
	string tmp="";
	double x,y,z;
	smi2gjf();
	if (atm!=NULL) delete [] atm;
	atm=NULL;
	vector< vector<string> > null(0);

	ifstream inf((para.smidir+"gg.txt").c_str());
	inf >> tmp >> ws;
	while (!inf.eof()) {
		tmp="";
		inf >> tmp >> ws;
		if (tmp=="V30") {
			null.resize(null.size()+1);
			null.at(null.size()-1).resize(7,"0");
			for (i=0;i<6;i++) {
				inf >> null.at(null.size()-1).at(i) >> ws;	
				//cout << null.at(null.size()-1).at(i) << " ";
			}
			do {
				inf >> tmp >> ws;
				if (tmp.length()>=4) {
					if (tmp[0]=='C' && tmp[1]=='H' && tmp[2]=='G' && tmp[3]=='=') {
						null.at(null.size()-1).at(6)=tmp.substr(4);
						//cout << null.at(null.size()-1).at(6) << " ";
					}	
				}
			} while (tmp!="M" && !inf.eof());
			//cout << endl;
		}
	}
	inf.close();

	natom=null.size();
	atm = new DEATOM [natom];

	if (natom>0) {
		for (k=0;k<natom;k++) {
			atm[k].name=null.at(k).at(1);
			atm[k].x[0]=atof(null.at(k).at(2).c_str());
			atm[k].x[1]=atof(null.at(k).at(3).c_str());
			atm[k].x[2]=atof(null.at(k).at(4).c_str());
			atm[k].chg=atoi(null.at(k).at(6).c_str());
			//atm[k].find_r();
		}

		/*
		   for (k=0;k<natom;k++) {
		   if (atm[k].name=="H") {
		   for (int k1=k+1;k1<natom;k1++) {
		   if (atm[k1].name!="H") swap(atm[k],atm[k1]);
		   }
		   }
		   }
		   */
	}
	null.clear();
	vector< vector<string> >().swap(null);

	dist = new double *[natom];
	connect = new int *[natom];
	order= new int *[natom];
	for (i=0;i<natom;i++) {
		order[i]=new int [natom];
		connect[i]=new int [natom];
		dist[i]=new double [natom];
	}

	for (i=0;i<natom;i++) {
		dist[i][i]=0.0;
	}
	k=0;
	for (i=0;i<natom;i++) {
		if (atm[i].name != "H") {
			k++;
		}
	}
	for (i=0;i<k;i++) {
		Cindex.push_back(i+1);
		Cyindex.push_back(0);
	}
	if_circle=0;
	//cal_r();
	check_bnd();
	return;
}


void MOLECULE::smi2gjf() {
	if (0) {
		ofstream out((para.smidir+smiles+".smi").c_str());
		out << smiles << endl;
		out.close();
		//system(("ssh cluster '  /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -ogjf -O "+para.smidir+"\""+smiles+".gjf\" --gen3D --conformer --converge 1000000 -c  --nconf 100 ' ").c_str());

		system(("ssh cluster ' /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D -x3 --conformer --converge 100000 --score energy --nconf 100 -c  2> /dev/null ' ").c_str());
	}
    if (1) {
        stringstream ss(smiles);
        ofstream out((para.smidir+smiles+".mol").c_str());

        OBConversion conv(&ss, &out);
        if(conv.SetInAndOutFormats("SMI","MOL"))
        {
            //conv.AddOption("h", OBConversion::GENOPTIONS);
            conv.AddOption("gen3D", OBConversion::GENOPTIONS);
            //conv.AddOption("canonical", OBConversion::GENOPTIONS);
			conv.AddOption("weighted", OBConversion::GENOPTIONS);
			conv.AddOption("3", OBConversion::OUTOPTIONS);
            conv.Convert();
        }
		out.close();
    }

	///home/akitainu/bin/openbabel-install/bin/obabel
	//system(("  /home/akitainu/bin/openbabel-install/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D -x3 --conformer --converge 100000 --score energy --nconf 100 -c  > /dev/null  ").c_str());
	system(("(grep -A10000 \"BEGIN ATOM\" "+para.smidir+"\""+smiles+".mol\" | grep -B10000 \"END ATOM\" |  tail -n +2 | tac | tail -n +2 | tac > "+para.smidir+"gg.txt) 2> /dev/null").c_str());
	return;
}


void MOLECULE::cal_r() {
	int i,j,k;
	double r;
	for (i=0;i<natom;i++) {
		for (j=i+1;j<natom;j++) {
			r = (atm[i].x[0]-atm[j].x[0])*(atm[i].x[0]-atm[j].x[0])+(atm[i].x[1]-atm[j].x[1])*(atm[i].x[1]-atm[j].x[1])+(atm[i].x[2]-atm[j].x[2])*(atm[i].x[2]-atm[j].x[2]);
			r = sqrt(r);
			dist[i][j] = dist[j][i] = r;
		}
	}
	return;
}


void MOLECULE::check_bnd(int aaa) {
	int i,j,k,sum[natom];
	double r,tol=1.15;
	for (i=0;i<natom;i++) {
		connect[i][i]=0;
		for (j=i+1;j<natom;j++) connect[i][j]=connect[j][i]=0;
	}

	if (1) {
		// /home/akitainu/bin/openbabel-install/bin/obabel

		if (0) {
			system(("ssh cluster ' /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D -x3 --conformer --converge 100000 --score energy --nconf 100 -c  2> /dev/null ' ").c_str());
		}
		if (1) {
        	stringstream ss(smiles);
        	ofstream out((para.smidir+smiles+".mol").c_str());

        	OBConversion conv(&ss, &out);
        	if(conv.SetInAndOutFormats("SMI","MOL"))
        	{
            	//conv.AddOption("h", OBConversion::GENOPTIONS);
            	conv.AddOption("gen3D", OBConversion::GENOPTIONS);
            	//conv.AddOption("canonical", OBConversion::GENOPTIONS);
            	conv.AddOption("weighted", OBConversion::GENOPTIONS);
            	conv.AddOption("3", OBConversion::OUTOPTIONS);
            	conv.Convert();
        	}
        	out.close();
		}

		//system(("  /home/akitainu/bin/openbabel-install/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D -x3 --conformer --converge 100000 --score energy --nconf 100 -c  > /dev/null  ").c_str());
		//grep -A10000 "BEGIN BOND" u1.mol | grep -B10000 "END BOND" | tail -n +2 | tac | tail -n +2 | tac | awk '{print " " $5 " " $6 " " $4}'
		system(("(grep -A10000 \"BEGIN BOND\" "+para.smidir+"\""+smiles+".mol\" | grep -B10000 \"END BOND\" |  tail -n +2 | tac | tail -n +2 | tac | awk '{print \" \" $5 \" \" $6 \" \" $4}' > "+para.smidir+"gg.txt) 2> /dev/null").c_str());

		//for (i=0;i<natom;i++) {
		//    k=-1;
		k=rd_nps(i);
		//    if (k!=-1 && k!=atm[i].nbnd) atm[i].nbnd=k;	
		//}
		for (k=0;k<natom;k++) {
			if (atm[k].name=="H") {
				for (int k1=k+1;k1<natom;k1++) {
					if (atm[k1].name!="H") {
						swap(atm[k],atm[k1]);
						for (int k2=0;k2<natom;k2++) {
							if (connect[k][k2]!=connect[k1][k2]) swap(connect[k][k2],connect[k1][k2]);
							if (connect[k2][k]!=connect[k2][k1]) swap(connect[k2][k],connect[k2][k1]);
						}
						//for (int k2=0;k2<natom;k2++) {
						//    swap(connect[k2][k],connect[k2][k1]);
						//}
					}
				}
			}
		}
		for (k=0;k<natom;k++) {
			if (atm[k].chg==0 && atm[k].name!="H") {
				for (int k1=k+1;k1<natom;k1++) {
					if (atm[k1].chg!=0 && atm[k1].name!="H") {
						//cout << "NAME " << atm[k].name << " " << atm[k1].name << endl;
						swap(atm[k],atm[k1]);
						for (int k2=0;k2<natom;k2++) {
							if (connect[k][k2]!=connect[k1][k2]) swap(connect[k][k2],connect[k1][k2]);
						}
						for (int k2=0;k2<natom;k2++) {
							if (connect[k2][k]!=connect[k2][k1]) swap(connect[k2][k],connect[k2][k1]);
						}
					}
				}
			}
		}

	}

	for (i=0;i<natom;i++) {
		sum[i]=0;
		for (j=0;j<natom;j++) sum[i]+=connect[i][j];
	}

	for (i=0;i<natom;i++) {
		for (j=0;j<natom;j++) order[i][j]=connect[i][j];
	}
	k=0;

	for (i=0;i<natom;i++) {
		for (j=0;j<natom;j++) {
			connect[i][j]=order[i][j];
		}
	}

	return;
}


void MOLECULE::check_bnd() {
	int i,j,k,sum[natom],chgg=0;;
	double r,tol=1.15;

	for (i=0;i<natom;i++) {
		connect[i][i]=0;
		for (j=i+1;j<natom;j++) connect[i][j]=connect[j][i]=0;
	}

	if (1) {
		// /home/akitainu/bin/openbabel-install/bin/obabel

		if (0) {
			system(("ssh cluster ' /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D -x3 --weighted --conformer --converge 100000 --score energy --nconf 100 -c  2> /dev/null ' ").c_str());
		}
        if (1) {
            stringstream ss(smiles);
            ofstream out((para.smidir+smiles+".mol").c_str());

            OBConversion conv(&ss, &out);
            if(conv.SetInAndOutFormats("SMI","MOL"))
            {
                //conv.AddOption("h", OBConversion::GENOPTIONS);
                conv.AddOption("gen3D", OBConversion::GENOPTIONS);
                //conv.AddOption("canonical", OBConversion::GENOPTIONS);
                conv.AddOption("weighted", OBConversion::GENOPTIONS);
                conv.AddOption("3", OBConversion::OUTOPTIONS);
                conv.Convert();
            }
            out.close();
        }


		//system(("  /home/akitainu/bin/openbabel-install/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D -x3 --weighted --conformer --converge 100000 --score energy --nconf 100 -c  > /dev/null  ").c_str());
		//grep -A10000 "BEGIN BOND" u1.mol | grep -B10000 "END BOND" | tail -n +2 | tac | tail -n +2 | tac | awk '{print " " $5 " " $6 " " $4}'
		system(("(grep -A10000 \"BEGIN BOND\" "+para.smidir+"\""+smiles+".mol\" | grep -B10000 \"END BOND\" |  tail -n +2 | tac | tail -n +2 | tac | awk '{print \" \" $5 \" \" $6 \" \" $4}' > "+para.smidir+"gg.txt) 2> /dev/null").c_str());

		//for (i=0;i<natom;i++) {
		//    k=-1;
		k=rd_nps(i);
		//	if (k!=-1 && k!=atm[i].nbnd) atm[i].nbnd=k;
		//}

		if (0) {
			for (k=0;k<natom;k++) {
				if (atm[k].name=="H") {
					for (int k1=k+1;k1<natom;k1++) {
						if (atm[k1].name!="H") {
							swap(atm[k],atm[k1]);
							for (int k2=0;k2<natom;k2++) {
								if (connect[k][k2]!=connect[k1][k2]) swap(connect[k][k2],connect[k1][k2]);
							}
							for (int k2=0;k2<natom;k2++) {
								if (connect[k2][k]!=connect[k2][k1]) swap(connect[k2][k],connect[k2][k1]);
							}
						}
					}
				}
			}
		}

		if (1) {
        	int nonH=0;
        	for (k=0;k<natom;k++) {
            	if (atm[k].name!="H") nonH++;
        	}
		
			//long double time1=time(NULL);

			for (k=0;k<nonH;k++) {
				if (atm[k].name=="H") {
					for (int k1=nonH;k1<natom;k1++) {
						if (atm[k1].name!="H") {
							swap(atm[k],atm[k1]);
                            for (int k2=0;k2<natom;k2++) {
                                if (connect[k][k2]!=connect[k1][k2]) swap(connect[k][k2],connect[k1][k2]);
                            }
                            for (int k2=0;k2<natom;k2++) {
                                if (connect[k2][k]!=connect[k2][k1]) swap(connect[k2][k],connect[k2][k1]);
                            }
							break; //20200708
						}
					}
				}
			}
			//time1=time(NULL)-time1;
            //cout << "CONT TABLE SORTING1 TIME: " << setprecision(10) << time1 << " sec." << endl;

			//time1=time(NULL);
			for (k=0;k<nonH;k++) { // k<natom
				if (atm[k].name!="H") {
					bool aall0=1;
					for (int k1=0;k1<=k;k1++) { 
						if (connect[k1][k]!=0 || connect[k-k1][k]!=0) { // orig connect[k1][k]!=0
							aall0=0;
							break;
						}
					}

					if (aall0) {
						for (int k1=k+1;k1<nonH;k1++) { //k1<natom
							if (atm[k1].name!="H") { 
								bool ball0=1;
								for (int k2=0;k2<=k;k2++) {
									if (connect[k2][k]!=0 || connect[k-k2][k]!=0) { // orig connect[k2][k]!=0
										ball0=0;
										break;
									}
								} 
								if (ball0) {
									//cout << "NAME " << atm[k].name << " " << atm[k1].name << endl;
									swap(atm[k],atm[k1]);
									for (int k2=0;k2<natom;k2++) {
										if (connect[k][k2]!=connect[k1][k2]) swap(connect[k][k2],connect[k1][k2]);
									}
									for (int k2=0;k2<natom;k2++) {
										if (connect[k2][k]!=connect[k2][k1]) swap(connect[k2][k],connect[k2][k1]);
									}
								}
							}
						}
					}
				}
			}
			//time1=time(NULL)-time1;
			//cout << "CONT TABLE SORTING2 TIME: " << setprecision(10) << time1 << " sec." << endl;
		}

		if (0) {
			cout << setw(2) << left << "NU" << " ";
			for (int k1=0;k1<natom;k1++) cout << setw(2) << left << atm[k1].name << " ";
			cout << endl;
			for (int k1=0;k1<natom;k1++) {
				cout << setw(2) << left << atm[k1].name << " ";
				for (int k2=0;k2<natom;k2++) {
					cout << setw(2) << left << connect[k1][k2] << " ";
				}
				cout << endl;
			}

		}


	}

	for (i=0;i<natom;i++) {
		sum[i]=0;
		for (j=0;j<natom;j++) sum[i]+=connect[i][j];
	}

	for (i=0;i<natom;i++) {
		for (j=0;j<natom;j++) order[i][j]=connect[i][j];
	}
	k=0;

	for (i=0;i<natom;i++) {
		for (j=0;j<natom;j++) {
			connect[i][j]=order[i][j];
		}
	}
	return;
}


void MOLECULE::smi2cod() {
	int i,j,n,t,m,u;
	vector<int> k;

	//cout << "SMI2COD" << endl;

	//empty();
	clear(); // 20190722

	for (i=0;i<natom;i++) {
		k.clear();
		if (atm[i].name == "H") continue; 
		else if (atm[i].name == "C") { 
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size() == 4) {
				if (atm[i].chg==0) Mindex.push_back(1);
			}
			else if (k.size() == 3) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==4 && atm[i].chg==0) Mindex.push_back(2);
				else if (t==3 && atm[i].chg==-1) Mindex.push_back(59);
			}
			else if (k.size()==2) {
				if (atm[i].chg==0) {
					for (n=0;n<k.size();n++) {
						if (k[n]==3 || k[n]==1) {
							Mindex.push_back(3);
							break;
						} 
						else if (k[n]==2) {
							Mindex.push_back(4);
							break;
						}
					}					
				}
			}
		}
		else if (atm[i].name=="O") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==1) {
				if (atm[i].chg==0) Mindex.push_back(6); // && k[0]==2
				else if (atm[i].chg==-1) Mindex.push_back(29);  // k[0]==1
			}
			else if (k.size()==0) {
				if (atm[i].chg==-1) Mindex.push_back(28);
			}
			else if (k.size()==2) {
				if (atm[i].chg==0) Mindex.push_back(5);
			}
		}
		else if (atm[i].name=="N") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==4) {
				if (atm[i].chg==1) Mindex.push_back(15);
			}
			if (k.size()==3) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==4 && atm[i].chg==1) Mindex.push_back(16);
				else if (t==3  && atm[i].chg==0) Mindex.push_back(7);
			}
			if (k.size()==2) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==3 && atm[i].chg==0) Mindex.push_back(8);
				else if (t==2 && atm[i].chg==-1) Mindex.push_back(35);
			}
			if (k.size()==1) {
				if (atm[i].chg==0) Mindex.push_back(9);
			}
		}
		else if (atm[i].name=="P") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==6) {
				if(atm[i].chg==-1) Mindex.push_back(30);
			}
			if (k.size()==5) {
				if (atm[i].chg==0) Mindex.push_back(31);
			}
			if (k.size()==4) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==4 && atm[i].chg==1) Mindex.push_back(17);
				else if (t==5 && atm[i].chg==0) Mindex.push_back(32);
			}
			if (k.size()==3) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==4 && atm[i].chg==1) Mindex.push_back(18);
				if (t==4 && atm[i].chg==0) Mindex.push_back(66);
				else if (t==3 && atm[i].chg==0) Mindex.push_back(21);
			}
			if (k.size()==2) {
				if (atm[i].chg==0) Mindex.push_back(22);
			}
			if (k.size()==1) {
				if (atm[i].chg==0) Mindex.push_back(23);
			}
		}
		else if (atm[i].name=="F") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==1) {
				if (atm[i].chg==0) Mindex.push_back(11);
				else if (atm[i].chg==-1) Mindex.push_back(24);
			}
			else if (k.size()==0) {
				if (atm[i].chg==-1) Mindex.push_back(24);
			}
		}
		else if (atm[i].name=="Cl") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==1) {
				if (atm[i].chg==0) Mindex.push_back(12);
				else if (atm[i].chg==-1) Mindex.push_back(25);
			}
			else if (k.size()==0) {
				if (atm[i].chg==-1) Mindex.push_back(25);
			}
			else if (k.size()==3) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==6 && atm[i].chg==0) Mindex.push_back(62);
			}
			else if (k.size()==4) {
				if (atm[i].chg==0) Mindex.push_back(62);
			}
		}
		else if (atm[i].name=="Br") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==1) {
				if (atm[i].chg==0) Mindex.push_back(13);
				else if (atm[i].chg==-1) Mindex.push_back(26);
			}
			else if (k.size()==0) {
				if (atm[i].chg==-1) Mindex.push_back(26);
			}
		}
		else if (atm[i].name=="I") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==1) {
				if (atm[i].chg==0) Mindex.push_back(14);
				else if (atm[i].chg==-1) Mindex.push_back(27);
			}
			else if (k.size()==0) {
				if (atm[i].chg==-1) Mindex.push_back(27);
			}
		}
		else if (atm[i].name=="S") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==2) {
				if (atm[i].chg==0) Mindex.push_back(19);
			}
			else if (k.size()==1) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==2 && atm[i].chg==0) Mindex.push_back(20);
				else if (t==1 && atm[i].chg==-1) Mindex.push_back(63);
			}
			else if (k.size()==3) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];				
				if (t==5 && atm[i].chg==0)  Mindex.push_back(61);
				else if (t==4 && atm[i].chg==0) Mindex.push_back(34);
				else if (t==3 && atm[i].chg==1) Mindex.push_back(65);
			}
			else if (k.size()==4) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==6 && atm[i].chg==0) Mindex.push_back(61);
			}
		}
		else if (atm[i].name=="B") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==4) {
				if (atm[i].chg==-1) Mindex.push_back(44);
			}
		}
		else if (atm[i].name=="In") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==4) {
				if (atm[i].chg==3) Mindex.push_back(57);
			}
		}
		else if (atm[i].name=="Ga") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (k.size()==4) {
				if (atm[i].chg==3) Mindex.push_back(64);
			}
		}
	}
	Rindex.push_back(-1);
	Pindex.push_back(0);

	if (0) {
		cout << "M: " << Mindex.size() << " | ";
		for (int s=0;s<Mindex.size();s++) cout << Mindex.at(s) << " ";
		cout << endl;
	}

	int *tmp=NULL;
	tmp=new int [natom];
	for (i=0;i<natom;i++) tmp[i]=-10;
	n=0;
	for (i=0;i<natom;i++) {
		if (atm[i].name=="H") {
			order[i][i]=0;
			for (j=0;j<natom;j++) { //j=0;j=i+1
				order[i][j]=order[j][i]=0;
			}
		}
		else {
			tmp[n]=i;
			n++;
		}
	}


	for (i=1;i<natom;i++) { //i=1
		t=0;
		m=0;
		if (atm[i].name!="H") {
			for (int c1=0;c1<3;c1++) {
				m=0;
				for (j=i;j>=0;j--) { // j=i;j>=0;j--
					if (order[i][j]>m) {
						m=order[i][j];
					}
				}

				for (j=i;j>=0;j--) { // j=i;j>=0;j--
					if (order[i][j]==m && m>0 && atm[j].name!="H") { //order[i][j]==m
						for (int s=0;s<n;s++) {
							if (tmp[s]==j && tmp[s]!=-10) {
								Rindex.push_back(order[i][j]);
								Pindex.push_back(Cindex.at(s));
								order[j][i]=order[i][j]=0;
								t=1;
								break;
							}
						}
					}
					if (t) break;
				}
			}
		}
	}

	delete [] tmp;
	tmp=NULL;

	if (0) {
		cout << "R: " << Rindex.size() << " | ";
		for (int s=0;s<Rindex.size();s++) cout << Rindex.at(s) << " ";
		cout << endl;
		cout << "P: " << Pindex.size() << " | ";
		for (int s=0;s<Pindex.size();s++) cout << Pindex.at(s) << " ";
		cout << endl;

		cout << "Cy(fir): " << Cyindex.size() << " | ";
		for (int s=0;s<Cyindex.size();s++) cout << Cyindex.at(s) << " ";
		cout << endl;

		cout << "order(before Cy): " << endl;
		cout << setw(2) << left << "NU" << " ";
		for (int k1=0;k1<natom;k1++) cout << setw(2) << left << atm[k1].name << " ";
		cout << endl;
		for (int k1=0;k1<natom;k1++) {
			cout << setw(2) << left << atm[k1].name << " ";
			for (int k2=0;k2<natom;k2++) {
				cout << setw(2) << left << order[k1][k2] << " ";
			}
			cout << endl;
		}

	}


	n=1;
	for (i=0;i<natom;i++) {
		if (atm[i].name!="H") {
			for (j=i+1;j<natom;j++) {
				if (atm[j].name!="H" && order[i][j]!=0) {
					Cyindex.at(i)=Cyindex.at(i)*10+n;
					//cout << j << " " << natom << " " << Cyindex.size() << endl;
					Cyindex.at(j)=Cyindex.at(j)*10+n;
					order[i][j]=order[j][i]=0;
					if_circle=n;
					n++;
				}
			}
		}
	}

	if (0) {
		cout << "Cy: " << Cyindex.size() << " | ";
		for (int s=0;s<Cyindex.size();s++) cout << Cyindex.at(s) << " ";
		cout << endl;
	}

	return;
}


void MOLECULE::report() {
	int i=0,j=0;
	if (para.protect) prct();

	ofstream outf((para.smidir+smiles+".mds").c_str());
	outf<<"natom "<<setw(4)<<left<<setfill(' ')<<Cindex.size();
	outf<<endl;
	outf<<"Mindex ";
	for (i=0;i<Mindex.size();i++) outf<<setw(4)<<left<<setfill(' ')<<Mindex.at(i)<<" ";
	outf<<endl;
	outf<<"Pindex ";
	for (i=0;i<Pindex.size();i++) outf<<setw(4)<<left<<setfill(' ')<<Pindex.at(i)<<" ";
	outf<<endl;
	outf<<"Cindex ";
	for (i=0;i<Cindex.size();i++) outf<<setw(4)<<left<<setfill(' ')<<Cindex.at(i)<<" ";
	outf<<endl;
	outf<<"Rindex ";
	for (i=0;i<Rindex.size();i++) outf<<setw(4)<<left<<setfill(' ')<<Rindex.at(i)<<" ";
	outf<<endl;
	outf<<"Cyindex ";
	for (i=0;i<Cyindex.size();i++) outf<<setw(4)<<left<<setfill(' ')<<Cyindex.at(i)<<" ";
	outf<<endl;
	outf<<"if_circle "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;
	if (para.protect) {
		outf<<"protection ";
		//for (i=0;i<Cindex.size();i++) outf<<j<<" ";
		for (i=0;i<protect.size();i++) outf<<setw(4)<<left<<setfill(' ')<<protect.at(i)<<" ";
		outf<<endl;
	}
	outf.close();

	//print(); //20190801

	/*
	   Cindex.clear();
	   Mindex.clear();
	   Pindex.clear();
	   Cyindex.clear();
	   Rindex.clear();
	   */
	empty(); //20200101
	if_circle=0;
	clean();
	return;
}


void MOLECULE::clean() {
	int i;
	if(dist==NULL && connect==NULL && order==NULL && atm==NULL) return;
	for (i=0;i<natom;i++) {
		if(dist[i]!=NULL) delete [] dist[i];
		if(connect[i]!=NULL) delete [] connect[i];
		if(order[i]!=NULL) delete [] order[i];
	}
	if(dist!=NULL) delete [] dist;
	if(connect!=NULL) delete [] connect;
	if(order!=NULL) delete [] order;
	if(atm!=NULL) delete [] atm;
	dist=NULL;
	connect=NULL;
	order=NULL;
	atm=NULL;
	/*
	   if(dist.size()) return;
	   for (i=0;i<natom;i++) {
	   dist.at(i).resize(0);
	   connect.at(i).resize(0);
	   order.at(i).resize(0);
	   }
	   dist.resize(0);
	   connect.resize(0);
	   order.resize(0);
	   atm.resize(0);
	   */
	return;
}


void MOLECULE::outmds(ostream &outs) {
	int i,j,k;
	outs<<"natom "<<setw(4)<<left<<setfill(' ')<<Cindex.size();
	outs<<endl;
	outs<<"Mindex ";
	for (i=0;i<Mindex.size();i++) outs<<setw(4)<<left<<setfill(' ')<<Mindex.at(i)<<" ";
	outs<<endl;
	outs<<"Pindex ";
	for (i=0;i<Pindex.size();i++) outs<<setw(4)<<left<<setfill(' ')<<Pindex.at(i)<<" ";
	outs<<endl;
	outs<<"Cindex ";
	for (i=0;i<Cindex.size();i++) outs<<setw(4)<<left<<setfill(' ')<<Cindex.at(i)<<" ";
	outs<<endl;
	outs<<"Rindex ";
	for (i=0;i<Rindex.size();i++) outs<<setw(4)<<left<<setfill(' ')<<Rindex.at(i)<<" ";
	outs<<endl;
	outs<<"Cyindex ";
	for (i=0;i<Cyindex.size();i++) outs<<setw(4)<<left<<setfill(' ')<<Cyindex.at(i)<<" ";

	outs<<endl;
	outs<<"if_circle "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl<<endl;
	return;
}


int MOLECULE::combine(MOLECULE &B, int k,int p) {
	int i,j,tmp,counter,x;
	x=1;

	//time_t start, end;
	//time(&start);
	//double time1=(double)time(NULL);//clock();//(double)time(NULL);

	for (i=0;i<data.a[Mindex[k]].norder;i++) {
		for (j=0;j<data.a[B.Mindex[p]].norder;j++) {
			if (Bindex[k].at(i)>=1 && Bindex[k].at(i)==B.Bindex[p].at(j)) x=0;
		}
	}
	if (x) {
		reset();
		mds2smi();
		//ofstream log("molecule.log",ios::app);
		//log<<"Warning : "<<smiles<<" cannot combine with "<<B.smiles<<" at point "<<k<<" and "<<p<<" because of different bond order."<<endl;
    	//time1=double(time(NULL)-time1);//double(clock()-time1)/double(CLOCKS_PER_SEC);//double(time(NULL)-time1);
    	//cout << "COMBINATION FINISHED: " << setw(15) << setprecision(10) << fixed << time1 << " sec." << endl;
		return 0;
	}
	if (p==0) ;
	else {
		B.recode(p);
		B.reset();
	}
	if (if_circle && B.if_circle) {
		for (i=0;i<B.Cindex.size();i++) {
			if (B.Cyindex.at(i)>0) B.Cyindex.at(i)+=if_circle;
		}
	}
	counter=0;
	while (1) {
		tmp=0;
		for (i=0;i<data.a[Mindex.at(k)].norder;i++) {
			if (Bindex[k].at(i)>0) {
				j=i; // the position
				tmp=Bindex[k].at(i); // the order
				break;
			}
		}
		x=1;
		for (i=0;i<data.a[B.Mindex.at(0)].norder;i++) {
			if (B.Bindex[0].at(i)==tmp) { //at
				x=0;
				break;
			}
		}
		if (x) {
			B.clear();
			B.empty();
			B.input();
    		//time1=double(time(NULL)-time1);//double(clock()-time1)/double(CLOCKS_PER_SEC);//(double)time(NULL)-time1;
    		//cout << "COMBINATION FINISHED: " << setw(15) << setprecision(10) << fixed << time1 << " sec." << endl;

			return 0;
		}
		if (tmp!=0) {
			x=1;
			for (i=0;i<data.a[B.Mindex.at(0)].norder;i++) {
				if (B.Bindex[0].at(i)==tmp) { //at
					x=0;
					break;
				}
			}
			if (1) {
				B.Rindex.at(0)=tmp; //at
				int tmpp=Cindex.size();
				for (i=0;i<B.Cindex.size();i++) {
					B.Cindex.at(i)+=tmpp;  //at
					if (i==0) B.Pindex.at(i)=Cindex.at(k); //at
					else B.Pindex.at(i)+=tmpp; //at
				}
				for (i=0;i<B.Cindex.size();i++) {
					Pindex.push_back(B.Pindex.at(i));
					Cindex.push_back(B.Cindex.at(i));
					Rindex.push_back(B.Rindex.at(i));
					Cyindex.push_back(B.Cyindex.at(i));
					Mindex.push_back(B.Mindex.at(i));
					if (para.protect) protect.push_back(B.protect.at(i));
				}
				B.clear();
				B.empty();
				B.input();
				if_circle=B.if_circle+if_circle;
				reset();
				mds2smi();

    			//time1=double(time(NULL)-time1);//double(clock()-time1)/double(CLOCKS_PER_SEC);//(double)time(NULL)-time1;
    			//cout << "COMBINATION FINISHED: " << setw(15) << setprecision(10) << fixed << time1 << " sec." << endl;
				return 1;
			}
		}
		else if (tmp==0) {
			if (p) {
				B.clear();
				B.empty();
				B.input();
			}
			reset();
			mds2smi();
    		//time1=double(time(NULL)-time1);//double(clock()-time1)/double(CLOCKS_PER_SEC);//(double)time(NULL)-time1;
    		//cout << "COMBINATION FINISHED: " << setw(15) << setprecision(10) << fixed << time1 << " sec." << endl;

			return 0;
		}
	}

    //time1=double(time(NULL)-time1);//double(clock()-time1)/double(CLOCKS_PER_SEC);//(double)time(NULL)-time1;
    //cout << "COMBINATION FINISHED: " << setw(15) << setprecision(10) << fixed << time1 << " sec." << endl;

	return 1;
}


void MOLECULE::recode(int k) {
	if (smiles=="" && molesmi!="null") smiles=molesmi;
	empty();
	//clear();
	clean();
	init(k);
	//check_bnd(k); 
	//check_bnd();
	smi2cod();
	report();
	read(smiles);
	reset();
	mds2smi();
	canonicalize_SMILES(); //20200706
	return;
}


void MOLECULE::wipe() {
	vector<int>().swap(Cindex);
	vector<int>().swap(Pindex);
	vector<int>().swap(Cyindex);
	vector<int>().swap(Rindex);
	vector<int>().swap(Mindex);
	vector<int>().swap(protect);
	int i;
	for (i=0;i<1024;i++) {
		vector<int>().swap(Bindex[i]);
		vector<int>().swap(atomsmi[i]);
	}
	delete [] atm;
	atm=NULL;
	//atm.resize(0);
	//vector<DEATOM>().swap(atm);
	molesmi.clear();
	smiles.clear();
	return;
}


void MOLECULE::init(int p) {
	int i=0,j=0,k=0,n=0,s=0,q=0;

	string trash="";
	//MOLECULE *tmp=NULL;
	//tmp=new MOLECULE [1];
	double x,y,z;
	smi2gjf();
	clean();
	atm=NULL;
	string temp="";
	//tmp[0].smiles=smiles;
	//tmp[0].init();
	vector< vector<string> > null(0);

	//vector<int> point(tmp[0].Cindex.size(),0);
	ifstream inf((para.smidir+"gg.txt").c_str());

	inf >> temp >> ws;
	while (!inf.eof()) {
		temp="";
		inf >> temp >> ws;
		if (temp=="V30") {
			null.resize(null.size()+1);
			null.at(null.size()-1).resize(7,"0");
			for (i=0;i<6;i++) inf >> null.at(null.size()-1).at(i) >> ws;
			do {
				inf >> temp >> ws;
				if (temp.length()>=4) {
					if (temp[0]=='C' && temp[1]=='H' && temp[2]=='G' && temp[3]=='=') {
						null.at(null.size()-1).at(6)=temp.substr(4);
					}	
				}
			} while (temp!="M" && !inf.eof());
		}
	}
	inf.close();

	natom=null.size();
	atm = new DEATOM [natom];

	if (natom>0) {
		for (k=0;k<natom;k++) {
			atm[k].name=null.at(k).at(1);
			atm[k].x[0]=atof(null.at(k).at(2).c_str());
			atm[k].x[1]=atof(null.at(k).at(3).c_str());
			atm[k].x[2]=atof(null.at(k).at(4).c_str());
			atm[k].chg=atoi(null.at(k).at(6).c_str());
			//atm[k].find_r();
		}
		/*
		   for (k=0;k<natom;k++) {
		   if (atm[k].name=="H") {
		   for (int k1=k+1;k1<natom;k1++) {
		   if (atm[k1].name!="H") swap(atm[k],atm[k1]);
		   }
		   }
		   }
		   */
	}

	//tmp[0].clean();
	//delete [] tmp;
	//tmp=NULL;

	vector< vector<string> >().swap(null);

	dist = new double *[natom];
	connect = new int *[natom];
	order= new int *[natom];
	for (i=0;i<natom;i++) {
		order[i]=new int [natom];
		connect[i]=new int [natom];
		dist[i]=new double [natom];
	}


	for (i=0;i<natom;i++) {
		dist[i][i]=0.0;
	}
	k=0;
	for (i=0;i<natom;i++) {
		if (atm[i].name != "H") {
			k++;
		}
	}
	for (i=0;i<k;i++) {
		Cindex.push_back(i+1);
		Cyindex.push_back(0);
	}
	if_circle=0;
	//cal_r();
	check_bnd();
	return;
}


void MOLECULE::input() {
	if (smiles=="" && molesmi!="null") smiles=molesmi;
	empty();
	/*
	if (atm!=NULL) {
	   delete [] atm;
	   atm=NULL;
	}
	*/
	//if (atm.size()) atm.resize(0);
	clean();
	init();
	//cout << "init() finished" << endl;
	//check_bnd();
	smi2cod();
	//cout << "smi2cod() finished" << endl;
	report();
	//cout << "report() finished" << endl;
	read(smiles);
	//cout << "read() finished" << endl;
	reset();
	mds2smi();
	canonicalize_SMILES(); //20200706
	system("rm ./*.smi ./*.mds ./*.mol 2> /dev/null");
	return;
}


void MOLECULE::empty() {
	Cindex.resize(0);
	Pindex.resize(0);
	Mindex.resize(0);
	Rindex.resize(0);
	Cyindex.resize(0);
	protect.resize(0);
	for (int i=0;i<1024;i++) {
		Bindex[i].resize(0);
		atomsmi[i].resize(0);
	}
	return;
}


int MOLECULE::rd_nps(int a) {
	//system(("ssh cluster ' /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D --weighted --conformer --converge 100000 --score energy --nconf 100 -c  2> /dev/null ' ").c_str());
	int i=0,j=0,k=0,l=0,s=0,m,n,p,q,x=0,y=0;

	if (1) {
		ifstream inf1((para.smidir+"gg.txt").c_str());
		a++;
		inf1 >> ws;
		for (i=0;i<natom;i++) atm[i].nbnd=0;
		while (inf1.is_open() && !inf1.eof()) {
			inf1 >> j >> ws >> k >> ws >> l >> ws;
			connect[j-1][k-1]=connect[k-1][j-1]=l;
			atm[j-1].nbnd+=l;
			atm[k-1].nbnd+=l;
			if (j==a || k==a) s+=l;
		}
		inf1.close();
	}


	if (0) {
		ifstream inf((para.smidir+smiles+".mol").c_str());
		string tmp;
		a++;

		if (1) {
			if (inf.is_open()) {
				inf>>ws;
				getline(inf,tmp);
				inf>>ws;
				inf>>x>>ws>>y>>ws;
				getline(inf,tmp);
				for (i=0;i<x;i++) getline(inf,tmp);
				for (i=0;i<y;i++) {
					inf>>j>>k>>l>>m>>n>>p>>q>>ws;
					if (j==a || k==a) s+=l;
				}
			}
		}

		if (0) {
			if (inf.is_open()) {
				a++;
				inf>>ws;
				getline(inf,tmp);
				inf>>ws;
				getline(inf,tmp);
				for (i=0;i<natom;i++) getline(inf,tmp);  // i<natom
				inf >> ws;
				//for (i=0;i<natom;i++) { // i<natom
				while (!inf.eof()) {
					string a1="nan",b1="nan";
					double g1=a+1,g2=a+2;
					int cur=inf.tellg();
					inf >> a1 >> ws >> b1 >> ws;
					g1=atof(a1.c_str());
					g2=atof(b1.c_str());
					inf.seekg(cur);
					if (!isnan(g1) && !isnan(g2)) {
						inf>>j>>ws>>k>>ws>>l>>ws>>m>>ws>>n>>ws>>p>>ws>>q>>ws;
						if (j==a || k==a) s+=l;
					}
					else {
						getline(inf,a1);
						inf >> ws;
					}
				}
				//}
			}
		}
		inf.close();
	}

	return s;
}



int MOLECULE::add(int pt, int id) {
	int i,j,n=0;
	for (i=0;i<data.a[Mindex.at(pt)].norder;i++) {
		for (j=0;j<data.a[id].norder;j++) {
			if (Bindex[pt].at(i)==data.a[id].order[j]) { //at 
				Pindex.push_back(Cindex.at(pt));
				Cindex.push_back(Cindex.size()+1);
				Mindex.push_back(id);
				Rindex.push_back(Bindex[pt].at(i)); //at
				Cyindex.push_back(0);
				if (para.protect) protect.push_back(0);
				n=1;
				break;
			}
		}
		if (n) break;
	}
	if (n==0) {
		//ofstream log("molecule.log",ios::app);
		//log<<"Warning : "<<"for "<<smiles<<", it can't excute addition at this point "<<pt<<" with this atom (id), "<<id<<"."<<endl;
		//log.close();
		reset();
		return 0;
	}
	reset();
	return 1;
}


int MOLECULE::subtract(int n,int mode) {
	// mode 1: subtract an atom
	// mode 2: subtract all atoms between Cindex=n and the end of the corresponding branch
	int i,j,k;
	vector<int> ref,pr;

	//long double time1=time(NULL);

	if (n>=Cindex.size()) return 0;
	if (n<=0) {
		//ofstream log("molecule.log",ios::app);
		//log<<"Warning : This molecule would be empty because of point, "<<n<<endl;
		//log.close();
		return 0;
	}

	ref.push_back(Cindex.at(n));
	//ref.push_back(n);  // wrong? 20190721

	if (0) {
		for (i=0;i<Cindex.size();i++) {
			//if (Cindex.at(i)>n) {
			j=ref.size();
			for (k=0;k<j;k++) {
				bool go=0;
				//if (Pindex.at(i)==ref.at(k) && mode==1) go=1;
				//if (Pindex.at(i)>=ref.at(k) && mode==2) go=1;
				if (Cindex.at(i)==ref.at(k) && mode==1) go=1;
				if (Cindex.at(i)>=ref.at(k) && mode==2) go=1;
				if (go) {
					ref.push_back(Cindex.at(i));
					if (para.protect) pr.push_back(protect.at(i));
					if (1) if (para.protect && protect.at(i)) return 0;   // terminate if it would subtract protected atom
					if (1) if (data.a[Mindex.at(i)].chg) return 0;  // terminate if it would subtract charged atom
				}
			}
			//}
		}
	}
	//for (i=0;i<Cindex.size();i++) {
		//if (Cindex.at(i)==Cindex.at(n) && mode==1) {
            if (para.protect) pr.push_back(protect.at(n));
            if (1) if (para.protect && protect.at(n)) return 0;   // terminate if it would subtract protected atom
            if (1) if (data.a[Mindex.at(n)].chg) return 0;  // terminate if it would subtract charged atom
		//}
	//}

	if (para.protect) {
		for (i=0;i<pr.size();i++) {
			if (pr.at(i)) return 0;
		}		
	}

	vector<int> Merge_C(0);
	if (mode==1) {  //20200131
		bool totbndchk=0,bndchk=0;
		for (i=0;i<ref.size();i++) {
			for (j=0;j<Cindex.size();j++) {
				if (ref.at(i)==Cindex.at(j)) { // atom j to be subtracted
					for (int k1=0;k1<Cindex.size();k1++) {
						if (Cindex.at(k1)==Pindex.at(j)) { // For the unique parental atom k1 of atom j 
							Merge_C.push_back(Cindex.at(k1));
							totbndchk=0;
							bndchk=0;

							int bndsum=Rindex.at(k1);  // # of bonds of atom k1 if atom j is deleted and k1 is connected with each of the descendant atom k
							vector<int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
							if (k1>=1) ordcount.at(Rindex.at(k1)-1)+=1;  //20200712
							for (k=0;k<Cindex.size();k++) { 
								if (Pindex.at(k)==Cindex.at(j)) { // For each descendant atom k of atom j // orig Pindex.at(k)==Cindex.at(j)
									Merge_C.push_back(Cindex.at(k));
									bndsum+=Rindex.at(k);
									ordcount.at(Rindex.at(k)-1)+=1;
								}
								if (Pindex.at(k)==Cindex.at(k1) && k!=j) {
                                    bndsum+=Rindex.at(k);
                                    ordcount.at(Rindex.at(k)-1)+=1;
								}
							}

							int maxbnd=0; // max available # of bonds for parental atom k1 (calc. from pool)
							for (int k2=0;k2<6;k2++) { 
								maxbnd+=data.a[Mindex.at(k1)].order.at(k2);
							}
							if (bndsum<=maxbnd) totbndchk=1;
							else totbndchk=0;

							for (int k2=0;k2<6;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
								if (ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
									bndchk=0;
									break;
								}
								else if (k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
									bndchk=1;
								}
							}

							if (totbndchk && bndchk) break;
							else {
								cout << "Failed subtraction: (totbndchk, bndchk) = (" << totbndchk << ", " << bndchk << ")" << endl;
								return 0;
							}
						}
					}

				}
			}
		}
	}

	int Csize=Cindex.size();
	for (i=0;i<ref.size();i++) {
		for (j=0;j<Csize;j++) {
			if (ref.at(i)==Cindex.at(j)) {
				//if (1) ring_no_chk(j);
                for (k=0;k<Csize;k++) {  // 20200131 prevent unpaired ring number
                    if (0) if (Cyindex.at(j) && Cyindex.at(k)==Cyindex.at(j) && k!=j) Cyindex.at(k)=0;
					if (1) { //20200509
						if (Cyindex.at(j) && Cyindex.at(k) && k!=j) {
                    		int digits=(int)log10(Cyindex.at(j))+1;
							int digits1=(int)log10(Cyindex.at(k))+1;
                    		for (int ap=1;ap<=digits;ap++) {
								int w=Cyindex.at(j)/(int)pow(10,ap-1);
                        		w=w%10;
								for (int bp=1;bp<=digits1;bp++) {
									int w1=Cyindex.at(k)/(int)pow(10,bp-1);
									w1=w1%10;
									if (w1==w) {
										Cyindex.at(k)=Cyindex.at(k)%(int)pow(10,bp-1)+(Cyindex.at(k)/(int)pow(10,bp))*(int)pow(10,bp-1);
									}
								}
                    		}
						}
						
					}
                }
                Mindex.erase(Mindex.begin()+j);
                Cindex.erase(Cindex.begin()+j);
                Pindex.erase(Pindex.begin()+j);
                Rindex.erase(Rindex.begin()+j);
				Cyindex.erase(Cyindex.begin()+j);	
				if (para.protect) protect.erase(protect.begin()+j);
				Csize=Cindex.size();
			}
		}
	}
	if (0) {
		for (i=0;i<Cindex.size();i++) {
			if (Cindex.at(i)!=(i+1)) {
				j=Cindex.at(i);
				Cindex.at(i)=i+1;
				for (k=0;k<Cindex.size();k++) { //k=i orig
					if (Pindex.at(k)==j) Pindex.at(k)=Cindex.at(i);
				}
			}
		}
	}
	if (1) {
		if (Merge_C.size()>=2) {
			for (k=1;k<Merge_C.size();k++) {
				for (i=0;i<Cindex.size();i++) {
					if (Cindex.at(i)==Merge_C.at(k)) {
						Pindex.at(i)=-2;
						break;
					}
				}
			}
		}

    	for (i=0;i<Cindex.size();i++) {
        	if (Cindex.at(i)!=(i+1)) {
				j=Cindex.at(i);
				Cindex.at(i)=i+1;
				if (Merge_C.size()>=2) {
					for (k=0;k<Merge_C.size();k++) {
						if (j==Merge_C.at(k)) {
							Merge_C.at(k)=Cindex.at(i);
							break;
						}
					}
				}
                for (k=0;k<Cindex.size();k++) { //k=i orig
                    if (Pindex.at(k)==j) Pindex.at(k)=Cindex.at(i);
                }
        	}
    	}

		if (Merge_C.size()>=2) {
			for (k=1;k<Merge_C.size();k++) {
				Pindex.at(Merge_C.at(k)-1)=Merge_C.at(0);
			}
		}
		
	}

	//if (1) ring_no_chk(-1);
	if (0) {
		if (mode==1 && if_circle) { //20200131 
			for (int k1=0;k1<Cindex.size();k1++) {
				if (Cyindex.at(k1)) {
					for (int k2=k1+1;k2<Cindex.size();k2++) {
						if (Cyindex.at(k1)==Cyindex.at(k2) && k1!=k2) {
							if ( fabs(Pindex.at(k1)-Pindex.at(k2))<=4 ) {
								Cyindex.at(k1)=Cyindex.at(k2)=0;
								if_circle--;
							}
							break;
						}
					}
				}
			}
		}
	}
	vector<int>().swap(ref);
	vector<int>().swap(pr);
	vector<int>().swap(Merge_C);

	del_unpaired_ring_no();
	decyc_small_ring(5);

	reset();

    //time1=time(NULL)-time1;
    //cout << "SUBTRACTION FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}



int MOLECULE::exchange(int n,int id,int id2,int bond) {
	int i,j,k,M,P,C,R,bd[3];
	reset();
	vector<int> tmp;

	//long double time1=time(NULL);

	if (!n) return 0;
	for (i=0;i<data.a[Mindex.at(n)].norder;i++) tmp.push_back(Bindex[n].at(i));
	M=Mindex.at(n);
	P=Pindex.at(n);
	C=Cindex.at(n);
	R=Rindex.at(n);
	if (para.protect) {
		if (protect.at(n)) return 0;
		if (P>0) if (protect.at(P-1)) return 0;
	}
	//if (data.a[Mindex.at(id)].chg) return 0;
	//if (data.a[Mindex.at(id2)].chg) return 0;
	if (R+bond<1 || R+bond>3) {
		//ofstream log("molecule.log",ios::app);
		//log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
		//log.close();
		reset();
		return 0;
	}
	if (R==2) {
		bd[0]=bd[1]=bd[2]=0;
		for (i=0;i<tmp.size();i++) {
			if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
		}
		tmp.clear();
		for (i=0;i<3;i++) bd[i]=data.a[M].bd[i]-bd[i];
		if (bond<0) { 
			bd[0]+=1;
			bd[1]-=1;
		}
		else if (bond>0) {
			bd[2]+=1;
			bd[1]-=1;
		}
		k=0;
		j=0;
		i=id;
		if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[M].chg==data.a[i].chg) {
			Mindex.at(n)=i;
			Rindex.at(n)=R+bond;
		}
		else {
			reset();
			//ofstream log("molecule.log",ios::app);
			//log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
			//log.close();
			return 0;
		}
		bd[0]=bd[1]=bd[2]=0;
		for (i=0;i<data.a[Mindex.at(P-1)].norder;i++) tmp.push_back(Bindex[P-1].at(i));
		for (i=0;i<tmp.size();i++) {
			if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
		}
		tmp.clear();
		for (i=0;i<3;i++) bd[i]=data.a[Mindex.at(P-1)].bd[i]-bd[i];
		if (bond<0) {
			bd[0]+=1;
			bd[1]-=1;
		}
		else if (bond>0) {
			bd[2]+=1;
			bd[1]-=1;
		}
		k=0;
		j=0;
		i=id2;
		if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[Mindex.at(P-1)].chg==data.a[i].chg) {
			Mindex.at(P-1)=i;
			Rindex.at(n)=R+bond;
		}
		else {
			//ofstream log("molecule.log",ios::app);
			//log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
			//log.close();
			Rindex.at(n)=R;
			Mindex.at(n)=M;
			reset();
			return 0;
		}
	}
	else if (R==3) {
		bd[0]=bd[1]=bd[2]=0;
		for (i=0;i<tmp.size();i++) {
			if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
		}
		tmp.clear();
		for (i=0;i<3;i++) bd[i]=data.a[M].bd[i]-bd[i];
		if (bond==-1) { 
			bd[1]+=1;
			bd[2]-=1;
		}
		else if (bond==-2) {
			bd[0]+=1;
			bd[2]-=1;
		}
		k=0;
		j=0;
		i=id;
		if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[M].chg==data.a[i].chg) {
			Mindex.at(n)=i;
			Rindex.at(n)=R+bond;
		}
		else {
			//ofstream log("molecule.log",ios::app);
			//log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
			//log.close();
			reset();
			return 0;
		}
		bd[0]=bd[1]=bd[2]=0;
		for (i=0;i<data.a[Mindex.at(P-1)].norder;i++) tmp.push_back(Bindex[P-1].at(i));
		for (i=0;i<tmp.size();i++) {
			if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
		}
		tmp.clear();
		for (i=0;i<3;i++) bd[i]=data.a[Mindex.at(P-1)].bd[i]-bd[i];
		if (bond==-1) {
			bd[1]+=1;
			bd[2]-=1;
		}
		else if (bond==-2) {
			bd[0]+=1;
			bd[2]-=1;
		}
		k=0;
		j=0;
		i=id2;
		if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[Mindex.at(P-1)].chg==data.a[i].chg) {
			Mindex.at(P-1)=i;
			Rindex.at(n)=R+bond;
		}
		else {
			//ofstream log("molecule.log",ios::app);
			//log<<"Warning : for "<<smiles<<", it can't be exchanged at point, "<<n<<" to become as "<<id<<" and "<<id2<<" with add "<<bond<<" order."<<endl;
			//log.close();
			Rindex.at(n)=R;
			Mindex.at(n)=M;
			reset();
			return 0;
		}
	}
	else if (R==1) {
		bd[0]=bd[1]=bd[2]=0;
		for (i=0;i<tmp.size();i++) {
			if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
		}
		tmp.clear();
		for (i=0;i<3;i++) bd[i]=data.a[M].bd[i]-bd[i];
		if (bond==1) { 
			bd[1]+=1;
			bd[0]-=1;
		}
		else if (bond==2) {
			bd[2]+=1;
			bd[0]-=1;
		}
		k=0;
		j=0;
		i=id;
		if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[M].chg==data.a[i].chg) {
			Mindex.at(n)=i;
			Rindex.at(n)=R+bond;
		}
		else {
			reset();
			return 0;
		}
		bd[0]=bd[1]=bd[2]=0;
		for (i=0;i<data.a[Mindex.at(P-1)].norder;i++) tmp.push_back(Bindex[P-1].at(i));
		for (i=0;i<tmp.size();i++) {
			if (tmp.at(i)) bd[tmp.at(i)-1]+=1;
		}
		tmp.clear();
		for (i=0;i<3;i++) bd[i]=data.a[Mindex.at(P-1)].bd[i]-bd[i];
		if (bond==1) {
			bd[1]+=1;
			bd[0]-=1;
		}
		else if (bond==2) {
			bd[2]+=1;
			bd[0]-=1;
		}
		k=0;
		j=0;
		i=id2;
		if (bd[0]<=data.a[i].bd[0] && bd[1]<=data.a[i].bd[1] && bd[2]<=data.a[i].bd[2] && data.a[Mindex.at(P-1)].chg==data.a[i].chg) {
			Mindex.at(P-1)=i;
			Rindex.at(n)=R+bond;
		}
		else {
			Rindex.at(n)=R;
			Mindex.at(n)=M;
			reset();
			return 0;
		}
	}

    del_unpaired_ring_no(); //20200627
    decyc_small_ring(5);

	reset();

    //time1=time(NULL)-time1;
    //cout << "EXCHANGE FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}


int MOLECULE::chk_chg(MOLECULE &mol) {
	rechg();
	mol.rechg();
	int num[2]={0,0};
	char dot='.';
	ionic="";
	num[0]=fabs(mol.chg);
	num[1]=fabs(chg);
	if (num[0]==0 && num[1]==0) ionic=molesmi;
	else if (num[0]==0) ionic=mol.molesmi;
	else if (num[1]==0) ionic=molesmi;
	else {
		int i;
		for (i=0;i<num[0];i++) {
			ionic+="(";
			ionic+=molesmi;
			ionic+=")";
		}
		ionic+=dot;
		for (i=0;i<num[1];i++) {
			ionic+="(";
			ionic+=mol.molesmi;
			ionic+=")";
		}
	}
	return 1;
}


int MOLECULE::rechg() {
	int i=0,j=0;
	chg=0;
	for (j=0;j<Cindex.size();j++) {
		chg+=data.a[Mindex.at(j)].chg;
	}
	return 1;
}


void MOLECULE::neutralize() {
	int i,j,k,m,n;
	chg=0;
	for (i=0;i<Cindex.size();i++) {
		chg+=data.a[Mindex.at(i)].chg;
		if (data.a[Mindex.at(i)].chg) j=i;
	}
	if (chg) {
		if (chg<0) {
			m=0;
			while (1) {
				m++;
				n=0;
				k=rand()%Cindex.size();
				for (i=0;i<data.a[Mindex.at(k)].nbond;i++) {
					if (Bindex[k].at(i) && Bindex[k].at(i)<3) n=Bindex[k].at(i);  //at
				}
				if (k!=j && n) break;
				if (m>5) return;
			}
			if (n==2) {
				if (rand()%2) {
					Cindex.push_back(Mindex.size());
					Pindex.push_back(Cindex.at(k));
					Mindex.push_back(16);
					Cyindex.push_back(0);
					Rindex.push_back(2);
				}
				else {
					Cindex.push_back(Mindex.size());
					Pindex.push_back(Cindex.at(k));
					Mindex.push_back(18);
					Cyindex.push_back(0);
					Rindex.push_back(2);
				}
			}
			else if (n==1){
				Cindex.push_back(Mindex.size());
				Pindex.push_back(Cindex.at(k));
				Mindex.push_back(rand()%3+15);
				Cyindex.push_back(0);
				Rindex.push_back(1);
			}
		}
		else if (chg>0) {
			m=0;
			while (1) {
				n=0;
				m++;
				k=rand()%Cindex.size();
				for (i=0;i<data.a[Mindex.at(k)].nbond;i++) {
					if (Bindex[k].at(i) && Bindex[k].at(i)==1) n=Bindex[k].at(i); //at
				}
				if (k!=j && n) break;
				if (m>5) return;
			}
			if (n) {
				int tmp[3]={30,29,35};
				Cindex.push_back(Mindex.size());
				Pindex.push_back(Cindex.at(k));
				Mindex.push_back(tmp[rand()%3]);
				Cyindex.push_back(0);
				Rindex.push_back(1);
			}
		}
	}
	reset();
	mds2smi();
	return;
}


void MOLECULE::mds23d(ostream &outs) {
	int i,j,k,n;
	double sum;
	Matrix vec[Mindex.size()][6];
	int **cont;
	cont=new int *[Mindex.size()];
	for (i=0;i<Mindex.size();i++) cont[i]=new int [6];
	vector<double> x[3];
	vector<int> mm;
	vector<OPT> opt;
	opt.resize(1);
	for (i=0;i<Mindex.size();i++) mm.push_back(Mindex.at(i));
	for (i=0;i<mm.size();i++) {
		if (mm.at(i)==10) mm.at(i)=5;
		else if (mm.at(i)==28) mm.at(i)=29;
	}
	for (i=0;i<mm.size();i++) {
		k=data.a[mm.at(i)].norder;
		n=0;
		for (j=0;j<data.a[mm.at(i)].norder;j++) {
			if (Bindex[i].at(j)>1) { //at
				k+=(Bindex[i].at(j)-1); //at
				n=1;
			}
		}
		if (k>data.a[Mindex.at(i)].norder) {
			for (j=1;j<data.num;j++) {
				if (n==1 && data.a[j].atm==data.a[mm.at(i)].atm && data.a[j].norder==k && data.a[j].chg==data.a[mm.at(i)].chg) {
					mm.at(i)=j;
					break;
				}
			}
		}
	}
	for (i=0;i<Mindex.size();i++) {
		for (j=0;j<6;j++) {
			vec[i][j].resize(1,1);
			cont[i][j]=0;
		}
	}

	for (i=0;i<mm.size();i++) {
		for (j=0;j<data.a[mm.at(i)].norder;j++) {
			vec[i][j].resize(3,1);
			cont[i][j]=1;
		}
	}
	/* initailize */
	if (data.a[mm.at(0)].type==1) {
		vec[0][0].Mtrx[0][0]=0.000000;
		vec[0][0].Mtrx[1][0]=0.000000;
		vec[0][0].Mtrx[2][0]=1.000000;
		vec[0][1].Mtrx[0][0]=-0.78335;
		vec[0][1].Mtrx[1][0]=-0.4523;
		vec[0][1].Mtrx[2][0]=-0.4264;
		vec[0][2].Mtrx[0][0]=0.000000;
		vec[0][2].Mtrx[1][0]=0.9045;
		vec[0][2].Mtrx[2][0]=-0.4264;
		vec[0][3].Mtrx[0][0]=0.78335;
		vec[0][3].Mtrx[1][0]=-0.4523;
		vec[0][3].Mtrx[2][0]=-0.4264;
		vec[0][0].rotate(vec[0][1],vec[0][0],15.0*3.1415926/180.0);
		vec[0][3].rotate(vec[0][1],vec[0][3],15.0*3.1415926/180.0);
		vec[0][2].rotate(vec[0][1],vec[0][2],15.0*3.1415926/180.0);
	}
	else if (data.a[mm.at(0)].type==2) {
		Matrix ref(3,1);
		ref[0][0]=ref[1][0]=ref[2][0]=1.0;
		vec[0][0].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][1].rotate(ref,vec[0][1],120.0*3.1415926/180.0);
		vec[0][2].rotate(ref,vec[0][2],-120.0*3.1415926/180.0);
	}
	else if (data.a[mm.at(0)].type==3) {
		Matrix ref(3,1);
		ref[0][0]=ref[1][0]=ref[2][0]=1.0;
		vec[0][0].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][2].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][1].rotate(ref,vec[0][1],120.0*3.1415926/180.0);
		vec[0][2].rotate(ref,vec[0][2],-120.0*3.1415926/180.0);

	}
	else if (data.a[mm.at(0)].type==4) {
		vec[0][0].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[2][0]=1.0/sqrt(3.0);
		vec[0][1].Mtrx[0][0]=-1.0/2.0;
		vec[0][1].Mtrx[1][0]=1.0/2.0;
		vec[0][1].Mtrx[2][0]=-1.0/sqrt(2.0);
	}
	else if (data.a[mm.at(0)].type==5) {
		vec[0][0].Mtrx[0][0]=0.3;
		vec[0][0].Mtrx[1][0]=0.1;
		vec[0][0].Mtrx[2][0]=1.0000;
		vec[0][1].Mtrx[0][0]=0.3;
		vec[0][1].Mtrx[1][0]=0.1;
		vec[0][1].Mtrx[2][0]=-1.0000;
		vec[0][2].Mtrx[0][0]=-0.50000*sqrt(3.0000);
		vec[0][2].Mtrx[1][0]=-sqrt(3.0000)/6.0000*sqrt(3.00000);
		vec[0][2].Mtrx[2][0]=0.100000;
		vec[0][3].Mtrx[0][0]=0.50000*sqrt(3.0000);
		vec[0][3].Mtrx[1][0]=-sqrt(3.0000)/6.0000*sqrt(3.00000);
		vec[0][3].Mtrx[2][0]=0.300000;
		vec[0][4].Mtrx[0][0]=0.1000000;
		vec[0][4].Mtrx[1][0]=1.0000;
		vec[0][4].Mtrx[2][0]=0.30000;
		for (i=0;i<5;i++) {
			sum=0.0;
			for (j=0;j<3;j++) sum+=vec[0][i].Mtrx[j][0]*vec[0][i].Mtrx[j][0];
			sum=sqrt(sum);
			for (j=0;j<3;j++) vec[0][i].Mtrx[j][0]=vec[0][i].Mtrx[j][0]/sum;
		}
	}
	else if (data.a[mm.at(0)].type==6) {
		vec[0][0].Mtrx[0][0]=0.2;
		vec[0][0].Mtrx[1][0]=0.1;
		vec[0][0].Mtrx[2][0]=1.0000;
		vec[0][1].Mtrx[0][0]=0.3;
		vec[0][1].Mtrx[1][0]=0.5;
		vec[0][1].Mtrx[2][0]=-1.0000;
		vec[0][2].Mtrx[0][0]=1.00000;
		vec[0][2].Mtrx[1][0]=0.1;
		vec[0][2].Mtrx[2][0]=0.3;
		vec[0][3].Mtrx[0][0]=-1.0000;
		vec[0][3].Mtrx[1][0]=0.5;
		vec[0][3].Mtrx[2][0]=0.3;
		vec[0][4].Mtrx[0][0]=0.3;
		vec[0][4].Mtrx[1][0]=1.000000;
		vec[0][4].Mtrx[2][0]=0.3;
		vec[0][5].Mtrx[0][0]=0.2;
		vec[0][5].Mtrx[1][0]=-1.0000;
		vec[0][5].Mtrx[2][0]=0.3;
		for (i=0;i<6;i++) {
			sum=0.0;
			for (j=0;j<3;j++) sum+=vec[0][i].Mtrx[j][0]*vec[0][i].Mtrx[j][0];
			sum=sqrt(sum);
			for (j=0;j<3;j++) vec[0][i].Mtrx[j][0]=vec[0][i].Mtrx[j][0]/sum;
		}
	}
	else if (data.a[mm.at(0)].type==7) {
		vec[0][0].Mtrx[0][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[1][0]=1.0/sqrt(3.0);
		vec[0][0].Mtrx[2][0]=1.0/sqrt(3.0);
	}
	x[0].push_back(0.0000000000);
	x[1].push_back(0.0000000000);
	x[2].push_back(0.0000000000);
	int P=0;
	double *pt=NULL;
	pt=new double [3];
	double *xx=NULL;
	xx=new double [3];	
	int memo;
	for (i=1;i<mm.size();i++) {
		P=Pindex.at(i)-1;
		for (j=0;j<data.a[mm.at(P)].norder;j++) {
			if (cont[P][j]) {
				cont[P][j]=0;
				memo=j;
				break;
			}
		}
		cont[i][0]=0;
		for (k=0;k<3;k++) x[k].push_back(x[k].at(P)+(data.a[mm.at(P)].rb+data.a[mm.at(i)].rb)*vec[P][j].Mtrx[k][0]);
		for (k=0;k<3;k++) vec[i][0].Mtrx[k][0]=-vec[P][j].Mtrx[k][0];
		if (data.a[mm.at(i)].type==1) { /* tetrahedral */
			for (k=0;k<3;k++) pt[k]=-0.5000*vec[i][0].Mtrx[k][0];
			int tmpp[3];
			for (k=0;k<3;k++) tmpp[3]=(fabs(vec[i][0].Mtrx[k][0])>1E-5);
			for (k=0;k<3;k++) {
				if (tmpp[k]) xx[k]=pt[k]+(double)rand()/RAND_MAX*3.0;
				else xx[k]=pt[k];
			}
			sum=0.0;
			if (tmpp[0]) {
				sum=vec[i][0].Mtrx[1][0]*(pt[1]-xx[1])+vec[i][0].Mtrx[2][0]*(pt[2]-xx[2]);
				xx[0]=sum/vec[i][0].Mtrx[0][0]+pt[0];
			}
			else if (tmpp[1]) {
				sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[2][0]*(pt[2]-xx[2]);
				xx[1]=sum/vec[i][0].Mtrx[1][0]+pt[1];
			}
			else {
				sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[1][0]*(pt[1]-xx[1]);
				xx[1]=sum/vec[i][0].Mtrx[2][0]+pt[2];
			}	
			sum=0;
			sum=(xx[0]-pt[0])*(xx[0]-pt[0])+(xx[1]-pt[1])*(xx[1]-pt[1])+(xx[2]-pt[2])*(xx[2]-pt[2]);
			sum=sqrt(sum);
			sum/=(sqrt(3.0)/2.0);
			xx[0]=(xx[0]-pt[0])/sum;
			xx[1]=(xx[1]-pt[1])/sum;
			xx[2]=(xx[2]-pt[2])/sum;
			for (k=1;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].Mtrx[0][0]=xx[0];
				vec[i][k].Mtrx[1][0]=xx[1];
				vec[i][k].Mtrx[2][0]=xx[2];
			}
			for (k=2;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].rotate(vec[i][0],vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
				sum=0;
				sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
				sum=sqrt(sum);
				vec[i][k].Mtrx[0][0]/=sum;
				vec[i][k].Mtrx[1][0]/=sum;
				vec[i][k].Mtrx[2][0]/=sum;
			}
		}
		else if (0) {
			for (k=0;k<3;k++) pt[k]=x[k].at(i)-0.5000*vec[i][0].Mtrx[k][0];
			for (k=0;k<2;k++) {
				xx[k]=pt[k]+(double)rand()/RAND_MAX*3.0;
			}
			sum=0;
			sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[1][0]*(pt[1]-xx[1]);
			if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[2]=sum/vec[i][0].Mtrx[2][0]+pt[2];
			else {
				xx[2]=pt[2];
				if (fabs(vec[i][0].Mtrx[1][0])>1E-5) xx[1]=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])/vec[i][0].Mtrx[1][0]+pt[1];
				else xx[1]=pt[1];
				if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[0]=vec[i][0].Mtrx[1][0]*(pt[1]-xx[1])/vec[i][0].Mtrx[0][0]+pt[0];
				else xx[0]=pt[0];
			}
			sum=0;
			sum=(xx[0]-pt[0])*(xx[0]-pt[0])+(xx[1]-pt[1])*(xx[1]-pt[1])+(xx[2]-pt[2])*(xx[2]-pt[2]);
			sum=sqrt(sum);
			sum=sum/(sqrt(3.0)/2.0);
			xx[0]=(xx[0]-pt[0])/sum;
			xx[1]=(xx[1]-pt[1])/sum;
			xx[2]=(xx[2]-pt[2])/sum;
			for (k=1;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].Mtrx[0][0]=xx[0]-0.50000*vec[i][0].Mtrx[0][0];
				vec[i][k].Mtrx[1][0]=xx[1]-0.50000*vec[i][0].Mtrx[1][0];
				vec[i][k].Mtrx[2][0]=xx[2]-0.50000*vec[i][0].Mtrx[2][0];
			}
			for (k=2;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].rotate(vec[i][0],vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
				sum=0;
				sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
				sum=sqrt(sum);
				vec[i][k].Mtrx[0][0]/=sum;
				vec[i][k].Mtrx[1][0]/=sum;
				vec[i][k].Mtrx[2][0]/=sum;
			}
		}
		else if (data.a[mm.at(i)].type==3||data.a[mm.at(i)].type==2) {
			Matrix a(3,1),aa(3,1);
			if (data.a[mm.at(P)].norder>2) {
				for (j=0;j<data.a[mm.at(P)].norder;j++) {
					if (j!=memo) break;
				}
				n=0;
				a[0][0]=vec[i][0].Mtrx[1][0]*vec[P][j].Mtrx[2][0]-vec[i][0].Mtrx[2][0]*vec[P][j].Mtrx[1][0];
				a[1][0]=vec[i][0].Mtrx[2][0]*vec[P][j].Mtrx[1][0]-vec[i][0].Mtrx[1][0]*vec[P][j].Mtrx[2][0];
				a[2][0]=vec[i][0].Mtrx[0][0]*vec[P][j].Mtrx[1][0]-vec[i][0].Mtrx[1][0]*vec[P][j].Mtrx[0][0];
				if (fabs(a[0][0])<1E-3 && fabs(a[1][0])<1E-3 && fabs(a[2][0])<1E-3) n=1;
				if (n) {
					while (1) {
						for (k=0;k<3;k++) aa[k][0]=(double)rand()/RAND_MAX;
						a[0][0]=aa[1][0]*vec[P][j].Mtrx[2][0]-aa[2][0]*vec[P][j].Mtrx[1][0];
						a[1][0]=aa[2][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[2][0];
						a[2][0]=aa[0][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[0][0];
						if (fabs(a[0][0])>1E-3 || fabs(a[1][0])>1E-3 || fabs(a[2][0])>1E-3) break;
					}
				}
				for (k=1;k<data.a[mm.at(i)].norder;k++) {
					vec[i][k].Mtrx[0][0]=vec[i][0].Mtrx[0][0];
					vec[i][k].Mtrx[1][0]=vec[i][0].Mtrx[1][0];
					vec[i][k].Mtrx[2][0]=vec[i][0].Mtrx[2][0];
				}
				for (k=1;k<data.a[mm.at(i)].norder;k++) {
					vec[i][k].rotate(a,vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
					sum=0;
					sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
					sum=sqrt(sum);
					vec[i][k].Mtrx[0][0]/=sum;
					vec[i][k].Mtrx[1][0]/=sum;
					vec[i][k].Mtrx[2][0]/=sum;
				}
			}
			else {
				while (1) {
					for (k=0;k<3;k++) aa[k][0]=(double)rand()/RAND_MAX;
					a[0][0]=aa[1][0]*vec[P][j].Mtrx[2][0]-aa[2][0]*vec[P][j].Mtrx[1][0];
					a[1][0]=aa[2][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[2][0];
					a[2][0]=aa[0][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[0][0];
					if (fabs(a[0][0])>1E-3 || fabs(a[1][0])>1E-3 || fabs(a[2][0])>1E-3) break;
				}
				for (k=1;k<data.a[mm.at(i)].norder;k++) {
					vec[i][k].Mtrx[0][0]=vec[i][0].Mtrx[0][0];
					vec[i][k].Mtrx[1][0]=vec[i][0].Mtrx[1][0];
					vec[i][k].Mtrx[2][0]=vec[i][0].Mtrx[2][0];
				}
				for (k=1;k<data.a[mm.at(i)].norder;k++) {
					vec[i][k].rotate(a,vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
					sum=0;
					sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
					sum=sqrt(sum);
					vec[i][k].Mtrx[0][0]/=sum;
					vec[i][k].Mtrx[1][0]/=sum;
					vec[i][k].Mtrx[2][0]/=sum;
				}
			}
		}
		else if (data.a[mm.at(i)].type==4) {
			Matrix a(3,1);
			Matrix aa(3,1);
			while (1) {
				for (k=0;k<3;k++) aa[k][0]=(double)rand()/RAND_MAX;
				a[0][0]=aa[1][0]*vec[P][j].Mtrx[2][0]-aa[2][0]*vec[P][j].Mtrx[1][0];
				a[1][0]=aa[2][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[2][0];
				a[2][0]=aa[0][0]*vec[P][j].Mtrx[1][0]-aa[1][0]*vec[P][j].Mtrx[0][0];
				if (fabs(a[0][0])>1E-3 || fabs(a[1][0])>1E-3 || fabs(a[2][0])>1E-3) break;
			}
			for (k=1;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].Mtrx[0][0]=vec[i][0].Mtrx[0][0];
				vec[i][k].Mtrx[1][0]=vec[i][0].Mtrx[1][0];
				vec[i][k].Mtrx[2][0]=vec[i][0].Mtrx[2][0];
			}
			for (k=1;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].rotate(a,vec[i][k],2.0000*3.14159/3.0000*(double)pow(-1,k));
				sum=0;
				sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
				sum=sqrt(sum);
				vec[i][k].Mtrx[0][0]/=sum;
				vec[i][k].Mtrx[1][0]/=sum;
				vec[i][k].Mtrx[2][0]/=sum;
			}
		}
		else if (data.a[mm.at(i)].type==5) {
			for (k=0;k<3;k++) pt[k]=x[k].at(i);
			for (k=0;k<2;k++) xx[k]=pt[k]+(double)rand()/RAND_MAX*3.0;
			sum=0;
			sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[1][0]*(pt[1]-xx[1]);
			if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[2]=sum/vec[i][0].Mtrx[2][0]+pt[2];
			else {
				xx[2]=pt[2];
				if (fabs(vec[i][0].Mtrx[1][0])>1E-5) xx[1]=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])/vec[i][0].Mtrx[1][0]+pt[1];
				else xx[1]=pt[1];
				if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[0]=vec[i][0].Mtrx[1][0]*(pt[1]-xx[1])/vec[i][0].Mtrx[0][0]+pt[0];
				else xx[0]=pt[0];
			}
			sum=0;
			sum=(xx[0]-pt[0])*(xx[0]-pt[0])+(xx[1]-pt[1])*(xx[1]-pt[1])+(xx[2]-pt[2])*(xx[2]-pt[2]);
			sum=sqrt(sum);
			xx[0]=(xx[0]-pt[0])/sum;
			xx[1]=(xx[1]-pt[1])/sum;
			xx[2]=(xx[2]-pt[2])/sum;
			for (k=1;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].Mtrx[0][0]=xx[0];
				vec[i][k].Mtrx[1][0]=xx[1];
				vec[i][k].Mtrx[2][0]=xx[2];
			}
			for (k=2;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].rotate(vec[i][0],vec[i][k],2.0000*3.1415926/4.00000*((double)k-1.00000));
				sum=0;
				sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
				sum=sqrt(sum);
				vec[i][k].Mtrx[0][0]/=sum;
				vec[i][k].Mtrx[1][0]/=sum;
				vec[i][k].Mtrx[2][0]/=sum;
			}
		}
		else if (data.a[mm.at(i)].type==6) {
			for (k=0;k<3;k++) pt[k]=x[k].at(i);
			for (k=0;k<2;k++) xx[k]=pt[k]+(double)rand()/RAND_MAX*3.0;
			sum=0;
			sum=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])+vec[i][0].Mtrx[1][0]*(pt[1]-xx[1]);
			if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[2]=sum/vec[i][0].Mtrx[2][0]+pt[2];
			else {
				xx[2]=pt[2];
				if (fabs(vec[i][0].Mtrx[1][0])>1E-5) xx[1]=vec[i][0].Mtrx[0][0]*(pt[0]-xx[0])/vec[i][0].Mtrx[1][0]+pt[1];
				else xx[1]=pt[1];
				if (fabs(vec[i][0].Mtrx[2][0])>1E-5) xx[0]=vec[i][0].Mtrx[1][0]*(pt[1]-xx[1])/vec[i][0].Mtrx[0][0]+pt[0];
				else xx[0]=pt[0];
			}
			sum=0;
			sum=(xx[0]-pt[0])*(xx[0]-pt[0])+(xx[1]-pt[1])*(xx[1]-pt[1])+(xx[2]-pt[2])*(xx[2]-pt[2]);
			sum=sqrt(sum);
			xx[0]=(xx[0]-pt[0])/sum;
			xx[1]=(xx[1]-pt[1])/sum;
			xx[2]=(xx[2]-pt[2])/sum;
			for (k=1;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].Mtrx[0][0]=xx[0];
				vec[i][k].Mtrx[1][0]=xx[1];
				vec[i][k].Mtrx[2][0]=xx[2];
			}
			for (k=2;k<data.a[mm.at(i)].norder;k++) {
				vec[i][k].rotate(vec[i][0],vec[i][k],108.5*3.1415926/180.0*((double)k-1.000));
				sum=0;
				sum=vec[i][k].Mtrx[0][0]*vec[i][k].Mtrx[0][0]+vec[i][k].Mtrx[1][0]*vec[i][k].Mtrx[1][0]+vec[i][k].Mtrx[2][0]*vec[i][k].Mtrx[2][0];
				sum=sqrt(sum);
				vec[i][k].Mtrx[0][0]/=sum;
				vec[i][k].Mtrx[1][0]/=sum;
				vec[i][k].Mtrx[2][0]/=sum;
			}
		}
		//else if (data.a[mm.at(i)].type==7) {}	
	}
	int q=0;
	for (i=0;i<mm.size();i++) {
		if (Cyindex.at(i)) {
			for (j=0;j<data.a[mm.at(i)].norder;j++) {
				if (cont[i][j]) {
					for (k=i+1;k<mm.size();k++) {
						if (Cyindex.at(i)==Cyindex.at(k)) {
							for (n=0;n<data.a[mm.at(k)].norder;n++) {
								if (cont[k][n]) {
									sum=0.0;
									sum=fabs(vec[i][j].Mtrx[0][0]+vec[k][n].Mtrx[0][0])+fabs(vec[i][j].Mtrx[1][0]+vec[k][n].Mtrx[1][0])+fabs(vec[i][j].Mtrx[2][0]+vec[k][n].Mtrx[2][0]);
									if (sum<1E-3) {
										cont[i][j]=cont[k][n]=0;
										q=1;
										break;
									}
								}
							}
						}
					}
					if (q) break;
				}
			}
			if (q) break;
		}
		if (q) break;
	}
	if (!q) {
		for (i=0;i<mm.size();i++) {
			if (Cyindex.at(i)) {
				for (j=0;j<data.a[mm.at(i)].norder;j++) {
					if (cont[i][j]) {
						cont[i][j]=0;
						break;
					}
				}
			}
		}
	}
	for (i=0;i<mm.size();i++) {
		for (j=0;j<data.a[mm.at(i)].norder;j++) {
			if (cont[i][j]) {
				for (k=0;k<3;k++) {
					x[k].push_back(x[k].at(i)+vec[i][j].Mtrx[k][0]*(data.a[mm.at(i)].rb+0.32));
				}
				cont[i][j]=0;
			}
		}
	}
	opt[0].count=if_circle;
	opt[0].num=x[0].size();
	opt[0].f.resize(opt[0].num*3,1);
	opt[0].df.resize(1,opt[0].num*3);
	opt[0].x0.resize(opt[0].num*3,1);
	opt[0].id=new int [opt[0].num];
	opt[0].table=new int *[opt[0].num];
	for (i=0;i<opt[0].num;i++) opt[0].table[i]=new int [opt[0].num];	

	for (i=0;i<opt[0].num;i++) {
		if (i<mm.size()) opt[0].id[i]=mm.at(i);
		else opt[0].id[i]=0;
	}
	for (i=0;i<opt[0].num;i++) {
		for (j=0;j<opt[0].num;j++) opt[0].table[i][j]=0;
	}
	for (i=1;i<mm.size();i++) {
		opt[0].table[i][Pindex.at(i)-1]=opt[0].table[Pindex.at(i)-1][i]=1;
	}
	if (if_circle) {
		for (i=0;i<Cyindex.size();i++) {
			if (Cyindex.at(i)) {
				for (j=i+1;j<Cyindex.size();j++) {
					if (Cyindex.at(i)==Cyindex.at(j)) {
						opt[0].table[i][j]=opt[0].table[j][i]=1;
					}
				}
			}
		}
	}
	int t=0;
	for (i=0;i<mm.size();i++) {
		k=0;
		for (j=0;j<opt[0].num;j++) k+=opt[0].table[i][j];
		if (k<data.a[mm.at(i)].norder) {
			for (n=0;n<(data.a[mm.at(i)].norder-k);n++) {
				opt[0].table[i][mm.size()+t]=opt[0].table[mm.size()+t][i]=1;
				t++;
			}
		}
	}
	for (i=0;i<opt[0].num;i++) {
		for (j=0;j<3;j++) {
			opt[0].x0[i+j*opt[0].num][0]=x[j].at(i);
		}
	}	
	opt[0].optimize();
	opt[0].output(outs);
	outs<<endl;
	delete [] opt[0].id;
	for (i=0;i<opt[0].num;i++) {
		delete [] opt[0].table[i];
	}
	delete [] opt[0].table;
	for (i=0;i<Mindex.size();i++) delete [] cont[i];
	delete [] cont;
	delete [] pt;
	delete [] xx;
	pt=NULL;
	xx=NULL;
	cont=NULL; 
	return;
}


int MOLECULE::print() {
	int i;
	//if (para.protect) prct();

	cout << setfill(' ');
    cout<<"natom "<<setw(4)<<left<<setfill(' ')<<Cindex.size();
    cout<<endl;
	cout<<"Pindex  ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Pindex.at(i)<<" ";
	cout<<endl;
	cout<<"Cindex  ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Cindex.at(i)<<" ";
	cout<<endl;
	cout<<"Rindex  ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Rindex.at(i)<<" ";
	cout<<endl;
	cout<<"Mindex  ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Mindex.at(i)<<" ";
	cout<<endl;
	cout<<"Cyindex ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Cyindex.at(i)<<" ";
	cout<<endl;
	cout<<"if_circle "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;
	if (para.protect) {
		cout<<"protect  ";
		for (i=0;i<protect.size();i++) cout<<setw(4)<<left<<protect.at(i)<<" ";
		cout<<endl;
	}
	return 1;
}		

int MOLECULE::print(ofstream &ouf) {
	int i;
	//if (para.protect) prct();

	ouf << setfill(' ');
    ouf<<"natom  "<<setw(4)<<left<<setfill(' ')<<Cindex.size()<<" ";
    ouf<<endl;
	ouf<<"Pindex  ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Pindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"Cindex  ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Cindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"Rindex  ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Rindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"Mindex  ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Mindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"Cyindex  ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Cyindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"if_circle "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;
	if (para.protect) {
		ouf<<"protect  ";
		for (i=0;i<protect.size();i++) ouf<<setw(4)<<left<<protect.at(i)<<" ";
		ouf<<endl;
	}
	return 1;
}

int MOLECULE::prct() {
	prct(0,Cindex.size());
	return 1;
}

int MOLECULE::prct(int a, int b) {

	//long double time1=time(NULL);

	if (a>Cindex.size() || b>Cindex.size() || a<0 || b<0 || Cindex.size()==0) return 0;
	int i=0,j=0;
	int count=0;
	if (para.protect) {
		for (i=0;i<protect.size();i++) {
			if (protect.at(i)) count++;
		}
	}
	bool go=1;
	if (count/(double)Cindex.size()>0.5) go=0;
	if (para.protect) { // para.protect && go
		if (0) { // double bonds
			protect.resize(0);
			protect.resize(Cindex.size(),0); //20190801

			for (i=a;i<b;i++) {
				if (Rindex.at(i)==2) { //e.g. CC=CC
					protect.at(i)=1; // i = C3
					for (int k1=a;k1<b;k1++) {
						if (Pindex.at(i)==Cindex.at(k1)) {
							protect.at(k1)=1; // k1 = C2
							for (int k2=a;k2<b;k2++) {
								if (Pindex.at(k1)==Cindex.at(k2)) protect.at(k2)=1;  // k2 = C1
							}
						}
						if (Pindex.at(k1)==Cindex.at(i)) protect.at(k1)=1;; // k1 = C4
					}
				}
			}
		}

		if (0) { // whole molecule
			protect.resize(0);
			protect.resize(Cindex.size(),0); //20190801
			for (i=a;i<b;i++) protect.at(i)=1;
		}

		if (1 && Cindex.size()>0) { // charged atom
			protect.resize(0);
			protect.resize(Cindex.size(),0); //20190801
			//for (i=a;i<b;i++) protect.at(i)=0;

			for (i=0;i<Cindex.size();i++) {
				if (data.a[Mindex.at(i)].chg) {
					protect.at(i)=1;
					for (j=a;j<b;j++) {
						//if ( abs(Cindex.at(i)-Cindex.at(j))<=1 ) protect.at(j)=1;
						if (Cindex.at(j)==Pindex.at(i) || Cindex.at(i)==Pindex.at(j)) protect.at(j)=1;
					}
				}
			}
		}

		if (1 && Cindex.size()>0) { // ring
			if (0) {
				protect.resize(0);
				protect.resize(Cindex.size(),0); //20190801
			}

			if (0) {
				for (i=a;i<b;i++) {
					int a=0,pos1=-1,pos2=-2;
					if (Cyindex.at(i)>0) {
						a=Cyindex.at(i);
						pos1=i;
						for (j=i+1;j<b;j++) {
							if (Cyindex.at(j)==a) {
								pos2=j;
								break;
							}
						}
						if (pos2>pos1) {
							for (int k=pos1;k<=pos2;k++) protect.at(k)=1;
						}
					}
				}
			}

			if (1) {
				vector<int> C_ringmember(0);
				for (int k1=if_circle;k1>=1;k1--) {
					C_ringmember.resize(0);
					int rep_quota=2;
					for (int k2=Cyindex.size()-1;k2>=0;k2--) {
						if (Cyindex.at(k2)) {
							int digits=(int)log10(Cyindex.at(k2))+1;
							for (int k3=1;k3<=digits;k3++) {
								int w=Cyindex.at(k2)/(int)pow(10,k3-1);
								w=w%10;
								if (w==k1) {
									int compa=Cindex.at(k2);
									int cont=0;
									while (compa>=1 && cont<=Cyindex.size()) {
										if (1) C_ringmember.push_back(compa);

										if (compa>=1) {
											if (Pindex.at(compa-1)>=1) compa=Cindex.at(Pindex.at(compa-1)-1);
											else break;
											cont++;
										}
										else break;
									} 
								}
							}
						}
					}
					if (C_ringmember.size()>0) {
						for (int k2=0;k2<C_ringmember.size()-1;k2++) {
							for (int k3=k2+1;k3<C_ringmember.size();k3++) {
								if (C_ringmember.at(k2)==C_ringmember.at(k3) && rep_quota<=0) {
									C_ringmember.at(k2)=C_ringmember.at(k3)=-1;
								}
								else if (C_ringmember.at(k2)==C_ringmember.at(k3) && rep_quota>0) rep_quota--;
							}
						}
						for (int k2=0;k2<C_ringmember.size();k2++) {
							if (C_ringmember.at(k2)!=-1) protect.at(C_ringmember.at(k2)-1)=1;
						}
					}
				}
				vector<int>().swap(C_ringmember);
			}

		}
	}

    //time1=time(NULL)-time1;
    //cout << "PROTECTION FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}

int MOLECULE::fcngroup2mds() {
	if (Mindex.size()>0) {
		bool if_fcn=1;
		if (1) {
			for (int i=0;i<Mindex.size();i++) {
				switch (Mindex.at(i)) {
					case 65:
					case 63:
					case 62:
					case 61:
					case 60:
					case 58:
					case 56:
					case 55:
					case 54:
					case 53:
					case 52:
					case 51:
					case 50:
					case 49:
					case 48:
					case 47:
					case 46:
					case 45:
					case 43:
					case 42:
					case 41:
					case 40:
					case 39:
					case 38:
					case 37:
					case 36:
					case 33:
						if_fcn=1;
						break;
					default:
						if_fcn=0;
						break;
				}
				if (if_fcn) break;
			}
			if (if_circle) if_fcn=1;
		}
		if (if_fcn) {
			smiles=molesmi;
			//clear();
			empty();
			clean();
			input();
			mds2smi();
			smiles=molesmi;
			if (para.protect) prct();
		}
	}
	else return 0;
	return 1;
}

int MOLECULE::decyc_small_ring(int size) {
    vector<int> C_ringmember(0);
	
    for (int k1=if_circle;k1>=1;k1--) {
        C_ringmember.resize(0);
        int rep_quota=1;
        for (int k2=Cyindex.size()-1;k2>=0;k2--) {
            if (Cyindex.at(k2)) {
                int digits=(int)log10(Cyindex.at(k2))+1;
                for (int k3=1;k3<=digits;k3++) {
                    int w=Cyindex.at(k2)/(int)pow(10,k3-1);
                    w=w%10;
                    if (w==k1) {
                    	int compa=Cindex.at(k2);
                        int cont=0;
                        while (compa>1 && cont<=Cyindex.size()) {
                            if (1) C_ringmember.push_back(compa);

                            if (compa>=1) {
                                if (Pindex.at(compa-1)>=1) compa=Cindex.at(Pindex.at(compa-1)-1);
                                else break;
                                cont++;
                            }
                            else break;
                        }
						if (compa==1) if (1) C_ringmember.push_back(compa);
                    }
                }
        	}
        }
        if (C_ringmember.size()>0) {
            for (int k2=0;k2<C_ringmember.size()-1;k2++) {
                for (int k3=k2+1;k3<C_ringmember.size();k3++) {
					if (C_ringmember.at(k2)!=-1 && C_ringmember.at(k3)!=-1) {
                    	if (C_ringmember.at(k2)==C_ringmember.at(k3) && rep_quota<=0) {
                        	C_ringmember.at(k2)=C_ringmember.at(k3)=-1;
                    	}
                    	else if (C_ringmember.at(k2)==C_ringmember.at(k3) && rep_quota>0) {
							C_ringmember.at(k3)=-1;
							rep_quota--;
						}
					}
                }
            }
        }
		int num_rmember=0;
		if (C_ringmember.size()>0) {
        	for (int k2=0;k2<C_ringmember.size();k2++) {
            	if (C_ringmember.at(k2)!=-1) num_rmember++;
        	}
			if (num_rmember<size) {
				for (int k2=0;k2<Cyindex.size();k2++) {
					if (Cyindex.at(k2)) {
						int digits=(int)log10(Cyindex.at(k2))+1;
                		for (int k3=1;k3<=digits;k3++) {
                    		int w=Cyindex.at(k2)/(int)pow(10,k3-1);
                    		w=w%10;
                    		if (w==k1) {
								Cyindex.at(k2)=Cyindex.at(k2)%(int)pow(10,k3-1)+(Cyindex.at(k2)/(int)pow(10,k3))*(int)pow(10,k3-1);						
							}
						}
					}
				}
				if_circle--;
			}
		}
    }
    vector<int>().swap(C_ringmember);
	

	return 1;

}


int MOLECULE::del_unpaired_ring_no(){
	int k;
    vector<int> freq(if_circle,0);
	if (if_circle) {
	    for (k=0;k<Cyindex.size();k++) {  // 20200131 prevent unpaired ring number
	        if (Cyindex.at(k)) {
	            int digits=(int)log10(Cyindex.at(k))+1;
	            for (int ap=1;ap<=digits;ap++) {
	                int w=Cyindex.at(k)/(int)pow(10,ap-1);
	                w=w%10;
	                if (w<=if_circle) freq.at(w-1)++;
					else Cyindex.at(k)=Cyindex.at(k)%(int)pow(10,ap-1)+(Cyindex.at(k)/(int)pow(10,ap))*(int)pow(10,ap-1);
	            }
	        }
	    }
	    for (k=freq.size()-1;k>=0;k--) {
	        if (freq.at(k)%2) {
	            for (int k1=0;k1<Cyindex.size();k1++) {
	                if (Cyindex.at(k1)==k+1) Cyindex.at(k1)=0;
	            }
				if_circle--;
	        }
	    }
	}
	vector<int>().swap(freq);

	return 1;
}


int MOLECULE::canonicalize_SMILES() {
	if (0) { //system code for obabel
		ofstream out((para.smidir+"tmp.smi").c_str());
		out << molesmi << endl;
		out.close();
		// /home/akitainu/bin/openbabel-install/bin/obabel
		system(("ssh cluster \" /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel  -ismi "+para.smidir+"tmp.smi -osmi -O "+para.smidir+"tmp.smi --gen3D --canonical \"").c_str());
		//system(("  /home/akitainu/bin/openbabel-install/bin/obabel  -ismi "+para.smidir+"tmp.smi  -osmi -O "+para.smidir+"tmp.smi  --gen3D  --canonical  > /dev/null ").c_str());
	
		molesmi="null";
		ifstream inf((para.smidir+"tmp.smi").c_str());
		inf >> molesmi >> ws;
		inf.close();
		smiles=molesmi;
	}
    
	if (1) {
		stringstream ss(molesmi);
    	stringstream ss_out("");

    	OBConversion conv(&ss, &ss_out);
    	if(conv.SetInAndOutFormats("SMI","SMI"))
    	{
        	//conv.AddOption("h", OBConversion::GENOPTIONS);
        	conv.AddOption("gen3D", OBConversion::GENOPTIONS);
        	conv.AddOption("canonical", OBConversion::GENOPTIONS);
        	conv.Convert();
    	}
		//smiles=molesmi=ss_out.str();
		ss_out >> molesmi >> ws;
		smiles=molesmi;
	}

	return 1;
}


