#include "MOLECULE.h"
#define p_circle 0.2
using namespace std;
using namespace OpenBabel;


extern PARAMETER para;


double MOLECULE::prob() {
	double p=(double)rand()/(double)RAND_MAX;
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
		//clean();
		clear();

		ctsisomer.resize(2,vector<string>(0));
		//ctsisomer.resize(mol.Cindex.size(),vector<string>(0));
		Cyindex.resize(mol.Cindex.size(),vector<int>(0));
		for (i=0;i<mol.Cindex.size();i++) {
			Cindex.push_back(mol.Cindex.at(i));
			Pindex.push_back(mol.Pindex.at(i));
			Rindex.push_back(mol.Rindex.at(i));
			//Cyindex.push_back(mol.Cyindex.at(i));
			for (j=0;j<mol.Cyindex.at(i).size();j++) Cyindex.at(i).push_back(mol.Cyindex.at(i).at(j));
			Mindex.push_back(mol.Mindex.at(i));
			if (para.protect) protect.push_back(mol.protect.at(i));
            //ctsisomer.push_back(mol.ctsisomer.at(i));
            ctsisomer.at(0).push_back(mol.ctsisomer.at(0).at(i));
            ctsisomer.at(1).push_back(mol.ctsisomer.at(1).at(i));
			//for (j=0;j<mol.ctsisomer.at(i).size();j++) ctsisomer.at(i).push_back(<mol.ctsisomer.at(i).at(j));
		}
		for (i=0;i<mol.Cybnd.size();i++) Cybnd.push_back(mol.Cybnd.at(i));
		if_circle=mol.if_circle;

		lnKow=mol.lnKow;
		Henry=mol.Henry;

		HOMO=mol.HOMO;
		LUMO=mol.LUMO;
		Egap=mol.Egap;

		IP=mol.IP;
		EA=mol.EA;
		hardness=mol.hardness;
		ephilic=mol.ephilic;

		fitness=mol.fitness;
		p_sel=mol.p_sel;
		front=mol.front;
		crowdist=mol.crowdist;
		natom=mol.natom;
		smiles=mol.smiles;
		molesmi=mol.molesmi;
		nsubcomp=1;
		comp_id=mol.comp_id;
		frac=mol.frac;
		chg=mol.chg;
		multiplicity=mol.multiplicity;
		
		data=mol.data;

		reset();
		//mds2smi();
	}
	return 1;
}

int MOLECULE::crossover(MOLECULE &aaa,int pp, int jj,bool cistrans) {
	if (Cindex.size()<2) return 0;
	if (aaa.Cindex.size()<2) return 0;
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
	//int delta;
	//int pos_i;
	//int pos_j;
	int lenA=0;
	int lenB=0;
	vector<int> t,f;
	vector<int> ci_ref, mi_ref,pi_ref,ri_ref,pri,cybndi; //,cyi
	vector<int> cj_ref, mj_ref,pj_ref,rj_ref,prj,cybndj; //,cyj
	vector< vector<int> > cyi(0,vector<int>(0));
	vector< vector<int> > cyj(0,vector<int>(0));
	//vector<string> cti_ref,ctj_ref;
	vector< vector<string> > cti_ref,ctj_ref;
	
	vector<bool> ifchg(4,0);
	if (para.ion && 0) {
		for (int k1=0;k1<Cindex.size();k1++) {
			if (Cindex.at(k1)<=Cindex.at(pp) && data->a[Mindex.at(k1)].chg) {
				ifchg.at(0)=1;
				break;
			}
		}
		for (int k1=0;k1<Cindex.size();k1++) {
			if (Cindex.at(k1)>=Cindex.at(pp) && data->a[Mindex.at(k1)].chg) {
				ifchg.at(1)=1;
				break;
			}
		}
		for (int k1=0;k1<aaa.Cindex.size();k1++) {
			if (aaa.Cindex.at(k1)<=aaa.Cindex.at(jj) && data->a[aaa.Mindex.at(k1)].chg) {
				ifchg.at(2)=1;
				break;
			}
		}
		for (int k1=0;k1<aaa.Cindex.size();k1++) {
			if (aaa.Cindex.at(k1)>=aaa.Cindex.at(jj) && data->a[aaa.Mindex.at(k1)].chg) {
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
    if (para.protect && if_circle && 0) {
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Cindex.at(k1)<=Cindex.at(pp) && Cyindex.at(k1).size()) { // && Cyindex.at(k1)
                ifcyc.at(0)=1;
                break;
            }
        }
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Cindex.at(k1)>=Cindex.at(pp) && Cyindex.at(k1).size()) { // && Cyindex.at(k1)
                ifcyc.at(1)=1;
                break;
            }
        }
        for (int k1=0;k1<aaa.Cindex.size();k1++) {
            if (aaa.Cindex.at(k1)<=aaa.Cindex.at(jj) && aaa.Cyindex.at(k1).size()) { // && aaa.Cyindex.a(k1)
                ifcyc.at(2)=1;
                break;
            }
        }
        for (int k1=0;k1<aaa.Cindex.size();k1++) {
            if (aaa.Cindex.at(k1)>=aaa.Cindex.at(jj) && aaa.Cyindex.at(k1).size()) { // && aaa.Cyindex.at(k1)
                ifcyc.at(3)=1;
                break;
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
	}

	//double time1=time(NULL);

	k=0;
	if (aaa.Rindex.at(jj)==Rindex.at(pp)) {
		k=1;
		m=jj;
		n=pp;
		//b=aaa.Pindex.at(m);
		//a=Pindex.at(n);
		//y=aaa.Cindex.at(m);
		//x=Cindex.at(n);
		ri_ref.push_back(Rindex.at(n));
		rj_ref.push_back(aaa.Rindex.at(m));
		ci_ref.push_back(Cindex.at(n));
		cj_ref.push_back(aaa.Cindex.at(m));
		pi_ref.push_back(Pindex.at(n));
		pj_ref.push_back(aaa.Pindex.at(m));
		mi_ref.push_back(Mindex.at(n));
		mj_ref.push_back(aaa.Mindex.at(m));
        //cti_ref.push_back(ctsisomer.at(n));

        cti_ref.resize(2,vector<string>(0));
        cti_ref.at(0).push_back(ctsisomer.at(0).at(n));
        cti_ref.at(1).push_back(ctsisomer.at(1).at(n));
		//cti_ref.resize(1,vector<string>(0));
		//for (int k1=0;k1<ctsisomer.at(n).size();k1++) cti_ref.at(0).push_back(ctsisomer.at(n).at(k1));

        //ctj_ref.push_back(aaa.ctsisomer.at(m));

        ctj_ref.resize(2,vector<string>(0));
        ctj_ref.at(0).push_back(aaa.ctsisomer.at(0).at(m));
        ctj_ref.at(1).push_back(aaa.ctsisomer.at(1).at(m));
        //ctj_ref.resize(1,vector<string>(0));
        //for (int k1=0;k1<aaa.ctsisomer.at(m).size();k1++) ctj_ref.at(0).push_back(aaa.ctsisomer.at(m).at(k1));

		//cyi.push_back(Cyindex.at(n));
		cyi.resize(1,vector<int>(0));
		for (int k1=0;k1<Cyindex.at(n).size();k1++) cyi.at(0).push_back(Cyindex.at(n).at(k1));

		//cyj.push_back(aaa.Cyindex.at(m));
		cyj.resize(1,vector<int>(0));
		for (int k1=0;k1<aaa.Cyindex.at(m).size();k1++) cyj.at(0).push_back(aaa.Cyindex.at(m).at(k1));

		if (para.protect) pri.push_back(protect.at(n));
		if (para.protect) prj.push_back(aaa.protect.at(m));
	}

	if (k==0) return 0;
	if (k==1) {
		lenA=Cindex.size();
		lenB=aaa.Cindex.size();
		/*
		for (n=1;n<Cindex.size();n++) {
			if (Cindex.at(n) > ci_ref.at(0)) {
				ref = ci_ref.size();
				for (m=0;m<ref;m++) {
					if (Pindex.at(n) == ci_ref.at(m)) {
						ci_ref.push_back(Cindex.at(n));
						mi_ref.push_back(Mindex.at(n));
						pi_ref.push_back(Pindex.at(n));
						ri_ref.push_back(Rindex.at(n));
                        //cti_ref.push_back(ctsisomer.at(n));
                        cti_ref.at(0).push_back(ctsisomer.at(0).at(n));
                        cti_ref.at(1).push_back(ctsisomer.at(1).at(n));
						//cti_ref.resize(cti_ref.size()+1,vector<string>(0));
						//for (int k1=0;k1<ctsisomer.at(n).size();k1++) cti_ref.at(cti_ref.size()-1).push_back(ctsisomer.at(n).at(k1));

						//cyi.push_back(Cyindex.at(n));
						cyi.resize(cyi.size()+1,vector<int>(0));
						for (int k1=0;k1<Cyindex.at(n).size();k1++) cyi.at(cyi.size()-1).push_back(Cyindex.at(n).at(k1));
						
						if (para.protect) pri.push_back(protect.at(n));
					}
				}
			}
		}
		*/

	    vector<int> curatm(0);

	    bool goout=0;
    	curatm.push_back(ci_ref.at(0));

    	do {
        	vector<int> tmpp(0);
        	for (int k3=0;k3<curatm.size();k3++) {
            	for (int k2=0;k2<Cindex.size();k2++) {
                	if (Pindex.at(k2)==curatm.at(k3)) {
	                    ci_ref.push_back(Cindex.at(k2));
    	                mi_ref.push_back(Mindex.at(k2));
        	            pi_ref.push_back(Pindex.at(k2));
            	        ri_ref.push_back(Rindex.at(k2));
                    	cti_ref.at(0).push_back(ctsisomer.at(0).at(k2));
                        cti_ref.at(1).push_back(ctsisomer.at(1).at(k2));

                        cyi.resize(cyi.size()+1,vector<int>(0));
						for (int k1=0;k1<Cyindex.at(k2).size();k1++) cyi.at(cyi.size()-1).push_back(Cyindex.at(k2).at(k1));

                        if (para.protect) pri.push_back(protect.at(k2));

						tmpp.push_back(Cindex.at(k2));
                	}
            	}
        	}
        	if (tmpp.size()>0) {
            	curatm.resize(0);
            	for (int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
        	}
        	else {
            	goout=1;
            	break;
        	}

    	} while (!goout);

		for (int i=0;i<ci_ref.size();i++) {
			for (int j=i+1;j<ci_ref.size();j++) {
				if (ci_ref.at(i)>ci_ref.at(j)) {
					swap(ci_ref.at(i),ci_ref.at(j));
					swap(mi_ref.at(i),mi_ref.at(j));
					swap(pi_ref.at(i),pi_ref.at(j));
					swap(ri_ref.at(i),ri_ref.at(j));
					swap(cti_ref.at(0).at(i),cti_ref.at(0).at(j));
					swap(cti_ref.at(1).at(i),cti_ref.at(1).at(j));
					swap(cyi.at(i),cyi.at(j));
					if (para.protect) swap(pri.at(i),pri.at(j));
				}
			}
		}

		/*
		for (m=1;m<aaa.Cindex.size();m++) {
			if (aaa.Cindex.at(m) > cj_ref.at(0) ) {
				ref = cj_ref.size();
				for (n=0;n<ref;n++) {
					if (aaa.Pindex.at(m) == cj_ref.at(n)) {
						cj_ref.push_back(aaa.Cindex.at(m));
						mj_ref.push_back(aaa.Mindex.at(m));
						pj_ref.push_back(aaa.Pindex.at(m));
						rj_ref.push_back(aaa.Rindex.at(m));
                        //ctj_ref.push_back(aaa.ctsisomer.at(m));
                        ctj_ref.at(0).push_back(aaa.ctsisomer.at(0).at(m));
                        ctj_ref.at(1).push_back(aaa.ctsisomer.at(1).at(m));
                        //ctj_ref.resize(ctj_ref.size()+1,vector<string>(0));
                        //for (int k1=0;k1<aaa.ctsisomer.at(m).size();k1++) ctj_ref.at(ctj_ref.size()-1).push_back(aaa.ctsisomer.at(m).at(k1));

						//cyj.push_back(aaa.Cyindex.at(m));
                        cyj.resize(cyj.size()+1,vector<int>(0));
                        for (int k1=0;k1<aaa.Cyindex.at(m).size();k1++) cyj.at(cyj.size()-1).push_back(aaa.Cyindex.at(m).at(k1));
						
						if (para.protect) prj.push_back(aaa.protect.at(m));
					}
				}
			}
		}
		*/

        goout=0;
		curatm.resize(0);
        curatm.push_back(cj_ref.at(0));

        do {
            vector<int> tmpp(0);
            for (int k3=0;k3<curatm.size();k3++) {
                for (int k2=0;k2<aaa.Cindex.size();k2++) {
                    if (aaa.Pindex.at(k2)==curatm.at(k3)) {
                        cj_ref.push_back(aaa.Cindex.at(k2));
                        mj_ref.push_back(aaa.Mindex.at(k2));
                        pj_ref.push_back(aaa.Pindex.at(k2));
                        rj_ref.push_back(aaa.Rindex.at(k2));
                        ctj_ref.at(0).push_back(aaa.ctsisomer.at(0).at(k2));
                        ctj_ref.at(1).push_back(aaa.ctsisomer.at(1).at(k2));

                        cyj.resize(cyj.size()+1,vector<int>(0));
                        for (int k1=0;k1<aaa.Cyindex.at(k2).size();k1++) cyj.at(cyj.size()-1).push_back(aaa.Cyindex.at(k2).at(k1));

                        if (para.protect) prj.push_back(aaa.protect.at(k2));
						
						tmpp.push_back(aaa.Cindex.at(k2));
                    }
                }
            }
            if (tmpp.size()>0) {
                curatm.resize(0);
                for (int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
            }
            else {
                goout=1;
                break;
            }

        } while (!goout);

        for (int i=0;i<cj_ref.size();i++) {
            for (int j=i+1;j<cj_ref.size();j++) {
                if (cj_ref.at(i)>cj_ref.at(j)) {
                    swap(cj_ref.at(i),cj_ref.at(j));
                    swap(mj_ref.at(i),mj_ref.at(j));
                    swap(pj_ref.at(i),pj_ref.at(j));
                    swap(rj_ref.at(i),rj_ref.at(j));
                    swap(ctj_ref.at(0).at(i),ctj_ref.at(0).at(j));
                    swap(ctj_ref.at(1).at(i),ctj_ref.at(1).at(j));
                    swap(cyj.at(i),cyj.at(j));
                    if (para.protect) swap(prj.at(i),prj.at(j));
                }
            }
        }

		if (1) {
			for (n=0;n<pri.size();n++) if (pri.at(n)) return 0;
			for (n=0;n<prj.size();n++) if (prj.at(n)) return 0;
    		if (para.ion) {
				bool michg=0,mjchg=0;
    			for (int k1=0;k1<mi_ref.size();k1++) {
					if (data->a[mi_ref.at(k1)].chg) {
						michg=1;
						break;
					}
				}
            	for (int k1=0;k1<mj_ref.size();k1++) {
                	if (data->a[mj_ref.at(k1)].chg) {
                    	mjchg=1;
                    	break;
                	}
            	}
				if (michg!=mjchg) return 0;
			}
		}

		//if (0) {
		//a=0;
		//for (n=0;n<cyi.size();n++) {
		//	if (cyi.at(n)) {
		//		cyi.at(n)+=if_circle;
		//		a++;
		//	}
		//}
		//if (a%2 && if_circle) {
		//	for (n=0;n<cyi.size();n++) if (cyi.at(n)) cyi.at(n)=if_circle;
		//}
		//else if (a%2) {
		//	for (n=0;n<cyi.size();n++) cyi.at(n)=0;
		//}
		//a=0;
		//for (n=0;n<cyj.size();n++) {
		//	if (cyj.at(n)) {
		//		cyj.at(n)+=aaa.if_circle;
		//		a++;
		//	}
		//}
		//if (a%2 && aaa.if_circle) {
		//	for (n=0;n<cyj.size();n++) if (cyj.at(n)) cyj.at(n)=aaa.if_circle;
		//}
		//if (a%2) {
		//	for (n=0;n<cyj.size();n++) cyj.at(n)=0;
		//}
		//}

		if (1) {
			cybndi.resize(if_circle+aaa.if_circle,0);
			for (n=0;n<aaa.if_circle;n++) cybndi.at(n)=aaa.Cybnd.at(n);
		
			cybndj.resize(if_circle+aaa.if_circle,0);
			for (n=0;n<if_circle;n++) cybndj.at(n)=Cybnd.at(n);
		
			int tm1=if_circle;
			for (n=0;n<cyi.size();n++) {
				for (int k1=0;k1<cyi.at(n).size();k1++) {
					cyi.at(n).at(k1)+=aaa.if_circle;
				}
				//if (cyi.at(n)) {
                //	int digits=(int)log10(cyi.at(n).at(k1))+1;
            	//	for (int ap=1;ap<=digits;ap++) {
                //		cyi.at(n)+=aaa.if_circle*(int)pow(10,ap-1);
                //		int w=cyi.at(n)/(int)pow(10,ap-1);
                //		w=w%10;
                //		cybndi.at(w-1)=Cybnd.at(w-1-aaa.if_circle);
                //	}
				//}
			}
			for (n=0;n<cyj.size();n++) {
                for (int k1=0;k1<cyj.at(n).size();k1++) {
                    cyj.at(n).at(k1)+=if_circle;
                }
				//if (cyj.at(n)) {
                //	int digits=(int)log10(cyj.at(n))+1;
                //	for (int ap=1;ap<=digits;ap++) {
				//		cyj.at(n)+=if_circle*(int)pow(10,ap-1);
                //		int w=cyj.at(n)/(int)pow(10,ap-1);
                //		w=w%10;
				//		cybndj.at(w-1)=aaa.Cybnd.at(w-1-if_circle);
                //	}
				//}
			}
            if_circle+=aaa.if_circle;
            aaa.if_circle+=tm1;
		}

		//x=ci_ref.at(0);
		//y=cj_ref.at(0);
		a=pi_ref.at(0);
		b=pj_ref.at(0);
		/*
		for (n=0;n<ci_ref.size();n++) {
			for (m=0;m<Cindex.size();m++) {
				if (Cindex.at(m) == ci_ref.at(n)) {
					Mindex.erase(Mindex.begin()+m);
					Cindex.erase(Cindex.begin()+m);
					Pindex.erase(Pindex.begin()+m);
					Rindex.erase(Rindex.begin()+m);
					//Cyindex.erase(Cyindex.begin()+m);
					Cyindex.erase(Cyindex.begin()+m);
                    //ctsisomer.erase(ctsisomer.begin()+m);
                    ctsisomer.at(0).erase(ctsisomer.at(0).begin()+m);
                    ctsisomer.at(1).erase(ctsisomer.at(1).begin()+m);
					//ctsisomer.erase(ctsisomer.begin()+m);
					if (para.protect) protect.erase(protect.begin()+m);
				}
			}
		}
		*/

		for (n=ci_ref.size()-1;n>=0;n--) {
			Mindex.erase(Mindex.begin()+ci_ref.at(n)-1);
			Cindex.erase(Cindex.begin()+ci_ref.at(n)-1);
			Pindex.erase(Pindex.begin()+ci_ref.at(n)-1);
			Rindex.erase(Rindex.begin()+ci_ref.at(n)-1);
			Cyindex.erase(Cyindex.begin()+ci_ref.at(n)-1);
			ctsisomer.at(0).erase(ctsisomer.at(0).begin()+ci_ref.at(n)-1);
			ctsisomer.at(1).erase(ctsisomer.at(1).begin()+ci_ref.at(n)-1);
			if (para.protect) protect.erase(protect.begin()+ci_ref.at(n)-1);
		}


		/*
		for (m=0;m<cj_ref.size();m++) {
			for (n=0;n<aaa.Cindex.size();n++) {
				if (aaa.Cindex.at(n) == cj_ref.at(m)) {
					aaa.Mindex.erase(aaa.Mindex.begin()+n);
					aaa.Cindex.erase(aaa.Cindex.begin()+n);
					aaa.Pindex.erase(aaa.Pindex.begin()+n);
					aaa.Rindex.erase(aaa.Rindex.begin()+n);
					//aaa.Cyindex.erase(aaa.Cyindex.begin()+n);
					aaa.Cyindex.erase(aaa.Cyindex.begin()+n);
                    //aaa.ctsisomer.erase(aaa.ctsisomer.begin()+n);
                    aaa.ctsisomer.at(0).erase(aaa.ctsisomer.at(0).begin()+n);
                    aaa.ctsisomer.at(1).erase(aaa.ctsisomer.at(1).begin()+n);
					//aaa.Cyindex.erase(aaa.Cyindex.begin()+n);
					if (para.protect) aaa.protect.erase(aaa.protect.begin()+n);
				}
			}
		}
		*/

        for (n=cj_ref.size()-1;n>=0;n--) {
            aaa.Mindex.erase(aaa.Mindex.begin()+cj_ref.at(n)-1);
            aaa.Cindex.erase(aaa.Cindex.begin()+cj_ref.at(n)-1);
            aaa.Pindex.erase(aaa.Pindex.begin()+cj_ref.at(n)-1);
            aaa.Rindex.erase(aaa.Rindex.begin()+cj_ref.at(n)-1);
            aaa.Cyindex.erase(aaa.Cyindex.begin()+cj_ref.at(n)-1);
            aaa.ctsisomer.at(0).erase(aaa.ctsisomer.at(0).begin()+cj_ref.at(n)-1);
            aaa.ctsisomer.at(1).erase(aaa.ctsisomer.at(1).begin()+cj_ref.at(n)-1);
            if (para.protect) aaa.protect.erase(aaa.protect.begin()+cj_ref.at(n)-1);
        }

		/*
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
		*/

		for (n=0;n<ci_ref.size();n++) {
			ci_ref.at(n)+=(lenB+2);
			pi_ref.at(n)+=(lenB+2);
		}
        for (n=0;n<cj_ref.size();n++) {
            cj_ref.at(n)+=(lenA+2);
            pj_ref.at(n)+=(lenA+2);
        }


		pi_ref.at(0) = b;
		pj_ref.at(0) = a;
		//ci_ref.at(0) = y;
		//cj_ref.at(0) = x;

		/*
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
		*/
	}

	for (n=0;n<cj_ref.size();n++) {
        Cindex.push_back(cj_ref.at(n));
        Pindex.push_back(pj_ref.at(n));
        Mindex.push_back(mj_ref.at(n));
        Rindex.push_back(rj_ref.at(n));
        //Cyindex.push_back(cyj.at(n));
        Cyindex.resize(Cyindex.size()+1,vector<int>(0));
        for (int k1=0;k1<cyj.at(n).size();k1++) Cyindex.at(Cyindex.size()-1).push_back(cyj.at(n).at(k1));
        //ctsisomer.push_back(ctj_ref.at(n));
        ctsisomer.at(0).push_back(ctj_ref.at(0).at(n));
        ctsisomer.at(1).push_back(ctj_ref.at(1).at(n));
        //ctsisomer.resize(ctsisomer.size()+1,vector<string>(0));
        //for (int k1=0;k1<ctj_ref.at(n).size();k1++)  ctsisomer.at(ctsisomer.size()-1).push_back(ctj_ref.at(n).at(k1));
        if (para.protect) protect.push_back(prj.at(n));
	}

    for (n=0;n<ci_ref.size();n++) {
        aaa.Cindex.push_back(ci_ref.at(n));
        aaa.Pindex.push_back(pi_ref.at(n));
        aaa.Mindex.push_back(mi_ref.at(n));
        aaa.Rindex.push_back(ri_ref.at(n));
        //aaa.Cyindex.push_back(cyi.at(n));
        aaa.Cyindex.resize(aaa.Cyindex.size()+1,vector<int>(0));
        for (int k1=0;k1<cyi.at(n).size();k1++) aaa.Cyindex.at(aaa.Cyindex.size()-1).push_back(cyi.at(n).at(k1));
        //aaa.ctsisomer.push_back(cti_ref.at(n));
        aaa.ctsisomer.at(0).push_back(cti_ref.at(0).at(n));
        aaa.ctsisomer.at(1).push_back(cti_ref.at(1).at(n));
        //aaa.ctsisomer.resize(aaa.ctsisomer.size()+1,vector<string>(0));
        //for (int k1=0;k1<cti_ref.at(n).size();k1++)  aaa.ctsisomer.at(aaa.ctsisomer.size()-1).push_back(cti_ref.at(n).at(k1));
        if (para.protect) aaa.protect.push_back(pri.at(n));
    }

    for (int i=0;i<Cindex.size();i++) {
        int j=Cindex.at(i);
        Cindex.at(i)=(i+1);
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Pindex.at(k1)==j) Pindex.at(k1)=Cindex.at(i);
        }
    }

    for (int i=0;i<aaa.Cindex.size();i++) {
        int j=aaa.Cindex.at(i);
        aaa.Cindex.at(i)=(i+1);
        for (int k1=0;k1<aaa.Cindex.size();k1++) {
            if (aaa.Pindex.at(k1)==j) aaa.Pindex.at(k1)=aaa.Cindex.at(i);
        }
    }

	/*
	n=0;
	ref=1;
	while(n<cj_ref.size()) {
		for (m=0;m<Cindex.size();m++) {
			if (Cindex.at(m) > cj_ref.at(n) && m >= 1) { // Cindex[m] > cj_ref[n] && Cindex[m-1] < cj_ref[n]  orig // 20191130
				if (Cindex.at(m-1) < cj_ref.at(n)) {
					Cindex.insert(Cindex.begin() + m, cj_ref.at(n));
					Pindex.insert(Pindex.begin() + m, pj_ref.at(n));
					Mindex.insert(Mindex.begin() + m, mj_ref.at(n));
					Rindex.insert(Rindex.begin() + m, rj_ref.at(n));
					//Cyindex.insert(Cyindex.begin() + m, cyj.at(n));
					Cyindex.insert(Cyindex.begin() + m, cyj.at(n));
                    //ctsisomer.insert(ctsisomer.begin() + m, ctj_ref.at(n));
                    ctsisomer.at(0).insert(ctsisomer.at(0).begin() + m, ctj_ref.at(0).at(n));
                    ctsisomer.at(1).insert(ctsisomer.at(1).begin() + m, ctj_ref.at(1).at(n));
					//ctsisomer.insert(ctsisomer.begin() + m, ctj_ref.at(n));
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
				//Cyindex.push_back(cyj.at(n));
				Cyindex.resize(Cyindex.size()+1,vector<int>(0));
				for (int k1=0;k1<cyj.at(n).size();k1++) Cyindex.at(Cyindex.size()-1).push_back(cyj.at(n).at(k1));
                //ctsisomer.push_back(ctj_ref.at(n));
                ctsisomer.at(0).push_back(ctj_ref.at(0).at(n));
                ctsisomer.at(1).push_back(ctj_ref.at(1).at(n));
				//ctsisomer.resize(ctsisomer.size()+1,vector<string>(0));
				//for (int k1=0;k1<ctj_ref.at(n).size();k1++)  ctsisomer.at(ctsisomer.size()-1).push_back(ctj_ref.at(n).at(k1));
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
			if (aaa.Cindex.at(n) > ci_ref.at(m) && n >= 1) { // aaa.Cindex[n] > ci_ref[m] && aaa.Cindex[n-1] < ci_ref[m]  orig  //20191130
				if (aaa.Cindex.at(n-1) < ci_ref.at(m)) {
					ref=n;
					aaa.Cindex.insert(aaa.Cindex.begin() + ref, ci_ref.at(m));
					aaa.Pindex.insert(aaa.Pindex.begin() + ref, pi_ref.at(m));
					aaa.Mindex.insert(aaa.Mindex.begin() + ref, mi_ref.at(m));
					aaa.Rindex.insert(aaa.Rindex.begin() + ref, ri_ref.at(m));
					//aaa.Cyindex.insert(aaa.Cyindex.begin() + ref, cyi.at(m));
					aaa.Cyindex.insert(aaa.Cyindex.begin() + ref, cyi.at(m));
                    //aaa.ctsisomer.insert(aaa.ctsisomer.begin() + ref, cti_ref.at(m));
                    aaa.ctsisomer.at(0).insert(aaa.ctsisomer.at(0).begin() + ref, cti_ref.at(0).at(m));
                    aaa.ctsisomer.at(1).insert(aaa.ctsisomer.at(1).begin() + ref, cti_ref.at(1).at(m));
					//aaa.Cyindex.insert(aaa.Cyindex.begin() + ref, cyi.at(m));
					//if (para.protect) aaa.protect.insert(aaa.protect.begin() + ref, pri.at(m));
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
				//aaa.Cyindex.push_back(cyi.at(m));
                aaa.Cyindex.resize(aaa.Cyindex.size()+1,vector<int>(0));
                for (int k1=0;k1<cyi.at(m).size();k1++) aaa.Cyindex.at(aaa.Cyindex.size()-1).push_back(cyi.at(m).at(k1));
                //aaa.ctsisomer.push_back(cti_ref.at(m));
                aaa.ctsisomer.at(0).push_back(cti_ref.at(0).at(m));
                aaa.ctsisomer.at(1).push_back(cti_ref.at(1).at(m));
                //aaa.ctsisomer.resize(aaa.ctsisomer.size()+1,vector<string>(0));
                //for (int k1=0;k1<cti_ref.at(m).size();k1++) aaa.ctsisomer.at(aaa.ctsisomer.size()-1).push_back(cti_ref.at(m).at(k1));
				if (para.protect) aaa.protect.push_back(pri.at(m));
				m++;
				break;
			}
		}
	}
	*/

	Cybnd.resize(0);
	for (n=0;n<cybndj.size();n++) Cybnd.push_back(cybndj.at(n));
    aaa.Cybnd.resize(0);
    for (n=0;n<cybndi.size();n++) aaa.Cybnd.push_back(cybndi.at(n));

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
    vector<int>().swap(cybndj);
    vector<int>().swap(cybndi);


	//if (1) ring_no_chk(-1);

	del_unpaired_ring_no();
	aaa.del_unpaired_ring_no();
	decyc_small_ring(5);
	aaa.decyc_small_ring(5);

	vector<int> PCofct(0);
    for (int i=0;i<Cindex.size();i++) {
		if (Rindex.at(i)==2) {
			PCofct.push_back(Pindex.at(i));
			PCofct.push_back(Cindex.at(i));
		}
    }
	for (int i=0;i<PCofct.size();i+=2) {
		if (PCofct.at(i)>0) {
			if (ctsisomer.at(0).at(PCofct.at(i)-1)!="" && ctsisomer.at(1).at(PCofct.at(i+1)-1)=="") {
				if (ctsisomer.at(0).at(PCofct.at(i)-1)=="/") {
					if (!cistrans) ctsisomer.at(1).at(PCofct.at(i+1)-1)="/";
					else ctsisomer.at(1).at(PCofct.at(i+1)-1)="\\";
				}
                if (ctsisomer.at(0).at(PCofct.at(i)-1)=="\\") {
                    if (!cistrans) ctsisomer.at(1).at(PCofct.at(i+1)-1)="\\";
                    else ctsisomer.at(1).at(PCofct.at(i+1)-1)="/";
                }
			}
			if (ctsisomer.at(0).at(PCofct.at(i)-1)=="" && ctsisomer.at(1).at(PCofct.at(i+1)-1)!="") {
				if (ctsisomer.at(1).at(PCofct.at(i+1)-1)=="/") {
					if (!cistrans) ctsisomer.at(0).at(PCofct.at(i)-1)="/";
					else ctsisomer.at(0).at(PCofct.at(i)-1)="\\";
				}
                if (ctsisomer.at(1).at(PCofct.at(i+1)-1)=="\\") {
                    if (!cistrans) ctsisomer.at(0).at(PCofct.at(i)-1)="\\";
                    else ctsisomer.at(0).at(PCofct.at(i)-1)="/";
                }
			}
			if (ctsisomer.at(0).at(PCofct.at(i)-1)=="" && ctsisomer.at(1).at(PCofct.at(i+1)-1)=="") {
				if (!cistrans) {
					ctsisomer.at(1).at(PCofct.at(i+1)-1)="/";
					ctsisomer.at(0).at(PCofct.at(i)-1)="/";
				}
				else {
                    ctsisomer.at(1).at(PCofct.at(i+1)-1)="\\";
                    ctsisomer.at(0).at(PCofct.at(i)-1)="/";
				}
			}
		}
	}

	//reset();
	//chk_cistrans(); //20200806
	mds2smi();

	//aaa.reset();
	//aaa.chk_cistrans(); //20200806
	aaa.mds2smi();

	//time1=time(NULL)-time1;
	//cout << "CROSSOVER FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}


void MOLECULE::read(string a) {
	for (int i=0;i<a.length();i++) {
		if (a[i]=='/') a[i]='u';
		if (a[i]=='\\') a[i]='d';
		if (a[i]=='*') a[i]='x';
	}
	
	ifstream inf((a+".mds").c_str());
	string b,b1;
	int i=0,j=0,k=0;
	inf >> ws;
	
	ctsisomer.resize(0,vector<string>(0));
	ctsisomer.resize(2,vector<string>(0));

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
			Cyindex.resize(j,vector<int>(0));
			for (i=0;i<j;i++) {
                b1="";
                inf>>b1>>ws;

                string buf="";
                for (int k1=0;k1<b1.length();k1++) {
                    if (b1[k1]!=',') buf+=b1[k1];
                    if (b1[k1]==',') {
                        if (buf!="") {
							Cyindex.at(i).push_back(atoi(buf.c_str()));
							//cout << "yy " << buf << ",";
						}
                        buf="";
                    }
					if (b1[k1]!=',' && k1>=b1.length()-1) {
						if (buf!="0") {
							Cyindex.at(i).push_back(atoi(buf.c_str()));
							//cout << "ll " << buf << endl;;
						}
						buf="";
					}
                }

				//inf>>k>>ws;
				//Cyindex.push_back(k);
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
        else if (b=="ctsisomer_start") {
        	for (i=0;i<j;i++) {
            	b1="";
        		inf>>b1>>ws;

        		if (b1=="_") ctsisomer.at(0).push_back("");
        		else ctsisomer.at(0).push_back(b1);
        	}
        }
        else if (b=="ctsisomer_end") {
            for (i=0;i<j;i++) {
                b1="";
                inf>>b1>>ws;

                if (b1=="_") ctsisomer.at(1).push_back("");
                else ctsisomer.at(1).push_back(b1);
            }
        }
		else if (para.protect && b=="protection") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				protect.push_back(k);
			}
		}
        else if (b=="Cybnd") {
			//Cybnd.resize(if_circle,vector<int>(0));
            for (i=0;i<if_circle;i++) {
                inf>>k>>ws;
                Cybnd.push_back(k);
            }
        }
		//else if (b=="SMILES") {
		//	b1="";
        //    inf>>b1>>ws;
		//	molesmi=smiles=b1;
		//}

	}
	inf.close();
	return;
}


int MOLECULE::cyclization(int pt1,int pt2,int pbnd) {
	//if (if_circle>=9) return 0; //20200510
	if (pt1==pt2) return 0;
	//if (Cyindex.at(pt1).size()) return 0;
	//if (Cyindex.at(pt2).size()) return 0;
	if (para.protect && protect.at(pt1)) return 0;
	if (para.protect && protect.at(pt2)) return 0;

	//int i,j,k,n=0,m=0,x,y;
	//int b_pos[2][2];

	//long double time1=time(NULL);
	
	
	if (1) {
	    vector<int> C_ringmember(0);
		
	    int rep_quota=1;
	    
		int compa=Cindex.at(pt1);
	    int cont=0;
	    while (compa>=1 && cont<=Cyindex.size()) {
	        C_ringmember.push_back(compa);
	
	        if (compa>=1) {
	            if (Pindex.at(compa-1)>=1) compa=Cindex.at(Pindex.at(compa-1)-1);
	            else break;
	            cont++;
	        }
	        else break;
	    }
		
		compa=Cindex.at(pt2);
		cont=0;
		while (compa>=1 && cont<=Cyindex.size()) {
	    	C_ringmember.push_back(compa);
	
	    	if (compa>=1) {
	        	if (Pindex.at(compa-1)>=1) compa=Cindex.at(Pindex.at(compa-1)-1);
	        	else break;
	        	cont++;
	    	}
	    	else break;
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
		}
		vector<int>().swap(C_ringmember);
		if (num_rmember<5) return 0;		
	}
	
	// Ensure the bonds of atom P(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
	if (1) {
        bool totbndchk=0;
        bool bndchk=0;

	    int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
		vector<int> ordcount(6,0);
		if (Rindex.at(pt1)>=1) {
			bndsum+=Rindex.at(pt1);
			ordcount.at(Rindex.at(pt1)-1)+=1;
		}
		if (Cyindex.at(pt1).size()) {
			for (int k1=0;k1<Cyindex.at(pt1).size();k1++) {
				bndsum+=Cybnd.at(Cyindex.at(pt1).at(k1)-1);
				ordcount.at(Cybnd.at(Cyindex.at(pt1).at(k1)-1)-1)++;
			}
		}
    	for (int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(pt1)) {
				bndsum+=Rindex.at(k1);
				ordcount.at(Rindex.at(k1)-1)+=1;
			}
		}
		//bndsum+=1;
		//ordcount.at(0)+=1;
        bndsum+=pbnd;
        ordcount.at(pbnd-1)+=1;


        int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(pt1)].order.at(k2);
        }
		/*
        if (id_bndsum<=id_maxbnd) totbndchk=1;
        else {
			totbndchk=0;
			return 0;
		}
		*/

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[Mindex.at(pt1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                bndchk=0;
				return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(pt1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

	}
	
	if (1) {
        bool totbndchk=0;
        bool bndchk=0;

	    int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
		vector<int> ordcount(6,0);
		if (Rindex.at(pt2)>=1) {
			bndsum+=Rindex.at(pt2);
			ordcount.at(Rindex.at(pt2)-1)+=1;
		}
		if (Cyindex.at(pt2).size()) {
            for (int k1=0;k1<Cyindex.at(pt2).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(pt2).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(pt2).at(k1)-1)-1)++;
            }
		}
    	for (int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(pt2)) {
				bndsum+=Rindex.at(k1);
				ordcount.at(Rindex.at(k1)-1)+=1;
			}
		}
		//bndsum+=1;
		//ordcount.at(0)+=1;
        bndsum+=pbnd;
        ordcount.at(pbnd-1)+=1;


        int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(pt2)].order.at(k2);
        }
		/*
        if (id_bndsum<=id_maxbnd) totbndchk=1;
        else {
			totbndchk=0;
			return 0;
		}
		*/

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[Mindex.at(pt2)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                bndchk=0;
				return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(pt2)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }
	}
	

	if_circle++;
	Cyindex.at(pt1).push_back(if_circle);
	Cyindex.at(pt2).push_back(if_circle);
	//Cyindex.at(pt1)=10*Cyindex.at(pt1)+if_circle;
	//Cyindex.at(pt2)=10*Cyindex.at(pt2)+if_circle;
	Cybnd.resize(if_circle,0);
	Cybnd.at(if_circle-1)=pbnd;

	//cout << "K1 " << endl;
	//print();
	
    del_unpaired_ring_no();
    decyc_small_ring(5);

	//cout << "K2 " << endl;
	//print();

	//reset();
	//chk_cistrans(); //20200806
	mds2smi();

	//time1=time(NULL)-time1;
    //cout << "CYCLIZATION FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}

void MOLECULE::mds2smi(bool ct_on)
{

	if (!ct_on) chk_cistrans(0,Cindex.size()-1);
	reset();

	int totalleng=0;
	//for (int i=0;i<Cindex.size();i++) totalleng+=data->a[Mindex.at(i)].name.length(); //20191011
	for (int i=0;i<atomsmi.size();i++) totalleng+=atomsmi.at(i).size();
	
	
	if (if_circle) totalleng+=6*if_circle; //2*if_circle
	//for (int i=0;i<ctsisomer.size();i++) {
	//	if (ctsisomer.at(i)!="") {
	//		string name=ctsisomer.at(i);
	//		totalleng+=name.length();
	//	}
	//}
	for (int j=0;j<2;j++) {
    	for (int i=0;i<ctsisomer.at(j).size();i++) {
    		if (ctsisomer.at(j).at(i)!="") {
    			string name=ctsisomer.at(j).at(i);
    			totalleng+=name.length();
    		}
    	}
	}
	totalleng+=10;
	//string x(totalleng,' ');
	//string x="";
	//for (int i=0;i<totalleng;i++) x+=" ";
	char *x=new char [totalleng];
	for (int i=0;i<totalleng;i++) x[i]=' ';
	
	vector<int> sdelay(Cindex.size(),0);
	
	bool cts=0;
	for (int i=0;i<ctsisomer.at(0).size();i++) {
		if (ctsisomer.at(0).at(i)!="")  {
			cts=1;
			break;
		}
	}
	if (cts) {
		vector<int> pos,m;
		vector<string> sgn;
		for (int i=0;i<ctsisomer.at(0).size();i++) { 
			if (ctsisomer.at(0).at(i)!="") {
				string n2=ctsisomer.at(0).at(i);
				for (int j=0;j<n2.length();j++) {
					m.push_back(i);
					pos.push_back(atomsmi.at(i).at(0));
					sdelay.at(i)+=1;
					//sgn.push_back(string(1,n2.at(j)));
					sgn.push_back(string(1,n2.at(n2.length()-j-1)));
				}		
			}
		}
		pos.reserve(pos.size());
		m.reserve(m.size());
		sgn.reserve(sgn.size());
		
		for (int i=0;i<pos.size();i++) {
			for (int j=i+1;j<pos.size();j++) {
				if (pos.at(i)>pos.at(j)) {  //at
					swap(pos.at(i),pos.at(j)); //at
					swap(m.at(i),m.at(j)); //at
					swap(sgn.at(i),sgn.at(j)); //at
				}
			}
		}
		
		for (int i=0;i<m.size();i++) {
			int posi=0,len=1;
			string sg=sgn.at(i);
			len=sg.length();

			posi=atomsmi.at(m.at(i)).at(0); //at
			for (int k=0;k<len;k++) {
				atomsmi.at(m.at(i)).insert(atomsmi.at(m.at(i)).begin(),-9);
			}
			
			for (int j=0;j<Cindex.size();j++) {
				for (int k=0;k<atomsmi.at(j).size();k++) {
					if (atomsmi.at(j).at(k)>=posi) atomsmi.at(j).at(k)+=len;
				}
			}
			
			for (int j=0;j<atomsmi.at(m.at(i)).size();j++) {
				if (atomsmi.at(m.at(i)).at(j)==-9) {
					//int g=data->a.at(Mindex.at(m.at(i))).index-1+sdelay.at(m.at(i));
					int g=atomsmi.at(m.at(i)).size();
					if (j>0) {
						for (int k1=1;k1<=g;k1++) { // atomsmi[m.at(i)].size()
							int left=-2;
							if (j-k1>=0) left=atomsmi.at(m.at(i)).at(j-k1);
							if (left!=-2 && left!=-9) {
								for (int k2=(j-k1);k2<=g;k2++) { //orig k2<atomsmi[m.at(i)].size()
									if (k2+1<g) { // orig k2+1<atomsmi[m.at(i)].size()
										if (atomsmi.at(m.at(i)).at(k2+1)==-9) atomsmi.at(m.at(i)).at(k2+1)=atomsmi.at(m.at(i)).at(k2)+1;
									}
								}
								break;
							}						
						}						
					}
					else if (j==0) {
						for (int k1=1;k1<=g;k1++) {
							int right=-1;
							if (j+k1<atomsmi.at(m.at(i)).size()) right=atomsmi.at(m.at(i)).at(j+k1);								
							if (right!=-1 && right!=-9) {
								for (int k2=(j+k1);k2>=0;k2--) {
									if (k2-1>=0) {
										if (atomsmi.at(m.at(i)).at(k2-1)==-9) atomsmi.at(m.at(i)).at(k2-1)=atomsmi.at(m.at(i)).at(k2)-1;
									}
								}
								break;
							}
						}						
					}
					
					//for (int k=0;k<len;k++) x.insert(atomsmi[m.at(i)].at(j),string(1,sg[len-k-1]));
					for (int k=totalleng-1;k-len>=atomsmi.at(m.at(i)).at(j);k--) {
						x[k]=x[k-len];
						x[k-len]=' ';
					}
					for (int k=0;k<len;k++) {
						x[atomsmi.at(m.at(i)).at(j)+k]=sg[k];
					}
						
					break;
				}
			}

		}
		pos.clear();
		m.clear();
		sgn.clear();
		vector<int>().swap(pos);
		vector<int>().swap(m);
		vector<string>().swap(sgn); 
	}

    vector<int> edelay(Cindex.size(),0);
    for (int i=0;i<ctsisomer.at(1).size();i++) {
        if (ctsisomer.at(1).at(i)!="")  {
            cts=1;
            break;
        }
    }
    if (cts) {
        vector<int> pos,m;
        vector<string> sgn;
        for (int i=0;i<ctsisomer.at(1).size();i++) {
            if (ctsisomer.at(1).at(i)!="") {
                string n2=ctsisomer.at(1).at(i);
                for (int j=0;j<n2.length();j++) {
                    m.push_back(i);
                    pos.push_back(atomsmi.at(i).at(0));
                    edelay.at(i)+=1;
                    //sgn.push_back(string(1,n2.at(j)));
                    sgn.push_back(string(1,n2.at(n2.length()-j-1)));
                }
            }
        }
        pos.reserve(pos.size());
        m.reserve(m.size());
        sgn.reserve(sgn.size());

        for (int i=0;i<pos.size();i++) {
            for (int j=i+1;j<pos.size();j++) {
                if (pos.at(i)>pos.at(j)) {  //at
                    swap(pos.at(i),pos.at(j)); //at
                    swap(m.at(i),m.at(j)); //at
                    swap(sgn.at(i),sgn.at(j)); //at
                }
            }
        }

        for (int i=0;i<m.size();i++) {
            int posi=0,len=1;
            string sg=sgn.at(i);
            len=sg.length();

            bool siz=0;

            if (atomsmi.at(m.at(i)).size()>data->a.at(Mindex.at(m.at(i))).index-1+sdelay.at(m.at(i))) siz=0;
            else siz=1;  //20200910

            if (!siz) posi=atomsmi.at(m.at(i)).at(data->a.at(Mindex.at(m.at(i))).index-1+sdelay.at(m.at(i)));
            else posi=atomsmi.at(m.at(i)).at(data->a.at(Mindex.at(m.at(i))).index-2+sdelay.at(m.at(i)));

            for (int k=0;k<len;k++) {
                if (!siz) atomsmi.at(m.at(i)).insert(atomsmi.at(m.at(i)).begin()+data->a.at(Mindex.at(m.at(i))).index-1+sdelay.at(m.at(i)),-9);
                else atomsmi.at(m.at(i)).push_back(-9);  //20200910
            }

            for (int j=0;j<Cindex.size();j++) {
                for (int k=0;k<atomsmi.at(j).size();k++) {
                    if (!siz) {
                        if (atomsmi.at(j).at(k)>=posi) atomsmi.at(j).at(k)+=len;
                    }
                    else { //20200910
                        if (atomsmi.at(j).at(k)>posi) atomsmi.at(j).at(k)+=len;
                    }
                }
            }

            for (int j=0;j<atomsmi.at(m.at(i)).size();j++) {
                if (atomsmi.at(m.at(i)).at(j)==-9) {
                    //int g=data->a.at(Mindex.at(m.at(i))).index-1+sdelay.at(m.at(i))+ringdelay.at(m.at(i));
                    int g=atomsmi.at(m.at(i)).size();
                    if (j>0) {
                        for (int k1=1;k1<=g;k1++) { //orig k1<atomsmi[m.at(i)].size()
                            int left=-2;
                            if (j-k1>=0) left=atomsmi.at(m.at(i)).at(j-k1);
                            if (left!=-2 && left!=-9) {
                                for (int k2=(j-k1);k2<=g;k2++) { //orig k2<atomsmi[m.at(i)].size()
                                    if (k2+1<atomsmi.at(m.at(i)).size()) {
                                        if (atomsmi.at(m.at(i)).at(k2+1)==-9) atomsmi.at(m.at(i)).at(k2+1)=atomsmi.at(m.at(i)).at(k2)+1;
                                    }
                                }
                                break;
                            }
                        }
                    }
                    else if (j==0) {
                        for (int k1=1;k1<=g;k1++) { // orig k1<atomsmi[m.at(i)].size()
                            int right=-1;
                            if (j+k1<atomsmi.at(m.at(i)).size()) right=atomsmi.at(m.at(i)).at(j+k1);

                            if (right!=-1 && right!=-9) {
                                for (int k2=(j+k1);k2>=0;k2--) {
                                    if (k2-1>=0) {
                                        if (atomsmi.at(m.at(i)).at(k2-1)==-9) atomsmi.at(m.at(i)).at(k2-1)=atomsmi.at(m.at(i)).at(k2)-1;
                                    }
                                }
                                break;
                            }
                        }
                    }
                    //for (int k=0;k<len;k++) x.insert(atomsmi[m.at(i)].at(j),string(1,sg[len-k-1]));
                    for (int k=totalleng-1;k-len>=atomsmi.at(m.at(i)).at(j);k--) {
                        x[k]=x[k-len];
                        x[k-len]=' ';
                    }
                    for (int k=0;k<len;k++) {
                        x[atomsmi.at(m.at(i)).at(j)+k]=sg[k];
                    }

                    break;
                }
            }

        }
        pos.clear();
        m.clear();
		sgn.clear();
        vector<int>().swap(pos);
        vector<int>().swap(m);
        vector<string>().swap(sgn);

	}

	vector<int> ringdelay(Cindex.size(),0);
	
	if (if_circle) {
		vector<int> pos,m;
		vector< vector<int> > num(0,vector<int>(0));
		vector< vector<int> > cyend(Cyindex.size(),vector<int>(0));
		vector<bool> iscyend(Cybnd.size(),0);
		for (int i=0;i<Cyindex.size();i++) { 
			if (Cyindex.at(i).size()) {
				m.push_back(i);
				pos.push_back(atomsmi.at(i).at(0));  //20200910
				num.resize(num.size()+1,vector<int>(0));

				int digits=0;
				for (int j=0;j<Cyindex.at(i).size();j++) {
					num.at(num.size()-1).push_back(Cyindex.at(i).at(j));

					if (Cyindex.at(i).at(j)<10) digits+=(int)(log10(Cyindex.at(i).at(j))+2);
					else if (Cyindex.at(i).at(j)>=10) digits+=(int)(log10(Cyindex.at(i).at(j))+1);

					if (!iscyend.at(Cyindex.at(i).at(j)-1)) iscyend.at(Cyindex.at(i).at(j)-1)=1;
					else cyend.at(i).push_back(Cyindex.at(i).at(j));
				}

				ringdelay.at(i)+=digits+1;
			}
		}
		
        for (int i=0;i<cyend.size();i++) {
            if (cyend.at(i).size()) {
				int ct=0;
				for (int j=0;j<cyend.at(i).size();j++) {
					if (Cybnd.at(cyend.at(i).at(j)-1)==2 || Cybnd.at(cyend.at(i).at(j)-1)==3) ct++;
				}
				ringdelay.at(i)+=3*ct;
			}	
		}

		pos.reserve(pos.size());
		m.reserve(m.size());
		num.reserve(num.size());
		
		for (int i=0;i<pos.size();i++) {
			for (int j=i+1;j<pos.size();j++) {
				if (pos.at(i)>pos.at(j)) {  //at
					swap(pos.at(i),pos.at(j)); //at
					swap(m.at(i),m.at(j)); //at
					swap(num.at(i),num.at(j)); //at
				}
			}
		}


		for (int i=0;i<m.size();i++) {
			int posi=0;
			int digits=0;
			for (int k=0;k<num.at(i).size();k++) {
				if (num.at(i).at(k)<10) digits+=2+(int)(log10(num.at(i).at(k))+1);
				else if (num.at(i).at(k)>=10) digits+=1+(int)(log10(num.at(i).at(k))+1);
			}

			stringstream nn("");

			bool siz=0;

			if (atomsmi.at(m.at(i)).size()>data->a.at(Mindex.at(m.at(i))).index-1+sdelay.at(m.at(i))) siz=0;
			else siz=1;  //20200910
			
			if (!siz) posi=atomsmi.at(m.at(i)).at(data->a.at(Mindex.at(m.at(i))).index-1+sdelay.at(m.at(i)));
			else posi=atomsmi.at(m.at(i)).at(data->a.at(Mindex.at(m.at(i))).index-2+sdelay.at(m.at(i)));

			for (int k=0;k<num.at(i).size();k++) {
				for (int k1=0;k1<cyend.at(m.at(i)).size();k1++) {
					if (num.at(i).at(k)==cyend.at(m.at(i)).at(k1)) {
						if (Cybnd.at(num.at(i).at(k)-1)==2) nn << "(=)";
						else if (Cybnd.at(num.at(i).at(k)-1)==3) nn << "(#)";

						break;
					}
				}
				if (num.at(i).at(k)<10) nn << "%" << num.at(i).at(k) << "%";
				else if (num.at(i).at(k)>=10) nn << "%" << num.at(i).at(k);
			}

			//cout << nn.str() << endl;

            int ct=0;
            if (cyend.at(m.at(i)).size()) {
                for (int j=0;j<cyend.at(m.at(i)).size();j++) {
                    if (Cybnd.at(cyend.at(m.at(i)).at(j)-1)==2 || Cybnd.at(cyend.at(m.at(i)).at(j)-1)==3) ct++;
                }
            }
            digits+=3*ct;

			for (int k=0;k<digits;k++) { // <digits
				if (!siz) atomsmi.at(m.at(i)).insert(atomsmi.at(m.at(i)).begin()+data->a.at(Mindex.at(m.at(i))).index-1+sdelay.at(m.at(i)),-9);
				else atomsmi.at(m.at(i)).push_back(-9);  //20200910
			}
			
			for (int j=0;j<Cindex.size();j++) {
				for (int k=0;k<atomsmi.at(j).size();k++) {
					if (!siz) {
						if (atomsmi.at(j).at(k)>=posi) atomsmi.at(j).at(k)+=digits; //+=digits
					}
					else { //20200910
						if (atomsmi.at(j).at(k)>posi) atomsmi.at(j).at(k)+=digits; //+=digits
					}	
				}
			}

			for (int j=0;j<atomsmi.at(m.at(i)).size();j++) {
				if (atomsmi.at(m.at(i)).at(j)==-9) {
					int g=atomsmi.at(m.at(i)).size();
					if (j>0) {
						for (int k1=1;k1<=g;k1++) { //orig k1<atomsmi[m.at(i)].size()
							int left=-2;
							if (j-k1>=0) left=atomsmi.at(m.at(i)).at(j-k1);
							if (left!=-2 && left!=-9) {
								for (int k2=(j-k1);k2<=g;k2++) { //orig k2<atomsmi[m.at(i)].size()
									if (k2+1<atomsmi.at(m.at(i)).size()) {
										if (atomsmi.at(m.at(i)).at(k2+1)==-9) atomsmi.at(m.at(i)).at(k2+1)=atomsmi.at(m.at(i)).at(k2)+1;
									}
								}
								break;
							}						
						}						
					}
					else if (j==0) {
						for (int k1=1;k1<=g;k1++) { // orig k1<atomsmi[m.at(i)].size()
							int right=-1;
							if (j+k1<atomsmi.at(m.at(i)).size()) right=atomsmi.at(m.at(i)).at(j+k1);
													
							if (right!=-1 && right!=-9) {
								for (int k2=(j+k1);k2>=0;k2--) {
									if (k2-1>=0) {
										if (atomsmi.at(m.at(i)).at(k2-1)==-9) atomsmi.at(m.at(i)).at(k2-1)=atomsmi.at(m.at(i)).at(k2)-1;
									}
								}
								break;
							}
						}						
					}
					
					string fp=nn.str();
                    for (int k=totalleng-1;k-digits>=atomsmi.at(m.at(i)).at(j);k--) { //k-digits>=
                        x[k]=x[k-digits]; //k-digits
                        x[k-digits]=' '; //k-digits
                    }
                    for (int k=0;k<digits;k++) { //k<digits
                        x[atomsmi.at(m.at(i)).at(j)+k]=fp[k]; //=fp[k]
                    }


					break;
				}
			}


		}
		pos.clear();
		m.clear();
		num.clear();
		vector<int>().swap(pos);
		vector<int>().swap(m);
		vector< vector<int> >().swap(num);

	    //for (int i=0;i<totalleng;i++) cout << x[i];
    	//cout << endl;

	}


	vector<int>().swap(sdelay);
	vector<int>().swap(ringdelay);

	chg=0;
	for (int i=0;i<Cindex.size();i++) chg+=data->a.at(Mindex.at(i)).chg;

	int implicit_ring_n=0;
	for (int i=0;i<Cindex.size();i++) {
		int M = Mindex.at(i);
		int R = Rindex.at(i);
		string kk=data->a.at(M).name; //20191011

		bool explicit_H=0;
		if (data->a.at(M).chg && kk.length()<=data->a.at(M).nbond && data->a.at(M).index!=data->a.at(M).nbond) explicit_H=1;
		if (M==67 || M==68) explicit_H=1;

		if (explicit_H) { // data.a[M].chg && kk.length()<=data.a[M].nbond && data.a[M].name[0]=='['
			int tmp=0;
			for (int j=0;j<data->a.at(M).norder;j++) {
				tmp+=Bindex.at(i).at(j);
			}
			if (1) {
				if (tmp>0 && tmp<10) kk[3]=char(tmp+48); //data->a[M].name[3]=char(tmp+48);
			}
		}
		
		//if (data->a.at(M).cypos.size()) {
		//	implicit_ring_n++;
		//	for (int j=data->a.at(M).cypos.at(0).size()-1;j>=0;j--) {
		//		if (if_circle+implicit_ring_n>0) {
		//			int dig=(int)log10(if_circle+implicit_ring_n)+1;
		//			stringstream buf("");
		//			buf << if_circle+implicit_ring_n;
		//			kk.erase(kk.begin()+data->a.at(M).cypos.at(1).at(j));
		//			kk.insert(kk.begin()+data->a.at(M).cypos.at(1).at(j),buf.str());
        //          kk.erase(kk.begin()+data->a.at(M).cypos.at(0).at(j));
        //          kk.insert(kk.begin()+data->a.at(M).cypos.at(0).at(j),buf.str());
		//		}
		//		//kk[data->a.at(M).cypos.at(0).at(j)]=char(if_circle+implicit_ring_n+48);
		//		//kk[data->a.at(M).cypos.at(1).at(j)]=char(if_circle+implicit_ring_n+48);
		//	}
		//}

        int num_3blanks=0;
        if (1 && data->a.at(M).norder>0) { //delete "(-)" , "(=)" , "(#)"
			if (1) {
	            for (int j=0;j<kk.size();j++) {
                    if (j>0 && j<kk.size()-1) {
						if (Rindex.at(i)==1) {
                        	if (kk[j-1]=='(' && kk[j]=='-' && kk[j+1]==')') {
                            	kk[j-1]=kk[j]=kk[j+1]=' ';
                            	num_3blanks++;
								break;
                        	}
						}
                        if (Rindex.at(i)==2) {
                            if (kk[j-1]=='(' && kk[j]=='=' && kk[j+1]==')') {
                                kk[j-1]=kk[j]=kk[j+1]=' ';
                                num_3blanks++;
                                break;
                            }
                        }
                        if (Rindex.at(i)==3) {
                            if (kk[j-1]=='(' && kk[j]=='#' && kk[j+1]==')') {
                                kk[j-1]=kk[j]=kk[j+1]=' ';
                                num_3blanks++;
                                break;
                            }
                        }

                    }
                }
				for (int j=0;j<Bindex.at(i).size();j++) {
					if (Bindex.at(i).at(j)) {
                		for (int k=0;k<kk.size();k++) {
                    		if (k>0 && k<kk.size()-1) {
                        		if (Bindex.at(i).at(j)==1) {
                            		if (kk[k-1]=='(' && kk[k]=='-' && kk[k+1]==')') {
                                		kk[k-1]=kk[k]=kk[k+1]=' ';
                                		num_3blanks++;
										break;
                            		}
                        		}
                        		if (Bindex.at(i).at(j)==2) {
                            		if (kk[k-1]=='(' && kk[k]=='=' && kk[k+1]==')') {
                                		kk[k-1]=kk[k]=kk[k+1]=' ';
                                		num_3blanks++;
										break;
                            		}
                        		}
                        		if (Bindex.at(i).at(j)==3) {
                            		if (kk[k-1]=='(' && kk[k]=='#' && kk[k+1]==')') {
                                		kk[k-1]=kk[k]=kk[k+1]=' ';
                                		num_3blanks++;
										break;
                            		}
                        		}

                    		}
                		}

					}
				}
                for (int j=0;j<kk.length();j++) {
                    if (kk[j]==' ') {
                        kk.erase(kk.begin()+j);
                        j--;
                    }
                }
			}
        }
		
        if (1 && Cyindex.at(i).size()) { //delete "(-)", "(=)", and "(#)" of cyclic flag
			for (int j=0;j<Cyindex.at(i).size();j++) {
				for (int k=0;k<kk.size();k++) {
					if (k>0 && k<kk.size()-1) {
						if (Cybnd.at(Cyindex.at(i).at(j)-1)==1) {
							if (kk[k-1]=='(' && kk[k]=='-' && kk[k+1]==')') {
								kk[k-1]=kk[k]=kk[k+1]=' ';
	                            num_3blanks++;
    	                        break;
							}
						}
                        if (Cybnd.at(Cyindex.at(i).at(j)-1)==2) {
                            if (kk[k-1]=='(' && kk[k]=='=' && kk[k+1]==')') {
                                kk[k-1]=kk[k]=kk[k+1]=' ';
                                num_3blanks++;
                                break;
                            }
                        }
                        if (Cybnd.at(Cyindex.at(i).at(j)-1)==3) {
                            if (kk[k-1]=='(' && kk[k]=='#' && kk[k+1]==')') {
                                kk[k-1]=kk[k]=kk[k+1]=' ';
                                num_3blanks++;
                                break;
                            }
                        }

					}
				}

            	for (int j=0;j<kk.length();j++) {
                	if (kk[j]==' ') {
                    	kk.erase(kk.begin()+j);
                    	j--;
                	}
            	}

			}

        }
		

		int num_1blank=0;
		if (1 && data->a.at(M).norder) { // //delete "-" of "....(-...."
			if (1) {
				for (int j=0;j<kk.size();j++) {
					if (j>0 && j<kk.size()-1) {
						if (kk[j-1]=='(' && kk[j]=='-' && kk[j+1]==')') {
							kk[j]=' ';
							num_1blank++;
						}
					}
				}
	            for (int j=0;j<kk.length();j++) {
    	            if (kk[j]==' ') {
        	            kk.erase(kk.begin()+j);
            	        j--;
                	}
            	}
			}
		}
		
		int ct=0;
		for (int j=0;j<atomsmi.at(i).size();j++) { // orig data->a[M].name.length() // data.a[M].nbond  //w<j  //data->a[M].name.length()
			int n=atomsmi.at(i).at(j);

			if (x[n]==' ') {
				//string nu=string(1,kk[ct]); //string(1,data->a[Mindex.at(i)].name[ct]);
				x[n]=kk[ct] ;//nu[0]; //data->a[Mindex.at(i)].name[j];
				//if (ct<data->a.at(M).nbond-3*num_3blanks-num_1blank) { // j<data->a[M].nbond
					//if (((ct-tmp)%3)==0 && ct>data->a[Mindex.at(i)].index && data->a[Mindex.at(i)].index<data->a[Mindex.at(i)].nbond) { 
						//x[n]=nu[0]; //data->a[Mindex.at(i)].name[j];
						//bool bnd=0;
						//if (n>=2) {
						//	if (x[n-1]=='-' || x[n-1]=='=' || x[n-1]=='#') bnd=1;
						//	if (x[n]==')' && x[n-2]=='(' && bnd) x[n-2]=x[n-1]=x[n]=' ';
						//}
					//}
					//else {
					//	x[n]=nu[0]; //data->a[Mindex.at(i)].name[j];
					//}
				//}
				//else if (ct>=data->a[M].nbond && ct<data->a[M].name.length()) { // orig j>=data->a[M].nbond && j<data->a[M].name.length()
					//x[n]=nu[0]; //data->a[Mindex.at(i)].name[j];
				//}
				ct++;

			}

		}

		//for (int g=0;g<totalleng;g++) cout << x[g];
		//cout << endl;
		
		//adata->a[M].name=kk;  //20191201
		//if (M==67 || M==68) Mindex.at(i)=1;  //20200822
	}

    atomsmi.resize(0);
    //atomsmi.shrink_to_fit();
    Bindex.resize(0);
    //Bindex.shrink_to_fit();

	//molesmi="";
	molesmi.clear();
	for (int i=0;i<totalleng;i++) {
		if (x[i]!=' ') {
			molesmi.push_back(x[i]);
		}
	}
    delete [] x;
    x=NULL;

	smiles=molesmi;

	//cout << "SMI: " << molesmi << endl;

    //if (cts) {
    //    ofstream chkct("./chkct",ios::app);
    //    chkct << molesmi << endl;
    //    chkct.close();
    //}

	return;
}


void MOLECULE::clear() {
	int n;
	molesmi.clear();
	//Bindex.clear();
	//atomsmi.clear();
	Bindex.resize(0);
	atomsmi.resize(0);
	//ctsisomer.clear();

	return;
}

void MOLECULE::chk_cistrans(int sposi, int lposi) {
	//vector<bool> posi(Cindex.size(),0);

	if (1) {
		if (ctsisomer.size()==2) {
			if (ctsisomer.at(0).size()==ctsisomer.at(1).size() && ctsisomer.at(0).size()==Cindex.size()) return;
		}
	}

    if (1) {
        //ctsisomer.clear();
        //ctsisomer.resize(Cindex.size(),"");
		ctsisomer.resize(0,vector<string> (0));
		ctsisomer.resize(2,vector<string>(Cindex.size(),""));
        //for (int i=0;i<posi.size();i++) posi.at(i)=1;
    }
	/*
	if (0) {
		if (ctsisomer.size()==Cindex.size()) {
			for (int i=0;i<ctsisomer.size();i++) {
				if (ctsisomer.at(i)!="/" && ctsisomer.at(i)!="\\" && ctsisomer.at(i)!="") ctsisomer.at(i)="";
			}
			for (int i=0;i<ctsisomer.size();i++) { //int i=Cindex.size()-1;i>=0;i-- //int i=0;i<ctsisomer.size();i++
				if (Rindex.at(i)==2 && Mindex.at(i)==2) {
					int P=Pindex.at(i);
					if (P>0) {
						if (Mindex.at(P-1)==2 && Cindex.at(P-1)>0) {
							vector<string> buf(0);
							
							if (Pindex.at(P-1)>0) {
								buf.push_back(ctsisomer.at(Pindex.at(P-1)-1));
								posi.at(Pindex.at(P-1)-1)=1;
							}

							for (int k1=0;k1<Cindex.size();k1++) {
								if (Cindex.at(P-1)==Pindex.at(k1) && Cindex.at(k1)!=Cindex.at(i)) {
									buf.push_back(ctsisomer.at(k1));
									posi.at(k1)=1;
									break;
								}
							}

							if (buf.size()==2) {
                            	if (buf.at(0)=="/" && buf.at(1)=="/") {
									if (Pindex.at(P-1)>0) ctsisomer.at(Pindex.at(P-1)-1)="\\";
                            	}
                            	if (buf.at(0)=="\\" && buf.at(1)=="\\") {
									if (Pindex.at(P-1)>0) ctsisomer.at(Pindex.at(P-1)-1)="/";
                            	}
							}

							//vector<string> buf(0);
							buf.clear();
							for (int k1=0;k1<Cindex.size();k1++) {
								if (Pindex.at(k1)==Cindex.at(i)) {
									buf.push_back(ctsisomer.at(k1));
									posi.at(k1)=1;
								}
							}
							if (buf.size()==2) {
								for (int k1=0;k1<Cindex.size();k1++) {
									if (buf.at(0)=="/" && buf.at(1)=="/") {
        		                	    if (Pindex.at(k1)==Cindex.at(i)) {
											if (ctsisomer.at(k1)=="/") {
												ctsisomer.at(k1)="\\";
												break;
											}
                            			}
                        			}
									if (buf.at(0)=="\\" && buf.at(1)=="\\") {
                        				for (int k1=0;k1<Cindex.size();k1++) {
                                    		if (Pindex.at(k1)==Cindex.at(i)) {
                                        		if (ctsisomer.at(k1)=="\\") {
                                            		ctsisomer.at(k1)="/";
                                            		break;
                                        		}
                                    		}
										}	
									}
								}
							}
						}
					}
				}
			}

			for (int i=0;i<ctsisomer.size();i++) {
				//if (ctsisomer.at(i)!="/" && ctsisomer.at(i)!="\\" && ctsisomer.at(i)!="") ctsisomer.at(i)="";
				if (!posi.at(i)) ctsisomer.at(i)="";
			}

			for (int i=0;i<posi.size();i++) {
				if (posi.at(i) && ctsisomer.at(i)!="") posi.at(i)=0;
			}
			
			//for (int i=0;i<posi.size();i++) {
            //    if (posi.at(i)) break;
			//	else if (!posi.at(i) && i>=posi.size()-1) return;
            //}
		}
		else ctsisomer.clear();

		if (!ctsisomer.size()) {
			ctsisomer.resize(Cindex.size(),"");
			for (int i=0;i<posi.size();i++) posi.at(i)=1;
		}
	}
	*/

    if (1) {
		int flip=0;
        for (int i=sposi;i<=lposi;i++) { // C\C=C\C  //i=lposi;i>=sposi;i-- //i=sposi;i<=lposi;i++
			bool go=0;
			if (Mindex.at(i)==2) go=1;
			if (Mindex.at(i)==8) go=1;
			if (Mindex.at(i)==16) go=1;

            if (Rindex.at(i)==2 && go) {
                int P=Pindex.at(i);

                if (P>0) {
					bool go1=0;
		            if (Mindex.at(P-1)==2) go1=1;
            		if (Mindex.at(P-1)==8) go1=1;
            		if (Mindex.at(P-1)==16) go1=1;

                    if (go1 && Cindex.at(P-1)>0) { // Mindex.at(P-1)==2 && Pindex.at(P-1)>0
						vector< vector<int> > chainsM(4,vector<int>(0));
                        vector< vector<int> > chainsR(4,vector<int>(0));

						bool cas=0;
						if (Pindex.at(P-1)>0) {
							cas=0;
						    if (Cindex.at(Pindex.at(P-1)-1)>0) {
        						int compa=Cindex.at(Pindex.at(P-1)-1);
        						int cont=0;
        						while (compa>=1 && cont<=Cindex.size()) {
            						if (compa>=1) {
                						if (Pindex.at(compa-1)>=1) {
											chainsM.at(0).push_back(Mindex.at(compa-1));
											chainsR.at(0).push_back(Rindex.at(compa-1));



                                    		vector<int> curatm(0);
                                    		curatm.push_back(Cindex.at(compa-1));

                                    		bool goout=0;

                                    		do {
												vector<int> tmpp(0);
                                        		for (int k3=0;k3<curatm.size();k3++) {
                                            		for (int k2=0;k2<Cindex.size();k2++) {
                                                		if (Pindex.at(k2)==curatm.at(k3) && Cindex.at(k2)!=Cindex.at(P-1)) {
                                                    		chainsM.at(0).push_back(Mindex.at(k2));
                                                    		chainsR.at(0).push_back(Rindex.at(k2));
                                                    		tmpp.push_back(Cindex.at(k2));
                                                		}
                                            		}
												}
                                            	if (tmpp.size()>0) {
                                                	curatm.resize(0);
                                                	for (int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
                                            	}
                                            	else {
                                               		goout=1;
                                                	break;
                                            	}
                                    		} while (!goout);


                    						compa=Cindex.at(Pindex.at(compa-1)-1);
                						}
                						else break;
                						cont++;
            						}
            						else break;
        						}
								if (compa==1) {
									chainsM.at(0).push_back(Mindex.at(0));
									chainsR.at(0).push_back(Rindex.at(0));
								}

    						}

                            for (int k1=0;k1<Cindex.size();k1++) {
                                if (Pindex.at(k1)==Cindex.at(P-1) && Cindex.at(k1)!=Cindex.at(i)) {
                                    chainsM.at(1).push_back(Mindex.at(k1));
                                    chainsR.at(1).push_back(Rindex.at(k1));

                                    vector<int> curatm(0);
                                    curatm.push_back(Cindex.at(k1));

                                    bool goout=0;

                                    do {
										vector<int> tmpp(0);
                                        for (int k3=0;k3<curatm.size();k3++) {
                                            for (int k2=0;k2<Cindex.size();k2++) {
                                                if (Pindex.at(k2)==curatm.at(k3)) {
                                                    chainsM.at(1).push_back(Mindex.at(k2));
                                                    chainsR.at(1).push_back(Rindex.at(k2));
                                                    tmpp.push_back(Cindex.at(k2));
                                                }
                                            }
										}
                                        if (tmpp.size()>0) {
                                            curatm.resize(0);
                                            for (int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
                                        }
                                        else {
                                            goout=1;
                                            break;
                                        }
                                    } while (!goout);
                                }
                            }
						}
						else {
							cas=1;
							int ct=0;
                            for (int k1=0;k1<Cindex.size();k1++) {
                                if (Pindex.at(k1)==Cindex.at(P-1) && Cindex.at(k1)!=Cindex.at(i)) {
                                    chainsM.at(ct).push_back(Mindex.at(k1));
                                    chainsR.at(ct).push_back(Rindex.at(k1));

                                    vector<int> curatm(0);
                                    curatm.push_back(Cindex.at(k1));

                                    bool goout=0;

                                    do {
										vector<int> tmpp(0);
                                        for (int k3=0;k3<curatm.size();k3++) {
                                            for (int k2=0;k2<Cindex.size();k2++) {
                                                if (Pindex.at(k2)==curatm.at(k3)) {
                                                    chainsM.at(ct).push_back(Mindex.at(k2));
                                                    chainsR.at(ct).push_back(Rindex.at(k2));
                                                    tmpp.push_back(Cindex.at(k2));
                                                }
                                            }
										}
                                        if (tmpp.size()>0) {
                                            curatm.resize(0);
                                            for (int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
                                        }
                                        else {
                                            goout=1;
                                            break;
                                        }
                                    } while (!goout);

									ct++;
                                }
                            }

						}

						int ct=2;
                        for (int k1=0;k1<Cindex.size();k1++) {
                            if (Pindex.at(k1)==Cindex.at(i)) {
                                chainsM.at(ct).push_back(Mindex.at(k1));
								chainsR.at(ct).push_back(Rindex.at(k1));

								vector<int> curatm(0);
								curatm.push_back(Cindex.at(k1));
									
								bool goout=0;
									
								do {
									vector<int> tmpp(0);
									for (int k3=0;k3<curatm.size();k3++) {
										for (int k2=0;k2<Cindex.size();k2++) {
                                			if (Pindex.at(k2)==curatm.at(k3)) {
                                    			chainsM.at(ct).push_back(Mindex.at(k2));
                                    			chainsR.at(ct).push_back(Rindex.at(k2));
												tmpp.push_back(Cindex.at(k2));
                                			}
										}
									}
									if (tmpp.size()>0) {
										curatm.resize(0);
										for (int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
									}
									else {
										goout=1;
										break;
									}
								} while (!goout);

								ct++;
                            }
                        }

						for (int g=0;g<chainsM.size();g++) {
							for (int g1=0;g1<chainsM.at(g).size();g1++) {
								for (int g2=g1+1;g2<chainsM.at(g).size();g2++) {
									if (chainsM.at(g).at(g1)>chainsM.at(g).at(g2)) swap(chainsM.at(g).at(g1),chainsM.at(g).at(g2));
								}
							}
						}

						if (0) {
                            for (int g=0;g<chainsM.size();g++) {
								cout << "F " << filenum << " | chg " << chg << " | " << g+1 << " | ";
                                for (int g1=0;g1<chainsM.at(g).size();g1++) {
                                    cout << chainsM.at(g).at(g1) << " ";
                                }
								cout << endl;
                            }
						}

						bool hascistrans=1;
						if (chainsM.at(0).size()==chainsM.at(1).size()) {
							if (chainsM.at(0).size()==0) hascistrans=0;
							else {
								for (int g=0;g<chainsM.at(0).size();g++) {
									if (chainsM.at(1).at(g)-chainsM.at(0).at(g)) {
										break;
									}
									else if (chainsM.at(1).at(g)==chainsM.at(0).at(g) && g>=chainsM.at(0).size()-1) {
										hascistrans=0;
									}
								}
							}
						}
						if (hascistrans) {
                            if (chainsM.at(2).size()==chainsM.at(3).size()) {
                                if (chainsM.at(2).size()==0) hascistrans=0;
                                else {
                                    for (int g=0;g<chainsM.at(2).size();g++) {
                                        if (chainsM.at(3).at(g)-chainsM.at(2).at(g)) {
                                            break;
                                        }
										else if (chainsM.at(3).at(g)==chainsM.at(2).at(g) && g>=chainsM.at(2).size()-1) {
											hascistrans=0;
										}
                                    }
                                }
                            }
						}

						bool y1=1,y2=1;
						if (hascistrans && if_circle) {
    						vector<int> C_ringmember(0);

    						for (int k1=if_circle;k1>=1;k1--) {
        						C_ringmember.resize(0);
        						int rep_quota=1;
        						for (int k2=Cyindex.size()-1;k2>=0;k2--) {
            						if (Cyindex.at(k2).size()) {
                						for (int k3=0;k3<Cyindex.at(k2).size();k3++) {
                    						if (Cyindex.at(k2).at(k3)==k1) {
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
        						if (C_ringmember.size()>0) {
                					for (int k2=0;k2<C_ringmember.size();k2++) {
                    					if (C_ringmember.at(k2)==-1) {
											C_ringmember.erase(C_ringmember.begin()+k2);
											k2--;
                    					}
                					}
        						}
								for (int g=0;g<C_ringmember.size();g++) {
									if (Cindex.at(i)==C_ringmember.at(g)) y1=0;
                                    if (Cindex.at(P-1)==C_ringmember.at(g)) y2=0;
								}

							}

						}

						bool z=0;
						if (y1 || y2) z=1;
						if (hascistrans && z) {
							if (0) {
                            	if (0) ctsisomer.at(1).at(Pindex.at(P-1)-1)+="/"; 

								if (1) ctsisomer.at(0).at(Cindex.at(P-1)-1)+="/";

								if (0) {
                           			for (int k1=0;k1<Cindex.size();k1++) {
                                		if (Cindex.at(P-1)==Pindex.at(k1) && Cindex.at(k1)!=Cindex.at(i)) {
                                    		ctsisomer.at(0).at(k1)+="\\"; 
											if (0) ctsisomer.at(0).at(Cindex.at(P-1)-1)+="\\";
                                    		break;
                                		}
                            		}
								}

								if (0) {
									int ct=0;
                            		for (int k1=0;k1<Cindex.size();k1++) {
                                		if (Pindex.at(k1)==Cindex.at(i)) {
                                    		ct=rand()%2;
                                    		string cis_trans="";
                                    		if (ct) cis_trans="\\";
                                    		else cis_trans="/";

                                    		ctsisomer.at(0).at(k1)+=cis_trans;

											if (0) ctsisomer.at(1).at(Cindex.at(i)-1)+=cis_trans;

                                    		ct++;
                                    		ct=ct%2;
                                		}
                            		}
								}

								if (1) ctsisomer.at(1).at(Cindex.at(i)-1)+='/';
							}
							if (1) {
								if (!cas) {
									if (0) {
										if (1) {
											if (!flip) ctsisomer.at(0).at(Cindex.at(P-1)-1)+="/";
											else ctsisomer.at(0).at(Cindex.at(P-1)-1)+="\\";
										}
										if (1) {
											if (!flip) ctsisomer.at(1).at(Cindex.at(i)-1)+='/';
											else ctsisomer.at(1).at(Cindex.at(P-1)-1)+="\\";
										}

										flip++;
										flip=flip%2;
									}
									if (1) {
										ctsisomer.at(0).at(Cindex.at(P-1)-1)+="/";
										ctsisomer.at(1).at(Cindex.at(i)-1)+='/';
									}
								}
								else {
	                                if (1) {
    	                                int ct=0;
                                    	for (int k1=0;k1<Cindex.size();k1++) {
                                        	if (Pindex.at(k1)==Cindex.at(P-1) && Cindex.at(k1)!=Cindex.at(i)) {
                                            	ct=rand()%2;
                                            	string cis_trans="";
                                            	if (ct) cis_trans="\\";
                                            	else cis_trans="/";

                                            	ctsisomer.at(0).at(k1)+=cis_trans;

                                            	ct++;
                                            	ct=ct%2;
                                        	}
                                    	}
                                	}
									
								}
							}

							if (0) cout << hascistrans << " | " << y1 << " | " << y2 << " | Append pos: " << Cindex.at(P-1)-1 << " " << Cindex.at(i)-1 << endl;

						}

					}
				}
			}
		}
	}

	/*
	vector<bool> posi(Cindex.size(),0);
	if (1) {
		if (ctsisomer.size()==Cindex.size()) {
			for (int i=Cindex.size()-1;i>=0;i--) {
				if (Rindex.at(i)==2 && Mindex.at(i)==2) {
					int P=Pindex.at(i);
					if (P>0) {
						if (Mindex.at(P-1)==2 && Cindex.at(P-1)>0) {
							//if (ctsisomer.at(Cindex.at(P-1)-1)=="") {
							//	ctsisomer.clear();
							//	break;
							//}
							//else posi.push_back(Cindex.at(P-1)-1);
							posi.at(Cindex.at(P-1)-1)=1;
						}
						for (int k1=0;k1<Cindex.size();k1++) {
							if (Pindex.at(k1)==Cindex.at(i)) {
								//if (ctsisomer.at(k1)=="") {
								//	ctsisomer.clear();
								//	break;
								//}
								//else posi.push_back(k1);
								posi.at(k1)=1;
							}
						}
					}
				}
			}
			for (int i=0;i<posi.size();i++) {
				if (posi.at(i) && ctsisomer.at(i)!="") posi.at(i)=0;
			}
			
			for (int i=0;i<posi.size();i++) {
                if (posi.at(i)) break;
				else if (!posi.at(i) && i>=posi.size()-1) return;
            }
		}
		else ctsisomer.clear();

		if (!ctsisomer.size()) {
			ctsisomer.resize(Cindex.size(),"");
			for (int i=0;i<posi.size();i++) posi.at(i)=1;
		}
	}
	if (0) {
		ctsisomer.clear();
		ctsisomer.resize(Cindex.size(),"");
	}
	
	for (int i=lposi;i>=sposi;i--) { // C\C=C\C
		if (Rindex.at(i)==2 && Mindex.at(i)==2) {
			int P=Pindex.at(i);

			if (P>0) {
				if (Mindex.at(P-1)==2 && Cindex.at(P-1)>0) { // Mindex.at(P-1)==2 && Pindex.at(P-1)>0
					string cis_trans="";
					if (prob()>=0.5) cis_trans="/";
					else cis_trans="\\";
					
					if (posi.at(Cindex.at(P-1)-1)) {
						ctsisomer.at(Cindex.at(P-1)-1)+=cis_trans;
						posi.at(Cindex.at(P-1)-1)=0;
					}

		            if (1) {
        		        int ct=rand()%2;
						string keep="";
                		for (int k1=0;k1<Cindex.size();k1++) {
                    		if (Pindex.at(k1)==Cindex.at(i)) {
                        		cis_trans="";
                        		if (ct) cis_trans="/";
                        		else cis_trans="\\";

                                if (!posi.at(k1)) keep=ctsisomer.at(k1);
                        		if (posi.at(k1)) {
									if (keep=="/") ctsisomer.at(k1)="\\";
									else if (keep=="\\") ctsisomer.at(k1)="/";
									else ctsisomer.at(k1)+=cis_trans;
									posi.at(k1)=0;
								}

                        		ct++;
                        		ct=ct%2;
                    		}
                		}
						
            		}
				}
			}		
		}
	}

	vector<bool>().swap(posi);
	*/
}


void MOLECULE::reset() {
	clear();
	Bindex.reserve(Cindex.size());
	Bindex.resize(Cindex.size(),vector<int> (0));
	//Bindex.shrink_to_fit();
	atomsmi.reserve(Cindex.size());
	atomsmi.resize(Cindex.size(),vector<int> (0));
	//atomsmi.shrink_to_fit();

	for (int i=0;i<Cindex.size();i++) {
		int b = Mindex.at(i);
		int c = Pindex.at(i);
		int x = Rindex.at(i);

		for (int j=0;j<data->a.at(b).norder;j++) {
			Bindex.at(i).push_back(data->a.at(b).order.at(j));
		}

		bool u=0;
		if (c==0) {
			for (int j=0;j<data->a.at(b).name.length();j++) { // data.a[b].nbond //20191230
				atomsmi.at(i).push_back(j);
			}
		}
		else if (c>0) {
			int M = Mindex.at(c-1);
			int j=0;
			for (j=0;j<data->a.at(M).norder;j++) {
				if (Bindex.at(c-1).at(j) == x) {
					for (int n=0;n<data->a.at(b).norder;n++) {
						if (Bindex.at(i).at(n) == x) {
							Bindex.at(i).at(n) = 0;
							Bindex.at(c-1).at(j) = 0;
							u=1;
							break;
						}
					}
					if (u) break;
				}
			}
			int smindex=0;
			if (u) smindex=data->a.at(M).index+j*3;

			for (j=0;j<data->a.at(b).name.length();j++) {  // data.a[b].nbond //20191230
				int a=j+atomsmi.at(c-1).at(smindex)+1;		
				atomsmi.at(i).push_back(a);
			}

			for (j=0;j<i;j++) {
				for (int n=0;n<data->a.at(Mindex.at(j)).name.length();n++) { // data.a[Mindex.at(j)].nbond
					if (atomsmi.at(j).at(n) > atomsmi.at(c-1).at(smindex)) {
						atomsmi.at(j).at(n) = atomsmi.at(j).at(n)+data->a.at(b).name.length();  // + data.a[b].nbond
					}  
				}
			}

		}

	}


    if (1) { //delete "(-)" , "(=)" , "(#)" used for parent atom, ring, H atoms.
        if (0) {
            for (int i=0;i<atomsmi.size();i++) {
                cout << "orig ";
                for (int j=0;j<atomsmi.at(i).size();j++) {
                    cout << atomsmi.at(i).at(j) << " ";
                }
                cout << endl;
            }
        }
        vector<int> cid(0);
        vector<int> rpos(0);
        vector<int> apos(0);
        for (int i=0;i<Bindex.size();i++) {
            int M=Mindex.at(i);
            int R=Rindex.at(i);
            int P_quota=1;
            for (int j=0;j<data->a.at(M).norder;j++) {
                if (Bindex.at(i).at(j)) {
                    int pos=data->a.at(M).index+3*j-1;
                    int tmp=atomsmi.at(i).at(pos); //pos-1
                    cid.push_back(Cindex.at(i));
                    rpos.push_back(pos);
                    apos.push_back(tmp);
                }
                else if (!Bindex.at(i).at(j) && data->a.at(M).order.at(j)==R && P_quota>0) {
                    int pos=data->a.at(M).index+3*j-1;
                    int tmp=atomsmi.at(i).at(pos); //pos-1
                    cid.push_back(Cindex.at(i));
                    rpos.push_back(pos);
                    apos.push_back(tmp);
                    P_quota--;
                }
            }
        }
		
        cid.reserve(cid.size());
        rpos.reserve(rpos.size());
        apos.reserve(apos.size());
		
        for (int i=0;i<cid.size();i++) {
            for (int j=i+1;j<cid.size();j++) {
                if (apos.at(j)<apos.at(i)) {
                    swap(cid.at(i),cid.at(j));
                    swap(apos.at(i),apos.at(j));
                    swap(rpos.at(i),rpos.at(j));
                }
            }
        }
        for (int i=0;i<cid.size();i++) {
            int rtmp=rpos.at(i);
            int atmp=atomsmi.at(cid.at(i)-1).at(rpos.at(i));
            for (int j=0;j<3;j++) {
                atomsmi.at(cid.at(i)-1).erase(atomsmi.at(cid.at(i)-1).begin()+rpos.at(i));
            }
            for (int k=0;k<atomsmi.size();k++) {
                for (int l=0;l<atomsmi.at(k).size();l++) {
                    if (atomsmi.at(k).at(l)>atmp) atomsmi.at(k).at(l)-=3;
                }
            }
            for (int j=i+1;j<cid.size();j++) {
                if (cid.at(j)==cid.at(i) && apos.at(j)>apos.at(i)) rpos.at(j)-=3;
                if (apos.at(j)>apos.at(i)) apos.at(j)-=3;
            }
            if (0) {
                for (int k=0;k<atomsmi.size();k++) {
                    cout << Cindex.at(k) << " " << Mindex.at(k) << " QQ ";
                    for (int j=0;j<atomsmi.at(k).size();j++) {
                        cout << atomsmi.at(k).at(j) << " ";
                    }
                    cout << endl;
                }
            }
        }

        vector<int>().swap(apos);
        vector<int>().swap(cid);
        vector<int>().swap(rpos);
    }


    if (1) { //delete "-" of "....(-....."
        vector<int> cid(0);
        vector<int> rpos(0);
        vector<int> apos(0);

        for (int i=0;i<Bindex.size();i++) {
            int M=Mindex.at(i);
            int R=Rindex.at(i);
            int P_quota=1;
			int ct=0;
            for (int j=0;j<data->a.at(M).norder;j++) {
                if (Bindex.at(i).at(j)) ;
                else if (!Bindex.at(i).at(j) && data->a.at(M).order.at(j)==R && P_quota>0) P_quota--;
				else {
					if (data->a.at(M).order.at(j)==1) {
                    	int pos=data->a.at(M).index+3*ct;
                    	int tmp=atomsmi.at(i).at(pos);
                    	cid.push_back(Cindex.at(i));
                    	rpos.push_back(pos);
                    	apos.push_back(tmp);
					}
					ct++;
				}
            }
        }

        cid.reserve(cid.size());
        rpos.reserve(rpos.size());
        apos.reserve(apos.size());

        for (int i=0;i<cid.size();i++) {
            for (int j=i+1;j<cid.size();j++) {
                if (apos.at(j)<apos.at(i)) {
                    swap(cid.at(i),cid.at(j));
                    swap(apos.at(i),apos.at(j));
                    swap(rpos.at(i),rpos.at(j));
                }
            }
        }
        for (int i=0;i<cid.size();i++) {
            int rtmp=rpos.at(i);
            int atmp=atomsmi.at(cid.at(i)-1).at(rpos.at(i));

            atomsmi.at(cid.at(i)-1).erase(atomsmi.at(cid.at(i)-1).begin()+rpos.at(i));

            for (int k=0;k<atomsmi.size();k++) {
                for (int l=0;l<atomsmi.at(k).size();l++) {
                    if (atomsmi.at(k).at(l)>atmp) atomsmi.at(k).at(l)--;
                }
            }
            for (int j=i+1;j<cid.size();j++) {
                if (cid.at(j)==cid.at(i) && apos.at(j)>apos.at(i)) rpos.at(j)--;
                if (apos.at(j)>apos.at(i)) apos.at(j)--;
            }
            if (0) {
                for (int k=0;k<atomsmi.size();k++) {
                    cout << Cindex.at(k) << " " << Mindex.at(k) << " QQ ";
                    for (int j=0;j<atomsmi.at(k).size();j++) {
                        cout << atomsmi.at(k).at(j) << " ";
                    }
                    cout << endl;
                }
            }
        }

        vector<int>().swap(apos);
        vector<int>().swap(cid);
        vector<int>().swap(rpos);
    }


	if (1) {
		int m=0;
		vector<int> t(0);
		vector<int> s(0);

		for (int i=0;i<Cyindex.size();i++) {
			if (Cyindex.at(i).size()) {
				for (int j=0;j<Cyindex.at(i).size();j++) {
					if (Cyindex.at(i).at(j)>m) m=Cyindex.at(i).at(j);
				}
			}
		}
		if (m==0) {
			if_circle=0;
			//for (int i=0;i<Cyindex.size();i++) Cyindex.at(i)=0;
			Cyindex.resize(0,vector<int>(0));
			Cyindex.resize(Cindex.size(),vector<int>(0));
		}


		if (m) {
			vector<int> cyclic;
			cyclic.resize(m,0);
			for (int i=0;i<Mindex.size();i++) {
				if (Cyindex.at(i).size())  {
					for (int j=0;j<Cyindex.at(i).size();j++) {
						cyclic.at(Cyindex.at(i).at(j)-1)+=1;
					}
				}
			}
			for (int i=cyclic.size()-1;i>=0;i--) { // orig i=0;i<cyclic.size();i++ //20200725
				if (cyclic.at(i)==2) {
					for (int j=0;j<Mindex.size();j++) {
						if (1) { //20200509
							t.resize(0);

							if (Cyindex.at(j).size()) { // t: each ring numbering
                    			for (int k=0;k<Cyindex.at(j).size();k++) {
									t.push_back(Cyindex.at(j).at(k));
                    			}
							}	
						}

						bool go1=0;
						if (1) {
							for (int j1=0;j1<t.size();j1++) {
								if (t.at(j1)==(i+1)) { // orig t.at(j1)==s.at(j2) && s.at(j2)==(i+1)
									go1=1;
									break;
								}
							}
						}


						if (go1) { // orig go1 && go2 20200725 // ( (s[0]==t[0] || s[0]==t[1] || s[1]==t[0] || s[1]==t[1]) && (s[0]==(i+1) || s[1]==(i+1)) )
							bool u=0;
							for (int m1=0;m1<data->a.at(Mindex.at(j)).norder;m1++) {
								if (Bindex.at(j).at(m1)==Cybnd.at(i)) { // orig Bindex[j].at(m1)==1 && Bindex[n].at(x1)==1
									Bindex.at(j).at(m1)=0;
									u=1;
									break;
								}
							}
							if (!u) {
								u=1;

                            	if (Cyindex.at(j).size()) {
                                	for (int k=0;k<Cyindex.at(j).size();k++) {
                                    	if (Cyindex.at(j).at(k)==(i+1)) {
											Cyindex.at(j).erase(Cyindex.at(j).begin()+k);
											k--;
										}
                                	}

        							for (int k=0;k<Cyindex.at(j).size();k++) {
                                        if (Cyindex.at(j).at(k)>(i+1)) {
                                            Cyindex.at(j).at(k)--;
                                        }
        							}

                                    if (Cybnd.size()>=(i+1)) Cybnd.erase(Cybnd.begin()+i);
                                    if_circle--;

                            	}
							}
						}
					}
				} 
				else {
					for (int j=0;j<Mindex.size();j++) {
						if (1) { //20200509
							if (Cyindex.at(j).size()) {
                                for (int k=0;k<Cyindex.at(j).size();k++) {
									if (Cyindex.at(j).at(k)==(i+1)) {
										Cyindex.at(j).erase(Cyindex.at(j).begin()+k);
										k--;
									}
                                }

                                for (int k=0;k<Cyindex.at(j).size();k++) {
                                    if (Cyindex.at(j).at(k)>(i+1)) {
                                        Cyindex.at(j).at(k)--;
                                    }
                                }

                                if (Cybnd.size()>=(i+1)) Cybnd.erase(Cybnd.begin()+i);
                                if_circle--;
							}
						}
					}
				}
			}
		}


		if (m) {
			m=0;
			for (int i=0;i<Cyindex.size();i++) {
				for (int j=0;j<Cyindex.at(i).size();j++) {
					if (Cyindex.at(i).at(j)>m) {
						m=Cyindex.at(i).at(j);
					}
				}
			}
			if_circle=m;
		}
	}

    if (0) {
        for (int k=0;k<Bindex.size();k++) {
            cout << Cindex.at(k) << " " << Mindex.at(k) << " QQ ";
            for (int j=0;j<Bindex.at(k).size();j++) {
                cout << Bindex.at(k).at(j) << " ";
            }
            cout << endl;
        }
    }


	//del_unpaired_ring_no();

	if (para.protect) prct(); //20200101
	//chk_cistrans();

	return;
}


void MOLECULE::init() {
	int i=0,j=0,k=0;
	string tmp="";
	double x,y,z;

	//OBMol mol;
	//smi2gjf(mol);
	smi2gjf();

	if (atm!=NULL) delete [] atm;
	atm=NULL;
    atm = new DEATOM [natom];

    ifstream inf((para.smidir+"gg.txt").c_str());
    inf >> tmp >> ws;
	k=0;
    while (!inf.eof()) {
        tmp="";
        inf >> tmp >> ws;
        if (tmp=="V30") {
			string nam="";
			int ch=0;
            inf >> tmp >> ws;
            inf >> nam >> ws; //inf >> atm[k].name >> ws;
            inf >> tmp >> ws;
            inf >> tmp >> ws;
            inf >> tmp >> ws;
            inf >> tmp >> ws;
            do {
                inf >> tmp >> ws;
                if (tmp.length()>=4) {
                    if (tmp[0]=='C' && tmp[1]=='H' && tmp[2]=='G' && tmp[3]=='=') {
						string buf=tmp.substr(4);
                        ch=atoi(buf.c_str()); //atm[k].chg=atoi(buf.c_str());
                    }
                }
            } while (tmp!="M" && !inf.eof());
			//k++;

			if (nam!="H" && nam!="Xx") {
				atm[k].name=nam;
				atm[k].chg=ch;
				k++;
			}
        }
    }
	inf.close();

	//chir_and_stereo(mol);

	//dist = new double *[natom];
	connect = new int *[natom];
	//order= new int *[natom];
	for (i=0;i<natom;i++) {
		//order[i]=new int [natom];
		connect[i]=new int [natom];
		//dist[i]=new double [natom];
	}

	//for (i=0;i<natom;i++) {
	//	dist[i][i]=0.0;
	//}
	
	/*
	k=0;
	for (i=0;i<natom;i++) {
		if (atm[i].name != "H") {
			k++;
		}
	}
	*/
	for (i=0;i<natom;i++) { //i<k
		Cindex.push_back(i+1);
		//Cyindex.push_back(0);
	}
	Cyindex.resize(natom,vector<int>(0));

	if_circle=0;
	//cal_r();
	check_bnd(); //mol
	return;
}


void MOLECULE::smi2gjf() { //OBMol &mol
	if (0) {
		ofstream out((para.smidir+smiles+".smi").c_str());
		out << smiles << endl;
		out.close();
		//system(("ssh cluster '  /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -ogjf -O "+para.smidir+"\""+smiles+".gjf\" --gen3D --conformer --converge 1000000 -c  --nconf 100 ' ").c_str());
		system(("ssh cluster ' /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D -x3 --conformer --converge 100000 --score energy --nconf 100 -c  2> /dev/null ' ").c_str());
	}
	if (1) {
        //string a=smiles;
		string a=molesmi;
        for (int i=0;i<a.length();i++) {
        	if (a[i]=='/') a[i]='u';
        	if (a[i]=='\\') a[i]='d';
			if (a[i]=='*') a[i]='x';
		}
        
		stringstream nu("");
		if (1) {
			//stringstream ss(smiles);
			stringstream ss(molesmi);
			stringstream ss1("");
        	//ofstream out((para.smidir+smiles+".mol").c_str());
        	ofstream out((para.smidir+a+".mol").c_str());

        	OBConversion conv(&ss,&ss1);
        	if(conv.SetInAndOutFormats("SMI","MOL")) {
            	conv.AddOption("gen3D", OBConversion::GENOPTIONS);
				conv.AddOption("3", OBConversion::OUTOPTIONS);
				//--minimize --ff --steps 0
				//conv.AddOption("canonical", OBConversion::GENOPTIONS);
				//conv.AddOption("minimize", OBConversion::OUTOPTIONS);
				//conv.AddOption("ff", OBConversion::OUTOPTIONS,"uff");
				//conv.AddOption("step", OBConversion::OUTOPTIONS,"1");
				//conv.AddOption("align", OBConversion::GENOPTIONS);
				//conv.AddOption("c", OBConversion::OUTOPTIONS);
				//a1=time(NULL);
				//obErrorLog.StopLogging();
				obErrorLog.SetOutputLevel(obMessageLevel::obError);
				//obErrorLog.SetOutputLevel(obMessageLevel::obWarning);
				obErrorLog.SetOutputStream(&cout);
            	conv.Convert();
        	}
			out << ss1.str();
			out.close();
			nu << ss1.str();
		}
	    //system(("(grep -A10000 \"BEGIN ATOM\" "+para.smidir+"\""+smiles+".mol\" | grep -B10000 \"END ATOM\" |  tail -n +2 | tac | tail -n +2 | tac > "+para.smidir+"gg.txt) 2> /dev/null").c_str());
    	//system(("(grep -A10000 \"BEGIN ATOM\" "+para.smidir+"\""+a+".mol\" | grep -B10000 \"END ATOM\" |  tail -n +2 | tac | tail -n +2 | tac > "+para.smidir+"gg.txt) 2> /dev/null").c_str());

        OBMol *mol=NULL;
        mol = new OBMol [1];
        OBConversion conv(&nu);
        conv.SetInFormat("MOL");
        conv.Read(&mol[0]);
        OBAtom *ar=NULL;
		int ct=0,Du=0;;
    	FOR_ATOMS_OF_MOL(ar,mol[0]) {
			ct++;
			if (ar->GetAtomicNum()==0) Du++;
        }
		//natom=ct;
		natom=mol[0].NumHvyAtoms()-Du;

        mol[0].Clear();

        if (mol!=NULL) {
            delete [] mol;
            mol=NULL;
        }


		stringstream ss("");
		//ss << natom;
		ss << ct;

		system(("(grep -A"+ss.str()+" \"BEGIN ATOM\" "+para.smidir+"\""+a+".mol\" | tail -n +2 > "+para.smidir+"gg.txt) 2> /dev/null").c_str());
	}
	
	return;
}

void MOLECULE::chir_and_stereo(OBMol &mol) {
	if (1) {
        vector<int> num(0);
        FOR_ATOMS_OF_MOL(ar, mol) {
            //cout << ar->GetIdx() << " " << ar->GetType() << " " << ar->GetFormalCharge() << " " << ar->GetSpinMultiplicity() << endl;
            if (ar->GetFormalCharge()!=0) atm[ar->GetIdx()-1].chg=ar->GetFormalCharge();
            if (ar->IsChiral()==1 && ar->GetAtomicNum()==6) num.push_back(ar->GetIdx());
        }
        for (int i=0;i<num.size();i++) {
			OBAtom *ar=NULL;
            ar = mol.GetAtom(num.at(i));
            if (ar->IsChiral() && ar->GetAtomicNum()==6) {
                OBChiralData* cd=(OBChiralData*)ar->GetData(OBGenericDataType::ChiralData);
                if (!cd) { //if no Chiral Data Set, need to make one!
                    cd=new OBChiralData;
                    ar->SetData(cd);
                }
                if (ar->GetHvyValence()==4) {
                    OBAtom *nbr=NULL;
                    vector<unsigned int> nbr_atms;
                    vector<OBEdgeBase*>::iterator i;
                    for (nbr=ar->BeginNbrAtom(i);nbr;nbr=ar->NextNbrAtom(i)) nbr_atms.push_back(nbr->GetIdx());
                    sort(nbr_atms.begin(),nbr_atms.end());
                    cd->SetAtom4Refs(nbr_atms,output);   // This saves the output atom4refs calculated above
                }
            //If it has co-ordinates then we can calculate the signed volume which stores the atom4refs used in the calculation and its result. This can then be compared to the output order and the corretchirality worked out. At the moment this overides input chirality, rather than checking for a conflict with input chirality (which might be a usefull check)*
                double vol=0;
                if (mol.HasNonZeroCoords()) {
                    vol=CalcSignedVolume(mol,ar);
                    if (vol>0.0) ar->SetClockwiseStereo();
                    else if(vol<0.0) ar->SetAntiClockwiseStereo();
                    CorrectChirality(mol,ar,calcvolume,output);
                }
                else { //If no co-ords, then use the atom4refs defined by the input format*
                    CorrectChirality(mol,ar); // will set the stereochem based on input/output atom4refs
                }

                if (ar->IsClockwise()) atm[ar->GetIdx()-1].chirality=2;
                else if (ar->IsAntiClockwise()) atm[ar->GetIdx()-1].chirality=1;
            }
			ar=NULL;
        }


	}
}


void MOLECULE::check_bnd() { //OBMol &mol
	int i,j,k,chgg=0; //sum[natom]
	double r,tol=1.15;

	for (i=0;i<natom;i++) {
		connect[i][i]=0;
		for (j=i+1;j<natom;j++) connect[i][j]=connect[j][i]=0;
	}

	if (1) {
		if (0) {
			system(("ssh cluster ' /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D -x3 --weighted --conformer --converge 100000 --score energy --nconf 100 -c  2> /dev/null ' ").c_str());	
		}

		stringstream nu("");
		
	    //string a=smiles;
		string a=molesmi;
	    for (int i=0;i<a.length();i++) {
	        if (a[i]=='/') a[i]='u';
	        if (a[i]=='\\') a[i]='d';
			if (a[i]=='*') a[i]='x';
		}
            
		//k=rd_nps(mol);
		k=rd_nps();

		if (0) {
			//int nonH=mol.NumHvyAtoms();
        	int nonH=0;
        	//for (k=0;k<natom;k++) {
            //	if (atm[k].name!="H") nonH++;
        	//}
			if (1) nonH=natom;
		
			if (0) {
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
							break;
						}
					}
				}
			}


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
			}

            if (1) {
                bool sw=0;
                for (k=1;k<nonH;k++) { // k<natom
                    if (atm[k].name!="H") {
                        bool aall0=1;
                        //bndpos1=999999999; //int bndpos1=999999999;
                        if (1) {
                            for (int k1=0;k1<=k;k1++) {
                                if (connect[k1][k]!=0 || connect[k-k1][k]!=0) { // orig connect[k1][k]!=0
                                    aall0=0;
                                    break;
                                }
                            }
                        }
                        if (aall0) {
                            for (int k1=k+1;k1<nonH;k1++) { //k1<natom
                                swap(atm[k],atm[k1]);

                                for (int k2=0;k2<natom;k2++) {
                                    if (connect[k][k2]!=connect[k1][k2]) swap(connect[k][k2],connect[k1][k2]);
                                }
                                for (int k2=0;k2<natom;k2++) {
                                    if (connect[k2][k]!=connect[k2][k1]) swap(connect[k2][k],connect[k2][k1]);
                                }
                            }

                            if (1) {
                                for (int k2=0;k2<nonH-1;k2++) {
                                    if (connect[k2][nonH-1]==1 && atm[k2].name=="C") {
                                        if (atm[k2].chirality==1) atm[k2].chirality=2;
                                        else if (atm[k2].chirality==2) atm[k2].chirality=1;
                                    }
                                }
                            }

							if (0) {
								bool g1=0,g2=0;
								for (int k2=0;k2<nonH-1;k2++) {
									g1=0;
									g2=0;
									if (connect[k2][nonH-1]==1) {
										g1=1;
										g2=0;
										for (int k3=0;k3<k2;k3++) {
											if (connect[k3][k2]==2) {
												g2=1;
												g1=1;
												if (1) {
		                                			if (atm[k2].cistrans[0]=="") {
        		                            			if (atm[k2].cistrans[1]=="/") atm[k2].cistrans[1]="\\";
                		                    			else if (atm[k2].cistrans[1]=="\\") atm[k2].cistrans[1]="/";
                        		        			}
												}
												if (0) {
                                					if (atm[k2].cistrans[1]=="") {
                                    					if (atm[k2].cistrans[0]=="/") atm[k2].cistrans[0]="\\";
                                    					else if (atm[k2].cistrans[0]=="\\") atm[k2].cistrans[0]="/";
                                					}
												}
												if (g1 && g2) break;
											}
										}
										if (g1 && g2) break;
									}
								}
							}

							k--;
                        }
                    }
                }
            }


		}

	}

    if (0) {
        cout << setfill(' ');
        cout << smiles << " " << molesmi << endl;
        cout << setw(2) << " " << " ";
        for (i=0;i<natom;i++) {
            cout << setw(2) << atm[i].name << " ";
        }
        cout << endl;
        for (i=0;i<natom;i++) {
            cout << setw(2) << atm[i].name << " ";
            for (j=0;j<natom;j++) {
                cout << setw(2) << connect[i][j] << " ";
            }
			cout << setw(2) << atm[i].nH << " ";
            cout << endl;
        }
    }


	//for (i=0;i<natom;i++) {
	//	sum[i]=0;
	//	for (j=0;j<natom;j++) sum[i]+=connect[i][j];
	//}

	//for (i=0;i<natom;i++) {
	//	for (j=0;j<natom;j++) order[i][j]=connect[i][j];
	//}
	//k=0;

	//for (i=0;i<natom;i++) {
	//	for (j=0;j<natom;j++) {
	//		connect[i][j]=order[i][j];
	//	}
	//}
	return;
}


void MOLECULE::smi2mds() {
	int i,j,n,t,m,u;
	vector<int> k;

	//cout << "SMI2COD" << endl;

	//empty();
	clear(); // 20190722
	ctsisomer.resize(0,vector<string>(0));
	ctsisomer.resize(2,vector<string>(0));

	for (i=0;i<natom;i++) {
		//cout << atm[i].name << " | " << atm[i].chg << " | " << atm[i].nH << endl;
		k.clear();
		if (atm[i].name=="H" && atm[i].chg==0) continue; 
		if (atm[i].name!="H") {  //20200823
			//ctsisomer.push_back(atm[i].cistrans);
			ctsisomer.at(0).push_back(atm[i].cistrans[0]);
			ctsisomer.at(1).push_back(atm[i].cistrans[1]);
		}

		if (atm[i].name == "C") { 
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

			if (k.size() == 4) {
				//cout << "Y" << endl;
				if (atm[i].chg==0 && atm[i].chirality==1) Mindex.push_back(68);
				else if (atm[i].chg==0 && atm[i].chirality==0) Mindex.push_back(1);
				else if (atm[i].chg==0 && atm[i].chirality==2) Mindex.push_back(67);
				//cout << Mindex.size() << endl;
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
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

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
			bool usepibnd=0;
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
				if (connect[i][j]==2) usepibnd=1;
			}
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

			if (k.size()==4) {
				if (atm[i].chg==1) Mindex.push_back(15);
			}
			if (k.size()==3) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==4 && atm[i].chg==1) Mindex.push_back(16);
				else if (t==3 && atm[i].chg==0) Mindex.push_back(7);
			}
			if (k.size()==2) {
				t=0;
				for (n=0;n<k.size();n++) t+=k[n];
				if (t==3 && atm[i].chg==0 && !usepibnd) Mindex.push_back(8);
				else if (t==3 && atm[i].chg==0 && usepibnd) Mindex.push_back(71);
				else if (t==2 && atm[i].chg==-1) Mindex.push_back(35);
				else if (t==2 && atm[i].chg==0) Mindex.push_back(7);
			}
			if (k.size()==1) {
				if (atm[i].chg==0) Mindex.push_back(9);
			}
		}
		else if (atm[i].name=="P") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

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
				else if (t==4 && atm[i].chg==0) Mindex.push_back(69);
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
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

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
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

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
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

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
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

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
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

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
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

			if (k.size()==4) {
				if (atm[i].chg==-1) Mindex.push_back(44);
			}
		}
		else if (atm[i].name=="In") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

			if (k.size()==4) {
				if (atm[i].chg==3) Mindex.push_back(57);
			}
		}
		else if (atm[i].name=="Ga") {
			for (j=0;j<natom;j++) {
				if (connect[i][j] != 0) k.push_back(connect[i][j]);
			}
			if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);

			if (k.size()==4) {
				if (atm[i].chg==3) Mindex.push_back(64);
			}
		}
        //else if (atm[i].name=="Xx") {
        //    for (j=0;j<natom;j++) {
        //        if (connect[i][j] != 0) k.push_back(connect[i][j]);
        //    }
        //    if (1) for (int u1=0;u1<atm[i].nH;u1++) k.push_back(1);
		//
        //    if (k.size()==1) {
        //        if (atm[i].chg==0) Mindex.push_back(70);
        //    }
        //}
	}
	//Rindex.push_back(-1);
	Rindex.push_back(0);
	Pindex.push_back(0);


	int nonH=natom;
	if (1) {
        for (i=1;i<nonH;i++) { //i=1
			t=0;
            if (atm[i].name!="H") {
				for (int c1=3;c1>=1;c1--) {
                	for (j=0;j<=i;j++) { // j=i;j>=0;j--
                    	if (connect[i][j]==c1 && atm[j].name!="H") { //order[i][j]==m
                        	Rindex.push_back(connect[i][j]); //order[i][j]
                        	Pindex.push_back(Cindex.at(j));
                        	connect[j][i]=connect[i][j]=0; //order[j][i]=order[i][j]=0;
							t=1;
                        	break;
                    	}
                	}
					if (t) break;
				}
            }
        }
	}

	
	
	n=1;
	for (i=0;i<nonH;i++) {
		if (atm[i].name!="H") {
			for (j=i+1;j<nonH;j++) {
				if (atm[j].name!="H" && connect[i][j]!=0) { //order[i][j]!=0
					//cout << j << " " << natom << " " << Cyindex.size() << endl;
					//Cyindex.at(i)=Cyindex.at(i)*10+n;
					//Cyindex.at(j)=Cyindex.at(j)*10+n;
					Cyindex.at(j).push_back(n);
					Cyindex.at(i).push_back(n);
					Cybnd.resize(n,0);
					Cybnd.at(n-1)=connect[i][j];

					connect[i][j]=connect[j][i]=0; // order[i][j]=order[j][i]=0;
					if_circle=n;
					n++;
				}
			}
		}
	}

	if (0) {
		for (i=0;i<Cindex.size();i++) {
            bool go=0;
            if (Mindex.at(i)==2) go=1;
            if (Mindex.at(i)==8) go=1;
            if (Mindex.at(i)==16) go=1;
			
			if (go && Rindex.at(i)==2) {
				int P=Pindex.at(i);

            	bool go1=0;
            	if (Mindex.at(P-1)==2) go1=1;
            	if (Mindex.at(P-1)==8) go1=1;
            	if (Mindex.at(P-1)==16) go1=1;

				if (go1) {
					//swap(ctsisomer.at(1).at(P-1),ctsisomer.at(1).at(i));  //20201105
					//swap(ctsisomer.at(0).at(P-1),ctsisomer.at(0).at(i));  //20201105

					if (P>0) {
						if (ctsisomer.at(1).at(P-1)!="" && ctsisomer.at(0).at(P-1)=="") {
							swap(ctsisomer.at(1).at(P-1),ctsisomer.at(0).at(P-1));

							if (0) { //20201105
								if (ctsisomer.at(0).at(P-1)=="/") ctsisomer.at(0).at(P-1)="\\";
								if (ctsisomer.at(0).at(P-1)=="\\") ctsisomer.at(0).at(P-1)="/";
							}
						}
					}
					if (ctsisomer.at(0).at(i)!="" && ctsisomer.at(1).at(i)=="") {
						swap(ctsisomer.at(1).at(i),ctsisomer.at(0).at(i));

                        if (0) { //20201105
                            if (ctsisomer.at(1).at(i)=="/") ctsisomer.at(1).at(i)="\\";
                            if (ctsisomer.at(1).at(i)=="\\") ctsisomer.at(1).at(i)="/";
                        }
					}

				}

			}
		}
	}

	if (Cindex.size()>2 && 1) {
		if (ctsisomer.at(0).at(0)!="" && ctsisomer.at(1).at(1)!="" && Rindex.at(1)==2 && Pindex.at(1)==1) {
			int ct1=0,ct2=0;
        	for (int k1=0;k1<Cindex.size();k1++) {
            	if (Pindex.at(k1)==Cindex.at(0) && Cindex.at(k1)!=Cindex.at(1)) {
					if (ct1==0) {
						//ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(0).at(0)+ctsisomer.at(0).at(Cindex.at(k1)-1);
						if (1) {
                        	if (ctsisomer.at(0).at(0)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
                        	else if (ctsisomer.at(0).at(0)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
						}
			
						if (0) {
                            if (ctsisomer.at(0).at(0)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
                            else if (ctsisomer.at(0).at(0)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
						}
					}
					else if (ct1==1) {
						ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(0).at(0)+ctsisomer.at(0).at(Cindex.at(k1)-1);
						//if (ctsisomer.at(0).at(0)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
						//else if (ctsisomer.at(0).at(0)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
					} 
					//if (ct1==1) ctsisomer.at(0).at(0)=ctsisomer.at(1).at(0)="";
					ct1++;
				}
				if (Pindex.at(k1)==Cindex.at(1)) {
					if (ct2==0) {
						//ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(1).at(1)+ctsisomer.at(0).at(Cindex.at(k1)-1);
						if (1) {
                        	if (ctsisomer.at(1).at(1)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
                        	else if (ctsisomer.at(1).at(1)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
						}
                        if (0) {
                            if (ctsisomer.at(1).at(1)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
                            else if (ctsisomer.at(1).at(1)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
                        }
					}
					else if (ct2==1) {
						ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(1).at(1)+ctsisomer.at(0).at(Cindex.at(k1)-1);
                    	//if (ctsisomer.at(1).at(1)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
                    	//else if (ctsisomer.at(1).at(1)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
					}
					//if (ct2==1) ctsisomer.at(1).at(1)=ctsisomer.at(0).at(1)="";
					ct2++;
				}
			}
			if (1) {
				ctsisomer.at(0).at(0)="";
				ctsisomer.at(1).at(1)="";
			}
		}
	}

	return;
}


void MOLECULE::report() {
	int i=0,j=0;
	if (para.protect) prct();

	//string a=smiles;
	string a=molesmi;
	for (int k=0;k<a.length();k++) {
		if (a[k]=='/') a[k]='u';
		if (a[k]=='\\') a[k]='d';
		if (a[k]=='*') a[k]='x';
	}

	//ofstream outf((para.smidir+smiles+".mds").c_str());
	ofstream outf((para.smidir+a+".mds").c_str());
	outf<<"SMILES "<<left<<molesmi<<" ";
    outf<<endl;
	outf<<"natom "<<setw(4)<<left<<setfill(' ')<<Cindex.size()<<" ";
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
	for (i=0;i<Cyindex.size();i++) {
		if (!Cyindex.at(i).size()) outf<<setw(4)<<left<<setfill(' ')<<"0"<<" ";
		else {
			stringstream buf("");
			for (j=0;j<Cyindex.at(i).size();j++) {
				if (j<Cyindex.at(i).size()-1) buf<<Cyindex.at(i).at(j)<<",";
				else buf<<Cyindex.at(i).at(j);

				//if (j<Cyindex.at(i).size()-1) outf<<left<<Cyindex.at(i).at(j)<<",";
				//else outf<<left<<Cyindex.at(i).at(j)<<" ";
			}
			outf<<setw(4)<<left<<setfill(' ')<<buf.str()<<" ";
		}
	}
	outf<<endl;
	//outf<<"ctsisomer ";
    //for (i=0;i<ctsisomer.size();i++) {
	//	if (ctsisomer.at(i)!="") outf<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(i)<<" ";
	//	else outf<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
	//}
    //outf<<endl;
    outf<<"ctsisomer_start ";
    for (i=0;i<ctsisomer.at(0).size();i++) {
      	if (ctsisomer.at(0).at(i)!="") outf<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(0).at(i)<<" ";
      	else outf<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    outf<<endl;
    outf<<"ctsisomer_end ";
    for (i=0;i<ctsisomer.at(1).size();i++) {
    	if (ctsisomer.at(1).at(i)!="") outf<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(1).at(i)<<" ";
    	else outf<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    outf<<endl;

	if (para.protect) {
		outf<<"protection ";
		//for (i=0;i<Cindex.size();i++) outf<<j<<" ";
		for (i=0;i<protect.size();i++) outf<<setw(4)<<left<<setfill(' ')<<protect.at(i)<<" ";
		outf<<endl;
	}

	outf<<"if_circle "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;

    outf<<"Cybnd ";
    for (i=0;i<Cybnd.size();i++) outf<<setw(4)<<left<<setfill(' ')<<Cybnd.at(i)<<" ";
	outf<<endl;

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
	//if(dist==NULL && connect==NULL && order==NULL && atm==NULL) return;
	if (connect==NULL && atm==NULL) return;
	for (i=0;i<natom;i++) {
		//if (dist!=NULL) delete [] dist[i]; //if(dist[i]!=NULL)
		if (connect!=NULL) delete [] connect[i]; //if(connect[i]!=NULL)
		//if (order!=NULL) delete [] order[i]; //if(order[i]!=NULL)
	}
	//if(dist!=NULL) delete [] dist;
	if(connect!=NULL) delete [] connect;
	//if(order!=NULL) delete [] order;
	if(atm!=NULL) delete [] atm;
	//dist=NULL;
	connect=NULL;
	//order=NULL;
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
    outs<<"SMILES "<<left<<molesmi<<" ";
    outs<<endl;
	outs<<"natom "<<setw(4)<<left<<setfill(' ')<<Cindex.size()<<" ";
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
	for (i=0;i<Cyindex.size();i++) {
        if (!Cyindex.at(i).size()) outs<<setw(4)<<left<<setfill(' ')<<"0"<<" ";
        else {
			stringstream buf("");
            for (j=0;j<Cyindex.at(i).size();j++) {
                if (j<Cyindex.at(i).size()-1) buf<<Cyindex.at(i).at(j)<<",";
                else buf<<Cyindex.at(i).at(j);

                //if (j<Cyindex.at(i).size()-1) outs<<left<<Cyindex.at(i).at(j)<<",";
                //else outs<<left<<Cyindex.at(i).at(j)<<" ";
            }
			outs<<setw(4)<<left<<setfill(' ')<<buf.str()<<" ";
        }
	}
	outs<<endl;
    //outs<<"ctsisomer ";
    //for (i=0;i<ctsisomer.size();i++) {
    //    if (ctsisomer.at(i)!="") outs<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(i)<<" ";
    //    else outs<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    //}
	//outs<<endl;
    outs<<"ctsisomer_start ";
    for (i=0;i<ctsisomer.at(0).size();i++) {
        if (ctsisomer.at(0).at(i)!="") outs<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(0).at(i)<<" ";
        else outs<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    outs<<endl;
    outs<<"ctsisomer_end ";
    for (i=0;i<ctsisomer.at(1).size();i++) {
        if (ctsisomer.at(1).at(i)!="") outs<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(1).at(i)<<" ";
        else outs<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    outs<<endl;


	outs<<"if_circle "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;


    outs<<"Cybnd ";
    for (i=0;i<Cybnd.size();i++) outs<<setw(4)<<left<<setfill(' ')<<Cybnd.at(i)<<" ";
	outs<<endl;


	return;
}

int MOLECULE::combination(MOLECULE &B, int k,int p,int bnd) {
	if (k>=Cindex.size() || k<0) return 0;
	if (p>=B.Cindex.size() || p<0) return 0;
	if (!data->a[Mindex.at(k)].bd[bnd-1]) return 0;
	if (!B.data->a[B.Mindex.at(p)].bd[bnd-1]) return 0;

	int nu=B.Cindex.at(p);

    if (1) {
        bool totbndchk=0;
        bool bndchk=0;

        int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        vector<int> ordcount(6,0);
        if (Rindex.at(k)>0) {
			bndsum+=Rindex.at(k);
			ordcount.at(Rindex.at(k)-1)+=1;
		}
        if (Cyindex.at(k).size()) {
            for (int k1=0;k1<Cyindex.at(k).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(k).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(k).at(k1)-1)-1)++;
            }

        }
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Pindex.at(k1)==Cindex.at(k)) {
                bndsum+=Rindex.at(k1);
                ordcount.at(Rindex.at(k1)-1)+=1;
            }
        }
		bndsum+=bnd;
		ordcount.at(bnd-1)+=1;

        int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(k)].order.at(k2);
        }
		/*
        if (bndsum<=maxbnd) totbndchk=1;
        else {
			cout << "1st molecule " << Cindex.at(k) << " | Required bonds: " << bndsum << " | Available bonds: " << maxbnd << endl;
            totbndchk=0;
            return 0;
        }
		*/

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[Mindex.at(k)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
				//cout << "1st molecule " << Cindex.at(k) << " | Bond order: " << (k2+1) << " | Required number: " << ordcount.at(k2) << " | Available number: " << data.a[Mindex.at(k)].bd[k2] << endl;
                bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(k)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

    }

    if (1) {
        bool totbndchk=0;
        bool bndchk=0;

        int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        vector<int> ordcount(6,0);
        if (B.Rindex.at(p)>0) {
			bndsum+=B.Rindex.at(p);
			ordcount.at(B.Rindex.at(p)-1)+=1;
		}
        if (B.Cyindex.at(p).size()) {
            for (int k1=0;k1<B.Cyindex.at(p).size();k1++) {
                bndsum+=B.Cybnd.at(B.Cyindex.at(p).at(k1)-1);
                ordcount.at(B.Cybnd.at(B.Cyindex.at(p).at(k1)-1)-1)++;
            }
        }
        for (int k1=0;k1<B.Cindex.size();k1++) {
            if (B.Pindex.at(k1)==B.Cindex.at(p)) {
                bndsum+=B.Rindex.at(k1);
                ordcount.at(B.Rindex.at(k1)-1)+=1;
            }
        }
        bndsum+=bnd;
        ordcount.at(bnd-1)+=1;

        int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[B.Mindex.at(p)].order.at(k2);
        }
		/*
        if (bndsum<=maxbnd) totbndchk=1;
        else {
			cout << "2nd molecule " << B.Cindex.at(p) << " | Required bonds: " << bndsum << " | Available bonds: " << maxbnd << endl;
            totbndchk=0;
            return 0;
        }
		*/

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[B.Mindex.at(p)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
				//cout << "2nd molecule " << B.Cindex.at(p) << " | Bond order: " << (k2+1) << " | Required number: " << ordcount.at(k2) << " | Available number: " << data.a[B.Mindex.at(p)].bd[k2]  << endl;
                bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[B.Mindex.at(p)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

    }


	vector<int> tmpP(B.Pindex.size(),0);
	vector<int> tmpC(B.Cindex.size(),0);
	vector<int> tmpM(B.Mindex.size(),0);
	vector<int> tmpR(B.Rindex.size(),0);
	//vector<int> tmpCy(B.Cyindex.size(),0);
	vector< vector<int> > tmpCy(B.Cyindex.size(),vector<int>(0));
	//vector<string> tmpCT(B.ctsisomer.size(),"");
	vector< vector<string> > tmpCT(2,vector<string>(B.ctsisomer.at(0).size(),""));
	vector<int> tmpPr(B.protect.size(),0);
	for (int k2=0;k2<B.Pindex.size();k2++) tmpP.at(k2)=B.Pindex.at(k2);
	for (int k2=0;k2<B.Cindex.size();k2++) tmpC.at(k2)=B.Cindex.at(k2);
	for (int k2=0;k2<B.Mindex.size();k2++) tmpM.at(k2)=B.Mindex.at(k2);
	for (int k2=0;k2<B.Rindex.size();k2++) tmpR.at(k2)=B.Rindex.at(k2);
	//for (int k2=0;k2<B.ctsisomer.size();k2++) tmpCT.at(k2)=B.ctsisomer.at(k2);
	for (int k2=0;k2<B.ctsisomer.at(0).size();k2++) tmpCT.at(0).at(k2)=B.ctsisomer.at(0).at(k2);
	for (int k2=0;k2<B.ctsisomer.at(1).size();k2++) tmpCT.at(1).at(k2)=B.ctsisomer.at(1).at(k2);
	//for (int k2=0;k2<B.Cyindex.size();k2++) tmpCy.at(k2)=B.Cyindex.at(k2);
	for (int k2=0;k2<B.Cyindex.size();k2++) {
		for (int k3=0;k3<B.Cyindex.at(k2).size();k3++) tmpCy.at(k2).push_back(B.Cyindex.at(k2).at(k3));
	}
	if (para.protect) for (int k2=0;k2<B.protect.size();k2++) tmpPr.at(k2)=B.protect.at(k2);

	if (p!=0) {
		tmpP.at(p)=0;
    	int compa=B.Cindex.at(p);
    	int cont=0;
    	while (compa>=1 && cont<=B.Cindex.size()) {
        	//if (1) C_ringmember.push_back(compa);

        	if (compa>=1) {
            	if (B.Pindex.at(compa-1)>=1) {
					tmpP.at(B.Pindex.at(compa-1)-1)=B.Cindex.at(compa-1);
					tmpR.at(B.Pindex.at(compa-1)-1)=B.Rindex.at(compa-1);
					/*
					for (int i=B.Pindex.at(compa-1)-1;i<tmpC.size()-1;i++) {
                		swap(tmpP.at(i),tmpP.at(i+1));
                		swap(tmpC.at(i),tmpC.at(i+1));
                		swap(tmpM.at(i),tmpM.at(i+1));
                		swap(tmpR.at(i),tmpR.at(i+1));
                		swap(tmpCy.at(i),tmpCy.at(i+1));
                		//swap(tmpCT.at(i),tmpCT.at(i+1));
                		swap(tmpCT.at(0).at(i),tmpCT.at(0).at(i+1));
                		swap(tmpCT.at(1).at(i),tmpCT.at(1).at(i+1));
                		if (para.protect) swap(tmpPr.at(i),tmpPr.at(i+1));
    				}
					*/

					compa=B.Cindex.at(B.Pindex.at(compa-1)-1);
				}
            	else break;
            	cont++;
        	}
        	else break;
    	}
		//tmpR.at(p)=-1;
		tmpR.at(p)=0;
	}

    if (0) {
        cout << "loop 2-1" << endl;
        for (int k2=0;k2<tmpP.size();k2++) cout << tmpP.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpC.size();k2++) cout << tmpC.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpM.size();k2++) cout << tmpM.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpR.size();k2++) cout << tmpR.at(k2) << " ";
        cout << endl;
    	for (int k2=0;k2<tmpCy.size();k2++) {
			if (!tmpCy.at(k2).size()) cout << "0" << " ";
            else {
				for (int k3=0;k3<tmpCy.at(k2).size();k3++) {
					if (k3<tmpCy.at(k2).size()-1) cout << tmpCy.at(k2).at(k3) << ","; 
					else cout << tmpCy.at(k2).at(k3) << " ";
				}
			}
		}
        cout << endl;
        if (para.protect) {
            for (int k2=0;k2<tmpPr.size();k2++) cout << tmpPr.at(k2) << " ";
            cout << endl;
        }
    }

	if (p!=0) {
		swap(tmpP.at(p),tmpP.at(0));
		swap(tmpC.at(p),tmpC.at(0));
		swap(tmpM.at(p),tmpM.at(0));
		swap(tmpR.at(p),tmpR.at(0));
		swap(tmpCy.at(p),tmpCy.at(0));
		//swap(tmpCT.at(p),tmpCT.at(0));
		swap(tmpCT.at(0).at(p),tmpCT.at(0).at(0));
		swap(tmpCT.at(1).at(p),tmpCT.at(1).at(0));
		if (para.protect) swap(tmpPr.at(p),tmpPr.at(0));
	}
	
    if (0) {
        cout << "loop 2-2" << endl;
        for (int k2=0;k2<tmpP.size();k2++) cout << tmpP.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpC.size();k2++) cout << tmpC.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpM.size();k2++) cout << tmpM.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpR.size();k2++) cout << tmpR.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpCy.size();k2++) {
            if (!tmpCy.at(k2).size()) cout << "0" << " ";
            else {
                for (int k3=0;k3<tmpCy.at(k2).size();k3++) {
                    if (k3<tmpCy.at(k2).size()-1) cout << tmpCy.at(k2).at(k3) << ",";
                    else cout << tmpCy.at(k2).at(k3) << " ";
                }
            }
		}
        cout << endl;
        if (para.protect) {
            for (int k2=0;k2<tmpPr.size();k2++) cout << tmpPr.at(k2) << " ";
            cout << endl;
        }
    }

	for (int i=0;i<tmpC.size();i++) {
		for (int j=i+1;j<tmpC.size();j++) {
			if (tmpP.at(i)==tmpC.at(j)) {
				swap(tmpP.at(j),tmpP.at(i));
                swap(tmpC.at(j),tmpC.at(i));
                swap(tmpM.at(j),tmpM.at(i));
                swap(tmpR.at(j),tmpR.at(i));
                swap(tmpCy.at(j),tmpCy.at(i));
				//swap(tmpCT.at(j),tmpCT.at(i));
				swap(tmpCT.at(0).at(j),tmpCT.at(0).at(i));
				swap(tmpCT.at(1).at(j),tmpCT.at(1).at(i));
				if (para.protect) swap(tmpPr.at(j),tmpPr.at(i));

				i--;
				break;
			}
		}
	}
	
    for (int i=0;i<tmpC.size();i++) {
        for (int j=i+1;j<tmpC.size();j++) {
			if (tmpP.at(i)==tmpC.at(j)) { // && tmpC.at(i)==tmpP.at(j)
			    vector<int>().swap(tmpP);
			    vector<int>().swap(tmpC);
			    vector<int>().swap(tmpM);
			    vector<int>().swap(tmpR);
			    vector< vector<int> >().swap(tmpCy);
			    vector<int>().swap(tmpPr);
				//vector<string>().swap(tmpCT);
				vector< vector<string> >().swap(tmpCT);

                return 0;
			}
        }
    }
	
    if (0) {
        cout << "loop 2-3" << endl;
        for (int k2=0;k2<tmpP.size();k2++) cout << tmpP.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpC.size();k2++) cout << tmpC.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpM.size();k2++) cout << tmpM.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpR.size();k2++) cout << tmpR.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpCy.size();k2++) {
            if (!tmpCy.at(k2).size()) cout << "0" << " ";
            else {
                for (int k3=0;k3<tmpCy.at(k2).size();k3++) {
                    if (k3<tmpCy.at(k2).size()-1) cout << tmpCy.at(k2).at(k3) << ",";
                    else cout << tmpCy.at(k2).at(k3) << " ";
                }
            }
        }
        cout << endl;
        if (para.protect) {
            for (int k2=0;k2<tmpPr.size();k2++) cout << tmpPr.at(k2) << " ";
            cout << endl;
        }
    }

	if (0) {
	for (int i=0;i<tmpC.size();i++) {
		tmpC.at(i)+=tmpC.size()+1;
		tmpP.at(i)+=tmpC.size()+1;
	}
	}

    if (0) {
        cout << "loop 2-4" << endl;
        for (int k2=0;k2<tmpP.size();k2++) cout << tmpP.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpC.size();k2++) cout << tmpC.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpM.size();k2++) cout << tmpM.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpR.size();k2++) cout << tmpR.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpCy.size();k2++) {
            if (!tmpCy.at(k2).size()) cout << "0" << " ";
            else {
                for (int k3=0;k3<tmpCy.at(k2).size();k3++) {
                    if (k3<tmpCy.at(k2).size()-1) cout << tmpCy.at(k2).at(k3) << ",";
                    else cout << tmpCy.at(k2).at(k3) << " ";
                }
            }
        }
        cout << endl;
        if (para.protect) {
            for (int k2=0;k2<tmpPr.size();k2++) cout << tmpPr.at(k2) << " ";
            cout << endl;
        }
    }

	if (0) {
    for (int i=0;i<tmpC.size();i++) {
        int j=tmpC.at(i);
        tmpC.at(i)=(i+1);
        for (int k1=0;k1<tmpP.size();k1++) {
            if (tmpP.at(k1)==j) tmpP.at(k1)=tmpC.at(i);
        }
    }
	if (tmpP.size()) tmpP.at(0)=0;
	}

    for (int i=0;i<tmpC.size();i++) {
        int j=tmpC.at(i);
        tmpC.at(i)=(i+1)+tmpC.size();
        for (int k1=0;k1<tmpC.size();k1++) { 
            if (tmpP.at(k1)==j) tmpP.at(k1)=tmpC.at(i);
        }
    }

	for (int i=0;i<tmpC.size();i++) {
		tmpC.at(i)-=tmpC.size();
		if (i!=0) tmpP.at(i)-=tmpC.size();
	}

    if (0) {
        cout << "loop 2-5" << endl;
        for (int k2=0;k2<tmpP.size();k2++) cout << tmpP.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpC.size();k2++) cout << tmpC.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpM.size();k2++) cout << tmpM.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpR.size();k2++) cout << tmpR.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpCy.size();k2++) {
            if (!tmpCy.at(k2).size()) cout << "0" << " ";
            else {
                for (int k3=0;k3<tmpCy.at(k2).size();k3++) {
                    if (k3<tmpCy.at(k2).size()-1) cout << tmpCy.at(k2).at(k3) << ",";
                    else cout << tmpCy.at(k2).at(k3) << " ";
                }
            }
        }
        cout << endl;
        if (para.protect) {
            for (int k2=0;k2<tmpPr.size();k2++) cout << tmpPr.at(k2) << " ";
            cout << endl;
        }
    }


    if (1) {
        for (int i=1;i<tmpC.size();i++) {
            if (tmpR.at(i)==2 && tmpP.at(i)>0) {
                if (tmpCT.at(1).at(tmpP.at(i)-1)!="" && tmpCT.at(0).at(i)!="") {
                    swap(tmpCT.at(0).at(i),tmpCT.at(1).at(i));
                    swap(tmpCT.at(0).at(tmpP.at(i)-1),tmpCT.at(1).at(tmpP.at(i)-1));
                }
            }
        }
    }


	int Cshift=Cindex.size();
    for (int i=0;i<tmpC.size();i++) {
        tmpC.at(i)+=Cshift;  //at
        if (i==0) { // orig i==p
        	tmpP.at(i)=Cindex.at(k); //at
        	tmpR.at(i)=bnd;
        }
        else tmpP.at(i)+=Cshift; //at
    }

    if (0) {
        cout << "loop 2-6" << endl;
        for (int k2=0;k2<tmpP.size();k2++) cout << tmpP.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpC.size();k2++) cout << tmpC.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpM.size();k2++) cout << tmpM.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpR.size();k2++) cout << tmpR.at(k2) << " ";
        cout << endl;
        for (int k2=0;k2<tmpCy.size();k2++) {
            if (!tmpCy.at(k2).size()) cout << "0" << " ";
            else {
                for (int k3=0;k3<tmpCy.at(k2).size();k3++) {
                    if (k3<tmpCy.at(k2).size()-1) cout << tmpCy.at(k2).at(k3) << ",";
                    else cout << tmpCy.at(k2).at(k3) << " ";
                }
            }
        }
        cout << endl;
        if (para.protect) {
            for (int k2=0;k2<tmpPr.size();k2++) cout << tmpPr.at(k2) << " ";
            cout << endl;
        }
    }

	Cybnd.resize(Cybnd.size()+B.if_circle,0);
    for (int k2=0;k2<tmpCy.size();k2++) {
        if (tmpCy.at(k2).size()) {
            for (int k3=0;k3<tmpCy.at(k2).size();k3++) {
				Cybnd.at(tmpCy.at(k2).at(k3)-1+if_circle)=B.Cybnd.at(tmpCy.at(k2).at(k3)-1);
				tmpCy.at(k2).at(k3)+=if_circle;
            }
        }
    }
	if_circle=B.if_circle+if_circle;

    for (int i=0;i<tmpC.size();i++) {
        Pindex.push_back(tmpP.at(i));
        Cindex.push_back(tmpC.at(i));
        Rindex.push_back(tmpR.at(i));
        //Cyindex.push_back(tmpCy.at(i));
		Cyindex.resize(Cyindex.size()+1,vector<int>(0));
		for (int j=0;j<tmpCy.at(i).size();j++) Cyindex.at(Cyindex.size()-1).push_back(tmpCy.at(i).at(j));
        Mindex.push_back(tmpM.at(i));
		//ctsisomer.push_back(tmpCT.at(i));
		ctsisomer.at(0).push_back(tmpCT.at(0).at(i));
		ctsisomer.at(1).push_back(tmpCT.at(1).at(i));
        if (para.protect) protect.push_back(tmpPr.at(i));
    }

	//print();

	vector<int>().swap(tmpP);
	vector<int>().swap(tmpC);
	vector<int>().swap(tmpM);
	vector<int>().swap(tmpR);
	vector< vector<int> >().swap(tmpCy);
	//vector<string>().swap(tmpCT);
	vector< vector<string> >().swap(tmpCT);
	vector<int>().swap(tmpPr);

	//reset();
	//chk_cistrans(); //20200806
    mds2smi();
	//if (0) cout << Cindex.at(k) << " " << nu << " " << molesmi << endl;


    return 1;

}

/*
int MOLECULE::change_stereo() {
    if (!Cindex.size()) return 0;
    vector< vector<int> > posi(2,vector<int>(0)); // 0 for transform to 67, 1 for transform to 68
    for (int i=0;i<Cindex.size();i++) {
        //if (Mindex.at(i)==67 && Mindex.at(i)==68) {
        //    int P=Pindex.at(i);
        //    if (P>0) {
        //        if (Mindex.at(P-1)==2 && ctsisomer.at(P-1)!="") posi.push_back(i);
        //    }
        //}
        if (Mindex.at(i)==67) {
            posi.at(1).push_back(i);
        }
        else if (Mindex.at(i)==68) {
            posi.at(0).push_back(i);
        }
        else if (Mindex.at(i)==1) {
            if (prob()<=0.5) posi.at(0).push_back(i);
            else posi.at(1).push_back(i);
        }

    }
    if (!posi.at(0).size() && !posi.at(1).size()) return 0;

    for (int i=0;i<posi.size();i++) {
        for (int j=0;j<posi.at(i).size();j++) {
            if (prob()<=0.33) change_stereo(posi.at(i).at(j),i);
        }
    }
    return 1;
}
*/

/*
int MOLECULE::change_stereo(int pos,int he) {
    //if (Mindex.at(pos)==67) Mindex.at(pos)=68;
    //else if (Mindex.at(pos)==68) Mindex.at(pos)=67;
    if (he==1) Mindex.at(pos)=68;
    else if (he==0) Mindex.at(pos)=67;

    return 1;
}
*/

int MOLECULE::change_chirality(int pos,int spec) {
	bool go=0;
	//if (Mindex.at(pos)==68 || Mindex.at(pos)==67 || Mindex.at(pos)==1) go=1;
	if (Mindex.at(pos)==68 || Mindex.at(pos)==67) go=1;

	if (go && spec!=Mindex.at(pos)) {
		if (spec==-1) {
    		if (Mindex.at(pos)==67) Mindex.at(pos)=68;
    		else if (Mindex.at(pos)==68) Mindex.at(pos)=67;
			//else if (Mindex.at(pos)==1) {
			//	if (prob()<=0.5) Mindex.at(pos)=68;
			//	else Mindex.at(pos)=67;
			//}
		}
		else {
			Mindex.at(pos)=spec;
		}
	}
	else return 0;
    return 1;
}

/*
int MOLECULE::change_cistrans() {
    if (!Cindex.size()) return 0;
    vector< vector<int> > posi(2,vector<int>(0));
    int ct=0;
    if (ctsisomer.size()==2) {
        for (int i=0;i<ctsisomer.size();i++) {
            for (int j=0;j<ctsisomer.at(i).size();j++) {
                if (ctsisomer.at(i).at(j)!="") {
                    posi.at(i).push_back(j);
                }
            }
        }
        if (!posi.at(0).size() && !posi.at(1).size()) return 0;

        for (int i=0;i<posi.size();i++) {
            for (int j=0;j<posi.at(i).size();j++) {
                if (prob()<=0.33) {
                    change_cistrans(posi.at(i).at(j),i);
                }
            }
        }
    }
    else return 0;
    return 1;
}
*/

int MOLECULE::change_cistrans(int pos,int he) {
    if (ctsisomer.size()==2) {
        if (ctsisomer.at(0).size()==Cindex.size() && ctsisomer.at(0).size()==ctsisomer.at(1).size()) {
            if (ctsisomer.at(he).at(pos)!="") {
                if (ctsisomer.at(he).at(pos)=="/") ctsisomer.at(he).at(pos)="\\";
                else if (ctsisomer.at(he).at(pos)=="\\") ctsisomer.at(he).at(pos)="/";
            }
            else return 0;
        }
        else return 0;
    }
    else return 0;
    return 1;
}

void MOLECULE::wipe() {
	vector<int>().swap(Cindex);
	vector<int>().swap(Pindex);
	vector< vector<int> >().swap(Cyindex);
	vector<int>().swap(Cybnd);
	vector<int>().swap(Rindex);
	vector<int>().swap(Mindex);
	vector<int>().swap(protect);
	
	

	//int i;
	//for (i=0;i<1024;i++) {
	//	vector<int>().swap(Bindex[i]);
	//	vector<int>().swap(atomsmi[i]);
	//}

	vector< vector<int> >().swap(Bindex);
	vector< vector<int> >().swap(atomsmi);
	//vector<string>().swap(ctsisomer);
	vector< vector<string> >().swap(ctsisomer);

	//delete [] atm;
	//atm=NULL;
	clean();

	//atm.resize(0);
	//vector<DEATOM>().swap(atm);
	molesmi.clear();
	smiles.clear();
	return;
}


void MOLECULE::input(bool ct_on) {
	if (smiles=="" && molesmi!="null") smiles=molesmi;
	int tmpchg=chg;
	empty();
	clean(); //20200701
	init(); //cout << "init end" << endl;
	smi2mds(); //cout << "smi2cod end" << endl;
	report(); //cout << "report end" << endl;
	//read(smiles); //cout << "read end" << endl;
	read(molesmi);  //cout << "read(molesmi) end" << endl;
	//reset();
	//chk_cistrans(); //20200806
	mds2smi(ct_on); //cout << "mds2smi end" << endl;
	canonicalize_SMILES(); //20200706

	//print();
	return;
}


void MOLECULE::empty() {
	Cindex.resize(0);
	Pindex.resize(0);
	Mindex.resize(0);
	Rindex.resize(0);
	Cyindex.resize(0);
	Cybnd.resize(0);
	protect.resize(0);
	//ctsisomer.resize(0);
	ctsisomer.resize(0,vector<string> (0));
	Bindex.resize(0,vector<int> (0));
	atomsmi.resize(0,vector<int> (0));
	//for (int i=0;i<1024;i++) {
	//	Bindex[i].resize(0);
	//	atomsmi[i].resize(0);
	//}
	return;
}


int MOLECULE::rd_nps() { //OBMol &mol
	//system(("ssh cluster ' /home/software/programs/openbabel-openbabel-2-4-1/build/bin/obabel -ismi "+para.smidir+"\""+smiles+".smi\" -omol -O "+para.smidir+"\""+smiles+".mol\" --gen3D --weighted --conformer --converge 100000 --score energy --nconf 100 -c  2> /dev/null ' ").c_str());
	//int i=0,j=0,k=0,l=0,s=0,m,n,p,q,x=0,y=0;

    stringstream nu("");

    if (1) {
        //string a=smiles;
		string a=molesmi;
        for (int i=0;i<a.length();i++) {
            if (a[i]=='/') a[i]='u';
            if (a[i]=='\\') a[i]='d';
			if (a[i]=='*') a[i]='x';
        }

        if (1) {
            //stringstream ss(smiles);
			stringstream ss(molesmi);
            stringstream ss1("");
            //ofstream out((para.smidir+smiles+".mol").c_str());
            ofstream out((para.smidir+a+".mol").c_str());

            OBConversion conv(&ss,&ss1);
            if(conv.SetInAndOutFormats("SMI","MOL")) {
                conv.AddOption("gen3D", OBConversion::GENOPTIONS);
                conv.AddOption("3", OBConversion::OUTOPTIONS);
                //--minimize --ff --steps 0
                //conv.AddOption("canonical", OBConversion::GENOPTIONS);
                //conv.AddOption("minimize", OBConversion::OUTOPTIONS);
                //conv.AddOption("ff", OBConversion::OUTOPTIONS,"uff");
                //conv.AddOption("step", OBConversion::OUTOPTIONS,"1");
                //conv.AddOption("align", OBConversion::GENOPTIONS);
                //obErrorLog.StopLogging();
                obErrorLog.SetOutputLevel(obMessageLevel::obError);
                //obErrorLog.SetOutputLevel(obMessageLevel::obWarning);
                obErrorLog.SetOutputStream(&cout);
                conv.Convert();
                //cout << smiles << endl;
            }
            out << ss1.str();
            out.close();
            nu << ss1.str();
        }
    }



	if (1) {
        OBMol *mol=NULL;
        mol = new OBMol [1];
        OBConversion conv(&nu);
        conv.SetInFormat("MOL");
        conv.Read(&mol[0]);
        OBStereoFacade facade(&mol[0]);
		
		vector<int> Hpos(0);
        FOR_ATOMS_OF_MOL(ar, mol[0]) {
            //atm[ar->GetId()].id=ar->GetId();
			if (ar->GetAtomicNum()==1 || ar->GetAtomicNum()==0) Hpos.push_back(ar->GetId());
        }

		FOR_BONDS_OF_MOL(bond, mol[0]) {
			OBAtom *bea=bond->GetBeginAtom();
			OBAtom *ena=bond->GetEndAtom();
			int bnd=bond->GetBondOrder();
			bool bea_not_HorDu=1;
			if (bea->GetAtomicNum()==1 || bea->GetAtomicNum()==0) bea_not_HorDu=0;
			bool ena_not_HorDu=1;
			if (ena->GetAtomicNum()==1 || ena->GetAtomicNum()==0) ena_not_HorDu=0;
			
			if (bea!=NULL && ena!=NULL && bea_not_HorDu && ena_not_HorDu) {
				int beanum=bea->GetId();
				int enanum=ena->GetId();
				for (int i=0;i<Hpos.size();i++) {
					if (bea->GetId()>Hpos.at(i)) beanum--;
					if (ena->GetId()>Hpos.at(i)) enanum--;
				}
				connect[beanum][enanum]=connect[enanum][beanum]=bnd;
			}
			else if (bea!=NULL && ena!=NULL && bea_not_HorDu && !ena_not_HorDu) {
                int beanum=bea->GetId();
                for (int i=0;i<Hpos.size();i++) {
                    if (bea->GetId()>Hpos.at(i)) beanum--;
                }
				if (ena->GetAtomicNum()!=0) atm[beanum].nH++;
			}
            else if (bea!=NULL && ena!=NULL && !bea_not_HorDu && ena_not_HorDu) {
                int enanum=ena->GetId();
                for (int i=0;i<Hpos.size();i++) {
                    if (ena->GetId()>Hpos.at(i)) enanum--;
                }
                if (bea->GetAtomicNum()!=0) atm[enanum].nH++;
            }
			bea==NULL;
			ena==NULL;
		}

		vector<int> ringmem(0),ctcenter(0);
		if (0) {
			FOR_RINGS_OF_MOL(ring, mol[0]) {
				for (int i=0;i<ring->_path.size();i++) ringmem.push_back(ring->_path.at(i));
			}
		}
        //cout << "ringmem ";
        //for (int i=0;i<ringmem.size();i++) cout << ringmem.at(i) << " ";
        //cout << endl;

		int flip=0;
		FOR_BONDS_OF_MOL(bond, mol[0]){
			if (bond->GetBondOrder()==2 && facade.HasCisTransStereo(bond->GetId())) {
				OBAtom *bea=bond->GetBeginAtom();
            	OBAtom *ena=bond->GetEndAtom();

				bool fine=1;
				if (ena->GetAtomicNum()==6 && bea->GetAtomicNum()==6) fine=1;
				if (ena->GetAtomicNum()==7 && bea->GetAtomicNum()==6) fine=1;
				if (ena->GetAtomicNum()==6 && bea->GetAtomicNum()==7) fine=1;
				if (ena->GetAtomicNum()==7 && bea->GetAtomicNum()==7) fine=1;

				if (fine) {
					bool go=1;
					if (0) {
						bool a1=0,a2=0;
						for (int i=0;i<ringmem.size();i++) {
							if (bea->GetIdx()==ringmem.at(i)) {
								a1=1;
								break;
							}				
						}
						for (int i=0;i<ringmem.size();i++) {
							if (ena->GetIdx()==ringmem.at(i)) {
								a2=1;
								break;
							}					
						}
						if (a1 && a2) go=0;
					}

					if (go) {
	                	OBCisTransStereo *ct1=NULL;
	                	ct1=facade.GetCisTransStereo(bond->GetId());
                		OBCisTransStereo::Config A=ct1->GetConfig(OBStereo::ShapeU);

						if (go) {
							int ct=0;
							if (1) {
								vector<string> tmpct(4,"");
                            	for (int g=0;g<A.refs.size();g++) {
									int tm=A.refs.at(g);
									for (int u1=0;u1<Hpos.size();u1++) {
										if (A.refs.at(g)>Hpos.at(u1)) tm--;
									}

                                	if (tm<natom) { // cis-trans (heavy atom)
										if (g==0) tmpct.at(g)="/";
										else if (g==1) tmpct.at(g)="\\";
										else if (g==2) tmpct.at(g)="/";
										else if (g==3) tmpct.at(g)="\\";
									}
                                	//else if (A.refs.at(g)>=natom && A.refs.at(g)<mol[0].NumHvyAtoms()) { // cis-trans (hydrogen atom)
                                    //	if (g==0) tmpct.at(g)="u";
                                    //	else if (g==1) tmpct.at(g)="d";
                                    //	else if (g==2) tmpct.at(g)="u";
                                    //	else if (g==3) tmpct.at(g)="d";
                                	//}
                                	else if (tm>=natom) { // lone pair
                                        OBAtom *ty=NULL;
                                        ty=mol[0].GetAtomById(A.refs.at(g));
										if (ty->GetAtomicNum()==0) {
                                    		//tmpct.at(g)="*";
                                    		if (g==0) tmpct.at(g)="/";
                                    		else if (g==1) tmpct.at(g)="\\";
                                    		else if (g==2) tmpct.at(g)="/";
                                    		else if (g==3) tmpct.at(g)="\\";
										}
										ty=NULL;
                                	}
								}

								if (tmpct.at(0)!="") tmpct.at(1)="";
								if (tmpct.at(1)!="") tmpct.at(0)="";
								if (tmpct.at(2)!="") tmpct.at(3)="";
								if (tmpct.at(3)!="") tmpct.at(2)="";
								
                                int beanum=bea->GetId();
                                for (int u1=0;u1<Hpos.size();u1++) {
                                    if (bea->GetId()>Hpos.at(u1)) beanum--;
                                }
                                int enanum=ena->GetId();
                                for (int u1=0;u1<Hpos.size();u1++) {
                                    if (ena->GetId()>Hpos.at(u1)) enanum--;
                                }


								if (tmpct.at(0)!="" && tmpct.at(1)=="") atm[beanum].cistrans[0]=tmpct.at(0);
								else if (tmpct.at(0)=="" && tmpct.at(1)!="") atm[beanum].cistrans[0]=tmpct.at(1);
								if (tmpct.at(2)!="" && tmpct.at(3)=="") atm[enanum].cistrans[1]=tmpct.at(2);
                            	else if (tmpct.at(2)=="" && tmpct.at(3)!="") atm[enanum].cistrans[1]=tmpct.at(3);
							}

							if (0) {
								for (int g=0;g<A.refs.size();g++) {	
									int tm=A.refs.at(g);
									for (int u1=0;u1<Hpos.size();u1++) {
										if (A.refs.at(g)>Hpos.at(u1)) tm--;
									}
									if (tm<natom) {
										if (ct) {
                                        	if (g==1 && atm[tm].name!="H" && atm[bea->GetId()].cistrans[0]=="" && atm[bea->GetId()].cistrans[1]=="") {
												OBAtom *ty=NULL;
												ty=mol[0].GetAtomById(A.refs.at(g));
												if (bea->IsConnected(ty)) {
													int beanum=bea->GetId();
													for (int u1=0;u1<Hpos.size();u1++) {
														if (bea->GetId()>Hpos.at(u1)) beanum--;
													}
													atm[beanum].cistrans[0]+="\\";
												}
												ty=NULL;
											}
                                        	if (g==3 && atm[tm].name!="H" && atm[ena->GetId()].cistrans[0]=="" && atm[ena->GetId()].cistrans[1]=="") {
												OBAtom *ty=NULL;
                                            	ty=mol[0].GetAtomById(A.refs.at(g));
												if (ena->IsConnected(ty)) {
                                                	int enanum=ena->GetId();
                                                	for (int u1=0;u1<Hpos.size();u1++) {
                                                    	if (ena->GetId()>Hpos.at(u1)) enanum--;
                                                	}
                                                	atm[enanum].cistrans[1]+="\\";
												}
												ty=NULL;
											}
										}
										else {
                                    		if (g==0 && atm[tm].name!="H" && atm[bea->GetId()].cistrans[0]=="" && atm[bea->GetId()].cistrans[1]=="") {
                                        		OBAtom *ty=NULL;
                                        		ty=mol[0].GetAtomById(A.refs.at(g));
												if (bea->IsConnected(ty)) {
                                            		int beanum=bea->GetId();
                                            		for (int u1=0;u1<Hpos.size();u1++) {
                                                		if (bea->GetId()>Hpos.at(u1)) beanum--;
                                            		}
                                            		atm[beanum].cistrans[0]+="/";
												}
												ty=NULL;
											}
                                    		if (g==2 && atm[tm].name!="H" && atm[ena->GetId()].cistrans[0]=="" && atm[ena->GetId()].cistrans[1]=="") {
                                        		OBAtom *ty=NULL;
                                        		ty=mol[0].GetAtomById(A.refs.at(g));
												if (ena->IsConnected(ty)) {
                                        			int enanum=ena->GetId();
                                            		for (int u1=0;u1<Hpos.size();u1++) {
                                                		if (ena->GetId()>Hpos.at(u1)) enanum--;
                                            		}
                                            		atm[enanum].cistrans[1]+="/";
												}
												ty=NULL;
											}
										}
									}
								}
							}
							
							ct++;
							ct=ct%2;
						}

						flip++;
						flip=flip%2;
					}
				}
			}
		}

    	if (1) {
        	FOR_ATOMS_OF_MOL(atom, mol[0]) {
            	if (atom->GetAtomicNum()==6 && facade.HasTetrahedralStereo(atom->GetId())) {
                	OBTetrahedralStereo *tr=NULL;
                	tr=facade.GetTetrahedralStereo(atom->GetId());
                	bool go=1;
                	OBTetrahedralStereo::Config A=tr->GetConfig(OBStereo::Clockwise,OBStereo::ViewFrom);

                	bool flip=0;
                	if (go) {
                    	for (int i=0;i<Hpos.size();i++) {
                        	if (A.refs.at(1)==Hpos.at(i)) {
                            	flip=1;
                            	break;
                        	}
                    	}
                	}

                	if (go) {
                    	int tmp=A.center;
                    	for (int i=0;i<Hpos.size();i++) {
                        	if (A.center>Hpos.at(i)) tmp--;
                    	}
                    	if (!flip) atm[tmp].chirality=2;
						else atm[tmp].chirality=1;
                	}
                	if (!go) {
                    	go=1;
                    	OBTetrahedralStereo::Config B=tr->GetConfig(OBStereo::AntiClockwise,OBStereo::ViewFrom);
                    	if (go) {
                        	for (int i=0;i<Hpos.size();i++) {
                            	if (B.refs.at(1)==Hpos.at(i)) {
                                	flip=1;
                                	break;
                            	}
                        	}
                    	}

                    	if (go) {
                        	int tmp=B.center;
                        	for (int i=0;i<Hpos.size();i++) {
                            	if (B.center>Hpos.at(i)) tmp--;
                        	}
                        	if (!flip) atm[tmp].chirality=1;
							else atm[tmp].chirality=2;
                    	}
                	}

            	}
        	}
    	}
	}


	return 1;
}


int MOLECULE::addition(int pt, int id, int bnd, bool cistrans) {
	if (pt<0 || pt>=Cindex.size()) return 0;
	if (bnd>3 || bnd<=0) return 0;
	if (!data->a[Mindex.at(pt)].bd[bnd-1]) return 0;
	if (!data->a[id].bd[bnd-1]) return 0;

	int chgg=0;
	for (int i=0;i<Cindex.size();i++) chgg+=data->a[Mindex.at(i)].chg;
	if (data->a[id].chg*chgg<0) return 0;

	bool isnotend=0;

    // Ensure the bonds of atom P(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
    if (1) {
        bool totbndchk=0;
        bool bndchk=0;

        int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        vector<int> ordcount(6,0);
        if (Rindex.at(pt)>=1) {
			bndsum+=Rindex.at(pt);
			ordcount.at(Rindex.at(pt)-1)+=1;
		}
        if (Cyindex.at(pt).size()) {
            for (int k1=0;k1<Cyindex.at(pt).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(pt).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(pt).at(k1)-1)-1)++;
            }
        }
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Pindex.at(k1)==Cindex.at(pt)) {
				isnotend=1;
                bndsum+=Rindex.at(k1);
                ordcount.at(Rindex.at(k1)-1)+=1;
            }
        }
		bndsum+=bnd;
		ordcount.at(bnd-1)+=1;

        int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(pt)].order.at(k2);
        }
        /*
        if (bndsum<=maxbnd) totbndchk=1;
        else {
            totbndchk=0;
            return 0;
        }
        */

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[Mindex.at(pt)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(pt)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

    }

	if (0) {
	if (pt>=Cindex.size()-1) {
		Pindex.push_back(Cindex.at(pt));
		Cindex.push_back(Cindex.size()+1);
		Mindex.push_back(id);
		Rindex.push_back(bnd);
		//Cyindex.push_back(0);
		Cyindex.resize(Cyindex.size()+1,vector<int>(0));
		//ctsisomer.push_back("");
		if (para.protect) protect.push_back(0);
		ctsisomer.at(0).push_back("");
        ctsisomer.at(1).push_back("");
	}
	else {
		print();

		for (int i=0;i<Cindex.size();i++) {
			int tmp=Cindex.at(i);
			Cindex.at(i)+=(Cindex.size()+2);
			for (int j=i+1;j<Cindex.size();j++) {
				if (Pindex.at(j)==tmp) Pindex.at(j)=Cindex.at(i);
			}
		}

		print();

        Pindex.insert(Pindex.begin()+pt+1,Cindex.at(pt));
        Cindex.insert(Cindex.begin()+pt+1,Cindex.size()+1);
        Mindex.insert(Mindex.begin()+pt+1,id);
        Rindex.insert(Rindex.begin()+pt+1,bnd);
        Cyindex.insert(Cyindex.begin()+pt+1,vector<int>(0));
		if (para.protect) protect.insert(protect.begin()+pt+1,0);
        ctsisomer.at(0).insert(ctsisomer.at(0).begin()+pt+1,"");
        ctsisomer.at(1).insert(ctsisomer.at(1).begin()+pt+1,"");

		print();

        for (int i=0;i<Cindex.size();i++) {
            int tmp=Cindex.at(i);
            Cindex.at(i)=(i+1);
            for (int j=i+1;j<Cindex.size();j++) {
                if (Pindex.at(j)==tmp) Pindex.at(j)=Cindex.at(i);
            }
        }

		print();
	}
	}


	Pindex.push_back(Cindex.at(pt));
    Cindex.push_back(Cindex.size()+1);
    Mindex.push_back(id);
    Rindex.push_back(bnd);
    //Cyindex.push_back(0);
    Cyindex.resize(Cyindex.size()+1,vector<int>(0));
    //ctsisomer.push_back("");
    ctsisomer.at(0).push_back("");
    ctsisomer.at(1).push_back("");

	bool fl1=0,fl2=0;
	if (0) {
    	if (Pindex.at(pt)>0) if (Mindex.at(Pindex.at(pt)-1)==2 || Mindex.at(Pindex.at(pt)-1)==7 || Mindex.at(Pindex.at(pt)-1)==16) fl1=1;
    	if (Mindex.at(pt)==2 || Mindex.at(pt)==7 || Mindex.at(pt)==16) fl2=1;
	}
	if (1) {
		if (Pindex.at(pt)>0) fl1=1;
		fl2=1;
	}

	if (Rindex.at(pt)==2 && bnd==1 && isnotend && Pindex.at(pt)>0 && fl1 && fl2) {
		if (ctsisomer.at(0).at(Pindex.at(pt)-1)=="" && ctsisomer.at(1).at(Pindex.at(pt)-1)=="") {
			ctsisomer.at(0).at(Pindex.at(pt)-1)="/";
			ctsisomer.at(1).at(Pindex.at(pt)-1)="";
		}
        if (ctsisomer.at(0).at(Cindex.at(pt)-1)=="" && ctsisomer.at(1).at(Cindex.at(pt)-1)=="") {
            ctsisomer.at(0).at(Cindex.at(pt)-1)="";
            if (!cistrans) ctsisomer.at(1).at(Cindex.at(pt)-1)="\\";
			else ctsisomer.at(1).at(Cindex.at(pt)-1)="/";
        }
	}
	else if (Rindex.at(pt)==2 && bnd==1 && !isnotend && Pindex.at(pt)>0 && fl1 && fl2) {
        if (ctsisomer.at(0).at(Pindex.at(pt)-1)=="" && ctsisomer.at(1).at(Pindex.at(pt)-1)=="") {
            ctsisomer.at(0).at(Pindex.at(pt)-1)="/";
            ctsisomer.at(1).at(Pindex.at(pt)-1)="";
        }
        if (ctsisomer.at(0).at(Cindex.at(pt)-1)=="" && ctsisomer.at(1).at(Cindex.at(pt)-1)=="") {
            ctsisomer.at(0).at(Cindex.at(pt)-1)="";
            if (!cistrans) ctsisomer.at(1).at(Cindex.at(pt)-1)="/";
			else ctsisomer.at(1).at(Cindex.at(pt)-1)="\\";
        }
	}

	if (0) {
	bool fl1=0,fl2=0;
	if (id==2 || id==7 || id==16) fl1=1;
	if (Mindex.at(pt)==2 || Mindex.at(pt)==7 || Mindex.at(pt)==16) fl2=1;

	if (bnd==2 && ctsisomer.at(0).at(pt)!="" && fl1 && fl2) {
		ctsisomer.at(0).push_back("");
		ctsisomer.at(1).push_back("/");
	}
	else if (bnd==2 && ctsisomer.at(0).at(pt)=="" && fl1 && fl2) {
		ctsisomer.at(0).at(pt)="/";
		ctsisomer.at(1).at(pt)="";
		ctsisomer.at(0).push_back("");
		ctsisomer.at(1).push_back("/");
	}
	else {
		ctsisomer.at(0).push_back("");
		ctsisomer.at(1).push_back("");
	}
	}

	if (para.protect) protect.push_back(0);
	

	//reset();
	//chk_cistrans(); //20200806
	mds2smi();

	/*
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
	*/
	return 1;
}

int MOLECULE::insertion(int n,int id,int bnd2par,int bnd2des,bool cistrans) {
    if (n<=0 || n>=Cindex.size()) return 0; // n: bnd pos
    int M,P,C,R;
    M=Mindex.at(n);
    P=Pindex.at(n);
    C=Cindex.at(n);
    R=Rindex.at(n);

    if (para.protect) {
        if (protect.at(n)) return 0;
        if (P>0) if (protect.at(P-1)) return 0;
    }
    if (bnd2par>3 || bnd2par<=0) return 0;
    if (bnd2des>3 || bnd2des<=0) return 0;
    //if (data->a[id].chg!=data->a[Mindex.at(n)].chg) return 0;
    //else if (data->a[id].chg==data->a[Mindex.at(n)].chg && data->a[id].chg!=0) {
    //    if (data->a[id].index<=2 && data->a[Mindex.at(n)].index>=7) return 0;
    //    if (data->a[id].index>=7 && data->a[Mindex.at(n)].index<=7) return 0;
    //}

    int chgg=0;
    for (int i=0;i<Cindex.size();i++) chgg+=data->a[Mindex.at(i)].chg;
    if (data->a[id].chg*chgg<0) return 0;
    else if (data->a[id].chg*chgg>=0 && data->a[id].chg!=0) {
        if (data->a[id].index<=2 && data->a[Mindex.at(n)].index>=7) return 0;
        if (data->a[id].index>=7 && data->a[Mindex.at(n)].index<=7) return 0;
    }

    /*
    if (P>0) if (!data->a[Mindex.at(P-1)].bd[bnd2par-1]) return 0;
    if (!data->a[id].bd[bnd2par-1]) return 0;
    if (!data->a[id].bd[bnd2des-1]) return 0;
    for (int k1=0;k1<Cindex.size();k1++) {
        if (Pindex.at(k1)==Cindex.at(n)) {
            if (!data->a[Mindex.at(k1)].bd[bnd2des-1]) return 0;
        }
    }
    */


    //reset();
    //vector<int> tmp;

    if (1) {  //desendent atm after new atom added, i.e. C(n)
        bool totbndchk=0;
        bool bndchk=0;

        int bndsum=0; // max available # of bonds of C(n) (calc. from pool)
        bndsum+=bnd2des;
        vector<int> ordcount(6,0);
        ordcount.at(bnd2des-1)+=1;
        if (Cyindex.at(n).size()) {
            for (int k1=0;k1<Cyindex.at(n).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(n).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(n).at(k1)-1)-1)++;
            }
        }
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Pindex.at(k1)==Cindex.at(n)) {
                bndsum+=Rindex.at(k1);
                ordcount.at(Rindex.at(k1)-1)+=1;
            }
        }

        int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[M].order.at(k2);
        }
        /*
        if (bndsum<=maxbnd) totbndchk=1;
            return 0;
        }
        */

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of C(n) after addition of atm(id)" match the definition of atom C(n) in pool
            if (ordcount.at(k2)>data->a[M].bd[k2]) { 
                bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[M].bd[k2]) { 
                bndchk=1;
            }
        }

    }
    // Ensure the bonds of atom P(n) (before insertion) are OK after atm(id) inserted and use a bnd2par bond to connect P(n)
    if (P>0) {
        bool totbndchk=0;
        bool bndchk=0;

        int bndsum=0;  // # of bonds of atom P(n) if atom M(n) is changed to atom(id) and a $bnd2par bond of atom(id) is used to connected with P(n).
        if (Rindex.at(P-1)>0) bndsum+=Rindex.at(P-1);
        vector<int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
        if (Rindex.at(P-1)>0) ordcount.at(Rindex.at(P-1)-1)+=1;  //20200712 // bnd P(n) to P(P(n))
        if (Cyindex.at(P-1).size()) {
            for (int k1=0;k1<Cyindex.at(P-1).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(P-1).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(P-1).at(k1)-1)-1)++;
            }
        }
        for (int k=0;k<Cindex.size();k++) {
            if (Pindex.at(k)==Cindex.at(P-1) && k!=n) {
                bndsum+=Rindex.at(k);
                ordcount.at(Rindex.at(k)-1)+=1;
            }
        }
        bndsum+=bnd2par;
        ordcount.at(bnd2par-1)+=1;

        int maxbnd=0; // max available # of bonds for atom P(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(P-1)].order.at(k2);
        }
        /*
        if (bndsum<=maxbnd) totbndchk=1;
        else {
            totbndchk=0;
            return 0;
        }
        */

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of P(n) after connecting with atom(id)" match the definition of atom P(n) in pool
            if (ordcount.at(k2)>data->a[Mindex.at(P-1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(P-1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

        //if (!totbndchk || !bndchk) {
            //cout << "Failed exchange: (totbndchk, bndchk) = (" << totbndchk << ", " << bndchk << ")" << endl;
            //return 0;
        //}
    }
    // Ensure the bonds of atm(id) are OK if atm(id) is inserted and use bnd2par & bnd2des bonds to connect P(n) & C(n), respectively.
    if (1) {
        bool totbndchk=0;
        bool bndchk=0;

        int id_bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        if (P>0) id_bndsum+=bnd2par;
        vector<int> id_ordcount(6,0);
        if (P>0) id_ordcount.at(bnd2par-1)+=1;

		id_bndsum+=bnd2des;
		id_ordcount.at(bnd2des-1)+=1;

        int id_maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            id_maxbnd+=data->a[id].order.at(k2);
        }
        /*
        if (id_bndsum<=id_maxbnd) totbndchk=1;
        else {
            totbndchk=0;
            return 0;
        }
        */

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (id_ordcount.at(k2)>data->a[id].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                bndchk=0;
                return 0;
            }
            else if (k2>=2 && id_ordcount.at(k2)<=data->a[id].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

    }


	for (int k3=0;k3<Cindex.size();k3++) {
		Cindex.at(k3)+=(Cindex.size()+2);
		Pindex.at(k3)+=(Cindex.size()+2);
	}

	vector<int> curatm(0);

	bool goout=0;
	Cindex.at(n)+=(Cindex.size()+2);
    Pindex.at(n)+=(Cindex.size()+2);
    Rindex.at(n)=bnd2des;

	curatm.push_back(C+Cindex.size()+2);

	do {
		vector<int> tmpp(0);
        for (int k3=0;k3<curatm.size();k3++) {
            for (int k2=0;k2<Cindex.size();k2++) {
                if (Pindex.at(k2)==curatm.at(k3)) {
					Pindex.at(k2)+=(Cindex.size()+2);
                    tmpp.push_back(Cindex.at(k2));
					Cindex.at(k2)+=(Cindex.size()+2);
                }
            }
		}
        if (tmpp.size()>0) {
            curatm.resize(0);
            for (int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
        }
        else {
            goout=1;
            break;
        }

	} while (!goout);

	if (1) {
		//int newC=C+2*Cindex.size()+4;
		//int newP=P+2*Cindex.size()+4;
		int newC=C+Cindex.size()+2;
		int newP=P+Cindex.size()+2;
        Pindex.insert(Pindex.begin()+C-1,newP);
        Cindex.insert(Cindex.begin()+C-1,newC);
        Mindex.insert(Mindex.begin()+C-1,id);
        Rindex.insert(Rindex.begin()+C-1,bnd2par);
        //Cyindex.insert(Cyindex.begin()+C-1,0);
		Cyindex.insert(Cyindex.begin()+C-1,vector<int>(0));
        //ctsisomer.insert("");
        ctsisomer.at(0).insert(ctsisomer.at(0).begin()+C-1,"");
        ctsisomer.at(1).insert(ctsisomer.at(1).begin()+C-1,"");
        if (para.protect) protect.insert(protect.begin()+C-1,0);		
	}

	Pindex.at(n+1)=Cindex.at(n);
	Pindex.at(0)=0;

    for (int i=0;i<Cindex.size();i++) {
        int j=Cindex.at(i);
        Cindex.at(i)=(i+1);
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Pindex.at(k1)==j) Pindex.at(k1)=Cindex.at(i);
        }
    }


	if (1) {
		if (R==2) {
			ctsisomer.at(0).at(C)=="";
			ctsisomer.at(1).at(C)=="";
			ctsisomer.at(0).at(P-1)=="";
        	ctsisomer.at(1).at(P-1)=="";
		}

        bool fl1=0,fl2=0,fl3=0;
		if (0) {
        	if (Pindex.at(C-1)>0) if (Mindex.at(Pindex.at(C-1)-1)==2 || Mindex.at(Pindex.at(C-1)-1)==7 || Mindex.at(Pindex.at(C-1)-1)==16) fl1=1;
			if (Mindex.at(C)==2 || Mindex.at(C)==7 || Mindex.at(C)==16) fl3=1;
        	if (id==2 || id==7 || id==16) fl2=1;
		}
		if (1) {
			if (Pindex.at(C-1)>0) fl1=1;
			fl2=1;
			fl3=1;
		}

		if (bnd2par==2 && bnd2des==1 && fl1 && fl2) {
			if (Rindex.at(Pindex.at(C-1)-1)!=2 && ctsisomer.at(0).at(Pindex.at(C-1)-1)=="" && ctsisomer.at(1).at(Pindex.at(C-1)-1)=="") {
				ctsisomer.at(0).at(Pindex.at(C-1)-1)="/";
				ctsisomer.at(1).at(Pindex.at(C-1)-1)="";
			}
			if (ctsisomer.at(0).at(C-1)=="" && ctsisomer.at(1).at(C-1)=="") {
				ctsisomer.at(0).at(C-1)="";
				if (!cistrans) ctsisomer.at(1).at(C-1)="/";
				else ctsisomer.at(1).at(C-1)="\\";
			}
		}
		else if (bnd2par==1 && bnd2des==2 && fl2 && fl3) {
			if (ctsisomer.at(0).at(C-1)=="" && ctsisomer.at(1).at(C-1)=="") {
				ctsisomer.at(0).at(C-1)="/";
				ctsisomer.at(1).at(C-1)="";
			}
			if (ctsisomer.at(0).at(C)=="" && ctsisomer.at(1).at(C)=="") {
				ctsisomer.at(0).at(C)="";
				if (!cistrans) ctsisomer.at(1).at(C)="/";
				else ctsisomer.at(1).at(C)="\\";
			}
		}
	}

    //reset();

    del_unpaired_ring_no(); //20200627
    decyc_small_ring(5);

    //reset();
    //chk_cistrans(); //20200806
    mds2smi();


	return 1;
}

int MOLECULE::subtraction(int n,int mode,bool cistrans) {
	// mode 1: subtract an atom
	// mode 2: subtract all atoms between Cindex=n and the end of the corresponding branch
	int i,j,k;
	//vector<int> ref,pr;

	//long double time1=time(NULL);

	if (Cindex.size()<2) return 0;
	if (n>=Cindex.size()) return 0;
	if (n<=0 && 0) {
		//ofstream log("molecule.log",ios::app);
		//log<<"Warning : This molecule would be empty because of point, "<<n<<endl;
		//log.close();
		return 0;
	}
	/*
    if (para.protect) pr.push_back(protect.at(n));
    if (para.protect) {
        for (i=0;i<pr.size();i++) {
            if (pr.at(i)) return 0;
        }
    }
	*/

	//ref.push_back(Cindex.at(n));
	//ref.push_back(n);  // wrong? 20190721

    if (1) if (para.protect && protect.at(n)) return 0;   // terminate if it would subtract protected atom
    if (1) if (data->a[Mindex.at(n)].chg) return 0;  // terminate if it would subtract charged atom

	bool usepibnd=0;
	bool isnotend=0;

	vector<int> Merge_C(0);
	if (mode==1) {  //20200131
		bool totbndchk=0,bndchk=0;
		int P=Pindex.at(n);
		if (P>0) {
			Merge_C.push_back(Cindex.at(P-1));
			totbndchk=0;
			bndchk=0;

			int bndsum=0; // # of bonds of atom P if atom j is deleted and P is connected with each of the descendant atom k
			if (Rindex.at(P-1)>0) bndsum+=Rindex.at(P-1);  
			vector<int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
			if (Rindex.at(P-1)>0) ordcount.at(Rindex.at(P-1)-1)+=1;  //20200712
	        if (Cyindex.at(P-1).size()) {
            	for (int k1=0;k1<Cyindex.at(P-1).size();k1++) {
                	bndsum+=Cybnd.at(Cyindex.at(P-1).at(k1)-1);
                	ordcount.at(Cybnd.at(Cyindex.at(P-1).at(k1)-1)-1)++;
            	}
        	}
            //if (Cyindex.at(n).size() && 1) {
            //    for (int k1=0;k1<Cyindex.at(n).size();k1++) {
            //        bndsum-=Cybnd.at(Cyindex.at(n).at(k1)-1);
            //        ordcount.at(Cybnd.at(Cyindex.at(n).at(k1)-1)-1)--;
            //    }
            //}
			for (k=0;k<Cindex.size();k++) { 
				if (Pindex.at(k)==Cindex.at(n)) { // For each descendant atom k of atom j // orig Pindex.at(k)==Cindex.at(j)
					isnotend=1;
					Merge_C.push_back(Cindex.at(k));
					bndsum+=Rindex.at(k);
					ordcount.at(Rindex.at(k)-1)++;
					if (Rindex.at(k)==2) usepibnd=1;
				}
				if (Pindex.at(k)==P && k!=n) {
                	bndsum+=Rindex.at(k);
                	ordcount.at(Rindex.at(k)-1)++;
				}
			}
			

			int maxbnd=0; // max available # of bonds for parental atom k1 (calc. from pool)
			for (int k2=0;k2<6;k2++) { 
				maxbnd+=data->a[Mindex.at(P-1)].order.at(k2);
			}
			/*
			if (bndsum<=maxbnd) totbndchk=1;
			else totbndchk=0;
			*/

			for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
				if (ordcount.at(k2)>data->a[Mindex.at(P-1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(P-1)].order.at(k2)
					bndchk=0;
					return 0;
				}
				else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(P-1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(P-1)].order.at(k2)
					bndchk=1;
				}
			}

			//if (bndchk) break; // orig totbndchk && bndchk
			//else {
				//cout << "Failed subtraction: (totbndchk, bndchk) = (" << totbndchk << ", " << bndchk << ")" << endl;
				//cout << "Failed subtraction: bndchk = " << bndchk << endl;
				//return 0;
			//}
		}
		else {
			Merge_C.push_back(0);
            for (k=0;k<Cindex.size();k++) {
                if (Pindex.at(k)==Cindex.at(n)) { // For each descendant atom k of atom j // orig Pindex.at(k)==Cindex.at(j)
                    Merge_C.push_back(Cindex.at(k));
                }
            }
			if (Merge_C.size()!=2) return 0;
		}
	}

	if (Cyindex.at(n).size()) {
		if (Cyindex.at(n).size()>1) sort(Cyindex.at(n).begin(),Cyindex.at(n).end());

		for (int g=Cyindex.at(n).size()-1;g>=0;g--) {
    		for (k=0;k<Cindex.size();k++) {  // 20200131 prevent unpaired ring number
				if (Cyindex.at(k).size() && k!=n) { 
					for (int k1=0;k1<Cyindex.at(k).size();k1++) {
						if (Cyindex.at(k).at(k1)==Cyindex.at(n).at(g)) { 
							Cyindex.at(k).erase(Cyindex.at(k).begin()+k1);
							k--;
							k1--;
						}
					}
                }
			}	
		}

        for (int g=Cyindex.at(n).size()-1;g>=0;g--) {
            for (k=0;k<Cindex.size();k++) {  // 20200131 prevent unpaired ring number
                if (Cyindex.at(k).size() && k!=n) { // orig Cyindex.at(j) && Cyindex.at(k) && k!=j
                    for (int k1=0;k1<Cyindex.at(k).size();k1++) {
                        if (Cyindex.at(k).at(k1)>Cyindex.at(n).at(g)) { //w1==w
							Cyindex.at(k).at(k1)--;
                        }
                    }
                }
            }
			if (Cybnd.size()>=Cyindex.at(n).at(g)) Cybnd.erase(Cybnd.begin()+Cyindex.at(n).at(g)-1);
			if_circle--;
        }

	}

    if (1) {
		if (Pindex.at(n)<=0 && isnotend) {
            for (int k1=0;k1<Cindex.size();k1++) {
                if (Pindex.at(k1)==Cindex.at(n) && Rindex.at(k1)==2) {
                    //ctsisomer.at(0).at(k1)="";
                    ctsisomer.at(1).at(k1)="";
                }
            }
		}
		if (Pindex.at(n)>0 && !isnotend) {
			ctsisomer.at(0).at(Pindex.at(n)-1)="";
			//ctsisomer.at(1).at(Pindex.at(n)-1)="";
		}
		else {
        	if (Rindex.at(n)==2) {
	            for (int k1=0;k1<Cindex.size();k1++) {
    	            if (Pindex.at(k1)==Cindex.at(n) && Rindex.at(k1)!=2) {
        	            //ctsisomer.at(0).at(k1)="";
            	        ctsisomer.at(1).at(k1)="";
                	}
            	}
        	}
        	if (Rindex.at(n)!=2 && Pindex.at(n)>1) {
				bool go=0;
                for (int k1=0;k1<Cindex.size();k1++) {
                    if (Pindex.at(k1)==Cindex.at(n) && Rindex.at(k1)==2) {
						go=1;
                        //ctsisomer.at(0).at(k1)="";
                        if (!cistrans) ctsisomer.at(1).at(k1)="/";
						else ctsisomer.at(1).at(k1)="\\";
                    }
                }
				if (go) ctsisomer.at(0).at(Pindex.at(n)-1)="/";
        	}
		}
    }


    Mindex.erase(Mindex.begin()+n);
	Cindex.erase(Cindex.begin()+n);
    Pindex.erase(Pindex.begin()+n);
	Rindex.erase(Rindex.begin()+n);
	Cyindex.erase(Cyindex.begin()+n);	
	//ctsisomer.erase(ctsisomer.begin()+n);
	ctsisomer.at(0).erase(ctsisomer.at(0).begin()+n);
	ctsisomer.at(1).erase(ctsisomer.at(1).begin()+n);
	if (para.protect) protect.erase(protect.begin()+n);

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
        	//if (Cindex.at(i)!=(i+1)) {
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
        	//}
    	}


		if (Merge_C.size()>=2) {
			for (k=1;k<Merge_C.size();k++) {
				Pindex.at(Merge_C.at(k)-1)=Merge_C.at(0);
				if (!Merge_C.at(0)) Rindex.at(Merge_C.at(k)-1)=-1;
			}
		}
		
	}

	//vector<int>().swap(ref);
	//vector<int>().swap(pr);
	vector<int>().swap(Merge_C);

	/*
    vector<int> CDofct(0);
    CDofct.push_back(Cindex.at(n));
    for (int i=n+1;i<Cindex.size();i++) {
        if (Rindex.at(i)==2 && Pindex.at(i)==Cindex.at(n)) CDofct.push_back(Cindex.at(i));
    }
    for (int i=0;i<CDofct.size();i+=2) {
        if (CDofct.at(i)>0) {
            if (ctsisomer.at(0).at(CDofct.at(i)-1)!="" && ctsisomer.at(1).at(CDofct.at(i+1)-1)=="") ctsisomer.at(1).at(CDofct.at(i+1)-1)="/";
            if (ctsisomer.at(0).at(CDofct.at(i)-1)=="" && ctsisomer.at(1).at(CDofct.at(i+1)-1)!="") ctsisomer.at(0).at(CDofct.at(i)-1)="/";
            if (ctsisomer.at(0).at(CDofct.at(i)-1)=="" && ctsisomer.at(1).at(CDofct.at(i+1)-1)=="") {
                ctsisomer.at(1).at(CDofct.at(i+1)-1)="/";
                ctsisomer.at(0).at(CDofct.at(i)-1)="/";
            }

        }
    }

    if (Rindex.at(n)==2 && Pindex.at(n)>0) {
        if (ctsisomer.at(0).at(Pindex.at(n)-1)!="" && ctsisomer.at(1).at(Cindex.at(n)-1)=="") ctsisomer.at(1).at(Cindex.at(n)-1)="/";
        if (ctsisomer.at(0).at(Pindex.at(n)-1)=="" && ctsisomer.at(1).at(Cindex.at(n)-1)!="") ctsisomer.at(0).at(Pindex.at(n)-1)="/";
        if (ctsisomer.at(0).at(Pindex.at(n)-1)=="" && ctsisomer.at(1).at(Cindex.at(n)-1)=="") {
            ctsisomer.at(1).at(Cindex.at(n)-1)="/";
            ctsisomer.at(0).at(Pindex.at(n)-1)="/";
        }
    }
	*/

	del_unpaired_ring_no();
	decyc_small_ring(5);

	//reset();
	//chk_cistrans(); //20200806
	mds2smi();

    //time1=time(NULL)-time1;
    //cout << "SUBTRACTION FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}

int MOLECULE::change_bnd(int n,int id1,int id2,int bnd,bool cistrans) {

    if (n<=0 || n>=Cindex.size()) return 0;
    int i,j,k,M,P,C,R;
    M=Mindex.at(n);
    P=Pindex.at(n);
    C=Cindex.at(n);
    R=Rindex.at(n);

    if (para.protect) {
        if (protect.at(n)) return 0;
        if (P>0) if (protect.at(P-1)) return 0;
    }
    if (bnd>3 || bnd<0) return 0;
    if (data->a[id2].chg!=data->a[Mindex.at(n)].chg) return 0;
    else if (data->a[id2].chg==data->a[Mindex.at(n)].chg && data->a[id2].chg!=0) {
        if (data->a[id2].index<=2 && data->a[Mindex.at(n)].index>=7) return 0;
        if (data->a[id2].index>=7 && data->a[Mindex.at(n)].index<=7) return 0;
    }
	if (P>0) {
    	if (data->a[id1].chg!=data->a[Mindex.at(P-1)].chg) return 0;
    	else if (data->a[id1].chg==data->a[Mindex.at(P-1)].chg && data->a[id1].chg!=0) {
        	if (data->a[id1].index<=2 && data->a[Mindex.at(P-1)].index>=7) return 0;
        	if (data->a[id1].index>=7 && data->a[Mindex.at(P-1)].index<=7) return 0;
    	}
	}

    if (!data->a[id1].bd[bnd-1]) return 0;
    if (!data->a[id2].bd[bnd-1]) return 0;


    //reset();
    //vector<int> tmp;


    // Ensure the bonds of atom C(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
    if (1) {
        bool totbndchk=0;
        bool bndchk=0;

        int id_bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        id_bndsum+=bnd;
        vector<int> id_ordcount(6,0);
        id_ordcount.at(bnd-1)+=1;
        if (Cyindex.at(n).size()) {
            for (int k1=0;k1<Cyindex.at(n).size();k1++) {
                id_bndsum+=Cybnd.at(Cyindex.at(n).at(k1)-1);
                id_ordcount.at(Cybnd.at(Cyindex.at(n).at(k1)-1)-1)++;
            }
        }
        for (int k1=0;k1<Cindex.size();k1++) {
            if (Pindex.at(k1)==Cindex.at(n)) {
                id_bndsum+=Rindex.at(k1);
                id_ordcount.at(Rindex.at(k1)-1)+=1;
            }
        }

        int id_maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            id_maxbnd+=data->a[id2].order.at(k2);
        }
        /*
        if (id_bndsum<=id_maxbnd) totbndchk=1;
            return 0;
        }
        */

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (id_ordcount.at(k2)>data->a[id2].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                bndchk=0;
                return 0;
            }
            else if (k2>=2 && id_ordcount.at(k2)<=data->a[id2].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

    }
    // Ensure the bonds of atom P(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
    if (P>0) {
        bool totbndchk=0;
        bool bndchk=0;

        int bndsum=0;  // # of bonds of atom P(n) if atom M(n) is changed to atom(id) and a $bnd2par bond of atom(id) is used to connected with P(n).
        if (Rindex.at(P-1)>0) bndsum+=Rindex.at(P-1);
        vector<int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
        if (Rindex.at(P-1)>0) ordcount.at(Rindex.at(P-1)-1)+=1;  //20200712 // bnd P(n) to P(P(n))
        if (Cyindex.at(P-1).size()) {
            for (int k1=0;k1<Cyindex.at(P-1).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(P-1).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(P-1).at(k1)-1)-1)++;
            }
        }
        for (k=0;k<Cindex.size();k++) {
            if (Pindex.at(k)==Cindex.at(P-1) && k!=n) {
                bndsum+=Rindex.at(k);
                ordcount.at(Rindex.at(k)-1)+=1;
            }
        }
        bndsum+=bnd;
        ordcount.at(bnd-1)+=1;

        int maxbnd=0; // max available # of bonds for atom P(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[id1].order.at(k2);
        }
        /*
        if (bndsum<=maxbnd) totbndchk=1;
        else {
            totbndchk=0;
            return 0;
        }
        */

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of P(n) after connecting with atom(id)" match the definition of atom P(n) in pool
            if (ordcount.at(k2)>data->a[id1].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[id1].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

        //if (!totbndchk || !bndchk) {
            //cout << "Failed exchange: (totbndchk, bndchk) = (" << totbndchk << ", " << bndchk << ")" << endl;
            //return 0;
        //}
    }


    Mindex.at(n)=id2;
	if (P>0) Mindex.at(P-1)=id1;
    Rindex.at(n)=bnd;

	if (1) {
		if (R==2 && bnd!=2) {
			if (P>0) {
				ctsisomer.at(0).at(n)="";
				ctsisomer.at(1).at(n)="";
				ctsisomer.at(0).at(P-1)="";
				ctsisomer.at(1).at(P-1)="";
			}
		}
		else if (R!=2 && bnd==2) {
			if (P>0) {
				if (ctsisomer.at(0).at(n)=="" && ctsisomer.at(1).at(n)=="") {
					ctsisomer.at(0).at(n)="";
					if (!cistrans) ctsisomer.at(1).at(n)="/";
					else ctsisomer.at(1).at(n)="\\";
				}
				if (ctsisomer.at(0).at(P-1)=="" && ctsisomer.at(1).at(P-1)=="") {
					ctsisomer.at(0).at(P-1)="/";
					ctsisomer.at(1).at(P-1)="";
				}
			}
		}
	}

    //reset();

    del_unpaired_ring_no(); //20200627
    decyc_small_ring(5);

    //reset();
    //chk_cistrans(); //20200806
    mds2smi();

    return 1;

}


int MOLECULE::change_ele(int n,int id,int bnd2par,int bnd2des,bool cistrans) {
	
    if (n<0 || n>=Cindex.size()) return 0;
    int i,j,k,M,P,C,R;
    M=Mindex.at(n);
    P=Pindex.at(n);
    C=Cindex.at(n);
    R=Rindex.at(n);

    if (para.protect) {
        if (protect.at(n)) return 0;
        if (P>0) if (protect.at(P-1)) return 0;
    }
    if (bnd2par>3 || bnd2par<=0) return 0;
    if (bnd2des>3 || bnd2des<=0) return 0;
    if (data->a[id].chg!=data->a[Mindex.at(n)].chg) return 0;
	else if (data->a[id].chg==data->a[Mindex.at(n)].chg && data->a[id].chg!=0) {
		if (data->a[id].index<=2 && data->a[Mindex.at(n)].index>=7) return 0;
		if (data->a[id].index>=7 && data->a[Mindex.at(n)].index<=7) return 0;
	}

	/*
	if (P>0) if (!data->a[Mindex.at(P-1)].bd[bnd2par-1]) return 0;
    if (!data->a[id].bd[bnd2par-1]) return 0;
    if (!data->a[id].bd[bnd2des-1]) return 0;
    for (int k1=0;k1<Cindex.size();k1++) {
        if (Pindex.at(k1)==Cindex.at(n)) {
			if (!data->a[Mindex.at(k1)].bd[bnd2des-1]) return 0;
        }
    }
	*/

    
    //reset();
    //vector<int> tmp;

	bool isnotend=0;
	bool usepibnd=0;

	// Ensure the bonds of atom C(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
	if (1) {
        bool totbndchk=0;
        bool bndchk=0;

	    int id_bndsum=0; // max available # of bonds of element(id) (calc. from pool)
		if (P>0) id_bndsum+=bnd2par;
		vector<int> id_ordcount(6,0);
		if (P>0) id_ordcount.at(bnd2par-1)+=1;
		if (Cyindex.at(n).size()) {
            for (int k1=0;k1<Cyindex.at(n).size();k1++) {
                id_bndsum+=Cybnd.at(Cyindex.at(n).at(k1)-1);
                id_ordcount.at(Cybnd.at(Cyindex.at(n).at(k1)-1)-1)++;
            }
		}
    	for (int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(n)) {
				isnotend=1;
				id_bndsum+=bnd2des;
				id_ordcount.at(bnd2des-1)+=1;
				if (Rindex.at(k1)==2) usepibnd=1;
			}
		}

        int id_maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            id_maxbnd+=data->a[id].order.at(k2); 
        }
		/*
        if (id_bndsum<=id_maxbnd) totbndchk=1;
        else {
			totbndchk=0;
			return 0;
		}
		*/

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (id_ordcount.at(k2)>data->a[id].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                bndchk=0;
				return 0;
            }
            else if (k2>=2 && id_ordcount.at(k2)<=data->a[id].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                bndchk=1;
            }
        }

	}


	// Ensure the bonds of atom P(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
	if (P>0) {
    	bool totbndchk=0;
        bool bndchk=0;

        int bndsum=0;  // # of bonds of atom P(n) if atom M(n) is changed to atom(id) and a $bnd2par bond of atom(id) is used to connected with P(n).
		if (Rindex.at(P-1)>0) bndsum+=Rindex.at(P-1); 
        vector<int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
        if (Rindex.at(P-1)>0) ordcount.at(Rindex.at(P-1)-1)+=1;  //20200712 // bnd P(n) to P(P(n))
        if (Cyindex.at(P-1).size()) {
            for (int k1=0;k1<Cyindex.at(P-1).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(P-1).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(P-1).at(k1)-1)-1)++;
            }
        }
        for (k=0;k<Cindex.size();k++) {
            if (Pindex.at(k)==Cindex.at(P-1) && k!=n) {
                bndsum+=Rindex.at(k);
                ordcount.at(Rindex.at(k)-1)+=1;
            }
        }
		bndsum+=bnd2par;
		ordcount.at(bnd2par-1)+=1;

        int maxbnd=0; // max available # of bonds for atom P(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(P-1)].order.at(k2);
        }
		/*
        if (bndsum<=maxbnd) totbndchk=1;
        else {
			totbndchk=0;
			return 0;
		}
		*/

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of P(n) after connecting with atom(id)" match the definition of atom P(n) in pool
        	if (ordcount.at(k2)>data->a[Mindex.at(P-1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
        		bndchk=0;
				return 0;
        	}
        	else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(P-1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
            	bndchk=1;
        	}
    	}

        //if (!totbndchk || !bndchk) {
            //cout << "Failed exchange: (totbndchk, bndchk) = (" << totbndchk << ", " << bndchk << ")" << endl;
            //return 0;
        //}
	}

    // Ensure the bonds of atom des(n) are OK if M(n) is changed to atom(id) and use a $bnd2des bond to connect des(n)
	if (isnotend) {
    	for (int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(n)) { // For the descendant atoms C(k1) of atom C(n)
            	bool totbndchk=0;
            	bool bndchk=0;

            	int bndsum=0; // # of bonds of atom C(k1) if atom C(n) is changed to atom(id) and is connected with C(k1) with a $bnd2des bond.
				bndsum+=bnd2des;  
            	vector<int> ordcount(6,0); // count # of bonds atom C(k1) will use after connecting with atom C(n) (in each bond order respectively)
        		ordcount.at(bnd2des-1)+=1;  //20200712
	        	if (Cyindex.at(k1).size()) {
		            for (int k2=0;k2<Cyindex.at(k1).size();k2++) {
                		bndsum+=Cybnd.at(Cyindex.at(k1).at(k2)-1);
                		ordcount.at(Cybnd.at(Cyindex.at(k1).at(k2)-1)-1)++;
            		}
        		}
            	for (k=0;k<Cindex.size();k++) { // For each des(k1) atom
                    if (Pindex.at(k)==Cindex.at(k1)) {
                        bndsum+=Rindex.at(k); // calc. the total # of bonds des(k1) atoms used to connect with C(k1).
                        ordcount.at(Rindex.at(k)-1)+=1;
                    }
                }

                int maxbnd=0; // max available # of bonds for C(k1) (calc. from pool)
                for (int k2=0;k2<6;k2++) {
                    maxbnd+=data->a[Mindex.at(k1)].order.at(k2);
                }
				/*
                if (bndsum<=maxbnd) totbndchk=1;
                else {
					totbndchk=0;
					return 0;
				}
				*/

                for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with des(k1) atoms" match the definition of atom k1 in pool
                    if (ordcount.at(k2)>data->a[Mindex.at(k1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                        bndchk=0;
						return 0;
                    }
                    else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(k1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                        bndchk=1;
                    }
                }

            }
        }

    }


	Mindex.at(n)=id;
	if (P>0) Rindex.at(n)=bnd2par;
	else Rindex.at(n)=-1;
	if (isnotend) {
    	for (int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(n)) {
				Rindex.at(k1)=bnd2des;
			}
		}
	}

	if (1) {
		if (R!=2 && bnd2par==2 && bnd2des!=2) {
			if (P>0) {
				if (ctsisomer.at(0).at(n)=="" && ctsisomer.at(1).at(n)=="") {
					ctsisomer.at(0).at(n)="";
					if (!cistrans) ctsisomer.at(1).at(n)="/";
					else ctsisomer.at(1).at(n)="\\";
				}
				if (ctsisomer.at(0).at(P-1)=="" && ctsisomer.at(1).at(P-1)=="") {
					ctsisomer.at(0).at(P-1)="/";
					ctsisomer.at(1).at(P-1)="";
				}
			}
		}
		else if (!usepibnd && bnd2par!=2 && bnd2des==2) {
			if (ctsisomer.at(0).at(n)=="" && ctsisomer.at(1).at(n)=="") {
				ctsisomer.at(0).at(n)="/";
				ctsisomer.at(1).at(n)="";
			}
			if (isnotend) {
				for (int k1=0;k1<Cindex.size();k1++) {
					if (Pindex.at(k1)==Cindex.at(n) && Rindex.at(k1)==2) {
						if (ctsisomer.at(0).at(k1)=="" && ctsisomer.at(1).at(k1)=="") {
							ctsisomer.at(0).at(k1)="";
							if (!cistrans) ctsisomer.at(1).at(k1)="/";
							else ctsisomer.at(1).at(k1)="\\";
						}
					}
				}
			}
		}
	}

	//reset();

    del_unpaired_ring_no(); //20200627
    decyc_small_ring(5);

    //reset();
	//chk_cistrans(); //20200806
	mds2smi();

	return 1;

}

int MOLECULE::decyclization(int ringnum) {

	if (!ringnum) return 0;
	if (ringnum>if_circle) return 0;

    for (int i=0;i<Cyindex.size();i++) {
        if (Cyindex.at(i).size()) {
            for (int k3=0;k3<Cyindex.at(i).size();k3++) {
                if (Cyindex.at(i).at(k3)==ringnum) {
					Cyindex.at(i).erase(Cyindex.at(i).begin()+k3);
					k3--;
                }
            }
        }
    }

    for (int i=0;i<Cyindex.size();i++) {
        if (Cyindex.at(i).size()) {
            for (int k3=0;k3<Cyindex.at(i).size();k3++) {
                if (Cyindex.at(i).at(k3)>ringnum) {
                    Cyindex.at(i).at(k3)--;
                }
            }
        }

    }


    if (Cybnd.size()>=ringnum) {
		Cybnd.erase(Cybnd.begin()+ringnum-1);
    }

	if_circle--;

    return 1;

}


int MOLECULE::rechg() {
	int i=0,j=0;
	chg=0;
	for (j=0;j<Cindex.size();j++) {
		chg+=data->a[Mindex.at(j)].chg;
	}
	return 1;
}


int MOLECULE::print() {
	int i;
	//if (para.protect) prct();

	cout << setfill(' ');
	cout<<"SMILES : "<<left<<molesmi<<endl;
	cout<<"natom : "<<setw(4)<<left<<Cindex.size()<<endl;
	cout<<"Pindex : ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Pindex.at(i)<<" ";
	cout<<endl;
	cout<<"Cindex : ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Cindex.at(i)<<" ";
	cout<<endl;
	cout<<"Rindex : ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Rindex.at(i)<<" ";
	cout<<endl;
	cout<<"Mindex : ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Mindex.at(i)<<" ";
	cout<<endl;
	cout<<"Cyindex : ";
	for (i=0;i<Cindex.size();i++) {
        if (!Cyindex.at(i).size()) cout<<setw(4)<<left<<setfill(' ')<<"0"<<" ";
        else {
			stringstream buf("");
            for (int j=0;j<Cyindex.at(i).size();j++) {
                if (j<Cyindex.at(i).size()-1) buf<<Cyindex.at(i).at(j)<<",";
                else buf<<Cyindex.at(i).at(j);

                //if (j<Cyindex.at(i).size()-1) cout<<left<<Cyindex.at(i).at(j)<<",";
                //else cout<<left<<Cyindex.at(i).at(j)<<" ";
            }
			cout<<setw(4)<<left<<setfill(' ')<<buf.str()<<" ";
        }
	}
	cout<<endl;
    cout<<"ctsisomer_start : ";
    for (i=0;i<ctsisomer.at(0).size();i++) {
        if (ctsisomer.at(0).at(i)!="") cout<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(0).at(i)<<" ";
        else cout<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    cout<<endl;
    cout<<"ctsisomer_end : ";
    for (i=0;i<ctsisomer.at(1).size();i++) {
        if (ctsisomer.at(1).at(i)!="") cout<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(1).at(i)<<" ";
        else cout<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    cout<<endl;


    cout<<"if_circle : "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;

    cout<<"Cybnd : ";
    for (i=0;i<Cybnd.size();i++) cout<<setw(4)<<left<<Cybnd.at(i)<<" ";
	cout<<endl;

	return 1;
}		

int MOLECULE::print(ofstream &ouf) {
	int i;
	//if (para.protect) prct();

	ouf << setfill(' ');
	ouf<<"SMILES : "<<left<<molesmi<<endl;
	ouf<<"natom : "<<setw(4)<<left<<Cindex.size()<<endl;
	ouf<<"Pindex : ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Pindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"Cindex : ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Cindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"Rindex : ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Rindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"Mindex : ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Mindex.at(i)<<" ";
	ouf<<endl;
	ouf<<"Cyindex : ";
	for (i=0;i<Cindex.size();i++) {
        if (!Cyindex.at(i).size()) ouf<<setw(4)<<left<<setfill(' ')<<"0"<<" ";
        else {
			stringstream buf("");
            for (int j=0;j<Cyindex.at(i).size();j++) {
                if (j<Cyindex.at(i).size()-1) buf<<Cyindex.at(i).at(j)<<",";
                else buf<<Cyindex.at(i).at(j);

                //if (j<Cyindex.at(i).size()-1) ouf<<left<<Cyindex.at(i).at(j)<<",";
                //else ouf<<left<<Cyindex.at(i).at(j)<<" ";
            }
			ouf<<setw(4)<<left<<setfill(' ')<<buf.str()<<" ";
        }
	}
	ouf<<endl;
	if (para.protect) {
		ouf<<"protect : ";
		for (i=0;i<protect.size();i++) ouf<<setw(4)<<left<<protect.at(i)<<" ";
		ouf<<endl;
	}
    ouf<<"ctsisomer_start : ";
    for (i=0;i<ctsisomer.at(0).size();i++) {
        if (ctsisomer.at(0).at(i)!="") ouf<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(0).at(i)<<" ";
        else ouf<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    ouf<<endl;
    ouf<<"ctsisomer_end : ";
    for (i=0;i<ctsisomer.at(1).size();i++) {
        if (ctsisomer.at(1).at(i)!="") ouf<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(1).at(i)<<" ";
        else ouf<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    ouf<<endl;


    ouf<<"if_circle : "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;

    ouf<<"Cybnd : ";
    for (i=0;i<Cybnd.size();i++) ouf<<setw(4)<<left<<Cybnd.at(i)<<" ";
	ouf<<endl;


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
				if (data->a[Mindex.at(i)].chg) {
					protect.at(i)=1;
					//for (j=a;j<b;j++) {
					//	//if ( abs(Cindex.at(i)-Cindex.at(j))<=1 ) protect.at(j)=1;
					//	if (Cindex.at(j)==Pindex.at(i) || Cindex.at(i)==Pindex.at(j)) protect.at(j)=1;
					//}
				}
			}
		}

		if (0 && Cindex.size()>0) { // ring
			if (0) {
				protect.resize(0);
				protect.resize(Cindex.size(),0); //20190801
			}

			if (1) {
				vector<int> C_ringmember(0);
				for (int k1=if_circle;k1>=1;k1--) {
					C_ringmember.resize(0);
					int rep_quota=1;
					for (int k2=Cyindex.size()-1;k2>=0;k2--) {
						if (Cyindex.at(k2).size()) {
							for (int k3=0;k3<Cyindex.at(k2).size();k3++) {
								if (Cyindex.at(k2).at(k3)==k1) {
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
								else if (C_ringmember.at(k2)==C_ringmember.at(k3) && rep_quota>0) {
									C_ringmember.at(k3)=-1;
									rep_quota--;
								}
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


int MOLECULE::decyc_small_ring(int size) {

    vector<int> C_ringmember(0);
	
    for (int k1=if_circle;k1>=1;k1--) { // k1>=ct+1
        C_ringmember.resize(0);
        int rep_quota=1;
        for (int k2=Cyindex.size()-1;k2>=0;k2--) {
            if (Cyindex.at(k2).size()) {
                for (int k3=0;k3<Cyindex.at(k2).size();k3++) {
                    if (Cyindex.at(k2).at(k3)==k1) {
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
		if (C_ringmember.size()>0) {
			int num_rmember=0;
        	for (int k2=0;k2<C_ringmember.size();k2++) {
            	if (C_ringmember.at(k2)!=-1) num_rmember++;
        	}
			if (num_rmember<size) {
				for (int k2=0;k2<Cyindex.size();k2++) {
					if (Cyindex.at(k2).size()) {
                		for (int k3=0;k3<Cyindex.at(k2).size();k3++) {
                    		if (Cyindex.at(k2).at(k3)==k1) {
								Cyindex.at(k2).erase(Cyindex.at(k2).begin()+k3);
								k3--;
							}
						}
					}
				}

                for (int k2=0;k2<Cyindex.size();k2++) {
                    if (Cyindex.at(k2).size()) {
                        for (int k3=0;k3<Cyindex.at(k2).size();k3++) {
                            if (Cyindex.at(k2).at(k3)>k1) {
                                Cyindex.at(k2).at(k3)--;
                            }
                        }
                    }
                }

				if (Cybnd.size()>=k1) Cybnd.erase(Cybnd.begin()+k1-1);
				if_circle--;
			}
		}
    }
    vector<int>().swap(C_ringmember);
	

	//for (int i=0;i<Cybnd.size();i++) { //i=ct
	//	if (!Cybnd.at(i)) {
	//		Cybnd.erase(Cybnd.begin()+i);
	//		i--;
	//	}
	//}

	return 1;

}


int MOLECULE::del_unpaired_ring_no(){
	int k;
    vector<int> freq(if_circle,0);
	if (if_circle) {
	    for (k=0;k<Cyindex.size();k++) {  // 20200131 prevent unpaired ring number
	        if (Cyindex.at(k).size()) {
	            for (int k1=0;k1<Cyindex.at(k).size();k1++) {
					freq.at(Cyindex.at(k).at(k1)-1)++;
	            }
	        }
	    }

	    for (k=freq.size()-1;k>=0;k--) { //k>=ct
			if (!freq.at(k) || freq.at(k)%2) {

                for (int k1=0;k1<Cyindex.size();k1++) {
                    if (Cyindex.at(k1).size()) {
                        for (int k2=0;k2<Cyindex.at(k1).size();k2++) {
                            if (Cyindex.at(k1).at(k2)==k+1) {
								Cyindex.at(k1).erase(Cyindex.at(k1).begin()+k2);
								k2--;
							}
                        }
                    }
                }

                for (int k1=0;k1<Cyindex.size();k1++) {
                    if (Cyindex.at(k1).size()) {
                        for (int k2=0;k2<Cyindex.at(k1).size();k2++) {
							if (Cyindex.at(k1).at(k2)>k+1) {
								Cyindex.at(k1).at(k2)--;
							}
                        }
                    }
                }

				if (Cybnd.size()>k) Cybnd.erase(Cybnd.begin()+k);

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

    stringstream nu("");
    if (0) {
        stringstream ss(molesmi);
        stringstream ss_out("");

        OBConversion conv(&ss, &ss_out);
        if(conv.SetInAndOutFormats("SMI","MOL"))
        {
            //conv.AddOption("h", OBConversion::GENOPTIONS);
            conv.AddOption("gen3D", OBConversion::GENOPTIONS);
            conv.AddOption("3", OBConversion::OUTOPTIONS);
            //conv.AddOption("canonical", OBConversion::GENOPTIONS);
            //conv.AddOption("minimize", OBConversion::GENOPTIONS);
            //conv.AddOption("ff", OBConversion::GENOPTIONS,"uff");
            //conv.AddOption("steps", OBConversion::GENOPTIONS,"1");
            //obErrorLog.StopLogging();
            obErrorLog.SetOutputLevel(obMessageLevel::obError);
            //obErrorLog.SetOutputLevel(obMessageLevel::obWarning);
            obErrorLog.SetOutputStream(&cout);
            conv.Convert();
            //cout << molesmi << endl;
        }
        //smiles=molesmi=ss_out.str();
        //ss_out >> molesmi >> ws;
        //smiles=molesmi;
        nu << ss_out.str();
    }

    if (0) {
        stringstream ss(nu.str());
        stringstream ss_out("");

        OBConversion conv(&ss, &ss_out);
        if(conv.SetInAndOutFormats("MOL","SMI"))
        {
            //conv.AddOption("h", OBConversion::GENOPTIONS);
            conv.AddOption("gen3D", OBConversion::GENOPTIONS);
            //conv.AddOption("3", OBConversion::OUTOPTIONS);
            //conv.AddOption("minimize", OBConversion::GENOPTIONS);
            //conv.AddOption("ff", OBConversion::GENOPTIONS,"uff");
            //conv.AddOption("steps", OBConversion::GENOPTIONS,"1");
            conv.AddOption("canonical", OBConversion::GENOPTIONS);
            //conv.AddOption("c", OBConversion::OUTOPTIONS);
            //obErrorLog.StopLogging();
            obErrorLog.SetOutputLevel(obMessageLevel::obError);
            //obErrorLog.SetOutputLevel(obMessageLevel::obWarning);
            obErrorLog.SetOutputStream(&cout);
            conv.Convert();
            //cout << molesmi << endl;
        }
        //smiles=molesmi=ss_out.str();
        ss_out >> molesmi >> ws;
        smiles=molesmi;
        //nu << ss_out.str();
    }


	if (1) {
		stringstream ss(molesmi);
    	stringstream ss_out("");

    	OBConversion conv(&ss, &ss_out);
    	if(conv.SetInAndOutFormats("SMI","SMI"))
    	{
        	//conv.AddOption("h", OBConversion::GENOPTIONS);
        	//conv.AddOption("gen3D", OBConversion::GENOPTIONS);
			//conv.AddOption("gen2D", OBConversion::GENOPTIONS);
        	conv.AddOption("canonical", OBConversion::GENOPTIONS);
            obErrorLog.SetOutputLevel(obMessageLevel::obError);
            //obErrorLog.SetOutputLevel(obMessageLevel::obWarning);
            obErrorLog.SetOutputStream(&cout);
        	conv.Convert();
    	}
		//smiles=molesmi=ss_out.str();
		ss_out >> molesmi >> ws;
		smiles=molesmi;
	}

	return 1;
}

