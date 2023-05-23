#include "MOLECULE.h"
//#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
using namespace std;
using namespace OpenBabel;

extern PARAMETER para;

unsigned int MOLECULE::replace(MOLECULE &mol) {
	if (mol.Cindex.size()==0) {
		cout<<"Error: It cannot be replaced by an empty molecule. smiles = " << mol.smiles << " | molesmi = " << mol.molesmi << endl;
		return 0;
	}
	else {
		unsigned int i=0,j=0;
		empty();
		clear();

		ctsisomer.resize(2,vector<string>(0));
		Cyindex.resize(mol.Cindex.size(),vector<unsigned int>(0));
		for (i=0;i<mol.Cindex.size();i++) {
			Cindex.push_back(mol.Cindex.at(i));
			Pindex.push_back(mol.Pindex.at(i));
			Rindex.push_back(mol.Rindex.at(i));
			chi.push_back(mol.chi.at(i));
			for (j=0;j<mol.Cyindex.at(i).size();j++) Cyindex.at(i).push_back(mol.Cyindex.at(i).at(j));
			Mindex.push_back(mol.Mindex.at(i));
			if (para.protect) protect.push_back(mol.protect.at(i));
            ctsisomer.at(0).push_back(mol.ctsisomer.at(0).at(i));
            ctsisomer.at(1).push_back(mol.ctsisomer.at(1).at(i));
		}
		for (i=0;i<mol.Cybnd.size();i++) Cybnd.push_back(mol.Cybnd.at(i));
		if_circle=mol.if_circle;

		natom=mol.natom;
		nsubcomp=1;
		comp_id=mol.comp_id;
		frac=mol.frac;
		chg=mol.chg;
		multiplicity=mol.multiplicity;

		data=mol.data;

		reset();

        smiles=mol.smiles;
        molesmi=mol.molesmi;
		canonicalize_SMILES();
	}
	return 1;
}

unsigned int MOLECULE::crossover(MOLECULE &aaa,unsigned int pp,unsigned int jj,bool cistrans) {
	if (para.protect) {
		if (protect.at(pp)) return 0; 
		if (aaa.protect.at(jj)) return 0; 
	}
    unsigned int Csize1=Cindex.size();
    for (unsigned int i=0;i<Cindex.size();i++) {
		if (Mindex.at(i)==70) Csize1--;
	}
    unsigned int Csize2=aaa.Cindex.size();
    for (unsigned int i=0;i<aaa.Cindex.size();i++) {
        if (aaa.Mindex.at(i)==70) Csize2--;
    }
    if (Csize1<2 || pp>=Cindex.size()) return 0;
    if (Csize2<2 || jj>=aaa.Cindex.size()) return 0;

    //if (Cindex.size()<2 || pp>=Cindex.size()) return 0;
    //if (aaa.Cindex.size()<2 || jj>=aaa.Cindex.size()) return 0;

	int lenA=0;
	int lenB=0;
	vector<unsigned int> t(0),f(0);
	vector<unsigned int> ci_ref(0),mi_ref(0),pi_ref(0),ri_ref(0),pri(0),cybndi(0),chii(0); 
	vector<unsigned int> cj_ref(0),mj_ref(0),pj_ref(0),rj_ref(0),prj(0),cybndj(0),chij(0); 
	vector< vector<unsigned int> > cyi(0,vector<unsigned int>(0));
	vector< vector<unsigned int> > cyj(0,vector<unsigned int>(0));
	vector< vector<string> > cti_ref(0),ctj_ref(0);
	
	vector<bool> ifchg(4,0);
	if (para.ion && 0) {
		for (unsigned int k1=0;k1<Cindex.size();k1++) {
			if (Cindex.at(k1)<=Cindex.at(pp) && data->a[Mindex.at(k1)].chg) {
				ifchg.at(0)=1;
				break;
			}
		}
		for (unsigned int k1=0;k1<Cindex.size();k1++) {
			if (Cindex.at(k1)>=Cindex.at(pp) && data->a[Mindex.at(k1)].chg) {
				ifchg.at(1)=1;
				break;
			}
		}
		for (unsigned int k1=0;k1<aaa.Cindex.size();k1++) {
			if (aaa.Cindex.at(k1)<=aaa.Cindex.at(jj) && data->a[aaa.Mindex.at(k1)].chg) {
				ifchg.at(2)=1;
				break;
			}
		}
		for (unsigned int k1=0;k1<aaa.Cindex.size();k1++) {
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
        for (unsigned int k1=0;k1<Cindex.size();k1++) {
            if (Cindex.at(k1)<=Cindex.at(pp) && Cyindex.at(k1).size()) { // && Cyindex.at(k1)
                ifcyc.at(0)=1;
                break;
            }
        }
        for (unsigned int k1=0;k1<Cindex.size();k1++) {
            if (Cindex.at(k1)>=Cindex.at(pp) && Cyindex.at(k1).size()) { // && Cyindex.at(k1)
                ifcyc.at(1)=1;
                break;
            }
        }
        for (unsigned int k1=0;k1<aaa.Cindex.size();k1++) {
            if (aaa.Cindex.at(k1)<=aaa.Cindex.at(jj) && aaa.Cyindex.at(k1).size()) { // && aaa.Cyindex.a(k1)
                ifcyc.at(2)=1;
                break;
            }
        }
        for (unsigned int k1=0;k1<aaa.Cindex.size();k1++) {
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

	if (aaa.Rindex.at(jj)==Rindex.at(pp)) {
		ri_ref.push_back(Rindex.at(pp));
		rj_ref.push_back(aaa.Rindex.at(jj));
		ci_ref.push_back(Cindex.at(pp));
		cj_ref.push_back(aaa.Cindex.at(jj));
		pi_ref.push_back(Pindex.at(pp));
		pj_ref.push_back(aaa.Pindex.at(jj));
		mi_ref.push_back(Mindex.at(pp));
		mj_ref.push_back(aaa.Mindex.at(jj));
        chii.push_back(chi.at(pp));
        chij.push_back(aaa.chi.at(jj));

        cti_ref.resize(2,vector<string>(0));
        cti_ref.at(0).push_back(ctsisomer.at(0).at(pp));
        cti_ref.at(1).push_back(ctsisomer.at(1).at(pp));

        ctj_ref.resize(2,vector<string>(0));
        ctj_ref.at(0).push_back(aaa.ctsisomer.at(0).at(jj));
        ctj_ref.at(1).push_back(aaa.ctsisomer.at(1).at(jj));

		cyi.resize(1,vector<unsigned int>(0));
		for (unsigned int k1=0;k1<Cyindex.at(pp).size();k1++) cyi.at(0).push_back(Cyindex.at(pp).at(k1));

		cyj.resize(1,vector<unsigned int>(0));
		for (unsigned int k1=0;k1<aaa.Cyindex.at(jj).size();k1++) cyj.at(0).push_back(aaa.Cyindex.at(jj).at(k1));

		if (para.protect) pri.push_back(protect.at(pp));
		if (para.protect) prj.push_back(aaa.protect.at(jj));

		lenA=Cindex.size();
		lenB=aaa.Cindex.size();
		if (1) {
			for (unsigned int n=1;n<Cindex.size();n++) {
				if (Cindex.at(n) > ci_ref.at(0)) {
					unsigned int ref = ci_ref.size();
					for (unsigned int m=0;m<ref;m++) {
						if (Pindex.at(n) == ci_ref.at(m)) {
							ci_ref.push_back(Cindex.at(n));
							mi_ref.push_back(Mindex.at(n));
							pi_ref.push_back(Pindex.at(n));
							ri_ref.push_back(Rindex.at(n));
							chii.push_back(chi.at(n));
                        	cti_ref.at(0).push_back(ctsisomer.at(0).at(n));
                        	cti_ref.at(1).push_back(ctsisomer.at(1).at(n));

							cyi.resize(cyi.size()+1,vector<unsigned int>(0));
							for (unsigned int k1=0;k1<Cyindex.at(n).size();k1++) cyi.at(cyi.size()-1).push_back(Cyindex.at(n).at(k1));
						
							if (para.protect) pri.push_back(protect.at(n));
						}
					}
				}
			}
		}
		else {
	    	vector<unsigned int> curatm(0);

	    	bool goout=0;
    		curatm.push_back(ci_ref.at(0));

    		do {
        		vector<unsigned int> tmpp(0);
        		for (unsigned int k3=0;k3<curatm.size();k3++) {
            		for (unsigned int k2=0;k2<Cindex.size();k2++) {
                		if (Pindex.at(k2)==curatm.at(k3)) {
	                    	ci_ref.push_back(Cindex.at(k2));
    	                	mi_ref.push_back(Mindex.at(k2));
        	            	pi_ref.push_back(Pindex.at(k2));
            	        	ri_ref.push_back(Rindex.at(k2));
							chii.push_back(chii.at(k2));
                    		cti_ref.at(0).push_back(ctsisomer.at(0).at(k2));
                        	cti_ref.at(1).push_back(ctsisomer.at(1).at(k2));

                        	cyi.resize(cyi.size()+1,vector<unsigned int>(0));
							for (unsigned int k1=0;k1<Cyindex.at(k2).size();k1++) cyi.at(cyi.size()-1).push_back(Cyindex.at(k2).at(k1));

                        	if (para.protect) pri.push_back(protect.at(k2));

							tmpp.push_back(Cindex.at(k2));
                		}
            		}
        		}
        		if (tmpp.size()>0) {
            		curatm.resize(0);
            		for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
        		}
        		else {
            		goout=1;
            		break;
        		}

    		} while (!goout);

			for (unsigned int i=0;i<ci_ref.size();i++) {
				for (unsigned int j=i+1;j<ci_ref.size();j++) {
					if (ci_ref.at(i)>ci_ref.at(j)) {
						swap(ci_ref.at(i),ci_ref.at(j));
						swap(mi_ref.at(i),mi_ref.at(j));
						swap(pi_ref.at(i),pi_ref.at(j));
						swap(ri_ref.at(i),ri_ref.at(j));
						swap(chii.at(i),chii.at(j));
						swap(cti_ref.at(0).at(i),cti_ref.at(0).at(j));
						swap(cti_ref.at(1).at(i),cti_ref.at(1).at(j));
						swap(cyi.at(i),cyi.at(j));
						if (para.protect) swap(pri.at(i),pri.at(j));
					}
				}
			}
		}

		if (1) {
			for (unsigned int m=1;m<aaa.Cindex.size();m++) {
				if (aaa.Cindex.at(m) > cj_ref.at(0) ) {
					unsigned int ref = cj_ref.size();
					for (unsigned int n=0;n<ref;n++) {
						if (aaa.Pindex.at(m) == cj_ref.at(n)) {
							cj_ref.push_back(aaa.Cindex.at(m));
							mj_ref.push_back(aaa.Mindex.at(m));
							pj_ref.push_back(aaa.Pindex.at(m));
							rj_ref.push_back(aaa.Rindex.at(m));
							chij.push_back(aaa.chi.at(m));
                        	ctj_ref.at(0).push_back(aaa.ctsisomer.at(0).at(m));
                        	ctj_ref.at(1).push_back(aaa.ctsisomer.at(1).at(m));

                        	cyj.resize(cyj.size()+1,vector<unsigned int>(0));
                        	for (unsigned int k1=0;k1<aaa.Cyindex.at(m).size();k1++) cyj.at(cyj.size()-1).push_back(aaa.Cyindex.at(m).at(k1));
						
							if (para.protect) prj.push_back(aaa.protect.at(m));
						}
					}
				}
			}
		}
		else {
            vector<unsigned int> curatm(0);
            bool goout=0;
            curatm.push_back(ci_ref.at(0));

        	do {
            	vector<unsigned int> tmpp(0);
            	for (unsigned int k3=0;k3<curatm.size();k3++) {
                	for (unsigned int k2=0;k2<aaa.Cindex.size();k2++) {
                    	if (aaa.Pindex.at(k2)==curatm.at(k3)) {
                        	cj_ref.push_back(aaa.Cindex.at(k2));
                        	mj_ref.push_back(aaa.Mindex.at(k2));
                        	pj_ref.push_back(aaa.Pindex.at(k2));
                        	rj_ref.push_back(aaa.Rindex.at(k2));
							chij.push_back(aaa.chi.at(k2));
                        	ctj_ref.at(0).push_back(aaa.ctsisomer.at(0).at(k2));
                        	ctj_ref.at(1).push_back(aaa.ctsisomer.at(1).at(k2));
                        	cyj.resize(cyj.size()+1,vector<unsigned int>(0));
                        	for (unsigned int k1=0;k1<aaa.Cyindex.at(k2).size();k1++) cyj.at(cyj.size()-1).push_back(aaa.Cyindex.at(k2).at(k1));

                        	if (para.protect) prj.push_back(aaa.protect.at(k2));
						
							tmpp.push_back(aaa.Cindex.at(k2));
                    	}
                	}
            	}
            	if (tmpp.size()>0) {
                	curatm.resize(0);
                	for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
            	}
            	else {
                	goout=1;
                	break;
            	}

        	} while (!goout);

        	for (unsigned int i=0;i<cj_ref.size();i++) {
            	for (unsigned int j=i+1;j<cj_ref.size();j++) {
                	if (cj_ref.at(i)>cj_ref.at(j)) {
                    	swap(cj_ref.at(i),cj_ref.at(j));
                    	swap(mj_ref.at(i),mj_ref.at(j));
                    	swap(pj_ref.at(i),pj_ref.at(j));
                    	swap(rj_ref.at(i),rj_ref.at(j));
						swap(chij.at(i),chij.at(j));
                    	swap(ctj_ref.at(0).at(i),ctj_ref.at(0).at(j));
                    	swap(ctj_ref.at(1).at(i),ctj_ref.at(1).at(j));
                    	swap(cyj.at(i),cyj.at(j));
                    	if (para.protect) swap(prj.at(i),prj.at(j));
                	}
            	}
        	}
		}

		if (1) {
			for (unsigned int n=0;n<pri.size();n++) if (pri.at(n)) return 0;
			for (unsigned int n=0;n<prj.size();n++) if (prj.at(n)) return 0;
    		if (para.ion) {
				/*
				bool michg=0,mjchg=0;
    			for (unsigned int k1=0;k1<mi_ref.size();k1++) {
					if (data->a[mi_ref.at(k1)].chg) {
						michg=1;
						break;
					}
				}
            	for (unsigned int k1=0;k1<mj_ref.size();k1++) {
                	if (data->a[mj_ref.at(k1)].chg) {
                    	mjchg=1;
                    	break;
                	}
            	}
				if (michg!=mjchg) return 0;
				*/
				double michg=0,mjchg=0;
				double mirefchg=0,mjrefchg=0;
				for (unsigned int k1=0;k1<Mindex.size();k1++) michg+=data->a.at(Mindex.at(k1)).chg;
				for (unsigned int k1=0;k1<aaa.Mindex.size();k1++) mjchg+=aaa.data->a.at(aaa.Mindex.at(k1)).chg;
				for (unsigned int k1=0;k1<mi_ref.size();k1++) mirefchg+=data->a.at(mi_ref.at(k1)).chg;
				for (unsigned int k1=0;k1<mj_ref.size();k1++) mjrefchg+=data->a.at(mj_ref.at(k1)).chg;

				if ((michg-mirefchg)*mjrefchg<0) return 0;
				if ((mjchg-mjrefchg)*mirefchg<0) return 0;
		
			}
		}

		if (1) {
			//cybndi.resize(if_circle+aaa.if_circle,0);
			cybndi.resize(aaa.if_circle,0);
			for (unsigned int n=0;n<aaa.if_circle;n++) cybndi.at(n)=aaa.Cybnd.at(n);
		
			//cybndj.resize(if_circle+aaa.if_circle,0);
			cybndj.resize(if_circle,0);
			for (unsigned int n=0;n<if_circle;n++) cybndj.at(n)=Cybnd.at(n);
		
			int tm1=if_circle;
			for (unsigned int n=0;n<cyi.size();n++) {
				for (unsigned int k1=0;k1<cyi.at(n).size();k1++) {
					cyi.at(n).at(k1)+=aaa.if_circle;
				}
			}
			for (unsigned int n=0;n<cyj.size();n++) {
                for (unsigned int k1=0;k1<cyj.at(n).size();k1++) {
                    cyj.at(n).at(k1)+=if_circle;
                }
			}
            if_circle+=aaa.if_circle;
            aaa.if_circle+=tm1;
		}

		unsigned int x=ci_ref.at(0);
		unsigned int y=cj_ref.at(0);
		unsigned int a=pi_ref.at(0);
		unsigned int b=pj_ref.at(0);

		if (1) {
			for (unsigned int n=0;n<ci_ref.size();n++) {
				for (unsigned int m=0;m<Cindex.size();m++) {
					if (Cindex.at(m) == ci_ref.at(n)) {
						Mindex.erase(Mindex.begin()+m);
						Cindex.erase(Cindex.begin()+m);
						Pindex.erase(Pindex.begin()+m);
						Rindex.erase(Rindex.begin()+m);
						chi.erase(chi.begin()+m);
						Cyindex.erase(Cyindex.begin()+m);
                    	ctsisomer.at(0).erase(ctsisomer.at(0).begin()+m);
                    	ctsisomer.at(1).erase(ctsisomer.at(1).begin()+m);
						if (para.protect) protect.erase(protect.begin()+m);

						break;
					}
				}
			}
		}
		else {
			for (int n=ci_ref.size()-1;n>=0;n--) {
				Mindex.erase(Mindex.begin()+ci_ref.at(n)-1);
				Cindex.erase(Cindex.begin()+ci_ref.at(n)-1);
				Pindex.erase(Pindex.begin()+ci_ref.at(n)-1);
				Rindex.erase(Rindex.begin()+ci_ref.at(n)-1);
				chi.erase(chi.begin()+ci_ref.at(n)-1);
				Cyindex.erase(Cyindex.begin()+ci_ref.at(n)-1);
				ctsisomer.at(0).erase(ctsisomer.at(0).begin()+ci_ref.at(n)-1);
				ctsisomer.at(1).erase(ctsisomer.at(1).begin()+ci_ref.at(n)-1);
				if (para.protect) protect.erase(protect.begin()+ci_ref.at(n)-1);
			}
		}


		if (1) {
			for (unsigned int m=0;m<cj_ref.size();m++) {
				for (unsigned int n=0;n<aaa.Cindex.size();n++) {
					if (aaa.Cindex.at(n) == cj_ref.at(m)) {
						aaa.Mindex.erase(aaa.Mindex.begin()+n);
						aaa.Cindex.erase(aaa.Cindex.begin()+n);
						aaa.Pindex.erase(aaa.Pindex.begin()+n);
						aaa.Rindex.erase(aaa.Rindex.begin()+n);
						aaa.chi.erase(aaa.chi.begin()+n);
						aaa.Cyindex.erase(aaa.Cyindex.begin()+n);
                    	aaa.ctsisomer.at(0).erase(aaa.ctsisomer.at(0).begin()+n);
                    	aaa.ctsisomer.at(1).erase(aaa.ctsisomer.at(1).begin()+n);
						if (para.protect) aaa.protect.erase(aaa.protect.begin()+n);

						break;
					}
				}
			}
		}
		else {
        	for (int n=cj_ref.size()-1;n>=0;n--) {
            	aaa.Mindex.erase(aaa.Mindex.begin()+cj_ref.at(n)-1);
            	aaa.Cindex.erase(aaa.Cindex.begin()+cj_ref.at(n)-1);
            	aaa.Pindex.erase(aaa.Pindex.begin()+cj_ref.at(n)-1);
            	aaa.Rindex.erase(aaa.Rindex.begin()+cj_ref.at(n)-1);
				aaa.chi.erase(aaa.chi.begin()+cj_ref.at(n)-1);
            	aaa.Cyindex.erase(aaa.Cyindex.begin()+cj_ref.at(n)-1);
            	aaa.ctsisomer.at(0).erase(aaa.ctsisomer.at(0).begin()+cj_ref.at(n)-1);
            	aaa.ctsisomer.at(1).erase(aaa.ctsisomer.at(1).begin()+cj_ref.at(n)-1);
            	if (para.protect) aaa.protect.erase(aaa.protect.begin()+cj_ref.at(n)-1);
        	}
		}

		if (1) {
			t.clear();
			for (unsigned int n=1;n<ci_ref.size();n++) {
				for (unsigned int m=0;m<ci_ref.size();m++) {
					if (pi_ref.at(n) == ci_ref.at(m)) {
						t.push_back(m);
					}
				}
			}
			f.clear();
			for (unsigned int n=1;n<cj_ref.size();n++) {
				for (unsigned int m=0;m<cj_ref.size();m++) {
					if (pj_ref.at(n) == cj_ref.at(m)) {
						f.push_back(m);
					}
				}
			}

			pi_ref.at(0) = b;
			pj_ref.at(0) = a;
			ci_ref.at(0) = y;
			cj_ref.at(0) = x;

			for (unsigned int n=1;n<ci_ref.size();n++) {
				ci_ref.at(n) = ci_ref.at(n-1) + 1;
				pi_ref.at(n) = ci_ref.at(t.at(n-1));
			}
			for (unsigned int n=1;n<cj_ref.size();n++) {
				cj_ref.at(n) = cj_ref.at(n-1) + 1;
				pj_ref.at(n) = cj_ref.at(f.at(n-1));
			}

			t.clear();
			f.clear();
			for (unsigned int n=1;n<Cindex.size();n++) {
				for (unsigned int m=0;m<Cindex.size();m++) {
					if (Pindex.at(n) == Cindex.at(m)) {
						t.push_back(m);
					}
				}
			}
			unsigned int ref=0;
			for (unsigned int n=1;n<Cindex.size();n++) {
				if (Cindex.at(n) != Cindex.at(n-1)+1) {
					ref=n;
					Cindex.at(ref) = Cindex.at(ref-1) + cj_ref.size() + 1;
					break;
				}
			}
			if (ref) {
				for (unsigned int n=ref+1;n<Cindex.size();n++) {
					Cindex.at(n) = Cindex.at(n-1) + 1;
				}

				for (unsigned int n=1;n<Cindex.size();n++) {
					Pindex.at(n) = Cindex.at(t.at(n-1));
				}
			}
			t.clear();
			f.clear();
			for (unsigned int n=1;n<aaa.Cindex.size();n++) {
				for (unsigned int m=0;m<aaa.Cindex.size();m++) {
					if (aaa.Pindex.at(n) == aaa.Cindex.at(m)) {
						t.push_back(m);
					}
				}
			}
			ref=0;
			for (unsigned int n=1;n<aaa.Cindex.size();n++) {
				if (aaa.Cindex.at(n) != aaa.Cindex.at(n-1)+1) {
					ref = n;
					aaa.Cindex.at(ref) = aaa.Cindex.at(ref-1) + ci_ref.size() + 1;
					break;
				}
			}
			if (ref) {
				for (unsigned int n=ref+1;n<aaa.Cindex.size();n++) {
					aaa.Cindex.at(n) = aaa.Cindex.at(n-1) + 1;
				}

				for (unsigned int n=1;n<aaa.Cindex.size();n++) {
					aaa.Pindex.at(n) = aaa.Cindex.at(t.at(n-1));
				}
			}
		}
		else {
	        for (unsigned int n=0;n<ci_ref.size();n++) {
    	        ci_ref.at(n)+=(lenB+2);
        	    pi_ref.at(n)+=(lenB+2);
        	}
        	for (unsigned int n=0;n<cj_ref.size();n++) {
            	cj_ref.at(n)+=(lenA+2);
            	pj_ref.at(n)+=(lenA+2);
        	}
		}
	}
	else return 0;


	if (1) {
		unsigned int n=0;
		while(n<cj_ref.size()) {
			for (unsigned int m=0;m<Cindex.size();m++) {
				if (Cindex.at(m) > cj_ref.at(n) && m >= 1) { // Cindex[m] > cj_ref[n] && Cindex[m-1] < cj_ref[n]  orig // 20191130
					if (Cindex.at(m-1) < cj_ref.at(n)) {
						Cindex.insert(Cindex.begin() + m, cj_ref.at(n));
						Pindex.insert(Pindex.begin() + m, pj_ref.at(n));
						Mindex.insert(Mindex.begin() + m, mj_ref.at(n));
						Rindex.insert(Rindex.begin() + m, rj_ref.at(n));
						chi.insert(chi.begin() + m, chij.at(n));
						Cyindex.insert(Cyindex.begin() + m, cyj.at(n));
                    	ctsisomer.at(0).insert(ctsisomer.at(0).begin() + m, ctj_ref.at(0).at(n));
                    	ctsisomer.at(1).insert(ctsisomer.at(1).begin() + m, ctj_ref.at(1).at(n));
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
					chi.push_back(chij.at(n));
					Cyindex.resize(Cyindex.size()+1,vector<unsigned int>(0));
					for (unsigned int k1=0;k1<cyj.at(n).size();k1++) Cyindex.at(Cyindex.size()-1).push_back(cyj.at(n).at(k1));
                	ctsisomer.at(0).push_back(ctj_ref.at(0).at(n));
                	ctsisomer.at(1).push_back(ctj_ref.at(1).at(n));
					if (para.protect) protect.push_back(prj.at(n));
					n++;
					break;
				}
			}
		}
		unsigned int m=0;
		unsigned int ref=1;
		while (m<ci_ref.size()) {
			for (unsigned int n=0;n<aaa.Cindex.size();n++) {
				if (aaa.Cindex.at(n) > ci_ref.at(m) && n >= 1) { // aaa.Cindex[n] > ci_ref[m] && aaa.Cindex[n-1] < ci_ref[m]  orig  //20191130
					if (aaa.Cindex.at(n-1) < ci_ref.at(m)) {
						ref=n;
						aaa.Cindex.insert(aaa.Cindex.begin() + ref, ci_ref.at(m));
						aaa.Pindex.insert(aaa.Pindex.begin() + ref, pi_ref.at(m));
						aaa.Mindex.insert(aaa.Mindex.begin() + ref, mi_ref.at(m));
						aaa.Rindex.insert(aaa.Rindex.begin() + ref, ri_ref.at(m));
						aaa.chi.insert(aaa.chi.begin() + ref, chii.at(m));
						aaa.Cyindex.insert(aaa.Cyindex.begin() + ref, cyi.at(m));
                    	aaa.ctsisomer.at(0).insert(aaa.ctsisomer.at(0).begin() + ref, cti_ref.at(0).at(m));
                    	aaa.ctsisomer.at(1).insert(aaa.ctsisomer.at(1).begin() + ref, cti_ref.at(1).at(m));
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
					aaa.chi.push_back(chii.at(m));
                	aaa.Cyindex.resize(aaa.Cyindex.size()+1,vector<unsigned int>(0));
                	for (unsigned int k1=0;k1<cyi.at(m).size();k1++) aaa.Cyindex.at(aaa.Cyindex.size()-1).push_back(cyi.at(m).at(k1));
                	aaa.ctsisomer.at(0).push_back(cti_ref.at(0).at(m));
                	aaa.ctsisomer.at(1).push_back(cti_ref.at(1).at(m));
					if (para.protect) aaa.protect.push_back(pri.at(m));
					m++;
					break;
				}
			}
		}
	}
	else {
    	for (unsigned int n=0;n<cj_ref.size();n++) {
        	Cindex.push_back(cj_ref.at(n));
        	Pindex.push_back(pj_ref.at(n));
        	Mindex.push_back(mj_ref.at(n));
        	Rindex.push_back(rj_ref.at(n));
			chi.push_back(chij.at(n));
        	Cyindex.resize(Cyindex.size()+1,vector<unsigned int>(0));
        	for (unsigned int k1=0;k1<cyj.at(n).size();k1++) Cyindex.at(Cyindex.size()-1).push_back(cyj.at(n).at(k1));
        	ctsisomer.at(0).push_back(ctj_ref.at(0).at(n));
        	ctsisomer.at(1).push_back(ctj_ref.at(1).at(n));
        	if (para.protect) protect.push_back(prj.at(n));
    	}

    	for (unsigned int n=0;n<ci_ref.size();n++) {
        	aaa.Cindex.push_back(ci_ref.at(n));
        	aaa.Pindex.push_back(pi_ref.at(n));
        	aaa.Mindex.push_back(mi_ref.at(n));
        	aaa.Rindex.push_back(ri_ref.at(n));
			aaa.chi.push_back(chii.at(n));
        	aaa.Cyindex.resize(aaa.Cyindex.size()+1,vector<unsigned int>(0));
        	for (unsigned int k1=0;k1<cyi.at(n).size();k1++) aaa.Cyindex.at(aaa.Cyindex.size()-1).push_back(cyi.at(n).at(k1));
        	aaa.ctsisomer.at(0).push_back(cti_ref.at(0).at(n));
        	aaa.ctsisomer.at(1).push_back(cti_ref.at(1).at(n));
        	if (para.protect) aaa.protect.push_back(pri.at(n));
    	}

    	for (unsigned int i=0;i<Cindex.size();i++) {
        	unsigned int j=Cindex.at(i);
        	Cindex.at(i)=(i+1);
        	for (unsigned int  k1=0;k1<Cindex.size();k1++) {
            	if (Pindex.at(k1)==j) Pindex.at(k1)=Cindex.at(i);
        	}
    	}

    	for (unsigned int i=0;i<aaa.Cindex.size();i++) {
        	unsigned int j=aaa.Cindex.at(i);
        	aaa.Cindex.at(i)=(i+1);
        	for (unsigned int  k1=0;k1<aaa.Cindex.size();k1++) {
            	if (aaa.Pindex.at(k1)==j) aaa.Pindex.at(k1)=aaa.Cindex.at(i);
        	}
    	}

	}

	//Cybnd.resize(0);
	//for (unsigned int n=0;n<cybndj.size();n++) Cybnd.push_back(cybndj.at(n));
    //aaa.Cybnd.resize(0);
    //for (unsigned int n=0;n<cybndi.size();n++) aaa.Cybnd.push_back(cybndi.at(n));

	for (unsigned int n=0;n<cybndi.size();n++) Cybnd.push_back(cybndi.at(n));
	for (unsigned int n=0;n<cybndj.size();n++) aaa.Cybnd.push_back(cybndj.at(n));

    vector<unsigned int>().swap(t);
    vector<unsigned int>().swap(f);
	vector<unsigned int>().swap(pi_ref);
    vector<unsigned int>().swap(ci_ref);
    vector<unsigned int>().swap(mi_ref);
    vector<unsigned int>().swap(ri_ref);
	vector<unsigned int>().swap(chii);
	vector<unsigned int>().swap(pj_ref);
    vector<unsigned int>().swap(cj_ref);
    vector<unsigned int>().swap(mj_ref);
    vector<unsigned int>().swap(rj_ref);
	vector<unsigned int>().swap(chij);
    vector<unsigned int>().swap(prj);
    vector<unsigned int>().swap(pri);
    vector< vector<unsigned int> >().swap(cyj);
    vector< vector<unsigned int> >().swap(cyi);
    vector<unsigned int>().swap(cybndj);
    vector<unsigned int>().swap(cybndi);

	del_unpaired_ring_no();
	aaa.del_unpaired_ring_no();
	//decyc_small_ring(5);
	//aaa.decyc_small_ring(5);

	vector<unsigned int> PCofct(0);
    for (unsigned int i=0;i<Cindex.size();i++) {
		if (Rindex.at(i)==2) {
			PCofct.push_back(Pindex.at(i));
			PCofct.push_back(Cindex.at(i));
		}
    }
	for (unsigned int i=0;i<PCofct.size();i+=2) {
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

	if (!chk_imine_ct(0,Cindex.size()-1,cistrans)) mds2smi();
	if (!aaa.chk_imine_ct(0,aaa.Cindex.size()-1,cistrans)) aaa.mds2smi();

	//chk_imine_ct(0,Cindex.size()-1,cistrans);
	//mds2smi();
	//aaa.chk_imine_ct(0,aaa.Cindex.size()-1,cistrans);
	//aaa.mds2smi();

	return 1;
}


unsigned int MOLECULE::read(string a) {
	for (unsigned int i=0;i<a.length();i++) {
		if (a[i]=='/') a[i]='u';
		if (a[i]=='\\') a[i]='d';
		if (a[i]=='*') a[i]='x';
	}
	
	ifstream inf((a+".mds").c_str());
	string b,b1;
	unsigned int i=0,j=0,k=0;
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
			Cyindex.resize(j,vector<unsigned int>(0));
			for (i=0;i<j;i++) {
                b1="";
                inf>>b1>>ws;

                string buf="";
                for (unsigned int k1=0;k1<b1.length();k1++) {
                    if (b1[k1]!=',') buf+=b1[k1];
                    if (b1[k1]==',') {
                        if (buf!="") {
							Cyindex.at(i).push_back(atoi(buf.c_str()));
						}
                        buf="";
                    }
					if (b1[k1]!=',' && k1>=b1.length()-1) {
						if (buf!="0") {
							Cyindex.at(i).push_back(atoi(buf.c_str()));
						}
						buf="";
					}
                }
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
            for (i=0;i<if_circle;i++) {
                inf>>k>>ws;
                Cybnd.push_back(k);
            }
        }
		else if (b=="Chirality") {
			for (i=0;i<j;i++) {
				inf>>k>>ws;
				chi.push_back(k);
			}
		}
		//else if (b=="SMILES") {
		//	b1="";
        //    inf>>b1>>ws;
		//	molesmi=smiles=b1;
		//}

	}
	inf.close();
	return 1;
}

unsigned int MOLECULE::cyclization(unsigned int pt1,unsigned int pt2,unsigned int pbnd) {
	if (pt1==pt2) return 0;
	if (para.protect) {
		if (protect.at(pt1)) return 0;
		if (protect.at(pt2)) return 0;
	}

	if (1) {
	    vector<int> C_ringmember(0);
		
	    int rep_quota=1;
	    
		unsigned int compa=Cindex.at(pt1);
	    unsigned int cont=0;
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
	        for (unsigned int k2=0;k2<C_ringmember.size()-1;k2++) {
	            for (unsigned int k3=k2+1;k3<C_ringmember.size();k3++) {
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
		unsigned int num_rmember=0;
		if (C_ringmember.size()>0) {
	    	for (unsigned int k2=0;k2<C_ringmember.size();k2++) {
	        	if (C_ringmember.at(k2)!=-1) num_rmember++;
	    	}
		}
		vector<int>().swap(C_ringmember);
		if (num_rmember<5) return 0;		
	}
	
	// Ensure the bonds of atom P(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
	if (1) {
        //bool bndchk=0;

	    unsigned int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
		vector<unsigned int> ordcount(6,0);
		if (Rindex.at(pt1)>=1) {
			bndsum+=Rindex.at(pt1);
			ordcount.at(Rindex.at(pt1)-1)+=1;
		}

		if (Cyindex.at(pt1).size()) {
			for (unsigned int k1=0;k1<Cyindex.at(pt1).size();k1++) {
				bndsum+=Cybnd.at(Cyindex.at(pt1).at(k1)-1);
				//cout << "Cybnd.at(Cyindex.at(pt1).at(k1)-1) " << Cybnd.at(Cyindex.at(pt1).at(k1)-1) << endl;
				ordcount.at(Cybnd.at(Cyindex.at(pt1).at(k1)-1)-1)+=1;
			}
		}

    	for (unsigned int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(pt1)) {
				bndsum+=Rindex.at(k1);
				ordcount.at(Rindex.at(k1)-1)+=1;
			}
		}
        bndsum+=pbnd;
        ordcount.at(pbnd-1)+=1;

        unsigned int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(pt1)].order.at(k2);
        }

        for (unsigned int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[Mindex.at(pt1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=0;
				return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(pt1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

	}
	

	if (1) {
        //bool bndchk=0;

	    unsigned int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
		vector<unsigned int> ordcount(6,0);
		if (Rindex.at(pt2)>=1) {
			bndsum+=Rindex.at(pt2);
			ordcount.at(Rindex.at(pt2)-1)+=1;
		}
		if (Cyindex.at(pt2).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(pt2).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(pt2).at(k1)-1);
				//cout << "Cybnd.at(Cyindex.at(pt2).at(k1)-1) " << Cybnd.at(Cyindex.at(pt2).at(k1)-1) << endl;
                ordcount.at(Cybnd.at(Cyindex.at(pt2).at(k1)-1)-1)++;
            }
		}
    	for (unsigned int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(pt2)) {
				bndsum+=Rindex.at(k1);
				ordcount.at(Rindex.at(k1)-1)+=1;
			}
		}
        bndsum+=pbnd;
        ordcount.at(pbnd-1)+=1;


        unsigned int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(pt2)].order.at(k2);
        }

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[Mindex.at(pt2)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=0;
				return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(pt2)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }
	}
	

	if_circle++;
	Cyindex.at(pt1).push_back(if_circle);
	Cyindex.at(pt2).push_back(if_circle);
	Cybnd.resize(if_circle,0);
	Cybnd.at(if_circle-1)=pbnd;

    del_unpaired_ring_no();
    //decyc_small_ring(5);

	//chk_cistrans(); //20200806
	if (!chk_imine_ct(0,Cindex.size()-1)) mds2smi();
	//chk_imine_ct(0,Cindex.size()-1);
	//mds2smi();

	return 1;
}


unsigned int MOLECULE::mds2smi(unsigned int chirality,bool cistrans,bool ct_on) {
	//if (!ct_on) chk_cistrans(0,Cindex.size()-1,cistrans);
	//chk_imine_ct(0,Cindex.size()-1,cistrans);
	chk_cistrans(0,Cindex.size()-1,cistrans);
	chk_chirality(0,Cindex.size()-1,chirality);
	reset();

	unsigned int totalleng=0;
	for (unsigned int i=0;i<atomsmi.size();i++) totalleng+=atomsmi.at(i).size();
	//for (unsigned int i=0;i<atomsmi.size();i++) totalleng+=data->a.at(Mindex.at(i)).name.length();
	
	if (if_circle) totalleng+=2*if_circle*(unsigned int)(log10(if_circle)+4);
	for (unsigned int i=0;i<Cybnd.size();i++) {
		if (Cybnd.at(i)==2 || Cybnd.at(i)==3) totalleng+=3;
	}
    for (unsigned int i=0;i<chi.size();i++) totalleng+=chi.at(i);
	for (unsigned int j=0;j<2;j++) {
    	for (unsigned int i=0;i<ctsisomer.at(j).size();i++) {
    		if (ctsisomer.at(j).at(i)!="") {
    			string name=ctsisomer.at(j).at(i);
    			totalleng+=name.length();
    		}
    	}
	}
	totalleng+=10;
	char *x=new char [totalleng];
	for (unsigned int i=0;i<totalleng;i++) x[i]=' ';

	vector<bool> no_bnd_and_suffix(Cindex.size(),0);
	for (unsigned int i=0;i<Cindex.size();i++) {
		if (atomsmi.at(i).size()>data->a.at(Mindex.at(i)).index-1) no_bnd_and_suffix.at(i)=0;
		else no_bnd_and_suffix.at(i)=1;
	}
    vector<unsigned int> m(Cindex.size(),0),chi1(Cindex.size(),0);
    vector<int> pos(atomsmi.size(),0);
    vector<string> ctf(ctsisomer.at(0).size(),""),cte(ctsisomer.at(1).size(),"");
    vector< vector<unsigned int> > cynum(Cyindex.size(),vector<unsigned int>(0));
    vector< vector<unsigned int> > cyend(Cyindex.size(),vector<unsigned int>(0));
    vector<bool> iscyend(Cybnd.size(),0);
    for (unsigned int i=0;i<Cindex.size();i++) {
        m.at(i)=i;
        pos.at(i)=atomsmi.at(i).at(0);
		ctf.at(i)=ctsisomer.at(0).at(i);
		cte.at(i)=ctsisomer.at(1).at(i);
		chi1.at(i)=chi.at(i);
		for (unsigned int j=0;j<Cyindex.at(i).size();j++) {
			cynum.at(i).push_back(Cyindex.at(i).at(j));
            if (!iscyend.at(Cyindex.at(i).at(j)-1)) iscyend.at(Cyindex.at(i).at(j)-1)=1;
            else cyend.at(i).push_back(Cyindex.at(i).at(j));
		}
    }

    pos.reserve(pos.size());
    m.reserve(m.size());
    ctf.reserve(ctf.size());
	cte.reserve(ctf.size());
	chi1.reserve(chi1.size());
	cynum.reserve(cynum.size());
	
	for (unsigned int i=0;i<pos.size();i++) {
        for (unsigned int j=i+1;j<pos.size();j++) {
            if (pos.at(i)>pos.at(j)) {
                swap(pos.at(i),pos.at(j));
                swap(m.at(i),m.at(j));
                swap(ctf.at(i),ctf.at(j));
				swap(cte.at(i),cte.at(j));
				swap(chi1.at(i),chi1.at(j));
				swap(cynum.at(i),cynum.at(j));
            }
        }
    }

    for (unsigned int i=0;i<m.size();i++) {
		if (ctf.at(i)!="") {
        	int posi=0;

        	posi=atomsmi.at(m.at(i)).at(0); //at
        	for (int k=0;k<(int)ctf.at(i).length();k++) {
            	atomsmi.at(m.at(i)).insert(atomsmi.at(m.at(i)).begin(),-9);
        	}

        	for (unsigned int j=0;j<Cindex.size();j++) {
            	for (unsigned int k=0;k<atomsmi.at(j).size();k++) {
                	if (atomsmi.at(j).at(k)>=posi) atomsmi.at(j).at(k)+=ctf.at(i).length();
            	}
        	}
            for (unsigned int j=0;j<atomsmi.at(m.at(i)).size();j++) {
                if (atomsmi.at(m.at(i)).at(j)==-9) {
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
                    for (int k=totalleng-1;(int)(k-ctf.at(i).length())>=atomsmi.at(m.at(i)).at(j);k--) {
                        x[k]=x[k-ctf.at(i).length()];
                        x[k-ctf.at(i).length()]=' ';
                    }
                    for (int k=0;k<(int)ctf.at(i).length();k++) {
                        x[atomsmi.at(m.at(i)).at(j)+k]=ctf.at(i).at(k);
                    }

                    break;
                }
            }

        }
        if (chi1.at(i)!=0) {
            int posi=0;

            posi=atomsmi.at(m.at(i)).at(data->a.at(Mindex.at(m.at(i))).index-2+ctf.at(i).length());

            for (int k=0;k<(int)chi1.at(i);k++) atomsmi.at(m.at(i)).insert(atomsmi.at(m.at(i)).begin()+data->a.at(Mindex.at(m.at(i))).index-2+ctf.at(i).length(),-9);

            for (unsigned int j=0;j<Cindex.size();j++) {
                for (unsigned int k=0;k<atomsmi.at(j).size();k++) {
                    if (atomsmi.at(j).at(k)>=posi) atomsmi.at(j).at(k)+=chi1.at(i);
                }
            }

            for (unsigned int j=0;j<atomsmi.at(m.at(i)).size();j++) {
                if (atomsmi.at(m.at(i)).at(j)==-9) {
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
                    for (int k=totalleng-1;(int)(k-chi1.at(i))>=atomsmi.at(m.at(i)).at(j);k--) {
                        x[k]=x[k-chi1.at(i)];
                        x[k-chi1.at(i)]=' ';
                    }
                    for (int k=0;k<(int)chi1.at(i);k++) {
                        x[atomsmi.at(m.at(i)).at(j)+k]='@';
                    }

                    break;
                }
            }

        }
        if (cynum.at(i).size() || cte.at(i)!="") {
            int posi=0;

            int digits=0;
            for (unsigned int k=0;k<cynum.at(i).size();k++) {
                if (0) { //obabel
                    if (cynum.at(i).at(k)<10) digits+=2+(int)(log10(cynum.at(i).at(k))+1);
                    else if (cynum.at(i).at(k)>=10) digits+=1+(int)(log10(cynum.at(i).at(k))+1);
                }
                if (1) { //rdkit
                    digits+=4+(int)(log10(cynum.at(i).at(k)));
                }
            }
            if (cyend.at(m.at(i)).size()) {
                for (unsigned int j=0;j<cyend.at(m.at(i)).size();j++) {
                    if (Cybnd.at(cyend.at(m.at(i)).at(j)-1)==2 || Cybnd.at(cyend.at(m.at(i)).at(j)-1)==3) {
						if (1) digits+=3; //obabel
						if (0) digits+=1; //rdkit
					}
                }
            }

            stringstream nn("");

            for (unsigned int k=0;k<cynum.at(i).size();k++) {
                for (unsigned int k1=0;k1<cyend.at(m.at(i)).size();k1++) {
                    if (cynum.at(i).at(k)==cyend.at(m.at(i)).at(k1)) {
                        if (1) { //obabel
                            if (Cybnd.at(cynum.at(i).at(k)-1)==2) nn << "(=)";
                            else if (Cybnd.at(cynum.at(i).at(k)-1)==3) nn << "(#)";
                        }
                        if (0) { //rdkit
                            if (Cybnd.at(cynum.at(i).at(k)-1)==2) nn << "=";
                            else if (Cybnd.at(cynum.at(i).at(k)-1)==3) nn << "#";
                        }

                        break;
                    }
                }
                if (0) { //obabel
                    if (cynum.at(i).at(k)<10) nn << "%" << cynum.at(i).at(k) << "%";
                    else if (cynum.at(i).at(k)>=10) nn << "%" << cynum.at(i).at(k);
                }
                if (1) { //rdkit
                    nn << "%(" << cynum.at(i).at(k) << ")";
                }
            }

			nn << cte.at(i);

			//cout << "nn.str(): " << nn.str() << endl;

            if (!no_bnd_and_suffix.at(m.at(i))) posi=atomsmi.at(m.at(i)).at(data->a.at(Mindex.at(m.at(i))).index-1+ctf.at(i).length()+chi1.at(i));
            else posi=atomsmi.at(m.at(i)).at(data->a.at(Mindex.at(m.at(i))).index-2+ctf.at(i).size()+chi1.at(i));

            for (int k=0;k<(int)(digits+cte.at(i).length());k++) {
                if (!no_bnd_and_suffix.at(m.at(i))) atomsmi.at(m.at(i)).insert(atomsmi.at(m.at(i)).begin()+data->a.at(Mindex.at(m.at(i))).index-1+ctf.at(i).length()+chi1.at(i),-9);
                else atomsmi.at(m.at(i)).push_back(-9);  //20200910
            }


            for (unsigned int j=0;j<Cindex.size();j++) {
                for (unsigned int k=0;k<atomsmi.at(j).size();k++) {
                    if (!no_bnd_and_suffix.at(m.at(i))) {
                        if (atomsmi.at(j).at(k)>=posi) atomsmi.at(j).at(k)+=(digits+cte.at(i).length());
                    }
                    else { //20200910
                        if (atomsmi.at(j).at(k)>posi) atomsmi.at(j).at(k)+=(digits+cte.at(i).length());
                    }
                }
            }
            for (unsigned int j=0;j<atomsmi.at(m.at(i)).size();j++) {
                if (atomsmi.at(m.at(i)).at(j)==-9) {
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
                    for (int k=totalleng-1;(int)(k-(digits+cte.at(i).length()))>=(int)atomsmi.at(m.at(i)).at(j);k--) {
                        x[k]=x[k-(digits+cte.at(i).length())];
                        x[k-(digits+cte.at(i).length())]=' ';
                    }
                    for (int k=0;k<(int)(digits+cte.at(i).length());k++) {
                        x[atomsmi.at(m.at(i)).at(j)+k]=nn.str().at(k);
                    }

                    break;
                }
            }

        }

        //for (unsigned int j=0;j<totalleng;j++) cout << x[j];
        //cout << endl;

	}

    pos.clear();
    m.clear();
    ctf.clear();
	cte.clear();
	chi1.clear();
	cynum.clear();
	cyend.clear();
	iscyend.clear();
    vector<int>().swap(pos);
    vector<unsigned int>().swap(m);
    vector<string>().swap(ctf);
	vector<string>().swap(cte);
	vector< vector<unsigned int> >().swap(cynum);
	vector< vector<unsigned int> >().swap(cyend);
	vector<bool>().swap(iscyend);

	chg=0;
	for (unsigned int i=0;i<Cindex.size();i++) chg+=data->a.at(Mindex.at(i)).chg;

	//int implicit_ring_n=0;
	for (unsigned int i=0;i<Cindex.size();i++) {
		unsigned int M = Mindex.at(i);
		//int R = Rindex.at(i);
		string kk=data->a.at(M).name; //20191011

		bool explicit_H=0;
		//if (data->a.at(M).chg && kk.length()<=data->a.at(M).nbond && data->a.at(M).index!=data->a.at(M).nbond) explicit_H=1; //kk.length()<=data->a.at(M).nbond
		if (data->a.at(M).index>=6 && data->a.at(M).index!=data->a.at(M).nbond) explicit_H=1;
		//if (data->a.at(M).atm=="P") explicit_H=1;
		unsigned int Hposi=0;
		if (explicit_H) {
			for (int g=0;g<(int)data->a.at(M).index-1;g++) {
				if ((int)kk[g]>=48 && (int)kk[g]<=57) {
					Hposi=g;
					break;
				}
				else if ((int)kk[g]<48 || (int)kk[g]>57) {
					if (g>=(int)data->a.at(M).index-2) {
						explicit_H=0;
						break;
					}
				}
			}
		}
		//if (M==67 || M==68) explicit_H=1;

		if (explicit_H) { // data.a[M].chg && kk.length()<=data.a[M].nbond && data.a[M].name[0]=='['
			int tmp=0;
			for (unsigned int j=0;j<data->a.at(M).norder;j++) {
				tmp+=Bindex.at(i).at(j);
				//if (Bindex.at(i).at(j)) tmp+=1;
			}
			if (1) {
				if (tmp>0 && tmp<10) kk[Hposi]=char(tmp+48); //data->a[M].name[3]=char(tmp+48);
			}
		}

        int num_3blanks=0;
        if (1 && data->a.at(M).norder>0) { //delete "(-)" , "(=)" , "(#)"
			if (1) {
	            for (unsigned int j=0;j<kk.size();j++) {
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
				for (unsigned int j=0;j<Bindex.at(i).size();j++) {
					if (Bindex.at(i).at(j)) {
                		for (unsigned int k=0;k<kk.size();k++) {
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
                for (int j=0;j<(int)kk.length();j++) {
                    if (kk[j]==' ') {
                        kk.erase(kk.begin()+j);
                        j--;
                    }
                }
			}
        }
		
        if (1 && Cyindex.at(i).size()) { //delete "(-)", "(=)", and "(#)" of cyclic flag
			for (unsigned int j=0;j<Cyindex.at(i).size();j++) {
				for (unsigned int k=0;k<kk.size();k++) {
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

            	for (int j=0;j<(int)kk.length();j++) {
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
				for (unsigned int j=0;j<kk.size();j++) {
					if (j>0 && j<kk.size()-1) {
						if (kk[j-1]=='(' && kk[j]=='-' && kk[j+1]==')') {
							kk[j]=' ';
							num_1blank++;
						}
					}
				}
	            for (int j=0;j<(int)kk.length();j++) {
    	            if (kk[j]==' ') {
        	            kk.erase(kk.begin()+j);
            	        j--;
                	}
            	}
			}
		}
		
		int ct=0;
		for (unsigned int j=0;j<atomsmi.at(i).size();j++) { // orig data->a[M].name.length() // data.a[M].nbond  //w<j  //data->a[M].name.length()
			int n=atomsmi.at(i).at(j);

			if (x[n]==' ') {
				x[n]=kk[ct];
				ct++;
			}

		}

		//for (unsigned int j=0;j<totalleng;j++) cout << x[j];
		//cout << endl;
		
	}

    atomsmi.resize(0,vector<int>(0));
    //atomsmi.shrink_to_fit();
    Bindex.resize(0,vector<unsigned int>(0));
    //Bindex.shrink_to_fit();

	//molesmi="";
	molesmi.clear();
	for (unsigned int i=0;i<totalleng;i++) {
		if (x[i]!=' ') {
			molesmi.push_back(x[i]);
		}
	}
    delete [] x;
    x=NULL;

	smiles=molesmi;
	canonicalize_SMILES();

	//cout << "SMI: " << molesmi << endl;

    //if (cts) {
    //    ofstream chkct("./chkct",ios::app);
    //    chkct << molesmi << endl;
    //    chkct.close();
    //}

	return 1;
}


unsigned int MOLECULE::clear() {
	//int n;
	molesmi.clear();
	//Bindex.clear();
	//atomsmi.clear();
	Bindex.resize(0,vector<unsigned int> (0));
	atomsmi.resize(0,vector<int> (0));
	//ctsisomer.clear();

	return 1;
}

unsigned int MOLECULE::chk_chirality(unsigned int sposi,unsigned int lposi,unsigned int chirality) {
	for (unsigned int i=sposi;i<=lposi;i++) { // C\C=C\C  //i=lposi;i>=sposi;i-- //i=sposi;i<=lposi;i++
		bool go=0;
		if (data->a.at(Mindex.at(i)).chigenic) go=1;
		if (go) {
			vector< vector<unsigned int> > chainsM(6,vector<unsigned int>(0));
			vector< vector<unsigned int> > chainsR(6,vector<unsigned int>(0));
			vector< vector<unsigned int> > branchpos(6,vector<unsigned int>(0));
			int P=Pindex.at(i);
			int C=Cindex.at(i);
			unsigned int ct=0;

			if (P>0) {
				unsigned int compa=Cindex.at(P-1),bnd1=Rindex.at(C-1);
				unsigned int cont=0;
				while (compa>=1 && cont<=Cindex.size()) {
					chainsM.at(0).push_back(Mindex.at(compa-1));
					chainsR.at(0).push_back(bnd1); //Rindex.at(compa-1)

					vector<unsigned int> curatm(0);
					curatm.push_back(Cindex.at(compa-1));

					bool goout=0;

					do {
						vector<unsigned int> tmpp(0);
						for (unsigned int k3=0;k3<curatm.size();k3++) {
							for (unsigned int k2=0;k2<Cindex.size();k2++) {
								if (Pindex.at(k2)==curatm.at(k3) && Cindex.at(k2)!=Cindex.at(P-1)) {
									chainsM.at(0).push_back(Mindex.at(k2));
									chainsR.at(0).push_back(Rindex.at(k2));
									branchpos.at(0).push_back(chainsM.at(0).size()-1);
									tmpp.push_back(Cindex.at(k2));
								}
							}
						}
						if (tmpp.size()>0) {
							curatm.resize(0);
							for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
						}
						else {
							goout=1;
							break;
						}
					} while (!goout);

					bnd1=Rindex.at(compa-1);

					if (compa<=1) break;

					compa=Cindex.at(Pindex.at(compa-1)-1);

					cont++;
				}
				ct++;
			}
		
			for (unsigned int k1=0;k1<Cindex.size();k1++) {
				if (Pindex.at(k1)==Cindex.at(i)) {
					chainsM.at(ct).push_back(Mindex.at(k1));
					chainsR.at(ct).push_back(Rindex.at(k1));

					vector<unsigned int> curatm(0);
					curatm.push_back(Cindex.at(k1));
								
					bool goout=0;
								
					do {
						vector<unsigned int> tmpp(0);
						for (unsigned int k3=0;k3<curatm.size();k3++) {
							for (unsigned int k2=0;k2<Cindex.size();k2++) {
								if (Pindex.at(k2)==curatm.at(k3)) {
									chainsM.at(ct).push_back(Mindex.at(k2));
									chainsR.at(ct).push_back(Rindex.at(k2));
									branchpos.at(ct).push_back(chainsM.at(ct).size()-1);
									tmpp.push_back(Cindex.at(k2));
								}
							}
						}
						if (tmpp.size()>0) {
							curatm.resize(0);
							for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
						}
						else {
							goout=1;
							break;
						}
					} while (!goout);

					ct++;
				}
			}

			if (chainsM.at(4).size() || chainsM.at(5).size()) {
	            vector< vector<unsigned int> >().swap(chainsM);
    	        vector< vector<unsigned int> >().swap(chainsR);
        	    vector< vector<unsigned int> >().swap(branchpos);
				return 0;
			}

			for (unsigned int g=0;g<chainsM.size();g++) {
				for (unsigned int g1=0;g1<chainsM.at(g).size();g1++) {
					for (unsigned int g2=g1+1;g2<chainsM.at(g).size();g2++) {
						if (chainsM.at(g).at(g1)>chainsM.at(g).at(g2)) {
							swap(chainsM.at(g).at(g1),chainsM.at(g).at(g2));
							swap(chainsR.at(g).at(g1),chainsR.at(g).at(g2));
						}
						else if (chainsM.at(g).at(g1)==chainsM.at(g).at(g2) && chainsR.at(g).at(g1)<chainsR.at(g).at(g2)) {
							swap(chainsM.at(g).at(g1),chainsM.at(g).at(g2));
							swap(chainsR.at(g).at(g1),chainsR.at(g).at(g2));
						}
					}
				}
			}

			bool haschirality=1;
			for (unsigned int g1=0;g1<4;g1++) {
				for (unsigned int g2=g1+1;g2<4;g2++) {
					if (chainsM.at(g1).size()==chainsM.at(g2).size()) {
						if (chainsM.at(g1).size()==0) {
							haschirality=0;
							break;
						}
						else {
							for (unsigned int g=0;g<chainsM.at(g1).size();g++) {
								if (chainsM.at(g2).at(g)-chainsM.at(g1).at(g)) {
									break;
								}
								else if (chainsM.at(g2).at(g)==chainsM.at(g1).at(g) && g>=chainsM.at(g1).size()-1) {
									haschirality=0;
								}
							}
						}
					}
					else haschirality=1;

					if (!haschirality) break;
				}
				if (!haschirality) break;
			}

			if (!haschirality) {
				for (unsigned int g1=0;g1<4;g1++) {
					for (unsigned int g2=g1+1;g2<4;g2++) {
						if (chainsR.at(g1).size()==chainsR.at(g2).size()) {
							if (chainsR.at(g1).size()==0) {
								haschirality=0;
								break;
							}
							else {
								for (unsigned int g=0;g<chainsR.at(g1).size();g++) {
									if (chainsR.at(g2).at(g)-chainsR.at(g1).at(g)) {
										haschirality=1;
										break;
									}
									else if (chainsR.at(g2).at(g)==chainsR.at(g1).at(g) && g>=chainsR.at(g1).size()-1) {
										haschirality=0;
									}
								}
							}
						}
						else haschirality=1;

						if (!haschirality) break;
					}
					if (!haschirality) break;
				}
			}

			if (!haschirality) {
				for (unsigned int g1=0;g1<4;g1++) {
					for (unsigned int g2=g1+1;g2<4;g2++) {
						if (branchpos.at(g1).size()==branchpos.at(g2).size()) {
							if (branchpos.at(g1).size()==0) {
								haschirality=0;
								break;
							}
							else {
								for (unsigned int g=0;g<branchpos.at(g1).size();g++) {
									if (branchpos.at(g2).at(g)-branchpos.at(g1).at(g)) {
										haschirality=1;
										break;
									}
									else if (branchpos.at(g2).at(g)==branchpos.at(g1).at(g) && g>=branchpos.at(g1).size()-1) {
										haschirality=0;
									}
								}
							}
						}
						else haschirality=1;

						if (!haschirality) break;
					}
					if (!haschirality) break;
				}
			}


			if (haschirality && !chi.at(i)) chi.at(i)=chirality;

			vector< vector<unsigned int> >().swap(chainsM);
			vector< vector<unsigned int> >().swap(chainsR);
			vector< vector<unsigned int> >().swap(branchpos);

		}
		//if (data->a.at(Mindex.at(i)).norder<4 && data->a.at(Mindex.at(i)).atm!="Group" && chi.at(i)) chi.at(i)=0;
		if (!data->a.at(Mindex.at(i)).chigenic && chi.at(i)) chi.at(i)=0;
	}


	return 1;

}

unsigned int MOLECULE::chk_cistrans(unsigned int sposi,unsigned int lposi,bool cistrans) {

	if (0) {
		if (ctsisomer.size()==2) {
			if (ctsisomer.at(0).size()==ctsisomer.at(1).size() && ctsisomer.at(0).size()==Cindex.size()) return 0;
		}
	}

    if (0) {
		ctsisomer.resize(0,vector<string> (0));
		ctsisomer.resize(2,vector<string>(Cindex.size(),""));
    }

    if (1) {
		vector< vector<bool> > shouldflag(2,vector<bool>(Cindex.size(),0));
        for (unsigned int i=sposi;i<=lposi;i++) { // C\C=C\C  //i=lposi;i>=sposi;i-- //i=sposi;i<=lposi;i++

			bool go=0;
			if (data->a.at(Mindex.at(i)).ctgenic) go=1;
			//bool go=1;

            if (Rindex.at(i)==2 && go) {
                int P=Pindex.at(i);
				int C=Cindex.at(i);

                if (P>0) {
					bool go1=0;
					if (data->a.at(Mindex.at(P-1)).ctgenic) go1=1;
					//bool go1=1;

                    if (go1 && Cindex.at(P-1)>0) { // Mindex.at(P-1)==2 && Pindex.at(P-1)>0
						vector< vector<unsigned int> > chainsM(6,vector<unsigned int>(0));
                        vector< vector<unsigned int> > chainsR(6,vector<unsigned int>(0));
						vector< vector<unsigned int> > branchpos(6,vector<unsigned int>(0));

						bool cas=0;
						if (Pindex.at(P-1)>0) {
							cas=0;
						    if (Cindex.at(Pindex.at(P-1)-1)>0) {
        						unsigned int compa=Cindex.at(Pindex.at(P-1)-1),bnd1=Rindex.at(P-1);
        						unsigned int cont=0;
        						while (compa>=1 && cont<=Cindex.size()) {
									chainsM.at(0).push_back(Mindex.at(compa-1));
									chainsR.at(0).push_back(bnd1); //Rindex.at(compa-1)

									vector<unsigned int> curatm(0);
									curatm.push_back(Cindex.at(compa-1));

									bool goout=0;

									do {
										vector<unsigned int> tmpp(0);
										for (unsigned int k3=0;k3<curatm.size();k3++) {
											for (unsigned int k2=0;k2<Cindex.size();k2++) {
												if (Pindex.at(k2)==curatm.at(k3) && Cindex.at(k2)!=Cindex.at(P-1)) {
													chainsM.at(0).push_back(Mindex.at(k2));
													chainsR.at(0).push_back(Rindex.at(k2));
													branchpos.at(0).push_back(chainsM.at(0).size()-1);
													tmpp.push_back(Cindex.at(k2));
												}
											}
										}
										if (tmpp.size()>0) {
											curatm.resize(0);
											for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
										}
										else {
											goout=1;
											break;
										}
									} while (!goout);

									bnd1=Rindex.at(compa-1);

									if (compa<=1) break;

									compa=Cindex.at(Pindex.at(compa-1)-1);

									cont++;
        						}
    						}

							int ct=1;
                            for (unsigned int k1=0;k1<Cindex.size();k1++) {
                                if (Pindex.at(k1)==Cindex.at(P-1) && Cindex.at(k1)!=Cindex.at(i)) {
                                    chainsM.at(ct).push_back(Mindex.at(k1));
                                    chainsR.at(ct).push_back(Rindex.at(k1));

                                    vector<unsigned int> curatm(0);
                                    curatm.push_back(Cindex.at(k1));

                                    bool goout=0;

                                    do {
										vector<unsigned int> tmpp(0);
                                        for (unsigned int k3=0;k3<curatm.size();k3++) {
                                            for (unsigned int k2=0;k2<Cindex.size();k2++) {
                                                if (Pindex.at(k2)==curatm.at(k3)) {
                                                    chainsM.at(ct).push_back(Mindex.at(k2));
                                                    chainsR.at(ct).push_back(Rindex.at(k2));
													branchpos.at(ct).push_back(chainsM.at(ct).size()-1);
                                                    tmpp.push_back(Cindex.at(k2));
                                                }
                                            }
										}
                                        if (tmpp.size()>0) {
                                            curatm.resize(0);
                                            for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
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
						else {
							cas=1;
							int ct=0;
                            for (unsigned int k1=0;k1<Cindex.size();k1++) {
                                if (Pindex.at(k1)==Cindex.at(P-1) && Cindex.at(k1)!=Cindex.at(i)) {
                                    chainsM.at(ct).push_back(Mindex.at(k1));
                                    chainsR.at(ct).push_back(Rindex.at(k1));

                                    vector<unsigned int> curatm(0);
                                    curatm.push_back(Cindex.at(k1));

                                    bool goout=0;

                                    do {
										vector<unsigned int> tmpp(0);
                                        for (unsigned int k3=0;k3<curatm.size();k3++) {
                                            for (unsigned int k2=0;k2<Cindex.size();k2++) {
                                                if (Pindex.at(k2)==curatm.at(k3)) {
                                                    chainsM.at(ct).push_back(Mindex.at(k2));
                                                    chainsR.at(ct).push_back(Rindex.at(k2));
													branchpos.at(ct).push_back(chainsM.at(ct).size()-1);
                                                    tmpp.push_back(Cindex.at(k2));
                                                }
                                            }
										}
                                        if (tmpp.size()>0) {
                                            curatm.resize(0);
                                            for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
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

						if (chainsM.at(2).size()) {
	                        vector< vector<unsigned int> >().swap(chainsM);
    	                    vector< vector<unsigned int> >().swap(chainsR);
        	                vector< vector<unsigned int> >().swap(branchpos);
							return 0;
						}

						int ct=3;
                        for (unsigned int k1=0;k1<Cindex.size();k1++) {
                            if (Pindex.at(k1)==Cindex.at(i)) {
                                chainsM.at(ct).push_back(Mindex.at(k1));
								chainsR.at(ct).push_back(Rindex.at(k1));

								vector<unsigned int> curatm(0);
								curatm.push_back(Cindex.at(k1));
									
								bool goout=0;
									
								do {
									vector<unsigned int> tmpp(0);
									for (unsigned int k3=0;k3<curatm.size();k3++) {
										for (unsigned int k2=0;k2<Cindex.size();k2++) {
                                			if (Pindex.at(k2)==curatm.at(k3)) {
                                    			chainsM.at(ct).push_back(Mindex.at(k2));
                                    			chainsR.at(ct).push_back(Rindex.at(k2));
												branchpos.at(ct).push_back(chainsM.at(ct).size()-1);
												tmpp.push_back(Cindex.at(k2));
                                			}
										}
									}
									if (tmpp.size()>0) {
										curatm.resize(0);
										for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
									}
									else {
										goout=1;
										break;
									}
								} while (!goout);

								ct++;
                            }
                        }

                        if (chainsM.at(5).size()) {
                            vector< vector<unsigned int> >().swap(chainsM);
                            vector< vector<unsigned int> >().swap(chainsR);
                            vector< vector<unsigned int> >().swap(branchpos);
                            return 0;
                        }

						for (unsigned int g=0;g<chainsM.size();g++) {
							for (unsigned int g1=0;g1<chainsM.at(g).size();g1++) {
								for (unsigned int g2=g1+1;g2<chainsM.at(g).size();g2++) {
									if (chainsM.at(g).at(g1)>chainsM.at(g).at(g2)) {
										swap(chainsM.at(g).at(g1),chainsM.at(g).at(g2));
										swap(chainsR.at(g).at(g1),chainsR.at(g).at(g2));
									}
									else if (chainsM.at(g).at(g1)==chainsM.at(g).at(g2) && chainsR.at(g).at(g1)<chainsR.at(g).at(g2)) {
                                        swap(chainsM.at(g).at(g1),chainsM.at(g).at(g2));
                                        swap(chainsR.at(g).at(g1),chainsR.at(g).at(g2));
									}
								}
							}
						}

						if (0) {
                            for (unsigned int g=0;g<chainsM.size();g++) {
								//cout << "F " << filenum << " | chg " << chg << " | " << g+1 << " | ";
                                for (unsigned int g1=0;g1<chainsM.at(g).size();g1++) {
                                    cout << chainsM.at(g).at(g1) << " ";
                                }
								cout << endl;
                            }
						}

						bool hascistrans=1;
						if (chainsM.at(0).size()==chainsM.at(1).size()) {
							if (chainsM.at(0).size()==0) hascistrans=0;
							else {
								for (unsigned int g=0;g<chainsM.at(0).size();g++) {
									if (chainsM.at(1).at(g)-chainsM.at(0).at(g)) {
										hascistrans=1;
										break;
									}
									else if (chainsM.at(1).at(g)==chainsM.at(0).at(g) && g>=chainsM.at(0).size()-1) {
										hascistrans=0;
									}
								}
							}
						}
						else hascistrans=1;

                        if (!hascistrans) {
                            if (chainsR.at(0).size()==chainsR.at(1).size()) {
                                if (chainsR.at(0).size()==0) hascistrans=0;
                                else {
                                    for (unsigned int g=0;g<chainsR.at(0).size();g++) {
                                        if (chainsR.at(1).at(g)-chainsR.at(0).at(g)) {
											hascistrans=1;
                                            break;
                                        }
                                        else if (chainsR.at(1).at(g)==chainsR.at(0).at(g) && g>=chainsR.at(0).size()-1) {
                                            hascistrans=0;
                                        }
                                    }
                                }
                            }
							else hascistrans=1;
                        }

                        if (!hascistrans) {
                            if (branchpos.at(0).size()==branchpos.at(1).size()) {
                                if (branchpos.at(0).size()==0) hascistrans=0;
                                else {
                                    for (unsigned int g=0;g<branchpos.at(0).size();g++) {
                                        if (branchpos.at(1).at(g)-branchpos.at(0).at(g)) {
                                            break;
                                        }
                                        else if (branchpos.at(1).at(g)==branchpos.at(0).at(g) && g>=branchpos.at(0).size()-1) {
                                            hascistrans=0;
                                        }
                                    }
                                }
                            }
							else hascistrans=1;
                        }


						if (hascistrans) {
                            if (chainsM.at(3).size()==chainsM.at(4).size()) {
                                if (chainsM.at(3).size()==0) hascistrans=0;
                                else {
                                    for (unsigned int g=0;g<chainsM.at(3).size();g++) {
                                        if (chainsM.at(4).at(g)-chainsM.at(3).at(g)) {
											hascistrans=1;
                                            break;
                                        }
										else if (chainsM.at(4).at(g)==chainsM.at(3).at(g) && g>=chainsM.at(3).size()-1) {
											hascistrans=0;
										}
                                    }
                                }
                            }
							else hascistrans=1;
						}

                        if (!hascistrans) {
                            if (chainsR.at(3).size()==chainsR.at(4).size()) {
                                if (chainsR.at(3).size()==0) hascistrans=0;
                                else {
                                    for (unsigned int g=0;g<chainsR.at(3).size();g++) {
                                        if (chainsR.at(4).at(g)-chainsR.at(3).at(g)) {
											hascistrans=1;
                                            break;
                                        }
                                        else if (chainsR.at(4).at(g)==chainsR.at(3).at(g) && g>=chainsR.at(3).size()-1) {
                                            hascistrans=0;
                                        }
                                    }
                                }
                            }
							else hascistrans=1;
                        }

                        if (!hascistrans) {
                            if (branchpos.at(3).size()==branchpos.at(4).size()) {
                                if (branchpos.at(3).size()==0) hascistrans=0;
                                else {
                                    for (unsigned int g=0;g<branchpos.at(3).size();g++) {
                                        if (branchpos.at(4).at(g)-branchpos.at(3).at(g)) {
                                            break;
                                        }
                                        else if (branchpos.at(4).at(g)==branchpos.at(3).at(g) && g>=branchpos.at(3).size()-1) {
                                            hascistrans=0;
                                        }
                                    }
                                }
                            }
							else hascistrans=1;
                        }



						bool y1=1,y2=1;
						if (hascistrans && if_circle) {
    						vector<int> C_ringmember(0);

    						for (int k1=if_circle;k1>=1;k1--) {
        						C_ringmember.resize(0);
        						int rep_quota=1;
        						for (int k2=Cyindex.size()-1;k2>=0;k2--) {
            						if (Cyindex.at(k2).size()) {
                						for (unsigned int k3=0;k3<Cyindex.at(k2).size();k3++) {
                    						if ((int)Cyindex.at(k2).at(k3)==k1) {
                        						unsigned int compa=Cindex.at(k2);
                        						unsigned int cont=0;
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
            						for (unsigned int k2=0;k2<C_ringmember.size()-1;k2++) {
                						for (unsigned int k3=k2+1;k3<C_ringmember.size();k3++) {
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
                					for (int k2=0;k2<(int)C_ringmember.size();k2++) {
                    					if (C_ringmember.at(k2)==-1) {
											C_ringmember.erase(C_ringmember.begin()+k2);
											k2--;
                    					}
                					}
        						}
								for (unsigned int g=0;g<C_ringmember.size();g++) {
									if ((int)Cindex.at(i)==C_ringmember.at(g)) y1=0;
                                    if ((int)Cindex.at(P-1)==C_ringmember.at(g)) y2=0;
								}

							}

						}

						bool z=0;
						if (y1 || y2) z=1;
						if (hascistrans && z) {
							if (!cas) {
								if (!cistrans) {
									//if (ctsisomer.at(0).at(P-1)=="" && ctsisomer.at(1).at(P-1)=="") ctsisomer.at(0).at(P-1)+="/"; //&& ctsisomer.at(1).at(P-1)==""
									//if (ctsisomer.at(1).at(C-1)=="" && ctsisomer.at(0).at(C-1)=="") ctsisomer.at(1).at(C-1)+="/"; //&& ctsisomer.at(0).at(C-1)==""
									if (ctsisomer.at(0).at(P-1)=="") ctsisomer.at(0).at(P-1)+="/";
									if (ctsisomer.at(1).at(C-1)=="") ctsisomer.at(1).at(C-1)+="/";
								}
								else {
									//if (ctsisomer.at(0).at(P-1)=="" && ctsisomer.at(1).at(P-1)=="") ctsisomer.at(0).at(P-1)+="/"; //&& ctsisomer.at(1).at(P-1)==""
									//if (ctsisomer.at(1).at(C-1)=="" && ctsisomer.at(0).at(C-1)=="") ctsisomer.at(1).at(C-1)+="\\"; //&& ctsisomer.at(0).at(C-1)==""
									if (ctsisomer.at(0).at(P-1)=="") ctsisomer.at(0).at(P-1)+="/";
									if (ctsisomer.at(1).at(C-1)=="") ctsisomer.at(1).at(C-1)+="\\";
								}
								shouldflag.at(0).at(P-1)=1;
								shouldflag.at(1).at(C-1)=1;
							}
							else {
							    if (Cindex.size()>2 && 0) {
        							if (ctsisomer.at(0).at(0)!="" && ctsisomer.at(1).at(1)!="" && Rindex.at(1)==2 && Pindex.at(1)==1) {
            							int ct1=0,ct2=0;
            							for (unsigned int k1=0;k1<Cindex.size();k1++) {
											if (Pindex.at(k1)==Cindex.at(0) && Cindex.at(k1)!=Cindex.at(1)) {
                                                if (ct1==0) {
                                                	if (ctsisomer.at(0).at(0)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
                                                	else if (ctsisomer.at(0).at(0)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
													//if (ctsisomer.at(0).at(0)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)+="/";
													//else if (ctsisomer.at(0).at(0)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)+="\\";

													shouldflag.at(0).at(Cindex.at(k1)-1)=1;
												}
												//if (ct1==0) {
												//	if (ctsisomer.at(0).at(0)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
												//	else if (ctsisomer.at(0).at(0)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
												//}
												//else if (ct1==1) {
												//  ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(0).at(0)+ctsisomer.at(0).at(Cindex.at(k1)-1);
												//}
												ct1++;
											}
											if (Pindex.at(k1)==Cindex.at(1)) {
												if (ct2==0) {
													ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(1).at(1)+ctsisomer.at(0).at(Cindex.at(k1)-1);
													
													shouldflag.at(0).at(Cindex.at(k1)-1)=1;
												}
												//if (ct2==0) {
												//    if (ctsisomer.at(1).at(1)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
												//    else if (ctsisomer.at(1).at(1)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
												//}
												//else if (ct2==1) {
												//  ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(1).at(1)+ctsisomer.at(0).at(Cindex.at(k1)-1);
                    							//}
                    							ct2++;
                							}
            							}
            							ctsisomer.at(0).at(0)="";
            							ctsisomer.at(1).at(1)="";
        							}
    							}

								if (1) {
									unsigned int ct1=0,ct2=0;
									for (unsigned int k1=0;k1<Cindex.size();k1++) {
										if (Pindex.at(k1)==Cindex.at(P-1) && Cindex.at(k1)!=Cindex.at(i)) {
											if (ct1==0) { // ctsisomer.at(0).at(k1)=="" &&
												if (!cistrans) ctsisomer.at(0).at(k1)+="/";
												else ctsisomer.at(0).at(k1)+="\\";
												shouldflag.at(0).at(k1)=1;
											}
											//if (ctsisomer.at(0).at(k1)=="" && ct1==0) {
											//	if (!cistrans) ctsisomer.at(0).at(k1)+="\\";
											//	else ctsisomer.at(0).at(k1)+="/";
											//}
											ct1++;
										}
										if (Pindex.at(k1)==Cindex.at(i)) {
											if (ct2==0) { //ctsisomer.at(0).at(k1)=="" && 
												ctsisomer.at(0).at(k1)+="/";
												
												shouldflag.at(0).at(k1)=1;
											}
											ct2++;
										}
									}
								}
							}
						}

						vector< vector<unsigned int> >().swap(chainsM);
						vector< vector<unsigned int> >().swap(chainsR);
						vector< vector<unsigned int> >().swap(branchpos);

					}
				}
			}
		}

		for (unsigned int i=0;i<2;i++) {
			for (unsigned int j=0;j<ctsisomer.at(i).size();j++) {
				if (!shouldflag.at(i).at(j) && ctsisomer.at(i).at(j)!="") ctsisomer.at(i).at(j)="";
				if (ctsisomer.at(i).at(j).length()>=2) {
					string null=ctsisomer.at(i).at(j);
					ctsisomer.at(i).at(j)=null.at(null.length()-1);
				}
			}
		}

		vector< vector<bool> >().swap(shouldflag);


		if (1) {
        	for (unsigned int i=0;i<Cindex.size();i++) {
            	if (Rindex.at(i)==2) {
                	int P=Pindex.at(i);
					if (P>0) {
						int PP=Pindex.at(P-1);
						if (PP>0) {
							if (ctsisomer.at(0).at(P-1)!="" && ctsisomer.at(1).at(PP-1)!="") {
								int PPP=Pindex.at(PP-1);
								if (PPP>0) {
									if (ctsisomer.at(1).at(PP-1)!=ctsisomer.at(0).at(P-1)) {
										//if (ctsisomer.at(1).at(PP-1)=="/") ctsisomer.at(1).at(PP-1)="\\";
										//else if (ctsisomer.at(1).at(PP-1)=="\\") ctsisomer.at(1).at(PP-1)="/";

                                    	if (ctsisomer.at(0).at(PPP-1)=="/") ctsisomer.at(0).at(PPP-1)="\\";
                                    	else if (ctsisomer.at(0).at(PPP-1)=="\\") ctsisomer.at(0).at(PPP-1)="/";

										//ctsisomer.at(1).at(PP-1)="";
										ctsisomer.at(1).at(PP-1)=ctsisomer.at(0).at(P-1);
									}
									//else ctsisomer.at(1).at(PP-1)="";
								}
							}
						}
					}
            	}
        	}
		}

	}

	return 1;

}


unsigned int MOLECULE::reset() {
	clear();
	Bindex.reserve(Cindex.size());
	Bindex.resize(Cindex.size(),vector<unsigned int> (0));
	//Bindex.shrink_to_fit();
	atomsmi.reserve(Cindex.size());
	atomsmi.resize(Cindex.size(),vector<int> (0));
	//atomsmi.shrink_to_fit();

	for (unsigned int i=0;i<Cindex.size();i++) {
		unsigned int b = Mindex.at(i);
		unsigned int c = Pindex.at(i);
		unsigned int x = Rindex.at(i);

		for (unsigned int j=0;j<data->a.at(b).norder;j++) {
			Bindex.at(i).push_back(data->a.at(b).order.at(j));
		}

		bool u=0;
		if (c==0) {
			for (unsigned int j=0;j<data->a.at(b).name.length();j++) { // data.a[b].nbond //20191230
				atomsmi.at(i).push_back(j);
			}
		}
		else if (c>0) {
			unsigned int M = Mindex.at(c-1);
			unsigned int j=0;
			for (j=0;j<data->a.at(M).norder;j++) {
				if (Bindex.at(c-1).at(j) == x) {
					for (unsigned int n=0;n<data->a.at(b).norder;n++) {
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
				for (unsigned int n=0;n<data->a.at(Mindex.at(j)).name.length();n++) { // data.a[Mindex.at(j)].nbond
					if (atomsmi.at(j).at(n) > atomsmi.at(c-1).at(smindex)) {
						atomsmi.at(j).at(n) = atomsmi.at(j).at(n)+data->a.at(b).name.length();  // + data.a[b].nbond
					}  
				}
			}

		}

	}


    if (1) { //delete "(-)" , "(=)" , "(#)" used for parent atom, ring, H atoms.
        if (0) {
            for (unsigned int i=0;i<atomsmi.size();i++) {
                cout << "orig ";
                for (unsigned int j=0;j<atomsmi.at(i).size();j++) {
                    cout << atomsmi.at(i).at(j) << " ";
                }
                cout << endl;
            }
        }
        vector<unsigned int> cid(0);
        vector<int> rpos(0);
        vector<int> apos(0);
        for (unsigned int i=0;i<Bindex.size();i++) {
            unsigned int M=Mindex.at(i);
            unsigned int R=Rindex.at(i);
            int P_quota=1;
            for (unsigned int j=0;j<data->a.at(M).norder;j++) {
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
		
        for (unsigned int i=0;i<cid.size();i++) {
            for (unsigned int j=i+1;j<cid.size();j++) {
                if (apos.at(j)<apos.at(i)) {
                    swap(cid.at(i),cid.at(j));
                    swap(apos.at(i),apos.at(j));
                    swap(rpos.at(i),rpos.at(j));
                }
            }
        }
        for (unsigned int i=0;i<cid.size();i++) {
            //int rtmp=rpos.at(i);
            int atmp=atomsmi.at(cid.at(i)-1).at(rpos.at(i));
            for (int j=0;j<3;j++) {
                atomsmi.at(cid.at(i)-1).erase(atomsmi.at(cid.at(i)-1).begin()+rpos.at(i));
            }
            for (unsigned int k=0;k<atomsmi.size();k++) {
                for (unsigned int l=0;l<atomsmi.at(k).size();l++) {
                    if (atomsmi.at(k).at(l)>atmp) atomsmi.at(k).at(l)-=3;
                }
            }
            for (unsigned int j=i+1;j<cid.size();j++) {
                if (cid.at(j)==cid.at(i) && apos.at(j)>apos.at(i)) rpos.at(j)-=3;
                if (apos.at(j)>apos.at(i)) apos.at(j)-=3;
            }
            if (0) {
                for (unsigned int k=0;k<atomsmi.size();k++) {
                    cout << Cindex.at(k) << " " << Mindex.at(k) << " QQ ";
                    for (unsigned int j=0;j<atomsmi.at(k).size();j++) {
                        cout << atomsmi.at(k).at(j) << " ";
                    }
                    cout << endl;
                }
            }
        }

        vector<int>().swap(apos);
        vector<unsigned int>().swap(cid);
        vector<int>().swap(rpos);
    }


    if (1) { //delete "-" of "....(-....."
        vector<int> cid(0);
        vector<int> rpos(0);
        vector<int> apos(0);

        for (unsigned int i=0;i<Bindex.size();i++) {
            unsigned int M=Mindex.at(i);
            unsigned int R=Rindex.at(i);
            int P_quota=1;
			unsigned int ct=0;
            for (unsigned int j=0;j<data->a.at(M).norder;j++) {
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

        for (unsigned int i=0;i<cid.size();i++) {
            for (unsigned int j=i+1;j<cid.size();j++) {
                if (apos.at(j)<apos.at(i)) {
                    swap(cid.at(i),cid.at(j));
                    swap(apos.at(i),apos.at(j));
                    swap(rpos.at(i),rpos.at(j));
                }
            }
        }
        for (unsigned int i=0;i<cid.size();i++) {
            //int rtmp=rpos.at(i);
            int atmp=atomsmi.at(cid.at(i)-1).at(rpos.at(i));

            atomsmi.at(cid.at(i)-1).erase(atomsmi.at(cid.at(i)-1).begin()+rpos.at(i));

            for (unsigned int k=0;k<atomsmi.size();k++) {
                for (unsigned int l=0;l<atomsmi.at(k).size();l++) {
                    if (atomsmi.at(k).at(l)>atmp) atomsmi.at(k).at(l)--;
                }
            }
            for (unsigned int j=i+1;j<cid.size();j++) {
                if (cid.at(j)==cid.at(i) && apos.at(j)>apos.at(i)) rpos.at(j)--;
                if (apos.at(j)>apos.at(i)) apos.at(j)--;
            }
            if (0) {
                for (unsigned int k=0;k<atomsmi.size();k++) {
                    cout << Cindex.at(k) << " " << Mindex.at(k) << " QQ ";
                    for (unsigned int j=0;j<atomsmi.at(k).size();j++) {
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
		unsigned int m=0;
		//vector<unsigned int> t(0);
		//vector<unsigned int> s(0);

		for (unsigned int i=0;i<Cyindex.size();i++) {
			if (Cyindex.at(i).size()) {
				for (unsigned int j=0;j<Cyindex.at(i).size();j++) {
					if (Cyindex.at(i).at(j)>m) m=Cyindex.at(i).at(j);
				}
			}
		}
		if (m==0) {
			if_circle=0;
			//for (int i=0;i<Cyindex.size();i++) Cyindex.at(i)=0;
			Cyindex.resize(0,vector<unsigned int>(0));
			Cyindex.resize(Cindex.size(),vector<unsigned int>(0));
		}


		if (m) {
			vector<int> cyclic;
			cyclic.resize(m,0);
			for (unsigned int i=0;i<Mindex.size();i++) {
				if (Cyindex.at(i).size())  {
					for (unsigned int j=0;j<Cyindex.at(i).size();j++) {
						cyclic.at(Cyindex.at(i).at(j)-1)+=1;
					}
				}
			}
			for (int i=cyclic.size()-1;i>=0;i--) { // orig i=0;i<cyclic.size();i++ //20200725
				if (cyclic.at(i)==2) {
					for (unsigned int j=0;j<Mindex.size();j++) {
						//if (1) { //20200509
						//	t.resize(0);
						//
						//	if (Cyindex.at(j).size()) { // t: each ring numbering
                    	//		for (unsigned int k=0;k<Cyindex.at(j).size();k++) {
						//			t.push_back(Cyindex.at(j).at(k));
                    	//		}
						//	}	
						//}
						//
						//bool go1=0;
						//if (1) {
						//	for (unsigned int j1=0;j1<t.size();j1++) {
						//		if ((int)t.at(j1)==(i+1)) { // orig t.at(j1)==s.at(j2) && s.at(j2)==(i+1)
						//			go1=1;
						//			break;
						//		}
						//	}
						//}
						bool go1=0;
                        if (Cyindex.at(j).size()) { // t: each ring numbering
                        	for (unsigned int k=0;k<Cyindex.at(j).size();k++) {
                                if ((int)Cyindex.at(j).at(k)==(i+1)) {
									go1=1;
									break;
								}
                            }
                        }

						if (go1) { // orig go1 && go2 20200725 // ( (s[0]==t[0] || s[0]==t[1] || s[1]==t[0] || s[1]==t[1]) && (s[0]==(i+1) || s[1]==(i+1)) )
							bool u=0;
							for (unsigned int m1=0;m1<data->a.at(Mindex.at(j)).norder;m1++) {
								if (Bindex.at(j).at(m1)==Cybnd.at(i)) { // orig Bindex[j].at(m1)==1 && Bindex[n].at(x1)==1
									Bindex.at(j).at(m1)=0;
									u=1;
									break;
								}
							}
							if (!u) {
								u=1;

                            	if (Cyindex.at(j).size()) {
                                	for (int k=0;k<(int)Cyindex.at(j).size();k++) {
                                    	if ((int)Cyindex.at(j).at(k)==(i+1)) {
											Cyindex.at(j).erase(Cyindex.at(j).begin()+k);
											k--;
										}
                                	}

        							for (unsigned int k=0;k<Cyindex.at(j).size();k++) {
                                        if ((int)Cyindex.at(j).at(k)>(i+1)) {
                                            Cyindex.at(j).at(k)--;
                                        }
        							}

                                    if ((int)Cybnd.size()>=(i+1)) Cybnd.erase(Cybnd.begin()+i);
                                    if_circle--;

                            	}
							}
						}
					}
				} 
				else {
					for (unsigned int j=0;j<Mindex.size();j++) {
						if (1) { //20200509
							if (Cyindex.at(j).size()) {
                                for (int k=0;k<(int)Cyindex.at(j).size();k++) {
									if ((int)Cyindex.at(j).at(k)==(i+1)) {
										Cyindex.at(j).erase(Cyindex.at(j).begin()+k);
										k--;
									}
                                }

                                for (unsigned int k=0;k<Cyindex.at(j).size();k++) {
                                    if ((int)Cyindex.at(j).at(k)>(i+1)) {
                                        Cyindex.at(j).at(k)--;
                                    }
                                }

                                if ((int)Cybnd.size()>=(i+1)) Cybnd.erase(Cybnd.begin()+i);
                                if_circle--;
							}
						}
					}
				}
			}
		}


		if (m) {
			m=0;
			for (unsigned int i=0;i<Cyindex.size();i++) {
				for (unsigned int j=0;j<Cyindex.at(i).size();j++) {
					if (Cyindex.at(i).at(j)>m) {
						m=Cyindex.at(i).at(j);
					}
				}
			}
			if_circle=m;
		}
	}

    if (0) {
        for (unsigned int k=0;k<Bindex.size();k++) {
            cout << Cindex.at(k) << " " << Mindex.at(k) << " QQ ";
            for (unsigned int j=0;j<Bindex.at(k).size();j++) {
                cout << Bindex.at(k).at(j) << " ";
            }
            cout << endl;
        }
    }


	//del_unpaired_ring_no();

	if (para.protect) prct(); //20200101
	//chk_cistrans();

	return 1;
}

unsigned int MOLECULE::smi2mds_OBabel() {
    stringstream smi(molesmi);
	stringstream mm("");

    OBMol *mol=NULL;
	mol = new OBMol [1];
	OpenBabel::OBConversion conv(&smi);
	if (1) {
    	conv.SetInFormat("SMI");
    	conv.Read(&mol[0]);

    	mol[0].AddHydrogens();
    	OpenBabel::PerceiveStereo(&mol[0]);
		if (1) {
  			OpenBabel::generateDiagram(&mol[0]);
  			mol[0].SetDimension(2);
		}

        if (1) {
            OpenBabel::OBBuilder builder;
            builder.Build(mol[0]);
            mol[0].SetDimension(3);
        }

		if (0) {
        	OBForceField* pFF = OBForceField::FindForceField("UFF");
        	pFF->SetLogFile(&cerr);
        	pFF->SetLogLevel(OBFF_LOGLVL_NONE);
			if (!pFF->Setup(mol[0])) cerr << "ERROR: could not setup force field." << endl;
			// Perform the actual minimization, maximum 1000 steps
			//pFF->ConjugateGradients(1000);
			// Perform the actual minimization, maximum 1000 steps
			pFF->SteepestDescent(150*mol[0].NumHvyAtoms());
			pFF->GetCoordinates(mol[0]);
        	pFF = NULL;
		}

		if (0) {
    		OpenBabel::OBBuilder builder;
    		builder.Build(mol[0]);
			//mol[0].SetDimension(3);
		}

    	OpenBabel::PerceiveStereo(&mol[0]);
	}
	OBStereoFacade facade(&mol[0]);


	if (0) {
    	string a=molesmi;
    	for (unsigned int i=0;i<a.length();i++) {
        	if (a[i]=='/') a[i]='u';
        	if (a[i]=='\\') a[i]='d';
        	if (a[i]=='*') a[i]='x';
    	}

    	conv.SetOutFormat("MOL");
		ofstream out((para.programdir+a+".mol").c_str());
    	conv.AddOption("3", OpenBabel::OBConversion::OUTOPTIONS);
    	conv.Write(&mol[0],&out);
    	out.close();
	}

	natom=0;
	vector<unsigned int> Hpos(0);
    unsigned int totatms=0;
    FOR_ATOMS_OF_MOL(atom, mol[0]) {
        if (atom->GetAtomicNum()==1) Hpos.push_back(atom->GetId());
        else natom++;
        totatms++;
    }

	vector<unsigned int> Hnum(natom,0);
	vector<unsigned int> chi1(natom,0);
    vector< vector<unsigned int> > ct(2,vector<unsigned int>(natom,0));
    vector< vector<unsigned int> > bnd(natom,vector<unsigned int> (0));

    FOR_BONDS_OF_MOL(bond, mol[0]) {
        OBAtom *bea=bond->GetBeginAtom();
        OBAtom *ena=bond->GetEndAtom();

        if (bea!=NULL && ena!=NULL && bea->GetAtomicNum()!=1 && ena->GetAtomicNum()!=1) {
            unsigned int beanum=bea->GetId();
            unsigned int enanum=ena->GetId();
            for (unsigned int i=0;i<Hpos.size();i++) {
                if (bea->GetId()>Hpos.at(i)) beanum--;
                if (ena->GetId()>Hpos.at(i)) enanum--;
            }
			bnd.at(beanum).push_back(bond->GetBondOrder());
			bnd.at(enanum).push_back(bond->GetBondOrder());
        }
        else if (bea!=NULL && ena!=NULL && bea->GetAtomicNum()!=1 && ena->GetAtomicNum()==1) {
            unsigned int beanum=bea->GetId();
            for (unsigned int i=0;i<Hpos.size();i++) {
                if (bea->GetId()>Hpos.at(i)) beanum--;
            }
            Hnum.at(beanum)++;
        }
        else if (bea!=NULL && ena!=NULL && bea->GetAtomicNum()==1 && ena->GetAtomicNum()!=1) {
            unsigned int enanum=ena->GetId();
            for (unsigned int i=0;i<Hpos.size();i++) {
                if (ena->GetId()>Hpos.at(i)) enanum--;
            }
			Hnum.at(enanum)++;
        }


		if (bond->GetBondOrder()==2 && facade.HasCisTransStereo(bond->GetId())) {
	        OBCisTransStereo *ct1=NULL;
	        ct1=facade.GetCisTransStereo(bond->GetId());
            OBCisTransStereo::Config A=ct1->GetConfig(OBStereo::ShapeU);

			vector<unsigned int> tmpct(4,0);
            for (unsigned int g=0;g<A.refs.size();g++) {
				unsigned int tm=A.refs.at(g);
				for (unsigned int u1=0;u1<Hpos.size();u1++) {
					if (A.refs.at(g)>Hpos.at(u1)) tm--;
				}

                if (tm<natom) { // cis-trans (heavy atom)
					if (g==0) tmpct.at(g)=1;
					else if (g==1) tmpct.at(g)=2;
					else if (g==2) tmpct.at(g)=1;
					else if (g==3) tmpct.at(g)=2;
				}
                else if (tm>=mol[0].NumHvyAtoms()) { // lone pair
                    if (g==0) tmpct.at(g)=1;
                    else if (g==1) tmpct.at(g)=2;
                    else if (g==2) tmpct.at(g)=1;
                    else if (g==3) tmpct.at(g)=2;
                }
			}

			if (tmpct.at(0)!=0) tmpct.at(1)=0;
			if (tmpct.at(1)!=0) tmpct.at(0)=0;
			if (tmpct.at(2)!=0) tmpct.at(3)=0;
			if (tmpct.at(3)!=0) tmpct.at(2)=0;

            unsigned int beanum=bea->GetId();
            for (unsigned int u1=0;u1<Hpos.size();u1++) {
                if (bea->GetId()>Hpos.at(u1)) beanum--;
            }
            unsigned int enanum=ena->GetId();
            for (unsigned int u1=0;u1<Hpos.size();u1++) {
                if (ena->GetId()>Hpos.at(u1)) enanum--;
            }

					
			if (tmpct.at(0)!=0 && tmpct.at(1)==0) ct.at(0).at(beanum)=tmpct.at(0);
			else if (tmpct.at(0)==0 && tmpct.at(1)!=0) ct.at(0).at(beanum)=tmpct.at(1);
			if (tmpct.at(2)!=0 && tmpct.at(3)==0) ct.at(1).at(enanum)=tmpct.at(2);
			else if (tmpct.at(2)==0 && tmpct.at(3)!=0) ct.at(1).at(enanum)=tmpct.at(3);

			ct1=NULL;
		}
		bea=NULL;
		ena=NULL;
    }


    FOR_ATOMS_OF_MOL(atom, mol[0]) {
        if (atom->GetAtomicNum()==6 && facade.HasTetrahedralStereo(atom->GetId())) {
            OBTetrahedralStereo *tr=NULL;
            tr=facade.GetTetrahedralStereo(atom->GetId());
            bool go=1;
            OBTetrahedralStereo::Config A=tr->GetConfig(OBStereo::Clockwise,OBStereo::ViewFrom);

            bool flip=0;
            if (go) {
                if (A.refs.at(0)<A.refs.at(1) && A.refs.at(1)<A.refs.at(2)) ;
                else if (A.refs.at(2)<A.refs.at(0) && A.refs.at(0)<A.refs.at(1)) ;
                else if (A.refs.at(1)<A.refs.at(2) && A.refs.at(2)<A.refs.at(0)) ;
                else flip=1;
            }

            if (go) {
                unsigned int tmp=A.center;
                for (unsigned int i=0;i<Hpos.size();i++) {
                    if (A.center>Hpos.at(i)) tmp--;
                }

                if (!flip) chi1.at(tmp)=2;
				else chi1.at(tmp)=1;
            }
			/*
            if (!go) {
                go=1;
                OBTetrahedralStereo::Config B=tr->GetConfig(OBStereo::AntiClockwise,OBStereo::ViewFrom);

                flip=0;
                if (go) {
                    if (B.refs.at(0)<B.refs.at(1) && B.refs.at(1)<B.refs.at(2)) ;
                    else if (B.refs.at(2)<B.refs.at(0) && B.refs.at(0)<B.refs.at(1)) ;
                    else if (B.refs.at(1)<B.refs.at(2) && B.refs.at(2)<B.refs.at(0)) ;
                    else flip=1;
                }

                if (go) {
	                unsigned int tmp=B.center;
    	            for (unsigned int i=0;i<Hpos.size();i++) {
        	            if (B.center>Hpos.at(i)) tmp--;
            	    }

                    if (!flip) chi1.at(tmp)=1;
					else chi1.at(tmp)=2;
                }
            }
			*/

        }
    }

	empty();
	Cyindex.resize(natom,vector<unsigned int>(0));
	if_circle=0;
	ctsisomer.resize(0,vector<string>(0));
	ctsisomer.resize(2,vector<string>(0));
	Cindex.resize(0);
	Mindex.resize(0);
	chi.resize(0);

	unsigned int Hvyatm_ct=0;
	//OpenBabel::OBElementTable ATMSYM;
	for (unsigned int i=0;i<totatms;i++) {
		OBAtom *atom=mol[0].GetAtomById(i);
		//string atmsym=ATMSYM.GetSymbol(atom->GetAtomicNum());
		string atmsym=OpenBabel::OBElements::GetSymbol(atom->GetAtomicNum());
		int atmchg=atom->GetFormalCharge();

		if (atmsym=="H" && atmchg==0) continue;
		else {
			Cindex.push_back(Hvyatm_ct+1);
			unsigned int totbnd=Hnum.at(Hvyatm_ct);
			for (unsigned int j=0;j<bnd.at(Hvyatm_ct).size();j++) totbnd+=bnd.at(Hvyatm_ct).at(j);
			unsigned int nbnd=Hnum.at(Hvyatm_ct)+bnd.at(Hvyatm_ct).size();

			if (ct.at(0).at(Hvyatm_ct)==1) ctsisomer.at(0).push_back("/");
			else if (ct.at(0).at(Hvyatm_ct)==2) ctsisomer.at(0).push_back("\\");
			else if (ct.at(0).at(Hvyatm_ct)==0) ctsisomer.at(0).push_back("");

            if (ct.at(1).at(Hvyatm_ct)==1) ctsisomer.at(1).push_back("/");
            else if (ct.at(1).at(Hvyatm_ct)==2) ctsisomer.at(1).push_back("\\");
            else if (ct.at(1).at(Hvyatm_ct)==0) ctsisomer.at(1).push_back("");

			if (chi1.at(Hvyatm_ct)==1) chi.push_back(1);
			else if (chi1.at(Hvyatm_ct)==2) chi.push_back(2);
			else if (chi1.at(Hvyatm_ct)==0) chi.push_back(0);


			if (atmsym=="C") {
				if (nbnd==4) {
					Mindex.push_back(1);
				}
				else if (nbnd==3) {
					if (totbnd==4 && atmchg==0) Mindex.push_back(2);
					else if (totbnd==3 && atmchg==-1) Mindex.push_back(59);
				}
				else if (nbnd==2) {
					if (atmchg==0) {
						for (unsigned int n=0;n<bnd.at(Hvyatm_ct).size();n++) {
							if (bnd.at(Hvyatm_ct).at(n)==3) {
								Mindex.push_back(3);
								break;
							}
							else if (bnd.at(Hvyatm_ct).at(n)==2) {
								Mindex.push_back(4);
								break;
							}
						}
					}
				}
			}
			else if (atmsym=="O") {
				if (nbnd==1) {
					if (atmchg==0) Mindex.push_back(6); 
					else if (atmchg==-1) Mindex.push_back(29);  
				}
				else if (nbnd==0) {
					if (atmchg==-1) Mindex.push_back(28);
				}
				else if (nbnd==2) {
					if (atmchg==0) Mindex.push_back(5);
				}
			}
			else if (atmsym=="N") {
				if (nbnd==4) {
					if (atmchg==1) Mindex.push_back(15);
				}
				if (nbnd==3) {
					if (totbnd==4 && atmchg==1) Mindex.push_back(16);
					else if (totbnd==3 && atmchg==0) Mindex.push_back(7);
				}
				if (nbnd==2) {
					if (totbnd==3 && atmchg==0) Mindex.push_back(8);
					else if (totbnd==2 && atmchg==-1) Mindex.push_back(35);
					else if (totbnd==2 && atmchg==0) Mindex.push_back(7);
				}
				if (nbnd==1) {
					if (atmchg==0) Mindex.push_back(9);
				}
			}
			else if (atmsym=="P") {
				if (nbnd==6) {
					if(atmchg==-1) Mindex.push_back(30);
				}
				if (nbnd==5) {
					if (atmchg==0) Mindex.push_back(31);
				}
				if (nbnd==4) {
					if (totbnd==4 && atmchg==1) Mindex.push_back(17);
					else if (totbnd==5 && atmchg==0) Mindex.push_back(32);
					else if (totbnd==4 && atmchg==0) Mindex.push_back(69);
				}
				if (nbnd==3) {
					if (totbnd==4 && atmchg==1) Mindex.push_back(18);
					if (totbnd==4 && atmchg==0) Mindex.push_back(66);
					else if (totbnd==3 && atmchg==0) Mindex.push_back(21);
				}
				if (nbnd==2) {
					if (atmchg==0) Mindex.push_back(22);
				}
				if (nbnd==1) {
					if (atmchg==0) Mindex.push_back(23);
				}
			}
			else if (atmsym=="F") {
				if (nbnd==1) {
					if (atmchg==0) Mindex.push_back(11);
					else if (atmchg==-1) Mindex.push_back(24);
				}
				else if (nbnd==0) {
					if (atmchg==-1) Mindex.push_back(24);
				}
			}
			else if (atmsym=="Cl") {
				if (nbnd==1) {
					if (atmchg==0) Mindex.push_back(12);
					else if (atmchg==-1) Mindex.push_back(25);
				}
				else if (nbnd==0) {
					if (atmchg==-1) Mindex.push_back(25);
				}
				else if (nbnd==3) {
					if (totbnd==6 && atmchg==0) Mindex.push_back(62);
				}
				else if (nbnd==4) {
					if (atmchg==0) Mindex.push_back(62);
				}
			}
			else if (atmsym=="Br") {
				if (nbnd==1) {
					if (atmchg==0) Mindex.push_back(13);
					else if (atmchg==-1) Mindex.push_back(26);
				}
				else if (nbnd==0) {
					if (atmchg==-1) Mindex.push_back(26);
				}
			}
			else if (atmsym=="I") {
				if (nbnd==1) {
					if (atmchg==0) Mindex.push_back(14);
					else if (atmchg==-1) Mindex.push_back(27);
				}
				else if (nbnd==0) {
					if (atmchg==-1) Mindex.push_back(27);
				}
			}
			else if (atmsym=="S") {
				if (nbnd==2) {
					if (atmchg==0) Mindex.push_back(19);
				}
				else if (nbnd==1) {
					if (totbnd==2 && atmchg==0) Mindex.push_back(20);
					else if (totbnd==1 && atmchg==-1) Mindex.push_back(63);
				}
				else if (nbnd==3) {
					if (totbnd==5 && atmchg==0)  Mindex.push_back(61);
					else if (totbnd==4 && atmchg==0) Mindex.push_back(34);
					else if (totbnd==3 && atmchg==1) Mindex.push_back(65);
				}
				else if (nbnd==4) {
					if (totbnd==6 && atmchg==0) Mindex.push_back(61);
				}
			}
			else if (atmsym=="B") {
				if (nbnd==4) {
					if (atmchg==-1) Mindex.push_back(44);
				}
			}
			else if (atmsym=="In") {
				if (nbnd==4) {
					if (atmchg==3) Mindex.push_back(57);
				}
			}
			else if (atmsym=="Ga") {
				if (nbnd==4) {
					if (atmchg==3) Mindex.push_back(64);
				}
			}
			else if (atmsym=="Xx") {
		    	if (nbnd==1) {
		        	if (atmchg==0) Mindex.push_back(70);
		    	}
			}

			Hvyatm_ct++;
		}
		atom=NULL;
	}

	vector<unsigned int>().swap(Hnum);
	vector<unsigned int>().swap(chi1);
	vector< vector<unsigned int> >().swap(ct);
	vector< vector<unsigned int> >().swap(bnd);

	Rindex.resize(natom,0);
	Pindex.resize(natom,0);

	vector<bool> hasP(natom,0);
	unsigned int ringnum=1;
    FOR_BONDS_OF_MOL(bond, mol[0]) {
        OBAtom *bea=bond->GetBeginAtom();
        OBAtom *ena=bond->GetEndAtom();

		if (bea->GetAtomicNum()!=1 && ena->GetAtomicNum()!=1) {
			unsigned int maxind=max(bea->GetId(),ena->GetId());
			unsigned int minind=min(bea->GetId(),ena->GetId());
        	for (unsigned int i=0;i<Hpos.size();i++) {
            	if (max(bea->GetId(),ena->GetId())>Hpos.at(i)) maxind--;
				if (min(bea->GetId(),ena->GetId())>Hpos.at(i)) minind--;
        	}

			if (!hasP.at(maxind)) {
				Rindex.at(maxind)=bond->GetBondOrder();
				Pindex.at(maxind)=Cindex.at(minind);
				hasP.at(maxind)=1;
			}
			else {
				Cyindex.at(minind).push_back(ringnum);
				Cyindex.at(maxind).push_back(ringnum);
				Cybnd.resize(ringnum,0);
				Cybnd.at(ringnum-1)=bond->GetBondOrder();
				ringnum++;
			}
		}
		bea=NULL;
		ena=NULL;
	}
	if_circle=Cybnd.size();

	vector<bool>().swap(hasP);
	vector<unsigned int>().swap(Hpos);

    mol[0].Clear();

	if (mol!=NULL) {
		delete [] mol;
		mol = NULL;
	}

	if (para.protect) prct();

	//mds2smi();
	//printmds();

	if (0) {
		for (unsigned int i=0;i<Cindex.size();i++) {
			if (Rindex.at(i)==2) {
				int P=Pindex.at(i);
				if (P>0) {
					if (ctsisomer.at(1).at(P-1)!="" && ctsisomer.at(0).at(P-1)=="") {
						swap(ctsisomer.at(1).at(P-1),ctsisomer.at(0).at(P-1));
					}
				}
				if (ctsisomer.at(0).at(i)!="" && ctsisomer.at(1).at(i)=="") {
					swap(ctsisomer.at(1).at(i),ctsisomer.at(0).at(i));
				}
			}
		}
	}

	if (Cindex.size()>2) {
		if (ctsisomer.at(0).at(0)!="" && ctsisomer.at(1).at(1)!="" && Rindex.at(1)==2 && Pindex.at(1)==1) {
			int ct1=0,ct2=0;
        	for (unsigned int k1=0;k1<Cindex.size();k1++) {
            	if (Pindex.at(k1)==Cindex.at(0) && Cindex.at(k1)!=Cindex.at(1)) {
					if (ct1==0) {
                        if (ctsisomer.at(0).at(0)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
                        else if (ctsisomer.at(0).at(0)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
					}
					ct1++;
				}
				if (Pindex.at(k1)==Cindex.at(1)) {
					if (ct2==0) {
						ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(1).at(1)+ctsisomer.at(0).at(Cindex.at(k1)-1);
					}
					ct2++;
				}
			}
			ctsisomer.at(0).at(0)="";
			ctsisomer.at(1).at(1)="";
		}
	}

	//mds2smi();

	return 1;

}

unsigned int MOLECULE::clean() {
	unsigned int i;
	if (connect==NULL && atm==NULL) return 0;
	if (connect!=NULL) {
		for (i=0;i<natom;i++) delete [] connect[i]; //if(connect[i]!=NULL)
		delete [] connect;
		connect=NULL;
	}
	if(atm!=NULL) {
		delete [] atm;
		atm=NULL;
	}
	return 1;
}

unsigned int MOLECULE::combination(MOLECULE &B,unsigned int k,unsigned int p,unsigned int bnd) {
	if (k>=Cindex.size() || k<0) return 0;
	if (p>=B.Cindex.size() || p<0) return 0;
	if (!data->a[Mindex.at(k)].bd[bnd-1]) return 0;
	if (!B.data->a[B.Mindex.at(p)].bd[bnd-1]) return 0;

	//int nu=B.Cindex.at(p);

    if (1) {
        //bool bndchk=0;

        int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        vector<unsigned int> ordcount(6,0);
        if (Rindex.at(k)>0) {
			bndsum+=Rindex.at(k);
			ordcount.at(Rindex.at(k)-1)+=1;
		}
        if (Cyindex.at(k).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(k).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(k).at(k1)-1);
				//cout << "Cybnd.at(Cyindex.at(k).at(k1)-1) " << Cybnd.at(Cyindex.at(k).at(k1)-1) << endl;
                ordcount.at(Cybnd.at(Cyindex.at(k).at(k1)-1)-1)++;
            }

        }
        for (unsigned int k1=0;k1<Cindex.size();k1++) {
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

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[Mindex.at(k)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(k)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

    }

    if (1) {
        //bool bndchk=0;

        int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        vector<unsigned int> ordcount(6,0);
        if (B.Rindex.at(p)>0) {
			bndsum+=B.Rindex.at(p);
			ordcount.at(B.Rindex.at(p)-1)+=1;
		}
        if (B.Cyindex.at(p).size()) {
            for (unsigned int k1=0;k1<B.Cyindex.at(p).size();k1++) {
                bndsum+=B.Cybnd.at(B.Cyindex.at(p).at(k1)-1);
				//cout << "B.Cybnd.at(B.Cyindex.at(p).at(k1)-1) " << B.Cybnd.at(B.Cyindex.at(p).at(k1)-1) << endl;
                ordcount.at(B.Cybnd.at(B.Cyindex.at(p).at(k1)-1)-1)++;
            }
        }
        for (unsigned int k1=0;k1<B.Cindex.size();k1++) {
            if (B.Pindex.at(k1)==B.Cindex.at(p)) {
                bndsum+=B.Rindex.at(k1);
                ordcount.at(B.Rindex.at(k1)-1)+=1;
            }
        }
        bndsum+=bnd;
        ordcount.at(bnd-1)+=1;

        unsigned int maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[B.Mindex.at(p)].order.at(k2);
        }

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a[B.Mindex.at(p)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
				//cout << "2nd molecule " << B.Cindex.at(p) << " | Bond order: " << (k2+1) << " | Required number: " << ordcount.at(k2) << " | Available number: " << data.a[B.Mindex.at(p)].bd[k2]  << endl;
                //bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[B.Mindex.at(p)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

    }


	vector<unsigned int> tmpP(B.Pindex.size(),0);
	vector<unsigned int> tmpC(B.Cindex.size(),0);
	vector<unsigned int> tmpM(B.Mindex.size(),0);
	vector<unsigned int> tmpR(B.Rindex.size(),0);
	vector<unsigned int> tmpchi(B.chi.size(),0);
	vector< vector<unsigned int> > tmpCy(B.Cyindex.size(),vector<unsigned int>(0));
	vector< vector<string> > tmpCT(2,vector<string>(B.ctsisomer.at(0).size(),""));
	vector<unsigned int> tmpPr(B.protect.size(),0);
	for (unsigned int k2=0;k2<B.Pindex.size();k2++) tmpP.at(k2)=B.Pindex.at(k2);
	for (unsigned int k2=0;k2<B.Cindex.size();k2++) tmpC.at(k2)=B.Cindex.at(k2);
	for (unsigned int k2=0;k2<B.Mindex.size();k2++) tmpM.at(k2)=B.Mindex.at(k2);
	for (unsigned int k2=0;k2<B.Rindex.size();k2++) tmpR.at(k2)=B.Rindex.at(k2);
	for (unsigned int k2=0;k2<B.chi.size();k2++) tmpchi.at(k2)=B.chi.at(k2);
	for (unsigned int k2=0;k2<B.ctsisomer.at(0).size();k2++) tmpCT.at(0).at(k2)=B.ctsisomer.at(0).at(k2);
	for (unsigned int k2=0;k2<B.ctsisomer.at(1).size();k2++) tmpCT.at(1).at(k2)=B.ctsisomer.at(1).at(k2);
	for (unsigned int k2=0;k2<B.Cyindex.size();k2++) {
		for (unsigned int k3=0;k3<B.Cyindex.at(k2).size();k3++) tmpCy.at(k2).push_back(B.Cyindex.at(k2).at(k3));
	}
	if (para.protect) for (unsigned int k2=0;k2<B.protect.size();k2++) tmpPr.at(k2)=B.protect.at(k2);

	if (p!=0) {
		tmpP.at(p)=0;
    	unsigned int compa=B.Cindex.at(p);
    	unsigned int cont=0;
    	while (compa>=1 && cont<=B.Cindex.size()) {
        	//if (1) C_ringmember.push_back(compa);

        	if (compa>=1) {
            	if (B.Pindex.at(compa-1)>=1) {
					tmpP.at(B.Pindex.at(compa-1)-1)=B.Cindex.at(compa-1);
					tmpR.at(B.Pindex.at(compa-1)-1)=B.Rindex.at(compa-1);

					compa=B.Cindex.at(B.Pindex.at(compa-1)-1);
				}
            	else break;
            	cont++;
        	}
        	else break;
    	}
		tmpR.at(p)=0;
	}


	if (p!=0) {
		swap(tmpP.at(p),tmpP.at(0));
		swap(tmpC.at(p),tmpC.at(0));
		swap(tmpM.at(p),tmpM.at(0));
		swap(tmpR.at(p),tmpR.at(0));
		swap(tmpchi.at(p),tmpchi.at(0));
		swap(tmpCy.at(p),tmpCy.at(0));
		swap(tmpCT.at(0).at(p),tmpCT.at(0).at(0));
		swap(tmpCT.at(1).at(p),tmpCT.at(1).at(0));
		if (para.protect) swap(tmpPr.at(p),tmpPr.at(0));
	}
	
	for (int i=0;i<(int)tmpC.size();i++) {
		for (int j=i+1;j<(int)tmpC.size();j++) {
			if (tmpP.at(i)==tmpC.at(j)) {
				swap(tmpP.at(j),tmpP.at(i));
                swap(tmpC.at(j),tmpC.at(i));
                swap(tmpM.at(j),tmpM.at(i));
                swap(tmpR.at(j),tmpR.at(i));
				swap(tmpchi.at(j),tmpchi.at(i));
                swap(tmpCy.at(j),tmpCy.at(i));
				swap(tmpCT.at(0).at(j),tmpCT.at(0).at(i));
				swap(tmpCT.at(1).at(j),tmpCT.at(1).at(i));
				if (para.protect) swap(tmpPr.at(j),tmpPr.at(i));

				i--;
				break;
			}
		}
	}
	
    for (unsigned int i=0;i<tmpC.size();i++) {
        for (unsigned int j=i+1;j<tmpC.size();j++) {
			if (tmpP.at(i)==tmpC.at(j)) { // && tmpC.at(i)==tmpP.at(j)
			    vector<unsigned int>().swap(tmpP);
			    vector<unsigned int>().swap(tmpC);
			    vector<unsigned int>().swap(tmpM);
			    vector<unsigned int>().swap(tmpR);
				vector<unsigned int>().swap(tmpchi);
			    vector< vector<unsigned int> >().swap(tmpCy);
			    vector<unsigned int>().swap(tmpPr);
				vector< vector<string> >().swap(tmpCT);

                return 0;
			}
        }
    }
	
    for (unsigned int i=0;i<tmpC.size();i++) {
        unsigned int j=tmpC.at(i);
        tmpC.at(i)=(i+1)+tmpC.size();
        for (unsigned int k1=0;k1<tmpC.size();k1++) { 
            if (tmpP.at(k1)==j) tmpP.at(k1)=tmpC.at(i);
        }
    }

	for (unsigned int i=0;i<tmpC.size();i++) {
		tmpC.at(i)-=tmpC.size();
		if (i!=0) tmpP.at(i)-=tmpC.size();
	}


    if (1) {
        for (unsigned int i=1;i<tmpC.size();i++) {
            if (tmpR.at(i)==2 && tmpP.at(i)>0) {
                if (tmpCT.at(1).at(tmpP.at(i)-1)!="" && tmpCT.at(0).at(i)!="") {
                    swap(tmpCT.at(0).at(i),tmpCT.at(1).at(i));
                    swap(tmpCT.at(0).at(tmpP.at(i)-1),tmpCT.at(1).at(tmpP.at(i)-1));
                }
            }
        }
    }


	int Cshift=Cindex.size();
    for (unsigned int i=0;i<tmpC.size();i++) {
        tmpC.at(i)+=Cshift;  //at
        if (i==0) { // orig i==p
        	tmpP.at(i)=Cindex.at(k); //at
        	tmpR.at(i)=bnd;
        }
        else tmpP.at(i)+=Cshift; //at
    }


	Cybnd.resize(Cybnd.size()+B.if_circle,0);
    for (unsigned int k2=0;k2<tmpCy.size();k2++) {
        if (tmpCy.at(k2).size()) {
            for (unsigned int k3=0;k3<tmpCy.at(k2).size();k3++) {
				Cybnd.at(tmpCy.at(k2).at(k3)-1+if_circle)=B.Cybnd.at(tmpCy.at(k2).at(k3)-1);
				tmpCy.at(k2).at(k3)+=if_circle;
            }
        }
    }
	if_circle=B.if_circle+if_circle;

    for (unsigned int i=0;i<tmpC.size();i++) {
        Pindex.push_back(tmpP.at(i));
        Cindex.push_back(tmpC.at(i));
        Rindex.push_back(tmpR.at(i));
		chi.push_back(tmpchi.at(i));
		Cyindex.resize(Cyindex.size()+1,vector<unsigned int>(0));
		for (unsigned int j=0;j<tmpCy.at(i).size();j++) Cyindex.at(Cyindex.size()-1).push_back(tmpCy.at(i).at(j));
        Mindex.push_back(tmpM.at(i));
		ctsisomer.at(0).push_back(tmpCT.at(0).at(i));
		ctsisomer.at(1).push_back(tmpCT.at(1).at(i));
        if (para.protect) protect.push_back(tmpPr.at(i));
    }

	//printmds();

	vector<unsigned int>().swap(tmpP);
	vector<unsigned int>().swap(tmpC);
	vector<unsigned int>().swap(tmpM);
	vector<unsigned int>().swap(tmpR);
	vector<unsigned int>().swap(tmpchi);
	vector< vector<unsigned int> >().swap(tmpCy);
	vector< vector<string> >().swap(tmpCT);
	vector<unsigned int>().swap(tmpPr);

	//chk_cistrans(); //20200806
	if (!chk_imine_ct(0,Cindex.size()-1)) mds2smi();
	//chk_imine_ct(0,Cindex.size()-1);
	//mds2smi();

    return 1;

}

unsigned int MOLECULE::change_chirality(unsigned int pos,unsigned int spec) {
	bool go=0;
	if (data->a.at(Mindex.at(pos)).index>=6 && data->a.at(Mindex.at(pos)).nbond!=data->a.at(Mindex.at(pos)).index && chi.at(pos)) go=1;
    if (go && spec!=chi.at(pos)) {
        if (spec==2 || spec==1) {
			chi.at(pos)=spec;

			if (!chk_imine_ct(0,Cindex.size()-1)) mds2smi();
			//chk_imine_ct(0,Cindex.size()-1);
			//mds2smi();
        }
        else return 0;
    }
    else return 0;
    return 1;

}


unsigned int MOLECULE::change_cistrans(unsigned int pos,unsigned int he) {
    if (ctsisomer.size()==2) {
        if (ctsisomer.at(0).size()==Cindex.size() && ctsisomer.at(0).size()==ctsisomer.at(1).size()) {
            if (ctsisomer.at(he).at(pos)!="") {
                if (ctsisomer.at(he).at(pos)=="/") ctsisomer.at(he).at(pos)="\\";
                else if (ctsisomer.at(he).at(pos)=="\\") ctsisomer.at(he).at(pos)="/";

				if (!chk_imine_ct(0,Cindex.size()-1)) mds2smi();
				//chk_imine_ct(0,Cindex.size()-1);
				//mds2smi();
            }
            else return 0;
        }
        else return 0;
    }
    else return 0;
    return 1;
}

unsigned int MOLECULE::wipe() {
	vector<unsigned int>().swap(Cindex);
	vector<unsigned int>().swap(Pindex);
	vector< vector<unsigned int> >().swap(Cyindex);
	vector<unsigned int>().swap(Cybnd);
	vector<unsigned int>().swap(Rindex);
	vector<unsigned int>().swap(chi);
	vector<unsigned int>().swap(Mindex);
	vector<unsigned int>().swap(protect);
	vector< vector<unsigned int> >().swap(Bindex);
	vector< vector<int> >().swap(atomsmi);
	vector< vector<string> >().swap(ctsisomer);

	clean();

	molesmi.clear();
	smiles.clear();
	return 1;
}


unsigned int MOLECULE::input(bool ct_on) {
	empty();
	clean(); 
	if (1) smi2mds_OBabel();
	clean();
	mds2smi(2,0,1); 
	canonicalize_SMILES();

	return 1;
}


unsigned int MOLECULE::empty() {
	Cindex.resize(0);
	Pindex.resize(0);
	Mindex.resize(0);
	Rindex.resize(0);
	chi.resize(0);
	Cyindex.resize(0,vector<unsigned int> (0));
	Cybnd.resize(0);
	protect.resize(0);
	ctsisomer.resize(0,vector<string> (0));
	Bindex.resize(0,vector<unsigned int> (0));
	atomsmi.resize(0,vector<int> (0));
	return 1;
}

unsigned int MOLECULE::addition(unsigned int pt,unsigned int id,unsigned int bnd, bool cistrans) {
	if (pt<0 || pt>=Cindex.size()) return 0;
	if (bnd>3 || bnd<=0) return 0;
	if (!data->a.at(Mindex.at(pt)).bd[bnd-1]) return 0;
	if (!data->a.at(id).bd[bnd-1]) return 0;

	int chgg=0;
	for (unsigned int i=0;i<Mindex.size();i++) chgg+=data->a.at(Mindex.at(i)).chg;
	if (para.ion) {
		if (data->a.at(id).chg*chgg<0) return 0;
	}
	else {
		if (data->a.at(id).chg!=0) return 0;
	}

	bool isnotend=0;

    // Ensure the bonds of atom P(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
    if (1) {
        //bool bndchk=0;

        int bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        vector<unsigned int> ordcount(6,0);
        if (Pindex.at(pt)>0) {
			bndsum+=Rindex.at(pt);
			ordcount.at(Rindex.at(pt)-1)+=1;
		}
        if (Cyindex.at(pt).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(pt).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(pt).at(k1)-1);
				//cout << "Cybnd.at(Cyindex.at(pt).at(k1)-1) " << Cybnd.at(Cyindex.at(pt).at(k1)-1) << endl;
                ordcount.at(Cybnd.at(Cyindex.at(pt).at(k1)-1)-1)+=1;
            }
        }
        for (unsigned int k1=0;k1<Cindex.size();k1++) {
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
            maxbnd+=data->a.at(Mindex.at(pt)).order.at(k2);
        }
        /*
        if (bndsum<=maxbnd) totbndchk=1;
        else {
            totbndchk=0;
            return 0;
        }
        */

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (ordcount.at(k2)>data->a.at(Mindex.at(pt)).bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a.at(Mindex.at(pt)).bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

    }

	if (0) {
		if (pt>=Cindex.size()-1) {
			Pindex.push_back(Cindex.at(pt));
			Cindex.push_back(Cindex.size()+1);
			Mindex.push_back(id);
			Rindex.push_back(bnd);
			chi.push_back(0);
			Cyindex.resize(Cyindex.size()+1,vector<unsigned int>(0));
			if (para.protect) protect.push_back(0);
			ctsisomer.at(0).push_back("");
			ctsisomer.at(1).push_back("");
		}
		else {
			for (unsigned int i=0;i<Cindex.size();i++) {
				unsigned int tmp=Cindex.at(i);
				Cindex.at(i)+=(Cindex.size()+2);
				for (unsigned int j=i+1;j<Cindex.size();j++) {
					if (Pindex.at(j)==tmp) Pindex.at(j)=Cindex.at(i);
				}
			}

			Pindex.insert(Pindex.begin()+pt+1,Cindex.at(pt));
			Cindex.insert(Cindex.begin()+pt+1,Cindex.size()+1);
			Mindex.insert(Mindex.begin()+pt+1,id);
			Rindex.insert(Rindex.begin()+pt+1,bnd);
			chi.insert(chi.begin()+pt+1,0);
			Cyindex.insert(Cyindex.begin()+pt+1,vector<unsigned int>(0));
			if (para.protect) protect.insert(protect.begin()+pt+1,0);
			ctsisomer.at(0).insert(ctsisomer.at(0).begin()+pt+1,"");
			ctsisomer.at(1).insert(ctsisomer.at(1).begin()+pt+1,"");

			for (unsigned int i=0;i<Cindex.size();i++) {
				unsigned int tmp=Cindex.at(i);
            	Cindex.at(i)=(i+1);
            	for (unsigned int j=i+1;j<Cindex.size();j++) {
                	if (Pindex.at(j)==tmp) Pindex.at(j)=Cindex.at(i);
            	}
        	}
		
		}
	}
	if (1) {
		unsigned int Csize=Cindex.size();
		Pindex.push_back(Cindex.at(pt));
		Cindex.push_back(Csize+1);
		Mindex.push_back(id);
		Rindex.push_back(bnd);
		chi.push_back(0);
		if (para.protect) protect.push_back(0);
		Cyindex.resize(Csize+1,vector<unsigned int>(0));
		ctsisomer.at(0).push_back("");
		ctsisomer.at(1).push_back("");
	}

	if (1) {
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

	//reset();
	//chk_cistrans(); //20200806
	if (!chk_imine_ct(0,Cindex.size()-1,cistrans)) mds2smi();
	//chk_imine_ct(0,Cindex.size()-1,cistrans);
	//mds2smi();

	return 1;
}

unsigned int MOLECULE::insertion(unsigned int n,unsigned int id,unsigned int bnd2par,unsigned int bnd2des,bool cistrans) {
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
    for (unsigned int i=0;i<Cindex.size();i++) chgg+=data->a[Mindex.at(i)].chg;
    if (para.ion) {
		if (data->a[id].chg*chgg<0) return 0;
    	//else if (data->a[id].chg*chgg>=0 && data->a[id].chg!=0) {
        //	if (data->a[id].index<=2 && data->a[Mindex.at(n)].index>=7) return 0;
       	//	if (data->a[id].index>=7 && data->a[Mindex.at(n)].index<=7) return 0;
    	//}
	}
	else {
		if (data->a[id].chg!=0) return 0;
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
        //bool bndchk=0;

        int bndsum=0; // max available # of bonds of C(n) (calc. from pool)
        bndsum+=bnd2des;
        vector<unsigned int> ordcount(6,0);
        ordcount.at(bnd2des-1)+=1;
        if (Cyindex.at(n).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(n).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(n).at(k1)-1);
				//cout << "Cybnd.at(Cyindex.at(n).at(k1)-1) " << Cybnd.at(Cyindex.at(n).at(k1)-1) << endl;
                ordcount.at(Cybnd.at(Cyindex.at(n).at(k1)-1)-1)++;
            }
        }
        for (unsigned int k1=0;k1<Cindex.size();k1++) {
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
                //bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[M].bd[k2]) { 
                //bndchk=1;
            }
        }

    }
    // Ensure the bonds of atom P(n) (before insertion) are OK after atm(id) inserted and use a bnd2par bond to connect P(n)
    if (P>0) {
        //bool bndchk=0;

        int bndsum=0;  // # of bonds of atom P(n) if atom M(n) is changed to atom(id) and a $bnd2par bond of atom(id) is used to connected with P(n).
        if (Rindex.at(P-1)>0) bndsum+=Rindex.at(P-1);
        vector<unsigned int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
        if (Rindex.at(P-1)>0) ordcount.at(Rindex.at(P-1)-1)+=1;  //20200712 // bnd P(n) to P(P(n))
        if (Cyindex.at(P-1).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(P-1).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(P-1).at(k1)-1);
				//cout << "Cyindex.at(P-1).at(k1) " << Cyindex.at(P-1).at(k1) << endl;
                ordcount.at(Cybnd.at(Cyindex.at(P-1).at(k1)-1)-1)++;
            }
        }
        for (unsigned int k=0;k<Cindex.size();k++) {
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
                //bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(P-1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

        //if (!totbndchk || !bndchk) {
            //cout << "Failed exchange: (totbndchk, bndchk) = (" << totbndchk << ", " << bndchk << ")" << endl;
            //return 0;
        //}
    }
    // Ensure the bonds of atm(id) are OK if atm(id) is inserted and use bnd2par & bnd2des bonds to connect P(n) & C(n), respectively.
    if (1) {
        //bool bndchk=0;

        int id_bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        if (P>0) id_bndsum+=bnd2par;
        vector<unsigned int> id_ordcount(6,0);
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
                //bndchk=0;
                return 0;
            }
            else if (k2>=2 && id_ordcount.at(k2)<=data->a[id].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

    }


	for (unsigned int k3=0;k3<Cindex.size();k3++) {
		Cindex.at(k3)+=(Cindex.size()+2);
		Pindex.at(k3)+=(Cindex.size()+2);
	}

	vector<unsigned int> curatm(0);

	bool goout=0;
	Cindex.at(n)+=(Cindex.size()+2);
    Pindex.at(n)+=(Cindex.size()+2);
    Rindex.at(n)=bnd2des;

	curatm.push_back(C+Cindex.size()+2);

	do {
		vector<unsigned int> tmpp(0);
        for (unsigned int k3=0;k3<curatm.size();k3++) {
            for (unsigned int k2=0;k2<Cindex.size();k2++) {
                if (Pindex.at(k2)==curatm.at(k3)) {
					Pindex.at(k2)+=(Cindex.size()+2);
                    tmpp.push_back(Cindex.at(k2));
					Cindex.at(k2)+=(Cindex.size()+2);
                }
            }
		}
        if (tmpp.size()>0) {
            curatm.resize(0);
            for (unsigned int k4=0;k4<tmpp.size();k4++) curatm.push_back(tmpp.at(k4));
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
		chi.insert(chi.begin()+C-1,0);
        //Cyindex.insert(Cyindex.begin()+C-1,0);
		Cyindex.insert(Cyindex.begin()+C-1,vector<unsigned int>(0));
        //ctsisomer.insert("");
        ctsisomer.at(0).insert(ctsisomer.at(0).begin()+C-1,"");
        ctsisomer.at(1).insert(ctsisomer.at(1).begin()+C-1,"");
        if (para.protect) protect.insert(protect.begin()+C-1,0);		
	}

	Pindex.at(n+1)=Cindex.at(n);
	Pindex.at(0)=0;

    for (unsigned int i=0;i<Cindex.size();i++) {
        unsigned int j=Cindex.at(i);
        Cindex.at(i)=(i+1);
        for (unsigned int k1=0;k1<Cindex.size();k1++) {
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
    //decyc_small_ring(5);

    //reset();
    //chk_cistrans(); //20200806
	if (!chk_imine_ct(0,Cindex.size()-1,cistrans)) mds2smi();
	//chk_imine_ct(0,Cindex.size()-1,cistrans);
	//mds2smi();

	return 1;
}

unsigned int MOLECULE::subtraction(unsigned int n,unsigned int bndfrm,bool cistrans) {
	// mode 1: subtract an atom
	// mode 2: subtract all atoms between Cindex=n and the end of the corresponding branch
	//vector<int> ref,pr;

	//long double time1=time(NULL);

	//if (Cindex.size()<2) return 0;
    //unsigned int Csize1=Cindex.size();
    //for (unsigned int i=0;i<Cindex.size();i++) {
    //    if (Mindex.at(i)==70) Csize1--;
    //}
	//if (Csize1<2) return 0;
	if (n>=Cindex.size() || Cindex.size()<=1) return 0;
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
	
	if (bndfrm && !Pindex.at(n)) return 0;

    if (1) if (para.protect && protect.at(n)) return 0;   // terminate if it would subtract protected atom
    if (1) if (data->a.at(Mindex.at(n)).chg) return 0;  // terminate if it would subtract charged atom

	//bool usepibnd=0;
	bool isnotend=0;

	vector<unsigned int> Merge_C(0);
	//bool bndchk=0;
	unsigned int P=Pindex.at(n);


	if (P>0) {
		Merge_C.push_back(Cindex.at(P-1));
		//bndchk=0;

		unsigned int bndsum=0; // # of bonds of atom P if atom j is deleted and P is connected with each of the descendant atom k
		if (Rindex.at(P-1)>0) bndsum+=Rindex.at(P-1);  
		vector<unsigned int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
		if (Rindex.at(P-1)>0) ordcount.at(Rindex.at(P-1)-1)+=1;  //20200712
		if (Cyindex.at(P-1).size()) {
			for (unsigned int k1=0;k1<Cyindex.at(P-1).size();k1++) {
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
		for (unsigned int k=0;k<Cindex.size();k++) { 
			if (Pindex.at(k)==Cindex.at(n)) { // For each descendant atom k of atom j // orig Pindex.at(k)==Cindex.at(j)
				isnotend=1;
				Merge_C.push_back(Cindex.at(k));
				if (bndfrm>0) {
					bndsum+=bndfrm; //bndsum+=Rindex.at(k);
					ordcount.at(bndfrm-1)++; //ordcount.at(Rindex.at(k)-1)++;
				}
				else {
					//bndsum+=Rindex.at(k);
					//ordcount.at(Rindex.at(k)-1)++;
					return 0;
				}

				//if (Rindex.at(k)==2) usepibnd=1;
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

		for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
			if (ordcount.at(k2)>data->a[Mindex.at(P-1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(P-1)].order.at(k2)
				//bndchk=0;
				return 0;
			}
			else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(P-1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(P-1)].order.at(k2)
				//bndchk=1;
			}
		}

	}
	else {
		Merge_C.push_back(0);
		for (unsigned int k=0;k<Cindex.size();k++) {
			if (Pindex.at(k)==Cindex.at(n)) { // For each descendant atom k of atom j // orig Pindex.at(k)==Cindex.at(j)
				Merge_C.push_back(Cindex.at(k));
			}
		}
		if (Merge_C.size()!=2) return 0;
	}
	if (!isnotend && bndfrm) {
		return 0;
	}

    // Ensure the bonds of atom des(n) are OK if bond order of des(n) is changed to $bndfrm.
	if (isnotend) {
    	for (unsigned int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(n) && Mindex.at(k1)!=70) { // For the descendant atoms C(k1) of atom C(n)
            	//bool bndchk=0;

				int bndsum=0; // # of bonds of atom C(k1) if its bond order R(k1) changes into .
            	vector<unsigned int> ordcount(6,0); // count # of bonds atom C(k1) will use after connecting with atom C(n) (in each bond order respectively)

                if (bndfrm>0) {
                    bndsum+=bndfrm; 
                    ordcount.at(bndfrm-1)++; 
                }
                else {
                    //bndsum+=Rindex.at(k1);
                    //ordcount.at(Rindex.at(k1)-1)++;
					if (P) return 0;
                }

	        	if (Cyindex.at(k1).size()) {
		            for (unsigned int k2=0;k2<Cyindex.at(k1).size();k2++) {
                		bndsum+=Cybnd.at(Cyindex.at(k1).at(k2)-1);
                		ordcount.at(Cybnd.at(Cyindex.at(k1).at(k2)-1)-1)++;
            		}
        		}
            	for (unsigned int k=0;k<Cindex.size();k++) { // For each des(k1) atom
                    if (Pindex.at(k)==Cindex.at(k1)) {
                        bndsum+=Rindex.at(k); // calc. the total # of bonds des(k1) atoms used to connect with C(k1).
                        ordcount.at(Rindex.at(k)-1)+=1;
                    }
                }

                int maxbnd=0; // max available # of bonds for C(k1) (calc. from pool)
                for (int k2=0;k2<6;k2++) {
                    maxbnd+=data->a[Mindex.at(k1)].order.at(k2);
                }

                for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with des(k1) atoms" match the definition of atom k1 in pool
                    if (ordcount.at(k2)>data->a[Mindex.at(k1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                        //bndchk=0;
						return 0;
                    }
                    else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(k1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                        //bndchk=1;
                    }
                }

            }
        }

    }


	if (Cyindex.at(n).size()) {
		if (Cyindex.at(n).size()>1) sort(Cyindex.at(n).begin(),Cyindex.at(n).end());

		for (int g=Cyindex.at(n).size()-1;g>=0;g--) {
    		for (int k=0;k<(int)Cyindex.size();k++) {  // 20200131 prevent unpaired ring number
				if (Cyindex.at(k).size() && k!=(int)n) { 
					for (int k1=0;k1<(int)Cyindex.at(k).size();k1++) {
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
            for (int k=0;k<(int)Cyindex.size();k++) {  // 20200131 prevent unpaired ring number
                if (Cyindex.at(k).size() && k!=(int)n) { // orig Cyindex.at(j) && Cyindex.at(k) && k!=j
                    for (int k1=0;k1<(int)Cyindex.at(k).size();k1++) {
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
            for (unsigned int k1=0;k1<Cindex.size();k1++) {
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
	            for (unsigned int k1=0;k1<Cindex.size();k1++) {
    	            if (Pindex.at(k1)==Cindex.at(n) && Rindex.at(k1)!=2) {
        	            //ctsisomer.at(0).at(k1)="";
            	        ctsisomer.at(1).at(k1)="";
                	}
            	}
        	}
        	if (Rindex.at(n)!=2 && Pindex.at(n)>1) {
				bool go=0;
                for (unsigned int k1=0;k1<Cindex.size();k1++) {
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
	chi.erase(chi.begin()+n);
	Cyindex.erase(Cyindex.begin()+n);	
	//ctsisomer.erase(ctsisomer.begin()+n);
	ctsisomer.at(0).erase(ctsisomer.at(0).begin()+n);
	ctsisomer.at(1).erase(ctsisomer.at(1).begin()+n);
	if (para.protect) protect.erase(protect.begin()+n);

	if (1) {
		if (Merge_C.size()>=2) {
			for (unsigned int k=1;k<Merge_C.size();k++) {
				for (unsigned int i=0;i<Cindex.size();i++) {
					if (Cindex.at(i)==Merge_C.at(k)) {
						Pindex.at(i)=2*Cindex.size()+3;
						break;
					}
				}
			}
		}

    	for (unsigned int i=0;i<Cindex.size();i++) {
        	//if (Cindex.at(i)!=(i+1)) {
				unsigned int j=Cindex.at(i);
				Cindex.at(i)=i+1;
				if (Merge_C.size()>=2) {
					for (unsigned int k=0;k<Merge_C.size();k++) {
						if (j==Merge_C.at(k)) {
							Merge_C.at(k)=Cindex.at(i);
							break;
						}
					}
				}
        		for (unsigned int k=0;k<Cindex.size();k++) { //k=i orig
        			if (Pindex.at(k)==j) Pindex.at(k)=Cindex.at(i);
        		}
        	//}
    	}


		if (Merge_C.size()>=2) {
			for (unsigned int k=1;k<Merge_C.size();k++) {
				Pindex.at(Merge_C.at(k)-1)=Merge_C.at(0);
				if (!Merge_C.at(0)) Rindex.at(Merge_C.at(k)-1)=0;//-1;
				else {
					if (bndfrm) Rindex.at(Merge_C.at(k)-1)=bndfrm;
				}
			}
		}
		
	}

	//vector<int>().swap(ref);
	//vector<int>().swap(pr);
	vector<unsigned int>().swap(Merge_C);

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
	//decyc_small_ring(5);

	//reset();
	//chk_cistrans(); //20200806
	if (!chk_imine_ct(0,Cindex.size()-1,cistrans)) mds2smi();
	//chk_imine_ct(0,Cindex.size()-1,cistrans);
	//mds2smi();

    //time1=time(NULL)-time1;
    //cout << "SUBTRACTION FINISHED: " << setprecision(10) << setw(15) << time1 << endl;

	return 1;
}

unsigned int MOLECULE::change_bnd(unsigned int n,unsigned int id1,unsigned int id2,unsigned int bnd,bool cistrans) {

    if (n<=0 || n>=Cindex.size()) return 0;
    unsigned int P,R;//M,P,C,R;
    //M=Mindex.at(n);
    P=Pindex.at(n);
    //C=Cindex.at(n);
    R=Rindex.at(n);

    if (para.protect) {
        if (protect.at(n)) return 0;
        if (P>0) if (protect.at(P-1)) return 0;
    }
    if (bnd>3 || bnd<0) return 0;
    if (data->a[id2].chg!=data->a[Mindex.at(n)].chg) return 0;



	if (P>0) {
		bool id1same=0,id2same=0;
		if (Mindex.at(P-1)==id1) id1same=1;
		if (data->a[id1].index<7 && data->a[Mindex.at(P-1)].index<7) {
			if (data->a[Mindex.at(P-1)].atm==data->a[id1].atm && data->a[id1].atm!="Group") id1same=1;
		}

		if (Mindex.at(n)==id2) id2same=1;
		if (data->a[id2].index<7 && data->a[Mindex.at(n)].index<7) {
			if (data->a[Mindex.at(n)].atm==data->a[id2].atm && data->a[id2].atm!="Group") id2same=1;
		}

		if (id1same && id2same && bnd==R) return 0;
	}

    //else if (data->a[id2].chg==data->a[Mindex.at(n)].chg && data->a[id2].chg!=0) {
    //    if (data->a[id2].index<=2 && data->a[Mindex.at(n)].index>=7) return 0;
    //    if (data->a[id2].index>=7 && data->a[Mindex.at(n)].index<=7) return 0;
    //}
	if (P>0) {
    	if (data->a[id1].chg!=data->a[Mindex.at(P-1)].chg) return 0;
    	//else if (data->a[id1].chg==data->a[Mindex.at(P-1)].chg && data->a[id1].chg!=0) {
        //	if (data->a[id1].index<=2 && data->a[Mindex.at(P-1)].index>=7) return 0;
        //	if (data->a[id1].index>=7 && data->a[Mindex.at(P-1)].index<=7) return 0;
    	//}
	}

    if (!data->a[id1].bd[bnd-1]) return 0;
    if (!data->a[id2].bd[bnd-1]) return 0;


    //reset();
    //vector<int> tmp;

    // Ensure the bonds of atom C(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
    if (1) {
        //bool bndchk=0;

        int id_bndsum=0; // max available # of bonds of element(id) (calc. from pool)
        id_bndsum+=bnd;
        vector<unsigned int> id_ordcount(6,0);
        id_ordcount.at(bnd-1)+=1;
        if (Cyindex.at(n).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(n).size();k1++) {
                id_bndsum+=Cybnd.at(Cyindex.at(n).at(k1)-1);
                id_ordcount.at(Cybnd.at(Cyindex.at(n).at(k1)-1)-1)++;
            }
        }
        for (unsigned int k1=0;k1<Cindex.size();k1++) {
            if (Pindex.at(k1)==Cindex.at(n)) {
                id_bndsum+=Rindex.at(k1);
                id_ordcount.at(Rindex.at(k1)-1)+=1;
            }
        }

        int id_maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            id_maxbnd+=data->a[id2].order.at(k2);
        }

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (id_ordcount.at(k2)>data->a[id2].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=0;
                return 0;
            }
            else if (k2>=2 && id_ordcount.at(k2)<=data->a[id2].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

    }
    // Ensure the bonds of atom P(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
    if (P>0) {
        //bool bndchk=0;

        int bndsum=0;  // # of bonds of atom P(n) if atom M(n) is changed to atom(id) and a $bnd2par bond of atom(id) is used to connected with P(n).
        if (Rindex.at(P-1)>0) bndsum+=Rindex.at(P-1);
        vector<unsigned int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
        if (Rindex.at(P-1)>0) ordcount.at(Rindex.at(P-1)-1)+=1;  //20200712 // bnd P(n) to P(P(n))
        if (Cyindex.at(P-1).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(P-1).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(P-1).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(P-1).at(k1)-1)-1)++;
            }
        }
        for (unsigned int k=0;k<Cindex.size();k++) {
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

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of P(n) after connecting with atom(id)" match the definition of atom P(n) in pool
            if (ordcount.at(k2)>data->a[id1].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=0;
                return 0;
            }
            else if (k2>=2 && ordcount.at(k2)<=data->a[id1].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

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
    //decyc_small_ring(5);

    //reset();
    //chk_cistrans(); //20200806
    if (!chk_imine_ct(0,Cindex.size()-1,cistrans)) mds2smi();
	//chk_imine_ct(0,Cindex.size()-1,cistrans);
	//mds2smi();

    return 1;

}

unsigned int MOLECULE::change_ele(unsigned int n,unsigned int id,unsigned int bnd2par,unsigned int bnd2des,bool cistrans) {
	
    if (n<0 || n>=Cindex.size()) return 0;
    unsigned int P,R;//M,P,C,R;
    //M=Mindex.at(n);
    P=Pindex.at(n);
    //C=Cindex.at(n);
    R=Rindex.at(n);

    if (para.protect) {
        if (protect.at(n)) return 0;
        //if (P>0) if (protect.at(P-1)) return 0;
    }
    //if (bnd2par>3 || bnd2par<=0) return 0;
    //if (bnd2des>3 || bnd2des<=0) return 0;
	if (bnd2par>3) return 0;
	if (bnd2des>3) return 0;
	if (bnd2par<=0 && Cindex.at(n)!=1) return 0;
    if (data->a[id].chg!=data->a[Mindex.at(n)].chg) return 0;
	if (!P && bnd2par>0) return 0;


	//else if (data->a[id].chg==data->a[Mindex.at(n)].chg && data->a[id].chg!=0) {
	//	if (data->a[id].index<=2 && data->a[Mindex.at(n)].index>=7) return 0;
	//	if (data->a[id].index>=7 && data->a[Mindex.at(n)].index<=7) return 0;
	//}

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

	vector<unsigned int> desbnds(0);
	bool isnotend=0;
	bool usepibnd=0;

	// Ensure the bonds of atom C(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
	if (1) {
        //bool bndchk=0;

	    int id_bndsum=0; // max available # of bonds of element(id) (calc. from pool)
		vector<unsigned int> id_ordcount(6,0);
		if (P>0 && bnd2par>0) {
			id_bndsum+=bnd2par;
			id_ordcount.at(bnd2par-1)+=1;
		}
		if (Cyindex.at(n).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(n).size();k1++) {
                id_bndsum+=Cybnd.at(Cyindex.at(n).at(k1)-1);
                id_ordcount.at(Cybnd.at(Cyindex.at(n).at(k1)-1)-1)++;
            }
		}
    	for (unsigned int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(n)) {
				desbnds.push_back(Rindex.at(k1));

				isnotend=1;
				if (bnd2des>0) {
					id_bndsum+=bnd2des;
					id_ordcount.at(bnd2des-1)+=1;
				}
				if (Rindex.at(k1)==2) usepibnd=1;
			}
		}
		if (isnotend && bnd2des<=0) return 0;

        int id_maxbnd=0; // max available # of bonds for atom C(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            id_maxbnd+=data->a[id].order.at(k2); 
        }

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with k" match the definition of atom k1 in pool
            if (id_ordcount.at(k2)>data->a[id].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=0;
				return 0;
            }
            else if (k2>=2 && id_ordcount.at(k2)<=data->a[id].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                //bndchk=1;
            }
        }

	}
	if (desbnds.size()) {
		for (unsigned int k1=0;k1<desbnds.size();k1++) {
			if (desbnds.at(k1)!=bnd2des) {
				vector<unsigned int>().swap(desbnds);
				break;
			}
			else if (k1>=desbnds.size()-1 && desbnds.at(k1)==bnd2des) {
				vector<unsigned int>().swap(desbnds);
				if (P>0) {
					if (data->a[Mindex.at(n)].atm==data->a[id].atm && bnd2par==R) return 0;
				}
				else {
					if (data->a[Mindex.at(n)].atm==data->a[id].atm) return 0;
				}
			}
		}
	}
	else {
		if (bnd2des>0) return 0;

        if (P>0) {
            if (data->a[Mindex.at(n)].atm==data->a[id].atm && bnd2par==R) return 0;
        }
        else {
            if (data->a[Mindex.at(n)].atm==data->a[id].atm) return 0;
        }
	}

	// Ensure the bonds of atom P(n) are OK if M(n) is changed to atom(id) and use a $bnd2par bond to connect P(n)
	if (P>0) {
        //bool bndchk=0;

        int bndsum=0;  // # of bonds of atom P(n) if atom M(n) is changed to atom(id) and a $bnd2par bond of atom(id) is used to connected with P(n).
		if (Rindex.at(P-1)>0) bndsum+=Rindex.at(P-1); 
        vector<unsigned int> ordcount(6,0); // count # of bonds atom k1 will use for the connection with atom k (in each bond order respectively)
        if (Rindex.at(P-1)>0) ordcount.at(Rindex.at(P-1)-1)+=1;  //20200712 // bnd P(n) to P(P(n))
        if (Cyindex.at(P-1).size()) {
            for (unsigned int k1=0;k1<Cyindex.at(P-1).size();k1++) {
                bndsum+=Cybnd.at(Cyindex.at(P-1).at(k1)-1);
                ordcount.at(Cybnd.at(Cyindex.at(P-1).at(k1)-1)-1)++;
            }
        }
        for (unsigned int k=0;k<Cindex.size();k++) {
            if (Pindex.at(k)==Cindex.at(P-1) && k!=n) {
                bndsum+=Rindex.at(k);
                ordcount.at(Rindex.at(k)-1)+=1;
            }
        }
		if (bnd2par>0) {
			bndsum+=bnd2par;
			ordcount.at(bnd2par-1)+=1;
		}

        int maxbnd=0; // max available # of bonds for atom P(n) (calc. from pool)
        for (int k2=0;k2<6;k2++) {
            maxbnd+=data->a[Mindex.at(P-1)].order.at(k2);
        }

        for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of P(n) after connecting with atom(id)" match the definition of atom P(n) in pool
        	if (ordcount.at(k2)>data->a[Mindex.at(P-1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
        		//bndchk=0;
				return 0;
        	}
        	else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(P-1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
            	//bndchk=1;
        	}
    	}

	}

    // Ensure the bonds of atom des(n) are OK if M(n) is changed to atom(id) and use a $bnd2des bond to connect des(n)
	if (isnotend) {
    	for (unsigned int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(n) && Mindex.at(k1)!=70) { // For the descendant atoms C(k1) of atom C(n)
            	//bool bndchk=0;

				int bndsum=0; // # of bonds of atom C(k1) if atom M(k1) is changed to atom(id) and use a $bnd2des bond to connect des(n).
            	vector<unsigned int> ordcount(6,0); // count # of bonds atom C(k1) will use after connecting with atom C(n) (in each bond order respectively)
				if (bnd2des>0) {
					bndsum+=bnd2des;
        			ordcount.at(bnd2des-1)+=1;  //20200712
				}
	        	if (Cyindex.at(k1).size()) {
		            for (unsigned int k2=0;k2<Cyindex.at(k1).size();k2++) {
                		bndsum+=Cybnd.at(Cyindex.at(k1).at(k2)-1);
                		ordcount.at(Cybnd.at(Cyindex.at(k1).at(k2)-1)-1)++;
            		}
        		}
            	for (unsigned int k=0;k<Cindex.size();k++) { // For each des(k1) atom
                    if (Pindex.at(k)==Cindex.at(k1)) {
                        bndsum+=Rindex.at(k); // calc. the total # of bonds des(k1) atoms used to connect with C(k1).
                        ordcount.at(Rindex.at(k)-1)+=1;
                    }
                }

                int maxbnd=0; // max available # of bonds for C(k1) (calc. from pool)
                for (int k2=0;k2<6;k2++) {
                    maxbnd+=data->a[Mindex.at(k1)].order.at(k2);
                }

                for (int k2=0;k2<3;k2++) {  // chk if "the bond orders of k1 after connecting with des(k1) atoms" match the definition of atom k1 in pool
                    if (ordcount.at(k2)>data->a[Mindex.at(k1)].bd[k2]) { // orig ordcount.at(k2)>data.a[Mindex.at(k1)].order.at(k2)
                        //bndchk=0;
						return 0;
                    }
                    else if (k2>=2 && ordcount.at(k2)<=data->a[Mindex.at(k1)].bd[k2]) { // orig k2>=5 && ordcount.at(k2)<=data.a[Mindex.at(k1)].order.at(k2)
                        //bndchk=1;
                    }
                }

            }
        }

    }


	Mindex.at(n)=id;
	if (P>0) Rindex.at(n)=bnd2par;
	else Rindex.at(n)=0;
	if (isnotend) {
    	for (unsigned int k1=0;k1<Cindex.size();k1++) {
        	if (Pindex.at(k1)==Cindex.at(n)) {
				Rindex.at(k1)=bnd2des;
			}
		}
	}

	if (1) {
		if (isnotend) {
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
					for (unsigned int k1=0;k1<Cindex.size();k1++) {
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
	}


    del_unpaired_ring_no(); //20200627
    //decyc_small_ring(5);

    //reset();
	//chk_cistrans(); //20200806
	if (!chk_imine_ct(0,Cindex.size()-1,cistrans)) mds2smi();
	//chk_imine_ct(0,Cindex.size()-1,cistrans);
	//mds2smi();

	return 1;

}

unsigned int MOLECULE::decyclization(unsigned int ringnum) {

	if (!ringnum) return 0;
	if (ringnum>if_circle) return 0;

    for (unsigned int i=0;i<Cyindex.size();i++) {
        if (Cyindex.at(i).size()) {
            for (int k3=0;k3<(int)Cyindex.at(i).size();k3++) {
                if (Cyindex.at(i).at(k3)==ringnum) {
					Cyindex.at(i).erase(Cyindex.at(i).begin()+k3);
					k3--;
                }
            }
        }
    }

    for (unsigned int i=0;i<Cyindex.size();i++) {
        if (Cyindex.at(i).size()) {
            for (unsigned int k3=0;k3<Cyindex.at(i).size();k3++) {
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

    del_unpaired_ring_no();
    //decyc_small_ring(5);

	if (!chk_imine_ct(0,Cindex.size()-1)) mds2smi();
	//chk_imine_ct(0,Cindex.size()-1);
	//mds2smi();

    return 1;

}


unsigned int MOLECULE::rechg() {
	unsigned int j=0;
	chg=0;
	for (j=0;j<Cindex.size();j++) {
		chg+=data->a[Mindex.at(j)].chg;
	}
	return 1;
}


unsigned int MOLECULE::printmds() {
	unsigned int i;
	//if (para.protect) prct();

	cout << setfill(' ');
	cout<<setw(8)<<left<<"SMILES"<<" "<<left<<molesmi<<endl;
	cout<<setw(8)<<left<<"orig_SMILES"<<" "<<left<<smiles<<endl;
	cout<<setw(8)<<left<<"natom"<<" "<<setw(4)<<left<<Cindex.size()<<endl;
	cout<<setw(8)<<left<<"Pindex"<<" ";
	for (i=0;i<Pindex.size();i++) cout<<setw(4)<<left<<Pindex.at(i)<<" ";
	cout<<endl;
	cout<<setw(8)<<left<<"Cindex"<<" ";
	for (i=0;i<Cindex.size();i++) cout<<setw(4)<<left<<Cindex.at(i)<<" ";
	cout<<endl;
	cout<<setw(8)<<left<<"Rindex"<<" ";
	for (i=0;i<Rindex.size();i++) cout<<setw(4)<<left<<Rindex.at(i)<<" ";
	cout<<endl;
	cout<<setw(8)<<left<<"Mindex"<<" ";
	for (i=0;i<Mindex.size();i++) cout<<setw(4)<<left<<Mindex.at(i)<<" ";
	cout<<endl;
	cout<<setw(8)<<left<<"Cyindex"<<" ";
	for (i=0;i<Cyindex.size();i++) {
        if (!Cyindex.at(i).size()) cout<<setw(4)<<left<<setfill(' ')<<"0"<<" ";
        else {
			stringstream buf("");
            for (unsigned int j=0;j<Cyindex.at(i).size();j++) {
                if (j<Cyindex.at(i).size()-1) buf<<Cyindex.at(i).at(j)<<",";
                else buf<<Cyindex.at(i).at(j);

                //if (j<Cyindex.at(i).size()-1) cout<<left<<Cyindex.at(i).at(j)<<",";
                //else cout<<left<<Cyindex.at(i).at(j)<<" ";
            }
			cout<<setw(4)<<left<<setfill(' ')<<buf.str()<<" ";
        }
	}
	cout<<endl;
    cout<<setw(8)<<left<<"Chirality"<<" ";
    for (i=0;i<chi.size();i++) cout<<setw(4)<<left<<chi.at(i)<<" ";
    cout<<endl;
    cout<<setw(8)<<left<<"ctsisomer_start"<<" ";
    for (i=0;i<ctsisomer.at(0).size();i++) {
        if (ctsisomer.at(0).at(i)!="") cout<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(0).at(i)<<" ";
        else cout<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    cout<<endl;
    cout<<setw(8)<<left<<"ctsisomer_end"<<" ";
    for (i=0;i<ctsisomer.at(1).size();i++) {
        if (ctsisomer.at(1).at(i)!="") cout<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(1).at(i)<<" ";
        else cout<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    cout<<endl;


    cout<<setw(8)<<left<<"if_circle"<<" "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;

    cout<<setw(8)<<left<<"Cybnd"<<" ";
    for (i=0;i<Cybnd.size();i++) cout<<setw(4)<<left<<Cybnd.at(i)<<" ";
	cout<<endl;

	return 1;
}		

unsigned int MOLECULE::printmds(ofstream &ouf) {
	unsigned int i;
	//if (para.protect) prct();

	ouf << setfill(' ');
	ouf<<setw(8)<<left<<"SMILES"<<" "<<left<<molesmi<<endl;
	ouf<<setw(8)<<left<<"orig_SMILES"<<" "<<left<<smiles<<endl;
	ouf<<setw(8)<<left<<"natom"<<" "<<setw(4)<<left<<Cindex.size()<<endl;
	ouf<<setw(8)<<left<<"Pindex"<<" ";
	for (i=0;i<Pindex.size();i++) ouf<<setw(4)<<left<<Pindex.at(i)<<" ";
	ouf<<endl;
	ouf<<setw(8)<<left<<"Cindex"<<" ";
	for (i=0;i<Cindex.size();i++) ouf<<setw(4)<<left<<Cindex.at(i)<<" ";
	ouf<<endl;
	ouf<<setw(8)<<left<<"Rindex"<<" ";
	for (i=0;i<Rindex.size();i++) ouf<<setw(4)<<left<<Rindex.at(i)<<" ";
	ouf<<endl;
	ouf<<setw(8)<<left<<"Mindex"<<" ";
	for (i=0;i<Mindex.size();i++) ouf<<setw(4)<<left<<Mindex.at(i)<<" ";
	ouf<<endl;
	ouf<<setw(8)<<left<<"Cyindex"<<" ";
	for (i=0;i<Cyindex.size();i++) {
        if (!Cyindex.at(i).size()) ouf<<setw(4)<<left<<setfill(' ')<<"0"<<" ";
        else {
			stringstream buf("");
            for (unsigned int j=0;j<Cyindex.at(i).size();j++) {
                if (j<Cyindex.at(i).size()-1) buf<<Cyindex.at(i).at(j)<<",";
                else buf<<Cyindex.at(i).at(j);

                //if (j<Cyindex.at(i).size()-1) ouf<<left<<Cyindex.at(i).at(j)<<",";
                //else ouf<<left<<Cyindex.at(i).at(j)<<" ";
            }
			ouf<<setw(4)<<left<<setfill(' ')<<buf.str()<<" ";
        }
	}
	ouf<<endl;
    ouf<<setw(8)<<left<<"Chirality"<<" ";
    for (i=0;i<chi.size();i++) ouf<<setw(4)<<left<<chi.at(i)<<" ";
    ouf<<endl;
	if (para.protect) {
		ouf<<setw(8)<<left<<"protect"<<" ";
		for (i=0;i<protect.size();i++) ouf<<setw(4)<<left<<protect.at(i)<<" ";
		ouf<<endl;
	}
    ouf<<setw(8)<<left<<"ctsisomer_start"<<" ";
    for (i=0;i<ctsisomer.at(0).size();i++) {
        if (ctsisomer.at(0).at(i)!="") ouf<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(0).at(i)<<" ";
        else ouf<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    ouf<<endl;
    ouf<<setw(8)<<left<<"ctsisomer_end"<<" ";
    for (i=0;i<ctsisomer.at(1).size();i++) {
        if (ctsisomer.at(1).at(i)!="") ouf<<setw(4)<<left<<setfill(' ')<<ctsisomer.at(1).at(i)<<" ";
        else ouf<<setw(4)<<left<<setfill(' ')<<"_"<<" ";
    }
    ouf<<endl;


    ouf<<setw(8)<<left<<"if_circle"<<" "<<setw(4)<<left<<setfill(' ')<<if_circle<<endl;

    ouf<<setw(8)<<left<<"Cybnd"<<" ";
    for (i=0;i<Cybnd.size();i++) ouf<<setw(4)<<left<<Cybnd.at(i)<<" ";
	ouf<<endl;


	return 1;
}

unsigned int MOLECULE::prct() {
	prct(0,Cindex.size());
	return 1;
}

unsigned int MOLECULE::prct(unsigned int a,unsigned int b) {

	//long double time1=time(NULL);

	if (a>Cindex.size() || b>Cindex.size() || a<0 || b<0 || Cindex.size()==0) return 0;
	unsigned int i=0;//,j=0;
	int count=0;
	if (para.protect) {
		for (i=0;i<protect.size();i++) {
			if (protect.at(i)) count++;
		}
	}

	//bool go=1;
	//if (count/(double)Cindex.size()>0.5) go=0;
	if (para.protect) { // para.protect && go
		if (0) { // double bonds
			protect.resize(0);
			protect.resize(Cindex.size(),0); //20190801

			for (i=a;i<b;i++) {
				if (Rindex.at(i)==2) { //e.g. CC=CC
					protect.at(i)=1; // i = C3
					for (unsigned int k1=a;k1<b;k1++) {
						if (Pindex.at(i)==Cindex.at(k1)) {
							protect.at(k1)=1; // k1 = C2
							for (unsigned int k2=a;k2<b;k2++) {
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
				if (data->a.at(Mindex.at(i)).chg) {
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
					unsigned int rep_quota=1;
					for (int k2=Cyindex.size()-1;k2>=0;k2--) {
						if (Cyindex.at(k2).size()) {
							for (unsigned int k3=0;k3<Cyindex.at(k2).size();k3++) {
								if ((int)Cyindex.at(k2).at(k3)==k1) {
									unsigned int compa=Cindex.at(k2);
									unsigned int cont=0;
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
						for (unsigned int k2=0;k2<C_ringmember.size()-1;k2++) {
							for (unsigned int k3=k2+1;k3<C_ringmember.size();k3++) {
								if (C_ringmember.at(k2)==C_ringmember.at(k3) && rep_quota<=0) {
									C_ringmember.at(k2)=C_ringmember.at(k3)=-1;
								}
								else if (C_ringmember.at(k2)==C_ringmember.at(k3) && rep_quota>0) {
									C_ringmember.at(k3)=-1;
									rep_quota--;
								}
							}
						}
						for (unsigned int k2=0;k2<C_ringmember.size();k2++) {
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


unsigned int MOLECULE::decyc_small_ring(unsigned int size) {

    vector<int> C_ringmember(0);
	
    for (int k1=if_circle;k1>=1;k1--) { // k1>=ct+1
        C_ringmember.resize(0);
        int rep_quota=1;
        for (int k2=Cyindex.size()-1;k2>=0;k2--) {
            if (Cyindex.at(k2).size()) {
                for (unsigned int k3=0;k3<Cyindex.at(k2).size();k3++) {
                    if ((int)Cyindex.at(k2).at(k3)==k1) {
                    	unsigned int compa=Cindex.at(k2);
                        unsigned int cont=0;
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
            for (unsigned int k2=0;k2<C_ringmember.size()-1;k2++) {
                for (unsigned int k3=k2+1;k3<C_ringmember.size();k3++) {
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
			unsigned int num_rmember=0;
        	for (unsigned int k2=0;k2<C_ringmember.size();k2++) {
            	if (C_ringmember.at(k2)!=-1) num_rmember++;
        	}
			if (num_rmember<size) {
				for (unsigned int k2=0;k2<Cyindex.size();k2++) {
					if (Cyindex.at(k2).size()) {
                		for (int k3=0;k3<(int)Cyindex.at(k2).size();k3++) {
                    		if ((int)Cyindex.at(k2).at(k3)==k1) {
								Cyindex.at(k2).erase(Cyindex.at(k2).begin()+k3);
								k3--;
							}
						}
					}
				}

                for (unsigned int k2=0;k2<Cyindex.size();k2++) {
                    if (Cyindex.at(k2).size()) {
                        for (unsigned int k3=0;k3<Cyindex.at(k2).size();k3++) {
                            if ((int)Cyindex.at(k2).at(k3)>k1) {
                                Cyindex.at(k2).at(k3)--;
                            }
                        }
                    }
                }

				if ((int)Cybnd.size()>=k1) Cybnd.erase(Cybnd.begin()+k1-1);
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


unsigned int MOLECULE::del_unpaired_ring_no(){
	int k;
    vector<unsigned int> freq(if_circle,0);
	if (if_circle) {
	    for (k=0;k<(int)Cyindex.size();k++) {  // 20200131 prevent unpaired ring number
	        if (Cyindex.at(k).size()) {
	            for (unsigned int k1=0;k1<Cyindex.at(k).size();k1++) {
					freq.at(Cyindex.at(k).at(k1)-1)++;
	            }
	        }
	    }

	    for (k=freq.size()-1;k>=0;k--) { //k>=ct
			if (!freq.at(k) || freq.at(k)%2) {

                for (unsigned int k1=0;k1<Cyindex.size();k1++) {
                    if (Cyindex.at(k1).size()) {
                        for (int k2=0;k2<(int)Cyindex.at(k1).size();k2++) {
                            if ((int)Cyindex.at(k1).at(k2)==k+1) {
								Cyindex.at(k1).erase(Cyindex.at(k1).begin()+k2);
								k2--;
							}
                        }
                    }
                }

                for (unsigned int k1=0;k1<Cyindex.size();k1++) {
                    if (Cyindex.at(k1).size()) {
                        for (int k2=0;k2<(int)Cyindex.at(k1).size();k2++) {
							if ((int)Cyindex.at(k1).at(k2)>k+1) {
								Cyindex.at(k1).at(k2)--;
							}
                        }
                    }
                }

				if ((int)Cybnd.size()>k) Cybnd.erase(Cybnd.begin()+k);

				if_circle--;

			}
	    }
	}
	vector<unsigned int>().swap(freq);

	return 1;
}


unsigned int MOLECULE::canonicalize_SMILES() {
	stringstream ss("");
    if (1) {
        stringstream ss1(molesmi);
        OBConversion conv(&ss1, &ss);
        if(conv.SetInAndOutFormats("SMI","SMI")) {
            //conv.AddOption("canonical", OBConversion::GENOPTIONS);
			conv.AddOption("h", OBConversion::OUTOPTIONS);
			//conv.AddOption("S", OBConversion::INOPTIONS);
			//conv.AddOption("S", OBConversion::INOPTIONS);
            //conv.AddOption("C", OBConversion::OUTOPTIONS);
            //conv.AddOption("gen3D", OBConversion::GENOPTIONS);
            conv.Convert();
        }
    }
	if (1) {
    	stringstream ss_out("");
    	OBConversion conv(&ss, &ss_out);
    	if(conv.SetInAndOutFormats("SMI","SMI"))
    	{
        	//conv.AddOption("h", OBConversion::GENOPTIONS);
        	//conv.AddOption("gen3D", OBConversion::GENOPTIONS);
			//conv.AddOption("gen2D", OBConversion::GENOPTIONS);
			//conv.AddOption("S", OBConversion::INOPTIONS);
			//conv.AddOption("k", OBConversion::OUTOPTIONS);
        	conv.AddOption("canonical", OBConversion::GENOPTIONS);
            obErrorLog.SetOutputLevel(obMessageLevel::obError);
            //obErrorLog.SetOutputLevel(obMessageLevel::obWarning);
            obErrorLog.SetOutputStream(&cout);
        	conv.Convert();
    	}
		ss_out >> molesmi >> ws;
		//smiles=molesmi;
	}

	return 1;
}

unsigned int MOLECULE::chk_imine_ct(unsigned int lpos,unsigned int hpos,bool cistrans) {
	bool change1=0,change2=0;
    for (int i=lpos;i<=(int)hpos;i++) {
        if (Mindex.at(i)==70) {
			bool do1=1;
			for (int j=i+1;j<(int)Cindex.size();j++) {
				if (Pindex.at(j)==Cindex.at(i) && Mindex.at(j)!=8) {
					if (subtraction(i)) {
						i--;
						j--;
						hpos=Cindex.size()-1;
						change1=1;
						do1=0;
						break;
					}
				}
				else if (Pindex.at(j)==Cindex.at(i) && Mindex.at(j)==8) {
					bool del=1;
					for (int g=j+1;g<(int)Cindex.size();g++) {
						if (Pindex.at(g)==Cindex.at(j) && Rindex.at(g)==2) {
							del=0;
							break;
						}
					}
					if (del) {
						if (subtraction(i)) {
							i--;
							j--;
							hpos=Cindex.size()-1;
							change1=1;
							do1=0;
							break;
						}
					}
				}
			}

            if (Pindex.at(i)>0 && do1) {
                if (Mindex.at(Pindex.at(i)-1)!=8) {
                    if (subtraction(i)) {
						i--;
						hpos=Cindex.size()-1;
						change1=1;
					}
                }
				else if (Mindex.at(Pindex.at(i)-1)==8 && Rindex.at(Pindex.at(i)-1)!=2) {
					bool del=1;
					for (int g=0;g<(int)Cindex.size();g++) {
						if (Pindex.at(i)==Pindex.at(g) && g!=i && Rindex.at(g)==2) {
							del=0;
							break;
						}
					}
                    if (del) {
                        if (subtraction(i)) {
                            i--;
                            hpos=Cindex.size()-1;
                            change1=1;
                        }
                    }
				}
            }
        }
    }
	for (int i=lpos;i<=(int)hpos;i++) {
		if (Mindex.at(i)==8) {
			if (Rindex.at(i)==2) {
				bool if_addDu=1;
				for (int j=i+1;j<(int)hpos;j++) {
					if (Pindex.at(j)==Cindex.at(i)) {
						if_addDu=0;
						break;
					}
				}
                if (Pindex.at(i)>0 && if_addDu) {
                    if (Rindex.at(Pindex.at(i)-1)!=1) if_addDu=0;
                    if (if_addDu) {
                        for (int j=Pindex.at(i);j<(int)hpos;j++) {
                            if (Pindex.at(j)==Pindex.at(i) && j!=i && Rindex.at(j)!=1) {
                                if_addDu=0;
                                break;
                            }
                        }
                    }
                }
				if (if_addDu) {
					if (addition(i,70,1,cistrans)) {
						change2=1;
						hpos=Cindex.size()-1;
					}
				}
			}
			//else if (Rindex.at(i)==0) {
			else {
				if (0) {
					bool usepibnd=0,hassigbnd=0;
					if (Rindex.at(i)==1) hassigbnd=1;
					for (int j=i+1;j<(int)hpos;j++) {
						if (Pindex.at(j)==Cindex.at(i)) {
							if (Rindex.at(j)==2) usepibnd=1;
							else if (Rindex.at(j)==1) hassigbnd=1;
						}
					}

					if (usepibnd && !hassigbnd) {
						if (addition(i,70,1,cistrans)) {
							change2=1;
							hpos=Cindex.size()-1;
						}
					}
				}
                if (1) {
                    bool do2=0;
                    if (Rindex.at(i)==1) {
                        for (int j=i+1;j<(int)hpos;j++) {
                            if (Pindex.at(j)==Cindex.at(i) && Rindex.at(j)==2) {
                                do2=1;
                                break;
                            }
                        }
                    }
                    if (do2) {
                        if (addition(i,70,1,cistrans)) {
                            change2=1;
                            hpos=Cindex.size()-1;
                        }
                    }
                }


				if (hpos>2) {
					if (ctsisomer.at(0).at(0)!="" && ctsisomer.at(1).at(1)!="" && Rindex.at(1)==2 && Pindex.at(1)==1) {
						unsigned int ct1=0,ct2=0;
						for (int k1=lpos;k1<(int)hpos;k1++) {
							if (Pindex.at(k1)==Cindex.at(0) && Cindex.at(k1)!=Cindex.at(1)) {
								if (ct1==0) {
									if (ctsisomer.at(0).at(0)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
									else if (ctsisomer.at(0).at(0)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
								}
								//else if (ct1==1) {
								//  ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(0).at(0)+ctsisomer.at(0).at(Cindex.at(k1)-1);
								//}
								ct1++;
							}
							if (Pindex.at(k1)==Cindex.at(1)) {
								if (ct2==0) {
									ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(1).at(1)+ctsisomer.at(0).at(Cindex.at(k1)-1);
								}
								//if (ct2==0) {
								//  if (ctsisomer.at(1).at(1)=="/") ctsisomer.at(0).at(Cindex.at(k1)-1)="\\"+ctsisomer.at(0).at(Cindex.at(k1)-1);
								//  else if (ctsisomer.at(1).at(1)=="\\") ctsisomer.at(0).at(Cindex.at(k1)-1)="/"+ctsisomer.at(0).at(Cindex.at(k1)-1);
								//}
								//else if (ct2==1) {
								//  ctsisomer.at(0).at(Cindex.at(k1)-1)=ctsisomer.at(1).at(1)+ctsisomer.at(0).at(Cindex.at(k1)-1);
								//}
								ct2++;
							}
						}
						ctsisomer.at(0).at(0)="";
						ctsisomer.at(1).at(1)="";
					}
				}
			}
		}
	}
	if (!change1 && !change2) return 0;
	else return 1;
}


unsigned int MOLECULE::readmds(string a) {
    //for (unsigned int i=0;i<a.length();i++) {
    //    if (a[i]=='/') a[i]='u';
    //    if (a[i]=='\\') a[i]='d';
    //    if (a[i]=='*') a[i]='x';
    //}

    ifstream inf(a.c_str());
    string b,b1;
    unsigned int i=0,j=0,k=0;
    inf >> ws;

	Pindex.resize(0);
	Cindex.resize(0);
	Mindex.resize(0);
	Rindex.resize(0);
	Cyindex.resize(0,vector<unsigned int>(0));
	protect.resize(0);
	chi.resize(0);
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
            Cyindex.resize(j,vector<unsigned int>(0));
            for (i=0;i<j;i++) {
                b1="";
                inf>>b1>>ws;

                string buf="";
                for (unsigned int k1=0;k1<b1.length();k1++) {
                    if (b1[k1]!=',') buf+=b1[k1];
                    if (b1[k1]==',') {
                        if (buf!="") {
                            Cyindex.at(i).push_back(atoi(buf.c_str()));
                        }
                        buf="";
                    }
                    if (b1[k1]!=',' && k1>=b1.length()-1) {
                        if (buf!="0") {
                            Cyindex.at(i).push_back(atoi(buf.c_str()));
                        }
                        buf="";
                    }
                }
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
            for (i=0;i<if_circle;i++) {
                inf>>k>>ws;
                Cybnd.push_back(k);
            }
        }
        else if (b=="Chirality") {
            for (i=0;i<j;i++) {
                inf>>k>>ws;
                chi.push_back(k);
            }
        }
        else if (b=="SMILES") {
        	b1="";
            inf>>b1>>ws;
        	molesmi=smiles=b1;
        }
    }
    inf.close();
    return 1;
}

