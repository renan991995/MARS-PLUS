#include "CASES_NEU.h"
using namespace std;

extern PARAMETER para;

void exhaustive_combination(MOLECULE &A, MOLECULE &B) {
	int count=0,i=0,j=0;
	string smi1=A.molesmi,smi2=B.molesmi;
	
	MOLECULE *buf=new MOLECULE [2];
	buf[0].replace(A);
	buf[1].replace(B);
	
	ofstream out("combination.txt");
	mark(out);
	for (i=0;i<A.Cindex.size();i++) {
    	for (j=0;j<B.Cindex.size();j++) {
			for (int b=1;b<4;b++) {
        	if (A.combine(B,i,j,b)) {
        		count++;
        		A.mds2smi();
            	A.smiles=A.molesmi;   
				         	
				A.input();
				
				out<<A.molesmi<<endl;
				cout << "EXHAUSTIVE COMBINATION: " << smi1 << " {"<<(i)<<"} + " << smi2 << " {"<<(j)<<"} , bnd " << b << " -> " << A.molesmi << " {No. " << count << "}" << endl;
				ofstream outs(("./mds/comb-"+A.molesmi+".enc").c_str());						
				A.print(outs);
				outs.close();
        
				if (0) {
					A.smiles=smi1;
					A.input();
				}
				if (1) A.replace(buf[0]);

				if (0) {
                	B.smiles=smi2;
                	B.input();
				}
				if (1) B.replace(buf[1]);
				
        	}
			}
    	}
	}

	out.close();

	delete [] buf;
	buf=NULL;
	
}

void exhaustive_crossover(MOLECULE &A, MOLECULE &B) {
	int count=0,i=0,j=0;
	A.canonicalize_SMILES();
	string smi1=A.molesmi,smi2=B.molesmi;
	
	MOLECULE *buf=new MOLECULE [2];
	buf[0].replace(A);	
	buf[1].replace(B);
	
	ofstream out("crossover.txt");
	mark(out);
	for (i=1;i<A.Cindex.size();i++) {
		for (j=1;j<B.Cindex.size();j++) {
			if (A.crossover(B,i,j)) {
				count++;
				A.mds2smi();
				A.smiles=A.molesmi;  
				A.input();
				
				B.mds2smi();
				B.smiles=B.molesmi;  
				B.input();
				
				out<<A.molesmi<<endl;
				ofstream outs(("./mds/cros-"+A.molesmi+".enc").c_str());
				A.print(outs);
				outs.close();
				
				out<<B.molesmi<<endl;
				outs.open(("./mds/cros-"+B.molesmi+".enc").c_str());
				B.print(outs);
				outs.close();

				cout << "EXHAUSTIVE CROSSOVER: " << smi1 << " {"<<(i)<<"} + " << smi2 << " {"<<(j)<<"} -> " << A.molesmi << " + " << B.molesmi << " {No. " << count << "}" << endl;

				if (0) {
					A.smiles=smi1;
					A.input();
				}
				if (1) A.replace(buf[0]);

				if (0) {
					B.smiles=smi2;
					B.input();
				}
				if (1) B.replace(buf[1]);

			}
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;	
	
}


void exhaustive_subtraction(MOLECULE &A) {
	int count=0,i=0,j=0;
	string smi=A.molesmi;
	
	MOLECULE *buf=new MOLECULE [1];
	buf[0].replace(A);
	
	ofstream out("subtraction.txt");
	mark(out);
	for (i=2;i<A.Cindex.size();i++) {
		if (A.subtract(i,1)) {
			count++;
			A.mds2smi();
			A.smiles=A.molesmi;  
			A.input();
			
			out<<A.molesmi<<endl;
			cout << "EXHAUSTIVE SUBTRACTION: " << smi << " {"<<(i)<<"} -> " << A.molesmi << " {No. " << count << "}" << endl;
	    	ofstream outs(("./mds/subtr-"+A.molesmi+".enc").c_str());
			A.print(outs);
			outs.close();
			
			if (0) {
				A.smiles=smi;
				A.input();
			}
			if (1) A.replace(buf[0]);
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;

}


void exhaustive_addition(MOLECULE &A) {
	int count=0,i=0,j=0;
	string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);

	
	ofstream out("addition.txt");
	mark(out);
	for (i=0;i<A.Cindex.size();i++) {
		for (j=1;j<A.data->num;j++) {
			if (A.data->a[j].probability>0) {
				for (int b=1;b<4;b++) {
				if (A.add(i,j,b)) {
					count++;
					A.mds2smi();
					A.smiles=A.molesmi;  
					A.input();
				
					out<<A.molesmi<<endl;
					cout << "EXHAUSTIVE ADDITION: " << smi << " {"<< j << " " << A.data->a[j].name << " on " <<(i)<<"} -> " << A.molesmi << " {No. " << count << "}" << endl;
					ofstream outs(("./mds/add-"+A.molesmi+".enc").c_str());
					A.print(outs);
					outs.close();
					
					if (0) {
						A.smiles=smi;
						A.input();
					}
					if (1) A.replace(buf[0]);
				}
				}
			}
		}
	}	
	
	out.close();

	delete [] buf;
	buf=NULL;

}

void exhaustive_exchange(MOLECULE &A) {
	int count=0,i=0,j=0,k=0,m=0;
	string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);
	
	ofstream out("exchange.txt");
	mark(out);
	for (i=1;i<A.Cindex.size();i++) {
		for (j=1;j<A.data->num;j++) {
			if (A.data->a[j].probability>0) {
				for (k=-2;k<3;k++) {
					for (m=1;m<A.data->num;m++) {
						if (A.data->a[m].probability>0) {
							if (A.exchange(i,j,m,k)) {
								count++;
								A.mds2smi();
								A.smiles=A.molesmi;  
								A.input();
						
								out<<A.molesmi<<endl;
								cout << "EXHAUSTIVE EXCHANGE: " << smi << " {id1: " << j << " " << A.data->a[j].name << " | id2: " << m << " " << A.data->a[m].name <<" | on " << i << "} -> " << A.molesmi << " {No. " << count << "}" << endl;
								ofstream outs(("./mds/exc-"+A.molesmi+".enc").c_str());
								A.print(outs);
								outs.close();
								
								if (0) {
									A.smiles=smi;
									A.input();
								}
								if (1) A.replace(buf[0]);
							}
						}
					}
				}
			}
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;

}

void exhaustive_cyclization(MOLECULE &A) {
	int count=0,i=0,j=0;
	string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);
	
	ofstream out("cyclization.txt");
	mark(out);
	for (i=0;i<A.Cindex.size();i++) {
		for (j=i+1;j<A.Cindex.size();j++) {
			if (A.ring(i,j)) {
				count++;
				A.mds2smi();
				A.smiles=A.molesmi;  
				A.input();
				
				out<<A.molesmi<<endl;
				cout << "EXHAUSTIVE CYCLIZATYION: " << smi << " {"<< j << " & " << i <<"} -> " << A.molesmi << " {No. " << count << "}" << endl;
				ofstream outs(("./mds/cyc-"+A.molesmi+".enc").c_str());
				A.print(outs);
				outs.close();

				if (0) {
					A.smiles=smi;
					A.input();
				}
				if (1) A.replace(buf[0]);
			}
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;
}

void CH4_to_Bz(MOLECULE &A) {

	ofstream out("benzene.txt");
	mark(out);	
	
	cout<<"Make benzene, c1ccccc1 from methane, C"<<endl;
	A.smiles="C";
	A.input();
	out<<A.molesmi<<endl;
	
	A.add(0,1,1);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.input();
	out<<A.molesmi<<endl;
	
	A.exchange(1,2,2,1);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.input();
	out<<A.molesmi<<endl;
	
	A.add(1,2,1);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.input();
	out<<A.molesmi<<endl;
	
	A.add(2,2,2);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.input();
	out<<A.molesmi<<endl;
	
	A.add(3,2,1);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.input();
	out<<A.molesmi<<endl;
	
	A.add(4,2,2);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.input();
	out<<A.molesmi<<endl;
	
	A.ring(0,5);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.input();
	out<<A.molesmi<<endl;
	
	out.close();
}


