#include "CASES.h"
using namespace std;

extern PARAMETER para;

void exhaustive_combination(MOLECULE &A, MOLECULE &B) {
	int count=0,i=0,j=0;
	string smi=A.molesmi;
	
	//MOLECULE *buf=new MOLECULE [1];
	//buf[0].replace(A);
	
	ofstream out("combination.txt");
	mark(out);
	for (i=0;i<A.Cindex.size();i++) {
    	for (j=0;j<B.Cindex.size();j++) {
        	if (A.combine(B,i,j)) {
        		count++;
        		A.mds2smi();
            	A.smiles=A.molesmi;   
				         	
				A.empty();
            	A.clean();
				A.input();
				A.mds2smi();
				if (para.protect) A.prct();
				
				A.canonicalize_SMILES();
				out<<A.molesmi<<endl;
				cout << "EXHAUSTIVE COMBINATION: " << smi << " ("<<(i+1)<<") + " << B.molesmi << " ("<<(j+1)<<") -> " << A.molesmi << "(No. " << count << ")" << endl;
				ofstream outs(("./mds/comb-"+A.molesmi+".enc").c_str());						
				A.print(outs);
				outs << "-----" << endl;
				outs.close();
        
				A.smiles=smi;
				A.empty();
				A.clean();
				A.input();
				A.mds2smi();
				if (para.protect) A.prct();	
				
				//A.replace(buf);
        	}
    	}
	}

	//delete [] buf;
	//buf=NULL;
	
	out.close();
}

void exhaustive_crossover(MOLECULE &A, MOLECULE &B) {
	int count=0,i=0,j=0;
	string smi1=A.molesmi,smi2=B.molesmi;
	
	//MOLECULE *buf=new MOLECULE [2];
	//buf[0].replace(A);	
	//buf[1].replace(B);
	
	ofstream out("crossover.txt");
	mark(out);
	for (i=1;i<A.Cindex.size();i++) {
		for (j=1;j<B.Cindex.size();j++) {
			if (A.crossover(B,i,j)) {
				A.mds2smi();
				A.smiles=A.molesmi;  
				A.empty();
            	A.clean();
				A.input();
				A.mds2smi();
				if (para.protect) A.prct();
				
				B.mds2smi();
				B.smiles=B.molesmi;  
				B.empty();
            	B.clean();
				B.input();
				B.mds2smi();
				if (para.protect) B.prct();
				
				A.canonicalize_SMILES();
				out<<A.molesmi<<endl;
				ofstream outs(("./mds/cros-"+A.molesmi+".enc").c_str());
				A.print(outs);
				outs.close();
				//A.mds23d(out);
				
				B.canonicalize_SMILES();
				out<<B.molesmi<<endl;
				outs.open(("./mds/cros-"+B.molesmi+".enc").c_str());
				B.print(outs);
				outs.close();
				//B.mds23d(out);

				A.smiles=smi1;
				A.empty();
				A.clean();
				A.input();
				A.mds2smi();
				if (para.protect) A.prct();	

				B.smiles=smi2;
				B.empty();
				B.clean();
				B.input();
				B.mds2smi();
				if (para.protect) B.prct();	

				//A.replace(buf[0]);
				//B.replace(buf[1]);
			}
		}
	}
	
	//delete [] buf;
	//buf=NULL;	
	
	out.close();
}


void exhaustive_subtraction(MOLECULE &A) {
	int i=0,j=0;
	string smi=A.molesmi;
	
	
	ofstream out("subtraction.txt");
	mark(out);
	for (i=1;i<A.Cindex.size();i++) {
		if (A.subtract(i,1)) {
			A.mds2smi();
			A.smiles=A.molesmi;  
			A.empty();
	    	A.clean();
			A.input();
			A.mds2smi();
			if (para.protect) A.prct();
			
			A.canonicalize_SMILES();
			out<<A.molesmi<<endl;
	    	ofstream outs(("./mds/subtr-"+A.molesmi+".enc").c_str());
			A.print(outs);
			outs.close();
	
			//A.mds23d(out);
			
			A.smiles=smi;
			A.empty();
			A.clean();
			A.input();
			A.mds2smi();
			if (para.protect) A.prct();	
		}
	}
	
	out.close();
}


void exhaustive_addition(MOLECULE &A) {
	int i=0,j=0;
	string smi=A.molesmi;
	
	ofstream out("addition.txt");
	mark(out);
	for (i=0;i<A.Cindex.size();i++) {
		for (j=1;j<66;j++) {
			if (A.add(i,j)) {
				A.mds2smi();
				A.smiles=A.molesmi;  
				A.empty();
		    	A.clean();
				A.input();
				A.mds2smi();
				if (para.protect) A.prct();
				
				A.canonicalize_SMILES();
				out<<A.molesmi<<endl;
				ofstream outs(("./mds/add-"+A.molesmi+".enc").c_str());
				A.print(outs);
				outs.close();
				
				//A.mds23d(out);
				
				A.smiles=smi;
				A.empty();
				A.clean();
				A.input();
				A.mds2smi();
				if (para.protect) A.prct();	
			}
		}
	}	
	
	out.close();
}

void exhaustive_exchange(MOLECULE &A) {
	int i=0,j=0,k=0,m=0;
	string smi=A.molesmi;
	
	ofstream out("exchange.txt");
	mark(out);
	for (i=1;i<A.Cindex.size();i++) {
		for (j=1;j<66;j++) {
			for (k=-2;k<3;k++) {
				for (m=1;m<66;m++) {
					if (A.exchange(i,j,m,k)) {
						A.mds2smi();
						A.smiles=A.molesmi;  
						A.empty();
				    	A.clean();
						A.input();
						A.mds2smi();
						if (para.protect) A.prct();
						
						A.canonicalize_SMILES();
						out<<A.molesmi<<endl;
						ofstream outs(("./mds/exc-"+A.molesmi+".enc").c_str());
						A.print(outs);
						outs.close();
						
						//A.mds23d(out);
						
						A.smiles=smi;
						A.empty();
						A.clean();
						A.input();
						A.mds2smi();
						if (para.protect) A.prct();	
					}
				}
			}
		}
	}
	
	out.close();
}

void exhaustive_cyclization(MOLECULE &A) {
	int i=0,j=0;
	string smi=A.molesmi;
	
	ofstream out("cyclization.txt");
	mark(out);
	for (i=0;i<A.Cindex.size();i++) {
		for (j=i+1;j<A.Cindex.size();j++) {
			if (A.ring(i,j)) {
				A.mds2smi();
				A.smiles=A.molesmi;  
				A.empty();
		    	A.clean();
				A.input();
				A.mds2smi();
				if (para.protect) A.prct();				
				
				A.canonicalize_SMILES();
				out<<A.molesmi<<endl;
				ofstream outs(("./mds/cyc-"+A.molesmi+".enc").c_str());
				A.print(outs);
				outs.close();
				
				//A.mds23d(out);

				A.smiles=smi;
				A.empty();
				A.clean();
				A.input();
				A.mds2smi();
				if (para.protect) A.prct();						
			}
		}
	}
	
	out.close();
}

void CH4_to_Bz(MOLECULE &A) {

	ofstream out("benzene.txt");
	mark(out);	
	
	cout<<"Make benzene, c1ccccc1 from methane, C"<<endl;
	A.smiles="C";
	A.empty();
	A.clean();
	A.input();
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<A.molesmi<<endl;
	//gs[0].mds23d(out);
	
	A.add(0,1);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.empty();
	A.clean();
	A.input();
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<A.molesmi<<endl;
	//gs[0].mds23d(out);
	
	A.exchange(1,2,2,1);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.empty();
	A.clean();
	A.input();
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<A.molesmi<<endl;
	//gs[0].mds23d(out);
	
	A.add(1,2);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.empty();
	A.clean();
	A.input();
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<A.molesmi<<endl;
	//gs[0].mds23d(out);
	
	A.add(2,2);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.empty();
	A.clean();
	A.input();
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<A.molesmi<<endl;
	//gs[0].mds23d(out);
	
	A.add(3,2);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.empty();
	A.clean();
	A.input();
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<A.molesmi<<endl;
	//gs[0].mds23d(out);
	
	A.add(4,2);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.empty();
	A.clean();
	A.input();
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<A.molesmi<<endl;
	//gs[0].mds23d(out);
	
	A.ring(0,5);
	A.mds2smi();
	A.smiles=A.molesmi;  
	A.empty();
	A.clean();
	A.input();
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<A.molesmi<<endl;
	//gs[0].mds23d(out);
	
	out.close();
}


