#include "CASES_IL_INDEPENDENT.h"
using namespace std;

extern PARAMETER para;

unsigned int exhaustive_combination1(IL &A, IL &B) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;
	string smi2_cat=B.ion[0].molesmi,smi2_an=B.ion[1].molesmi;
	
	IL *buf=new IL [2];
	buf[0].replace(A);
	buf[1].replace(B);
	
	ofstream out((para.logdir+"combination_IL.txt").c_str(),ios::app);
	//mark(out);
	
	out << "IL1: " << (smi1_cat+"."+smi1_an)  << endl;
	A.printmds(out);
	out << endl;
	out << "IL2: " << (smi2_cat+"."+smi2_an)  << endl;
	B.printmds(out);
	out << endl;
		
	for (unsigned int s=0;s<2;s++) {
		if (s==0) {
			count=0;
			out << setw(7) << left << "No." << " "
				<< setw(65) << left << "IL_cat_comb" << " "
				<< setw(17) << left << "IL1_cat_comb_pt" << " "
	        	<< setw(17) << left << "IL2_cat_comb_pt" << " "
				<< setw(17) << left << "bnd" << endl;
	    }    	
	    if (s==1) {
			count=0;
			out << setw(7) << left << "No." << " "
				<< setw(65) << left << "IL_an_comb" << " "
				<< setw(17) << left << "IL1_an_comb_pt" << " "
	        	<< setw(17) << left << "IL2_an_comb_pt" << " "
				<< setw(17) << left << "bnd" << endl;
	    }    	

		for (unsigned int i=0;i<A.ion[s].Cindex.size();i++) {
			for (unsigned int j=0;j<B.ion[s].Cindex.size();j++) {
				unsigned int ACi=A.ion[s].Cindex.at(i);
				unsigned int BCj=B.ion[s].Cindex.at(j);

				for (unsigned int b=1;b<4;b++) {
	        		if (A.ion[s].combination(B.ion[s],i,j,b)) {
					
	        			//A.ion[s].mds2smi();
						if (0) {
							A.ion[s].smiles=A.ion[s].molesmi;
							A.ion[s].input();
						}
						if (1) A.ion[s].canonicalize_SMILES();

						A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
						//out<<A.molesmi<<endl;

                        string mm=A.ion[s].molesmi;
                        for (unsigned int q=0;q<mm.length();q++) {
                            if (mm[q]=='/') mm[q]='u';
                            if (mm[q]=='\\') mm[q]='d';
							if (mm[q]=='*') mm[q]='x';
                        }
		
						ofstream of1("./tmp1");
						of1.close();
						if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str()); //A.ion[s].molesmi
						ifstream inf("./tmp1");
						inf >> ws;

						if (inf.eof()) {
							count++;
							if (1) {
								if (s==0) cout << "EXHAUSTIVE COMBINATION: " << (smi1_cat) << " { "<<ACi<<" } + " << (smi2_cat) << " { "<<BCj<<" } , bnd " << b << " -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
								if (s==1) cout << "EXHAUSTIVE COMBINATION: " << (smi1_an) << " { "<<ACi<<" } + " << (smi2_an) << " { "<<BCj<<" } , bnd " << b << " -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							}
						
						
							out << setw(7) << left << count << " "
								<< setw(65) << left << A.ion[s].molesmi << " "
								<< setw(17) << left << ACi << " "
				        		<< setw(17) << left << BCj << " "
								<< setw(17) << left << b << endl;					    	
							
							if (para.ifwritemds) {
								ofstream outs((para.mdsdir+mm+"_comb_IL.enc").c_str());
								A.ion[s].printmds(outs);
								outs.close();
							}

							ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
							outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_comb_IL.enc") << endl;
							outs.close();

							if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_comb_IL.enc") << endl; //A.ion[s].molesmi
							if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_comb_IL.enc") << endl; //A.ion[s].molesmi

							if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
						}
						inf.close();
	        
						if (0) {
							if (s==0) A.ion[s].smiles=smi1_cat;
							if (s==1) A.ion[s].smiles=smi1_an;
							A.ion[s].input();
						}
						if (1) {
							A.ion[s].replace(buf[0].ion[s]);
							//A.ion[s].reset();
						}
	
						if (0) {
							if (s==0) B.ion[s].smiles=smi2_cat;
							if (s==1) B.ion[s].smiles=smi2_an;
							B.ion[s].input();
						}
						if (1) {
							B.ion[s].replace(buf[1].ion[s]);
							//B.ion[s].reset();
						}
	        		}    
				}

			}
		}

		if (s==0) para.stat << "CAT | COMBINATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
		if (s==1) para.stat << "AN | COMBINATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
	}

	
	out.close();

	delete [] buf;
	buf=NULL;

	return 1;
}

unsigned int exhaustive_crossover1(IL &A, IL &B) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;
	string smi2_cat=B.ion[0].molesmi,smi2_an=B.ion[1].molesmi;
	
    IL *buf=new IL [2];
    buf[0].replace(A);
    buf[1].replace(B);

	ofstream out((para.logdir+"crossover_IL.txt").c_str(),ios::app);
	//mark(out);

    out << "IL1: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;
    out << "IL2: " << (smi2_cat+"."+smi2_an)  << endl;
    B.printmds(out);
    out << endl;
    
    for (unsigned int s=0;s<2;s++) {
		if (s==0) {
			count=0;
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL1_cat_cros" << " "
		        << setw(65) << left << "IL2_cat_cros" << " "
		        << setw(17) << left << "IL1_cat_cros_pt" << " "
		        << setw(17) << left << "IL2_cat_cros_pt" << endl;			
		}
		if (s==1) {
			count=0;
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL1_an_cros" << " "
		        << setw(65) << left << "IL2_an_cros" << " "
		        << setw(17) << left << "IL1_an_cros_pt" << " "
		        << setw(17) << left << "IL2_an_cros_pt" << endl;			
		}
	
		for (i=1;i<A.ion[s].Cindex.size();i++) {
			for (j=1;j<B.ion[s].Cindex.size();j++) {

                unsigned int ACi=A.ion[s].Cindex.at(i);
                unsigned int BCj=B.ion[s].Cindex.at(j);


				if (A.ion[s].crossover(B.ion[s],i,j)) {

					//A.ion[s].mds2smi();
					if (0) {
	            		A.ion[s].smiles=A.ion[s].molesmi;
						A.ion[s].input();
					}
					if (1) A.ion[s].canonicalize_SMILES();

					//B.ion[s].mds2smi();
					if (0) {
	            		B.ion[s].smiles=B.ion[s].molesmi;
						B.ion[s].input();
					}
					if (1) B.ion[s].canonicalize_SMILES();
					
					bool go=0;

					A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
					//out<<A.molesmi<<endl;

                    string mm=A.ion[s].molesmi;
                    for (unsigned int q=0;q<mm.length();q++) {
                        if (mm[q]=='/') mm[q]='u';
                        if (mm[q]=='\\') mm[q]='d';
						if (mm[q]=='*') mm[q]='x';
                    }

					ofstream of1("./tmp1");
					of1.close();
                    if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    ifstream inf("./tmp1");
                    inf >> ws;

					if (inf.eof()) {
						if (para.ifwritemds) {
							ofstream outs((para.mdsdir+mm+"_cros_IL.enc").c_str());
							A.ion[s].printmds(outs);
							outs.close();
						}
						go=1;

                        ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
                        outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_cros_IL.enc") << endl;
                        outs.close();

                        if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_cros_IL.enc") << endl; //A.ion[s].molesmi
                        if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_cros_IL.enc") << endl; //A.ion[s].molesmi

						if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
					}
					inf.close();
					
					B.molesmi=B.ion[0].molesmi+"."+B.ion[1].molesmi;
					//out<<B.molesmi<<endl;

                    mm=B.ion[s].molesmi;
                    for (unsigned int q=0;q<mm.length();q++) {
                        if (mm[q]=='/') mm[q]='u';
                        if (mm[q]=='\\') mm[q]='d';
						if (mm[q]=='*') mm[q]='x';
                    }

					of1.open("./tmp1");
					of1.close();
                    if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    inf.open("./tmp1");
                    inf >> ws;

					if (inf.eof()) {
						if (para.ifwritemds) {
							ofstream outs((para.mdsdir+mm+"_cros_IL.enc").c_str());
							B.ion[s].printmds(outs);
							outs.close();
						}
						go=1;

						ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
						outs << setw(50) << left << B.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_cros_IL.enc") << endl;
                        outs.close();

                        if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_cros_IL.enc") << endl; //A.ion[s].molesmi
                        if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_cros_IL.enc") << endl; //A.ion[s].molesmi

						if (para.enumeration) SMI_Enumerator(B.ion[s].molesmi);
					}
					inf.close();
	
					if (go) {
						count++;
						if (1) {
							if (s==0) cout << "EXHAUSTIVE CROSSOVER: " << (smi1_cat) << " { "<<ACi<<" } + " << (smi2_cat) << " { "<<BCj<<" } -> " << A.ion[s].molesmi << " + " << B.ion[s].molesmi << " {No. " << count << " }" << endl;
							if (s==1) cout << "EXHAUSTIVE CROSSOVER: " << (smi1_an) << " { "<<ACi<<" } + " << (smi2_an) << " { "<<BCj<<" } -> " << A.ion[s].molesmi << " + " << B.ion[s].molesmi << " {No. " << count << " }" << endl;
						}
					    out << setw(7) << left << count << " "
					        << setw(65) << left << A.ion[s].molesmi << " "
					        << setw(65) << left << B.ion[s].molesmi << " "
					        << setw(17) << left << ACi << " "
					        << setw(17) << left << BCj << endl;	
					}
	
					if (0) {
						if (s==0) A.ion[s].smiles=smi1_cat;
						if (s==1) A.ion[s].smiles=smi1_an;
						A.ion[s].input();	
					}
					if (1) {
						A.ion[s].replace(buf[0].ion[s]);
						//A.ion[s].reset();
					}
	
					if (0) {
						if (s==0) B.ion[s].smiles=smi2_cat;
						if (s==1) B.ion[s].smiles=smi2_an;
						B.ion[s].input();
					}
					if (1) {
						B.ion[s].replace(buf[1].ion[s]);
						//B.ion[s].reset();
					}
					
				}				
			}
		}


        if (s==0) para.stat << "CAT | CROSSOVER | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | CROSSOVER | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

	}
	
	out.close();

	delete [] buf;
	buf=NULL;

	return 1;
}


unsigned int exhaustive_subtraction1(IL &A) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);
	
	ofstream out((para.logdir+"subtraction_IL.txt").c_str(),ios::app);
	//mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;

	for (unsigned int s=0;s<2;s++) {
		if (s==0) {
			count=0;
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_cat_subtr" << " "
		        << setw(17) << left << "cat_subtr_pt" << " "
				<< setw(17) << left << "cat_subtr_bndfrm" << endl;			
		}
		if (s==1) {
			count=0;
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_an_subtr" << " "
		        << setw(17) << left << "an_subtr_pt" << " " 
				<< setw(17) << left << "an_subtr_bndfrm" << endl;			
		}
	
		for (i=0;i<A.ion[s].Cindex.size();i++) {
			for (unsigned int bd=0;bd<4;bd++){
				unsigned int ACi=A.ion[s].Cindex.at(i);

				if (A.ion[s].subtraction(i,bd,0)) {
					//A.ion[s].mds2smi();
					if (0) {
        	    		A.ion[s].smiles=A.ion[s].molesmi;
						A.ion[s].input();
					}
					if (1) A.ion[s].canonicalize_SMILES();
				
					A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
					//out<<A.molesmi<<endl;

                	string mm=A.ion[s].molesmi;
                	for (unsigned int q=0;q<mm.length();q++) {
                    	if (mm[q]=='/') mm[q]='u';
                    	if (mm[q]=='\\') mm[q]='d';
						if (mm[q]=='*') mm[q]='x';
                	}

					ofstream of1("./tmp1");
                    of1.close();
                    if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                	ifstream inf("./tmp1");
                	inf >> ws;

					if (inf.eof()) {
						count++;
						if (1) {
							if (s==0) cout << "EXHAUSTIVE SUBTRACTION: " << (smi1_cat) << " { "<<ACi<<" } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							if (s==1) cout << "EXHAUSTIVE SUBTRACTION: " << (smi1_an) << " { "<<ACi<<" } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
						}
 
            			out << setw(7) << left << count << " "
                    		<< setw(65) << left << A.ion[s].molesmi << " "
                    		<< setw(17) << left << ACi << " "
							<< setw(17) << left << bd << endl;

						if (para.ifwritemds) {
							ofstream outs((para.mdsdir+mm+"_subtr_IL.enc").c_str());
							A.ion[s].printmds(outs);
							outs.close();
						}

						ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
						outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_subtr_IL.enc") << endl;
                    	outs.close();

                    	if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_subtr_IL.enc") << endl;   //A.ion[s].molesmi
                    	if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_subtr_IL.enc") << endl; //A.ion[s].molesmi

						if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
					}
					inf.close();
	
					if (0) {			
						if (s==0) A.ion[s].smiles=smi1_cat;
						if (s==1) A.ion[s].smiles=smi1_an;
						A.ion[s].input();
					}
					if (1) {
						A.ion[s].replace(buf[0].ion[s]);
						//A.ion[s].reset();
					}
				}
			}
		}	

        if (s==0) para.stat << "CAT | SUBTRACTION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | SUBTRACTION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

	}

	
	out.close();

	delete [] buf;
	buf=NULL;

	return 1;
}


unsigned int exhaustive_addition1(IL &A) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

	ofstream out((para.logdir+"addition_IL.txt").c_str(),ios::app);
	//mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;
    
    for (unsigned int s=0;s<2;s++) {
    	if (s==0) {
			count=0;
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_cat_add" << " "
		        << setw(17) << left << "cat_add_pt" << " "
				<< setw(17) << left << "cat_add_id" << " "
				<< setw(7) << left << "cat_add_bnd" << " "
				<< setw(30) << left << "cat_add_name" << endl;    		
		}
    	if (s==1) {
			count=0;
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_an_add" << " "
		        << setw(17) << left << "an_add_pt" << " "
				<< setw(17) << left << "an_add_id" << " "
				<< setw(7) << left << "an_add_bnd" << " "
				<< setw(30) << left << "an_add_name" << endl;    		
		}
	
		for (i=0;i<A.ion[s].Cindex.size();i++) {
			for (j=1;j<A.ion[s].data->num;j++) {
				if (A.ion[s].data->a[j].probability>0){
					unsigned int ACi=A.ion[s].Cindex.at(i);
					for (unsigned int b=1;b<4;b++) {
						if (A.ion[s].addition(i,j,b)) {
							//A.ion[s].mds2smi();
							if (0) {
								A.ion[s].smiles=A.ion[s].molesmi;
								A.ion[s].input();
							}
							if (1) A.ion[s].canonicalize_SMILES();
						
							A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
							//out<<A.molesmi<<endl;

							string mm=A.ion[s].molesmi;
							for (unsigned int q=0;q<mm.length();q++) {
								if (mm[q]=='/') mm[q]='u';
								if (mm[q]=='\\') mm[q]='d';
								if (mm[q]=='*') mm[q]='x';
							}

							ofstream of1("./tmp1");
                    		of1.close();
                    		if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
							ifstream inf("./tmp1");
							inf >> ws;

							if (inf.eof()) {
								count++;
								if (1) {
									if (s==0) cout << "EXHAUSTIVE ADDITION: " << (smi1_cat) << " { "<< A.ion[s].data->a[j].id << " " << A.ion[s].data->a[j].name << " on " <<ACi<<" , bnd " << b << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
									if (s==1) cout << "EXHAUSTIVE ADDITION: " << (smi1_an) << " { "<< A.ion[s].data->a[j].id << " " << A.ion[s].data->a[j].name << " on " <<ACi<<" , bnd " << b << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
								}
								out << setw(7) << left << count << " "
									<< setw(65) << left << A.ion[s].molesmi << " "
									<< setw(17) << left << ACi << " "
									<< setw(17) << left << j << " "
									<< setw(7) << left << b << " "
									<< setw(30) << left << A.ion[s].data->a[j].name << endl;    

								if (para.ifwritemds) {
									ofstream outs((para.mdsdir+mm+"_add_IL.enc").c_str());
									A.ion[s].printmds(outs);
									outs.close();
								}

								ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
								outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_add_IL.enc") << endl;
								outs.close();

								if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_add_IL.enc") << endl; //A.ion[s].molesmi
								if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_add_IL.enc") << endl; //A.ion[s].molesmi

								if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
							}
							inf.close();
						
							if (0) {
								if (s==0) A.ion[s].smiles=smi1_cat;
								if (s==1) A.ion[s].smiles=smi1_an;
								A.ion[s].input();	
							}
							if (1) {
								A.ion[s].replace(buf[0].ion[s]);
								//A.ion[s].reset();
							}
						}
					}

				}
			}
		}
        if (s==0) para.stat << "CAT | ADDITION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | ADDITION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

	}

	
	out.close();

	delete [] buf;
	buf=NULL;

	return 1;
}

unsigned int exhaustive_insertion1(IL &A) {
    unsigned int count=0,i=0,j=0,k=0,m=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.logdir+"insertion_IL.txt").c_str(),ios::app);
    //mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;


    for (unsigned int s=0;s<2;s++) {
        if (s==0) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_ins" << " "
				<< setw(17) << left << "cat_ins_oC" << " "
                << setw(17) << left << "cat_ins_oP" << " "
                << setw(17) << left << "cat_ins_id" << " "
                << setw(30) << left << "cat_ins_name" << " "
                << setw(17) << left << "cat_bnd2par" << " "
                << setw(17) << left << "cat_bnd2des" << endl;
        }
        if (s==1) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_ins" << " "
				<< setw(17) << left << "an_ins_oC" << " "
                << setw(17) << left << "an_ins_oP" << " "
                << setw(17) << left << "an_ins_id" << " "
                << setw(30) << left << "an_ins_name" << " "
                << setw(17) << left << "an_bnd2par" << " "
                << setw(17) << left << "an_bnd2des" << endl;
        }
        for (i=1;i<A.ion[s].Cindex.size();i++) {
            for (j=1;j<A.ion[s].data->num;j++) {
                if (A.ion[s].data->a[j].probability>0) {
                    for (k=1;k<4;k++) {
                        for (m=1;m<4;m++) {
                            unsigned int ACi=A.ion[s].Cindex.at(i);
                            unsigned int APi=A.ion[s].Pindex.at(i);

                            if (A.ion[s].insertion(i,j,k,m)) {
                                //A.mds2smi();
                                if (0) {
                                    A.ion[s].smiles=A.ion[s].molesmi;
                                    A.ion[s].input();
                                }
                                if (1) A.ion[s].canonicalize_SMILES();

                                A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;

                                string mm=A.ion[s].molesmi;
                                for (unsigned int q=0;q<mm.length();q++) {
                                    if (mm[q]=='/') mm[q]='u';
                                    if (mm[q]=='\\') mm[q]='d';
									if (mm[q]=='*') mm[q]='x';
                                }

								ofstream of1("./tmp1");
                    			of1.close();
                    			if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                                ifstream inf("./tmp1");
                                inf >> ws;

                                if (inf.eof()) {
                                    count++;
									if (1) {
                                    	if (s==0) {
                                        	cout << "EXHAUSTIVE INSERTION: " << (smi1_cat) << " {CAT: " << A.ion[s].data->a[j].id << " " << A.ion[s].data->a[j].name << " on " << ACi << " , bnd2par " << k << " , bnd2des " << m << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
                                    	}
                                    	if (s==1) {
                                        	cout << "EXHAUSTIVE INSERTION: " << (smi1_an) << " {AN: " << A.ion[s].data->a[j].id << " " << A.ion[s].data->a[j].name << " on " << ACi << ", bnd2par " << k << " , bnd2des " << m << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
                                    	}
									}

                                    out << setw(7) << left << count << " "
                                        << setw(65) << left << A.ion[s].molesmi << " "
										<< setw(17) << left << APi << " "
                                        << setw(17) << left << ACi << " "
                                        << setw(17) << left << A.ion[s].data->a[j].id << " "
                                        << setw(30) << left << A.ion[s].data->a[j].name << " "
                                        << setw(17) << left << k << " "
                                        << setw(17) << left << m << endl;

									if (para.ifwritemds) {
										ofstream outs((para.mdsdir+mm+"_ins_IL.enc").c_str());
                                    	A.ion[s].printmds(outs);
                                    	outs.close();
									}

									ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
									outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_ins_IL.enc") << endl;
                                    outs.close();

                            		if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_ins_IL.enc") << endl; //A.ion[s].molesmi
                            		if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_ins_IL.enc") << endl; //A.ion[s].molesmi

									if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
                                }
                                inf.close();

                                if (0) {
                                    if (s==0) A.ion[s].smiles=smi1_cat;
                                    if (s==1) A.ion[s].smiles=smi1_an;
                                    A.ion[s].input();
                                }
                                if (1) {
                                    A.ion[s].replace(buf[0].ion[s]);
                                    //A.ion[s].reset();
                                }


                            }

                        }
                    }
                }
            }
        }
        if (s==0) para.stat << "CAT | INSERTION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | INSERTION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

    }



    out.close();

    delete [] buf;
    buf=NULL;

	return 1;
}


unsigned int exhaustive_change_bnd1(IL &A) {
    unsigned int count=0,i=0,j=0,k=0,m=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.logdir+"change_bnd_IL.txt").c_str(),ios::app);
    //mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;


    for (unsigned int s=0;s<2;s++) {
        if (s==0) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_chbnd" << " "
                << setw(17) << left << "cat_chbnd_oP" << " "
                << setw(17) << left << "cat_chbnd_oC" << " "
                << setw(17) << left << "cat_chbnd_obnd" << " "
				<< setw(17) << left << "cat_chbnd_nid_P" << " "
				<< setw(30) << left << "cat_chbnd_nPname" << " "
				<< setw(17) << left << "cat_chbnd_nid_C" << " "
				<< setw(30) << left << "cat_chbnd_nCname" << " "
				<< setw(17) << left << "cat_chbnd_nbnd" << endl;
        }
        if (s==1) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_chbnd" << " "
                << setw(17) << left << "an_chbnd_oP" << " "
                << setw(17) << left << "an_chbnd_oC" << " "
                << setw(17) << left << "an_chbnd_obnd" << " "
                << setw(17) << left << "an_chbnd_nid_P" << " "
				<< setw(30) << left << "an_chbnd_nPname" << " "
                << setw(17) << left << "an_chbnd_nid_C" << " "
				<< setw(30) << left << "an_chbnd_nCname" << " "
				<< setw(17) << left << "an_chbnd_nbnd" << endl;
        }

        for (i=1;i<A.ion[s].Cindex.size();i++) {
            for (j=1;j<A.ion[s].data->num;j++) {
                if (A.ion[s].data->a[j].probability>0) {
                    for (k=1;k<A.ion[s].data->num;k++) {
						if (A.ion[s].data->a[k].probability>0) {
                        	for (m=1;m<4;m++) {
                            	unsigned int ACi=A.ion[s].Cindex.at(i);
                            	unsigned int APi=A.ion[s].Pindex.at(i);
								unsigned int obnd=A.ion[s].Rindex.at(i);

                            	if (A.ion[s].change_bnd(i,j,k,m)) {
                                	//A.mds2smi();
                                	if (0) {
                                    	A.ion[s].smiles=A.ion[s].molesmi;
                                    	A.ion[s].input();
                                	}
                                	if (1) A.ion[s].canonicalize_SMILES();

                                	A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;

                                	string mm=A.ion[s].molesmi;
                                	for (unsigned int q=0;q<mm.length();q++) {
                                    	if (mm[q]=='/') mm[q]='u';
                                    	if (mm[q]=='\\') mm[q]='d';
										if (mm[q]=='*') mm[q]='x';
                                	}

									ofstream of1("./tmp1");
                    				of1.close();
                    				if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                                	ifstream inf("./tmp1");
                                	inf >> ws;

                                	if (inf.eof()) {
                                    	count++;
										if (1) {
                                    		if (s==0) {
                                        		cout << "EXHAUSTIVE CHANGE_BND: " << (smi1_cat) << " {CAT: " << A.ion[s].data->a[j].id << " " << A.ion[s].data->a[j].name << " on " << APi << " | " << A.ion[s].data->a[k].id << " " << A.ion[s].data->a[k].name << " on " << ACi << " | bnd " << obnd << " to " << m << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
                                    		}
                                    		if (s==1) {
                                        		cout << "EXHAUSTIVE CHANGE_BND: " << (smi1_an) << " {AN: " << A.ion[s].data->a[j].id << " " << A.ion[s].data->a[j].name << " on " << APi << " | " << A.ion[s].data->a[k].id << " " << A.ion[s].data->a[k].name << " on " << ACi << " | bnd " << obnd << " to " << m << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
                                    		}
										}

                                    	out << setw(7) << left << count << " "
                                        	<< setw(65) << left << A.ion[s].molesmi << " "
                                        	<< setw(17) << left << APi << " "
											<< setw(17) << left << ACi << " "
											<< setw(17) << left << obnd << " "
                                        	<< setw(17) << left << A.ion[s].data->a[j].id << " "
                                        	<< setw(30) << left << A.ion[s].data->a[j].name << " "
											<< setw(17) << left << A.ion[s].data->a[k].id << " "
                                            << setw(30) << left << A.ion[s].data->a[k].name << " "
                                        	<< setw(17) << left << m << endl;

										if (para.ifwritemds) {
											ofstream outs((para.mdsdir+mm+"_chbnd_IL.enc").c_str());
                                    		A.ion[s].printmds(outs);
                                    		outs.close();
										}

										ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
										outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_chbnd_IL.enc") << endl;
                                    	outs.close();

			                            if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_chbnd_IL.enc") << endl; //A.ion[s].molesmi
            			                if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_chbnd_IL.enc") << endl; //A.ion[s].molesmi

                                    	if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
                                	}
                                	inf.close();

                               		if (0) {
                                    	if (s==0) A.ion[s].smiles=smi1_cat;
                                    	if (s==1) A.ion[s].smiles=smi1_an;
                                    	A.ion[s].input();
                                	}
                                	if (1) {
                                    	A.ion[s].replace(buf[0].ion[s]);
                                    	//A.ion[s].reset();
                                	}
                            	}
                        	}
                    	}
					}
                }
            }
        }
        if (s==0) para.stat << "CAT | CHANGE_BND | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | CHANGE_BND | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

    }



    out.close();

    delete [] buf;
    buf=NULL;

	return 1;
}



unsigned int exhaustive_change_ele1(IL &A) {
    unsigned int count=0,i=0,j=0,k=0,m=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.logdir+"change_ele_IL.txt").c_str(),ios::app);
    //mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;


	for (unsigned int s=0;s<2;s++) {
		if (s==0) {
			count=0;
    		out << setw(7) << left << "No." << " "
        		<< setw(65) << left << "IL_chele" << " "
        		<< setw(17) << left << "cat_chele_pt" << " "
        		<< setw(17) << left << "cat_chele_id" << " "
        		<< setw(30) << left << "cat_chele_name" << " "
        		<< setw(17) << left << "cat_bnd2par" << " "
        		<< setw(17) << left << "cat_bnd2des" << endl;
		}
        if (s==1) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_chele" << " "
                << setw(17) << left << "an_chele_pt" << " "
                << setw(17) << left << "an_chele_id" << " "
                << setw(30) << left << "an_chele_name" << " "
                << setw(17) << left << "an_bnd2par" << " "
                << setw(17) << left << "an_bnd2des" << endl;
        }

    	for (i=0;i<A.ion[s].Cindex.size();i++) {
        	for (j=1;j<A.ion[s].data->num;j++) {
            	if (A.ion[s].data->a[j].probability>0) {
                	for (k=0;k<4;k++) {
						for (m=0;m<4;m++) {
                        	unsigned int ACi=A.ion[s].Cindex.at(i);

                        	if (A.ion[s].change_ele(i,j,k,m)) {
                            	//A.mds2smi();
                            	if (0) {
                                	A.ion[s].smiles=A.ion[s].molesmi;
                                	A.ion[s].input();
                            	}
                            	if (1) A.ion[s].canonicalize_SMILES();

                            	A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;

                                string mm=A.ion[s].molesmi;
                                for (unsigned int q=0;q<mm.length();q++) {
                                    if (mm[q]=='/') mm[q]='u';
                                    if (mm[q]=='\\') mm[q]='d';
									if (mm[q]=='*') mm[q]='x';
                                }

								ofstream of1("./tmp1");
                                of1.close();
                                if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                            	ifstream inf("./tmp1");
                            	inf >> ws;

                            	if (inf.eof()) {
                            		count++;
									if (1) {
										if (s==0) {
                                			cout << "EXHAUSTIVE CHANGE_ELE: " << (smi1_cat) << " {CAT: " << A.ion[s].data->a[j].id << " " << A.ion[s].data->a[j].name << " on " << ACi << " | bnd2par " << k << " | bnd2des " << m << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
										}
										if (s==1) {
                                        	cout << "EXHAUSTIVE CHANGE_ELE: " << (smi1_an) << " {AN: " << A.ion[s].data->a[j].id << " " << A.ion[s].data->a[j].name << " on " << ACi << " | bnd2par " << k << " | bnd2des " << m << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
										}
									}

                                	out << setw(7) << left << count << " "
                                    	<< setw(65) << left << A.ion[s].molesmi << " "
                                    	<< setw(17) << left << ACi << " "
                                    	<< setw(17) << left << A.ion[s].data->a[j].id << " "
                                    	<< setw(30) << left << A.ion[s].data->a[j].name << " "
                                    	<< setw(17) << left << k << " "
                                    	<< setw(17) << left << m << endl;

									if (para.ifwritemds) {
										ofstream outs((para.mdsdir+mm+"_chele_IL.enc").c_str());
                                		A.ion[s].printmds(outs);
                                		outs.close();
									}

									ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
									outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_chele_IL.enc") << endl;
                                	outs.close();

		                            if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_chele_IL.enc") << endl; //A.ion[s].molesmi
        		                    if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_chele_IL.enc") << endl; ////A.ion[s].molesmi

									if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
                            	}
                            	inf.close();

                            	if (0) {
                                	if (s==0) A.ion[s].smiles=smi1_cat;
                                	if (s==1) A.ion[s].smiles=smi1_an;
                                	A.ion[s].input();
                            	}
                            	if (1) {
                                	A.ion[s].replace(buf[0].ion[s]);
                                	//A.ion[s].reset();
                            	}

						
							}

						}
					}
				}
			}
		}
        if (s==0) para.stat << "CAT | CHANGE_ELE | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | CHANGE_ELE | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

	}



    out.close();

    delete [] buf;
    buf=NULL;

	return 1;
}


unsigned int exhaustive_cyclization1(IL &A) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);
	
	ofstream out((para.logdir+"cyclization_IL.txt").c_str(),ios::app);
	//mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;

	for (unsigned int s=0;s<2;s++) {
		if (s==0) {
			count=0;
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_cat_cyc" << " "
		        << setw(17) << left << "cat_cyc_pt1" << " "
		        << setw(17) << left << "cat_cyc_pt2" << " "
				<< setw(7) << left << "bnd" << endl;		
		}
		if (s==1) {
			count=0;
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_an_cyc" << " "
		        << setw(17) << left << "an_cyc_pt1" << " "
		        << setw(17) << left << "an_cyc_pt2" << " "
				<< setw(7) << left << "bnd" << endl;		
		}
	
	
		for (i=0;i<A.ion[s].Cindex.size();i++) {
			for (j=i+1;j<A.ion[s].Cindex.size();j++) {
				for (unsigned int k=1;k<=3;k++) {
				unsigned int ACi=A.ion[s].Cindex.at(i);
				unsigned int ACj=A.ion[s].Cindex.at(j);

				if (A.ion[s].cyclization(i,j,k)) {
					//A.ion[s].mds2smi();
					if (0) {
        				A.ion[s].smiles=A.ion[s].molesmi;
						A.ion[s].input();			
					}
					if (1) A.ion[s].canonicalize_SMILES();
					
					A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
					//out<<A.molesmi<<endl;

                    string mm=A.ion[s].molesmi;
                    for (unsigned int q=0;q<mm.length();q++) {
                        if (mm[q]=='/') mm[q]='u';
                        if (mm[q]=='\\') mm[q]='d';
						if (mm[q]=='*') mm[q]='x';
                    }

					ofstream of1("./tmp1");
                    of1.close();
                    if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    ifstream inf("./tmp1");
                    inf >> ws;

					if (inf.eof()) {
						count++;
						if (1) {
							if (s==0) cout << "EXHAUSTIVE CYCLIZATION: " << (smi1_cat) << " {CAT: "<< ACi << " , " << ACj << " | bnd " << k << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							if (s==1) cout << "EXHAUSTIVE CYCLIZATION: " << (smi1_an) << " {AN: "<< ACi << " , " << ACj << " | bnd " << k << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
						}
						
					    out << setw(7) << left << count << " "
					        << setw(65) << left << A.ion[s].molesmi << " "
					        << setw(17) << left << ACi << " "
					        << setw(17) << left << ACj << " "
							<< setw(7) << left << k << endl;

						if (para.ifwritemds) {
							ofstream outs((para.mdsdir+mm+"_cyc_IL.enc").c_str());
							A.ion[s].printmds(outs);
							outs.close();
						}

						ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
						outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_cyc_IL.enc") << endl;
                        outs.close();

                        if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_cyc_IL.enc") << endl; //A.ion[s].molesmi
                        if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_cyc_IL.enc") << endl; //A.ion[s].molesmi

						if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
					}
					inf.close();
					
					if (0) {
						if (s==0) A.ion[s].smiles=smi1_cat;
						if (s==1) A.ion[s].smiles=smi1_an;
						A.ion[s].input();	
					}	
					if (1) {
						A.ion[s].replace(buf[0].ion[s]);			
						//A.ion[s].reset();
					}
				}
				}
			}
		}
        if (s==0) para.stat << "CAT | CYCLIZATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | CYCLIZATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

	}

	
	out.close();

	delete [] buf;
	buf=NULL;

	return 1;
}

unsigned int exhaustive_change_cistrans1(IL &A) {
    unsigned int count=0,i=0,j=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.logdir+"change_cistrans_IL.txt").c_str(),ios::app);
    //mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;
    for (unsigned int s=0;s<2;s++) {	
        if (s==0) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_cat_ct" << " "
                << setw(17) << left << "cat_ct_pt" << " "
                << setw(17) << left << "cat_ct_h_ori" << " " 
				<< setw(17) << left << "cat_ct_h_after" << " " 
				<< setw(17) << left << "cat_ct_e_ori" << " "
				<< setw(17) << left << "cat_ct_e_after" << endl;
        }
        if (s==1) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_an_ct" << " "
                << setw(17) << left << "an_ct_pt" << " "
				<< setw(17) << left << "an_ct_h_ori" << " "
                << setw(17) << left << "an_ct_h_after" << " "
                << setw(17) << left << "an_ct_e_ori" << " "
                << setw(17) << left << "an_ct_e_after" << endl;
        }

		for (unsigned int w=0;w<2;w++) {
        	for (i=0;i<A.ion[s].ctsisomer.at(w).size();i++) {
                unsigned int ACi=A.ion[s].Cindex.at(i);
                string ACTh=A.ion[s].ctsisomer.at(0).at(i);
				string ACTe=A.ion[s].ctsisomer.at(1).at(i);

                if (A.ion[s].change_cistrans(i,w)) {
                    //A.ion[s].mds2smi();
                    if (0) {
                        A.ion[s].smiles=A.ion[s].molesmi;
                        A.ion[s].input();
                    }
                    if (1) A.ion[s].canonicalize_SMILES();

                    A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
                    //out<<A.molesmi<<endl;

                    string mm=A.ion[s].molesmi;
                    for (unsigned int q=0;q<mm.length();q++) {
                        if (mm[q]=='/') mm[q]='u';
                        if (mm[q]=='\\') mm[q]='d';
						if (mm[q]=='*') mm[q]='x';
                    }

                    ofstream of1("./tmp1");
                    of1.close();
                    if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    ifstream inf("./tmp1");
                    inf >> ws;
                    if (inf.eof()) {
                        count++;
						if (1) {
                        	if (s==0) {
								if (w==0) cout << "EXHAUSTIVE CHANGE_CISTRANS: " << (smi1_cat) << " {CAT: "<< ACi << " | front | from " << ACTh << " to " << A.ion[s].ctsisomer.at(0).at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
								else if (w==1) cout << "EXHAUSTIVE CHANGE_CISTRANS: " << (smi1_cat) << " {CAT: "<< ACi << " | end | from " << ACTe << " to " << A.ion[s].ctsisomer.at(1).at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							}
                        	if (s==1) {
								if (w==0) cout << "EXHAUSTIVE CHANGE_CISTRANS: " << (smi1_an) << " {AN: "<< ACi << " | front | from " << ACTh << " to " << A.ion[s].ctsisomer.at(0).at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
								else if (w==1) cout << "EXHAUSTIVE CHANGE_CISTRANS: " << (smi1_an) << " {AN: "<< ACi << " | end | from " << ACTe << " to " << A.ion[s].ctsisomer.at(1).at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							}
						}

                        out << setw(7) << left << count << " "
                            << setw(65) << left << A.ion[s].molesmi << " "
                            << setw(17) << left << ACi << " "
                            << setw(17) << left << ACTh << " "
							<< setw(17) << left << A.ion[s].ctsisomer.at(0).at(i) << " " 
							<< setw(17) << left << ACTe << " " 
							<< setw(17) << left << A.ion[s].ctsisomer.at(1).at(i) << endl;

                        //}

						if (para.ifwritemds) {
							ofstream outs((para.mdsdir+mm+"_chct_IL.enc").c_str());
                        	A.ion[s].printmds(outs);
                        	outs.close();
						}

						ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
						outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_chct_IL.enc") << endl;
                        outs.close();

                        if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_chct_IL.enc") << endl; //A.ion[s].molesmi
                        if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_chct_IL.enc") << endl; //A.ion[s].molesmi

						if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
                    }
                    inf.close();

                    if (0) {
                        if (s==0) A.ion[s].smiles=smi1_cat;
                        if (s==1) A.ion[s].smiles=smi1_an;
                        A.ion[s].input();
                    }
                    if (1) {
                        A.ion[s].replace(buf[0].ion[s]);
                        //A.ion[s].reset();
                    }
                }

        	}
		}
        if (s==0) para.stat << "CAT | CHANGE_CISTRANS | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | CHANGE_CISTRANS | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

    }


    out.close();

    delete [] buf;
    buf=NULL;

	return 1;
}


unsigned int exhaustive_change_chirality1(IL &A) {
    unsigned int count=0,i=0,j=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

	vector<unsigned int> stereoC(2);
	stereoC.at(0)=1;
	stereoC.at(1)=2;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.logdir+"change_chirality_IL.txt").c_str(),ios::app);
    //mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;
    for (unsigned int s=0;s<2;s++) {
        if (s==0) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_cat_chir" << " "
                << setw(17) << left << "cat_chir_pt" << " "
                << setw(17) << left << "cat_chir_ori" << " "
                << setw(17) << left << "cat_chir_after" << endl;
        }
        if (s==1) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_an_chir" << " "
                << setw(17) << left << "an_chir_pt" << " "
                << setw(17) << left << "an_chir_ori" << " "
                << setw(17) << left << "an_chir_after" << endl;
        }

        for (i=0;i<A.ion[s].Cindex.size();i++) {
            unsigned int ACi=A.ion[s].Cindex.at(i);
			unsigned int Achi=A.ion[s].chi.at(i);

			for (j=0;j<stereoC.size();j++) {
            	if (A.ion[s].change_chirality(i,stereoC.at(j))) {
                	//A.ion[s].mds2smi();
                	if (0) {
                    	A.ion[s].smiles=A.ion[s].molesmi;
                    	A.ion[s].input();
                	}
                	if (1) A.ion[s].canonicalize_SMILES();

                	A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
                	//out<<A.molesmi<<endl;

                    string mm=A.ion[s].molesmi;
                    for (unsigned int q=0;q<mm.length();q++) {
                        if (mm[q]=='/') mm[q]='u';
                        if (mm[q]=='\\') mm[q]='d';
						if (mm[q]=='*') mm[q]='x';
                    }

                    ofstream of1("./tmp1");
                    of1.close();
                    if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                	ifstream inf("./tmp1");
                	inf >> ws;
                	if (inf.eof()) {
                    	count++;
						if (1) {
                    		if (s==0) {
                        		cout << "EXHAUSTIVE CHANGE_CHIRALITY: " << (smi1_cat) << " {CAT: "<< ACi << " | " << Achi << " to " << A.ion[s].chi.at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
                    		}
                    		if (s==1) {
                        		cout << "EXHAUSTIVE CHANGE_CHIRALITY: " << (smi1_an) << " {AN: "<< ACi << " | " << Achi << " to " << A.ion[s].chi.at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
                    		}
						}

                    	out << setw(7) << left << count << " "
                        	<< setw(65) << left << A.ion[s].molesmi << " "
                        	<< setw(17) << left << ACi << " "
                        	<< setw(17) << left << Achi << " "
                        	<< setw(17) << left << A.ion[s].chi.at(i) << endl;

						if (para.ifwritemds) {
							ofstream outs((para.mdsdir+mm+"_chchir_IL.enc").c_str());
                    		A.ion[s].printmds(outs);
                    		outs.close();
						}

						ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
						outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_chchir_IL.enc") << endl;
                    	outs.close();

                        if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_chchir_IL.enc") << endl; //A.ion[s].molesmi
                        if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_chchir_IL.enc") << endl; //A.ion[s].molesmi

						if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
                	}
                	inf.close();

                	if (0) {
                    	if (s==0) A.ion[s].smiles=smi1_cat;
                    	if (s==1) A.ion[s].smiles=smi1_an;
                    	A.ion[s].input();
                	}
                	if (1) {
                    	A.ion[s].replace(buf[0].ion[s]);
                    	//A.ion[s].reset();
                	}
            	}
        	}
		}
        if (s==0) para.stat << "CAT | CHANGE_CHIRALITY | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | CHANGE_CHIRALITY | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

    }


    out.close();

    delete [] buf;
    buf=NULL;

	return 1;
}
 

unsigned int exhaustive_decyclization1(IL &A) {
    unsigned int count=0,i=0,j=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.logdir+"decyclization_IL.txt").c_str(),ios::app);
    //mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an) << endl;
    A.printmds(out);
    out << endl;

    for (unsigned int s=0;s<2;s++) {
        if (s==0) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_cat_decyc" << " "
                << setw(17) << left << "cat_ring_no" << " "
                << setw(17) << left << "cat_cyc_bnd" << endl;
        }
        if (s==1) {
			count=0;
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_an_decyc" << " "
                << setw(17) << left << "an_ring_no" << " "
                << setw(17) << left << "an_cyc_bnd" << endl;
        }

    	for (i=1;i<=A.ion[s].if_circle;i++) {
        	unsigned int cycnum=i;
        	unsigned int cybnd=A.ion[s].Cybnd.at(cycnum-1);
        	if (A.ion[s].decyclization(i)) {
                //A.ion[s].mds2smi();
                if (0) {
                    A.ion[s].smiles=A.ion[s].molesmi;
                    A.ion[s].input();
                }
                if (1) A.ion[s].canonicalize_SMILES();

                A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
                //out<<A.molesmi<<endl;

                string mm=A.ion[s].molesmi;
                for (unsigned int q=0;q<mm.length();q++) {
                    if (mm[q]=='/') mm[q]='u';
                    if (mm[q]=='\\') mm[q]='d';
					if (mm[q]=='*') mm[q]='x';
                }

                ofstream of1("./tmp1");
                of1.close();
                if (para.redu_duplicates) system(("fgrep -F \""+para.mdsdir+mm+"_\" "+para.mdsdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                ifstream inf("./tmp1");
                inf >> ws;
                if (inf.eof()) {
                    count++;
					if (1) {
                    	if (s==0) {
                        	cout << "EXHAUSTIVE DECYCLIZATION: " << (smi1_cat) << " {CAT: ring no. " << cycnum << " | bnd " << cybnd << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
                    	}
                    	if (s==1) {
                        	cout << "EXHAUSTIVE DECYCLIZATION: " << (smi1_an) << " {AN: ring no. " << cycnum << " | bnd " << cybnd << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
                    	}
					}

                    out << setw(7) << left << count << " "
                        << setw(65) << left << A.ion[s].molesmi << " "
                        << setw(17) << left << cycnum << " "
                        << setw(17) << left << cybnd << endl;

					if (para.ifwritemds) {
						ofstream outs((para.mdsdir+mm+"_decyc_IL.enc").c_str());
                    	A.ion[s].printmds(outs);
                    	outs.close();
					}

					ofstream outs((para.mdsdir+"DATLIST.txt").c_str(),ios::app);
					outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.mdsdir+mm+"_decyc_IL.enc") << endl;
                    outs.close();

                    if (para.glbouf.is_open() && s==0) para.glbouf << (para.mdsdir+mm+"_decyc_IL.enc") << endl; //A.ion[s].molesmi
                    if (para.glbouf1.is_open() && s==1) para.glbouf1 << (para.mdsdir+mm+"_decyc_IL.enc") << endl;  //A.ion[s].molesmi

                    if (para.enumeration) SMI_Enumerator(A.ion[s].molesmi);
                }
                inf.close();

                if (0) {
                    if (s==0) A.ion[s].smiles=smi1_cat;
                    if (s==1) A.ion[s].smiles=smi1_an;
                    A.ion[s].input();
                }
                if (1) {
                    A.ion[s].replace(buf[0].ion[s]);
                    //A.ion[s].reset();
                }
			}
        }
        if (s==0) para.stat << "CAT | DECYCLIZATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
        if (s==1) para.stat << "AN | DECYCLIZATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;

    }

    out.close();

    delete [] buf;
    buf=NULL;

	return 1;
}

