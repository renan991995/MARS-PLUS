#include "CASES_IL_INDEPENDENT.h"
using namespace std;

extern PARAMETER para;

void exhaustive_combination1(IL &A, IL &B) {
	int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;
	string smi2_cat=B.ion[0].molesmi,smi2_an=B.ion[1].molesmi;
	
	IL *buf=new IL [2];
	buf[0].replace(A);
	buf[1].replace(B);
	
	ofstream out((para.smidir+"combination_IL.txt").c_str());
	mark(out);
	
	out << "IL1: " << (smi1_cat+"."+smi1_an)  << endl;
	A.print(out);
	out << endl;
	out << "IL2: " << (smi2_cat+"."+smi2_an)  << endl;
	B.print(out);
	out << endl;
		
	for (int s=0;s<2;s++) {
		if (s==0) {
			out << setw(7) << left << "No." << " "
				<< setw(65) << left << "IL_cat_comb" << " "
				<< setw(17) << left << "IL1_cat_comb_pt" << " "
	        	<< setw(17) << left << "IL2_cat_comb_pt" << " "
				<< setw(17) << left << "bnd" << endl;
	    }    	
	    if (s==1) {
			out << setw(7) << left << "No." << " "
				<< setw(65) << left << "IL_an_comb" << " "
				<< setw(17) << left << "IL1_an_comb_pt" << " "
	        	<< setw(17) << left << "IL2_an_comb_pt" << " "
				<< setw(17) << left << "bnd" << endl;
	    }    	

		for (int i=0;i<A.ion[s].Cindex.size();i++) {
			for (int j=0;j<B.ion[s].Cindex.size();j++) {
				int ACi=A.ion[s].Cindex.at(i);
				int BCj=B.ion[s].Cindex.at(j);

				for (int b=1;b<4;b++) {
	        		if (A.ion[s].combine(B.ion[s],i,j,b)) {
					
	        			A.ion[s].mds2smi();
						if (0) {
							A.ion[s].smiles=A.ion[s].molesmi;
							A.ion[s].input();
						}
						if (1) A.ion[s].canonicalize_SMILES();

						A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
						//out<<A.molesmi<<endl;

                            string mm=A.ion[s].molesmi;
                            for (int q=0;q<mm.length();q++) {
                                if (mm[q]=='/') mm[q]='u';
                                if (mm[q]=='\\') mm[q]='d';
                            }
		
						//system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
						system(("grep -Fw \""+mm+"\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str()); //A.ion[s].molesmi
						ifstream inf("./tmp1");
						inf >> ws;

						if (inf.eof()) {
							count++;
							if (s==0) cout << "EXHAUSTIVE COMBINATION: " << (smi1_cat) << " { "<<i<<" } + " << (smi2_cat) << " { "<<j<<" } , bnd " << b << " -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							if (s==1) cout << "EXHAUSTIVE COMBINATION: " << (smi1_an) << " { "<<i<<" } + " << (smi2_an) << " { "<<j<<" } , bnd " << b << " -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
						
						
							out << setw(7) << left << count << " "
								<< setw(65) << left << A.ion[s].molesmi << " "
								<< setw(17) << left << ACi << " "
				        		<< setw(17) << left << BCj << " "
								<< setw(17) << left << b << endl;					    	
							/*
							string mm=A.ion[s].molesmi;
							for (int q=0;q<mm.length();q++) {
								if (mm[q]=='/') mm[q]='u';
								if (mm[q]=='\\') mm[q]='d';
							}
							*/

							//ofstream outs((para.smidir+"mds/"+A.ion[s].molesmi+"_comb_IL.enc").c_str());						
							ofstream outs((para.smidir+"mds/"+mm+"_comb_IL.enc").c_str());
							A.ion[s].print(outs);
							outs.close();

							outs.open((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
							//outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+A.ion[s].molesmi+"_comb_IL.enc") << endl;
							outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+mm+"_comb_IL.enc") << endl;
							outs.close();
						}
						inf.close();
	        
						if (0) {
							if (s==0) A.ion[s].smiles=smi1_cat;
							if (s==1) A.ion[s].smiles=smi1_an;
							A.ion[s].input();
						}
						if (1) {
							A.ion[s].replace(buf[0].ion[s]);
							A.ion[s].reset();
						}
	
						if (0) {
							if (s==0) B.ion[s].smiles=smi2_cat;
							if (s==1) B.ion[s].smiles=smi2_an;
							B.ion[s].input();
						}
						if (1) {
							B.ion[s].replace(buf[1].ion[s]);
							B.ion[s].reset();
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

void exhaustive_crossover1(IL &A, IL &B) {
	int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;
	string smi2_cat=B.ion[0].molesmi,smi2_an=B.ion[1].molesmi;
	
    IL *buf=new IL [2];
    buf[0].replace(A);
    buf[1].replace(B);

	ofstream out((para.smidir+"crossover_IL.txt").c_str());
	mark(out);

    out << "IL1: " << (smi1_cat+"."+smi1_an)  << endl;
    A.print(out);
    out << endl;
    out << "IL2: " << (smi2_cat+"."+smi2_an)  << endl;
    B.print(out);
    out << endl;
    
    for (int s=0;s<2;s++) {
		if (s==0) {
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL1_cat_cros" << " "
		        << setw(65) << left << "IL2_cat_cros" << " "
		        << setw(17) << left << "IL1_cat_cros_pt" << " "
		        << setw(17) << left << "IL2_cat_cros_pt" << endl;			
		}
		if (s==1) {
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL1_an_cros" << " "
		        << setw(65) << left << "IL2_an_cros" << " "
		        << setw(17) << left << "IL1_an_cros_pt" << " "
		        << setw(17) << left << "IL2_an_cros_pt" << endl;			
		}
	
		for (i=1;i<A.ion[s].Cindex.size();i++) {
			for (j=1;j<B.ion[s].Cindex.size();j++) {

                int ACi=A.ion[s].Cindex.at(i);
                int BCj=B.ion[s].Cindex.at(j);


				if (A.ion[s].crossover(B.ion[s],i,j)) {

					A.ion[s].mds2smi();
					if (0) {
	            		A.ion[s].smiles=A.ion[s].molesmi;
						A.ion[s].input();
					}
					if (1) A.ion[s].canonicalize_SMILES();

					B.ion[s].mds2smi();
					if (0) {
	            		B.ion[s].smiles=B.ion[s].molesmi;
						B.ion[s].input();
					}
					if (1) B.ion[s].canonicalize_SMILES();
					
					bool go=0;

					A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
					//out<<A.molesmi<<endl;

                    //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                    system(("grep -Fw \""+A.ion[s].molesmi+"\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    ifstream inf("./tmp1");
                    inf >> ws;

					if (inf.eof()) {
                        string mm=A.ion[s].molesmi;
                        for (int q=0;q<mm.length();q++) {
                            if (mm[q]=='/') mm[q]='u';
                            if (mm[q]=='\\') mm[q]='d';
                        }

						ofstream outs((para.smidir+"mds/"+mm+"_cros_IL.enc").c_str());
						A.ion[s].print(outs);
						outs.close();
						go=1;

                        outs.open((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
                        outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+mm+"_cros_IL.enc") << endl;
                        outs.close();
					}
					inf.close();
					
					B.molesmi=B.ion[0].molesmi+"."+B.ion[1].molesmi;
					//out<<B.molesmi<<endl;

                    //system(("ls "+para.smidir+"mds/ | grep -Fw \""+B.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                    system(("grep -Fw \""+B.ion[s].molesmi+"\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    inf.open("./tmp1");
                    inf >> ws;

					if (inf.eof()) {
                        string mm=B.ion[s].molesmi;
                        for (int q=0;q<mm.length();q++) {
                            if (mm[q]=='/') mm[q]='u';
                            if (mm[q]=='\\') mm[q]='d';
                        }

						ofstream outs((para.smidir+"mds/"+mm+"_cros_IL.enc").c_str());
						B.ion[s].print(outs);
						outs.close();
						go=1;

                        outs.open((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
                        outs << setw(50) << left << B.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+mm+"_cros_IL.enc") << endl;
                        outs.close();
					}
					inf.close();
	
					if (go) {
						count++;
						if (s==0) cout << "EXHAUSTIVE CROSSOVER: " << (smi1_cat) << " { "<<i<<" } + " << (smi2_cat) << " { "<<j<<" } -> " << A.ion[s].molesmi << " + " << B.ion[s].molesmi << " {No. " << count << " }" << endl;
						if (s==1) cout << "EXHAUSTIVE CROSSOVER: " << (smi1_an) << " { "<<i<<" } + " << (smi2_an) << " { "<<j<<" } -> " << A.ion[s].molesmi << " + " << B.ion[s].molesmi << " {No. " << count << " }" << endl;
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
						A.ion[s].reset();
					}
	
					if (0) {
						if (s==0) B.ion[s].smiles=smi2_cat;
						if (s==1) B.ion[s].smiles=smi2_an;
						B.ion[s].input();
					}
					if (1) {
						B.ion[s].replace(buf[1].ion[s]);
						B.ion[s].reset();
					}
					
				}				
			}
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;

}


void exhaustive_subtraction1(IL &A) {
	int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);
	
	ofstream out((para.smidir+"subtraction_IL.txt").c_str());
	mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.print(out);
    out << endl;

	for (int s=0;s<2;s++) {
		if (s==0) {
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_cat_subtr" << " "
		        << setw(17) << left << "cat_subtr_pt" <<  endl;			
		}
		if (s==1) {
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_an_subtr" << " "
		        << setw(17) << left << "an_subtr_pt" <<  endl;			
		}
	
		for (i=1;i<A.ion[s].Cindex.size();i++) {
			int ACi=A.ion[s].Cindex.at(i);

			if (A.ion[s].subtract(i,1)) {
				A.ion[s].mds2smi();
				if (0) {
            		A.ion[s].smiles=A.ion[s].molesmi;
					A.ion[s].input();
				}
				if (1) A.ion[s].canonicalize_SMILES();
				
				A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
				//out<<A.molesmi<<endl;

                //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                system(("grep -Fw \""+A.ion[s].molesmi+"\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                ifstream inf("./tmp1");
                inf >> ws;

				if (inf.eof()) {
					count++;
					if (s==0) cout << "EXHAUSTIVE SUBTRACTION: " << (smi1_cat) << " { "<<i<<" } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
					if (s==1) cout << "EXHAUSTIVE SUBTRACTION: " << (smi1_an) << " { "<<i<<" } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
 
            		out << setw(7) << left << count << " "
                    	<< setw(65) << left << A.ion[s].molesmi << " "
                    	<< setw(17) << left << ACi << endl;

                    string mm=A.ion[s].molesmi;
                    for (int q=0;q<mm.length();q++) {
                        if (mm[q]=='/') mm[q]='u';
                        if (mm[q]=='\\') mm[q]='d';
                    }

		    		ofstream outs((para.smidir+"mds/"+mm+"_subtr_IL.enc").c_str());
					A.ion[s].print(outs);
					outs.close();

                    outs.open((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
                    outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+mm+"_subtr_IL.enc") << endl;
                    outs.close();
				}
				inf.close();
	
				if (0) {			
					if (s==0) A.ion[s].smiles=smi1_cat;
					if (s==1) A.ion[s].smiles=smi1_an;
					A.ion[s].input();
				}
				if (1) {
					A.ion[s].replace(buf[0].ion[s]);
					A.ion[s].reset();
				}
			}
		}	
	}

	
	out.close();

	delete [] buf;
	buf=NULL;


}


void exhaustive_addition1(IL &A) {
	int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

	ofstream out((para.smidir+"addition_IL.txt").c_str());
	mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.print(out);
    out << endl;
    
    for (int s=0;s<2;s++) {
    	if (s==0) {
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_cat_add" << " "
		        << setw(17) << left << "cat_add_pt" << " "
				<< setw(17) << left << "cat_add_id" << " "
				<< setw(30) << left << "cat_add_name" << endl;    		
		}
    	if (s==1) {
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_an_add" << " "
		        << setw(17) << left << "an_add_pt" << " "
				<< setw(17) << left << "an_add_id" << " "
				<< setw(30) << left << "an_add_name" << endl;    		
		}
	
		for (i=0;i<A.ion[s].Cindex.size();i++) {
			for (j=1;j<A.ion[s].data->num;j++) {
				if (A.ion[s].data->a[j].probability>0){

    				int ACi=A.ion[s].Cindex.at(i);

					for (int b=1;b<4;b++) {
					if (A.ion[s].add(i,j,b)) {
						A.ion[s].mds2smi();
						if (0) {
		            		A.ion[s].smiles=A.ion[s].molesmi;
							A.ion[s].input();
						}
						if (1) A.ion[s].canonicalize_SMILES();
						
						A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
						//out<<A.molesmi<<endl;

                		//system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                		system(("grep -Fw \""+A.ion[s].molesmi+"\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
        				ifstream inf("./tmp1");
        				inf >> ws;

						if (inf.eof()) {
							count++;
							if (s==0) cout << "EXHAUSTIVE ADDITION: " << (smi1_cat) << " { "<< j << " " << A.ion[s].data->a[j].name << " on " <<i<<" } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							if (s==1) cout << "EXHAUSTIVE ADDITION: " << (smi1_an) << " { "<< j << " " << A.ion[s].data->a[j].name << " on " <<i<<" } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;

						    out << setw(7) << left << count << " "
						        << setw(65) << left << A.ion[s].molesmi << " "
						        << setw(17) << left << ACi << " "
								<< setw(17) << left << j << " "
								<< setw(30) << left << A.ion[s].data->a[j].name << endl;    

                            string mm=A.ion[s].molesmi;
                            for (int q=0;q<mm.length();q++) {
                                if (mm[q]=='/') mm[q]='u';
                                if (mm[q]=='\\') mm[q]='d';
                            }

							ofstream outs((para.smidir+"mds/"+mm+"_add_IL.enc").c_str());
							A.ion[s].print(outs);
							outs.close();

                            outs.open((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
		                    outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+mm+"_add_IL.enc") << endl;
        		            outs.close();

						}
						inf.close();
						
						if (0) {
							if (s==0) A.ion[s].smiles=smi1_cat;
							if (s==1) A.ion[s].smiles=smi1_an;
							A.ion[s].input();	
						}
						if (1) {
							A.ion[s].replace(buf[0].ion[s]);
							A.ion[s].reset();
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

void exhaustive_exchange1(IL &A) {
    int count=0,i=0,j=0,k=0,m=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.smidir+"exchange_IL.txt").c_str());
    mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.print(out);
    out << endl;


	for (int s=0;s<2;s++) {
		if (s==0) {
    		out << setw(7) << left << "No." << " "
        		<< setw(65) << left << "IL_exc" << " "
        		<< setw(17) << left << "cat_exc_pt" << " "
        		<< setw(17) << left << "cat_exc_id" << " "
        		<< setw(30) << left << "cat_exc_name" << " "
        		<< setw(17) << left << "cat_bnd2par" << " "
        		<< setw(17) << left << "cat_bnd2des" << endl;
		}
        if (s==1) {
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_exc" << " "
                << setw(17) << left << "an_exc_pt" << " "
                << setw(17) << left << "an_exc_id" << " "
                << setw(30) << left << "an_exc_name" << " "
                << setw(17) << left << "an_bnd2par" << " "
                << setw(17) << left << "an_bnd2des" << endl;
        }

    	for (i=1;i<A.ion[s].Cindex.size();i++) {
        	for (j=1;j<A.ion[s].data->num;j++) {
            	if (A.ion[s].data->a[j].probability>0) {
                	for (k=1;k<4;k++) {
						for (m=1;m<4;m++) {
                        	int ACi=A.ion[s].Cindex.at(i);
                        	//int APi=A.ion[s].Pindex.at(i);

                        	if (A.ion[s].exchange(i,j,k,m)) {
                            	A.mds2smi();
                            	if (0) {
                                	A.ion[s].smiles=A.ion[s].molesmi;
                                	A.ion[s].input();
                            	}
                            	if (1) A.ion[s].canonicalize_SMILES();

                            	A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
								system(("grep -Fw \""+A.ion[s].molesmi+"\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                            	ifstream inf("./tmp1");
                            	inf >> ws;

                            	if (inf.eof()) {
                            		count++;
									if (s==0) {
                                		cout << "EXHAUSTIVE EXCHANGE: " << (smi1_cat) << " {CAT: " << j << " " << A.ion[s].data->a[j].name << " on " << ACi << " , bnd2par " << k << " bnd2des " << m << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
									}
									if (s==1) {
                                        cout << "EXHAUSTIVE EXCHANGE: " << (smi1_an) << " {AN: " << j << " " << A.ion[s].data->a[j].name << " on " << ACi << ", bnd2par " << k << " bnd2des " << m << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
									}

                                	out << setw(7) << left << count << " "
                                    	<< setw(65) << left << A.ion[s].molesmi << " "
                                    	<< setw(17) << left << ACi << " "
                                    	<< setw(17) << left << A.ion[s].data->a[j].id << " "
                                    	<< setw(30) << left << A.ion[s].data->a[j].name << " "
                                    	<< setw(17) << left << k << " "
                                    	<< setw(17) << left << m << endl;

		                            string mm=A.ion[s].molesmi;
        		                    for (int q=0;q<mm.length();q++) {
                		                if (mm[q]=='/') mm[q]='u';
                        		        if (mm[q]=='\\') mm[q]='d';
                            		}

                                	ofstream outs((para.smidir+"mds/"+mm+"_exc_IL.enc").c_str());
                                	A.ion[s].print(outs);
                                	outs.close();

                                	outs.open((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
									outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+mm+"_exc_IL.enc") << endl;
                                	outs.close();
                            	}
                            	inf.close();

                            	if (0) {
                                	if (s==0) A.ion[s].smiles=smi1_cat;
                                	if (s==1) A.ion[s].smiles=smi1_an;
                                	A.ion[s].input();
                            	}
                            	if (1) {
                                	A.ion[s].replace(buf[0].ion[s]);
                                	A.ion[s].reset();
                            	}

						
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


void exhaustive_cyclization1(IL &A) {
	int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);
	
	ofstream out((para.smidir+"cyclization_IL.txt").c_str());
	mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.print(out);
    out << endl;

	for (int s=0;s<2;s++) {
		if (s==0) {
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_cat_cyc" << " "
		        << setw(17) << left << "cat_cyc_pt1" << " "
		        << setw(17) << left << "cat_cyc_pt2" << endl;		
		}
		if (s==1) {
		    out << setw(7) << left << "No." << " "
		        << setw(65) << left << "IL_an_cyc" << " "
		        << setw(17) << left << "an_cyc_pt1" << " "
		        << setw(17) << left << "an_cyc_pt2" << endl;		
		}
	
	
		for (i=0;i<A.ion[s].Cindex.size();i++) {
			for (j=i+1;j<A.ion[s].Cindex.size();j++) {
				int ACi=A.ion[s].Cindex.at(i);
				int ACj=A.ion[s].Cindex.at(j);

				if (A.ion[s].ring(i,j)) {
					A.ion[s].mds2smi();
					if (0) {
        				A.ion[s].smiles=A.ion[s].molesmi;
						A.ion[s].input();			
					}
					if (1) A.ion[s].canonicalize_SMILES();
					
					A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
					//out<<A.molesmi<<endl;

                    //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                    system(("grep -Fw \""+A.ion[s].molesmi+"\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    ifstream inf("./tmp1");
                    inf >> ws;

					if (inf.eof()) {
						count++;
						if (s==0) cout << "EXHAUSTIVE CYCLIZATION: " << (smi1_cat) << " {CAT: "<< i << " , " << j << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
						if (s==1) cout << "EXHAUSTIVE CYCLIZATION: " << (smi1_an) << " {AN: "<< i << " , " << j << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
						
					    out << setw(7) << left << count << " "
					        << setw(65) << left << A.ion[s].molesmi << " "
					        << setw(17) << left << ACi << " "
					        << setw(17) << left << ACj << endl;

                            string mm=A.ion[s].molesmi;
                            for (int q=0;q<mm.length();q++) {
                                if (mm[q]=='/') mm[q]='u';
                                if (mm[q]=='\\') mm[q]='d';
                            }
						ofstream outs((para.smidir+"mds/"+mm+"_cyc_IL.enc").c_str());
						A.ion[s].print(outs);
						outs.close();

                        outs.open((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
                        outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+mm+"_cyc_IL.enc") << endl;
                        outs.close();

					}
					inf.close();
					
					if (0) {
						if (s==0) A.ion[s].smiles=smi1_cat;
						if (s==1) A.ion[s].smiles=smi1_an;
						A.ion[s].input();	
					}	
					if (1) {
						A.ion[s].replace(buf[0].ion[s]);			
						A.ion[s].reset();
					}
				}

			}
		}		
	}

	
	out.close();

	delete [] buf;
	buf=NULL;

}

void exhaustive_changect1(IL &A) {
    int count=0,i=0,j=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.smidir+"changect_IL.txt").c_str());
    mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.print(out);
    out << endl;
    for (int s=0;s<2;s++) {
        if (s==0) {
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_cat_ct" << " "
                << setw(17) << left << "cat_ct_pt" << " "
                << setw(17) << left << "cat_ct_h_ori" << " " 
				<< setw(17) << left << "cat_ct_h_after" << " " 
				<< setw(17) << left << "cat_ct_e_ori" << " "
				<< setw(17) << left << "cat_ct_e_after" << endl;
        }
        if (s==1) {
            out << setw(7) << left << "No." << " "
                << setw(65) << left << "IL_an_ct" << " "
                << setw(17) << left << "an_ct_pt" << " "
				<< setw(17) << left << "an_ct_h_ori" << " "
                << setw(17) << left << "an_ct_h_after" << " "
                << setw(17) << left << "an_ct_e_ori" << " "
                << setw(17) << left << "an_ct_e_after" << endl;
        }

		for (int w=0;w<2;w++) {
        	for (i=0;i<A.ion[s].ctsisomer.at(w).size();i++) {
                int ACi=A.ion[s].Cindex.at(i);
                string ACTh=A.ion[s].ctsisomer.at(0).at(i);
				string ACTe=A.ion[s].ctsisomer.at(1).at(i);

                if (A.ion[s].change_cistrans(i,w)) {
                    A.ion[s].mds2smi();
                    if (0) {
                        A.ion[s].smiles=A.ion[s].molesmi;
                        A.ion[s].input();
                    }
                    if (1) A.ion[s].canonicalize_SMILES();

                    A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
                    //out<<A.molesmi<<endl;

                    //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                    system(("grep -Fw \""+A.ion[s].molesmi+"\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    ifstream inf("./tmp1");
                    inf >> ws;
                    if (inf.eof()) {
                        count++;
                        if (s==0) {
							if (w==0) cout << "EXHAUSTIVE CHANGE OF CIS/TRANS: " << (smi1_cat) << " {CAT: "<< i << " | h | from " << ACTh << " to " << A.ion[s].ctsisomer.at(0).at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							else if (w==1) cout << "EXHAUSTIVE CHANGE OF CIS/TRANS: " << (smi1_cat) << " {CAT: "<< i << " | e | from " << ACTe << " to " << A.ion[s].ctsisomer.at(1).at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
						}
                        if (s==1) {
							if (w==0) cout << "EXHAUSTIVE CHANGE OF CIS/TRANS: " << (smi1_an) << " {AN: "<< i << " | h | from " << ACTh << " to " << A.ion[s].ctsisomer.at(0).at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
							else if (w==1) cout << "EXHAUSTIVE CHANGE OF CIS/TRANS: " << (smi1_an) << " {AN: "<< i << " | e | from " << ACTe << " to " << A.ion[s].ctsisomer.at(1).at(i) << " } -> " << A.ion[s].molesmi << " {No. " << count << " }" << endl;
						}

                        out << setw(7) << left << count << " "
                            << setw(65) << left << A.ion[s].molesmi << " "
                            << setw(17) << left << ACi << " "
                            << setw(17) << left << ACTh << " "
							<< setw(17) << left << A.ion[s].ctsisomer.at(0).at(i) << " " 
							<< setw(17) << left << ACTe << " " 
							<< setw(17) << left << A.ion[s].ctsisomer.at(1).at(i) << endl;

                            string mm=A.ion[s].molesmi;
                            for (int q=0;q<mm.length();q++) {
                                if (mm[q]=='/') mm[q]='u';
                                if (mm[q]=='\\') mm[q]='d';
                            }
                        ofstream outs((para.smidir+"mds/"+mm+"_ct_IL.enc").c_str());
                        A.ion[s].print(outs);
                        outs.close();

                        outs.open((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
                        outs << setw(50) << left << A.ion[s].molesmi << "   " << left << (para.smidir+"mds/"+mm+"_ct_IL.enc") << endl;
                        outs.close();

                    }
                    inf.close();

                    if (0) {
                        if (s==0) A.ion[s].smiles=smi1_cat;
                        if (s==1) A.ion[s].smiles=smi1_an;
                        A.ion[s].input();
                    }
                    if (1) {
                        A.ion[s].replace(buf[0].ion[s]);
                        A.ion[s].reset();
                    }
                }

        	}
		}
    }


    out.close();

    delete [] buf;
    buf=NULL;



} 
