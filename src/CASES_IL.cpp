#include "CASES_IL.h"
using namespace std;

extern PARAMETER para;

unsigned int exhaustive_combination(IL &A, IL &B) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;
	string smi2_cat=B.ion[0].molesmi,smi2_an=B.ion[1].molesmi;
	
	IL *buf=new IL [2];
	buf[0].replace(A);
	buf[1].replace(B);
	
	ofstream out((para.programdir+"combination_IL.txt").c_str());
	mark(out);
	
	out << "IL1: " << (smi1_cat+"."+smi1_an)  << endl;
	A.printmds(out);
	out << endl;
	out << "IL2: " << (smi2_cat+"."+smi2_an)  << endl;
	B.printmds(out);
	out << endl;

	out << setw(7) << left << "No." << " "
		<< setw(65) << left << "IL_comb" << " "
		<< setw(17) << left << "IL1_cat_comb_pt" << " "
		<< setw(17) << left << "IL1_an_comb_pt" << " "
		<< setw(17) << left << "bnd_cat" << " "
        << setw(17) << left << "IL2_cat_comb_pt" << " "
        << setw(17) << left << "IL2_an_comb_pt" << " "
		<< setw(17) << left << "bnd_an" << endl;

	vector< vector<bool> > cat_taboo(A.ion[0].Cindex.size(),vector<bool> (B.ion[0].Cindex.size(),0));
	vector< vector<bool> > an_taboo(A.ion[1].Cindex.size(),vector<bool> (B.ion[1].Cindex.size(),0));

	for (i=0;i<A.ion[0].Cindex.size();i++) {
    	for (j=0;j<B.ion[0].Cindex.size();j++) {
    		for (unsigned int i1=0;i1<A.ion[1].Cindex.size();i1++) {
    			for (unsigned int j1=0;j1<B.ion[1].Cindex.size();j1++) {
					unsigned int ACi=A.ion[0].Cindex.at(i);
					unsigned int ACi1=A.ion[1].Cindex.at(i1);
					unsigned int BCj=B.ion[0].Cindex.at(j);
					unsigned int BCj1=B.ion[1].Cindex.at(j1);

					for (unsigned int b=1;b<4;b++) {
					for (unsigned int b1=1;b1<4;b1++) {
					bool go1=0;
					if (!cat_taboo.at(i).at(j)) {
						if (A.ion[0].combination(B.ion[0],i,j,b)) go1=1;
						else cat_taboo.at(i).at(j)=1;
					}
					//else cat_taboo.at(i).at(j)=1;

					if (!an_taboo.at(i1).at(j1)) {
						if (A.ion[1].combination(B.ion[1],i1,j1,b1)) go1=1;
						else an_taboo.at(i1).at(j1)=1;
					}
					//else an_taboo.at(i1).at(j1)=1;

		        	if (go1) {
						
		        		A.mds2smi();
						if (0) {
		            		A.ion[0].smiles=A.ion[0].molesmi;
		            		A.ion[1].smiles=A.ion[1].molesmi;
							A.input();
						}
						if (1) A.canonicalize_SMILES();

						A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
						//out<<A.molesmi<<endl;
			
						//system(("ls "+para.programdir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
						system(("grep -Fw \""+A.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
						ifstream inf("./tmp1");
						inf >> ws;

						if (inf.eof()) {
							count++;
							cout << "EXHAUSTIVE COMBINATION: " << (smi1_cat+"."+smi1_an) << " { "<<i<<" , "<<i1<<" } + " << (smi2_cat+"."+smi2_an) << " { "<<j<<" , "<<j1<<" } , bnd { " << b << " , " << b1 << " } -> " << A.molesmi << " {No. " << count << " }" << endl;
							out << setw(7) << left << count << " "
								<< setw(65) << left << A.molesmi << " "
								<< setw(17) << left << ACi << " "
								<< setw(17) << left << ACi1 << " "
								<< setw(17) << left << b << " "
								<< setw(17) << left << BCj << " "
								<< setw(17) << left << BCj1 << " "
								<< setw(17) << left << b1 << endl;
							ofstream outs((para.programdir+"mds/"+A.molesmi+"_comb_IL.enc").c_str());						
							A.printmds(outs);
							outs.close();

							outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
							outs << setw(50) << left << A.molesmi << "   " << left << (para.programdir+"mds/"+A.molesmi+"_comb_IL.enc") << endl;
							outs.close();
						}
						inf.close();
		        
						if (0) {
							A.ion[0].smiles=smi1_cat;
							A.ion[1].smiles=smi1_an;
							A.input();
						}
						if (1) {
							A.replace(buf[0]);
							A.reset();
						}
		
						if (0) {
							B.ion[0].smiles=smi2_cat;
							B.ion[1].smiles=smi2_an;
		                	B.input();
						}
						if (1) {
							B.replace(buf[1]);
							B.reset();
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

	vector< vector<bool> >().swap(cat_taboo);
	vector< vector<bool> >().swap(an_taboo);

	return 1;
}

unsigned int exhaustive_crossover(IL &A, IL &B) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;
	string smi2_cat=B.ion[0].molesmi,smi2_an=B.ion[1].molesmi;
	
    IL *buf=new IL [2];
    buf[0].replace(A);
    buf[1].replace(B);

	ofstream out((para.programdir+"crossover_IL.txt").c_str());
	mark(out);

    out << "IL1: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;
    out << "IL2: " << (smi2_cat+"."+smi2_an)  << endl;
    B.printmds(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "IL_cros" << " "
        << setw(17) << left << "IL1_cat_cros_pt" << " "
        << setw(17) << left << "IL1_an_cros_pt" << " "
        << setw(17) << left << "IL2_cat_cros_pt" << " "
        << setw(17) << left << "IL2_an_cros_pt" << endl;

    vector< vector<bool> > cat_taboo(A.ion[0].Cindex.size(),vector<bool> (B.ion[0].Cindex.size(),0));
    vector< vector<bool> > an_taboo(A.ion[1].Cindex.size(),vector<bool> (B.ion[1].Cindex.size(),0));

	for (i=1;i<A.ion[0].Cindex.size();i++) {
		for (j=1;j<B.ion[0].Cindex.size();j++) {
			for (unsigned int i1=1;i1<A.ion[1].Cindex.size();i1++) {
				for (unsigned int j1=1;j1<B.ion[1].Cindex.size();j1++) {
                    unsigned int ACi=A.ion[0].Cindex.at(i);
                    unsigned int ACi1=A.ion[1].Cindex.at(i1);
                    unsigned int BCj=B.ion[0].Cindex.at(j);
                    unsigned int BCj1=B.ion[1].Cindex.at(j1);

                    bool go1=0;

					if (!cat_taboo.at(i).at(j)) {
                    	if (A.ion[0].crossover(B.ion[0],i,j)) go1=1;
						else cat_taboo.at(i).at(j)=1;
					}
					//else cat_taboo.at(i).at(j)=1;

					if (!an_taboo.at(i1).at(j1)) {
                    	if (A.ion[1].crossover(B.ion[1],i1,j1)) go1=1;
						else an_taboo.at(i1).at(j1)=1;
					}
					//else an_taboo.at(i1).at(j1)=1;

					if (go1) {

						A.mds2smi();
						if (0) {
		            		A.ion[0].smiles=A.ion[0].molesmi;
		            		A.ion[1].smiles=A.ion[1].molesmi; 
							A.input();
						}
						if (1) A.canonicalize_SMILES();

						B.mds2smi();
						if (0) {
		            		B.ion[0].smiles=B.ion[0].molesmi;
		            		B.ion[1].smiles=B.ion[1].molesmi; 
							B.input();
						}
						if (1) B.canonicalize_SMILES();
						
						bool go=0;

						A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
						//out<<A.molesmi<<endl;

                        //system(("ls "+para.programdir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                        system(("grep -Fw \""+A.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                        ifstream inf("./tmp1");
                        inf >> ws;

						if (inf.eof()) {
							ofstream outs((para.programdir+"mds/"+A.molesmi+"_cros_IL.enc").c_str());
							A.printmds(outs);
							outs.close();
							go=1;

                            outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
                            outs << setw(50) << left << A.molesmi << "   " << left << (para.programdir+"mds/"+A.molesmi+"_cros_IL.enc") << endl;
                            outs.close();
						}
						inf.close();
						
						B.molesmi=B.ion[0].molesmi+"."+B.ion[1].molesmi;
						//out<<B.molesmi<<endl;

                        //system(("ls "+para.programdir+"mds/ | grep -Fw \""+B.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                        system(("grep -Fw \""+B.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                        inf.open("./tmp1");
                        inf >> ws;

						if (inf.eof()) {
							ofstream outs((para.programdir+"mds/"+B.molesmi+"_cros_IL.enc").c_str());
							B.printmds(outs);
							outs.close();
							go=1;

                            outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
                            outs << setw(50) << left << B.molesmi << "   " << left << (para.programdir+"mds/"+B.molesmi+"_cros_IL.enc") << endl;
                            outs.close();
						}
						inf.close();
		
						if (go) {
							count++;
							cout << "EXHAUSTIVE CROSSOVER: " << (smi1_cat+"."+smi1_an) << " { "<<i<<" , "<<i1<<" } + " << (smi2_cat+"."+smi2_an) << " { "<<j<<" , "<<j1<<" } -> " << A.molesmi << " + " << B.molesmi << " {No. " << count << " }" << endl;
                            out << setw(7) << left << count << " "
                                << setw(65) << left << A.molesmi << " "
                                << setw(17) << left << ACi << " "
                                << setw(17) << left << ACi1 << " "
                                << setw(17) << left << BCj << " "
                                << setw(17) << left << BCj1 << endl;
						}
		
						if (0) {
							A.ion[0].smiles=smi1_cat;
							A.ion[1].smiles=smi1_an;
							A.input();	
						}
						if (1) {
							A.replace(buf[0]);
							A.reset();
						}
		
						if (0) {
							B.ion[0].smiles=smi2_cat;
							B.ion[1].smiles=smi2_an;
							B.input();
						}
						if (1) {
							B.replace(buf[1]);
							B.reset();
						}
						
					}				
				}
			}
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;

    vector< vector<bool> >().swap(cat_taboo);
    vector< vector<bool> >().swap(an_taboo);

	return 1;
}


unsigned int exhaustive_subtraction(IL &A) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);
	
	ofstream out((para.programdir+"subtraction_IL.txt").c_str());
	mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "IL_subtr" << " "
        << setw(17) << left << "cat_subtr_pt" << " "
		<< setw(17) << left << "cat_subtr_bndfrm" << " "
        << setw(17) << left << "an_subtr_pt" << " " 
		<< setw(17) << left << "an_subtr_bndfrm" << endl;

	for (i=0;i<A.ion[0].Cindex.size();i++) {
		for (unsigned int bd=0;bd<4;bd++) {
			for (j=0;j<A.ion[1].Cindex.size();j++) {
				for (unsigned int bd1=0;bd1<4;bd1++) {
					unsigned int ACi=A.ion[0].Cindex.at(i);
					unsigned int ACj=A.ion[1].Cindex.at(j);

		            bool go1=0;
		            if (A.ion[0].subtraction(i,bd,0)) go1=1;
		            if (A.ion[1].subtraction(j,bd1,0)) go1=1;

					if (go1) {
						A.mds2smi();
						if (0) {
            				A.ion[0].smiles=A.ion[0].molesmi;
            				A.ion[1].smiles=A.ion[1].molesmi; 
							A.input();
						}
						if (1) A.canonicalize_SMILES();
				
						A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
						//out<<A.molesmi<<endl;

                		//system(("ls "+para.programdir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                		system(("grep -Fw \""+A.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                		ifstream inf("./tmp1");
                		inf >> ws;

						if (inf.eof()) {
							count++;
							cout << "EXHAUSTIVE SUBTRACTION: " << (smi1_cat+"."+smi1_an) << " { "<<i<<" , "<<j<<" } -> " << A.molesmi << " {No. " << count << " }" << endl;
                			out << setw(7) << left << count << " "
                    			<< setw(65) << left << A.molesmi << " "
                    			<< setw(17) << left << ACi << " "
								<< setw(17) << left << bd << " "
                    			<< setw(17) << left << ACj << " "
								<< setw(17) << left << bd1 << endl;

		    				ofstream outs((para.programdir+"mds/"+A.molesmi+"_subtr_IL.enc").c_str());
							A.printmds(outs);
							outs.close();

                    		outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
                    		outs << setw(50) << left << A.molesmi << "   " << left << (para.programdir+"mds/"+A.molesmi+"_subtr_IL.enc") << endl;
                    		outs.close();
						}
						inf.close();
	
						if (0) {			
							A.ion[0].smiles=smi1_cat;
							A.ion[1].smiles=smi1_an;
							A.input();
						}
						if (1) {
							A.replace(buf[0]);
							A.reset();
						}
					}
				}
			}
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;

	return 1;
}


unsigned int exhaustive_addition(IL &A) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

	ofstream out((para.programdir+"addition_IL.txt").c_str());
	mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "IL_add" << " "
        << setw(17) << left << "cat_add_pt" << " "
		<< setw(17) << left << "cat_add_id" << " "
		<< setw(30) << left << "cat_add_name" << " "
        << setw(17) << left << "an_add_pt" << " "
		<< setw(17) << left << "an_add_id" << " "
        << setw(30) << left << "an_add_name" << endl;

    vector< vector< vector<bool> > > cat_taboo(A.ion[0].Cindex.size(),vector< vector<bool> > (A.ion[0].data->num, vector<bool> (4,0)));
    vector< vector< vector<bool> > > an_taboo(A.ion[1].Cindex.size(),vector< vector<bool> > (A.ion[1].data->num, vector<bool> (4,0)));

	for (i=0;i<A.ion[0].Cindex.size();i++) {
		for (unsigned int i1=0;i1<A.ion[1].Cindex.size();i1++) {
			for (j=1;j<A.ion[0].data->num;j++) {
				if (A.ion[0].data->a[j].probability>0){
					for (unsigned int j1=1;j1<A.ion[1].data->num;j1++) {
						if (A.ion[1].data->a[j1].probability>0) {
            				unsigned int ACi=A.ion[0].Cindex.at(i);
            				unsigned int ACi1=A.ion[1].Cindex.at(i1);

							for (unsigned int b=1;b<4;b++) {
							for (unsigned int b1=1;b1<4;b1++) {
            				bool go1=0;
							if (!cat_taboo.at(i).at(j).at(b)) {
            					if (A.ion[0].addition(i,j,b)) go1=1;
								else cat_taboo.at(i).at(j).at(b)=1;
							}
							//else cat_taboo.at(i).at(j).at(b)=1;

							if (!an_taboo.at(i1).at(j1).at(b1)) {
            					if (A.ion[1].addition(i1,j1,b1)) go1=1;
								else an_taboo.at(i1).at(j1).at(b1)=1;
                            }
                            //else an_taboo.at(i1).at(j1).at(b1)=1;


							if (go1) {
								A.mds2smi();
								if (0) {
				            		A.ion[0].smiles=A.ion[0].molesmi;
				            		A.ion[1].smiles=A.ion[1].molesmi; 
									A.input();
								}
								if (1) A.canonicalize_SMILES();
								
								A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
								//out<<A.molesmi<<endl;

                        		//system(("ls "+para.programdir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                        		system(("grep -Fw \""+A.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                				ifstream inf("./tmp1");
                				inf >> ws;

								if (inf.eof()) {
									count++;
									cout << "EXHAUSTIVE ADDITION: " << (smi1_cat+"."+smi1_an) << " { "<< j << " " << A.ion[0].data->a[j].name << " , " << j1 << " " << A.ion[1].data->a[j1].name << " on " <<i<<" , "<<j<<" } -> " << A.molesmi << " {No. " << count << " }" << endl;

    								out << setw(7) << left << count << " "
        								<< setw(65) << left << A.molesmi << " "
								        << setw(17) << left << ACi << " "
								        << setw(17) << left << A.ion[0].data->a[j].id << " "
								        << setw(30) << left << A.ion[0].data->a[j].name << " "
								        << setw(17) << left << ACi1 << " "
								        << setw(17) << left << A.ion[1].data->a[j1].id << " "
								        << setw(30) << left << A.ion[1].data->a[j1].name << endl;

									ofstream outs((para.programdir+"mds/"+A.molesmi+"_add_IL.enc").c_str());
									A.printmds(outs);
									outs.close();

		                            outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
        		                    outs << setw(50) << left << A.molesmi << "   " << left << (para.programdir+"mds/"+A.molesmi+"_add_IL.enc") << endl;
                		            outs.close();

								}
								inf.close();
								
								if (0) {
									A.ion[0].smiles=smi1_cat;
									A.ion[1].smiles=smi1_an;
									A.input();	
								}
								if (1) {
									A.replace(buf[0]);
									A.reset();
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
	
	out.close();

	delete [] buf;
	buf=NULL;

    vector< vector< vector<bool> > >().swap(cat_taboo);
    vector< vector< vector<bool> > >().swap(an_taboo);

	return 1;
}

unsigned int exhaustive_insertion(IL &A) {
    unsigned int count=0,i=0,j=0,k=0,m=0;
    string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);

    ofstream out((para.programdir+"insertion_IL.txt").c_str());
    mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "IL_ins" << " "
        << setw(17) << left << "cat_ins_pt" << " "
        << setw(17) << left << "cat_ins_id" << " "
        << setw(30) << left << "cat_ins_name" << " "
        << setw(17) << left << "cat_bnd2par" << " "
        << setw(17) << left << "cat_bnd2des" << " "
        << setw(17) << left << "an_ins_pt" << " "
        << setw(17) << left << "an_ins_id" << " "
        << setw(30) << left << "an_ins_name" << " "
        << setw(17) << left << "an_bnd2par" << " "
        << setw(17) << left << "an_bnd2des" << endl;

    vector< vector< vector<bool> > > cat_taboo(A.ion[0].Cindex.size(),vector< vector<bool> > (A.ion[0].data->num,vector<bool> (3,0)));
    vector< vector< vector<bool> > > an_taboo(A.ion[1].Cindex.size(),vector< vector<bool> > (A.ion[1].data->num,vector<bool> (3,0)));

    for (i=1;i<A.ion[0].Cindex.size();i++) {
        for (j=1;j<A.ion[0].data->num;j++) {
            if (A.ion[0].data->a[j].probability>0) {
                for (k=1;k<4;k++) {
                    for (m=1;m<4;m++) {
                        if (A.ion[0].data->a[m].probability>0) {
                            for (unsigned int i1=1;i1<A.ion[1].Cindex.size();i1++) {
                                for (unsigned int j1=1;j<A.ion[1].data->num;j1++) {
                                    if (A.ion[1].data->a[j1].probability>0) {
                                        for (unsigned int k1=1;k1<4;k1++) {
                                            for (unsigned int m1=1;m1<4;m++) {
                                                if (A.ion[1].data->a[m1].probability>0) {
                                                    unsigned int ACi=A.ion[0].Cindex.at(i);
                                                    unsigned int ACi1=A.ion[1].Cindex.at(i1);

                                                    bool go1=0;
                                                    if (!cat_taboo.at(i).at(j).at(m-1)) {
                                                        if (A.ion[0].insertion(i,j,k,m)) go1=1;
                                                        else cat_taboo.at(i).at(j).at(m-1)=1;
                                                    }

                                                    if (!an_taboo.at(i1).at(j1).at(m1-1)) {
                                                        if (A.ion[1].insertion(i1,j1,k1,m1)) go1=1;
                                                        else an_taboo.at(i1).at(j1).at(m1-1)=1;
                                                    }

                                                    if (go1) {
                                                        A.mds2smi();
                                                        if (0) {
                                                            A.ion[0].smiles=A.ion[0].molesmi;
                                                            A.ion[1].smiles=A.ion[1].molesmi;
                                                            A.input();
                                                        }
                                                        if (1) A.canonicalize_SMILES();

                                                        A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
                                                        //out<<A.molesmi<<endl;

                                                        //system(("ls "+para.programdir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                                                        system(("grep -Fw \""+A.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                                                        ifstream inf("./tmp1");
                                                        inf >> ws;

                                                        if (inf.eof()) {
                                                            count++;
                                                            cout << "EXHAUSTIVE EXCHANGE: " << (smi1_cat+"."+smi1_an) << " {CAT: " << j << " " << A.ion[0].data->a[j].name << " , bnd2par " << k << " , bnd2des " << m << " ; "
                                                                                                                        << " AN: " << j1 << " " << A.ion[1].data->a[j1].name << " , bnd2par " << k << " , bnd2des " << m1 << " } -> " << A.molesmi << " {No. " << count << " }" << endl;

                                                            out << setw(7) << left << count << " "
                                                                << setw(65) << left << A.molesmi << " "
                                                                << setw(17) << left << ACi << " "
                                                                << setw(17) << left << A.ion[0].data->a[j].id << " "
                                                                << setw(30) << left << A.ion[0].data->a[j].name << " "
                                                                << setw(17) << left << k << " "
                                                                << setw(17) << left << m << " "
                                                                << setw(17) << left << ACi1 << " "
                                                                << setw(17) << left << A.ion[1].data->a[j1].id << " "
                                                                << setw(30) << left << A.ion[1].data->a[j1].name << " "
                                                                << setw(17) << left << k1 << " "
                                                                << setw(17) << left << m1 << endl;

                                                            ofstream outs((para.programdir+"mds/"+A.molesmi+"_ins_IL.enc").c_str());
                                                            A.printmds(outs);
                                                            outs.close();

                                                            outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
                                                            outs << setw(50) << left << A.molesmi << "   " << left << (para.programdir+"mds/"+A.molesmi+"_ins_IL.enc") << endl;
                                                            outs.close();
                                                        }
                                                        inf.close();
                                                        if (0) {
                                                            A.ion[0].smiles=smi1_cat;
                                                            A.ion[1].smiles=smi1_an;
                                                            A.input();
                                                        }
                                                        if (1) {
                                                            A.replace(buf[0]);
                                                            A.reset();
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
                }
            }
        }
    }

    out.close();

    delete [] buf;
    buf=NULL;

    vector< vector< vector<bool> > >().swap(cat_taboo);
    vector< vector< vector<bool> > >().swap(an_taboo);

	return 1;
}



unsigned int exhaustive_change_ele(IL &A) {
	unsigned int count=0,i=0,j=0,k=0,m=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;
	
    IL *buf=new IL [1];
    buf[0].replace(A);

	ofstream out((para.programdir+"change_ele_IL.txt").c_str());
	mark(out);

    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;

	out << setw(7) << left << "No." << " "
		<< setw(65) << left << "IL_chele" << " "
		<< setw(17) << left << "cat_chele_pt" << " "
		<< setw(17) << left << "cat_chele_id" << " "
		<< setw(30) << left << "cat_chele_name" << " "
		<< setw(17) << left << "cat_bnd2par" << " "
		<< setw(17) << left << "cat_bnd2des" << " "
        << setw(17) << left << "an_chele_pt" << " "
        << setw(17) << left << "an_chele_id" << " "
        << setw(30) << left << "an_chele_name" << " "
        << setw(17) << left << "an_bnd2par" << " "
        << setw(17) << left << "an_bnd2des" << endl;

    vector< vector< vector<bool> > > cat_taboo(A.ion[0].Cindex.size(),vector< vector<bool> > (A.ion[0].data->num,vector<bool> (3,0)));
    vector< vector< vector<bool> > > an_taboo(A.ion[1].Cindex.size(),vector< vector<bool> > (A.ion[1].data->num,vector<bool> (3,0)));

	for (i=1;i<A.ion[0].Cindex.size();i++) {
		for (j=1;j<A.ion[0].data->num;j++) {
			if (A.ion[0].data->a[j].probability>0) {
				for (k=1;k<4;k++) {
					for (m=1;m<4;m++) {
						if (A.ion[0].data->a[m].probability>0) {
							for (unsigned int i1=1;i1<A.ion[1].Cindex.size();i1++) {
								for (unsigned int j1=1;j<A.ion[1].data->num;j1++) {
									if (A.ion[1].data->a[j1].probability>0) {
										for (unsigned int k1=1;k1<4;k1++) {
											for (unsigned int m1=1;m1<4;m++) {
												if (A.ion[1].data->a[m1].probability>0) {
													unsigned int ACi=A.ion[0].Cindex.at(i);
													unsigned int ACi1=A.ion[1].Cindex.at(i1);

                            						bool go1=0;
													if (!cat_taboo.at(i).at(j).at(m-1)) {
                            							if (A.ion[0].change_ele(i,j,k,m)) go1=1;
														else cat_taboo.at(i).at(j).at(m-1)=1;
													}

													if (!an_taboo.at(i1).at(j1).at(m1-1)) {
	                            						if (A.ion[1].change_ele(i1,j1,k1,m1)) go1=1;
                                                        else an_taboo.at(i1).at(j1).at(m1-1)=1;
                                                    }

													if (go1) {
														A.mds2smi();
														if (0) {
										            		A.ion[0].smiles=A.ion[0].molesmi;
										            		A.ion[1].smiles=A.ion[1].molesmi; 
															A.input();
														}
														if (1) A.canonicalize_SMILES();
														
														A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
														//out<<A.molesmi<<endl;

                        								//system(("ls "+para.programdir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                        								system(("grep -Fw \""+A.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                            							ifstream inf("./tmp1");
                                						inf >> ws;

														if (inf.eof()) {
															count++;
															cout << "EXHAUSTIVE CHANGE_ELE: " << (smi1_cat+"."+smi1_an) << " {CAT: " << j << " " << A.ion[0].data->a[j].name << " , bnd2par " << k << " , bnd2des " << m << " ; " 
																														<< " AN: " << j1 << " " << A.ion[1].data->a[j1].name << " , bnd2par " << k << " , bnd2des " << m1 << " } -> " << A.molesmi << " {No. " << count << " }" << endl;

						                                	out << setw(7) << left << count << " "
						                                    	<< setw(65) << left << A.molesmi << " "
						                                    	<< setw(17) << left << ACi << " "
						                                    	<< setw(17) << left << A.ion[0].data->a[j].id << " "
						                                    	<< setw(30) << left << A.ion[0].data->a[j].name << " "
						                                    	<< setw(17) << left << k << " "
						                                    	<< setw(17) << left << m << " "
						                                    	<< setw(17) << left << ACi1 << " "
						                                    	<< setw(17) << left << A.ion[1].data->a[j1].id << " "
						                                    	<< setw(30) << left << A.ion[1].data->a[j1].name << " "
						                                    	<< setw(17) << left << k1 << " "
						                                    	<< setw(17) << left << m1 << endl;

															ofstream outs((para.programdir+"mds/"+A.molesmi+"_chele_IL.enc").c_str());
															A.printmds(outs);
															outs.close();

								                            outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
								                            outs << setw(50) << left << A.molesmi << "   " << left << (para.programdir+"mds/"+A.molesmi+"_chele_IL.enc") << endl;
								                            outs.close();
														}
														inf.close();
					
														if (0) {									
															A.ion[0].smiles=smi1_cat;
															A.ion[1].smiles=smi1_an;
															A.input();	
														}
														if (1) {
															A.replace(buf[0]);
															A.reset();
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
				}
			}
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;

	vector< vector< vector<bool> > >().swap(cat_taboo);
	vector< vector< vector<bool> > >().swap(an_taboo);

	return 1;
}

unsigned int exhaustive_cyclization(IL &A) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;

    IL *buf=new IL [1];
    buf[0].replace(A);
	
	ofstream out((para.programdir+"cyclization_IL.txt").c_str());
	mark(out);


    out << "IL: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "IL_cyc" << " "
        << setw(17) << left << "cat_cyc_pt1" << " "
        << setw(17) << left << "cat_cyc_pt2" << " "
        << setw(17) << left << "an_cyc_pt1" << " "
        << setw(17) << left << "an_cyc_pt2" << endl;

    vector< vector<bool> > cat_taboo(A.ion[0].Cindex.size(),vector<bool> (A.ion[0].Cindex.size(),0));
    vector< vector<bool> > an_taboo(A.ion[1].Cindex.size(),vector<bool> (A.ion[1].Cindex.size(),0));

	for (i=0;i<A.ion[0].Cindex.size();i++) {
		for (j=i+1;j<A.ion[0].Cindex.size();j++) {
			for (unsigned int k=1;k<=3;k++) {
				for (unsigned int i1=0;i<A.ion[1].Cindex.size();i1++) {
					for (unsigned int j1=i1+1;j1<A.ion[1].Cindex.size();j1++) {
						for (unsigned int k1=1;k1<=3;k1++) {
							unsigned int ACi=A.ion[0].Cindex.at(i);
							unsigned int ACj=A.ion[0].Cindex.at(j);
							unsigned int ACi1=A.ion[1].Cindex.at(i1);
							unsigned int ACj1=A.ion[1].Cindex.at(j1);

							bool go1=0;
							if (!cat_taboo.at(i).at(j)) {
								if (A.ion[0].cyclization(i,j,k)) go1=1;
								else cat_taboo.at(i).at(j)=1;
							}
							//else cat_taboo.at(i).at(j)=1;

							if (!an_taboo.at(i1).at(j1)) {
								if (A.ion[1].cyclization(i1,j1,k1)) go1=1;
								else an_taboo.at(i1).at(j1)=1;
							}
							//else an_taboo.at(i1).at(j1)=1;

							if (go1) {
								A.mds2smi();
								if (0) {
									A.ion[0].smiles=A.ion[0].molesmi;
									A.ion[1].smiles=A.ion[1].molesmi; 
									A.input();			
								}
								if (1) A.canonicalize_SMILES();
								
								A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
								//out<<A.molesmi<<endl;

								//system(("ls "+para.programdir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
								system(("grep -Fw \""+A.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
								ifstream inf("./tmp1");
								inf >> ws;

								if (inf.eof()) {
									count++;
									cout << "EXHAUSTIVE CYCLIZATION: " << (smi1_cat+"."+smi1_an) << " {CAT: "<< j << " , " << i <<" | AN: "<< j1 << " , " << i1 <<" } -> " << A.molesmi << " {No. " << count << " }" << endl;

									out << setw(7) << left << count << " "
										<< setw(65) << left << A.molesmi << " "
										<< setw(17) << left << ACi << " "
										<< setw(17) << left << ACj << " "
										<< setw(17) << left << ACi1 << " "
										<< setw(17) << left << ACj1 << endl;

									ofstream outs((para.programdir+"mds/"+A.molesmi+"_cyc_IL.enc").c_str());
									A.printmds(outs);
									outs.close();

									outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
									outs << setw(50) << left << A.molesmi << "   " << left << (para.programdir+"mds/"+A.molesmi+"_cyc_IL.enc") << endl;
									outs.close();

								}
								inf.close();
								
								if (0) {
									A.ion[0].smiles=smi1_cat;
									A.ion[1].smiles=smi1_an;
									A.input();	
								}	
								if (1) {
									A.replace(buf[0]);			
									A.reset();
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

	vector< vector<bool> >().swap(cat_taboo);
	vector< vector<bool> >().swap(an_taboo);

	return 1;
}


unsigned int comp_swap(IL &A,IL &B) {
	unsigned int count=0,i=0,j=0;
	string smi1_cat=A.ion[0].molesmi,smi1_an=A.ion[1].molesmi;
	string smi2_cat=B.ion[0].molesmi,smi2_an=B.ion[1].molesmi;

    IL *buf=new IL [2];
    buf[0].replace(A);
	buf[1].replace(B);
	
	ofstream out((para.programdir+"switch_IL.txt").c_str());
	mark(out);

    out << "IL1: " << (smi1_cat+"."+smi1_an)  << endl;
    A.printmds(out);
    out << endl;
    out << "IL2: " << (smi2_cat+"."+smi2_an)  << endl;
    B.printmds(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "IL1_swap" << " "
		<< setw(65) << left << "IL2_swap"<< endl;

	swap(A.ion[1],B.ion[1]);
	A.reset();
	B.reset();
	A.mds2smi();
	B.mds2smi();

	bool go=0;
	
	if (0) {
		A.ion[0].smiles=A.ion[0].molesmi;
		A.ion[1].smiles=A.ion[1].molesmi; 
		A.input();	
	}
	if (1) A.canonicalize_SMILES();
	A.molesmi=A.ion[0].molesmi+"."+A.ion[1].molesmi;
	//out<<A.molesmi<<endl;
	
	if (0) {
		B.ion[0].smiles=B.ion[0].molesmi;
		B.ion[1].smiles=B.ion[1].molesmi; 
		B.input();
	}	
	if (1) B.canonicalize_SMILES();			
	B.molesmi=B.ion[0].molesmi+"."+B.ion[1].molesmi;	
	//out<<B.molesmi<<endl;

    //system(("ls "+para.programdir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
    system(("grep -Fw \""+A.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());	
    ifstream inf("./tmp1");
    inf >> ws;

    if (inf.eof()) {
        ofstream outs((para.programdir+"mds/"+A.molesmi+"_swap_IL.enc").c_str());
    	A.printmds(outs);
        outs.close();
        go=1;

        outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
        outs << setw(50) << left << A.molesmi << "   " << left << (para.programdir+"mds/"+A.molesmi+"_swap_IL.enc") << endl;
        outs.close();
    }
    inf.close();


    //system(("ls "+para.programdir+"mds/ | grep -Fw \""+B.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
    system(("grep -Fw \""+B.molesmi+"\" "+para.programdir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
    inf.open("./tmp1");
    inf >> ws;

    if (inf.eof()) {
        ofstream outs((para.programdir+"mds/"+B.molesmi+"_swap_IL.enc").c_str());
        B.printmds(outs);
        outs.close();
        go=1;

        outs.open((para.programdir+"mds/DATLIST.txt").c_str(),ios::app);
        outs << setw(50) << left << B.molesmi << "   " << left << (para.programdir+"mds/"+B.molesmi+"_swap_IL.enc") << endl;
        outs.close();
    }
    inf.close();


	if (go) {
		count++;
		cout << "SWITCH: " << (smi1_cat+"."+smi1_an) << " + " << (smi2_cat+"."+smi2_an) <<" -> " << A.molesmi << " + " << B.molesmi << " {No. " << count << " }" << endl;
    	out << setw(7) << left << count << " "
        	<< setw(65) << left << A.molesmi << " "
        	<< setw(65) << left << B.molesmi << endl;
	}

	if (0) {
		A.ion[0].smiles=smi1_cat;
		A.ion[1].smiles=smi1_an;
		A.input();						
	}
	if (1) {
		A.replace(buf[0]);
		A.reset();
	}

	if (0) {
		B.ion[0].smiles=smi2_cat;
		B.ion[1].smiles=smi2_an;
		B.input();		
	}
	if (1) {
		B.replace(buf[1]);
		B.reset();
	}
	
	out.close();

	delete [] buf;
	buf=NULL;

	return 1;
}

