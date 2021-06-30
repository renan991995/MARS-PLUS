#include "CASES_NEU.h"
using namespace std;

extern PARAMETER para;

void exhaustive_combination(MOLECULE &A, MOLECULE &B) {
	int count=0,i=0,j=0;
	string smi1=A.molesmi,smi2=B.molesmi;
	
	MOLECULE *buf=new MOLECULE [2];
	buf[0].replace(A);
	buf[1].replace(B);
	
	ofstream out("combination.txt",ios::app);
	//mark(out);

    out << "MOL1: " << smi1 << endl;
    A.print(out);
    out << endl;
    out << "MOL2: " << smi2 << endl;
    B.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_comb" << " "
        << setw(17) << left << "MOL1_comb_pt" << " "
        << setw(17) << left << "MOL2_comb_pt" << " "
        << setw(17) << left << "bnd" << endl;


	for (i=0;i<A.Cindex.size();i++) {
    	for (j=0;j<B.Cindex.size();j++) {
            int ACi=A.Cindex.at(i);
            int BCj=B.Cindex.at(j);
			for (int b=1;b<4;b++) {
        		if (A.combination(B,i,j,b)) {
        			//A.mds2smi();
					A.canonicalize_SMILES();
            		A.smiles=A.molesmi;   
				         	
					//A.input();
				
					//out<<A.molesmi<<endl;

                	string mm=A.molesmi;
                	for (int q=0;q<mm.length();q++) {
                    	if (mm[q]=='/') mm[q]='u';
                    	if (mm[q]=='\\') mm[q]='d';
						if (mm[q]=='*') mm[q]='x';
                	}

                	//system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                	system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str()); //A.ion[s].molesmi
					//system(("grep -Fw \""+mm+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                	ifstream inf("./tmp1");
                	inf >> ws;

                	if (inf.eof()) {
                    	count++;
						if (1) cout << "EXHAUSTIVE COMBINATION: " << smi1 << " { " << ACi << " } + " << smi2 << " { " << BCj << " } , bnd " << b << " -> " << A.molesmi << " {No. " << count << " }" << endl;

                    	out << setw(7) << left << count << " "
                        	<< setw(65) << left << A.molesmi << " "
                        	<< setw(17) << left << ACi << " "
                        	<< setw(17) << left << BCj << " "
                        	<< setw(17) << left << b << endl;

						if (para.ifwritemds) {
                			ofstream outs((para.smidir+"mds/"+mm+"_comb.enc").c_str());
							//ofstream outs((para.outdir+mm+"_comb.enc").c_str());
                			A.print(outs);
                			outs.close();
						}

            			ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
						//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
                		outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_comb.enc") << endl;
						//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_comb.enc") << endl;
                		outs.close();

						if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

						if (para.enumeration) SMI_Enumerator(A.molesmi);
					}
        
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

	para.stat << "NEU | COMBINATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
	
}

void exhaustive_crossover(MOLECULE &A, MOLECULE &B) {
	int count=0,i=0,j=0;
	A.canonicalize_SMILES();
	string smi1=A.molesmi,smi2=B.molesmi;
	
	MOLECULE *buf=new MOLECULE [2];
	buf[0].replace(A);	
	buf[1].replace(B);
	
	ofstream out("crossover.txt",ios::app);
	//mark(out);

    out << "MOL1: " << smi1 << endl;
    A.print(out);
    out << endl;
    out << "MOL2: " << smi2 << endl;
    B.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL1_cros" << " "
        << setw(65) << left << "MOL2_cros" << " "
        << setw(17) << left << "MOL1_cros_pt" << " "
        << setw(17) << left << "MOL2_cros_pt" << endl;

	for (i=1;i<A.Cindex.size();i++) {
		for (j=1;j<B.Cindex.size();j++) {
            int ACi=A.Cindex.at(i);
            int BCj=B.Cindex.at(j);

			if (A.crossover(B,i,j)) {
				//count++;
				//A.mds2smi();
				A.canonicalize_SMILES();
				A.smiles=A.molesmi;  
				//A.input();
				
				//B.mds2smi();
				B.canonicalize_SMILES();
				B.smiles=B.molesmi;  
				//B.input();
				
				bool go=0;

				//out<<A.molesmi<<endl;
				//out<<B.molesmi<<endl;

                string mm=A.molesmi;
                for (int q=0;q<mm.length();q++) {
                    if (mm[q]=='/') mm[q]='u';
                    if (mm[q]=='\\') mm[q]='d';
                    if (mm[q]=='*') mm[q]='x';
                }

                //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                system(("grep -Fw \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
				//system(("grep -Fw \""+A.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                ifstream inf("./tmp1");
                inf >> ws;

                if (inf.eof()) {

					if (para.ifwritemds) {
                		ofstream outs((para.smidir+"mds/"+mm+"_cros.enc").c_str());
						//ofstream outs((para.outdir+mm+"_cros.enc").c_str());
                		A.print(outs);
                		outs.close();
					}
					go=1;

                	ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
					//ofstream outs((para.outdir+"mds/DATLIST.txt").c_str(),ios::app);
                	outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_cros.enc") << endl;
					//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_cros.enc") << endl;
                	outs.close();

					if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

					if (para.enumeration) SMI_Enumerator(A.molesmi);
				}
				inf.close();

                mm=B.molesmi;
                for (int q=0;q<mm.length();q++) {
                    if (mm[q]=='/') mm[q]='u';
                    if (mm[q]=='\\') mm[q]='d';
                    if (mm[q]=='*') mm[q]='x';
                }

                //system(("ls "+para.smidir+"mds/ | grep -Fw \""+B.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                system(("grep -Fw \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
				//system(("grep -Fw \""+B.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                inf.open("./tmp1");
                inf >> ws;

                if (inf.eof()) {

					if (para.ifwritemds) {
                		ofstream outs((para.smidir+"mds/"+mm+"_cros.enc").c_str());
						//ofstream outs((para.outdir+mm+"_cros.enc").c_str());
                		B.print(outs);
                		outs.close();
					}

                	ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
					//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
                	outs << setw(50) << left << B.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_cros.enc") << endl;
					//outs << setw(50) << left << B.molesmi << "   " << left << (para.outdir+mm+"_cros.enc") << endl;
                	outs.close();

					if (para.glbouf.is_open()) para.glbouf << B.molesmi << endl;

					if (para.enumeration) SMI_Enumerator(B.molesmi);
                }
                inf.close();

                if (go) {
                    count++;
                    if (1) cout << "EXHAUSTIVE CROSSOVER: " << smi1 << " { " << ACi << " } + " << smi2 << " { " << BCj << " } -> " << A.molesmi << " + " << B.molesmi << " {No. " << count << " }" << endl;
                    out << setw(7) << left << count << " "
                        << setw(65) << left << A.molesmi << " "
                        << setw(65) << left << B.molesmi << " "
                        << setw(17) << left << ACi << " "
                        << setw(17) << left << BCj << endl;
                }


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
	
	para.stat << "NEU | CROSSOVER | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}


void exhaustive_subtraction(MOLECULE &A) {
	int count=0,i=0,j=0;
	string smi=A.molesmi;
	
	MOLECULE *buf=new MOLECULE [1];
	buf[0].replace(A);
	
	ofstream out("subtraction.txt",ios::app);
	//mark(out);

    out << "MOL: " << smi << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_subtr" << " "
        << setw(17) << left << "MOL_subtr_pt" <<  endl;


	for (i=0;i<A.Cindex.size();i++) {
		int ACi=A.Cindex.at(i);
		if (A.subtraction(i,1)) {
			//count++;
			//A.mds2smi();
			A.canonicalize_SMILES();
			A.smiles=A.molesmi;  
			//A.input();
			
			//out<<A.molesmi<<endl;

            string mm=A.molesmi;
            for (int q=0;q<mm.length();q++) {
                if (mm[q]=='/') mm[q]='u';
                if (mm[q]=='\\') mm[q]='d';
                if (mm[q]=='*') mm[q]='x';
            }

            //system(("ls "+para.smidir+"mds/ | grep -F \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
            system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
			//system(("grep -Fw \""+A.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
            ifstream inf("./tmp1");
            inf >> ws;

            if (inf.eof()) {
                count++;
                if (1) cout << "EXHAUSTIVE SUBTRACTION: " << smi << " { " << ACi << " } -> " << A.molesmi << " {No. " << count << " }" << endl;

                out << setw(7) << left << count << " "
                    << setw(65) << left << A.molesmi << " "
                    << setw(17) << left << ACi << endl;
				if (para.ifwritemds) {
            		ofstream outs((para.smidir+"mds/"+mm+"_subtr.enc").c_str());
					//ofstream outs((para.outdir+mm+"_subtr.enc").c_str());
           			A.print(outs);
            		outs.close();
				}

           		ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
				//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
            	outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_subtr.enc") << endl;
				//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_subtr.enc") << endl;
            	outs.close();

				if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

				if (para.enumeration) SMI_Enumerator(A.molesmi);
            }
            inf.close();

			
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

	para.stat << "NEU | SUBTRACTION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}


void exhaustive_addition(MOLECULE &A) {
	int count=0,i=0,j=0;
	string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);

	
	ofstream out("addition.txt",ios::app);
	//mark(out);

    out << "MOL: " << smi << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_add" << " "
        << setw(17) << left << "MOL_add_pt" << " "
        << setw(17) << left << "Add_id" << " "
        << setw(30) << left << "Add_name" << endl;

	for (i=0;i<A.Cindex.size();i++) {
		int ACi=A.Cindex.at(i);
		for (j=1;j<A.data->num;j++) {
			if (A.data->a[j].probability>0) {
				for (int b=1;b<4;b++) {
					if (A.addition(i,j,b)) {
						//count++;
						//A.mds2smi();
						A.canonicalize_SMILES();
						A.smiles=A.molesmi;  
						//A.input();
				
						//out<<A.molesmi<<endl;

                        string mm=A.molesmi;
                        for (int q=0;q<mm.length();q++) {
                            if (mm[q]=='/') mm[q]='u';
                            if (mm[q]=='\\') mm[q]='d';
                            if (mm[q]=='*') mm[q]='x';
                        }

                    	//system(("ls "+para.smidir+"mds/ | grep -F \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                    	system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
						//system(("grep -Fw \""+A.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    	ifstream inf("./tmp1");
                    	inf >> ws;

                    	if (inf.eof()) {
                        	count++;
							if (1) cout << "EXHAUSTIVE ADDITION: " << smi << " { "<< j << " " << A.data->a[j].name << " on " << ACi <<" } -> " << A.molesmi << " {No. " << count << " }" << endl;

                        	out << setw(7) << left << count << " "
                            	<< setw(65) << left << A.molesmi << " "
                            	<< setw(17) << left << ACi << " "
                            	<< setw(17) << left << j << " "
                            	<< setw(30) << left << A.data->a[j].name << endl;

							if (para.ifwritemds) {
                				ofstream outs((para.smidir+"mds/"+mm+"_add.enc").c_str());
								//ofstream outs((para.outdir+mm+"_add.enc").c_str());
                				A.print(outs);
                				outs.close();
							}

                			ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
							//ofstream outs((para.outdir+"mds/DATLIST.txt").c_str(),ios::app);
                			outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_add.enc") << endl;
							//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_add.enc") << endl;
                			outs.close();

							if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

							if (para.enumeration) SMI_Enumerator(A.molesmi);
                    	}
                    	inf.close();
					
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

	para.stat << "NEU | ADDITION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}

void exhaustive_insertion(MOLECULE &A) {
    int count=0,i=0,j=0,k=0,m=0;
    string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);

    ofstream out("insertion.txt",ios::app);
    //mark(out);

    out << "MOL: " << smi << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_ins" << " "
        << setw(17) << left << "MOL_ins_pt" << " "
        << setw(17) << left << "Ins_id" << " "
        << setw(30) << left << "Ins_name" << " "
        << setw(17) << left << "Ins_bnd2par" << " "
        << setw(17) << left << "Ins_bnd2des" << endl;

    for (i=1;i<A.Cindex.size();i++) {
		int ACi=A.Cindex.at(i);
        for (j=1;j<A.data->num;j++) {
            if (A.data->a[j].probability>0) {
				for (m=1;m<=3;m++) {
                	for (k=1;k<=3;k++) {
                        if (A.insertion(i,j,m,k)) {
                            //count++;
                            //A.mds2smi();
                            A.canonicalize_SMILES();
                            A.smiles=A.molesmi;
                            //A.input();

                            //out<<A.molesmi<<endl;

                            string mm=A.molesmi;
                            for (int q=0;q<mm.length();q++) {
                                if (mm[q]=='/') mm[q]='u';
                                if (mm[q]=='\\') mm[q]='d';
                                if (mm[q]=='*') mm[q]='x';
                            }

                            system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
							//system(("grep -Fw \""+A.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                            ifstream inf("./tmp1");
                            inf >> ws;

                            if (inf.eof()) {
                                count++;

                            	if (1) cout << "EXHAUSTIVE INSERTION: " << smi << " {id: " << j << " " << A.data->a[j].name << " | on " << ACi << " | bnd2par " << m << " | bnd2des " << k << " } -> " << A.molesmi << " {No. " << count << " }" << endl;

                                out << setw(7) << left << count << " "
                                    << setw(65) << left << A.molesmi << " "
                                    << setw(17) << left << ACi << " "
                                    << setw(17) << left << A.data->a[j].id << " "
                                    << setw(30) << left << A.data->a[j].name << " "
                                    << setw(17) << left << k << " "
                                    << setw(17) << left << m << endl;

                            	string mm=A.molesmi;
                            	for (int q=0;q<mm.length();q++) {
                                	if (mm[q]=='/') mm[q]='u';
                                	if (mm[q]=='\\') mm[q]='d';
									if (mm[q]=='*') mm[q]='x';
                            	}

								if (para.ifwritemds) {
                            		ofstream outs((para.smidir+"mds/"+mm+"_ins.enc").c_str());
									//ofstream outs((para.outdir+"mds/"+mm+"_ins.enc").c_str());
                            		A.print(outs);
                            		outs.close();
								}

                            	ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
								//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
                            	outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_ins.enc") << endl;
								//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_ins.enc") << endl;
                            	outs.close();

								if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

								if (para.enumeration) SMI_Enumerator(A.molesmi);
                            }
                            inf.close();


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

    out.close();

    delete [] buf;
    buf=NULL;

	para.stat << "NEU | INSERTION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}

void exhaustive_change_bnd(MOLECULE &A) {
    int count=0,i=0,j=0,k=0,m=0;
    string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);

    ofstream out("change_bnd.txt",ios::app);
    //mark(out);

    out << "MOL: " << smi << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_chbnd" << " "
        << setw(17) << left << "MOL_chbnd_oatm_P" << " "
        << setw(17) << left << "MOL_chbnd_oatm_C" << " "
        << setw(17) << left << "MOL_chbnd_obnd" << " "
        << setw(17) << left << "MOL_chbnd_nid_P" << " "
        << setw(30) << left << "MOL_chbnd_nPname" << " "
        << setw(17) << left << "MOL_chbnd_nid_C" << " "
        << setw(30) << left << "MOL_chbnd_nCname" << " "
        << setw(17) << left << "MOL_chbnd_nbnd" << endl;


    for (i=1;i<A.Cindex.size();i++) {
		int APi=A.Pindex.at(i);
        int ACi=A.Cindex.at(i);
		int obnd=A.Rindex.at(i);
        for (j=1;j<A.data->num;j++) {
            if (A.data->a[j].probability>0) {
                for (k=1;k<A.data->num;k++) {
					if (A.data->a[k].probability>0) {
                    	for (m=1;m<=3;m++) {
                            if (A.change_bnd(i,j,k,m)) {
                                //count++;
                                //A.mds2smi();
                                A.canonicalize_SMILES();
                                A.smiles=A.molesmi;
                                //A.input();

                                //out<<A.molesmi<<endl;

                                string mm=A.molesmi;
                                for (int q=0;q<mm.length();q++) {
                                    if (mm[q]=='/') mm[q]='u';
                                    if (mm[q]=='\\') mm[q]='d';
                                    if (mm[q]=='*') mm[q]='x';
                                }

                                system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
								//system(("grep -Fw \""+A.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                                ifstream inf("./tmp1");
                                inf >> ws;

                                if (inf.eof()) {
                                    count++;

                                    if (1) cout << "EXHAUSTIVE CHANGE_BND: " << smi << " {id1: " << j << " " << A.data->a[j].name << " for " << APi << " | id2: " << k << " " << A.data->a[k].name << " for " << ACi << " | bnd "<< m <<" } -> " << A.molesmi << " {No. " << count << " }" << endl;


                                    out << setw(7) << left << count << " "
                                        << setw(65) << left << A.molesmi << " "
										<< setw(17) << left << APi << " "
                                        << setw(17) << left << ACi << " "
										<< setw(17) << left << obnd << " "
                                        << setw(17) << left << A.data->a[j].id << " "
                                        << setw(30) << left << A.data->a[j].name << " "
                                        << setw(17) << left << A.data->a[k].id << " "
                                        << setw(30) << left << A.data->a[k].name << " "
                                        << setw(17) << left << m << endl;

									if (para.ifwritemds) {
                                    	ofstream outs((para.smidir+"mds/"+mm+"_chbnd.enc").c_str());
										//ofstream outs((para.outdir+mm+"_chbnd.enc").c_str());
                                    	A.print(outs);
                                    	outs.close();
									}

                                    ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
									//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
                                    outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_chbnd.enc") << endl;
									//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_chbnd.enc") << endl;
                                    outs.close();

                                    if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

                                    if (para.enumeration) SMI_Enumerator(A.molesmi);
                                }
                                inf.close();


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

	para.stat << "NEU | CHANGE_BND | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}


void exhaustive_change_ele(MOLECULE &A) {
	int count=0,i=0,j=0,k=0,m=0;
	string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);
	
	ofstream out("change_ele.txt",ios::app);
	//mark(out);

    out << "MOL: " << smi << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_chele" << " "
        << setw(17) << left << "MOL_chele_oele_C" << " "
        << setw(17) << left << "MOL_chele_id" << " "
        << setw(30) << left << "MOL_chele_name" << " "
        << setw(17) << left << "MOL_chele_bnd2par" << " "
        << setw(17) << left << "MOL_chele_bnd2des" << endl;

	for (i=1;i<A.Cindex.size();i++) {
		int ACi=A.Cindex.at(i);
		for (j=1;j<A.data->num;j++) {
			if (A.data->a[j].probability>0) {
            	for (k=1;k<=3;k++) {
					for (m=1;m<=3;m++) {
						if (A.change_ele(i,j,k,m)) {
							//count++;
							//A.mds2smi();
							A.canonicalize_SMILES();
							A.smiles=A.molesmi;  
							//A.print();

							//A.input();
						
							//out<<A.molesmi<<endl;

                            string mm=A.molesmi;
                            for (int q=0;q<mm.length();q++) {
                                if (mm[q]=='/') mm[q]='u';
                                if (mm[q]=='\\') mm[q]='d';
                                if (mm[q]=='*') mm[q]='x';
                            }

                            system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
							//system(("grep -Fw \""+A.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                            ifstream inf("./tmp1");
                            inf >> ws;

                            if (inf.eof()) {
                                count++;

								if (1) cout << "EXHAUSTIVE CHANGE_ELE: " << smi << " {id1: " << j << " " << A.data->a[j].name << " on " << ACi << " | bnd2par " << k << " bnd2des " << m << " } -> " << A.molesmi << " {No. " << count << " }" << endl;

                                out << setw(7) << left << count << " "
                                    << setw(65) << left << A.molesmi << " "
                                    << setw(17) << left << ACi << " "
                                    << setw(17) << left << A.data->a[j].id << " "
                                    << setw(30) << left << A.data->a[j].name << " "
                                    << setw(17) << left << k << " "
                                    << setw(17) << left << m << endl;

								if (para.ifwritemds) {
                					ofstream outs((para.smidir+"mds/"+mm+"_chele.enc").c_str());
									//ofstream outs((para.outdir+mm+"_chele.enc").c_str());
                					A.print(outs);
                					outs.close();
								}

                				ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
								//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
                				outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_chele.enc") << endl;
								//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_chele.enc") << endl;
                				outs.close();

								if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

								if (para.enumeration) SMI_Enumerator(A.molesmi);
                            }
                            inf.close();

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
	
	out.close();

	delete [] buf;
	buf=NULL;

	para.stat << "NEU | CHANGE_ELE | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}

void exhaustive_cyclization(MOLECULE &A) {
	int count=0,i=0,j=0;
	string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);
	
	ofstream out("cyclization.txt",ios::app);
	//mark(out);

    out << "MOL: " << smi  << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_cyc" << " "
        << setw(17) << left << "MOL_cyc_pt1" << " "
        << setw(17) << left << "MOL_cyc_pt2" << " "
		<< setw(17) << left << "Bnd" << endl;


	for (i=0;i<A.Cindex.size();i++) {
		for (j=i+1;j<A.Cindex.size();j++) {
            int ACi=A.Cindex.at(i);
            int ACj=A.Cindex.at(j);
			for (int k=1;k<=3;k++) {
				if (A.cyclization(i,j,k)) {
					//count++;
					//A.mds2smi();
					A.canonicalize_SMILES();
					A.smiles=A.molesmi;  
					//A.input();
				
					//out<<A.molesmi<<endl;

                    string mm=A.molesmi;
                    for (int q=0;q<mm.length();q++) {
                        if (mm[q]=='/') mm[q]='u';
                        if (mm[q]=='\\') mm[q]='d';
                        if (mm[q]=='*') mm[q]='x';
                    }

                    //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                    system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
					//system(("grep -Fw \""+A.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                    ifstream inf("./tmp1");
                    inf >> ws;

                    if (inf.eof()) {
                        count++;

						if (1) cout << "EXHAUSTIVE CYCLIZATION: " << smi << " { " << ACi << " & " << ACj <<" | bnd " << k << " } -> " << A.molesmi << " {No. " << count << " }" << endl;

                        out << setw(7) << left << count << " "
                            << setw(65) << left << A.molesmi << " "
                            << setw(17) << left << ACi << " "
                            << setw(17) << left << ACj << " "
							<< setw(17) << left << k << endl;

						if (para.ifwritemds) {
                			ofstream outs((para.smidir+"mds/"+mm+"_cyc.enc").c_str());
							//ofstream outs((para.outdir+mm+"_cyc.enc").c_str());
                			A.print(outs);
                			outs.close();
						}

                		ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
						//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
                		outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_cyc.enc") << endl;
						//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_cyc.enc") << endl;
                		outs.close();

						if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

						if (para.enumeration) SMI_Enumerator(A.molesmi);
                    }
                    inf.close();


					if (0) {
						A.smiles=smi;
						A.input();
					}
					if (1) A.replace(buf[0]);
				}
			}
		}
	}
	
	out.close();

	delete [] buf;
	buf=NULL;

	para.stat << "NEU | CYCLIZATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}

void exhaustive_decyclization(MOLECULE &A) {
    int count=0,i=0,j=0;
    string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);

    ofstream out("decyclization.txt",ios::app);
    //mark(out);

    out << "MOL: " << smi << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_decyc" << " "
        << setw(17) << left << "MOL_cyc_num" << " "
        << setw(17) << left << "MOL_bnd" << endl;

    for (i=1;i<=A.if_circle;i++) {
        int cycnum=i;
		int cybnd=A.Cybnd.at(cycnum-1);
        if (A.decyclization(i)) {
            //count++;
            //A.mds2smi();
            A.canonicalize_SMILES();
            A.smiles=A.molesmi;
            //A.input();

            //out<<A.molesmi<<endl;

            string mm=A.molesmi;
            for (int q=0;q<mm.length();q++) {
                if (mm[q]=='/') mm[q]='u';
                if (mm[q]=='\\') mm[q]='d';
                if (mm[q]=='*') mm[q]='x';
            }

            //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
            system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
			//system(("grep -Fw \""+A.molesmi+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
            ifstream inf("./tmp1");
            inf >> ws;

            if (inf.eof()) {
                count++;

                if (1) cout << "EXHAUSTIVE DECYCLIZATION: " << smi << " { ring no. "<< cycnum << " | bnd " << cybnd << " } -> " << A.molesmi << " {No. " << count << " }" << endl;

                out << setw(7) << left << count << " "
                    << setw(65) << left << A.molesmi << " "
                    << setw(17) << left << cycnum << " "
                    << setw(17) << left << cybnd << endl;

				if (para.ifwritemds) {
                	ofstream outs((para.smidir+"mds/"+mm+"_decyc.enc").c_str());
					//ofstream outs((para.outdir+mm+"_decyc.enc").c_str());
                	A.print(outs);
                	outs.close();
				}

                ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
				//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
                outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_decyc.enc") << endl;
				//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_decyc.enc") << endl;
                outs.close();

				if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

				if (para.enumeration) SMI_Enumerator(A.molesmi);
            }
            inf.close();


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

	para.stat << "NEU | DECYCLIZATION | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}

void exhaustive_change_cistrans(MOLECULE &A) {
    int count=0,i=0,j=0;
    string smi=A.molesmi;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);

    ofstream out("change_cistrans.txt",ios::app);
    //mark(out);

    out << "MOL: " << smi  << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_ct" << " "
        << setw(17) << left << "MOL_ct_pt" << " "
        << setw(17) << left << "MOL_ct_h_ori" << " "
        << setw(17) << left << "MOL_ct_h_after" << " "
        << setw(17) << left << "MOL_ct_e_ori" << " "
        << setw(17) << left << "MOL_ct_e_after" << endl;

    for (int w=0;w<2;w++) {
        for (i=0;i<A.ctsisomer.at(w).size();i++) {
            int ACi=A.Cindex.at(i);
            string ACTh=A.ctsisomer.at(0).at(i);
            string ACTe=A.ctsisomer.at(1).at(i);

            if (A.change_cistrans(i,w)) {
                //A.mds2smi();
                if (0) {
                    A.smiles=A.molesmi;
                    A.input();
                }
                if (1) A.canonicalize_SMILES();

                //out<<A.molesmi<<endl;

                string mm=A.molesmi;
                for (int q=0;q<mm.length();q++) {
                    if (mm[q]=='/') mm[q]='u';
                    if (mm[q]=='\\') mm[q]='d';
					if (mm[q]=='*') mm[q]='x';
                }

                //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
				//system(("grep -Fw \""+mm+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                ifstream inf("./tmp1");
                inf >> ws;
                if (inf.eof()) {
                    count++;
					if (1) {
                    	if (w==0) cout << "EXHAUSTIVE CHANGE_CISTRANS: " << smi << " { "<< ACi << " | front | from " << ACTh << " to " << A.ctsisomer.at(0).at(i) << " } -> " << A.molesmi << " {No. " << count << " }" << endl;
                    	else if (w==1) cout << "EXHAUSTIVE CHANGE_CISTRANS: " << smi << " { "<< ACi << " | end | from " << ACTe << " to " << A.ctsisomer.at(1).at(i) << " } -> " << A.molesmi << " {No. " << count << " }" << endl;
					}

                    out << setw(7) << left << count << " "
                        << setw(65) << left << A.molesmi << " "
                        << setw(17) << left << ACi << " "
                        << setw(17) << left << ACTh << " "
                        << setw(17) << left << A.ctsisomer.at(0).at(i) << " "
                        << setw(17) << left << ACTe << " "
                        << setw(17) << left << A.ctsisomer.at(1).at(i) << endl;

                    //string mm=A.ion[s].molesmi;
                    //for (int q=0;q<mm.length();q++) {
                    //    if (mm[q]=='/') mm[q]='u';
                    //    if (mm[q]=='\\') mm[q]='d';
                    //}

					if (para.ifwritemds) {
                    	ofstream outs((para.smidir+"mds/"+mm+"_chct.enc").c_str());
						//ofstream outs((para.outdir+mm+"_chct.enc").c_str());
                    	A.print(outs);
                    	outs.close();
					}

                    ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
					//ofstream outs((para.outdir+"/DATLIST.txt").c_str(),ios::app);
                    outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_chct.enc") << endl;
					//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_chct.enc") << endl;
                    outs.close();

					if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

                    if (para.enumeration) SMI_Enumerator(A.molesmi);
                }
                inf.close();

                if (0) {
                    A.smiles=smi;
                    A.input();
                }
                if (1) {
                    A.replace(buf[0]);
                }
            }

        }
    }


    out.close();

    delete [] buf;
    buf=NULL;

	para.stat << "NEU | CHANGE_CISTRANS | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}


void exhaustive_change_chirality(MOLECULE &A) {
    int count=0,i=0,j=0;
    string smi=A.molesmi;

    vector<int> stereoC(2);
    stereoC.at(0)=67;
    stereoC.at(1)=68;

    MOLECULE *buf=new MOLECULE [1];
    buf[0].replace(A);

    ofstream out("change_chirality.txt",ios::app);
    //mark(out);

    out << "MOL: " << smi << endl;
    A.print(out);
    out << endl;

    out << setw(7) << left << "No." << " "
        << setw(65) << left << "MOL_chir" << " "
        << setw(17) << left << "MOL_chir_pt" << " "
        << setw(17) << left << "MOL_chir_ori" << " "
        << setw(17) << left << "MOL_chir_after" << endl;

    for (i=0;i<A.Cindex.size();i++) {
        int ACi=A.Cindex.at(i);
        int AMi=A.Mindex.at(i);

        for (j=0;j<stereoC.size();j++) {
            if (A.change_chirality(i,stereoC.at(j))) {
                //A.mds2smi();
                if (0) {
                    A.smiles=A.molesmi;
                    A.input();
                }
                if (1) A.canonicalize_SMILES();

                //out<<A.molesmi<<endl;

                string mm=A.molesmi;
                for (int q=0;q<mm.length();q++) {
                    if (mm[q]=='/') mm[q]='u';
                    if (mm[q]=='\\') mm[q]='d';
					if (mm[q]=='*') mm[q]='x';
                }

                //system(("ls "+para.smidir+"mds/ | grep -Fw \""+A.molesmi+"_\" > ./tmp1 2> /dev/null").c_str());
                system(("grep -F \"/"+mm+"_\" "+para.smidir+"mds/DATLIST.txt > ./tmp1 2> /dev/null").c_str());
				//system(("grep -Fw \""+mm+"\" "+para.outdir+"DATLIST.txt > ./tmp1 2> /dev/null").c_str());
                ifstream inf("./tmp1");
                inf >> ws;

                if (inf.eof()) {
                    count++;
                    if (0) cout << "EXHAUSTIVE CHANGE_CHIRALITY: " << smi << " {CAT: "<< ACi << " | " << AMi << " to " << A.Mindex.at(i) << " } -> " << A.molesmi << " {No. " << count << " }" << endl;

                    out << setw(7) << left << count << " "
                        << setw(65) << left << A.molesmi << " "
                        << setw(17) << left << ACi << " "
                        << setw(17) << left << AMi << " "
                        << setw(17) << left << A.Mindex.at(i) << endl;

                    //string mm=A.molesmi;
                    //for (int q=0;q<mm.length();q++) {
                    //  if (mm[q]=='/') mm[q]='u';
                    //  if (mm[q]=='\\') mm[q]='d';
                    //}

					if (para.ifwritemds) {
                    	ofstream outs((para.smidir+"mds/"+mm+"_chchir.enc").c_str());
						//ofstream outs((para.outdir+mm+"_chchir.enc").c_str());
                		A.print(outs);
                    	outs.close();
					}

                    ofstream outs((para.smidir+"mds/DATLIST.txt").c_str(),ios::app);
					//ofstream outs((para.outdir+"DATLIST.txt").c_str(),ios::app);
                    outs << setw(50) << left << A.molesmi << "   " << left << (para.smidir+"mds/"+mm+"_chchir.enc") << endl;
					//outs << setw(50) << left << A.molesmi << "   " << left << (para.outdir+mm+"_chchir.enc") << endl;
                    outs.close();

					if (para.glbouf.is_open()) para.glbouf << A.molesmi << endl;

                    if (para.enumeration) SMI_Enumerator(A.molesmi);
                }
                inf.close();

                if (0) {
                    A.smiles=smi;
                    A.input();
                }
                if (1) {
                    A.replace(buf[0]);
                    //A.reset();
                }
            }
        }
    }

    out.close();

    delete [] buf;
    buf=NULL;

	para.stat << "NEU | CHANGE_CHIRALITY | ROUND: " << setw(4) << para.round << " " << " | COUNTS: " << count << endl;
}



void CH4_to_Bz_path1(MOLECULE &A) { //MOLECULE &A

	ofstream out("benzene_path1.txt");
	//mark(out);	
	
	//cout<<"Make benzene, c1ccccc1 from methane, C"<<endl; //C
	//A.smiles="C";
	//A.input();
	A.empty();
	A.clear();
	A.Pindex.push_back(0);
	A.Cindex.push_back(1);
	A.Mindex.push_back(1);
	A.Rindex.push_back(0);
	A.Cyindex.push_back(vector<int>(0));
	if (para.protect) A.protect.push_back(1);
	A.ctsisomer.resize(2,vector<string>(0));
	A.ctsisomer.at(0).push_back("");
	A.ctsisomer.at(1).push_back("");
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<"Initialization "<<A.molesmi<<endl;

	A.addition(0,2,1); //CC
	//A.mds2smi();
	A.canonicalize_SMILES();
	A.smiles=A.molesmi;  
	//A.input();
	out<<"addition(0,2,1) "<<A.molesmi<<endl;

    A.change_ele(0,2,1,2); 
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out<<"change_ele(0,2,1,2) "<<A.molesmi<<endl;	

	/*
    A.change_bnd(1,2,2,2); //C=C
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out<<"change_bnd(1,2,2,2) "<<A.molesmi<<endl;
	*/

	A.addition(1,2,1); //C=CC
	//A.mds2smi();
	A.canonicalize_SMILES();
	A.smiles=A.molesmi;  
	//A.input();
	out<<"addition(1,2,1) "<<A.molesmi<<endl;
	
	A.addition(2,2,2); //C=CC=C
	//A.mds2smi();
	A.canonicalize_SMILES();
	A.smiles=A.molesmi;  
	//A.input();
	out<<"addition(2,2,2) "<<A.molesmi<<endl;
	
	A.addition(3,2,1); //C=CC=CC
	//A.mds2smi();
	A.canonicalize_SMILES();
	A.smiles=A.molesmi;  
	//A.input();
	out<<"addition(3,2,1) "<<A.molesmi<<endl;
	
	A.addition(4,2,2); //C=CC=CC=C
	//A.mds2smi();
	A.canonicalize_SMILES();
	A.smiles=A.molesmi;  
	//A.input();
	out<<"addition(4,2,2) "<<A.molesmi<<endl;
	
	A.cyclization(0,5,1); //C1=CC=CC=C1
	//A.mds2smi();
	A.canonicalize_SMILES();
	A.smiles=A.molesmi;  
	//A.input();
	out<<"cyclization(0,5,1) "<<A.molesmi<<endl;

	out.close();
}

void CH4_to_Bz_path2(MOLECULE &A) { //MOLECULE &A

    ofstream out("benzene_path2.txt");
    //mark(out);

    //cout<<"Make benzene, c1ccccc1 from methane, C"<<endl; //C
    //A.smiles="C";
    //A.input();
	A.empty();
	A.clear();
    A.Pindex.push_back(0);
    A.Cindex.push_back(1);
    A.Mindex.push_back(1);
    A.Rindex.push_back(0);
    A.Cyindex.push_back(vector<int>(0));
    if (para.protect) A.protect.push_back(1);
	A.ctsisomer.resize(2,vector<string>(0));
    A.ctsisomer.at(0).push_back("");
    A.ctsisomer.at(1).push_back("");
    A.mds2smi();
	A.canonicalize_SMILES();
    out<<"Initialization "<<A.molesmi<<endl;

    A.addition(0,2,1); //CC
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
	out<<"addition(0,2,1) "<<A.molesmi<<endl;

	/*
    A.change_ele(0,2,1,2); //C=C
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out<<"change_ele(0,2,1,2) "<<A.molesmi<<endl;
	*/
    
    A.change_bnd(1,2,2,2); //C=C
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out<<"change_bnd(1,2,2,2) "<<A.molesmi<<endl;

    A.addition(1,2,1); //C=CC
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
	out<<"addition(1,2,1) "<<A.molesmi<<endl;

    A.addition(2,2,2); //C=CC=C
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out<<"addition(2,2,2) "<<A.molesmi<<endl;

    A.addition(3,2,1); //C=CC=CC
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out<<"addition(3,2,1) "<<A.molesmi<<endl;

    A.addition(4,2,2); //C=CC=CC=C
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out<<"addition(4,2,2) "<<A.molesmi<<endl;

    A.cyclization(0,5,1); //C1=CC=CC=C1
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out<<"cyclization(0,5,1) "<<A.molesmi<<endl;

    out.close();
}


void CH4_to_aspirin(MOLECULE &A) {
	ofstream out("aspirin.txt");
	//mark(out);

    //cout << "Creating aspirin from methane, C" << endl;

    //A.smiles = "C";
    //A.input();
	A.empty();
	A.clear();
    A.Pindex.push_back(0);
    A.Cindex.push_back(1);
    A.Mindex.push_back(1);
    A.Rindex.push_back(0);
    A.Cyindex.push_back(vector<int>(0));
    if (para.protect) A.protect.push_back(1);
	A.ctsisomer.resize(2,vector<string>(0));
    A.ctsisomer.at(0).push_back("");
    A.ctsisomer.at(1).push_back("");
	A.mds2smi();
	A.canonicalize_SMILES();
	out<<"Initialization "<<A.molesmi<<endl;

    A.addition(0,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
	//A.input();
    out << "addition(0,2,1) " << A.molesmi << endl;

    //A.addition(1,6,2);
    ////A.mds2smi();
    //A.canonicalize_SMILES();
    //A.smiles=A.molesmi;
    ////A.input();
    //out << "addition(1,6,2) " << A.molesmi << endl;

    A.addition(1,8,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(1,8,2) " << A.molesmi << endl;

    A.change_ele(2,6,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "change_ele(2,6,2,1) " << A.molesmi << endl;

    //A.addition(1,5,1);
    ////A.mds2smi();
    //A.canonicalize_SMILES();
    //A.smiles=A.molesmi;
    ////A.input();
    //out << "addition(1,5,1) " << A.molesmi << endl;

    //A.addition(3,2,1);
    ////A.mds2smi();
    //A.canonicalize_SMILES();
    //A.smiles=A.molesmi;
    ////A.input();
    //out << "addition(3,2,1) " << A.molesmi << endl;

    A.addition(1,2,1);
    A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(1,2,1) " << A.molesmi << endl;

    A.insertion(3,5,1,1);
    A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "insertion(3,5,1,1) " << A.molesmi << endl;


    A.addition(4,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(4,2,2) " << A.molesmi << endl;

    A.addition(5,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(5,2,1) " << A.molesmi << endl;

    //A.addition(6,6,2);
    ////A.mds2smi();
    //A.canonicalize_SMILES();
    //A.smiles=A.molesmi;
    ////A.input();
    //out << "addition(6,6,2) " << A.molesmi << endl;

    A.addition(6,11,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(6,11,1) " << A.molesmi << endl;

    A.change_bnd(7,2,6,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "change_bnd(7,2,6,2) " << A.molesmi << endl;


    A.addition(6,5,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(6,5,1) " << A.molesmi << endl;

    A.addition(5,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(5,2,1) " << A.molesmi << endl;

    A.addition(9,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(9,2,2) " << A.molesmi << endl;

    A.addition(10,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(10,2,1) " << A.molesmi << endl;

    A.addition(11,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(11,2,2) " << A.molesmi << endl;

    A.cyclization(4,12,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "cyclization(4,12,1) " << A.molesmi << endl;

	out.close();

}

void CH4_to_alpha_carotene_path3(MOLECULE &A) {
    ofstream out("carotene_path3.txt");
    //mark(out);

    //cout << "Creating aspirin from methane, C" << endl;

    //A.smiles = "C";
    //A.input();
	A.empty();
	A.clear();
    A.Pindex.push_back(0);
    A.Cindex.push_back(1);
    A.Mindex.push_back(1);
    A.Rindex.push_back(0);
    A.Cyindex.push_back(vector<int>(0));
    if (para.protect) A.protect.push_back(1);
	A.ctsisomer.resize(2,vector<string>(0));
    A.ctsisomer.at(0).push_back("");
    A.ctsisomer.at(1).push_back("");
	A.mds2smi();
	A.canonicalize_SMILES();
    out<<"Initialization "<<A.molesmi<<endl;

    A.addition(0,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(0,1,1) " << A.molesmi << endl;

    A.addition(0,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(0,1,1) " << A.molesmi << endl;

    A.addition(0,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(0,1,1) " << A.molesmi << endl;

    A.addition(3,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(3,1,1) " << A.molesmi << endl;

    A.addition(4,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(4,1,1) " << A.molesmi << endl;

    A.addition(5,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(5,2,1) " << A.molesmi << endl;

    A.addition(6,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
	out << "addition(6,2,2) " << A.molesmi << endl;

    A.addition(6,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(6,1,1) " << A.molesmi << endl;

    A.addition(7,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(7,2,1) " << A.molesmi << endl;

    A.cyclization(0,7,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "cyclization(0,7,1) " << A.molesmi << endl;

    A.addition(9,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(9,2,2) " << A.molesmi << endl;

    A.addition(10,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(10,2,1) " << A.molesmi << endl;

    A.addition(11,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(11,2,2) " << A.molesmi << endl;

    A.addition(12,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(12,2,1) " << A.molesmi << endl;

    A.addition(11,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(11,1,1) " << A.molesmi << endl;

    A.addition(13,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(13,2,2) " << A.molesmi << endl;

    A.addition(15,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(15,2,1) " << A.molesmi << endl;

    A.addition(16,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(16,2,2) " << A.molesmi << endl;

    A.addition(17,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(17,2,1) " << A.molesmi << endl;

    A.addition(16,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(16,1,1) " << A.molesmi << endl;

    A.addition(18,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(18,2,2) " << A.molesmi << endl;

    A.addition(20,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(20,2,1) " << A.molesmi << endl;

    A.addition(21,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(21,2,2) " << A.molesmi << endl;

    A.addition(22,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(22,2,1) " << A.molesmi << endl;

    A.addition(22,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(22,1,1) " << A.molesmi << endl;

    A.addition(23,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(23,2,2) " << A.molesmi << endl;

    A.addition(25,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(25,2,1) " << A.molesmi << endl;

    A.addition(26,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(26,2,2) " << A.molesmi << endl;

    A.addition(27,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(27,2,1) " << A.molesmi << endl;

    A.addition(27,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(27,1,1) " << A.molesmi << endl;

    A.addition(28,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(28,2,2) " << A.molesmi << endl;

    A.addition(30,1,1);
    A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(30,1,1) " << A.molesmi << endl;

    A.addition(31,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(31,2,1) " << A.molesmi << endl;

    A.addition(32,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(32,1,1) " << A.molesmi << endl;

    A.addition(32,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(32,2,2) " << A.molesmi << endl;

    A.addition(34,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(34,1,1) " << A.molesmi << endl;

    //A.change_cistrans(32,0);
    ////A.mds2smi();
    //A.canonicalize_SMILES();
    //A.smiles=A.molesmi;
    ////A.input();
    //out << "change_cistrans(32,0) " << A.molesmi << endl;

    A.addition(35,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(35,1,1) " << A.molesmi << endl;

    A.addition(36,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(36,1,1) " << A.molesmi << endl;

    A.addition(37,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(37,1,1) " << A.molesmi << endl;

    A.addition(37,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(37,1,1) " << A.molesmi << endl;

    A.cyclization(31,37,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "cyclization(31,37,1) " << A.molesmi << endl;

}

void CH4_to_alpha_carotene_path1(MOLECULE &A) {
    ofstream out("carotene_path1.txt");
    //mark(out);

    //cout << "Creating aspirin from methane, C" << endl;

    //A.smiles = "C";
    //A.input();
	A.empty();
	A.clear();
    A.Pindex.push_back(0);
    A.Cindex.push_back(1);
    A.Mindex.push_back(1);
    A.Rindex.push_back(0);
    A.Cyindex.push_back(vector<int>(0));
    if (para.protect) A.protect.push_back(1);
	A.ctsisomer.resize(2,vector<string>(0));
    A.ctsisomer.at(0).push_back("");
    A.ctsisomer.at(1).push_back("");
	A.mds2smi();
	A.canonicalize_SMILES();
    out<<"Initialization "<<A.molesmi<<endl;

    A.addition(0,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(0,1,1) " << A.molesmi << endl;

    A.addition(0,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(0,1,1) " << A.molesmi << endl;

    A.addition(0,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(0,1,1) " << A.molesmi << endl;

    A.addition(3,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(3,1,1) " << A.molesmi << endl;

    A.addition(4,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(4,1,1) " << A.molesmi << endl;

    A.addition(5,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(5,2,1) " << A.molesmi << endl;

    A.addition(6,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(6,1,1) " << A.molesmi << endl;

    A.addition(6,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
	out << "addition(6,2,2) " << A.molesmi << endl;

    A.cyclization(0,8,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "cyclization(0,8,1) " << A.molesmi << endl;

    A.addition(8,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(8,2,1) " << A.molesmi << endl;

    A.addition(9,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(9,2,2) " << A.molesmi << endl;

    A.addition(10,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(10,2,1) " << A.molesmi << endl;

    A.addition(11,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(11,1,1) " << A.molesmi << endl;

    A.addition(11,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(11,2,2) " << A.molesmi << endl;

    A.addition(13,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(13,2,1) " << A.molesmi << endl;

    A.addition(14,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(14,2,2) " << A.molesmi << endl;

    A.addition(15,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(15,2,1) " << A.molesmi << endl;

    A.addition(16,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(16,1,1) " << A.molesmi << endl;

    A.addition(16,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(16,2,2) " << A.molesmi << endl;

    A.addition(18,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(18,2,1) " << A.molesmi << endl;

    A.addition(19,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(19,2,2) " << A.molesmi << endl;

    A.addition(20,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(20,2,1) " << A.molesmi << endl;

    A.addition(21,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(21,2,2) " << A.molesmi << endl;

    A.addition(22,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(22,1,1) " << A.molesmi << endl;

    A.addition(22,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(22,2,1) " << A.molesmi << endl;

    A.addition(24,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(24,2,2) " << A.molesmi << endl;

    A.addition(25,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(25,2,1) " << A.molesmi << endl;

    A.addition(26,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(26,2,2) " << A.molesmi << endl;

    A.addition(27,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(27,1,1) " << A.molesmi << endl;

    A.addition(27,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(27,2,1) " << A.molesmi << endl;

    A.addition(29,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(29,2,2) " << A.molesmi << endl;

    A.addition(30,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(30,1,1) " << A.molesmi << endl;

    A.addition(31,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(31,2,1) " << A.molesmi << endl;

    A.addition(32,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(32,1,1) " << A.molesmi << endl;

    A.addition(32,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(32,2,2) " << A.molesmi << endl;

    A.addition(34,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(34,1,1) " << A.molesmi << endl;

    //A.change_cistrans(32,0);
    ////A.mds2smi();
    //A.canonicalize_SMILES();
    //A.smiles=A.molesmi;
    ////A.input();
    //out << "change_cistrans(32,0) " << A.molesmi << endl;

    A.addition(35,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(35,1,1) " << A.molesmi << endl;

    A.addition(36,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(36,1,1) " << A.molesmi << endl;

    A.addition(37,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(37,1,1) " << A.molesmi << endl;

    A.addition(37,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(37,1,1) " << A.molesmi << endl;

    A.cyclization(31,37,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "cyclization(31,37,1) " << A.molesmi << endl;

    A.change_cistrans(22,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "change_cistrans(22,1) " << A.molesmi << endl;


}


void CH4_to_alpha_carotene_path2(MOLECULE &A) {
    ofstream out("carotene_path2.txt");
    //mark(out);

    //cout << "Creating aspirin from methane, C" << endl;

    //A.smiles = "C";
    //A.input();
	A.empty();
	A.clear();
    A.Pindex.push_back(0);
    A.Cindex.push_back(1);
    A.Mindex.push_back(1);
    A.Rindex.push_back(0);
    A.Cyindex.push_back(vector<int>(0));
    if (para.protect) A.protect.push_back(1);
	A.ctsisomer.resize(2,vector<string>(0));
    A.ctsisomer.at(0).push_back("");
    A.ctsisomer.at(1).push_back("");
	A.mds2smi();
	A.canonicalize_SMILES();
    out<<"Initialization "<<A.molesmi<<endl;

    A.addition(0,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(0,2,1) " << A.molesmi << endl;

    A.addition(1,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(1,1,1) " << A.molesmi << endl;

    A.addition(2,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(2,1,1) " << A.molesmi << endl;

    A.addition(3,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(3,1,1) " << A.molesmi << endl;

    A.addition(4,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(4,1,1) " << A.molesmi << endl;

    A.addition(5,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(5,1,1) " << A.molesmi << endl;

    A.addition(5,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(5,1,1) " << A.molesmi << endl;

    A.addition(5,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(5,2,1) " << A.molesmi << endl;

    A.cyclization(1,8,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "cyclization(1,8,2) " << A.molesmi << endl;

    A.addition(8,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(8,2,1) " << A.molesmi << endl;

    A.addition(9,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(9,2,2) " << A.molesmi << endl;

    A.addition(10,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(10,2,1) " << A.molesmi << endl;

    A.addition(11,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(11,1,1) " << A.molesmi << endl;

    A.addition(11,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(11,2,2) " << A.molesmi << endl;

    A.addition(13,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(13,2,1) " << A.molesmi << endl;

    A.addition(14,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(14,2,2) " << A.molesmi << endl;

    A.addition(15,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(15,2,1) " << A.molesmi << endl;

    A.addition(16,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(16,1,1) " << A.molesmi << endl;

    A.addition(16,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(16,2,2) " << A.molesmi << endl;

    A.addition(18,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(18,2,1) " << A.molesmi << endl;

    A.addition(19,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(19,2,2) " << A.molesmi << endl;

    A.addition(20,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(20,2,1) " << A.molesmi << endl;

    A.addition(21,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(21,2,2) " << A.molesmi << endl;

    A.addition(22,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(22,1,1) " << A.molesmi << endl;

    A.addition(22,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(22,2,1) " << A.molesmi << endl;

    A.addition(24,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(24,2,2) " << A.molesmi << endl;

    A.addition(25,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(25,2,1) " << A.molesmi << endl;

    A.addition(26,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(26,2,2) " << A.molesmi << endl;

    A.addition(27,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(27,1,1) " << A.molesmi << endl;

    A.addition(27,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(27,2,1) " << A.molesmi << endl;

    A.addition(29,2,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(29,2,2) " << A.molesmi << endl;

    A.addition(30,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(30,1,1) " << A.molesmi << endl;

    A.addition(31,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(31,2,1) " << A.molesmi << endl;

    A.addition(32,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(32,1,1) " << A.molesmi << endl;

    A.addition(31,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(31,1,1) " << A.molesmi << endl;

    A.addition(34,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(34,1,1) " << A.molesmi << endl;

    A.addition(34,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(34,1,1) " << A.molesmi << endl;

    A.addition(34,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(34,1,1) " << A.molesmi << endl;

    A.addition(37,1,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(37,1,1) " << A.molesmi << endl;

    A.addition(38,2,1);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "addition(38,2,1) " << A.molesmi << endl;

    A.cyclization(32,39,2);
    //A.mds2smi();
    A.canonicalize_SMILES();
    A.smiles=A.molesmi;
    //A.input();
    out << "cyclization(32,39,2) " << A.molesmi << endl;

}
