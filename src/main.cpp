#include "ELEMENTS.h"
#include "MOLECULE.h"
//#include "MATRIX.h"
//#include "OPT.h"
#include "IL.h"
#include "PARAMETER.h"
#include "CASES_NEU.h"
#include "CASES_IL.h"
#include "CASES_IL_INDEPENDENT.h"
#include "UTILITY.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
using namespace std;


PARAMETER para;


int main(int argc, char **argv) {
	mark(cout);
	if (argc<=1) {
		cout << "Please specify at least a control file (.in file)" << endl;
		exit(1);
	}
	rd_para(argv);
	output_para();
	clearlog();

	MOLECULE *gs=NULL;
	IL *gsion=NULL;

	POOL data;
	if (!data.read_in()) data.set_up();

	mk_datlist(); 

	para.stat.open((para.logdir+"statistics.txt").c_str());

	if (!para.ion) {
		cout<<"Read the molecules."<<endl;    


		stringstream num(""),num1("");

		while (para.round<=para.epoch) {

			int nmol=0; //=para.gssize;
			if (gs!=NULL) {
				delete [] gs;
				gs=NULL;
			}
			gs=new MOLECULE [2];
			gs[0].data=&data;
			gs[1].data=&data;

			num.str("");
			num << para.round-1;

			if (para.round==1) nmol=cal_nmol(para.guess);
			else nmol=cal_nmol(para.logdir+"mol-"+num.str()+"gen.txt");

			num1.str("");
			num1 << para.round;
			para.glbouf.open((para.logdir+"mol-"+num1.str()+"gen.txt").c_str());
			for (int i=0;i<nmol;i++) {
				string smi="";//gs[0].molesmi;

				if (para.round==1) smi=rd_1molsmi(para.guess,i);
				else smi=rd_1molsmi((para.logdir+"mol-"+num.str()+"gen.txt"),i);

				if (para.operation=="Change_ele" || para.operation=="All") {
					if (para.round==1) {
						gs[0].molesmi=gs[0].smiles=smi;   
						gs[0].input();
					}
					else gs[0].readmds(smi);

					cout<<"Do exhaustive element change"<<endl;
					exhaustive_change_ele(gs[0]);
				}

				if (para.operation=="Change_bnd" || para.operation=="All") {
                    if (para.round==1) {
                        gs[0].molesmi=gs[0].smiles=smi;
                        gs[0].input();
                    }
                    else gs[0].readmds(smi);

					cout<<"Do exhaustive bond change"<<endl;
					exhaustive_change_bnd(gs[0]);
				}

				if (para.operation=="Addition" || para.operation=="All") {
                    if (para.round==1) {
                        gs[0].molesmi=gs[0].smiles=smi;
                        gs[0].input();
                    }
                    else gs[0].readmds(smi);
					
					cout<<"Do exhaustive addition"<<endl;
					exhaustive_addition(gs[0]);
				}

				if (para.operation=="Insertion" || para.operation=="All") {
                    if (para.round==1) {
                        gs[0].molesmi=gs[0].smiles=smi;
                        gs[0].input();
                    }
                    else gs[0].readmds(smi);

					cout<<"Do exhaustive insertion"<<endl;
					exhaustive_insertion(gs[0]);
				}

				if (para.operation=="Subtraction" || para.operation=="All") {
                    if (para.round==1) {
                        gs[0].molesmi=gs[0].smiles=smi;
                        gs[0].input();
                    }
                    else gs[0].readmds(smi);

					cout<<"Do exhaustive subtraction"<<endl;
					exhaustive_subtraction(gs[0]);
				}

				if (para.operation=="Cistrans_inversion" || para.operation=="All") {
                    if (para.round==1) {
                        gs[0].molesmi=gs[0].smiles=smi;
                        gs[0].input();
                    }
                    else gs[0].readmds(smi);

					cout<<"Do exhaustive cis-trans inversion"<<endl;
					exhaustive_change_cistrans(gs[0]);
				}

				if (para.operation=="Chirality_inversion" || para.operation=="All") {
                    if (para.round==1) {
                        gs[0].molesmi=gs[0].smiles=smi;
                        gs[0].input();
                    }
                    else gs[0].readmds(smi);

					cout<<"Do exhaustive chirality inversion"<<endl;
					exhaustive_change_chirality(gs[0]);
				}

				if (para.operation=="Cyclization" || para.operation=="All") {
                    if (para.round==1) {
                        gs[0].molesmi=gs[0].smiles=smi;
                        gs[0].input();
                    }
                    else gs[0].readmds(smi);

					cout<<"Do exhaustive cyclization"<<endl;
					exhaustive_cyclization(gs[0]);
				}

				if (para.operation=="Decyclization" || para.operation=="All") {
                    if (para.round==1) {
                        gs[0].molesmi=gs[0].smiles=smi;
                        gs[0].input();
                    }
                    else gs[0].readmds(smi);

					cout<<"Do exhaustive decyclization"<<endl;
					exhaustive_decyclization(gs[0]);
				}

				for (int j=i;j<nmol;j++) {
					string smi2=""; //gs[1].molesmi;

					if (para.round==1) smi2=rd_1molsmi(para.guess,j);
					else smi2=rd_1molsmi((para.logdir+"mol-"+num.str()+"gen.txt"),j);

					if (para.operation=="Crossover" || para.operation=="All") {
	                    if (para.round==1) {
    	                    gs[0].molesmi=gs[0].smiles=smi;   gs[0].input();
							gs[1].molesmi=gs[1].smiles=smi2;   gs[1].input();
            	        }
                	    else {
							gs[0].readmds(smi);
							gs[1].readmds(smi2);
						}

						cout<<"Do exhaustive crossover"<<endl;
						exhaustive_crossover(gs[0],gs[1]);
					}

					if (para.operation=="Combination" || para.operation=="All") {
                        if (para.round==1) {
                            gs[0].molesmi=gs[0].smiles=smi;   gs[0].input();
                            gs[1].molesmi=gs[1].smiles=smi2;   gs[1].input();
                        }
                        else {
                            gs[0].readmds(smi);
                            gs[1].readmds(smi2);
                        }

						cout<<"Do exhaustive combination"<<endl;
						exhaustive_combination(gs[0],gs[1]);
					}
				}

			}
			para.glbouf.close();

			para.round++;

		}

		cout<<"See the result for each txt file. (combination.txt, crossover.txt, addition.txt, subtraction.txt, change_ele.txt, change_bnd.txt, cyclization.txt, decyclization.txt, insertion.txt, change_cistrans.txt, change_chirality.txt)"<<endl;		
	}


	if (para.ion) {

		cout<<"Read the ILs."<<endl;


		stringstream num(""),num1("");
		while (para.round<=para.epoch) {

			int ncat=0,nan=0;
			if (gsion!=NULL) {
				delete [] gsion;
				gsion=NULL;
			}
			gsion=new IL [2];
			gsion[0].ion[0].data=&data;
			gsion[0].ion[1].data=&data;
			gsion[1].ion[0].data=&data;
			gsion[1].ion[1].data=&data;

			num.str("");
			num << para.round-1;

			if (para.round==1) ncat=cal_nIL(para.guess);
			else ncat=cal_nmol((para.logdir+"cat-"+num.str()+"gen.txt"));

			if (para.round==1) nan=cal_nIL(para.guess);
			else nan=cal_nmol((para.logdir+"an-"+num.str()+"gen.txt"));

			num1.str("");
			num1 << para.round;
			para.glbouf.open((para.logdir+"cat-"+num1.str()+"gen.txt").c_str());
			para.glbouf1.open((para.logdir+"an-"+num1.str()+"gen.txt").c_str());

			for (int i=0;i<ncat;i++) {
				for (int j=0;j<nan;j++) {
					string smicat="",smian=""; //gsion[0].ion[0].molesmi //gsion[0].ion[1].molesmi

					if (para.round==1) rd_1ILsmi(para.guess,i,j,smicat,smian);
					else {
						smicat=rd_1molsmi((para.logdir+"cat-"+num.str()+"gen.txt"),i);
						smian=rd_1molsmi((para.logdir+"an-"+num.str()+"gen.txt"),j);
					}

					if (para.operation=="Cistrans_inversion" || para.operation=="All") {
						if (para.round==1) {
							gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
							gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
							gsion[0].input();
						}
						else {
							gsion[0].ion[0].readmds(smicat);
							gsion[0].ion[1].readmds(smian);
						}

						cout<<"Do exhaustive change_cistrans"<<endl;
						exhaustive_change_cistrans1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					if (para.operation=="Chirality_inversion" || para.operation=="All") {
                        if (para.round==1) {
                            gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
                            gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                            gsion[0].input();
                        }
                        else {
                            gsion[0].ion[0].readmds(smicat);
                            gsion[0].ion[1].readmds(smian);
                        }

						cout<<"Do exhaustive change_chirality"<<endl;
						exhaustive_change_chirality1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					if (para.operation=="Subtraction" || para.operation=="All") {
                        if (para.round==1) {
                            gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
                            gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                            gsion[0].input();
                        }
                        else {
                            gsion[0].ion[0].readmds(smicat);
                            gsion[0].ion[1].readmds(smian);
                        }

						cout<<"Do exhaustive subtraction"<<endl;
						exhaustive_subtraction1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					if (para.operation=="Combination" || para.operation=="All") {
                        if (para.round==1) {
                            gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
                            gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                            gsion[0].input();
                        }
                        else {
                            gsion[0].ion[0].readmds(smicat);
                            gsion[0].ion[1].readmds(smian);
                        }

						cout<<"Do exhaustive addition"<<endl;
						exhaustive_addition1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					if (para.operation=="Insertion" || para.operation=="All") {
                        if (para.round==1) {
                            gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
                            gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                            gsion[0].input();
                        }
                        else {
                            gsion[0].ion[0].readmds(smicat);
                            gsion[0].ion[1].readmds(smian);
                        }

						cout<<"Do exhaustive insertion"<<endl;
						exhaustive_insertion1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					if (para.operation=="Change_ele" || para.operation=="All") {
                        if (para.round==1) {
                            gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
                            gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                            gsion[0].input();
                        }
                        else {
                            gsion[0].ion[0].readmds(smicat);
                            gsion[0].ion[1].readmds(smian);
                        }

						cout<<"Do exhaustive change_ele"<<endl;
						exhaustive_change_ele1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					if (para.operation=="Change_bnd" || para.operation=="All") {
                        if (para.round==1) {
                            gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
                            gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                            gsion[0].input();
                        }
                        else {
                            gsion[0].ion[0].readmds(smicat);
                            gsion[0].ion[1].readmds(smian);
                        }

						cout<<"Do exhaustive change_bnd"<<endl;
						exhaustive_change_bnd1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					if (para.operation=="Cyclization" || para.operation=="All") {
                        if (para.round==1) {
                            gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
                            gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                            gsion[0].input();
                        }
                        else {
                            gsion[0].ion[0].readmds(smicat);
                            gsion[0].ion[1].readmds(smian);
                        }

						cout<<"Do exhaustive cyclization"<<endl;
						exhaustive_cyclization1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					if (para.operation=="Decyclization" || para.operation=="All") {
                        if (para.round==1) {
                            gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;
                            gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                            gsion[0].input();
                        }
                        else {
                            gsion[0].ion[0].readmds(smicat);
                            gsion[0].ion[1].readmds(smian);
                        }

						cout<<"Do exhaustive decyclization"<<endl;
						exhaustive_decyclization1(gsion[0]);
						//system("rm ./mds/* 2> /dev/null");
					}

					for (int k=i+1;k<ncat;k++) {
						for (int l=j+1;l<nan;l++) {

							string smicat2="",smian2=""; //gsion[1].ion[0].molesmi //gsion[1].ion[0].molesmi

							if (para.round==1) rd_1ILsmi(para.guess,k,l,smicat2,smian2);
							else {
								smicat2=rd_1molsmi((para.logdir+"cat-"+num.str()+"gen.txt"),k);
								smian2=rd_1molsmi((para.logdir+"an-"+num.str()+"gen.txt"),l);
							}

							if (para.operation=="Crossover" || para.operation=="All") {
								if (para.round==1) {
									gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;    gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
									gsion[1].ion[0].molesmi=gsion[1].ion[0].smiles=smicat2;   gsion[1].ion[1].molesmi=gsion[1].ion[1].smiles=smian2;
									gsion[0].input();
									gsion[1].input();
								}
								else {
									gsion[0].ion[0].readmds(smicat);    gsion[0].ion[1].readmds(smian);
									gsion[1].ion[0].readmds(smicat2);   gsion[1].ion[1].readmds(smian2);
								}

								cout<<"Do exhaustive crossover"<<endl;
								exhaustive_crossover1(gsion[0],gsion[1]);
							}

							if (para.operation=="Combination" || para.operation=="All") {
                                if (para.round==1) {
                                    gsion[0].ion[0].molesmi=gsion[0].ion[0].smiles=smicat;    gsion[0].ion[1].molesmi=gsion[0].ion[1].smiles=smian;
                                    gsion[1].ion[0].molesmi=gsion[1].ion[0].smiles=smicat2;   gsion[1].ion[1].molesmi=gsion[1].ion[1].smiles=smian2;
                                    gsion[0].input();
                                    gsion[1].input();
                                }
                                else {
                                    gsion[0].ion[0].readmds(smicat);    gsion[0].ion[1].readmds(smian);
                                    gsion[1].ion[0].readmds(smicat2);   gsion[1].ion[1].readmds(smian2);
                                }

								cout<<"Do exhaustive combination"<<endl;
								exhaustive_combination1(gsion[0],gsion[1]);
							}
						}
					}
				}
			}
			para.glbouf.close();
			para.glbouf1.close();

			para.round++;

		}


		cout<<"See the result for each txt file. (combination_IL.txt, crossover_IL.txt, addition_IL.txt, insertion_IL.txt, subtraction_IL.txt, cyclization_IL.txt, decyclization_IL.txt, change_ele_IL, change_bnd_IL, change_cistrans_IL.txt, and change_chirality_IL.txt)"<<endl;
	}


	para.stat.close();


	if (1) {
        if (gs!=NULL) {
            delete [] gs;
            gs=NULL;
        }
        //gs=new MOLECULE [1];
        //gs[0].data=&data;
        gs=new MOLECULE [2];
        gs[1].data=gs[0].data=&data;

		cout<<"Make benzene from methane (path 1)"<<endl;
		CH4_to_Bz_path1(gs[0]);

		cout<<"Make benzene from methane (path 2)"<<endl;
		CH4_to_Bz_path2(gs[0]);

		cout<<"Make aspirin from methane"<<endl;
		CH4_to_aspirin(gs[0]);

		cout<<"Make alpha-carotene from methane (path 1)"<<endl;
		CH4_to_alpha_carotene_path1(gs[0]);

		cout<<"Make alpha-carotene from methane (path 2)"<<endl;
		CH4_to_alpha_carotene_path2(gs[0]);

		cout<<"Make alpha-carotene from methane (path 3)"<<endl;
		CH4_to_alpha_carotene_path3(gs[0]);

        cout<<"Tamiflu total synthesis"<<endl;
        Tamiflu_Corey(gs[0],gs[1]);


        cout<<"See the result for each txt file. (benzene_path1.txt, benzene_path2.txt, aspirin.txt, carotene_path1.txt, carotene_path2.txt, and carotene_path3.txt)"<<endl;

	}




	//gs[0].data=gs[1].data=NULL;   
	delete [] gs;
	gs=NULL;

	//gsion[0].ion[0].data=gsion[0].ion[1].data=NULL;
	//gsion[1].ion[0].data=gsion[1].ion[1].data=NULL;
	delete [] gsion;
	gsion=NULL;

	return 1;
}


