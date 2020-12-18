#include "ATOM.h"
#include "MOLECULE.h"
#include "MATRIX.h"
#include "OPT.h"
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
	
    MOLECULE *gs=NULL;
	IL *gsion=NULL;

	POOL data;
	data.set_up();
    
	mk_datlist(); 
   
	if (!para.ion) {
		cout<<"Read the molecules."<<endl;    
		rd_mols(gs,&data);
		string smi1=gs[0].molesmi,smi2=gs[1].molesmi;

        cout<<"Do exhaustive exchange"<<endl;
        exhaustive_exchange(gs[0]);

		gs[0].smiles=smi1;  gs[0].input();
		gs[1].smiles=smi2;  gs[1].input();
        cout<<"Do exhaustive addition"<<endl;
        exhaustive_addition(gs[0]);
	    
        gs[0].smiles=smi1;  gs[0].input();
        gs[1].smiles=smi2;  gs[1].input();
		cout<<"Do exhaustive combination"<<endl;
		exhaustive_combination(gs[0],gs[1]);
    	
        gs[0].smiles=smi1;  gs[0].input();
        gs[1].smiles=smi2;  gs[1].input();
		cout<<"Do exhaustive crossover"<<endl;
		exhaustive_crossover(gs[0],gs[1]);
		
        gs[0].smiles=smi1;  gs[0].input();
        gs[1].smiles=smi2;  gs[1].input();
		cout<<"Do exhaustive subtraction"<<endl;
		exhaustive_subtraction(gs[0]);
		
		
		//cout<<"Do exhaustive addition"<<endl;
		//exhaustive_addition(gs[0]);
		
		
		//cout<<"Do exhaustive exchange"<<endl;
		//exhaustive_exchange(gs[0]);
				
        gs[0].smiles=smi1;  gs[0].input();
        gs[1].smiles=smi2;  gs[1].input();
		cout<<"Do exhaustive cyclization"<<endl;
		exhaustive_cyclization(gs[0]);
		
	
        gs[0].smiles=smi1;  gs[0].input();
        gs[1].smiles=smi2;  gs[1].input();
		cout<<"Make benzene, c1ccccc1 from methane, C"<<endl;
		CH4_to_Bz(gs[0]);
	
	
		cout<<"Mission completed!"<<endl;
		cout<<"See the result for each txt file. (combination.txt, crossover.txt, addition.txt, subtraction.txt, circle.txt, and benzene.txt)"<<endl;		
	}
	
	if (para.ion) {
		cout<<"Read the ILs."<<endl;    
		rd_ILs(gsion,&data);
		string smi1_cat=gsion[0].ion[0].molesmi,smi1_an=gsion[0].ion[1].molesmi;    
		string smi2_cat=gsion[1].ion[0].molesmi,smi2_an=gsion[1].ion[1].molesmi;

        cout<<"Do exhaustive combination"<<endl;
        exhaustive_combination1(gsion[0],gsion[1]);

        gsion[0].ion[0].smiles=smi1_cat;    gsion[0].ion[1].smiles=smi1_an;
        gsion[1].ion[0].smiles=smi2_cat;    gsion[1].ion[1].smiles=smi2_an;
        gsion[0].input();                    gsion[1].input();
        cout<<"Do exhaustive change of cis/trans"<<endl;
        exhaustive_changect1(gsion[0]);

		gsion[0].ion[0].smiles=smi1_cat;    gsion[0].ion[1].smiles=smi1_an;
		gsion[1].ion[0].smiles=smi2_cat;    gsion[1].ion[1].smiles=smi2_an;
		gsion[0].input();					 gsion[1].input();
		cout<<"Do switch"<<endl;
		swit(gsion[0],gsion[1]);

        gsion[0].ion[0].smiles=smi1_cat;    gsion[0].ion[1].smiles=smi1_an;
        gsion[1].ion[0].smiles=smi2_cat;    gsion[1].ion[1].smiles=smi2_an;
        gsion[0].input();                    gsion[1].input();
        cout<<"Do exhaustive crossover"<<endl;
        exhaustive_crossover1(gsion[0],gsion[1]);

        gsion[0].ion[0].smiles=smi1_cat;    gsion[0].ion[1].smiles=smi1_an;
        gsion[1].ion[0].smiles=smi2_cat;    gsion[1].ion[1].smiles=smi2_an;
        gsion[0].input();                    gsion[1].input();
        cout<<"Do exhaustive subtraction"<<endl;
        exhaustive_subtraction1(gsion[0]);

        //gsion[0].ion[0].smiles=smi1_cat;    gsion[0].ion[1].smiles=smi1_an;
        //gsion[1].ion[0].smiles=smi2_cat;    gsion[1].ion[1].smiles=smi2_an;
        //gsion[0].input();                    gsion[1].input();
        //cout<<"Do exhaustive combination"<<endl;
        //exhaustive_combination1(gsion[0],gsion[1]);

        gsion[0].ion[0].smiles=smi1_cat;    gsion[0].ion[1].smiles=smi1_an;
        gsion[1].ion[0].smiles=smi2_cat;    gsion[1].ion[1].smiles=smi2_an;
        gsion[0].input();                    gsion[1].input();		
		cout<<"Do exhaustive addition"<<endl;
		exhaustive_addition1(gsion[0]);
		
		
        gsion[0].ion[0].smiles=smi1_cat;    gsion[0].ion[1].smiles=smi1_an;
        gsion[1].ion[0].smiles=smi2_cat;    gsion[1].ion[1].smiles=smi2_an;
        gsion[0].input();                    gsion[1].input();		
		cout<<"Do exhaustive exchange"<<endl;
		exhaustive_exchange1(gsion[0]);
		
		
        gsion[0].ion[0].smiles=smi1_cat;    gsion[0].ion[1].smiles=smi1_an;
        gsion[1].ion[0].smiles=smi2_cat;    gsion[1].ion[1].smiles=smi2_an;
        gsion[0].input();                    gsion[1].input();		
		cout<<"Do exhaustive cyclization"<<endl;
		exhaustive_cyclization1(gsion[0]);
		

		//cout<<"Make benzene, c1ccccc1 from methane, C"<<endl;
		//CH4_to_Bz(gs[0]);
	
	
		cout<<"Mission completed!"<<endl;
		cout<<"See the result for each txt file. (combination-IL.txt, crossover-IL.txt, addition-IL.txt, subtraction-IL.txt, circle-IL.txt, and benzene-IL.txt)"<<endl;		

		for (int i=0;i<para.gssize;i++) {
			for (int j=0;j<gsion[i].nsubcomp;j++) {
				if (gsion[i].ion[j].data!=NULL) {
					delete gsion[i].ion[j].data;
					gsion[i].ion[j].data=NULL;
				}
			}
		}
	}

    
	delete [] gs;
    gs=NULL;
    
    delete [] gsion;
    gsion=NULL;
    
    return 1;
}


