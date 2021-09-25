#include "UTILITY.h"
using namespace std;

extern PARAMETER para;

unsigned int rd_para(char **&argv){
    ifstream inf(argv[1],ios::in);
    inf >> ws;
    while (!inf.eof()) {
        string tmp="null";
        unsigned int cursor=inf.tellg();
        inf >> tmp >> ws;
        if (tmp[0]=='#') {
            inf.seekg(cursor);
            getline(inf,tmp);
            inf >> ws;
            tmp="null";
        }
        if (tmp=="CHEMICAL_IUPUTLIST") inf >> para.guess >> ws;
		else if (tmp=="OPERATIONS") inf >> para.operation >> ws;
        else if (tmp=="LOG_DIR") inf >> para.logdir >> ws;
		else if (tmp=="PROGRAM_DIR") inf >> para.programdir >> ws;
        else if (tmp=="IF_ION") inf >> para.ion >> ws;
        else if (tmp=="IF_PROTECT") inf >> para.protect >> ws;
		else if (tmp=="IF_ENUMERATION") inf >> para.enumeration >> ws;
		else if (tmp=="MDS_DIR") inf >> para.mdsdir >> ws;
		else if (tmp=="IF_OUTPUT_MDS") inf >> para.ifwritemds >> ws;
		else if (tmp=="EPOCH") inf >> para.epoch >> ws;
		else if (tmp=="ELEMENT_LIST") inf >> para.element_list >> ws;
    }
    inf.close();

    return 1;
}

unsigned int output_para(){
    cout << "=====================================================================" << endl;
    cout << left << setw(25) << "CHEMICAL_IUPUTLIST" << " " << left << setw(45) << para.guess << endl
        << left << setw(25) << "OPERATIONS" << " " << left << setw(45) << para.operation << endl
        << left << setw(25) << "LOG_DIR" << " " << left << setw(45) << para.logdir << endl
        << left << setw(25) << "PROGRAM_DIR" << " " << left << setw(45) << para.programdir << endl
        << left << setw(25) << "MDS_DIR" << " " << left << setw(45) << para.mdsdir << endl
        << left << setw(25) << "ELEMENT_LIST" << " " << left << setw(45) << para.element_list << endl
		<< left << setw(25) << "EPOCH" << " " << left << setw(45) << para.epoch << " " << endl
        << left << setw(25) << "IF_ION" << " " << left << setw(45) << para.ion << endl
        << left << setw(25) << "IF_PROTECT" << " " << left << setw(45) << para.protect << " " << endl
		<< left << setw(25) << "IF_OUTPUT_MDS" << " " << left << setw(45) << para.ifwritemds << " " << endl
        << left << setw(25) << "IF_ENUMERATION" << " " << left << setw(45) << para.enumeration << " " << endl;
    cout << "=====================================================================" << endl;

    return 1;
}


unsigned int rd_ILs(IL *&gsion,POOL *pol){
	if (para.ion) {
		ifstream inf(para.guess.c_str());
        string anion="";
        string cation="";		
        
        unsigned int i=0,cur=0;
        inf >> ws;
        while (!inf.eof()) {
            anion="";
            cation="";
            cur=inf.tellg();
            inf >> cation >> ws >> anion >> ws;
            if (anion!="" && anion[0]!='#' && cation!="" && cation[0]!='#') i++;
            else {
            	inf.seekg(cur);
            	getline(inf,anion);
            	inf >> ws;
			}
        }
        inf.close();
        
        gsion=new IL [i];
        para.gssize=i;

		for (unsigned int k=0;k<para.gssize;k++) {
			for (unsigned int k1=0;k1<gsion[k].nsubcomp;k1++) {
				gsion[k].ion[k1].data=pol;
			}
		}
        
        inf.open(para.guess.c_str());
		i=0;
        inf >> ws;
        while (!inf.eof()) {
            anion="";
            cation="";
            cur=inf.tellg();
            inf >> cation >> ws >> anion >> ws;
            if (anion!="" && anion[0]!='#' && cation!="" && cation[0]!='#') {
                cout << "Current step : read ionic liquid " << (i+1) << " , " << cation << " " << anion << endl;
                gsion[i].ion[0].smiles=gsion[i].ion[0].molesmi=cation;
                gsion[i].ion[1].smiles=gsion[i].ion[1].molesmi=anion;
                gsion[i].ion[0].input();
                gsion[i].ion[1].input();
                if (para.protect) gsion[i].prct();
                i++;
                if (i>=para.gssize) break;
            }
            else {
            	inf.seekg(cur);
            	getline(inf,anion);
            	inf >> ws;            	
			}
        }
        inf.close();		
	}

	return 1;
}

unsigned int rd_mols(MOLECULE *&gs,POOL *pol){
	if (!para.ion) {
		ifstream inf(para.guess.c_str());
        string smi="";
        
        unsigned int i=0,cur=0;
        inf >> ws;
        while (!inf.eof()) {
            smi="";
            cur=inf.tellg();
            inf >> smi >> ws;
            if (smi!="" && smi[0]!='#') i++;
            else {
            	inf.seekg(cur);
            	getline(inf,smi);
            	inf >> ws;
			}
        }
        inf.close();
        
        gs=new MOLECULE [i];
        para.gssize=i;

        for (unsigned int k=0;k<para.gssize;k++) {
            gs[k].data=pol;
        }
        
        inf.open(para.guess.c_str());
		i=0;
        inf >> ws;
        while (!inf.eof()) {
            smi="";
            cur=inf.tellg();
            inf >> smi >> ws;
            if (smi!="" && smi[0]!='#') {
                cout << "Current step : read molecule " << (i+1) << " , " << smi << endl;
                gs[i].smiles=gs[i].molesmi=smi;
                gs[i].input();
                if (para.protect) gs[i].prct();
                i++;
                if (i>=para.gssize) break;
            }
            else {
            	inf.seekg(cur);
            	getline(inf,smi);
            	inf >> ws;            	
			}
        }
        inf.close();		
	}

	return 1;
}

unsigned int cal_nmol(string r_inputlist) {
    string smi="";
    ifstream inf((r_inputlist).c_str());
    unsigned int nmol=0;
    inf >> ws;
    while (inf.is_open() && !inf.eof()) {
    	smi="";
		unsigned int cur=inf.tellg();
        inf >> smi >> ws;
        if (smi!="" && smi[0]!='#') nmol++;
		else {
            inf.seekg(cur);
            getline(inf,smi);
            inf >> ws;
		}
    }
    inf.close();

	return nmol;
}

unsigned int cal_nIL(string r_inputlist) {
    string smi="",smi1="";
    ifstream inf((r_inputlist).c_str());
    unsigned int nmol=0;
    inf >> ws;
    while (inf.is_open() && !inf.eof()) {
        smi="";
        unsigned int cur=inf.tellg();
        inf >> smi >> ws >> smi1 >> ws;
        if (smi!="") {
			if (smi[0]!='#') {
				inf.seekg(cur);
				inf >> smi >> ws >> smi1 >> ws;
				nmol++;
			}
	        else {
    	        inf.seekg(cur);
        	    getline(inf,smi);
            	inf >> ws;
        	}
		}
        else {
            inf.seekg(cur);
            getline(inf,smi);
            inf >> ws;
        }
    }
    inf.close();

    return nmol;
}


string rd_1molsmi(string r_inputlist,unsigned int num) {
	string smi="";
    ifstream inf((r_inputlist).c_str());
	inf >> ws;
    for (unsigned int i=0;i<=num;i++) {
		if (inf.eof()) break;
        smi="";
		unsigned int cur=inf.tellg();
        inf >> smi >> ws;
        if (smi=="" || smi[0]=='#') {
            inf.seekg(cur);
            getline(inf,smi);
            inf >> ws;
			i--;
			smi="";
        }
    }
    inf.close();

	cout << "Current step : read molecule " << (num+1) << " , " << smi << endl;

	return smi;

}

unsigned int rd_1ILsmi(string r_inputlist,unsigned int numcat,unsigned int numan,string &catsmi,string &ansmi) {
    string smi="",smi1="";
    ifstream inf((r_inputlist).c_str());
    inf >> ws;
    for (unsigned int i=0;i<=numcat;i++) {
        if (inf.eof()) break;
        smi="";
		smi1="";
        unsigned int cur=inf.tellg();
        inf >> smi >> ws >> smi1 >> ws;
        if (smi=="" || smi[0]=='#') {
            inf.seekg(cur);
            getline(inf,smi);
            inf >> ws;
            i--;
            smi="";
			smi1="";
        }
    }
    inf.close();

	catsmi=smi;

	cout << "Current step : read cation " << (numcat+1) << " , " << catsmi << endl;

	smi="";
    smi1="";

    inf.open((r_inputlist).c_str());
    inf >> ws;
    for (unsigned int i=0;i<=numan;i++) {
        if (inf.eof()) break;
        smi="";
        smi1="";
        unsigned int cur=inf.tellg();
        inf >> smi >> ws >> smi1 >> ws;
        if (smi=="" || smi[0]=='#') {
            inf.seekg(cur);
            getline(inf,smi);
            inf >> ws;
            i--;
            smi="";
            smi1="";
        }
    }
    inf.close();

    ansmi=smi1;

	cout << "Current step : read anion " << (numan+1) << " , " << ansmi << endl;

    return 1;

}

unsigned int clearlog() {
	system(("rm "+para.logdir+"combination.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"crossover.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"addition.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"subtraction.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"insertion.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"cyclization.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"decyclization.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"change_ele.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"change_bnd.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"change_cistrans.txt 2> /dev/null").c_str());
	system(("rm "+para.logdir+"change_chirality.txt 2> /dev/null").c_str());

    system(("rm "+para.logdir+"combination_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"crossover_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"addition_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"subtraction_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"insertion_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"cyclization_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"decyclization_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"change_ele_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"change_bnd_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"change_cistrans_IL.txt 2> /dev/null").c_str());
    system(("rm "+para.logdir+"change_chirality_IL.txt 2> /dev/null").c_str());

	return 1;
}


unsigned int mark(ostream& out) {
	out <<"***************************************************************************************"<<endl
        <<"               Molecular Data Structure for Computer-aided Molecular Design            "<<endl
        <<"                                      Hsuan Hao Hsu                                    "<<endl
        <<"                                    Chen Hsuan Huang                                   "<<endl
        <<"                            Shiang Tai Lin (stlin@ntu.edu.tw)                          "<<endl
        <<"                      Computational Molecular Engineering Laboratory                   "<<endl
        <<"      Department of Chemical Engineering, National Taiwan University, Taipei, Taiwan   "<<endl
        <<"                                   Copyright (c) 2019                                  "<<endl
        <<"                                   All rights Reserved                                 "<<endl
        <<"***************************************************************************************"<<endl;
	return 1;
}


//get filename extensions of files
string getExt(string pathname) {
    int period = pathname.find_last_of(".");   //Finds the last period character of the string
    string ext = pathname.substr(period+1);   //use +1 because the period is not really needed
    return ext;
}


//get filenames in a directory
unsigned int getdir(string dir, vector<string> &filenames, string fileextention) {
    DIR *dp=NULL;     //declare a pointer for a directory
    struct dirent *dirp=NULL;
    if ((dp = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }
    while ((dirp = readdir(dp)) != NULL) {     //If dirent pointer is not NULL
        string name = (dirp->d_name) ;
        if (getExt(name) == fileextention && name!=("."+fileextention) && name!=(fileextention)) {
            filenames.push_back(name);     //store filenames into vector
        }
    }
	dirp=NULL;
	dp=NULL;
    closedir(dp);     //close directory pointer
    return 1;
}

//make dat list
unsigned int mk_datlist() {
	string Ldir=para.mdsdir+"DATLIST.txt";

    ofstream LIS(Ldir.c_str());
    DIR *dp=NULL;     //declare a pointer for a directory
    struct dirent *dirp=NULL;
    if ((dp = opendir((para.mdsdir).c_str())) == NULL) { //para.programdir+"mds/"  //para.mdsdir
        cout << "Error(" << errno << ") opening " << (para.mdsdir) << endl; //(para.programdir+"mds/")  //para.mdsdir
        return errno;
    }
    while ((dirp = readdir(dp)) != NULL) {     //If dirent pointer is not NULL
        string name = (dirp->d_name) ;
        if (getExt(name)=="enc" && name!=".enc" && name!="enc") {
			string SMI="";
			for (unsigned int i=0;i<name.find_first_of('_');i++) SMI+=name[i];
            if (SMI!="" && SMI!="null") LIS << setw(50) << left << SMI << "   " << left << (para.mdsdir+name) << endl; //(para.programdir+"mds/"+name)  //para.mdsdir+name
        }
    }
    dirp=NULL;
    dp=NULL;
    closedir(dp);     //close directory pointer
    LIS.close();

    cout << "DATLIST.txt established." << endl;

    return 1;
}


unsigned int smi2mol(string SMI,OBMol &mol) {
    stringstream ss1("");
    if (1) {
        stringstream ss("");
        ss << SMI;
        OpenBabel::OBConversion conv(&ss,&ss1);
        if (conv.SetInAndOutFormats("SMI","MOL")) {
            conv.AddOption("gen3D", OpenBabel::OBConversion::GENOPTIONS);
            //conv.AddOption("align", OBConversion::GENOPTIONS);
            //conv.AddOption("canonical", OpenBabel::OBConversion::GENOPTIONS);
            conv.AddOption("3", OpenBabel::OBConversion::OUTOPTIONS);
            conv.Convert();
            //cout <<"ss1: " << endl
            //  << ss1.str() << endl;
        }
    }
    if (1) {
        OpenBabel::OBConversion conv(&ss1);
        conv.SetInFormat("MOL");
        conv.Read(&mol);
    }

    return 1;
}


string mol2smi(OBMol &mol) {
    stringstream ss1("");
    stringstream ss("");
    stringstream buf("");
    for (unsigned int i=0;i<mol.NumHvyAtoms();i++) {
        if (i<=mol.NumHvyAtoms()-2) buf << (i+1) << "-";
        else buf << (i+1);
    }

    if (1) {
        OpenBabel::OBConversion conv;
        //conv.SetInFormat("MOL");
        //conv.Read(&mol);
        conv.SetOutFormat("MOL");
        conv.AddOption("3", OpenBabel::OBConversion::OUTOPTIONS);
        conv.Write(&mol,&ss1);
        //cout <<"ss1: " << endl
        //   << ss1.str() << endl;
    }
    if (1) {
        OpenBabel::OBConversion conv(&ss1,&ss);
        if (conv.SetInAndOutFormats("MOL","SMI")) {
            conv.AddOption("o", OpenBabel::OBConversion::OUTOPTIONS , buf.str().c_str());
            //conv.AddOption("gen3D", OpenBabel::OBConversion::GENOPTIONS);
            //conv.AddOption("align", OBConversion::GENOPTIONS);
            //conv.AddOption("canonical", OpenBabel::OBConversion::GENOPTIONS);
            //conv.AddOption("3", OpenBabel::OBConversion::OUTOPTIONS);
            conv.Convert();
            //cout <<"ss: " << ss.str() << endl;

        }
    }
    //if (1) {
    //    OBConversion conv(&ss1);
    //    conv.SetInFormat("MOL");
    //    conv.Read(&mol);
    //}

    return ss.str();
}

unsigned int SMI_Enumerator(string SMI) {
    OBMol mol;
    smi2mol(SMI,mol);

    string mm=SMI;
    for (unsigned int q=0;q<mm.length();q++) {
        if (mm[q]=='/') mm[q]='u';
        if (mm[q]=='\\') mm[q]='d';
		if (mm[q]=='*') mm[q]='x';
    }

	//ofstream ouf1(para.programdir+"mds/"+mm+".synonyms");
	ofstream ouf1(para.mdsdir+mm+".synonyms");

	ouf1 << "Original SMILES: " << SMI << endl;

    vector<int> nonH_ord(0);
    vector<int> H_ord(0);
    vector<int> ord(0);
    vector<string> memo(0);
    memo.push_back(SMI);
    unsigned int natom=0;

    FOR_ATOMS_OF_MOL(atom, mol) {
        if (atom->GetAtomicNum()!=1 && atom->GetAtomicNum()!=0) {
            nonH_ord.push_back(atom->GetIdx()); //atom->GetIdx()
            natom++;
        }
        else if (atom->GetAtomicNum()==1 || atom->GetAtomicNum()==0) {
            H_ord.push_back(atom->GetIdx()); //atom->GetIdx()
        }
    }
    nonH_ord.reserve(nonH_ord.size());
    H_ord.reserve(H_ord.size());

    ord.resize(0);
    ord.reserve(nonH_ord.size()+H_ord.size());
    if (0) {
        if (1) {
            sort(nonH_ord.begin(),nonH_ord.end());
            do {
                ord.resize(0);
                for (unsigned int i=0;i<nonH_ord.size();i++) ord.push_back(nonH_ord.at(i));
                for (unsigned int i=0;i<H_ord.size();i++) ord.push_back(H_ord.at(i));

                mol.RenumberAtoms(ord);
                string buf=mol2smi(mol);

                //cout << buf << endl;
                ouf1 << buf << endl;
            } while ( next_permutation(nonH_ord.begin(),nonH_ord.end()) ); //ord.begin(),ord.end()
        }
        if (0) {
            for (unsigned int i=0;i<(nonH_ord.size()+H_ord.size());i++) ord.push_back(i+1);
            do {
                mol.RenumberAtoms(ord);

                cout << mol2smi(mol) << endl;
            } while ( next_permutation(nonH_ord.begin(),nonH_ord.end()) ); //ord.begin(),ord.end()

        }
    }
    if (1) {
        long long unsigned int N_arrange=tgamma(nonH_ord.size()+1),ct=0,n_repeat=0;
        unsigned int digits = 0;
        if (nonH_ord.size()>=2) digits = log10(N_arrange)+1;
        if (0) {
            for (unsigned int i=0;i<(nonH_ord.size()+H_ord.size());i++) ord.push_back(i+1);
        }

        for (long long unsigned int i=1;i<=N_arrange;i++) {
            if (1) {
                random_shuffle(nonH_ord.begin(),nonH_ord.end());
                ord.resize(0);
                for (unsigned int j=0;j<nonH_ord.size();j++) ord.push_back(nonH_ord.at(j));
                for (unsigned int j=0;j<H_ord.size();j++) ord.push_back(H_ord.at(j));
            }
            if (0) {
                random_shuffle(ord.begin(),ord.end());
            }

            mol.RenumberAtoms(ord);
            string buf=mol2smi(mol);

            for (unsigned int j=0;j<memo.size();j++) {
                if (memo.at(j)==buf) {
                    //i--;
                    //i=memo.size()-1;
                    n_repeat++;
                    break;
                }
                else if (memo.at(j)!=buf && j>=memo.size()-1) {
                    memo.push_back(buf);
                    //cout << setw(digits+1) << (ct+1) << " " << setw(80) << buf ;
                    ouf1 << setw(digits+1) << (ct+1) << " " << setw(80) << buf ;
                    ct++;
                    //mol.RenumberAtoms(ori_ord);
                    n_repeat=0;

                    if (0) for (unsigned int k=0;k<memo.size();k++) cout << memo.at(k) << endl;
                }
            }

            if (n_repeat>=pow(nonH_ord.size(),3)) break;
        }
    }

    ouf1.close();


	return 1;
}
