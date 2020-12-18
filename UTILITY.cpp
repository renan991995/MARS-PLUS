#include "UTILITY.h"
using namespace std;

extern PARAMETER para;

void rd_para(char **&argv){
    ifstream inf(argv[1],ios::in);
    inf >> ws;
    while (!inf.eof()) {
        string tmp="null";
        int cursor=inf.tellg();
        inf >> tmp >> ws;
        if (tmp[0]=='#') {
            inf.seekg(cursor);
            getline(inf,tmp);
            inf >> ws;
            tmp="null";
        }
        if (tmp=="OBABEL_ROUTE") inf >> para.obabel >> ws;
        else if (tmp=="CHEMICAL_IUPUTLIST") inf >> para.guess >> ws;
        else if (tmp=="EXHAUSTIVE_OPERATIONS") inf >> para.operation >> ws;
        else if (tmp=="ION") inf >> para.ion >> ws;
        else if (tmp=="OUTFILE") inf >> para.out_route >> ws;
        else if (tmp=="PROTECT") inf >> para.protect >> ws;
    }
    inf.close();

    return;
}


void rd_ILs(IL *&gsion,POOL *pol){
	if (para.ion) {
		ifstream inf(para.guess.c_str());
        string anion="";
        string cation="";		
        
        int i=0,cur=0;
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

		for (int k=0;k<para.gssize;k++) {
			for (int k1=0;k1<gsion[k].nsubcomp;k1++) {
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
                gsion[i].ion[0].smiles=cation;
                gsion[i].ion[1].smiles=anion;
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
}

void rd_mols(MOLECULE *&gs,POOL *pol){
	if (!para.ion) {
		ifstream inf(para.guess.c_str());
        string smi="";
        
        int i=0,cur=0;
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

        for (int k=0;k<para.gssize;k++) {
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
                gs[i].smiles=smi;
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
}


void mark(ostream& out) {
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
	return;
}

/*
void Obabel_conv() {
	// Create molecule from SMILES string
	OBMol mol;
	stringstream ss("c1ccccc1"),sout;
	OBConversion conv(&ss,&sout);
	//if (conv.SetInFormat("smi") && conv.Read(&mol))
	// Conversion without manipulation
	//OBConversion conv(&is, &os);
	if (conv.SetInAndOutFormats("SMI", "MOL"))
	{
 		// Option "h" adds explicit hydrogens
 		//conv.AddOption("h", OBConversion::GENOPTIONS);
		conv.AddOption("gen3D", OBConversion::GENOPTIONS);
		conv.AddOption("canonical", OBConversion::GENOPTIONS);
		conv.AddOption("3", OBConversion::OUTOPTIONS);
 		conv.Convert();
	}
	cout << sout.str() << endl;

}
*/


//get filename extensions of files
string getExt(string pathname) {
    int period = pathname.find_last_of(".");   //Finds the last period character of the string
    string ext = pathname.substr(period+1);   //use +1 because the period is not really needed
    return ext;
}


//get filenames in a directory
int getdir(string dir, vector<string> &filenames, string fileextention) {
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
    return 0;
}

//make dat list
int mk_datlist() {
	string Ldir=para.smidir+"mds/DATLIST.txt";

    ofstream LIS(Ldir.c_str());
    DIR *dp=NULL;     //declare a pointer for a directory
    struct dirent *dirp=NULL;
    if ((dp = opendir((para.smidir+"mds/").c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << (para.smidir+"mds/") << endl;
        return errno;
    }
    while ((dirp = readdir(dp)) != NULL) {     //If dirent pointer is not NULL
        string name = (dirp->d_name) ;
        if (getExt(name)=="enc" && name!=".enc" && name!="enc") {
			string SMI="";
			for (int i=0;i<name.find_first_of('_');i++) SMI+=name[i];
            if (SMI!="" && SMI!="null") LIS << setw(50) << left << SMI << "   " << left << (para.smidir+"mds/"+name) << endl;
        }
    }
    dirp=NULL;
    dp=NULL;
    closedir(dp);     //close directory pointer
    LIS.close();

    cout << "DATLIST.txt established." << endl;

    return 0;
}

