#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/bond.h>
#include <openbabel/bondtyper.h>
#include <openbabel/ring.h>
#include <openbabel/atom.h>
#include <openbabel/data.h>
#include <openbabel/data_utilities.h>
#include <openbabel/obiter.h>
//#include <openbabel/chiral.h>
#include <openbabel/distgeom.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/stereo/tetranonplanar.h>
#include <openbabel/stereo/tetraplanar.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/builder.h>
#include <openbabel/typer.h>
#include <openbabel/tokenst.h>
#include <openbabel/oberror.h>
//#include <openbabel/elements.h>
#include <openbabel/obmolecformat.h>
//#include <openbabel/stereo/bindings.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
//#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>


#include <bits/stdc++.h>
#include <sstream>
#include <cerrno>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>
#include <iterator>
#include <ctime>
#include <cstdlib>
using namespace std;
//using namespace OpenBabel;

unsigned int rd_para(char **&argv,bool &IF_OBABEL, bool &IF_RDKIT, string &SYN_DIR, string &SMI_LIST, unsigned long long int &TOLERANCE);
unsigned int print_para(bool IF_OBABEL, bool IF_RDKIT, string SYN_DIR, string SMI_LIST, unsigned long long int TOLERANCE);
unsigned int smi2mol(string SMI,OpenBabel::OBMol &mol);
string mol2smi(OpenBabel::OBMol &mol);
string canonicalize_smi(string smi);
int myrandom(int i);

int main(int argc, char **argv)
{
	if (argc!=2) {
		cout << "Usage: [excutable_file] [control_file]" << endl
			 << "For example: ./SMI_Enumerator ./control.in" << endl;
		return 0;
	}
	srand(time(0));
	bool IF_OBABEL=1;
	bool IF_RDKIT=1;
	string SYN_DIR="./SYNONYMS/";
	string SMI_LIST="./SMI.txt";
	unsigned long long int TOLERANCE=1e+15;
	rd_para(argv,IF_OBABEL,IF_RDKIT,SYN_DIR,SMI_LIST,TOLERANCE);

	bool TOLERANCE_set=1;
    if (TOLERANCE==1e+15) TOLERANCE_set=0;

	print_para(IF_OBABEL,IF_RDKIT,SYN_DIR,SMI_LIST,TOLERANCE);

	system(("mkdir "+SYN_DIR+" 2> /dev/null").c_str());
	srand(time(NULL));
	ifstream inf(SMI_LIST.c_str());
	inf >> ws;
	while (!inf.eof() && inf.is_open()) {
		string SMI="";
		unsigned int cur=inf.tellg();
		string null="#";
		inf >> null >> ws;
		if (null[0]!='#') SMI=null;
		else {
			inf.seekg(cur);
			getline(inf,null);
			inf >> ws;
			continue;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////

		vector<string> memo(0);
		unsigned int stereoH=0;

		if (IF_OBABEL) {
			std::cout << SMI << " -----------OBABEL-----------" << std::endl;
			OpenBabel::OBMol mol;

			smi2mol(SMI,mol);
			string can_SMI1=canonicalize_smi(mol2smi(mol));

			vector<unsigned int> nonH_ord(0);
			//vector<unsigned int> chiralH(0);
			vector<unsigned int> H_ord(0);
			vector<int> ord(0);
			memo.push_back(SMI);

			OpenBabel::OBStereoFacade facade(&mol);

			FOR_ATOMS_OF_MOL(atom, mol) {
				if (atom->GetAtomicNum()!=1) { // && atom->GetAtomicNum()!=0
					nonH_ord.push_back(atom->GetIdx()); //atom->GetIdx()
					if (1) {
						if (facade.HasTetrahedralStereo(atom->GetId())) {
							OpenBabel::OBTetrahedralStereo *tr=facade.GetTetrahedralStereo(atom->GetId());
							OpenBabel::OBTetrahedralStereo::Config A=tr->GetConfig(OpenBabel::OBStereo::Clockwise,OpenBabel::OBStereo::ViewFrom);
							for (unsigned int i=0;i<A.refs.size();i++) {
								OpenBabel::OBAtom *atm=mol.GetAtomById(A.refs.at(i));
								if (A.refs.at(i)<mol.NumAtoms()) {
									if (atm->GetAtomicNum()==1) {
										nonH_ord.push_back(atm->GetIdx());
										stereoH++;

										//chiralH.push_back(atom->GetIdx());
										//chiralH.push_back(atm->GetIdx());
									}
								}
								atm=NULL;
							}
							tr=NULL;
						}
					}
				}
			}

			if (1) {
				FOR_BONDS_OF_MOL(bond, mol){
					if (bond->GetBondOrder()==2 && facade.HasCisTransStereo(bond->GetId())) {
						OpenBabel::OBCisTransStereo *ct1=facade.GetCisTransStereo(bond->GetId());
						OpenBabel::OBCisTransStereo::Config A=ct1->GetConfig(OpenBabel::OBStereo::ShapeU);
						for (unsigned int g=0;g<A.refs.size();g++) {
							OpenBabel::OBAtom *atm=mol.GetAtomById(A.refs.at(g));
							if (A.refs.at(g)<mol.NumAtoms()) {
								if (atm->GetAtomicNum()==1) {
									nonH_ord.push_back(atm->GetIdx());
									stereoH++;
								}
							}
							atm=NULL;
						}
						ct1=NULL;
					}
				}
			}

			FOR_ATOMS_OF_MOL(atom, mol) {
				if (atom->GetAtomicNum()==1) { // || atom->GetAtomicNum()==0
					if (0) H_ord.push_back(atom->GetIdx());
					if (1) {
						if (nonH_ord.size()) {
							for (unsigned int i=0;i<nonH_ord.size();i++) {
								if (nonH_ord.at(i)==atom->GetIdx()) break;
								else if (nonH_ord.at(i)!=atom->GetIdx() && i>=nonH_ord.size()-1) {
									H_ord.push_back(atom->GetIdx()); //atom->GetIdx()
								}
							}
						}
						else H_ord.push_back(atom->GetIdx());
					}

					//if (chiralH.size()) {
					//    for (unsigned int i=0;i<chiralH.size();i+=2) {
					//        if (chiralH.at(i+1)==atom->GetIdx()) break;
					//        else if (chiralH.at(i+1)!=atom->GetIdx() && (i+1)>=chiralH.size()-1) {
					//            H_ord.push_back(atom->GetIdx()); //atom->GetIdx()
					//        }
					//    }
					//}
					//else H_ord.push_back(atom->GetIdx());
				}
			}


            string a=can_SMI1;
            for (unsigned int i=0;i<a.length();i++) {
                if (a[i]=='/') a[i]='u';
                else if (a[i]=='\\') a[i]='d';
                else if (a[i]=='*') a[i]='x';
            }

			if (!TOLERANCE_set) {
            	if (nonH_ord.size()>2) TOLERANCE=nonH_ord.size()*(nonH_ord.size()-1)*(nonH_ord.size()-2);
            	else TOLERANCE=2;
			}

            ofstream ouf1((SYN_DIR+"Obabel_"+a+".txt").c_str());
            ouf1 << "Original SMILES: " << SMI << endl;
			ouf1 << "TOLERANCE: " << TOLERANCE << endl;
            ouf1 << "Canonical SMILES: " << canonicalize_smi(SMI) << endl;
            ouf1 << "Canonical SMILES (3Dbuilder->smi): " << can_SMI1 << endl;

			nonH_ord.reserve(nonH_ord.size()); //+chiralH.size()/2
			H_ord.reserve(H_ord.size());

			//sort(nonH_ord.begin(),nonH_ord.end());

			ord.resize(0);
			ord.reserve(nonH_ord.size()+H_ord.size()); //+chiralH.size()/2
			if (0) {
				long long unsigned int N_arrange=tgamma(nonH_ord.size()+1),ct=0,n_repeat=0;
				unsigned int digits = 0;
				if (nonH_ord.size()>=2) digits = log10(N_arrange)+1;

				//sort(nonH_ord.begin(),nonH_ord.end());
				do {
					ord.resize(0);
					for (unsigned int i=0;i<nonH_ord.size();i++) ord.push_back(nonH_ord.at(i));
					for (unsigned int i=0;i<H_ord.size();i++) ord.push_back(H_ord.at(i));

					mol.RenumberAtoms(ord);
					string buf=mol2smi(mol); //
					string can_SMI2=canonicalize_smi(buf);

					if (can_SMI1==can_SMI2) {
						for (unsigned int j=0;j<memo.size();j++) {
							if (memo.at(j)==buf) {
								//i--;
								//i=memo.size()-1;
								n_repeat++;
								break;
							}
							else if (memo.at(j)!=buf && j>=memo.size()-1) {
								memo.push_back(buf);
								//stringstream o1("");
								//for (unsigned int k1=0;k1<ord.size();k1++) {
								//    if (k1<=ord.size()-2) o1 << ord.at(k1) << "-";
								//    else o1 << ord.at(k1);
								//}

								cout << setw(digits+1) << (ct+1) << " " << setw(80) << buf << endl; //<< " " << right << setw(80) << o1.str() << endl;
								ouf1 << setw(digits+1) << (ct+1) << " " << setw(80) << buf << endl; //" " << right << setw(80) << o1.str() << endl;
								ct++;
								//mol.RenumberAtoms(ori_ord);
								n_repeat=0;

								if (0) for (unsigned int k=0;k<memo.size();k++) cout << memo.at(k) << endl;
							}
						}
					}
					else n_repeat++;

				} while ( next_permutation(nonH_ord.begin(),nonH_ord.end()) ); //ord.begin(),ord.end()
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
						random_shuffle(nonH_ord.begin(),nonH_ord.end(),myrandom);
						//shuffle(nonH_ord.begin(),nonH_ord.end(),default_random_engine(123)); //time(0)
						//for (unsigned int j=0;j<chiralH.size();j+=2) {
						//	for (unsigned int k=0;k<nonH_ord.size();k++) {
						//		if (chiralH.at(j)==nonH_ord.at(k)) {
						//			if (k<=nonH_ord.size()-2) {
						//				nonH_ord.insert(nonH_ord.begin()+(k+1),chiralH.at(j+1));
						//			}
						//			else {
						//				nonH_ord.push_back(chiralH.at(j+1));
						//			}
						//			break;
						//		}
						//	}
						//}
						ord.resize(0);
						for (unsigned int j=0;j<nonH_ord.size();j++) ord.push_back(nonH_ord.at(j));
						for (unsigned int j=0;j<H_ord.size();j++) ord.push_back(H_ord.at(j));
					}
					if (0) {
						random_shuffle(ord.begin(),ord.end());
					}

					mol.RenumberAtoms(ord);
					string buf=mol2smi(mol);
					string can_SMI2=canonicalize_smi(buf);

					//if (!(i/100)) cout << "HHH " << can_SMI1 << " | " << can_SMI2 << endl;

					if (can_SMI1==can_SMI2) {
						for (unsigned int j=0;j<memo.size();j++) {
							if (memo.at(j)==buf) {
								//i--;
								//i=memo.size()-1;
								sort(nonH_ord.begin(),nonH_ord.end());
								n_repeat++;
								break;
							}
							else if (memo.at(j)!=buf && j>=memo.size()-1) {
								memo.push_back(buf);
								//stringstream o1("");
								//for (unsigned int k1=0;k1<ord.size();k1++) {
								//	if (k1<=ord.size()-2) o1 << ord.at(k1) << "-";
								//	else o1 << ord.at(k1);
								//}

								cout << setw(digits+1) << (ct+1) << " " << setw(80) << buf << endl; //<< " " << right << setw(80) << o1.str() << endl;
								ouf1 << setw(digits+1) << (ct+1) << " " << setw(80) << buf << endl; //" " << right << setw(80) << o1.str() << endl;
								ct++;
								//mol.RenumberAtoms(ori_ord);
								n_repeat=0;

								if (0) for (unsigned int k=0;k<memo.size();k++) cout << memo.at(k) << endl;
							}
						}
					}
					else {
						sort(nonH_ord.begin(),nonH_ord.end());
						//next_permutation(nonH_ord.begin(),nonH_ord.end());
						n_repeat++;
					}

					//if (n_repeat>=pow(nonH_ord.size(),2)) break;
					//if (nonH_ord.size()>1) {
					//	//if (n_repeat>=nonH_ord.size()*(ceil(log(nonH_ord.size()))-1)) break;
					//	if (n_repeat>=nonH_ord.size()*(nonH_ord.size()-1)*(nonH_ord.size()-2)) break;
					//}
					//else break;
					if (n_repeat>=TOLERANCE) break;
				}
			}

			ouf1.close();
			mol.Clear();
		}

		memo.resize(0);

		if (IF_RDKIT) {
			std::cout << SMI << " -----------RDKIT-----------" << std::endl;

			RDKit::ROMol *mmol = RDKit::SmilesToMol(SMI);
			mmol = RDKit::MolOps::addHs( *mmol );
			RDKit::DGeomHelpers::EmbedMolecule( *mmol );
			//RDKit::MMFF::MMFFOptimizeMolecule( *mmol , 1000 , "MMFF94s" );
			RDKit::UFF::UFFOptimizeMolecule( *mmol , 1000 );

			mmol = RDKit::MolOps::removeHs( *mmol );
			RDKit::DGeomHelpers::EmbedMolecule( *mmol );

			string can_SMIR=RDKit::MolToSmiles( *mmol );
			//cout << "can " <<  can_SMIR << endl;

			//std::cout << "Number of atoms : " << mmol->getNumAtoms() << std::endl;
			//std::cout << "Original order: " << RDKit::MolToMolBlock( *mmol, true, -1, true, true ) << std::endl;
			//std::cout << "Original SMI: " << RDKit::MolToSmiles( *mmol, true, false, -1, false, false, false, false ) << std::endl;
			//std::cout << "Original order: " << RDKit::MolToMolBlock( *mmol ) << std::endl;
			//std::cout << "Original SMI: " << RDKit::MolToSmiles( *mmol ) << std::endl;


			if (1) {
				memo.push_back(SMI);
				std::vector<unsigned int> ord1(mmol->getNumAtoms(),0);
				for (unsigned int i=0;i<ord1.size();i++) ord1.at(i)=i;
				long long unsigned int N_arrange=tgamma(mmol->getNumAtoms()+1),ct=0,n_repeat=0;
				unsigned int digits = 0;
				if (mmol->getNumAtoms()>=2) digits = log10(N_arrange)+1;

				if (!TOLERANCE_set) {
					if (ord1.size()>2) TOLERANCE=ord1.size()*(ord1.size()-1)*(ord1.size()-2);
    				else TOLERANCE=2;
				}

	            string a=can_SMIR;
    	        for (unsigned int i=0;i<a.length();i++) {
        	        if (a[i]=='/') a[i]='u';
            	    else if (a[i]=='\\') a[i]='d';
                	else if (a[i]=='*') a[i]='x';
            	}

            	ofstream ouf2((SYN_DIR+"RDKit_"+a+".txt").c_str());
            	ouf2 << "Original SMILES: " << SMI << endl;
				ouf2 << "TOLERANCE: " << TOLERANCE << endl;
            	ouf2 << "Canonical SMILES: " << can_SMIR << endl;

				for (long long unsigned int i=1;i<=N_arrange;i++) {
					random_shuffle(ord1.begin(),ord1.end());

					RDKit::ROMol *mmol1=RDKit::MolOps::renumberAtoms( *mmol , ord1 );
					//RDKit::DGeomHelpers::EmbedMolecule( *mmol1 );
					//RDKit::UFF::UFFOptimizeMolecule( *mmol1 , 1000 );
					string buf=RDKit::MolToSmiles( *mmol1, true, false, -1, false, false, false, false );
					string can_SMIR1=RDKit::MolToSmiles( *mmol1 );

					//if (!(i/100)) cout << "GGG " << can_SMIR << " | " << can_SMIR1 << endl;

					if (can_SMIR==can_SMIR1) {
						for (unsigned int j=0;j<memo.size();j++) {
							if (memo.at(j)==buf) {
								//i--;
								//i=memo.size()-1;
								n_repeat++;
								//sort(ord1.begin(),ord1.end());
								break;
							}
							else if (memo.at(j)!=buf && j>=memo.size()-1) {
								memo.push_back(buf);
								cout << setw(digits+1) << (ct+1) << " " << setw(80) << buf << endl;
								ouf2 << setw(digits+1) << (ct+1) << " " << setw(80) << buf << endl;
								ct++;
								//mol.RenumberAtoms(ori_ord);
								n_repeat=0;

								if (0) for (unsigned int k=0;k<memo.size();k++) cout << memo.at(k) << endl;
							}
						}
					}
					else {
						n_repeat++;
						sort(ord1.begin(),ord1.end());
					}

					mmol1=NULL;

					//if (n_repeat>=pow(mmol->getNumAtoms(),3)) break;
					if (n_repeat>=TOLERANCE) break;
				}

				ouf2.close();
			}

			mmol = NULL;

		}

	}
	return 1;
}

unsigned int smi2mol(string SMI,OpenBabel::OBMol &mol) {
	stringstream smi(SMI);
	OpenBabel::OBConversion conv(&smi);
	conv.SetInFormat("SMI");
	//conv.AddOption("S",OpenBabel::OBConversion::INOPTIONS);
	conv.Read(&mol);
	OpenBabel::PerceiveStereo(&mol);
	mol.AddHydrogens();
	OpenBabel::PerceiveStereo(&mol);

	return 1;
}

string mol2smi(OpenBabel::OBMol &mol) {
	stringstream ss1("");
	stringstream ss("");
	stringstream buf("");
	for (unsigned int i=0;i<mol.NumAtoms();i++) { //mol.NumHvyAtoms()
		if (i<=mol.NumAtoms()-2) buf << (i+1) << "-"; //mol.NumHvyAtoms()
		else buf << (i+1);
	}

	if (1) {
        OpenBabel::OBBuilder builder;
        builder.Build(mol);

		OpenBabel::OBConversion conv;
		//conv.SetInFormat("MOL");
		//conv.Read(&mol);
		conv.SetOutFormat("MOL");
		conv.AddOption("3", OpenBabel::OBConversion::OUTOPTIONS);
		//conv.AddOption("v", OpenBabel::OBConversion::OUTOPTIONS);
		//conv.AddOption("S",OpenBabel::OBConversion::INOPTIONS);
		//conv.AddOption("S",OpenBabel::OBConversion::OUTOPTIONS);
		//conv.AddOption("s",OpenBabel::OBConversion::INOPTIONS);
		//conv.AddOption("minimize", OpenBabel::OBConversion::GENOPTIONS);
		//conv.AddOption("ff", OpenBabel::OBConversion::GENOPTIONS, "uff");
		//conv.AddOption("steps", OpenBabel::OBConversion::GENOPTIONS, "1000");
		//conv.AddOption("c", OpenBabel::OBConversion::OUTOPTIONS);
		OpenBabel::PerceiveStereo(&mol);
		conv.Write(&mol,&ss1);
		//cout <<"ss1: " << endl
		//	 << ss1.str() << endl;
		//unsigned df=0;
		//OpenBabel::OBElementTable ATMSYM;
		//FOR_ATOMS_OF_MOL(atom, mol) {
		//	cout << df << " " << ATMSYM.GetSymbol(atom->GetAtomicNum()) << " | AnC: " << atom->IsAntiClockwise() << " | C: " << atom->IsClockwise() << endl;
		//}
	}
	if (1) {
		OpenBabel::OBConversion conv(&ss1,&ss);
		if (conv.SetInAndOutFormats("MOL","SMI")) {
			conv.AddOption("o",OpenBabel::OBConversion::OUTOPTIONS,buf.str().c_str());
			//conv.AddOption("s",OpenBabel::OBConversion::INOPTIONS);
			//conv.AddOption("U", OpenBabel::OBConversion::OUTOPTIONS);
			//conv.AddOption("S",OpenBabel::OBConversion::INOPTIONS);
			//conv.AddOption("h", OpenBabel::OBConversion::OUTOPTIONS);
			//conv.AddOption("U", OpenBabel::OBConversion::OUTOPTIONS);
			//conv.AddOption("C", OpenBabel::OBConversion::OUTOPTIONS);
			//conv.AddOption("f", OpenBabel::OBConversion::OUTOPTIONS, "1");
			//conv.AddOption("gen3D", OpenBabel::OBConversion::GENOPTIONS);
			//conv.AddOption("3", OpenBabel::OBConversion::OUTOPTIONS);
			conv.Convert();
			//cout <<"ss: " << ss.str() << endl;

		}
	}

	string null="";
	ss >> null >> ws;

	return null;
}

string canonicalize_smi(string smi) {
	stringstream ss_out("");
	stringstream ss_out1("");

	if (1) {
		stringstream ss(smi);
		OpenBabel::OBConversion conv(&ss, &ss_out);
		if(conv.SetInAndOutFormats("SMI","SMI")) {
			conv.AddOption("h", OpenBabel::OBConversion::OUTOPTIONS);
			conv.AddOption("S",OpenBabel::OBConversion::INOPTIONS);
			//conv.AddOption("U", OpenBabel::OBConversion::OUTOPTIONS);
			//conv.AddOption("S", OpenBabel::OBConversion::INOPTIONS);
			//conv.AddOption("canonical", OpenBabel::OBConversion::GENOPTIONS);
			//conv.AddOption("C", OpenBabel::OBConversion::OUTOPTIONS);
			//conv.AddOption("l", OpenBabel::OBConversion::OUTOPTIONS,"1");
			//conv.AddOption("gen2D", OpenBabel::OBConversion::GENOPTIONS);
			//conv.AddOption("minimize", OpenBabel::OBConversion::GENOPTIONS);
			//conv.AddOption("ff", OpenBabel::OBConversion::GENOPTIONS, "uff");
			//conv.AddOption("steps", OpenBabel::OBConversion::GENOPTIONS, "1000");
			conv.Convert();
		}
	}
	if (1) {
		OpenBabel::OBConversion conv(&ss_out, &ss_out1);
		if(conv.SetInAndOutFormats("SMI","SMI")) { //"SMI"
			conv.AddOption("S", OpenBabel::OBConversion::INOPTIONS);
			//conv.AddOption("U", OpenBabel::OBConversion::OUTOPTIONS);
			conv.AddOption("canonical", OpenBabel::OBConversion::GENOPTIONS);
			//conv.AddOption("C", OpenBabel::OBConversion::OUTOPTIONS);
			//conv.AddOption("f", OpenBabel::OBConversion::OUTOPTIONS, "1");
			//conv.AddOption("gen2D", OpenBabel::OBConversion::GENOPTIONS);
			//conv.AddOption("minimize", OpenBabel::OBConversion::GENOPTIONS);
			//conv.AddOption("ff", OpenBabel::OBConversion::GENOPTIONS, "uff");
			//conv.AddOption("steps", OpenBabel::OBConversion::GENOPTIONS, "1000");
			conv.Convert();
		}
	}
	string null="";
	ss_out1 >> null >> ws;

	return null;
}

int myrandom(int i) { 
	return rand()%i;
}

unsigned int rd_para(char **&argv,bool &IF_OBABEL, bool &IF_RDKIT, string &SYN_DIR, string &SMI_LIST, unsigned long long int &TOLERANCE) {
	ifstream inf(argv[1]);
	inf >> ws;
	if (!inf.is_open() || inf.eof()) {
		inf.close();
		return 0;
	}

    while (!inf.eof() && inf.is_open()) {
        unsigned int cur=inf.tellg();
        string null="#";
        inf >> null >> ws;
		if (null=="IF_OBABEL") {
			inf >> null >> ws;
			IF_OBABEL=(bool)atoi(null.c_str());
		}
		else if (null=="IF_RDKIT") {
			inf >> null >> ws;
			IF_RDKIT=(bool)atoi(null.c_str());
		}
		else if (null=="SYN_DIR") {
			inf >> null >> ws;
			SYN_DIR=null;
		}
        else if (null=="SMI_LIST") {
            inf >> null >> ws;
            SMI_LIST=null;
        }
        else if (null=="TOLERANCE") {
            inf >> null >> ws;
            TOLERANCE=atoi(null.c_str());;
        }
        else if (null[0]=='#') {
            inf.seekg(cur);
            getline(inf,null);
            inf >> ws;
        }
	}	

	return 1;
}

unsigned int print_para(bool IF_OBABEL, bool IF_RDKIT, string SYN_DIR, string SMI_LIST, unsigned long long int TOLERANCE) {
    cout << "=====================================================================" << endl;
    cout << left << setw(15) << "SYN_DIR" << " " << left << setw(15) << SYN_DIR << endl
		<< left << setw(15) << "SMI_LIST" << " " << left << setw(15) << SMI_LIST << endl
        << left << setw(15) << "IF_OBABEL" << " " << left << setw(15) << IF_OBABEL << endl
        << left << setw(15) << "IF_RDKIT" << " " << left << setw(15) << IF_RDKIT << endl;
	if (TOLERANCE==1e+15) {
    	cout << left << setw(15) << "TOLERANCE" << " " << left << setw(15) << "nonH_ord.size()*(nonH_ord.size()-1)*(nonH_ord.size()-2)" << endl;
	}
	else {
		cout << left << setw(15) << "TOLERANCE" << " " << left << setw(15) << TOLERANCE << endl;
	}
	cout << "=====================================================================" << endl;
	return 1;
}
