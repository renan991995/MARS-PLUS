#include "ELEMENTS.h"
using namespace std;

extern PARAMETER para;

unsigned int POOL::set_up() {
	num=71;
	unsigned int i,j;//,k;
	
	a.resize(num);
	for (i=0;i<num;i++) {
		a.at(i).order.resize(6);
	}
	for (i=0;i<num;i++) {
		a.at(i).id=i;
		for (j=0;j<6;j++) a.at(i).order.at(j)=0;
		for (j=0;j<3;j++) a.at(i).bd[j]=0;
		a.at(i).index=2; // for neutral atom
		a.at(i).chg=0;	// for neutral atom
		a.at(i).probability=1;
	}
	
	a.at(0).name="[H](-)";
	a.at(0).nbond=6;
	a.at(0).order.at(0)=1;
	a.at(0).bd[0]=1;
	a.at(0).norder=1;
	a.at(0).index=4;
	a.at(0).atm="H";
	a.at(0).probability=0;

	a.at(1).name="C(-)(-)(-)(-)";
	a.at(1).order.at(0)=a.at(1).order.at(1)=a.at(1).order.at(2)=a.at(1).order.at(3)=1;
	a.at(1).nbond=13;
	a.at(1).norder=4;
	a.at(1).bd[0]=4;
	a.at(1).atm="C";

	a.at(2).name="C(=)(-)(-)";
	a.at(2).order.at(0)=2;
	a.at(2).order.at(1)=a.at(2).order.at(2)=1;
	a.at(2).nbond=10;
	a.at(2).norder=3;
	a.at(2).bd[0]=2;
	a.at(2).bd[1]=1;
	a.at(2).atm="C";

	a.at(3).name="C(#)(-)";
	a.at(3).order.at(0)=3;
	a.at(3).order.at(1)=1;
	a.at(3).nbond=7;
	a.at(3).norder=2;
	a.at(3).bd[0]=1;
	a.at(3).bd[2]=1;
	a.at(3).atm="C";

	a.at(4).name="C(=)(=)";
	a.at(4).order.at(0)=a.at(4).order.at(1)=2;
	a.at(4).nbond=7;
	a.at(4).norder=2;
	a.at(4).bd[1]=2;
	a.at(4).atm="C";

	a.at(5).name="O(-)(-)";
	a.at(5).order.at(0)=a.at(5).order.at(1)=1;
	a.at(5).nbond=7;
	a.at(5).norder=2;
	a.at(5).bd[0]=2;
	a.at(5).atm="O";

	a.at(6).name="O(=)";
	a.at(6).order.at(0)=2;
	a.at(6).nbond=4;
	a.at(6).norder=1;
	a.at(6).bd[1]=1;
	a.at(6).atm="O";

	a.at(7).name="N(-)(-)(-)";
	a.at(7).order.at(0)=a.at(7).order.at(1)=a.at(7).order.at(2)=1;
	a.at(7).nbond=10;
	a.at(7).norder=3;
	a.at(7).bd[0]=3;
	a.at(7).atm="N";

	a.at(8).name="N(=)(-)";
	a.at(8).order.at(0)=2;
	a.at(8).order.at(1)=1;
	a.at(8).nbond=7;
	a.at(8).norder=2;
	a.at(8).bd[0]=1;
	a.at(8).bd[1]=1;
	a.at(8).atm="N";

	a.at(9).name="N(#)";
	a.at(9).order.at(0)=3;
	a.at(9).nbond=4;
	a.at(9).norder=1;
	a.at(9).bd[2]=1;
	a.at(9).atm="N";

	a.at(10).name="O(-)";
	a.at(10).order.at(0) = 1;
	a.at(10).nbond=4;
	a.at(10).norder=1;
	a.at(10).bd[0]=1;
	a.at(10).atm="O";

	a.at(11).name = "F(-)";
	a.at(11).order.at(0) = 1;
	a.at(11).nbond=4;
	a.at(11).norder=1;
	a.at(11).bd[0]=1;
	a.at(11).atm="F";

	a.at(12).name="Cl(-)";
	a.at(12).order.at(0) = 1;
	a.at(12).nbond=5;
	a.at(12).norder=1;
	a.at(12).index=3;
	a.at(12).bd[0]=1;
	a.at(12).atm="Cl";	

    a.at(13).name="Br(-)";
    a.at(13).order.at(0) = 1;
    a.at(13).nbond = 5;
    a.at(13).norder=1;
	a.at(13).bd[0]=1;
	a.at(13).index=3;
	a.at(13).atm="Br";

	a.at(14).name ="I(-)";
	a.at(14).order.at(0)=1;
	a.at(14).nbond = 4;
	a.at(14).norder=1;
	a.at(14).bd[0]=1;
	a.at(14).atm="I";
	
	a.at(15).name="[NH0+](-)(-)(-)(-)";
	for (i=0;i<4;i++) a.at(15).order.at(i)=1;
	a.at(15).nbond=18;
	a.at(15).norder=4;
	a.at(15).bd[0]=4;
	a.at(15).chg=1;
	a.at(15).index=7;
	a.at(15).atm="N";

	a.at(16).name="[NH0+](=)(-)(-)";
	a.at(16).order.at(0)=2;
	a.at(16).order.at(1)=a.at(16).order.at(2)=1;
	a.at(16).nbond=15;
	a.at(16).norder=3;
	a.at(16).bd[0]=2;
	a.at(16).bd[1]=1;
	a.at(16).chg=1;
	a.at(16).index=7;
	a.at(16).atm="N";

	a.at(17).name="[PH0+](-)(-)(-)(-)";
	for (i=0;i<4;i++) a.at(17).order.at(i)=1;
	a.at(17).norder=4;
	a.at(17).nbond=18;
	a.at(17).bd[0]=4;
	a.at(17).chg=1;
	a.at(17).index=7;
	a.at(17).atm="P";

	a.at(18).name="[PH0+](=)(-)(-)";
	a.at(18).order.at(0)=2;
	a.at(18).order.at(1)=a.at(18).order.at(2)=1;
	a.at(18).norder=3;
	a.at(18).nbond=15;
	a.at(18).bd[0]=2;
	a.at(18).bd[1]=1;
	a.at(18).chg=1;
	a.at(18).index=7;
	a.at(18).atm="P";
	a.at(18).probability=0; 

	a.at(19).name="S(-)(-)";
	a.at(19).order.at(0)=a.at(19).order.at(1)=1;
	a.at(19).norder=2;
	a.at(19).nbond=7;
	a.at(19).bd[0]=2;
	a.at(19).atm="S";

	a.at(20).name="S(=)";
	a.at(20).order.at(0)=2;
	a.at(20).nbond=4;
	a.at(20).norder=1;
	a.at(20).bd[1]=1;
	a.at(20).atm="S";

	a.at(21).name="P(-)(-)(-)";
	a.at(21).order.at(0)=a.at(21).order.at(1)=a.at(21).order.at(2)=1;
	a.at(21).nbond=10;
	a.at(21).norder=3;
	a.at(21).bd[0]=3;
	a.at(21).atm="P";

	a.at(22).name="P(=)(-)";
	a.at(22).order.at(0)=2;
	a.at(22).order.at(1)=1;
	a.at(22).nbond=7;
	a.at(22).norder=2;
	a.at(22).bd[0]=1;
	a.at(22).bd[1]=1;
	a.at(22).atm="P";
	
	a.at(23).name="P(#)";
	a.at(23).order.at(0)=3;
	a.at(23).nbond=4;
	a.at(23).norder=1;
	a.at(23).bd[2]=1;
	a.at(23).atm="P";

	a.at(24).name="[F-]";
	a.at(24).nbond=4;
	a.at(24).norder=0;
	a.at(24).chg=-1;
	a.at(24).atm="F";
	a.at(24).index=4; 

	a.at(25).name="[Cl-]";
	a.at(25).nbond=5;
	a.at(25).norder=0;
	a.at(25).chg=-1;
	a.at(25).atm="Cl";
	a.at(25).index=5; 

	a.at(26).name="[Br-]";
	a.at(26).nbond=5;
	a.at(26).norder=0;
	a.at(26).chg=-1;
	a.at(26).atm="Br";
	a.at(26).index=5; 

	a.at(27).name="[I-]";
	a.at(27).nbond=4;
	a.at(27).norder=0;
	a.at(27).chg=-1;
	a.at(27).atm="I";
	a.at(27).index=4; 

	a.at(28).name="[OH1-]";
	a.at(28).nbond=6;
	a.at(28).norder=0;
	a.at(28).chg=-1;
	a.at(28).atm="O";
	a.at(28).index=6; 

	a.at(29).name="[OH0-](-)";
	a.at(29).order.at(0)=1;
	a.at(29).nbond=9;
	a.at(29).norder=1;
	a.at(29).chg=-1;
	a.at(29).bd[0]=1;
	a.at(29).index=7;	
	a.at(29).atm="O";

	a.at(30).name="[PH0-](-)(-)(-)(-)(-)(-)";
	for (i=0;i<6;i++) a.at(30).order.at(i)=1;
	a.at(30).index=7;
	a.at(30).nbond=24;
	a.at(30).norder=6;
	a.at(30).chg=-1;
	a.at(30).bd[0]=6;
	a.at(30).atm="P";

	a.at(31).name="P(-)(-)(-)(-)(-)";
	for (i=0;i<5;i++) a.at(31).order.at(i)=1;
	a.at(31).nbond=16;
	a.at(31).norder=5;
	a.at(31).bd[0]=5;
	a.at(31).atm="P";

	a.at(32).name="P(=)(-)(-)(-)";
	a.at(32).order.at(0)=2;
	for (i=1;i<4;i++) a.at(32).order.at(i)=1;
	a.at(32).nbond=13;
	a.at(32).norder=4;
	a.at(32).bd[0]=3;
	a.at(32).bd[1]=1;
	a.at(32).atm="P";
	a.at(32).probability=0;

	a.at(33).name="S(-)(=O)(=O)([O-])";
	a.at(33).order.at(0)=1;
	a.at(33).nbond=4;
	a.at(33).norder=1;
	a.at(33).chg=-1;
	a.at(33).bd[0]=1;
	a.at(33).atm="S";

	a.at(34).name="S(=)(-)(-)";
	a.at(34).order.at(0)=2;
	a.at(34).order.at(1)=a.at(34).order.at(2)=1;
	a.at(34).norder=3;
	a.at(34).nbond=10;
	a.at(34).bd[0]=2;
	a.at(34).bd[1]=1;
	a.at(34).atm="S";
	a.at(34).probability=0;

	a.at(35).name="[NH0-](-)(-)";
	a.at(35).order.at(0)=a.at(35).order.at(1)=1;
	a.at(35).nbond=12;
	a.at(35).norder=2;
	a.at(35).chg=-1;
	a.at(35).bd[0]=2;
	a.at(35).index=7;
	a.at(35).atm="N";	

	a.at(36).name="C(-)(-)(-)([N+]%99999C=CN(C)C=%99999)";
    a.at(36).order.at(0)=a.at(36).order.at(1)=a.at(36).order.at(2)=1;
    a.at(36).nbond=10;
    a.at(36).norder=3;
    a.at(36).chg=1;
    a.at(36).bd[0]=3;
    a.at(36).atm="C";

    a.at(37).name="C(-)(-)(-)([N+]%99999C=CN(C)C(C)=%99999)";
    a.at(37).order.at(0)=a.at(37).order.at(1)=a.at(37).order.at(2)=1;
    a.at(37).nbond=10;
    a.at(37).norder=3;
    a.at(37).chg=1;
    a.at(37).bd[0]=3;
    a.at(37).atm="C";

    a.at(38).name="C(-)(-)(-)(N%99999C=C[N+](C)=C%99999)";
    a.at(38).order.at(0)=a.at(38).order.at(1)=a.at(38).order.at(2)=1;
    a.at(38).nbond=10;
    a.at(38).norder=3;
    a.at(38).chg=1;
    a.at(38).bd[0]=3;
    a.at(38).atm="C";

    a.at(39).name="C(-)(-)(-)([N+]%99999=CC=CC(C)=C%99999)";
    a.at(39).order.at(0)=a.at(39).order.at(1)=a.at(39).order.at(2)=1;
    a.at(39).nbond=10;
    a.at(39).norder=3;
    a.at(39).chg=1;
    a.at(39).bd[0]=3;
    a.at(39).atm="C";

    a.at(40).name="C(-)(-)(-)(C%99999=C[N+](C)=CC=C%99999)";
    a.at(40).order.at(0)=a.at(40).order.at(1)=a.at(40).order.at(2)=1;
    a.at(40).nbond=10;
    a.at(40).norder=3;
    a.at(40).chg=1;
    a.at(40).bd[0]=3;
    a.at(40).atm="C";

    a.at(41).name="C(-)(%99999=[NH+]C=CC=C%99999)";
    a.at(41).order.at(0)=1;
    a.at(41).nbond=4;
    a.at(41).norder=1;
    a.at(41).chg=1;
    a.at(41).bd[0]=1;
    a.at(41).atm="C";

    a.at(42).name="[NH0+](-)(%99999=CC=CC=C%99999)";
    a.at(42).order.at(0)=1;
    a.at(42).nbond=9;
    a.at(42).norder=1;
    a.at(42).chg=1;
    a.at(42).bd[0]=1;
    a.at(42).index=7;
    a.at(42).atm="N";

    a.at(43).name="[NH0-](S(=O)(=O)C(F)(F)(F))(S(=O)(=O)C(F)(F)(F))";
    a.at(43).nbond=48; 
    a.at(43).norder=0;
    a.at(43).chg=-1;
    a.at(43).index=48;  
    a.at(43).atm="N";

	a.at(44).name="[BH0-](-)(-)(-)(-)";
    a.at(44).order.at(0)=a.at(44).order.at(1)=a.at(44).order.at(2)=a.at(44).order.at(3)=1;
    a.at(44).nbond=18;
    a.at(44).norder=4;
    a.at(44).chg=-1;
    a.at(44).bd[0]=4;
    a.at(44).index=7;
    a.at(44).atm="B";

	a.at(45).name="C(-)(-)(-)(C(=O)([O-]))";
    a.at(45).order.at(0)=a.at(45).order.at(1)=a.at(45).order.at(2)=1;
    a.at(45).nbond=10;
    a.at(45).norder=3;
    a.at(45).chg=-1;
    a.at(45).bd[0]=3;
    a.at(45).atm="C";

	a.at(46).name="C(#N)([S-])";
    a.at(46).nbond=11;
    a.at(46).norder=0;
    a.at(46).chg=-1;
    a.at(46).index=11;
    a.at(46).atm="C";

    a.at(47).name="C(-)(-)(-)(OP(=O)(OC)([O-]))";
    a.at(47).order.at(0)=a.at(47).order.at(1)=a.at(47).order.at(2)=1;
    a.at(47).nbond=10;
    a.at(47).norder=3;
    a.at(47).chg=-1;
    a.at(47).bd[0]=3;
    a.at(47).atm="C";

    a.at(48).name="C(#N)([N-]C#N)";
    a.at(48).nbond=14; 
    a.at(48).norder=0;
    a.at(48).chg=-1;
    a.at(48).index=14; 
    a.at(48).atm="C";

    a.at(49).name="[BH0-](C#N)(C#N)(C#N)(C#N)"; 
    a.at(49).nbond=26; 
    a.at(49).norder=0;
    a.at(49).chg=-1;
    a.at(49).index=26; 
    a.at(49).atm="B";

    a.at(50).name="S(OCCOCCOC)(=O)(=O)([O-])"; 
    a.at(50).nbond=25; 
    a.at(50).norder=0;
    a.at(50).chg=-1;
    a.at(50).index=25;
    a.at(50).atm="S";

    a.at(51).name="S(c(cc%99999)ccc%99999C)(=O)(=O)([O-])"; 
    a.at(51).nbond=28; 
    a.at(51).norder=0;
    a.at(51).chg=-1;
    a.at(51).index=28; 
    a.at(51).atm="S";

    a.at(52).name="[PH0-](F)(F)(F)(C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)"; 
    a.at(52).nbond=66; 
    a.at(52).norder=0;
    a.at(52).chg=-1;
    a.at(52).atm="P";
    a.at(52).index=66; 

    a.at(53).name="[NH0-](S(=O)(=O)C(C(F)(F)F)(F)F)(S(=O)(=O)C(C(F)(F)F)(F)F)"; 
    a.at(53).nbond=58; 
    a.at(53).norder=0;
    a.at(53).chg=-1;
    a.at(53).atm="N";
    a.at(53).index=58; 

    a.at(54).name="[CH0-](S(C(F)(F)(F))(=O)(=O))(S(C(F)(F)(F))(=O)(=O))(S(C(F)(F)(F))(=O)(=O))"; 
    a.at(54).nbond=75; 
    a.at(54).norder=0;
    a.at(54).chg=-1;
    a.at(54).atm="C";
    a.at(54).index=75; 

    a.at(55).name="[PH0-](F)(F)(F)(F)(F)(F)"; 
    a.at(55).nbond=24; 
    a.at(55).norder=0;
    a.at(55).chg=-1;
    a.at(55).atm="C";
    a.at(55).index=24; 

    a.at(56).name="[In+3]([Cl-])([Cl-])([Cl-])([Cl-])"; 
    a.at(56).nbond=34; 
    a.at(56).norder=0;
    a.at(56).chg=-1;
    a.at(56).atm="In";
    a.at(56).index=34; 
    a.at(56).probability=0;

    a.at(57).name="[In+3](-)(-)(-)(-)";
    a.at(57).order.at(0)=a.at(57).order.at(1)=a.at(57).order.at(2)=a.at(57).order.at(3)=1;
    a.at(57).nbond=18;
    a.at(57).norder=4;
    a.at(57).chg=3;
    a.at(57).bd[0]=4;
    a.at(57).index=7;
    a.at(57).atm="In";
    a.at(57).probability=0;

    a.at(58).name="Cl(=O)(=O)([O-])(=O)"; 
    a.at(58).nbond=20; 
    a.at(58).norder=0;
    a.at(58).chg=-1;
    a.at(58).atm="Cl";
    a.at(58).index=20; 

    a.at(59).name="[CH0-](-)(-)(-)"; 
	a.at(59).order.at(0)=a.at(59).order.at(1)=a.at(59).order.at(2)=1;
    a.at(59).nbond=15;
    a.at(59).norder=3;
    a.at(59).chg=-1;
	a.at(59).bd[0]=3;
    a.at(59).atm="C";
    a.at(59).index=7;

    a.at(60).name="[NH0+](=O)([O-])([O-])"; 
    a.at(60).nbond=22;
    a.at(60).norder=0;
    a.at(60).chg=-1;
    a.at(60).atm="N";
    a.at(60).index=22; 

    a.at(61).name="S(=)(=)(-)(-)";
    a.at(61).order.at(0)=a.at(61).order.at(1)=2;
	a.at(61).order.at(2)=a.at(61).order.at(3)=1;
    a.at(61).nbond=13;
    a.at(61).norder=4;
    a.at(61).bd[0]=a.at(61).bd[1]=2;
    a.at(61).atm="S";
    a.at(61).probability=0;

    a.at(62).name="Cl(=)(=)(=)(-)"; 
    a.at(62).order.at(0)=a.at(62).order.at(1)=a.at(62).order.at(2)=2;
	a.at(62).order.at(3)=1;
    a.at(62).nbond=14;
    a.at(62).index=3;
    a.at(62).norder=4;
	a.at(62).bd[0]=1;
    a.at(62).bd[1]=3;
    a.at(62).atm="Cl";
    a.at(62).index=3;
    a.at(62).probability=0;

    a.at(63).name="[SH0-](-)";
	a.at(63).order.at(0)=1;
    a.at(63).nbond=9;
    a.at(63).norder=1;
    a.at(63).chg=-1;
	a.at(63).bd[0]=1;
    a.at(63).index=7;
    a.at(63).atm="S";

    a.at(64).id=64;
    a.at(64).name="[Ga+3](-)(-)(-)(-)";
    a.at(64).order.at(0)=a.at(64).order.at(1)=a.at(64).order.at(2)=a.at(64).order.at(3)=1;
    a.at(64).nbond=18;
    a.at(64).norder=4;
    a.at(64).chg=3;
    a.at(64).bd[0]=4;
    a.at(64).index=7;
    a.at(64).atm="Ga";
    a.at(64).probability=0;

    a.at(65).name="[SH0+](-)(-)(-)";
    a.at(65).order.at(0)=a.at(65).order.at(1)=a.at(65).order.at(2)=1;
    a.at(65).nbond=15;
    a.at(65).norder=3;
    a.at(65).chg=1;
    a.at(65).bd[0]=3;
    a.at(65).index=7;
    a.at(65).atm="S";
    a.at(65).probability=0;

    a.at(66).name="P(=)(-)(-)";
    a.at(66).order.at(0)=2;
    for (i=1;i<3;i++) a.at(66).order.at(i)=1;
    a.at(66).nbond=10;
    a.at(66).norder=3;
    a.at(66).bd[0]=2;
    a.at(66).bd[1]=1;
    a.at(66).atm="P";
    a.at(66).probability=0;

    a.at(67).name="[CH0@@](-)(-)(-)(-)";
    a.at(67).order.at(0)=a.at(67).order.at(1)=a.at(67).order.at(2)=a.at(67).order.at(3)=1;
    a.at(67).nbond=19;
    a.at(67).norder=4;
    a.at(67).bd[0]=4;
    a.at(67).atm="C";
    a.at(67).probability=0;
    a.at(67).index=8;

    a.at(68).name="[CH0@](-)(-)(-)(-)";
    a.at(68).order.at(0)=a.at(68).order.at(1)=a.at(68).order.at(2)=a.at(68).order.at(3)=1;
    a.at(68).nbond=18;
    a.at(68).norder=4;
    a.at(68).bd[0]=4;
    a.at(68).atm="C";
    a.at(68).probability=0;
    a.at(68).index=7;

    a.at(69).name="P(-)(-)(-)(-)";
    for (i=0;i<4;i++) a.at(69).order.at(i)=1;
    a.at(69).nbond=13;
    a.at(69).norder=4;
    a.at(69).bd[0]=4;
    a.at(69).atm="P";

    a.at(70).name="*(-)"; //lone pair or pseudo-atom
    a.at(70).nbond=4;
    a.at(70).order.at(0)=1;
    a.at(70).bd[0]=1;
    a.at(70).norder=1;
    a.at(70).atm="Xx";
    a.at(70).probability=0;

	/*
    a.at(71).name="[NH0+@@](-)(-)(-)(-)";
    for (i=0;i<4;i++) a.at(71).order.at(i)=1;
    a.at(71).nbond=20;
    a.at(71).norder=4;
    a.at(71).bd[0]=4;
    a.at(71).chg=1;
    a.at(71).index=9;
    a.at(71).atm="N";
	a.at(71).probability=0;

    a.at(72).name="[NH0+@](-)(-)(-)(-)";
    for (i=0;i<4;i++) a.at(72).order.at(i)=1;
    a.at(72).nbond=19;
    a.at(72).norder=4;
    a.at(72).bd[0]=4;
    a.at(72).chg=1;
    a.at(72).index=8;
    a.at(72).atm="N";
    a.at(72).probability=0;

    a.at(73).name="[PH0+@@](-)(-)(-)(-)";
    for (i=0;i<4;i++) a.at(73).order.at(i)=1;
    a.at(73).norder=4;
    a.at(73).nbond=20;
    a.at(73).bd[0]=4;
    a.at(73).chg=1;
    a.at(73).index=9;
    a.at(73).atm="P";
	a.at(73).probability=0;

    a.at(74).name="[PH0+@](-)(-)(-)(-)";
    for (i=0;i<4;i++) a.at(74).order.at(i)=1;
    a.at(74).norder=4;
    a.at(74).nbond=19;
    a.at(74).bd[0]=4;
    a.at(74).chg=1;
    a.at(74).index=8;
    a.at(74).atm="P";
    a.at(74).probability=0;

    a.at(75).name="[PH0@@](=)(-)(-)(-)";
    a.at(75).order.at(0)=2;
    for (i=1;i<4;i++) a.at(75).order.at(i)=1;
    a.at(75).nbond=19;
    a.at(75).norder=4;
    a.at(75).bd[0]=3;
    a.at(75).bd[1]=1;
	a.at(75).index=8;
    a.at(75).atm="P";
    a.at(75).probability=0;

    a.at(76).name="[PH0@](=)(-)(-)(-)";
    a.at(76).order.at(0)=2;
    for (i=1;i<4;i++) a.at(76).order.at(i)=1;
    a.at(76).nbond=18;
    a.at(76).norder=4;
    a.at(76).bd[0]=3;
    a.at(76).bd[1]=1;
    a.at(76).index=7;
    a.at(76).atm="P";
    a.at(76).probability=0;

    a.at(77).name="[CH0@@](-)(-)(-)([N+]%99999C=CN(C)C=%99999)";
    a.at(77).order.at(0)=a.at(77).order.at(1)=a.at(77).order.at(2)=1;
    a.at(77).nbond=16;
	a.at(77).index=8;
    a.at(77).norder=3;
    a.at(77).chg=1;
    a.at(77).bd[0]=3;
    a.at(77).atm="C";
	a.at(77).probability=0;

    a.at(78).name="[CH0@](-)(-)(-)([N+]%99999C=CN(C)C=%99999)";
    a.at(78).order.at(0)=a.at(78).order.at(1)=a.at(78).order.at(2)=1;
    a.at(78).nbond=15;
    a.at(78).index=7;
    a.at(78).norder=3;
    a.at(78).chg=1;
    a.at(78).bd[0]=3;
    a.at(78).atm="C";
	a.at(78).probability=0;

    a.at(79).name="[CH0@@](-)(-)(-)([N+]%99999C=CN(C)C(C)=%99999)";
    a.at(79).order.at(0)=a.at(79).order.at(1)=a.at(79).order.at(2)=1;
    a.at(79).nbond=16;
	a.at(79).index=8;
    a.at(79).norder=3;
    a.at(79).chg=1;
    a.at(79).bd[0]=3;
    a.at(79).atm="C";
	a.at(79).probability=0;

    a.at(80).name="[CH0@](-)(-)(-)([N+]%99999C=CN(C)C(C)=%99999)";
    a.at(80).order.at(0)=a.at(80).order.at(1)=a.at(80).order.at(2)=1;
    a.at(80).nbond=15;
    a.at(80).index=7;
    a.at(80).norder=3;
    a.at(80).chg=1;
    a.at(80).bd[0]=3;
    a.at(80).atm="C";
    a.at(80).probability=0;

    a.at(81).name="[CH0@@](-)(-)(-)(N%99999C=C[N+](C)=C%99999)";
    a.at(81).order.at(0)=a.at(81).order.at(1)=a.at(81).order.at(2)=1;
    a.at(81).nbond=16;
	a.at(81).index=8;
    a.at(81).norder=3;
    a.at(81).chg=1;
    a.at(81).bd[0]=3;
    a.at(81).atm="C";
	a.at(81).probability=0;

    a.at(82).name="[CH0@](-)(-)(-)(N%99999C=C[N+](C)=C%99999)";
    a.at(82).order.at(0)=a.at(82).order.at(1)=a.at(82).order.at(2)=1;
    a.at(82).nbond=15;
    a.at(82).index=7;
    a.at(82).norder=3;
    a.at(82).chg=1;
    a.at(82).bd[0]=3;
    a.at(82).atm="C";
    a.at(82).probability=0;

    a.at(83).name="[CH0@@](-)(-)(-)([N+]%99999=CC=CC(C)=C%99999)";
    a.at(83).order.at(0)=a.at(83).order.at(1)=a.at(83).order.at(2)=1;
    a.at(83).nbond=16;
	a.at(83).index=8;
    a.at(83).norder=3;
    a.at(83).chg=1;
    a.at(83).bd[0]=3;
    a.at(83).atm="C";
	a.at(83).probability=0;

    a.at(84).name="[CH0@](-)(-)(-)([N+]%99999=CC=CC(C)=C%99999)";
    a.at(84).order.at(0)=a.at(84).order.at(1)=a.at(84).order.at(2)=1;
    a.at(84).nbond=15;
    a.at(84).index=7;
    a.at(84).norder=3;
    a.at(84).chg=1;
    a.at(84).bd[0]=3;
    a.at(84).atm="C";
    a.at(84).probability=0;

    a.at(85).name="[CH0@@](-)(-)(-)(C%99999=C[N+](C)=CC=C%99999)";
    a.at(85).order.at(0)=a.at(85).order.at(1)=a.at(85).order.at(2)=1;
    a.at(85).nbond=16;
	a.at(85).index=8;
    a.at(85).norder=3;
    a.at(85).chg=1;
    a.at(85).bd[0]=3;
    a.at(85).atm="C";
	a.at(85).probability=0;

    a.at(86).name="[CH0@](-)(-)(-)(C%99999=C[N+](C)=CC=C%99999)";
    a.at(86).order.at(0)=a.at(86).order.at(1)=a.at(86).order.at(2)=1;
    a.at(86).nbond=15;
    a.at(86).index=7;
    a.at(85).norder=3;
    a.at(86).chg=1;
    a.at(86).bd[0]=3;
    a.at(86).atm="C";
    a.at(86).probability=0;

    a.at(87).name="[BH0-@@](-)(-)(-)(-)";
    a.at(87).order.at(0)=a.at(87).order.at(1)=a.at(87).order.at(2)=a.at(87).order.at(3)=1;
	a.at(87).nbond=20;
	a.at(87).index=9;
    a.at(87).norder=4;
    a.at(87).chg=-1;
    a.at(87).bd[0]=4;
    a.at(87).atm="B";
	a.at(87).probability=0;

    a.at(88).name="[BH0-@](-)(-)(-)(-)";
    a.at(88).order.at(0)=a.at(88).order.at(1)=a.at(88).order.at(2)=a.at(88).order.at(3)=1;
    a.at(88).nbond=19;
    a.at(88).index=8;
    a.at(88).norder=4;
    a.at(88).chg=-1;
    a.at(88).bd[0]=4;
    a.at(88).atm="B";
    a.at(88).probability=0;

    a.at(89).name="[CH0@@](-)(-)(-)(C(=O)([O-]))";
    a.at(89).order.at(0)=a.at(89).order.at(1)=a.at(89).order.at(2)=1;
    a.at(89).nbond=16;
	a.at(89).index=8;
    a.at(89).norder=3;
    a.at(89).chg=-1;
    a.at(89).bd[0]=3;
    a.at(89).atm="C";
	a.at(89).probability=0;

    a.at(90).name="[CH0@](-)(-)(-)(C(=O)([O-]))";
    a.at(90).order.at(0)=a.at(90).order.at(1)=a.at(90).order.at(2)=1;
    a.at(90).nbond=15;
    a.at(90).index=7;
    a.at(90).norder=3;
    a.at(90).chg=-1;
    a.at(90).bd[0]=3;
    a.at(90).atm="C";
    a.at(90).probability=0;

    a.at(91).name="[CH0@@](-)(-)(-)(OP(=O)(OC)([O-]))";
    a.at(91).order.at(0)=a.at(91).order.at(1)=a.at(91).order.at(2)=1;
    a.at(91).nbond=18;
	a.at(91).index=8;
    a.at(91).norder=3;
    a.at(91).chg=-1;
    a.at(91).bd[0]=3;
    a.at(91).atm="C";
	a.at(91).probability=0;

    a.at(92).name="[CH0@](-)(-)(-)(OP(=O)(OC)([O-]))";
    a.at(92).order.at(0)=a.at(92).order.at(1)=a.at(92).order.at(2)=1;
    a.at(92).nbond=15;
    a.at(92).index=7;
    a.at(92).norder=3;
    a.at(92).chg=-1;
    a.at(92).bd[0]=3;
    a.at(92).atm="C";
    a.at(92).probability=0;

    a.at(93).name="[PH0@@](-)(-)(-)(-)";
    for (i=0;i<4;i++) a.at(93).order.at(i)=1;
    a.at(93).nbond=19;
    a.at(93).index=8;
    a.at(93).norder=4;
    a.at(93).bd[0]=4;
    a.at(93).atm="P";
	a.at(93).probability=0;

    a.at(94).name="[PH0@](-)(-)(-)(-)";
    for (i=0;i<4;i++) a.at(94).order.at(i)=1;
    a.at(94).nbond=18;
    a.at(94).index=7;
    a.at(94).norder=4;
    a.at(94).bd[0]=4;
    a.at(94).atm="P";
    a.at(94).probability=0;
	*/

	//stringstream ss("");
    //for (i=0;i<num;i++) {
	//	cout << left << setw(4) << a.at(i).id << " "
	//		 << left << setw(4) << a.at(i).atm << " "
	//		 << left << setw(65) << a.at(i).name << " "
	//		 << left << setw(4) << a.at(i).index << " "
	//		 << left << setw(4) << a.at(i).nbond << " ";
	//	ss.str("");
    //    for (j=0;j<6;j++) {
	//		if (j<=4) ss << a.at(i).order.at(j) << " ";
	//		else ss << a.at(i).order.at(j);
	//	}
	//	cout << left << setw(15) << ss.str() << " "
	//		 << left << setw(4) << a.at(i).norder << " ";
	//	ss.str("");
    //    for (j=0;j<3;j++) {
	//		if (j<=2) ss << a.at(i).bd[j] << " ";
	//		else ss << a.at(i).bd[j];
	//	}
	//	cout << left << setw(9) << ss.str() << " "
    //    	 << left << setw(4) << a.at(i).chg << " "
    //    	 << left << setw(4) << a.at(i).probability << endl;
    //}
	//exit(0);

	return 1;
}

unsigned int POOL::read_in() {
	if (para.element_list!="none" && para.element_list!="") {
		ifstream inf(para.element_list.c_str());
		inf >> ws;
		if (!inf.is_open() || inf.eof()) {
			inf.close();
			return 0;
		}
		
		while (!inf.eof()) {
			string null="#";
			int cur=inf.tellg();
			inf >> null >> ws;
			inf.seekg(cur);
			if (null[0]!='#') {
				a.resize(a.size()+1);

				inf >> null >> ws;  a.at(a.size()-1).id=(unsigned int)atoi(null.c_str());
				inf >> null >> ws;  a.at(a.size()-1).atm=null;
				inf >> null >> ws;  a.at(a.size()-1).name=null;

                inf >> null >> ws;  a.at(a.size()-1).index=(unsigned int)atoi(null.c_str());
                inf >> null >> ws;  a.at(a.size()-1).nbond=(unsigned int)atoi(null.c_str());
                inf >> null >> ws;  a.at(a.size()-1).chg=atoi(null.c_str());
                inf >> null >> ws;  a.at(a.size()-1).probability=atof(null.c_str());
                inf >> null >> ws;  a.at(a.size()-1).chigenic=(bool)atoi(null.c_str());
                inf >> null >> ws;  a.at(a.size()-1).ctgenic=(bool)atoi(null.c_str());

				for (unsigned int i=0;i<6;i++) {
                	inf >> null >> ws;  a.at(a.size()-1).order.push_back((unsigned int)atoi(null.c_str()));
				}
				a.at(a.size()-1).norder=0;
				a.at(a.size()-1).bd[0]=a.at(a.size()-1).bd[1]=a.at(a.size()-1).bd[2]=0;
				for (unsigned int i=0;i<6;i++) {
					if (a.at(a.size()-1).order.at(i)) {
						a.at(a.size()-1).norder+=1;
						a.at(a.size()-1).bd[a.at(a.size()-1).order.at(i)-1]+=1;
					}
				}
				//inf >> null >> ws;  a.at(a.size()-1).chg=atoi(null.c_str());
				//inf >> null >> ws;  
				//a.at(a.size()-1).probability=atof(null.c_str()); 
				//inf >> null >> ws;
                //if (null=="1") a.at(a.size()-1).chigenic=1;
                //else if (null=="0") a.at(a.size()-1).chigenic=0;
                //inf >> null >> ws;
                //if (null=="1") a.at(a.size()-1).ctgenic=1;
                //else if (null=="0") a.at(a.size()-1).ctgenic=0;
			}
			else {
				getline(inf,null);
				inf >> ws;
			}
		}
		inf.close();

		num=a.size();

		return 1;
	}
	else return 0;
}
