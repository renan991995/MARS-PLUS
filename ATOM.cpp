#include "ATOM.h"
using namespace std;

void DEATOM::find_r() {
	if (name=="C") {
		r_bnd=0.77;
		nbnd=4;
	}
	else if (name=="N") {
		r_bnd=0.74;
		nbnd=3;
	}
	else if (name=="O") {
		r_bnd=0.74;
		nbnd=2;
	}
	else if (name=="F") {
		r_bnd=0.72;
		nbnd=1;
	}
	else if (name=="Cl") {
		r_bnd=0.99;
		nbnd=1;
	}
	else if (name=="Br") {
		r_bnd=1.20; 
		nbnd=1;
	}
	else if (name=="I") {
		r_bnd=1.39;
		nbnd=1;
	}
	else if (name=="H") {
		r_bnd=0.37;
		nbnd=1;
	}
	else if (name=="S") {
		r_bnd=1.05;
		nbnd=2;
	}
	else if (name=="P") {
		r_bnd=1.07;
		nbnd=5;
	}
	
	return;
}

void POOL::set_up() {
	num=70;
	int i,j,k;
	
	a.resize(num);
	for (i=0;i<num;i++) a.at(i).order.resize(6);
	for (i=0;i<num;i++) {
		a.at(i).nh=0;
		for (j=0;j<6;j++) a.at(i).order.at(j)=0;
		for (j=0;j<3;j++) a.at(i).bd[j]=0;
		a.at(i).index=2; // for neutral atom
		a.at(i).chg=0;	// for neutral atom
	}
	
	/*
	if (a!=NULL) {
		for (i=0;i<num;i++) {
			if (a[i].order!=NULL) {
				delete [] a[i].order;
				a[i].order=NULL;
			}		
		}	
		delete [] a;
		a=NULL;
	}
	if (a==NULL) {
		a = new ATOM [num];
		for (i=0;i<num;i++) a[i].order=new int [6];
		for (i=0;i<num;i++) {
			a[i].nh=0;
			for (j=0;j<6;j++) a[i].order[j]=0;
			for (j=0;j<3;j++) a[i].bd[j]=0;
			a[i].index=2; // for neutral atom
			a[i].chg=0;	// for neutral atom
			a[i].ang0=9999;
		}
	}
	*/
	double f=M_PI/180.0;
	a[0].id=0;
	a[0].name="H(-)";
	a[0].nbond=4;
	a[0].order[0]=1;
	a[0].order[1]=0;
	a[0].order[2]=0;
	a[0].order[3]=0;
	a[0].bd[0]=1;
	a[0].norder=1;
	a[0].rb=0.37;
	a[0].ang0=0.0;
	a[0].atm="H";
	a[0].probability=0.0;

	a[1].id=1;
	a[1].name="C(-)(-)(-)(-)";
	a[1].order[0]=a[1].order[1]=a[1].order[2]=a[1].order[3]=1;
	a[1].nbond=13;
	a[1].norder=4;
	a[1].rb=0.77;
	a[1].type=1;
	a[1].bd[0]=4;
	a[1].ang0=109.5*f;
	a[1].atm="C";
	a[1].probability=0.8*0.3;
	a[1].index=2;

	a[2].id=2;
	a[2].name="C(=)(-)(-)";
	a[2].order[0]=2;
	a[2].order[1]=a[2].order[2]=1;
	a[2].order[3]=0;
	a[2].nbond=10;
	a[2].norder=3;
	a[2].rb=0.77;
	a[2].type=3;
	a[2].bd[0]=2;
	a[2].bd[1]=1;
	a[2].ang0=120.0*f;
	a[2].atm="C";
	a[2].probability=0.8*0.3;

	a[3].id=3;
	a[3].name="C(#)(-)";
	a[3].order[0]=3;
	a[3].order[1]=1;
	a[3].order[2]=a[3].order[3]=0;
	a[3].nbond=7;
	a[3].norder=2;
	a[3].rb=0.77;
	a[3].type=4;
	a[3].bd[0]=1;
	a[3].bd[2]=1;
	a[3].ang0=180*f;
	a[3].atm="C";
	a[3].probability=0.2*0.8;

	a[4].id=4;
	a[4].name="C(=)(=)";
	a[4].order[0]=a[4].order[1]=2;
	a[4].order[2]=a[4].order[3]=0;
	a[4].nbond=7;
	a[4].norder=2;
	a[4].rb=0.77;
	a[4].type=4;
	a[4].bd[1]=2;
	a[4].ang0=180*f;
	a[4].atm="C";
	a[4].probability=0.01*0.8;

	a[5].id=5;
	a[5].name="O(-)(-)";
	a[5].order[0]=a[5].order[1]=1;
	a[5].order[2]=a[5].order[3]=0;
	a[5].nbond=7;
	a[5].norder=2;
	a[5].rb=0.74;
	a[5].type=4;
	a[5].bd[0]=2;
	a[5].ang0=104.5*f;
	a[5].atm="O";
	a[5].probability=0.4*0.3;

	a[6].id=6;
	a[6].name="O(=)";
	a[6].order[0]=2;
	a[6].order[1]=a[6].order[2]=a[6].order[3]=0;
	a[6].nbond=4;
	a[6].norder=1;
	a[6].rb=0.74;
	a[6].type=7;
	a[6].bd[1]=1;
	a[6].ang0=0;
	a[6].atm="O";
	a[6].probability=0.4*0.3;

	a[7].id=7;
	a[7].name="N(-)(-)(-)";
	a[7].order[0]=a[7].order[1]=a[7].order[2]=1;
	a[7].order[3]=0;
	a[7].nbond=10;
	a[7].norder=3;
	a[7].rb=0.74;
	a[7].type=2;
	a[7].bd[0]=3;
	a[7].ang0=107.5*f;
	a[7].atm="N";
	a[7].probability=0.25*0.2;

	a[8].id=8;
	a[8].name="N(=)(-)";
	a[8].order[0]=2;
	a[8].order[1]=1;
	a[8].order[2]=a[8].order[3]=0;
	a[8].nbond=7;
	a[8].norder=2;
	a[8].rb=0.74;
	a[8].type=4;
	a[8].bd[0]=1;
	a[8].bd[1]=1;
	a[8].ang0=115.0*f;
	a[8].atm="N";
	a[8].probability=0.2*0.25;

	a[9].id=9;
	a[9].name="N(#)";
	a[9].order[0]=3;
	a[9].order[1]=a[9].order[2]=a[9].order[3]=0;
	a[9].nbond=4;
	a[9].norder=1;
	a[9].rb=0.74;
	a[9].type=7;
	a[9].bd[2]=1;
	a[9].ang0=0;
	a[9].atm="N";
	a[9].probability=0.2*0.1;

	a[10].id=10;
	a[10].name="O(-)";
	a[10].order[0] = 1;
	a[10].order[1] = a[10].order[2] =a[10].order[3] =0;
	a[10].nbond=4;
	a[10].norder=1;
	a[10].nh=1;
	a[10].rb=0.74;
	a[10].type=4;
	a[10].bd[0]=1;
	a[10].ang0=0;
	a[10].atm="O";
	a[10].probability=0.3*0.3;

	a[11].id = 11;
	a[11].name = "F(-)";
	a[11].order[0] = 1;
	a[11].order[1] = a[11].order[2] =a[11].order[3] =0;
	a[11].nbond=4;
	a[11].norder=1;
	a[11].rb=0.72;
	a[11].type=7;
	a[11].bd[0]=1;
	a[11].ang0=0;
	a[11].atm="F";
	a[11].probability=0.1*0.05;

	a[12].id = 12;
	a[12].name="Cl(-)";
	a[12].order[0] = 1;
	a[12].order[1]=a[12].order[2]=a[12].order[3] = 0;
	a[12].nbond=5;
	a[12].norder=1;
	a[12].type=7;
	a[12].index=3;
	a[12].rb=0.99;
	a[12].bd[0]=1;
	a[12].ang0=0.0;
	a[12].atm="Cl";	
	a[12].probability=0.1*0.05;

	a[13].id = 13;
    a[13].name="Br(-)";
    a[13].order[0] = 1;
    a[13].order[1] = a[13].order[2] = a[13].order[3] =0;
    a[13].nbond = 5;
    a[13].norder=1;
	a[13].rb=1.200000;
	a[13].type=7;
	a[13].bd[0]=1;
	a[13].index=3;
	a[13].ang0=0;
	a[13].atm="Br";
	a[13].probability=0.1*0.05;

	a[14].id = 14;
	a[14].name ="I(-)";
	a[14].order[0]=1;
	a[14].order[1] = a[14].order[2] = a[14].order[3] = 0;
	a[14].nbond = 4;
	a[14].norder=1;
	a[14].rb=1.3900000;
	a[14].type=7;
	a[14].bd[0]=1;
	a[14].ang0=0;
	a[14].atm="I";
	a[14].probability=0;
	
	a[15].id=15;
	a[15].name="[NH0+](-)(-)(-)(-)";
	for (i=0;i<4;i++) a[15].order[i]=1;
	a[15].nbond=18;
	a[15].norder=4;
	a[15].rb=0.74;
	a[15].type=1;
	a[15].bd[0]=4;
	a[15].chg=1;
	a[15].index=7;
	a[15].ang0=109.5*f;
	a[15].atm="N";
	a[15].probability=0.005;

	a[16].id=16;
	a[16].name="[NH0+](=)(-)(-)";
	a[16].order[0]=2;
	a[16].order[1]=a[16].order[2]=1;
	a[16].order[3]=0;
	a[16].nbond=15;
	a[16].norder=3;
	a[16].rb=0.74;
	a[16].type=3;
	a[16].bd[0]=2;
	a[16].bd[1]=1;
	a[16].chg=1;
	a[16].index=7;
	a[16].ang0=120.0*f;
	a[16].atm="N";
	a[16].probability=0.005;

	a[17].id=17;
	a[17].name="[PH0+](-)(-)(-)(-)";
	for (i=0;i<4;i++) a[17].order[i]=1;
	a[17].norder=4;
	a[17].nbond=18;
	a[17].rb=1.07;
	a[17].type=1;
	a[17].bd[0]=4;
	a[17].chg=1;
	a[17].index=7;
	a[17].ang0=109.5*f;
	a[17].atm="P";
	a[17].probability=0.005;

	a[18].id=18;
	a[18].name="[PH0+](=)(-)(-)";
	a[18].order[0]=2;
	a[18].order[1]=a[18].order[2]=1;
	a[18].order[3]=0;
	a[18].norder=3;
	a[18].nbond=15;
	a[18].rb=1.07;
	a[18].bd[0]=2;
	a[18].bd[1]=1;
	a[18].chg=1;
	a[18].index=7;
	a[18].ang0=120.0*f;
	a[18].atm="P";
	a[18].probability=-0.001; //0.005

	a[19].id=19;
	a[19].name="S(-)(-)";
	a[19].order[0]=a[19].order[1]=1;
	a[19].order[2]=a[19].order[3]=0;
	a[19].norder=2;
	a[19].nbond=7;
	a[19].rb=1.05;
	a[19].type=4;
	a[19].bd[0]=2;
	a[19].ang0=104.5*f;
	a[19].atm="S";
	a[19].probability=0.01*0.2;

	a[20].id=20;
	a[20].name="S(=)";
	a[20].order[0]=2;
	a[20].order[1]=a[20].order[2]=a[20].order[3]=0;
	a[20].nbond=4;
	a[20].norder=1;
	a[20].rb=1.05;
	a[20].type=7;
	a[20].bd[1]=1;
	a[20].ang0=0;
	a[20].atm="S";
	a[20].probability=0.01*0.2;

	a[21].id=21;
	a[21].name="P(-)(-)(-)";
	a[21].order[0]=a[21].order[1]=a[21].order[2]=1;
	a[21].order[3]=0;
	a[21].nbond=10;
	a[21].norder=3;
	a[21].rb=1.07;
	a[21].type=2;
	a[21].bd[0]=3;
	a[21].ang0=107.5*f;
	a[21].atm="P";
	a[21].probability=0.01*0.1;

	a[22].id=22;
	a[22].name="P(=)(-)";
	a[22].order[0]=2;
	a[22].order[1]=1;
	a[22].order[2]=a[22].order[3]=0;
	a[22].nbond=7;
	a[22].norder=2;
	a[22].rb=1.07;
	a[22].type=4;
	a[22].bd[0]=1;
	a[22].bd[1]=1;
	a[22].ang0=115.0*f;
	a[22].atm="P";
	a[22].probability=0.001*0.1;
	
	a[23].id=23;
	a[23].name="P(#)";
	a[23].order[0]=3;
	a[23].order[1]=a[23].order[2]=a[23].order[3]=0;
	a[23].nbond=4;
	a[23].norder=1;
	a[23].rb=1.07;
	a[23].type=7;
	a[23].bd[2]=1;
	a[23].ang0=0;
	a[23].atm="P";
	a[23].probability=0.0001*0.1;

	a[24].id=24;
	a[24].name="[F-]";
	for (i=0;i<4;i++) a[24].order[i]=0;
	a[24].nbond=4;
	a[24].norder=0;
	a[24].rb=0.72;
	a[24].type=0;
	a[24].chg=-1;
	a[24].atm="F";
	a[24].index=4; //0
	a[24].probability=0.005;

	a[25].id=25;
	a[25].name="[Cl-]";
	for (i=0;i<4;i++) a[25].order[i]=0;
	a[25].nbond=5;
	a[25].norder=0;
	a[25].rb=0.99;
	a[25].type=0;
	a[25].chg=-1;
	a[25].atm="Cl";
	a[25].probability=0.005;
	a[25].index=5; //0

	a[26].id=26;
	a[26].name="[Br-]";
	for (i=0;i<4;i++) a[26].order[i]=0;
	a[26].nbond=5;
	a[26].norder=0;
	a[26].rb=1.20;
	a[26].type=0;
	a[26].chg=-1;
	a[26].atm="Br";
	a[26].probability=0.005;
	a[26].index=5; //0

	a[27].id=27;
	a[27].name="[I-]";
	for (i=0;i<4;i++) a[27].order[i]=0;
	a[27].nbond=4;
	a[27].norder=0;
	a[27].rb=1.39;
	a[27].type=0;
	a[27].chg=-1;
	a[27].atm="I";
	a[27].index=4; //0
	a[27].probability=0.000;

	a[28].id=28;
	a[28].name="[OH1-]";
	for (i=0;i<4;i++) a[28].order[i]=0;
	a[28].nbond=6;
	a[28].norder=0;
	a[28].chg=-1;
	a[28].nh=1;
	a[28].rb=0.74;
	a[28].type=0;	
	a[28].atm="O";
	a[28].index=6; //0
	a[28].probability=0.005;

	a[29].id=29;
	a[29].name="[OH0-](-)";
	a[29].order[0]=1;
	a[29].nbond=9;
	a[29].norder=1;
	a[29].chg=-1;
	a[29].rb=0.74;
	a[29].type=7;
	a[29].bd[0]=1;
	a[29].index=7;	
	a[29].atm="O";
	a[29].probability=0.005;

	a[30].id=30;
	a[30].name="[PH0-](-)(-)(-)(-)(-)(-)";
	for (i=0;i<6;i++) a[30].order[i]=1;
	a[30].index=7;
	a[30].nbond=24;
	a[30].norder=6;
	a[30].chg=-1;
	a[30].rb=1.07;
	a[30].type=6;
	a[30].bd[0]=6;
	a[30].ang0=1.82557;
	a[30].atm="P";
	a[30].probability=0.005;

	a[31].id=31;
	a[31].name="P(-)(-)(-)(-)(-)";
	for (i=0;i<5;i++) a[31].order[i]=1;
	a[31].nbond=16;
	a[31].norder=5;
	a[31].chg=0;
	a[31].rb=1.07;
	a[31].type=5;
	a[31].bd[0]=5;
	a[31].ang0=90.0*f;
	a[31].atm="P";
	a[31].probability=0.00001;

	a[32].id=32;
	a[32].name="P(=)(-)(-)(-)";
	a[32].order[0]=2;
	for (i=1;i<4;i++) a[32].order[i]=1;
	a[32].nbond=13;
	a[32].norder=4;
	a[32].chg=0;
	a[32].rb=1.07;
	a[32].type=1;
	a[32].bd[0]=3;
	a[32].bd[1]=1;
	a[32].ang0=109.5*f;
	a[32].atm="P";
	a[32].probability=0.0;

	a[33].id=33;
	a[33].name="S(-)(=O)(=O)([O-])";
	a[33].order[0]=1;
	a[33].index=2;
	a[33].nbond=4;
	a[33].norder=1;
	a[33].chg=-1;
	a[33].rb=1.05;	
	a[33].type=3;
	a[33].bd[0]=1;
	a[33].ang0=109.5*f;
	a[33].atm="S";
	a[33].probability=0.005;

	a[34].id=34;
	a[34].name="S(=)(-)(-)";
	a[34].order[0]=2;
	a[34].order[1]=a[34].order[2]=1;
	a[34].norder=3;
	a[34].nbond=10;
	a[34].chg=0;
	a[34].rb=1.05;
	a[34].type=3;
	a[34].bd[0]=2;
	a[34].bd[1]=1;
	a[34].ang0=120.0*f;
	a[34].atm="S";
	a[34].probability=0.0;

	a[35].id=35;
	a[35].name="[NH0-](-)(-)";
	a[35].order[0]=a[35].order[1]=1;
	a[35].nbond=12;
	a[35].norder=2;
	a[35].chg=-1;
	a[35].rb=0.74;
	a[35].type=4;
	a[35].bd[0]=2;
	a[35].index=7;
	a[35].ang0=104.5*f;
	a[35].atm="N";	
	a[35].probability=0.005;

	a[36].id=36; //GG
	a[36].name="C(-)(-)(-)([N+]1C=CN(C)C=1)";
    a[36].order[0]=a[36].order[1]=a[36].order[2]=1;
    a[36].nbond=10;
    a[36].norder=3;
    a[36].chg=1;
    a[36].rb=0.77;
    a[36].type=1;
    a[36].bd[0]=3;
    a[36].index=2;
    a[36].ang0=109.5*f;
    a[36].atm="C";
    a[36].probability=0.005;

    a[37].id=37; //GG
    a[37].name="C(-)(-)(-)([N+]1C=CN(C)C(C)=1)";
    a[37].order[0]=a[37].order[1]=a[37].order[2]=1;
    a[37].nbond=10;
    a[37].norder=3;
    a[37].chg=1;
    a[37].rb=0.77;
    a[37].type=1;
    a[37].bd[0]=3;
    a[37].index=2;
    a[37].ang0=109.5*f;
    a[37].atm="C";
    a[37].probability=0.005;

    a[38].id=38;
    a[38].name="C(-)(-)(-)(N1C=C[N+](C)=C1)";
    a[38].order[0]=a[38].order[1]=a[38].order[2]=1;
    a[38].nbond=10;
    a[38].norder=3;
    a[38].chg=1;
    a[38].rb=0.77;
    a[38].type=1;
    a[38].bd[0]=3;
    a[38].index=2;
    a[38].ang0=109.5*f;
    a[38].atm="C";
    a[38].probability=0.005;

    a[39].id=39; //GG
    a[39].name="C(-)(-)(-)([N+]1=CC=CC(C)=C1)";
    a[39].order[0]=a[39].order[1]=a[39].order[2]=1;
    a[39].nbond=10;
    a[39].norder=3;
    a[39].chg=1;
    a[39].rb=0.77;
    a[39].type=1;
    a[39].bd[0]=3;
    a[39].index=2;
    a[39].ang0=109.5*f;
    a[39].atm="C";
    a[39].probability=0.005;

    a[40].id=40; //GG
    a[40].name="C(-)(-)(-)(C1=C[N+](C)=CC=C1)";
    a[40].order[0]=a[40].order[1]=a[40].order[2]=1;
    a[40].nbond=10;
    a[40].norder=3;
    a[40].chg=1;
    a[40].rb=0.77;
    a[40].type=1;
    a[40].bd[0]=3;
    a[40].index=2;
    a[40].ang0=109.5*f;
    a[40].atm="C";
    a[40].probability=0.005;

    a[41].id=41;
    a[41].name="C(-)(1=[NH+]C=CC=C1)";
    a[41].order[0]=1;
    a[41].nbond=4;
    a[41].norder=1;
    a[41].chg=1;
    a[41].rb=0.77;
    a[41].type=3;
    a[41].bd[0]=1;
    a[41].index=2;
    a[41].ang0=120.0*f;
    a[41].atm="C";
    a[41].probability=0.005;

    a[42].id=42; 
    a[42].name="[NH0+](-)(1=CC=CC=C1)";
    a[42].order[0]=1;
    a[42].nbond=9;
    a[42].norder=1;
    a[42].chg=1;
    a[42].rb=0.74;
    a[42].type=3;
    a[42].bd[0]=1;
    a[42].index=7;
    a[42].ang0=120.0*f;
    a[42].atm="N";
    a[42].probability=0.005;

    a[43].id=43;
    a[43].name="[NH0-](S(=O)(=O)C(F)(F)(F))(S(=O)(=O)C(F)(F)(F))";
    a[43].nbond=48; //6
    a[43].norder=0;
    a[43].chg=-1;
    a[43].rb=0.74;
    a[43].type=4;
    a[43].index=48;  //0
    a[43].ang0=104.5*f;
    a[43].atm="N";
    a[43].probability=0.005;

	a[44].id=44;
	a[44].name="[BH0-](-)(-)(-)(-)";
    a[44].order[0]=a[44].order[1]=a[44].order[2]=a[44].order[3]=1;
    a[44].nbond=18;
    a[44].norder=4;
    a[44].chg=-1;
    a[44].rb=0.82;  // Wikipedia
    a[44].type=1;
    a[44].bd[0]=4;
    a[44].index=7;
    a[44].ang0=109.5*f;
    a[44].atm="B";
    a[44].probability=0.005;

	a[45].id=45;
	a[45].name="C(-)(-)(-)(C(=O)([O-]))";
    a[45].order[0]=a[45].order[1]=a[45].order[2]=1;
    a[45].nbond=10;
    a[45].norder=3;
    a[45].chg=-1;
    a[45].rb=0.77;  
    a[45].type=1;
    a[45].bd[0]=3;
    a[45].index=2;
    a[45].ang0=109.5*f;
    a[45].atm="C";
    a[45].probability=0.005;

	a[46].id=46;
	a[46].name="C(#N)([S-])";
    a[46].nbond=1;
    a[46].norder=0;
    a[46].chg=-1;
    a[46].rb=0.77;  
    a[46].type=4;
    a[46].index=0;
    a[46].ang0=180.0*f;
    a[46].atm="C";
    a[46].probability=0.005;

    a[47].id=47;
    a[47].name="C(-)(-)(-)(OP(=O)(OC)([O-]))";
    a[47].order[0]=a[47].order[1]=a[47].order[2]=1;
    a[47].nbond=10;
    a[47].norder=3;
    a[47].chg=-1;
    a[47].rb=0.77;  
    a[47].type=1;
    a[47].bd[0]=3;
    a[47].index=2;
    a[47].ang0=109.5*f;
    a[47].atm="C";
    a[47].probability=0.005;

    a[48].id=48;
    a[48].name="C(#N)([N-]C#N)";
    a[48].nbond=14; //1
    a[48].norder=0;
    a[48].chg=-1;
    a[48].rb=0.77;  
    a[48].type=4;
    a[48].index=14; //0
    a[48].ang0=180.0*f;
    a[48].atm="C";
    a[48].probability=0.005;

// 20191109

    a[49].id=49;
    a[49].name="[BH0-](C#N)(C#N)(C#N)(C#N)"; // N#C[B-](C#N)(C#N)C#N
    a[49].nbond=26;  //4
    a[49].norder=0;
    a[49].chg=-1;
    a[49].rb=0.82;
    a[49].type=1;
    a[49].index=26; //0
    a[49].ang0=109.5*f;
    a[49].atm="B";
    a[49].probability=0.005;

    a[50].id=50;
    a[50].name="S(OCCOCCOC)(=O)(=O)([O-])"; // COCCOCCOS(=O)(=O)[O-]
    a[50].nbond=25; //1
    a[50].norder=0;
    a[50].chg=-1;
    a[50].rb=0.82;
    a[50].type=1;
    a[50].index=25; //0
    a[50].ang0=109.5*f;
    a[50].atm="S";
    a[50].probability=0.005;

    a[51].id=51;
    a[51].name="S(c(cc1)ccc1C)(=O)(=O)([O-])"; // Cc1ccc(cc1)S(=O)(=O)[O-]
    a[51].nbond=28; //1
    a[51].norder=0;
    a[51].chg=-1;
    a[51].rb=0.82;
    a[51].type=1;
    a[51].index=28; //0
    a[51].ang0=109.5*f;
    a[51].atm="S";
    a[51].probability=0.005;

    a[52].id=52;
    a[52].name="[PH0-](F)(F)(F)(C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)"; // F[P-](C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F
    a[52].nbond=66; //6
    a[52].norder=0;
    a[52].chg=-1;
    a[52].rb=1.07;
    a[52].type=6;
    a[52].ang0=1.82557;
    a[52].atm="P";
    a[52].index=66; //0
    a[52].probability=0.005;

    a[53].id=53;
    a[53].name="[NH0-](S(=O)(=O)C(C(F)(F)F)(F)F)(S(=O)(=O)C(C(F)(F)F)(F)F)"; // FC(S(=O)(=O)[N-]S(=O)(=O)C(C(F)(F)F)(F)F)(C(F)(F)F)F
    a[53].nbond=58; //6
    a[53].norder=0;
    a[53].chg=-1;
    a[53].rb=0.74;
    a[53].type=4;
    a[53].ang0=104.5*f;
    a[53].atm="N";
    a[53].index=58; //0
    a[53].probability=0.005;

    a[54].id=54;
    a[54].name="[CH0-](S(C(F)(F)(F))(=O)(=O))(S(C(F)(F)(F))(=O)(=O))(S(C(F)(F)(F))(=O)(=O))"; // FC(S(=O)(=O)[C-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F)(F)F
    a[54].nbond=75; //6
    a[54].norder=0;
    a[54].chg=-1;
    a[54].rb=0.77;
    a[54].type=3;
    a[54].ang0=120.0*f;
    a[54].atm="C";
    a[54].index=75; //0
    a[54].probability=0.005;

    a[55].id=55;
    a[55].name="[PH0-](F)(F)(F)(F)(F)(F)"; // [PH0-](F)(F)(F)(F)(F)(F)
    a[55].nbond=24; //6
    a[55].norder=0;
    a[55].chg=-1;
    a[55].rb=0.77;
    a[55].type=3;
    a[55].ang0=120.0*f;
    a[55].atm="C";
    a[55].index=24; //0
    a[55].probability=0.005;

    a[56].id=56;
    a[56].name="[In+3]([Cl-])([Cl-])([Cl-])([Cl-])"; // [Cl-][In+3]([Cl-])([Cl-])[Cl-]
    a[56].nbond=34; //6
    a[56].norder=0;
    a[56].chg=-1;
    a[56].rb=1.44;
    a[56].type=1;
    a[56].ang0=90.0*f;
    a[56].atm="In";
    a[56].index=34; //0
    a[56].probability=-1;

    a[57].id=57;
    a[57].name="[In+3](-)(-)(-)(-)";
    a[57].order[0]=a[57].order[1]=a[57].order[2]=a[57].order[3]=1;
    a[57].nbond=18;
    a[57].norder=4;
    a[57].chg=3;
    a[57].rb=1.44;  // Wikipedia
    a[57].type=1;
    a[57].bd[0]=4;
    a[57].index=7;
    a[57].ang0=109.5*f;
    a[57].atm="In";
    a[57].probability=-1;

    a[58].id=58;
    a[58].name="Cl(=O)(=O)([O-])(=O)"; // [O-][Cl](=O)(=O)=O
    a[58].nbond=20; //2
    a[58].norder=0;
    a[58].chg=-1;
    a[58].rb=0.99;
    a[58].type=1;
    a[58].ang0=90.0*f;
    a[58].atm="Cl";
    a[58].index=20; //0
    a[58].probability=0.005;

    a[59].id=59;
    a[59].name="[CH0-](-)(-)(-)"; 
	a[59].order[0]=a[59].order[1]=a[59].order[2]=1;
    a[59].nbond=15;
    a[59].norder=3;
    a[59].chg=-1;
    a[59].rb=0.74;
    a[59].type=3;
	a[59].bd[0]=3;
    a[59].ang0=120.0*f;
    a[59].atm="C";
    a[59].index=7;
    a[59].probability=0.005;

    a[60].id=60;
    a[60].name="[NH0+](=O)([O-])([O-])"; // [O-][N+](=O)[O-]
    a[60].nbond=22; //6
    a[60].norder=0;
    a[60].chg=-1;
    a[60].rb=0.74;
    a[60].type=3;
    a[60].ang0=90.0*f;
    a[60].atm="N";
    a[60].index=22; //0
    a[60].probability=0.005;

    a[61].id=61;
    a[61].name="S(=)(=)(-)(-)";
    a[61].order[0]=a[61].order[1]=2;
	a[61].order[2]=a[61].order[3]=1;
    a[61].index=2;
    a[61].nbond=13;
    a[61].norder=4;
    a[61].chg=0;
    a[61].rb=1.05;
    a[61].type=3;
    a[61].bd[0]=a[61].bd[1]=2;
    a[61].ang0=109.5*f;
    a[61].atm="S";
    a[61].probability=-0.001;

    a[62].id=62;
    a[62].name="Cl(=)(=)(=)(-)"; // [O-][Cl](=O)(=O)=O
    a[62].order[0]=a[62].order[1]=a[62].order[2]=2;
	a[62].order[3]=1;
    a[62].nbond=14;
    a[62].norder=4;
    a[62].chg=0;
    a[62].rb=0.99;
    a[62].type=1;
	a[62].bd[0]=1;
    a[62].bd[1]=3;
    a[62].ang0=90.0*f;
    a[62].atm="Cl";
    a[62].index=3;
    a[62].probability=-0.001;

    a[63].id=63;
    a[63].name="[SH0-](-)";
	a[63].order[0]=1;
    a[63].nbond=9;
    a[63].norder=1;
    a[63].chg=-1;
    a[63].rb=1.05;
    a[63].type=7;
	a[63].bd[0]=1;
    a[63].index=7;
    a[63].ang0=180.0*f;
    a[63].atm="S";
    a[63].probability=0.005;

    a[64].id=64;
    a[64].name="[Ga+3](-)(-)(-)(-)";
    a[64].order[0]=a[64].order[1]=a[64].order[2]=a[64].order[3]=1;
    a[64].nbond=18;
    a[64].norder=4;
    a[64].chg=3;
    a[64].rb=1.44;  // Wikipedia
    a[64].type=1;
    a[64].bd[0]=4;
    a[64].index=7;
    a[64].ang0=109.5*f;
    a[64].atm="Ga";
    a[64].probability=-1;

    a[65].id=65;
    a[65].name="[SH0+](-)(-)(-)";
    a[65].order[0]=a[65].order[1]=a[65].order[2]=1;
    a[65].nbond=15;
    a[65].norder=3;
    a[65].chg=1;
    a[65].rb=1.05;
    a[65].type=7;
    a[65].bd[0]=3;
    a[65].index=7;
    a[65].ang0=180.0*f;
    a[65].atm="S";
    a[65].probability=-1;

    a[66].id=66;
    a[66].name="P(=)(-)(-)";
    a[66].order[0]=2;
    for (i=1;i<3;i++) a[66].order[i]=1;
    a[66].nbond=10;
    a[66].norder=3;
    a[66].chg=0;
    a[66].rb=1.07;
    a[66].type=1;
    a[66].bd[0]=2;
    a[66].bd[1]=1;
    a[66].ang0=109.5*f;
    a[66].atm="P";
    a[66].probability=0.001;

    a[67].id=67;
    a[67].name="[CH0@@](-)(-)(-)(-)";
    a[67].order[0]=a[67].order[1]=a[67].order[2]=a[67].order[3]=1;
    a[67].nbond=19;
    a[67].norder=4;
    a[67].rb=0.77;
    a[67].type=1;
    a[67].bd[0]=4;
    a[67].ang0=109.5*f;
    a[67].atm="C";
    a[67].probability=-1;
    a[67].index=8;

    a[68].id=68;
    a[68].name="[CH0@](-)(-)(-)(-)";
    a[68].order[0]=a[68].order[1]=a[68].order[2]=a[68].order[3]=1;
    a[68].nbond=18;
    a[68].norder=4;
    a[68].rb=0.77;
    a[68].type=1;
    a[68].bd[0]=4;
    a[68].ang0=109.5*f;
    a[68].atm="C";
    a[68].probability=-1;
    a[68].index=7;

    a[69].id=69;
    a[69].name="P(-)(-)(-)(-)";
    for (i=0;i<4;i++) a[69].order[i]=1;
    a[69].nbond=13;
    a[69].norder=4;
    a[69].chg=0;
    a[69].rb=1.07;
    a[69].type=5;
    a[69].bd[0]=4;
    a[69].ang0=90.0*f;
    a[69].atm="P";
    a[69].probability=0.00001;


	return;
}
