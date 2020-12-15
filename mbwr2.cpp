#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <cstdio>

using namespace std;

void MBWR(int pset, double rho, double Tr, double *Pcf, double *Gcf, double *B2, double *Ucf, double *Acf, double *Scf, double *Cvcf);
void IG(double rho, double Tr, double mass, double epslon, double sigma, double *Pig, double *Eig, double *Sig, double *Aig, double *Gig, double *Cvig);

int main()
{
	
	double rho=0.0,mass=0.0,epslon=0.0,sigma=0.0;
	char null[128];
	int  displaysize=10;

	cerr<<"Thermodynamics and Transport Properties of Lennard-Jones Fluid"<<endl;
	cerr<<"Modified Benedict-Webb-Rubin (MBWR) Equation of State for Lennard-Jones Fluid"<<endl;
	cerr<<"Ref: Mol. Phys., 78 (3), 591-618, 1993"<<endl<<endl;

	ifstream inf("input_Ar.txt",ios::in);
	if(!inf) {cout<<"File: input.txt could not be found!"<<endl; exit(1);}
	ofstream outf("output.txt",ios::out);

	int pset;
	inf>>null>>pset;
	inf>>null>>mass>>null>>epslon>>null>>sigma>>null>>null>>rho;
	outf<<"Total Properties"<<endl;
	outf<<"   rho*     T*       P*       E*       S*        A*        G*        Cv*"<<endl;
	outf<<setiosflags(ios::fixed)<<setprecision(displaysize-5);

	while(rho>0)
	{
		double Tr=0.0,Pcf=0.0,Gcf=0.0,B2=0.0,Ucf=0.0,Acf=0.0,Scf=0.0,Cvcf=0.0;
		double Pig=0.0,Eig=0.0,Sig=0.0,Aig=0.0,Gig=0.0,D=0.0,Cvig=0.0;

		inf>>Tr;
		MBWR(pset,rho,Tr,&Pcf,&Gcf,&B2,&Ucf,&Acf,&Scf,&Cvcf);
		IG(rho,Tr,mass,epslon,sigma,&Pig,&Eig,&Sig,&Aig,&Gig,&Cvig);

		outf<<setw(displaysize)<<rho<<setw(displaysize)<<Tr<<setw(displaysize)<<Pcf+Pig<<setw(displaysize)<<Ucf+Eig<<setw(displaysize)<<Scf+Sig<<setw(displaysize+1)<<Acf+Aig<<setw(displaysize+1)<<Gcf+Gig<<setw(displaysize+1)<<Cvcf+Cvig<<endl;
		inf>>rho;
	}

	inf.seekg(0);
	inf.clear();
	inf>>null>>pset;
	inf>>null>>mass>>null>>epslon>>null>>sigma>>null>>null>>rho;
	outf<<endl<<endl<<endl;
	outf<<"Residual Properties"<<endl;
	outf<<"   rho*     T*       P*       E*       S*        A*        G*        Cv*"<<endl;
	while(rho>0)
	{
		double Tr=0.0,Pcf=0.0,Gcf=0.0,B2=0.0,Ucf=0.0,Acf=0.0,Scf=0.0,Cvcf=0.0;
		double Pig=0.0,Eig=0.0,Sig=0.0,Aig=0.0,Gig=0.0,D=0.0,Cvig=0.0;

		inf>>Tr;
		MBWR(pset,rho,Tr,&Pcf,&Gcf,&B2,&Ucf,&Acf,&Scf,&Cvcf);

		outf<<setw(displaysize)<<rho<<setw(displaysize)<<Tr<<setw(displaysize)<<Pcf<<setw(displaysize)<<Ucf<<setw(displaysize)<<Scf<<setw(displaysize+1)<<Acf<<setw(displaysize+1)<<Gcf<<setw(displaysize+1)<<Cvcf<<endl;
		inf>>rho;
	}

	inf.close();
	outf.close();
       return 0;
}

void MBWR(int pset, double rho, double Tr, double *Pcf, double *Gcf, double *B2, double *Ucf, double *Acf, double *Scf, double *Cvcf)
{
	double x[32];
	if( pset ==1979 )
	{
		//1979 parameter
		x[ 0]= -0.44480725e-1;
		x[ 1]= +0.72738221e+1;
		x[ 2]= -0.14343368e+2;
		x[ 3]= +0.38397096e+1;
		x[ 4]= -0.20057745e+1;
		x[ 5]= +0.19084472e+1;
		x[ 6]= -0.57441787e+1;
		x[ 7]= +0.25110073e+2;
		x[ 8]= -0.45232787e+4;
		x[ 9]= +0.89327162e-2;
		x[10]= +0.98163358e+1;
		x[11]= -0.61434572e+2;
		x[12]= +0.14161454e+2;
		x[13]= +0.43353841e+2;
		x[14]= +0.11078327e+4;
		x[15]= -0.35429519e+2;
		x[16]= +0.10591298e+2;
		x[17]= +0.49770046e+3;
		x[18]= -0.35338542e+3;
		x[19]= +0.45036093e+4;
		x[20]= +0.77805296e+1;
		x[21]= +0.13567114e+5;
		x[22]= -0.85818023e+1;
		x[23]= +0.16646578e+5;
		x[24]= -0.14092234e+2;
		x[25]= +0.19386911e+5;
		x[26]= +0.38585868e+2;
		x[27]= +0.33800371e+4;
		x[28]= -0.18567754e+3;
		x[29]= +0.84874693e+4;
		x[30]= +0.97508689e+2;
		x[31]= -0.14483060e+2;
	} else 
	{
		//1993 parameter
		x[ 0]= +0.8623085097507421e+0;
		x[ 1]= +2.976218765822098e+0;
		x[ 2]= -8.402230115796038e+0;
		x[ 3]= +0.1054136629203355e+0;
		x[ 4]= -0.8564583828174598e+0;
		x[ 5]= +1.582759470107601e+0;
		x[ 6]= +0.7639421948305453e+0;
		x[ 7]= +1.753173414312048e+0;
		x[ 8]= +2.798291772190376e+3;
		x[ 9]= -4.8394220260857657e-2;
		x[10]= +0.9963265197721935e+0;
		x[11]= -3.698000291272493e+1;
		x[12]= +2.084012299434647e+1;
		x[13]= +8.305402124717285e+1;
		x[14]= -9.574799715203068e+2;
		x[15]= -1.477746229234994e+2;
		x[16]= +6.398607852471505e+1;
		x[17]= +1.603993673294834e+1;
		x[18]= +6.805916615864377e+1;
		x[19]= -2.791293578795945e+3;
		x[20]= -6.245128304568454e+0;
		x[21]= -8.116836104958410e+3;
		x[22]= +1.488735559561229e+1;
		x[23]= -1.059346754655084e+4;
		x[24]= -1.131607632802822e+2;
		x[25]= -8.867771540418822e+3;
		x[26]= -3.986982844450543e+1;
		x[27]= -4.689270299917261e+3;
		x[28]= +2.593535277438717e+2;
		x[29]= -2.694523589434903e+3;
		x[30]= -7.218487631550215e+2;
		x[31]= +1.721802063863269e+2; 
	}

	double A[9];
	{
		A[ 1]= (x[0]*Tr + x[1]*sqrt(Tr) + x[2] + x[3]/Tr + x[4]/Tr/Tr);
		A[ 2]= (x[5]*Tr + x[6] +x[7]/Tr + x[8]/Tr/Tr);
		A[ 3]= (x[9]*Tr + x[10] + x[11]/Tr);
		A[ 4]= x[12];
		A[ 5]= (x[13]/Tr + x[14]/Tr/Tr);
		A[ 6]= (x[15]/Tr);
		A[ 7]= (x[16]/Tr + x[17]/Tr/Tr);
		A[ 8]= (x[18]/Tr/Tr);
	}

	double B[7];
	{
		B[ 1]= (x[19]/Tr/Tr + x[20]/Tr/Tr/Tr);
		B[ 2]= (x[21]/Tr/Tr + x[22]/Tr/Tr/Tr/Tr);
		B[ 3]= (x[23]/Tr/Tr + x[24]/Tr/Tr/Tr);
		B[ 4]= (x[25]/Tr/Tr + x[26]/Tr/Tr/Tr/Tr);
		B[ 5]= (x[27]/Tr/Tr + x[28]/Tr/Tr/Tr);
		B[ 6]= (x[29]/Tr/Tr + x[30]/Tr/Tr/Tr + x[31]/Tr/Tr/Tr/Tr);
	}

	double r=3.0;
	double F=exp(-r*rho*rho);

	double G[7];
	{
		G[ 1]=(1-F)/(2*r);
		G[ 2]=-(F*pow(rho,2) - 2*G[1])/(2*r);
		G[ 3]=-(F*pow(rho,4) - 4*G[2])/(2*r);
		G[ 4]=-(F*pow(rho,6) - 6*G[3])/(2*r);
		G[ 5]=-(F*pow(rho,8) - 8*G[4])/(2*r);
		G[ 6]=-(F*pow(rho,10)-10*G[5])/(2*r);
	}

	double C[9];
	{
		C[ 1]=x[1]*sqrt(Tr)/2+x[2]+2*x[3]/Tr+3*x[4]/Tr/Tr;
		C[ 2]=x[6]+2*x[7]/Tr+3*x[8]/Tr/Tr;
		C[ 3]=x[10]+2*x[11]/Tr;
		C[ 4]=x[12];
		C[ 5]=2*x[13]/Tr+3*x[14]/Tr/Tr;
		C[ 6]=2*x[15]/Tr;
		C[ 7]=2*x[16]/Tr+3*x[17]/Tr/Tr;
		C[ 8]=3*x[18]/Tr/Tr;
	}

	double D[9];
	{
		D[ 1]=3*x[19]/Tr/Tr+4*x[20]/Tr/Tr/Tr;
		D[ 2]=3*x[21]/Tr/Tr+5*x[22]/Tr/Tr/Tr/Tr;
		D[ 3]=3*x[23]/Tr/Tr+4*x[24]/Tr/Tr/Tr;
		D[ 4]=3*x[25]/Tr/Tr+5*x[26]/Tr/Tr/Tr/Tr;
		D[ 5]=3*x[27]/Tr/Tr+4*x[28]/Tr/Tr/Tr;
		D[ 6]=3*x[29]/Tr/Tr+4*x[30]/Tr/Tr/Tr+5*x[31]/Tr/Tr/Tr/Tr;
	}

	//second virial coefficient
	*B2=3.0/2.0/3.14159265359*A[1]/Tr;

	int i;

	//Energy
	for(i=1;i<=8;i++)
	{
		if(i<=6)	*Ucf += C[i]*pow(rho,i)/(i) + D[i]*G[i];
		else *Ucf += C[i]*pow(rho,i)/(i);
	}

	//Pressure
	for(i=1;i<=8;i++)
	{
		if(i<=6)	*Pcf += A[i]*pow(rho,i+1) + F*B[i]*pow(rho,2*i+1);
		else *Pcf += A[i]*pow(rho,i+1);
	}

	//Helmholtz free energy
	for(i=1;i<=8;i++)
	{
		if(i<=6)	*Acf += A[i]*pow(rho,i)/(i) + B[i]*G[i];
		else *Acf += A[i]*pow(rho,i)/i;
	}

	//Chemical potential or Gibbs free energy
	*Gcf= *Acf+ *Pcf/rho;

	//Entropy
	*Scf= (*Ucf-*Acf)/Tr;

	//Constant volume heat capacity
	double CdT[9];
	{
		CdT[ 1]=x[1]/sqrt(Tr)/4-2*x[3]/Tr/Tr-6*x[4]/Tr/Tr/Tr;
		CdT[ 2]=-2*x[7]/Tr/Tr-6*x[8]/Tr/Tr/Tr;
		CdT[ 3]=-2*x[11]/Tr/Tr;
		CdT[ 4]=0;
		CdT[ 5]=-2*x[13]/Tr/Tr-6*x[14]/Tr/Tr/Tr;
		CdT[ 6]=-2*x[15]/Tr/Tr;
		CdT[ 7]=-2*x[16]/Tr/Tr-6*x[17]/Tr/Tr/Tr;
		CdT[ 8]=-6*x[18]/Tr/Tr/Tr;
	}

	double DdT[9];
	{
		DdT[ 1]=-6*x[19]/Tr/Tr/Tr-12*x[20]/Tr/Tr/Tr/Tr;
		DdT[ 2]=-6*x[21]/Tr/Tr/Tr-20*x[22]/Tr/Tr/Tr/Tr/Tr;
		DdT[ 3]=-6*x[23]/Tr/Tr/Tr-12*x[24]/Tr/Tr/Tr/Tr;
		DdT[ 4]=-6*x[25]/Tr/Tr/Tr-20*x[26]/Tr/Tr/Tr/Tr/Tr;
		DdT[ 5]=-6*x[27]/Tr/Tr/Tr-12*x[28]/Tr/Tr/Tr/Tr;
		DdT[ 6]=-6*x[29]/Tr/Tr/Tr-12*x[30]/Tr/Tr/Tr/Tr-20*x[31]/Tr/Tr/Tr/Tr/Tr;
	}

	for(i=1;i<=8;i++)
	{
		if(i<=6)	*Cvcf += CdT[i]*pow(rho,i)/(i) + DdT[i]*G[i];
		else *Cvcf += CdT[i]*pow(rho,i)/(i);
	}

}

void IG(double rho, double Tr, double mass, double epslon, double sigma, double *Pig, double *Eig, double *Sig, double *Aig, double *Gig, double *Cvig)
{
	double Rig=8.314;
	double pi=3.141592653;
	double na=6.0221367e23;
	double h=6.6260755e-34;
	*Pig=rho*Tr;
	*Eig=3.0/2.0*Tr;
	*Sig=(5.0/2.0+log( pow((2.0*pi*(mass*Tr*4.184*epslon)/(h*na*h*na)),(3.0/2.0)) *1e-30/rho*sigma*sigma*sigma));
	*Aig=(*Eig)-Tr*(*Sig);
	*Gig=(*Eig)+(*Pig)/rho-Tr*(*Sig);
	*Cvig=3.0/2.0;
}

