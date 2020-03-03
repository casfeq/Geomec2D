/*
	This header is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The class defined 
	here contains the functions for determination of the poroelasticity problem's parameters.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.
*/

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class problemParameters
{
public:
	// Class variables
	double K; // [m^2]
	double mu_f; // [Pa.s]
	double rho_f; // [kg/m^3]
	double rho_s; // [kg/m^3]
	double rho; // [kg/m^3]
	double c_s; // [Pa^-1]
	double c_f; // [Pa^-1]
	double alpha; // [adim]
	double G; // [Pa]
	double lambda; // [Pa]
	double phi; // [adim]
	double sigmab; // [Pa]
	double F; // [N/m]
	double g; // [m/s^2]
	double Q;
	double c;
	double M;
	double dx;
	double dy;
	double dt_vv;
	double dt_carlos;
	double u0;
	double v0;
	double e0;
	double P0;
	double Lx;
	double Ly;
	vector<vector<double>> uDisplacementField;
	vector<vector<double>> vDisplacementField;
	vector<vector<double>> pressureField;
	int Nu;
	int Nv;
	int NP;
	vector<vector<int>> uDisplacementFVCoordinates;
	vector<vector<int>> vDisplacementFVCoordinates;
	vector<vector<int>> pressureFVCoordinates;
	vector<vector<int>> uDisplacementFVIndex;
	vector<vector<int>> vDisplacementFVIndex;
	vector<vector<int>> pressureFVIndex;

	// Class functions
	int getUDisplacementFVPosition(int,int);
	int getVDisplacementFVPosition(int,int);
	int getPressureFVPosition(int,int);
	void computeProblemParameters();
	void applyTerzaghiInitialConditions();
	void applySealedColumnInitialConditions();
	void applyMandelInitialConditions();

	// Constructor
	problemParameters(double,double,double,double,double,double,double,double,double,double,double,
		double,double,double,vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,
		vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,
		vector<vector<int>>,vector<vector<int>>,double);

	// Destructor
	~problemParameters();
};

problemParameters::problemParameters(double deltax, double deltay, double porousMediumPermeability,
	double porousMediumPorosity, double solidMatrixDensity, double solidMatrixCompressibility,
	double fluidViscosity, double fluidDensity, double fluidCompressibility, double shearModulus,
	double lames1stParameter, double prescribedNormalStress, double b, double H,
	vector<vector<double>> uField, vector<vector<double>> vField, vector<vector<double>> pField,
	vector<vector<int>> cooU, vector<vector<int>> cooV, vector<vector<int>> cooP,
	vector<vector<int>> idU, vector<vector<int>> idV, vector<vector<int>> idP, double gravity)
{
	K=porousMediumPermeability;
	phi=porousMediumPorosity;
	rho_s=solidMatrixDensity;
	c_s=solidMatrixCompressibility;
	mu_f=fluidViscosity;
	rho_f=fluidDensity;
	c_f=fluidCompressibility;
	G=shearModulus;
	lambda=lames1stParameter;
	sigmab=prescribedNormalStress;
	dx=deltax;
	dy=deltay;
	Lx=b;
	Ly=H;
	F=sigmab*Lx;
	uDisplacementField=uField;
	vDisplacementField=vField;
	pressureField=pField;
	Nu=uField.size();
	Nv=vField.size();
	NP=pField.size();
	uDisplacementFVCoordinates=cooU;
	vDisplacementFVCoordinates=cooV;
	pressureFVCoordinates=cooP;
	uDisplacementFVIndex=idU;
	vDisplacementFVIndex=idV;
	pressureFVIndex=idP;
	rho=rho_f*phi+rho_s*(1-phi);
	g=gravity;

	computeProblemParameters();
}

problemParameters::~problemParameters(){}

int problemParameters::getUDisplacementFVPosition(int x, int y)
{
	int uDisplacementFVPosition;

	uDisplacementFVPosition=uDisplacementFVIndex[x][y]-1;

	return uDisplacementFVPosition;
}

int problemParameters::getVDisplacementFVPosition(int x, int y)
{
	int vDisplacementFVPosition;

	vDisplacementFVPosition=vDisplacementFVIndex[x][y]+Nu-1;

	return vDisplacementFVPosition;
}

int problemParameters::getPressureFVPosition(int x, int y)
{
	int pressureFVPosition;

	pressureFVPosition=pressureFVIndex[x][y]+Nu+Nv-1;

	return pressureFVPosition;
}

void problemParameters::computeProblemParameters()
{
	alpha=1-c_s*(lambda+2/3*G);
	Q=1/(c_s*alpha+(c_f-c_s)*phi);
	M=2*G+lambda;
	c=(K*Q*M)/(mu_f*(M+Q*alpha*alpha));
	dt_vv=(dx*dy)/(6*c);
	dt_carlos=(dx*dy)/c;

	return;
}

void problemParameters::applySealedColumnInitialConditions()
{

	for(int FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		uDisplacementField[FVCounter][0]=0;
	}

	for(int FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		vDisplacementField[FVCounter][0]=0;
	}

	for(int FVCounter=0; FVCounter<NP; FVCounter++)
	{
		pressureField[FVCounter][0]=0;
	}

	return;
}

void problemParameters::applyTerzaghiInitialConditions()
{
	int i, j, v_P, P_P;
	double As, yValue, vValue, pValue;
	
	As=sigmab/(M+alpha*alpha*Q)+(rho*g*Ly)/(2*(M+alpha*alpha*Q))+0.5*(rho-alpha*rho_f)*(g*Ly/M);
	P0=0.5*rho_f*g*Ly+alpha*Q*(rho-alpha*rho_f)*g*Ly/(2*M)-alpha*Q*As;

	for(int FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		uDisplacementField[FVCounter][0]=0;
	}

	for(int FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		yValue=Ly-i*dy;
		vValue=-((rho-alpha*rho_f)*g*yValue*yValue)/(2*M)+As*yValue;
		vDisplacementField[v_P-Nu][0]=vValue;
	}

	for(int FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		yValue=(Ly-dy/2)-i*dy;
		pValue=P0-rho_f*g*(Ly-yValue);
		pressureField[P_P-Nu-Nv][0]=pValue;
	}

	return;
}

void problemParameters::applyMandelInitialConditions()
{
	int Nx=Lx/dx;
	int Ny=Ly/dy;
	int i, j;
	int u_P, v_P;

	e0=F/(Lx*(2*alpha*alpha*Q+M+lambda));
	P0=-(alpha*Q*F)/(Lx*(2*alpha*alpha*Q+M+lambda));

	for(int FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		u0=-((alpha*alpha*Q+lambda)/(2*G))*e0*j*dx;
		uDisplacementField[u_P][0]=u0;
	}

	for(int FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		v0=(M+alpha*alpha*Q)/(2*G)*e0*(Ly-i*dy);
		vDisplacementField[v_P-Nu][0]=v0;
	}

	for(int FVCounter=0; FVCounter<NP; FVCounter++)
	{
		pressureField[FVCounter][0]=P0;
	}

	return;
}