/*
	The class defined here contains the functions for determination of the poroelasticity problem's parameters for materials with double porosity.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.
*/

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class problemDoubleParameters
{
public:
	// Class variables
	double phiPore, phiFrac;
	double c_f, c_s, c_b, G, lambda;
	double mu_f, KPore,KFrac;
	double alpha, S11, S12, S22, consolCoef;
	double psiPore, psiFrac;
	double leakTerm;
	double sigmab;
	vector<vector<int>> vDisplacementFVCoordinates;
	vector<vector<int>> vDisplacementFVIndex;
	vector<vector<double>> vDisplacementField;
	vector<vector<double>> pressurePoreField;
	vector<vector<double>> pressureFracField;

	// Class functions
	int getVDisplacementFVPosition(int,int);
	void computeProblemParameters();
	double computeLeakTerm(double);
	void applyDoublePoreInitialConditions(double,double,vector<vector<double>>,
		vector<vector<double>>,vector<vector<double>>);

	// Constructor
	problemDoubleParameters(double,double,double,double,double,double,double,double,double,double,
		vector<vector<int>>,vector<vector<int>>);

	// Destructor
	~problemDoubleParameters();
};

problemDoubleParameters::problemDoubleParameters(double porosityPore, double porosityFrac,
	double fluidCompressibility, double solidCompressibility, double shearModulus,
	double lame1stParameter, double fluidViscosity, double permeabilityPore,
	double permeabilityFrac, double prescribedStress, vector<vector<int>> cooV,
	vector<vector<int>> idV)
{
	phiPore=porosityPore;
	phiFrac=porosityFrac;
	c_f=fluidCompressibility;
	c_s=solidCompressibility;
	G=shearModulus;
	lambda=lame1stParameter;
	c_b=1/(lambda+2*G/3);
	mu_f=fluidViscosity;
	KPore=permeabilityPore;
	KFrac=permeabilityFrac;
	sigmab=prescribedStress;
	vDisplacementFVCoordinates=cooV;
	vDisplacementFVIndex=idV;

	computeProblemParameters();
}

problemDoubleParameters::~problemDoubleParameters(){}

int problemDoubleParameters::getVDisplacementFVPosition(int x, int y)
{
	int vDisplacementFVPosition;

	vDisplacementFVPosition=vDisplacementFVIndex[x][y]-1;

	return vDisplacementFVPosition;
}

void problemDoubleParameters::computeProblemParameters()
{
	psiPore=phiPore/(phiPore+phiFrac);
	psiFrac=phiFrac/(phiPore+phiFrac);
	alpha=1-(c_s/c_b);
	S11=phiPore*(c_f-c_s)+psiPore*alpha*(1-psiPore*alpha-phiFrac)*c_b;
	S12=psiPore*alpha*(psiFrac*alpha-phiFrac)*c_b;
	S22=phiFrac*(c_f-c_s)+psiFrac*alpha*(1-psiFrac*alpha)*c_b-phiFrac*psiPore*alpha*c_b;

	double M=2*G+lambda;
	double consolCoefPore=(KPore/mu_f)/(S11+alpha*psiPore*alpha*psiPore/M);
	double consolCoefFrac=(KFrac/mu_f)/(S22+alpha*psiFrac*alpha*psiFrac/M);
	consolCoef=min(consolCoefPore,consolCoefFrac);

	return;
}

double problemDoubleParameters::computeLeakTerm(double beta)
{
	leakTerm=beta*0.4/min(phiPore*phiPore,phiFrac*phiFrac);
	leakTerm=leakTerm*KPore/mu_f;
	if(phiPore==0 || phiFrac==0) leakTerm=0;
	
	return leakTerm;
}

void problemDoubleParameters::applyDoublePoreInitialConditions(double dy, double Ly,
	vector<vector<double>> vField, vector<vector<double>> pPoreField,
	vector<vector<double>> pFracField)
{
	double M=2*G+lambda;
	double p0=-(alpha/(M*(S11+2*S12+S22)+alpha*alpha))*sigmab;
	int v_P;
	double yValue;

	for(int i=0; i<pPoreField.size(); i++)
	{
		pPoreField[i][0]=p0;
		pFracField[i][0]=p0;
	}

	for(int i=0; i<vDisplacementFVIndex.size(); i++)
	{
		for(int j=0; j<vDisplacementFVIndex[0].size(); j++)
		{
			v_P=getVDisplacementFVPosition(i,j);
			yValue=Ly-i*dy;
			vField[v_P][0]=((sigmab+alpha*p0)/M)*yValue;
		}
	}

	vDisplacementField=vField;
	pressurePoreField=pPoreField;
	pressureFracField=pFracField;

	return;
}