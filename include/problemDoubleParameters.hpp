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
	double alpha, A11, A12, A22, consolCoef;
	double psiPore, psiFrac;
	double leakTerm;

	// Class functions
	void computeProblemParameters();
	double computeLeakTerm(double);

	// Constructor
	problemDoubleParameters(double,double,double,double,double,double,double,double,double);

	// Destructor
	~problemDoubleParameters();
};

problemDoubleParameters::problemDoubleParameters(double porosityPore, double porosityFrac,
	double fluidCompressibility, double solidCompressibility, double shearModulus,
	double lame1stParameter, double fluidViscosity, double permeabilityPore,
	double permeabilityFrac)
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

	computeProblemParameters();
}

problemDoubleParameters::~problemDoubleParameters(){}

void problemDoubleParameters::computeProblemParameters()
{
	psiPore=phiPore/(phiPore+phiFrac);
	psiFrac=phiFrac/(phiPore+phiFrac);
	alpha=1-(c_s/c_b);
	A11=phiPore*(c_f-c_s)+psiPore*alpha*(1-psiPore*alpha-phiFrac)*c_b;
	A12=psiPore*alpha*(psiFrac*alpha-phiFrac)*c_b;
	A22=phiFrac*(c_f-c_s)+psiFrac*alpha*(1-psiFrac*alpha)*c_b-phiFrac*psiPore*alpha*c_b;

	double M=2*G+lambda;
	double consolCoefPore=(KPore/mu_f)/(A11+alpha*psiPore*alpha*psiPore/M);
	double consolCoefFrac=(KFrac/mu_f)/(A22+alpha*psiFrac*alpha*psiFrac/M);
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