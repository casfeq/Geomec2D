/*
	This header is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The class defined
	here contains the functions for assembly of the coefficients matrix of the linear system which
	represents the discretized problem of poroelasticity.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.
*/

#include <iostream>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

class coefficientsAssembly
{
public:
	// Class variables
	vector<vector<double>> coefficientsMatrix;
	vector<vector<int>> boundaryConditionType;
	int Nu, Nv, NP, NPM;
	vector<vector<int>> uDisplacementFVIndex;
	vector<vector<int>> vDisplacementFVIndex;
	vector<vector<int>> pressureFVIndex;
	vector<vector<int>> uDisplacementFVCoordinates;
	vector<vector<int>> vDisplacementFVCoordinates;
	vector<vector<int>> pressureFVCoordinates;
	vector<vector<int>> horizontalFacesStatus;
	vector<vector<int>> verticalFacesStatus;
	string gridType;
	string interpScheme;
	vector<double> sparseCoefficientsRow;
	vector<double> sparseCoefficientsColumn;
	vector<double> sparseCoefficientsValue;

	// Class functions
	void resizeLinearProblem();
	int getUDisplacementFVPosition(int,int);
	int getVDisplacementFVPosition(int,int);
	int getPressureFVPosition(int,int);
	int getMacroPressureFVPosition(int,int);
	void assemblyCoefficientsMatrix(double,double,double,double,double,double,double,double,double,
		double,double);
	void assemblyXMomentum(double,double,double,double,double);
	void assemblyYMomentum(double,double,double,double,double);
	void assemblyContinuity(double,double,double,double,double,double,double,double,double);
	void addUDisplacementToXMomentum(double,double,double,double);
	void addVDisplacementToXMomentum(double,double,double,double);
	void addPressureToXMomentum(double,double);
	void addBCToXMomentum(double,double,double);
	void addUDisplacementToYMomentum(double,double,double,double);
	void addVDisplacementToYMomentum(double,double,double,double);
	void addPressureToYMomentum(double,double);
	void addBCToYMomentum(double,double,double);
	void addTransientToContinuity(double,double,double,double);
	void addFluidFlowToContinuity(double,double,double,double,double,double,double,double);
	void addDisplacementToContinuity(double,double,double,double,double,double);
	void addBCToContinuity();
	void addStaggeredVDisplacementToXMomentum(double,double,double,double);
	void addStaggeredPressureToXMomentum(double,double);
	void addStaggeredUDisplacementToYMomentum(double,double,double,double);
	void addStaggeredPressureToYMomentum(double,double);
	void addStaggeredFluidFlowToContinuity(double,double,double,double);
	void addStaggeredDisplacementToContinuity(double,double,double,double);
	void addCollocatedFluidFlowToContinuity(double,double,double,double);
	void addCDSVDisplacementToXMomentum(double,double,double,double);
	void addCDSPressureToXMomentum(double,double);
	void addCDSUDisplacementToYMomentum(double,double,double,double);
	void addCDSPressureToYMomentum(double,double);
	void addCDSDisplacementToContinuity(double,double,double,double);
	void addDirichletBCToXMomentum(double,double,double,int);
	void addDirichletBCToYMomentum(double,double,double,int);
	void addDirichletBCToContinuity(int);
	void add1DPISFluidFlowToContinuity(double,double,double,double,double,double);
	void addI2DPISFluidFlowToContinuity(double,double,double,double,double);
	void addI2DPISDisplacementToContinuity(double,double,double,double);
	void addC2DPISFluidFlowToContinuity(double,double,double,double,double,double);
	void addC2DPISDisplacementToContinuity(double,double,double,double,double,double);
	void assemblySparseMatrix(vector<vector<double>>);
	void assemblyMandelCoefficientsMatrix(double,double,double,double,double);
	void addMandelRigidMotion();
	void increaseMandelCoefficientsMatrixSize();
	void addMandelStaggeredStressToVDisplacement(double);
	void addMandelStaggeredStress(double,double,double,double,double);
	void addMandelCollocatedStressToVDisplacement(double);
	void addMandelCollocatedStress(double,double,double,double,double);
	void assemblyMacroPorosityMatrix(double,double,double,double,double,double,double,double,double,
		double,double,double,double);
	void increaseMacroPorosityCoefficientsMatrixSize();
	void addMacroPressureToXMomentum(double,double);
	void addStaggeredMacroPressureToXMomentum(double,double);
	void addCDSMacroPressureToXMomentum(double,double);
	void addMacroPressureToYMomentum(double,double);
	void addStaggeredMacroPressureToYMomentum(double,double);
	void addCDSMacroPressureToYMomentum(double,double);
	void addMacroTransientToContinuity(double,double,double,double);
	void addMacroDisplacementToContinuity(double,double,double,double,double,double);
	void addStaggeredMacroDisplacementToContinuity(double,double,double,double);
	void addCDSMacroDisplacementToContinuity(double,double,double,double);
	void addI2DPISMacroDisplacementToContinuity(double,double,double,double);
	void addMacroFluidFlowToContinuity(double,double,double,double,double,double,double,double,
		double);
	void addStaggeredMacroFluidFlowToContinuity(double,double,double,double);
	void addCollocatedMacroFluidFlowToContinuity(double,double,double,double);
	void addI2DPISMacroFluidFlowToMicroContinuity(double,double,double,double,double);
	void addI2DPISMacroFluidFlowToMacroContinuity(double,double,double,double,double);
	void addI2DPISMicroFluidFlowToMacroContinuity(double,double,double,double,double);
	void addMacroPressuretoContinuity(double,double,double);
	void addMacroBCToContinuity();
	void addMacroDirichletBCToContinuity(int);

	// Constructor
	coefficientsAssembly(vector<vector<int>>,int,int,int,vector<vector<int>>,vector<vector<int>>,
		vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,
		vector<vector<int>>,vector<vector<int>>,string,string);

	// Destructor
	~coefficientsAssembly();
};

coefficientsAssembly::coefficientsAssembly(vector<vector<int>> bcType, 
	int numberOfActiveUDisplacementFV, int numberOfActiveVDisplacementFV,
	int numberOfActivePressureFV, vector<vector<int>> idU, vector<vector<int>> idV,
	vector<vector<int>> idP, vector<vector<int>> cooU, vector<vector<int>> cooV,
	vector<vector<int>> cooP, vector<vector<int>> horFStatus, vector<vector<int>> verFStatus,
	string myGridType, string myInterpScheme)
{
	boundaryConditionType=bcType;
	Nu=numberOfActiveUDisplacementFV;
	Nv=numberOfActiveVDisplacementFV;
	NP=numberOfActivePressureFV;
	uDisplacementFVIndex=idU;
	vDisplacementFVIndex=idV;
	pressureFVIndex=idP;
	uDisplacementFVCoordinates=cooU;
	vDisplacementFVCoordinates=cooV;
	pressureFVCoordinates=cooP;
	horizontalFacesStatus=horFStatus;
	verticalFacesStatus=verFStatus;
	gridType=myGridType;
	interpScheme=myInterpScheme;

	resizeLinearProblem();
}

coefficientsAssembly::~coefficientsAssembly(){}

void coefficientsAssembly::resizeLinearProblem()
{
	int rowNo, colNo;

	// Resize coefficientsMatrix
	rowNo=NP+Nv+Nu;
	colNo=rowNo;
	coefficientsMatrix.resize(rowNo);
	for(int i=0; i<rowNo; i++)
	{
		coefficientsMatrix[i].resize(colNo);
		for(int j=0; j<colNo; j++)
		{
			coefficientsMatrix[i][j]=0;
		}
	}

	return;
}

int coefficientsAssembly::getUDisplacementFVPosition(int x, int y)
{
	int uDisplacementFVPosition;

	uDisplacementFVPosition=uDisplacementFVIndex[x][y]-1;

	return uDisplacementFVPosition;
}

int coefficientsAssembly::getVDisplacementFVPosition(int x, int y)
{
	int vDisplacementFVPosition;

	vDisplacementFVPosition=vDisplacementFVIndex[x][y]+Nu-1;

	return vDisplacementFVPosition;
}

int coefficientsAssembly::getPressureFVPosition(int x, int y)
{
	int pressureFVPosition;

	pressureFVPosition=pressureFVIndex[x][y]+Nu+Nv-1;

	return pressureFVPosition;
}

int coefficientsAssembly::getMacroPressureFVPosition(int x, int y)
{
	int pressureFVPosition;

	pressureFVPosition=pressureFVIndex[x][y]+Nu+Nv+NP-1;

	return pressureFVPosition;
}

void coefficientsAssembly::assemblyCoefficientsMatrix(double dx, double dy, double dt, double G,
	double lambda, double alpha, double K, double mu_f, double Q, double rho, double g)
{
	assemblyXMomentum(dx,dy,G,lambda,alpha);
	assemblyYMomentum(dx,dy,G,lambda,alpha);
	assemblyContinuity(dx,dy,dt,alpha,K,mu_f,Q,G,lambda);
	assemblySparseMatrix(coefficientsMatrix);

	return;
}

void coefficientsAssembly::assemblyXMomentum(double dx, double dy, double G, double lambda,
	double alpha)
{	
	addUDisplacementToXMomentum(dx,dy,G,lambda);
	addVDisplacementToXMomentum(dx,dy,G,lambda);
	addPressureToXMomentum(dy,alpha);
	addBCToXMomentum(dx,dy,G);

	return;
}

void coefficientsAssembly::assemblyYMomentum(double dx, double dy, double G, double lambda,
	double alpha)
{	
	addUDisplacementToYMomentum(dx,dy,G,lambda);
	addVDisplacementToYMomentum(dx,dy,G,lambda);
	addPressureToYMomentum(dx,alpha);
	addBCToYMomentum(dx,dy,G);

	return;
}

void coefficientsAssembly::assemblyContinuity(double dx, double dy, double dt, double alpha,
	double K, double mu_f, double Q, double G, double lambda)
{
	addTransientToContinuity(dx,dy,dt,Q);
	addFluidFlowToContinuity(dx,dy,dt,K,mu_f,alpha,G,lambda);
	addDisplacementToContinuity(dx,dy,dt,alpha,G,lambda);
	addBCToContinuity();

	return;
}

void coefficientsAssembly::addUDisplacementToXMomentum(double dx, double dy, double G,
	double lambda)
{
	int u_P, u_E, u_W, u_N, u_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(i==0) // Northern border
		{
			u_S=getUDisplacementFVPosition(i+1,j);

			if(j==0 || j==uDisplacementFVIndex[0].size()-1)	value=0.5;

			coefficientsMatrix[u_P][u_S]-=G*(dx/dy)*value;
			coefficientsMatrix[u_P][u_P]+=G*(dx/dy)*value;

			value=1;
		}
		else if(i==uDisplacementFVIndex.size()-1) // Southern border
		{
			u_N=getUDisplacementFVPosition(i-1,j);

			if(j==0 || j==uDisplacementFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[u_P][u_N]-=G*(dx/dy)*value;
			coefficientsMatrix[u_P][u_P]+=G*(dx/dy)*value;

			value=1;
		}
		else
		{
			u_N=getUDisplacementFVPosition(i-1,j);
			u_S=getUDisplacementFVPosition(i+1,j);

			if(j==0 || j==uDisplacementFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[u_P][u_N]-=G*(dx/dy)*value;
			coefficientsMatrix[u_P][u_S]-=G*(dx/dy)*value;
			coefficientsMatrix[u_P][u_P]+=2*G*(dx/dy)*value;

			value=1;
		}

		if(j==0) // Western border
		{
			u_E=getUDisplacementFVPosition(i,j+1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) if(gridType!="staggered") value=0.5;

			coefficientsMatrix[u_P][u_E]-=(2*G+lambda)*(dy/dx)*value;
			coefficientsMatrix[u_P][u_P]+=(2*G+lambda)*(dy/dx)*value;

			value=1;
		}
		else if(j==uDisplacementFVIndex[0].size()-1) // Eastern border
		{
			u_W=getUDisplacementFVPosition(i,j-1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) if(gridType!="staggered") value=0.5;

			coefficientsMatrix[u_P][u_W]-=(2*G+lambda)*(dy/dx)*value;
			coefficientsMatrix[u_P][u_P]+=(2*G+lambda)*(dy/dx)*value;

			value=1;
		}
		else
		{
			u_E=getUDisplacementFVPosition(i,j+1);
			u_W=getUDisplacementFVPosition(i,j-1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) if(gridType!="staggered") value=0.5;

			coefficientsMatrix[u_P][u_E]-=(2*G+lambda)*(dy/dx)*value;
			coefficientsMatrix[u_P][u_W]-=(2*G+lambda)*(dy/dx)*value;
			coefficientsMatrix[u_P][u_P]+=2*(2*G+lambda)*(dy/dx)*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addVDisplacementToXMomentum(double dx, double dy, double G,
	double lambda)
{
	if(gridType=="staggered") addStaggeredVDisplacementToXMomentum(dx,dy,G,lambda);
	else if(gridType=="collocated")
	{
		if(interpScheme=="CDS") addCDSVDisplacementToXMomentum(dx,dy,G,lambda);
		if(interpScheme=="1DPIS") addCDSVDisplacementToXMomentum(dx,dy,G,lambda);
		if(interpScheme=="I2DPIS") addCDSVDisplacementToXMomentum(dx,dy,G,lambda);
		if(interpScheme=="C2DPIS") addCDSVDisplacementToXMomentum(dx,dy,G,lambda);
	}

	return;
}

void coefficientsAssembly::addPressureToXMomentum(double dy, double alpha)
{
	if(gridType=="staggered") addStaggeredPressureToXMomentum(dy,alpha);
	else if(gridType=="collocated")
	{
		if(interpScheme=="CDS") addCDSPressureToXMomentum(dy,alpha);
		if(interpScheme=="1DPIS") addCDSPressureToXMomentum(dy,alpha);
		if(interpScheme=="I2DPIS") addCDSPressureToXMomentum(dy,alpha);
		if(interpScheme=="C2DPIS") addCDSPressureToXMomentum(dy,alpha);
	}

	return;
}

void coefficientsAssembly::addBCToXMomentum(double dx, double dy, double G)
{
	addDirichletBCToXMomentum(dx,dy,G,0);
	addDirichletBCToXMomentum(dx,dy,G,2);
	addDirichletBCToXMomentum(dx,dy,G,1);
	addDirichletBCToXMomentum(dx,dy,G,3);

	return;
}

void coefficientsAssembly::addUDisplacementToYMomentum(double dx, double dy, double G,
	double lambda)
{
	if(gridType=="staggered") addStaggeredUDisplacementToYMomentum(dx,dy,G,lambda);
	else if(gridType=="collocated")
	{
		if(interpScheme=="CDS") addCDSUDisplacementToYMomentum(dx,dy,G,lambda);
		if(interpScheme=="1DPIS") addCDSUDisplacementToYMomentum(dx,dy,G,lambda);
		if(interpScheme=="I2DPIS") addCDSUDisplacementToYMomentum(dx,dy,G,lambda);
		if(interpScheme=="C2DPIS") addCDSUDisplacementToYMomentum(dx,dy,G,lambda);
	}

	return;
}

void coefficientsAssembly::addVDisplacementToYMomentum(double dx, double dy, double G,
	double lambda)
{
	int v_P, v_E, v_W, v_N, v_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(i==0) // Northern border
		{
			v_S=getVDisplacementFVPosition(i+1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) if(gridType!="staggered") value=0.5;

			coefficientsMatrix[v_P][v_S]-=(2*G+lambda)*(dx/dy)*value;
			coefficientsMatrix[v_P][v_P]+=(2*G+lambda)*(dx/dy)*value;

			value=1;
		}
		else if(i==vDisplacementFVIndex.size()-1) // Southern border
		{
			v_N=getVDisplacementFVPosition(i-1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) if(gridType!="staggered") value=0.5;

			coefficientsMatrix[v_P][v_N]-=(2*G+lambda)*(dx/dy)*value;
			coefficientsMatrix[v_P][v_P]+=(2*G+lambda)*(dx/dy)*value;

			value=1;
		}
		else
		{
			v_N=getVDisplacementFVPosition(i-1,j);
			v_S=getVDisplacementFVPosition(i+1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) if(gridType!="staggered") value=0.5;

			coefficientsMatrix[v_P][v_N]-=(2*G+lambda)*(dx/dy)*value;
			coefficientsMatrix[v_P][v_S]-=(2*G+lambda)*(dx/dy)*value;
			coefficientsMatrix[v_P][v_P]+=2*(2*G+lambda)*(dx/dy)*value;

			value=1;
		}

		if(j==0) // Western border
		{
			v_E=getVDisplacementFVPosition(i,j+1);

			if(i==0 || i==vDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[v_P][v_E]-=G*(dy/dx)*value;
			coefficientsMatrix[v_P][v_P]+=G*(dy/dx)*value;

			value=1;
		}
		else if(j==vDisplacementFVIndex[0].size()-1) // Eastern border
		{
			v_W=getVDisplacementFVPosition(i,j-1);

			if(i==0 || i==vDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[v_P][v_W]-=G*(dy/dx)*value;
			coefficientsMatrix[v_P][v_P]+=G*(dy/dx)*value;

			value=1;
		}
		else
		{
			v_E=getVDisplacementFVPosition(i,j+1);
			v_W=getVDisplacementFVPosition(i,j-1);

			if(i==0 || i==vDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[v_P][v_E]-=G*(dy/dx)*value;
			coefficientsMatrix[v_P][v_W]-=G*(dy/dx)*value;
			coefficientsMatrix[v_P][v_P]+=2*G*(dy/dx)*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addPressureToYMomentum(double dx, double alpha)
{
	if(gridType=="staggered") addStaggeredPressureToYMomentum(dx,alpha);
	else if(gridType=="collocated")
	{
		if(interpScheme=="CDS") addCDSPressureToYMomentum(dx,alpha);
		if(interpScheme=="1DPIS") addCDSPressureToYMomentum(dx,alpha);
		if(interpScheme=="I2DPIS") addCDSPressureToYMomentum(dx,alpha);
		if(interpScheme=="C2DPIS") addCDSPressureToYMomentum(dx,alpha);
	}

	return;
}

void coefficientsAssembly::addBCToYMomentum(double dx, double dy, double G)
{
	addDirichletBCToYMomentum(dx,dy,G,1);
	addDirichletBCToYMomentum(dx,dy,G,3);
	addDirichletBCToYMomentum(dx,dy,G,0);
	addDirichletBCToYMomentum(dx,dy,G,2);

	return;
}

void coefficientsAssembly::addTransientToContinuity(double dx, double dy, double dt, double Q)
{	
	int P_P;
	int FVCounter;
	int i, j;
	double Mp=(1/Q)*(dx*dy/dt);
	int borderCounter=0;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(gridType=="staggered")
		{
			coefficientsMatrix[P_P][P_P]+=Mp;
		}
		else if(gridType=="collocated")
		{
			if(i==0 || i==pressureFVIndex.size()-1) borderCounter++;
			if(j==0 || j==pressureFVIndex[0].size()-1) borderCounter++;

			coefficientsMatrix[P_P][P_P]+=Mp/pow(2,borderCounter);

			borderCounter=0;
		}
	}
}

void coefficientsAssembly::addFluidFlowToContinuity(double dx, double dy, double dt, double K,
	double mu_f, double alpha, double G, double lambda)
{
	if(gridType=="staggered") addStaggeredFluidFlowToContinuity(dx,dy,K,mu_f);
	else if(gridType=="collocated") 
	{
		addCollocatedFluidFlowToContinuity(dx,dy,K,mu_f);

		if(interpScheme=="1DPIS") add1DPISFluidFlowToContinuity(dx,dy,dt,alpha,G,lambda);
		if(interpScheme=="I2DPIS") addI2DPISFluidFlowToContinuity(dx,dy,dt,alpha,G);
		if(interpScheme=="C2DPIS") addC2DPISFluidFlowToContinuity(dx,dy,dt,alpha,G,lambda);
	}

	return;
}

void coefficientsAssembly::addDisplacementToContinuity(double dx, double dy, double dt,
	double alpha, double G, double lambda)
{
	if(gridType=="staggered") addStaggeredDisplacementToContinuity(dx,dy,dt,alpha);
	else if(gridType=="collocated")
	{
		if(interpScheme=="CDS") addCDSDisplacementToContinuity(dx,dy,dt,alpha);
		if(interpScheme=="1DPIS") addCDSDisplacementToContinuity(dx,dy,dt,alpha);
		if(interpScheme=="I2DPIS") addI2DPISDisplacementToContinuity(dx,dy,dt,alpha);
		if(interpScheme=="C2DPIS") addC2DPISDisplacementToContinuity(dx,dy,dt,alpha,G,lambda);
	}
	
	return;
}

void coefficientsAssembly::addBCToContinuity()
{
	if(gridType=="collocated")
		for(int counter=0; counter<4; counter++)
			addDirichletBCToContinuity(counter);

	return;
}

void coefficientsAssembly::addStaggeredVDisplacementToXMomentum(double dx, double dy, double G,
	double lambda)
{
	int u_P;
	int v_P, v_W, v_S, v_SW;
	int FVCounter;
	int i, j;

	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(j==0) // Western border
		{
			v_P=getVDisplacementFVPosition(i,j);
			v_S=getVDisplacementFVPosition(i+1,j);

			coefficientsMatrix[u_P][v_P]-=lambda;
			coefficientsMatrix[u_P][v_S]-=-lambda;
		}
		else if(j==uDisplacementFVIndex[0].size()-1) // Eastern border
		{
			v_W=getVDisplacementFVPosition(i,j-1);
			v_SW=getVDisplacementFVPosition(i+1,j-1);

			coefficientsMatrix[u_P][v_W]-=-lambda;
			coefficientsMatrix[u_P][v_SW]-=lambda;
		}
		else
		{
			if(i==0) // Northern border
			{
				v_P=getVDisplacementFVPosition(i,j);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[u_P][v_P]-=lambda;
				coefficientsMatrix[u_P][v_W]-=-lambda;
				coefficientsMatrix[u_P][v_S]-=-G-lambda;
				coefficientsMatrix[u_P][v_SW]-=G+lambda;
			}
			else if(i==uDisplacementFVIndex.size()-1) // Southern border
			{
				v_P=getVDisplacementFVPosition(i,j);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[u_P][v_P]-=G+lambda;
				coefficientsMatrix[u_P][v_W]-=-G-lambda;
				coefficientsMatrix[u_P][v_S]-=-lambda;
				coefficientsMatrix[u_P][v_SW]-=lambda;
			}	
			else
			{
				v_P=getVDisplacementFVPosition(i,j);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[u_P][v_P]-=G+lambda;
				coefficientsMatrix[u_P][v_W]-=-G-lambda;
				coefficientsMatrix[u_P][v_S]-=-G-lambda;
				coefficientsMatrix[u_P][v_SW]-=G+lambda;
			}
		}
	}

	return;
}

void coefficientsAssembly::addStaggeredPressureToXMomentum(double dy, double alpha)
{
	int u_P;
	int P_P, P_W;
	int FVCounter;
	int i, j;

	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(j==0) // FV on the western border
		{
			P_P=getPressureFVPosition(i,j);

			coefficientsMatrix[u_P][P_P]-=-alpha*dy;
		}
		else if(j==uDisplacementFVIndex[0].size()-1) // FV on the eastern border
		{
			P_W=getPressureFVPosition(i,j-1);

			coefficientsMatrix[u_P][P_W]-=alpha*dy;
		}
		else // FV not on the western or eastern border
		{
			P_P=getPressureFVPosition(i,j);
			P_W=getPressureFVPosition(i,j-1);

			coefficientsMatrix[u_P][P_P]-=-alpha*dy;
			coefficientsMatrix[u_P][P_W]-=alpha*dy;
		}
	}

	return;
}

void coefficientsAssembly::addStaggeredUDisplacementToYMomentum(double dx, double dy, double G,
	double lambda)
{
	int u_P, u_E, u_N, u_NE;
	int v_P;
	int FVCounter;
	int i, j;

	for(FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(i==0) // Northern border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_E=getUDisplacementFVPosition(i,j+1);

			coefficientsMatrix[v_P][u_P]-=lambda;
			coefficientsMatrix[v_P][u_E]-=-lambda;
		}
		else if(i==vDisplacementFVIndex.size()-1) // Southern border
		{
			u_N=getUDisplacementFVPosition(i-1,j);
			u_NE=getUDisplacementFVPosition(i-1,j+1);

			coefficientsMatrix[v_P][u_N]-=-lambda;
			coefficientsMatrix[v_P][u_NE]-=lambda;
		}
		else
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_NE=getUDisplacementFVPosition(i-1,j+1);

				coefficientsMatrix[v_P][u_P]-=lambda;
				coefficientsMatrix[v_P][u_E]-=-G-lambda;
				coefficientsMatrix[v_P][u_N]-=-lambda;
				coefficientsMatrix[v_P][u_NE]-=G+lambda;
			}
			else if(j==vDisplacementFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_NE=getUDisplacementFVPosition(i-1,j+1);

				coefficientsMatrix[v_P][u_P]-=G+lambda;
				coefficientsMatrix[v_P][u_E]-=-lambda;
				coefficientsMatrix[v_P][u_N]-=-G-lambda;
				coefficientsMatrix[v_P][u_NE]-=lambda;
			}
			else
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_NE=getUDisplacementFVPosition(i-1,j+1);

				coefficientsMatrix[v_P][u_P]-=G+lambda;
				coefficientsMatrix[v_P][u_E]-=-G-lambda;
				coefficientsMatrix[v_P][u_N]-=-G-lambda;
				coefficientsMatrix[v_P][u_NE]-=G+lambda;
			}
		}
	}

	return;
}

void coefficientsAssembly::addStaggeredPressureToYMomentum(double dx, double alpha)
{
	int v_P;
	int P_P, P_N;
	int FVCounter;
	int i, j;
	
	for(FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(i==0) // FV on the northern border
		{
			P_P=getPressureFVPosition(i,j);

			coefficientsMatrix[v_P][P_P]-=alpha*dx;
		}
		else if(i==vDisplacementFVIndex.size()-1) // FV on the southern border
		{
			P_N=getPressureFVPosition(i-1,j);

			coefficientsMatrix[v_P][P_N]-=-alpha*dx;
		}
		else // FV not on the northern or southern border
		{
			P_P=getPressureFVPosition(i,j);
			P_N=getPressureFVPosition(i-1,j);

			coefficientsMatrix[v_P][P_P]-=alpha*dx;
			coefficientsMatrix[v_P][P_N]-=-alpha*dx;
		}
	}

	return;
}

void coefficientsAssembly::addStaggeredFluidFlowToContinuity(double dx, double dy, double K,
	double mu_f)
{
	int P_P, P_E, P_W, P_N, P_S;
	int FVCounter;
	int i, j;
	int bcType;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			P_S=getPressureFVPosition(i+1,j);

			coefficientsMatrix[P_P][P_S]-=(K/mu_f)*(dx/dy);
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dx/dy);

			bcType=boundaryConditionType[0][2];
			if(bcType==1) coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dx/dy);
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			P_N=getPressureFVPosition(i-1,j);

			coefficientsMatrix[P_P][P_N]-=(K/mu_f)*(dx/dy);
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dx/dy);

			bcType=boundaryConditionType[2][2];
			if(bcType==1) coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dx/dy);
		}
		else
		{
			P_N=getPressureFVPosition(i-1,j);
			P_S=getPressureFVPosition(i+1,j);

			coefficientsMatrix[P_P][P_N]-=(K/mu_f)*(dx/dy);
			coefficientsMatrix[P_P][P_S]-=(K/mu_f)*(dx/dy);
			coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dx/dy);
		}

		if(j==0) // Western border
		{
			P_E=getPressureFVPosition(i,j+1);

			coefficientsMatrix[P_P][P_E]-=(K/mu_f)*(dy/dx);
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dy/dx);

			bcType=boundaryConditionType[1][2];
			if(bcType==1) coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dy/dx);
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			P_W=getPressureFVPosition(i,j-1);

			coefficientsMatrix[P_P][P_W]-=(K/mu_f)*(dy/dx);
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dy/dx);

			bcType=boundaryConditionType[3][2];
			if(bcType==1) coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dy/dx);
		}
		else
		{
			P_E=getPressureFVPosition(i,j+1);
			P_W=getPressureFVPosition(i,j-1);

			coefficientsMatrix[P_P][P_E]-=(K/mu_f)*(dy/dx);
			coefficientsMatrix[P_P][P_W]-=(K/mu_f)*(dy/dx);
			coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dy/dx);
		}
	}

	return;
}

void coefficientsAssembly::addStaggeredDisplacementToContinuity(double dx, double dy, double dt,
	double alpha)
{
	int u_P, u_E;
	int v_P, v_S;
	int P_P;
	int FVCounter;
	int i, j;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);
		u_E=getUDisplacementFVPosition(i,j+1);
		v_P=getVDisplacementFVPosition(i,j);
		v_S=getVDisplacementFVPosition(i+1,j);
		P_P=getPressureFVPosition(i,j);

		coefficientsMatrix[P_P][u_P]-=alpha*(dy/dt);
		coefficientsMatrix[P_P][v_P]-=-alpha*(dx/dt);
		coefficientsMatrix[P_P][u_E]-=-alpha*(dy/dt);
		coefficientsMatrix[P_P][v_S]-=alpha*(dx/dt);
	}

	return;
}

void coefficientsAssembly::addCollocatedFluidFlowToContinuity(double dx, double dy, double K,
	double mu_f)
{
	int P_P, P_E, P_W, P_N, P_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			P_S=getPressureFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_S]-=(K/mu_f)*(dx/dy)*value;
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dx/dy)*value;

			value=1;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			P_N=getPressureFVPosition(i-1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_N]-=(K/mu_f)*(dx/dy)*value;
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dx/dy)*value;

			value=1;
		}
		else
		{
			P_N=getPressureFVPosition(i-1,j);
			P_S=getPressureFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_N]-=(K/mu_f)*(dx/dy)*value;
			coefficientsMatrix[P_P][P_S]-=(K/mu_f)*(dx/dy)*value;
			coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dx/dy)*value;

			value=1;
		}

		if(j==0) // Western border
		{
			P_E=getPressureFVPosition(i,j+1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_E]-=(K/mu_f)*(dy/dx)*value;
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dy/dx)*value;

			value=1;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			P_W=getPressureFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_W]-=(K/mu_f)*(dy/dx)*value;
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dy/dx)*value;

			value=1;
		}
		else
		{
			P_E=getPressureFVPosition(i,j+1);
			P_W=getPressureFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_E]-=(K/mu_f)*(dy/dx)*value;
			coefficientsMatrix[P_P][P_W]-=(K/mu_f)*(dy/dx)*value;
			coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dy/dx)*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addCDSVDisplacementToXMomentum(double dx, double dy, double G,
	double lambda)
{
	int u_P;
	int v_P, v_E, v_W, v_N, v_S, v_NE, v_NW, v_SE, v_SW;
	int FVCounter;
	int i, j;
	int nodesCounter=0;

	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				v_P=getVDisplacementFVPosition(i,j);
				v_E=getVDisplacementFVPosition(i,j+1);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_SE=getVDisplacementFVPosition(i+1,j+1);

				coefficientsMatrix[u_P][v_P]-=0.25*(G+lambda);
				coefficientsMatrix[u_P][v_E]-=-0.25*(G-lambda);
				coefficientsMatrix[u_P][v_S]-=0.25*(G-lambda);
				coefficientsMatrix[u_P][v_SE]-=-0.25*(G+lambda);
			}
			else if(j==uDisplacementFVIndex[0].size()-1) // Eastern border
			{
				v_P=getVDisplacementFVPosition(i,j);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[u_P][v_P]-=-0.25*(G+lambda);
				coefficientsMatrix[u_P][v_W]-=0.25*(G-lambda);
				coefficientsMatrix[u_P][v_S]-=-0.25*(G-lambda);
				coefficientsMatrix[u_P][v_SW]-=0.25*(G+lambda);
			}
			else
			{
				v_E=getVDisplacementFVPosition(i,j+1);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_SE=getVDisplacementFVPosition(i+1,j+1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[u_P][v_E]-=-0.25*(G-lambda);
				coefficientsMatrix[u_P][v_W]-=0.25*(G-lambda);
				coefficientsMatrix[u_P][v_SE]-=-0.25*(G+lambda);
				coefficientsMatrix[u_P][v_SW]-=0.25*(G+lambda);
			}
		}
		else if(i==uDisplacementFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				v_P=getVDisplacementFVPosition(i,j);
				v_E=getVDisplacementFVPosition(i,j+1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);

				coefficientsMatrix[u_P][v_P]-=-0.25*(G+lambda);
				coefficientsMatrix[u_P][v_E]-=0.25*(G-lambda);
				coefficientsMatrix[u_P][v_N]-=-0.25*(G-lambda);
				coefficientsMatrix[u_P][v_NE]-=0.25*(G+lambda);
			}
			else if(j==uDisplacementFVIndex[0].size()-1) // Eastern border
			{
				v_P=getVDisplacementFVPosition(i,j);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_NW=getVDisplacementFVPosition(i-1,j-1);

				coefficientsMatrix[u_P][v_P]-=0.25*(G+lambda);
				coefficientsMatrix[u_P][v_W]-=-0.25*(G-lambda);
				coefficientsMatrix[u_P][v_N]-=0.25*(G-lambda);
				coefficientsMatrix[u_P][v_NW]-=-0.25*(G+lambda);
			}
			else
			{
				v_E=getVDisplacementFVPosition(i,j+1);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_NW=getVDisplacementFVPosition(i-1,j-1);

				coefficientsMatrix[u_P][v_E]-=0.25*(G-lambda);
				coefficientsMatrix[u_P][v_W]-=-0.25*(G-lambda);
				coefficientsMatrix[u_P][v_NE]-=0.25*(G+lambda);
				coefficientsMatrix[u_P][v_NW]-=-0.25*(G+lambda);
			}
		}
		else
		{
			if(j==0) // Western border
			{
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_SE=getVDisplacementFVPosition(i+1,j+1);

				coefficientsMatrix[u_P][v_N]-=-0.25*(G-lambda);
				coefficientsMatrix[u_P][v_S]-=0.25*(G-lambda);
				coefficientsMatrix[u_P][v_NE]-=0.25*(G+lambda);
				coefficientsMatrix[u_P][v_SE]-=-0.25*(G+lambda);
			}
			else if(j==uDisplacementFVIndex[0].size()-1) // Eastern border
			{
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_NW=getVDisplacementFVPosition(i-1,j-1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[u_P][v_N]-=0.25*(G-lambda);
				coefficientsMatrix[u_P][v_S]-=-0.25*(G-lambda);
				coefficientsMatrix[u_P][v_NW]-=-0.25*(G+lambda);
				coefficientsMatrix[u_P][v_SW]-=0.25*(G+lambda);
			}
			else
			{
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_NW=getVDisplacementFVPosition(i-1,j-1);
				v_SE=getVDisplacementFVPosition(i+1,j+1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[u_P][v_NE]-=0.25*(G+lambda);
				coefficientsMatrix[u_P][v_NW]-=-0.25*(G+lambda);
				coefficientsMatrix[u_P][v_SE]-=-0.25*(G+lambda);
				coefficientsMatrix[u_P][v_SW]-=0.25*(G+lambda);
			}
		}
	}

	return;
}

void coefficientsAssembly::addCDSPressureToXMomentum(double dy, double alpha)
{
	int u_P;
	int P_P, P_E, P_W;
	int FVCounter;
	int i, j;
	double value=1;
	
	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{	
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(j==0) // FV on the western border
		{
			P_P=getPressureFVPosition(i,j);
			P_E=getPressureFVPosition(i,j+1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[u_P][P_P]-=-0.5*alpha*dy*value;
			coefficientsMatrix[u_P][P_E]-=-0.5*alpha*dy*value;

			value=1;
		}
		else if(j==uDisplacementFVIndex[0].size()-1) // FV on the eastern border
		{
			P_P=getPressureFVPosition(i,j);
			P_W=getPressureFVPosition(i,j-1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[u_P][P_P]-=0.5*alpha*dy*value;
			coefficientsMatrix[u_P][P_W]-=0.5*alpha*dy*value;

			value=1;
		}
		else // FV not on the western or eastern border
		{
			P_E=getPressureFVPosition(i,j+1);
			P_W=getPressureFVPosition(i,j-1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[u_P][P_E]-=-0.5*alpha*dy*value;
			coefficientsMatrix[u_P][P_W]-=0.5*alpha*dy*value;

			value=1;
		}
	}
	
	return;
}

void coefficientsAssembly::addCDSUDisplacementToYMomentum(double dx, double dy, double G,
	double lambda)
{
	int v_P;
	int u_P, u_E, u_W, u_N, u_S, u_NE, u_NW, u_SE, u_SW;
	int FVCounter;
	int i, j;
	int nodesCounter=0;

	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_SE=getUDisplacementFVPosition(i+1,j+1);

				coefficientsMatrix[v_P][u_P]-=0.25*(G+lambda);
				coefficientsMatrix[v_P][u_E]-=0.25*(G-lambda);
				coefficientsMatrix[v_P][u_S]-=-0.25*(G-lambda);
				coefficientsMatrix[v_P][u_SE]-=-0.25*(G+lambda);
			}
			else if(j==vDisplacementFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_SW=getUDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[v_P][u_P]-=-0.25*(G+lambda);
				coefficientsMatrix[v_P][u_W]-=-0.25*(G-lambda);
				coefficientsMatrix[v_P][u_S]-=0.25*(G-lambda);
				coefficientsMatrix[v_P][u_SW]-=0.25*(G+lambda);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[v_P][u_E]-=0.25*(G-lambda);
				coefficientsMatrix[v_P][u_W]-=-0.25*(G-lambda);
				coefficientsMatrix[v_P][u_SE]-=-0.25*(G+lambda);
				coefficientsMatrix[v_P][u_SW]-=0.25*(G+lambda);
			}
		}
		else if(i==vDisplacementFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_NE=getUDisplacementFVPosition(i-1,j+1);

				coefficientsMatrix[v_P][u_P]-=-0.25*(G+lambda);
				coefficientsMatrix[v_P][u_E]-=-0.25*(G-lambda);
				coefficientsMatrix[v_P][u_N]-=0.25*(G-lambda);
				coefficientsMatrix[v_P][u_NE]-=0.25*(G+lambda);
			}
			else if(j==vDisplacementFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_NW=getUDisplacementFVPosition(i-1,j-1);

				coefficientsMatrix[v_P][u_P]-=0.25*(G+lambda);
				coefficientsMatrix[v_P][u_W]-=0.25*(G-lambda);
				coefficientsMatrix[v_P][u_N]-=-0.25*(G-lambda);
				coefficientsMatrix[v_P][u_NW]-=-0.25*(G+lambda);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_NW=getUDisplacementFVPosition(i-1,j-1);

				coefficientsMatrix[v_P][u_E]-=-0.25*(G-lambda);
				coefficientsMatrix[v_P][u_W]-=0.25*(G-lambda);
				coefficientsMatrix[v_P][u_NE]-=0.25*(G+lambda);
				coefficientsMatrix[v_P][u_NW]-=-0.25*(G+lambda);
			}
		}
		else
		{
			if(j==0) // Western border
			{
				u_N=getUDisplacementFVPosition(i-1,j);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);

				coefficientsMatrix[v_P][u_N]-=0.25*(G-lambda);
				coefficientsMatrix[v_P][u_S]-=-0.25*(G-lambda);
				coefficientsMatrix[v_P][u_NE]-=0.25*(G+lambda);
				coefficientsMatrix[v_P][u_SE]-=-0.25*(G+lambda);
			}
			else if(j==vDisplacementFVIndex[0].size()-1) // Eastern border
			{
				u_N=getUDisplacementFVPosition(i-1,j);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[v_P][u_N]-=-0.25*(G-lambda);
				coefficientsMatrix[v_P][u_S]-=0.25*(G-lambda);
				coefficientsMatrix[v_P][u_NW]-=-0.25*(G+lambda);
				coefficientsMatrix[v_P][u_SW]-=0.25*(G+lambda);
			}
			else
			{
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[v_P][u_NE]-=0.25*(G+lambda);
				coefficientsMatrix[v_P][u_NW]-=-0.25*(G+lambda);
				coefficientsMatrix[v_P][u_SE]-=-0.25*(G+lambda);
				coefficientsMatrix[v_P][u_SW]-=0.25*(G+lambda);
			}
		}
	}

	return;
}

void coefficientsAssembly::addCDSPressureToYMomentum(double dx, double alpha)
{
	int v_P;
	int P_P, P_N, P_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(i==0) // FV on the northern border
		{
			P_P=getPressureFVPosition(i,j);
			P_S=getPressureFVPosition(i+1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[v_P][P_P]-=0.5*alpha*dx*value;
			coefficientsMatrix[v_P][P_S]-=0.5*alpha*dx*value;

			value=1;
		}
		else if(i==vDisplacementFVIndex.size()-1) // FV on the southern border
		{
			P_P=getPressureFVPosition(i,j);
			P_N=getPressureFVPosition(i-1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[v_P][P_P]-=-0.5*alpha*dx*value;
			coefficientsMatrix[v_P][P_N]-=-0.5*alpha*dx*value;

			value=1;
		}
		else // FV not on the northern or southern border
		{
			P_N=getPressureFVPosition(i-1,j);
			P_S=getPressureFVPosition(i+1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[v_P][P_N]-=-0.5*alpha*dx*value;
			coefficientsMatrix[v_P][P_S]-=0.5*alpha*dx*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addCDSDisplacementToContinuity(double dx, double dy, double dt,
	double alpha)
{
	int P_P;
	int u_P, u_E, u_W;
	int v_P, v_N, v_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			v_P=getVDisplacementFVPosition(i,j);
			v_S=getVDisplacementFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][v_P]-=-0.5*alpha*(dx/dt)*value;
			coefficientsMatrix[P_P][v_S]-=0.5*alpha*(dx/dt)*value;

			value=1;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			v_P=getVDisplacementFVPosition(i,j);
			v_N=getVDisplacementFVPosition(i-1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][v_P]-=0.5*alpha*(dx/dt)*value;
			coefficientsMatrix[P_P][v_N]-=-0.5*alpha*(dx/dt)*value;

			value=1;
		}
		else
		{
			v_N=getVDisplacementFVPosition(i-1,j);
			v_S=getVDisplacementFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][v_N]-=-0.5*alpha*(dx/dt)*value;
			coefficientsMatrix[P_P][v_S]-=0.5*alpha*(dx/dt)*value;

			value=1;
		}

		if(j==0) // Western border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_E=getUDisplacementFVPosition(i,j+1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][u_P]-=0.5*alpha*(dy/dt)*value;
			coefficientsMatrix[P_P][u_E]-=-0.5*alpha*(dy/dt)*value;

			value=1;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_W=getUDisplacementFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][u_P]-=-0.5*alpha*(dy/dt)*value;
			coefficientsMatrix[P_P][u_W]-=0.5*alpha*(dy/dt)*value;

			value=1;
		}
		else
		{
			u_E=getUDisplacementFVPosition(i,j+1);
			u_W=getUDisplacementFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][u_E]-=-0.5*alpha*(dy/dt)*value;
			coefficientsMatrix[P_P][u_W]-=0.5*alpha*(dy/dt)*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addDirichletBCToXMomentum(double dx, double dy, double G, int counter)
{
	int u_P;
	int bcType;
	int i, j;

	bcType=boundaryConditionType[counter][0];
	
	if(bcType==1)
	{
		switch(counter)
		{
			case 0:
				for(j=0; j<uDisplacementFVIndex[0].size(); j++)
				{
					u_P=getUDisplacementFVPosition(0,j);

					if(gridType=="staggered")
					{
						coefficientsMatrix[u_P][u_P]+=2*G*(dx/dy);
					}
					else if(gridType=="collocated")
					{
						coefficientsMatrix[u_P].assign(Nu+Nv+NP+NPM,0);
						coefficientsMatrix[u_P][u_P]+=1;
					}
				}
				break;

			case 1:
				for(i=0; i<uDisplacementFVIndex.size(); i++)
				{
					u_P=getUDisplacementFVPosition(i,0);

					coefficientsMatrix[u_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[u_P][u_P]+=1;
				}
				break;

			case 2:
				for(j=0; j<uDisplacementFVIndex[0].size(); j++)
				{
					u_P=getUDisplacementFVPosition(uDisplacementFVIndex.size()-1,j);

					if(gridType=="staggered")
					{
						coefficientsMatrix[u_P][u_P]+=2*G*(dx/dy);
					}					
					else if(gridType=="collocated")
					{
						coefficientsMatrix[u_P].assign(Nu+Nv+NP+NPM,0);
						coefficientsMatrix[u_P][u_P]+=1;
					}
				}
				break;

			case 3:
				for(i=0; i<uDisplacementFVIndex.size(); i++)
				{
					u_P=getUDisplacementFVPosition(i,uDisplacementFVIndex[0].size()-1);

					coefficientsMatrix[u_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[u_P][u_P]+=1;
				}
				break;
		}
	}
	/**/
	return;
}

void coefficientsAssembly::addDirichletBCToYMomentum(double dx, double dy, double G, int counter)
{
	int v_P;
	int bcType;
	int i, j;

	bcType=boundaryConditionType[counter][1];

	if(bcType==1)
	{
		switch(counter)
		{
			case 0:
				for(j=0; j<vDisplacementFVIndex[0].size(); j++)
				{
					v_P=getVDisplacementFVPosition(0,j);

					coefficientsMatrix[v_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[v_P][v_P]+=1;
				}
				break;

			case 1:
				for(i=0; i<vDisplacementFVIndex.size(); i++)
				{
					v_P=getVDisplacementFVPosition(i,0);

					if(gridType=="staggered")
					{
						coefficientsMatrix[v_P][v_P]+=2*G*(dy/dx);
					}
					else if(gridType=="collocated")
					{
						coefficientsMatrix[v_P].assign(Nu+Nv+NP+NPM,0);
						coefficientsMatrix[v_P][v_P]+=1;	
					}
				}
				break;

			case 2:
				for(j=0; j<vDisplacementFVIndex[0].size(); j++)
				{
					v_P=getVDisplacementFVPosition(vDisplacementFVIndex.size()-1,j);

					coefficientsMatrix[v_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[v_P][v_P]+=1;
				}
				break;

			case 3:
				for(i=0; i<vDisplacementFVIndex.size(); i++)
				{
					v_P=getVDisplacementFVPosition(i,vDisplacementFVIndex[0].size()-1);

					if(gridType=="staggered")
					{
						coefficientsMatrix[v_P][v_P]+=2*G*(dy/dx);
					}
					else if(gridType=="collocated")
					{
						coefficientsMatrix[v_P].assign(Nu+Nv+NP+NPM,0);
						coefficientsMatrix[v_P][v_P]+=1;
					}
				}
				break;
		}
	}

	return;
}

void coefficientsAssembly::addDirichletBCToContinuity(int counter)
{
	int P_P;
	int bcType;
	int i, j;

	bcType=boundaryConditionType[counter][2];

	if(bcType==1)
	{
		switch(counter)
		{
			case 0:
				for(j=0; j<pressureFVIndex[0].size(); j++)
				{
					P_P=getPressureFVPosition(0,j);

					coefficientsMatrix[P_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[P_P][P_P]+=1;
				}
				break;

			case 1:
				for(i=0; i<pressureFVIndex.size(); i++)
				{
					P_P=getPressureFVPosition(i,0);

					coefficientsMatrix[P_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[P_P][P_P]+=1;
				}
				break;

			case 2:
				for(j=0; j<pressureFVIndex[0].size(); j++)
				{
					P_P=getPressureFVPosition(pressureFVIndex.size()-1,j);

					coefficientsMatrix[P_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[P_P][P_P]+=1;
				}
				break;

			case 3:
				for(i=0; i<pressureFVIndex.size(); i++)
				{
					P_P=getPressureFVPosition(i,pressureFVIndex[0].size()-1);

					coefficientsMatrix[P_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[P_P][P_P]+=1;
				}
				break;
		}
	}

	return;
}

void coefficientsAssembly::add1DPISFluidFlowToContinuity(double dx, double dy, double dt,
	double alpha, double G, double lambda)
{
	int P_P, P_E, P_W, P_N, P_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			P_S=getPressureFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*value;
			coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*value;

			value=1;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			P_N=getPressureFVPosition(i-1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*value;
			coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*value;

			value=1;
		}
		else
		{
			P_N=getPressureFVPosition(i-1,j);
			P_S=getPressureFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*value;
			coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*value;
			coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*value;

			value=1;
		}

		if(j==0) // Western border
		{
			P_E=getPressureFVPosition(i,j+1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*value;
			coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*value;

			value=1;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			P_W=getPressureFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*value;
			coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*value;

			value=1;
		}
		else
		{
			P_E=getPressureFVPosition(i,j+1);
			P_W=getPressureFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*value;
			coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*value;
			coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addI2DPISFluidFlowToContinuity(double dx, double dy, double dt,
	double alpha, double G)
{
	int P_P, P_E, P_W, P_N, P_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_S=getPressureFVPosition(i+1,j);
	
				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;		
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);				
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
			}
		}
	}
	
	return;
}

void coefficientsAssembly::addI2DPISDisplacementToContinuity(double dx, double dy, double dt,
	double alpha)
{
	int P_P;
	int u_P, u_E, u_W, u_N, u_S, u_NE, u_NW, u_SE, u_SW;
	int v_P, v_E, v_W, v_N, v_S, v_NE, v_NW, v_SE, v_SW;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{	
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_E=getVDisplacementFVPosition(i,j+1);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_SE=getVDisplacementFVPosition(i+1,j+1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_P]-=-alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_E]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_W]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_SE]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_SW]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				coefficientsMatrix[P_P][u_P]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				coefficientsMatrix[P_P][u_P]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_E=getVDisplacementFVPosition(i,j+1);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_NW=getVDisplacementFVPosition(i-1,j-1);

				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=-(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_P]-=alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_E]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_W]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_NE]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_NW]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
			}
		}
		else
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=-(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_P]-=alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_N]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_NE]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_S]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_SE]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_P]-=-alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_N]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_S]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_NW]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_SW]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_NW=getVDisplacementFVPosition(i-1,j-1);
				v_SE=getVDisplacementFVPosition(i+1,j+1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[P_P][u_E]-=-(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_W]-=(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_NE]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_NW]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_SE]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_SW]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_N]-=-(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_S]-=(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_NE]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_NW]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_SE]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_SW]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
			}
		}
	}

	return;
}

void coefficientsAssembly::addC2DPISFluidFlowToContinuity(double dx, double dy, double dt,
	double alpha, double G, double lambda)
{
	int P_P, P_E, P_W, P_N, P_S;
	int FVCounter;
	int i, j;
	double value=1;
	double Hx, Hy;

	Hx=alpha/(2*(G*dx*dx+(2*G+lambda)*dy*dy)*dt);
	Hy=alpha/(2*((2*G+lambda)*dx*dx+G*dy*dy)*dt);

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_S=getPressureFVPosition(i+1,j);
	
				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_S]-=0.25*alpha*dx*dx*dx*dy*Hy;
				coefficientsMatrix[P_P][P_P]+=0.25*alpha*dx*dx*dx*dy*Hy;
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;		
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][P_N]-=0.25*alpha*dx*dx*dx*dy*Hy;
				coefficientsMatrix[P_P][P_P]+=0.25*alpha*dx*dx*dx*dy*Hy;
				coefficientsMatrix[P_P][P_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_E]-=0.25*alpha*dx*dy*dy*dy*Hx;
				coefficientsMatrix[P_P][P_P]+=0.25*alpha*dx*dy*dy*dy*Hx;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_W]-=0.25*alpha*dx*dy*dy*dy*Hx;
				coefficientsMatrix[P_P][P_P]+=0.25*alpha*dx*dy*dy*dy*Hx;				
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][P_N]-=0.25*alpha*dx*dx*dx*dy*Hy;
				coefficientsMatrix[P_P][P_S]-=0.25*alpha*dx*dx*dx*dy*Hy;
				coefficientsMatrix[P_P][P_P]+=2*0.25*alpha*dx*dx*dx*dy*Hy;
				coefficientsMatrix[P_P][P_E]-=0.25*alpha*dx*dy*dy*dy*Hx;
				coefficientsMatrix[P_P][P_W]-=0.25*alpha*dx*dy*dy*dy*Hx;
				coefficientsMatrix[P_P][P_P]+=2*0.25*alpha*dx*dy*dy*dy*Hx;
			}
		}
	}
	
	return;
}

void coefficientsAssembly::addC2DPISDisplacementToContinuity(double dx, double dy, double dt,
	double alpha, double G, double lambda)
{
	int P_P;
	int u_P, u_E, u_W, u_N, u_S, u_NE, u_NW, u_SE, u_SW;
	int v_P, v_E, v_W, v_N, v_S, v_NE, v_NW, v_SE, v_SW;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{	
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_E=getVDisplacementFVPosition(i,j+1);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_SE=getVDisplacementFVPosition(i+1,j+1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][u_W]-=((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][u_SE]-=-((G+lambda)*dx*dy)
					/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)*alpha*(dx/dt);
				coefficientsMatrix[P_P][u_SW]-=-((G+lambda)*dx*dy)
					/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_P]-=((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_P]-=-alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_E]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_W]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_SE]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_SW]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				coefficientsMatrix[P_P][u_P]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				coefficientsMatrix[P_P][u_P]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_E=getVDisplacementFVPosition(i,j+1);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_NW=getVDisplacementFVPosition(i-1,j-1);

				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][u_W]-=((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][u_NE]-=-((G+lambda)*dx*dy)
					/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)*alpha*(dx/dt);
				coefficientsMatrix[P_P][u_NW]-=-((G+lambda)*dx*dy)
					/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_P]-=((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_P]-=-alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_E]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_W]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_NE]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_NW]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
			}
		}
		else
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_SE=getVDisplacementFVPosition(i+1,j+1);

				coefficientsMatrix[P_P][u_P]-=-((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_P]-=-alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_N]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_NE]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_S]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_SE]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_S]-=((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_NE]-=-((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_SE]-=-((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_NW=getVDisplacementFVPosition(i-1,j-1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[P_P][u_P]-=((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_P]-=-alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_N]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_S]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_NW]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_SW]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_S]-=((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_NW]-=-((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_SW]-=-((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_NW=getVDisplacementFVPosition(i-1,j-1);
				v_SE=getVDisplacementFVPosition(i+1,j+1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[P_P][u_E]-=-((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_NE]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_NW]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_SE]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][u_SW]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt);
				coefficientsMatrix[P_P][v_N]-=-((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_NE]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_NW]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_SE]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_SW]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt);
			}
		}
	}

	return;
}

void coefficientsAssembly::assemblySparseMatrix(vector<vector<double>> myCoefficientsMatrix)
{
	int rowNo, colNo;
	double i, j;
	double value;

	rowNo=myCoefficientsMatrix.size();
	colNo=rowNo;

	sparseCoefficientsRow.clear();
	sparseCoefficientsColumn.clear();
	sparseCoefficientsValue.clear();

	for(i=0; i<rowNo; i++)
	{
		for(j=0; j<colNo; j++)
		{
			value=myCoefficientsMatrix[i][j];
			if(value!=0)
			{
				sparseCoefficientsRow.push_back(i);
				sparseCoefficientsColumn.push_back(j);
				sparseCoefficientsValue.push_back(value);
			}
		}
	}

	return;
}

void coefficientsAssembly::assemblyMandelCoefficientsMatrix(double dx, double dy, double G,
	double lambda, double alpha)
{
	addMandelRigidMotion();
	increaseMandelCoefficientsMatrixSize();

	if(gridType=="staggered")
	{
		addMandelStaggeredStressToVDisplacement(dx);
		addMandelStaggeredStress(dx,dy,G,lambda,alpha);
	}
	else if(gridType=="collocated")
	{
		addMandelCollocatedStressToVDisplacement(dx);
		addMandelCollocatedStress(dx,dy,G,lambda,alpha);
	}

	return;
}

void coefficientsAssembly::addMandelRigidMotion()
{
	int i, j;
	int v_P, v_ref;

	for(int FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(j!=vDisplacementFVIndex[0].size()-1 && i==0) // North but not eastern border
		{
			v_ref=getVDisplacementFVPosition(i,vDisplacementFVIndex[0].size()-1);
			coefficientsMatrix[v_P].assign(Nu+Nv+NP,0);
			coefficientsMatrix[v_P][v_P]+=1;
			coefficientsMatrix[v_P][v_ref]-=1;
		}

	}

	return;
}

void coefficientsAssembly::increaseMandelCoefficientsMatrixSize()
{
	int rowNo=coefficientsMatrix.size()+1;
	int colNo=coefficientsMatrix[0].size()+1;

	coefficientsMatrix.resize(rowNo);
	for(int i=0; i<rowNo; i++)
		coefficientsMatrix[i].resize(colNo);

	return;
}

void coefficientsAssembly::addMandelStaggeredStressToVDisplacement(double dx)
{
	int v_P=getVDisplacementFVPosition(0,vDisplacementFVIndex[0].size()-1);
	int sigma_P=coefficientsMatrix[0].size()-1;

	coefficientsMatrix[v_P][sigma_P]-=dx;

	return;
}

void coefficientsAssembly::addMandelStaggeredStress(double dx, double dy, double G, double lambda,
	double alpha)
{
	int i, j;
	int u_P, u_E, v_P, v_S, P_P;
	int sigma_P=coefficientsMatrix[0].size()-1;

	for(int FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		if(j!=vDisplacementFVIndex[0].size()-1 && i==0) // North but not eastern border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_E=getUDisplacementFVPosition(i,j+1);
			v_P=getVDisplacementFVPosition(i,j);
			v_S=getVDisplacementFVPosition(i+1,j);
			P_P=getPressureFVPosition(i,j);

			coefficientsMatrix[sigma_P][u_P]+=-lambda;
			coefficientsMatrix[sigma_P][u_E]+=lambda;
			coefficientsMatrix[sigma_P][v_P]+=(2*G+lambda)*(dx/dy);
			coefficientsMatrix[sigma_P][v_S]+=-(2*G+lambda)*(dx/dy);
			coefficientsMatrix[sigma_P][P_P]+=-alpha*dx;
		}
	}

	coefficientsMatrix[sigma_P][sigma_P]+=dx;

	return;
}

void coefficientsAssembly::addMandelCollocatedStressToVDisplacement(double dx)
{
	int v_P=getVDisplacementFVPosition(0,vDisplacementFVIndex[0].size()-1);
	int sigma_P=coefficientsMatrix[0].size()-1;

	coefficientsMatrix[v_P][sigma_P]-=dx*0.5;

	return;
}

void coefficientsAssembly::addMandelCollocatedStress(double dx, double dy, double G, double lambda,
	double alpha)
{
	int i, j;
	int u_P, u_E, u_W, v_P, v_S, P_P;
	int sigma_P=coefficientsMatrix[0].size()-1;

	for(int FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		if(j<vDisplacementFVIndex[0].size()-1 && i==0 && j>0) // N but not NE nor NW
		{
			u_E=getUDisplacementFVPosition(i,j+1);
			u_W=getUDisplacementFVPosition(i,j-1);
			v_P=getVDisplacementFVPosition(i,j);
			v_S=getVDisplacementFVPosition(i+1,j);
			P_P=getPressureFVPosition(i,j);

			coefficientsMatrix[sigma_P][u_E]+=lambda/2;
			coefficientsMatrix[sigma_P][u_W]+=-lambda/2;
			coefficientsMatrix[sigma_P][v_P]+=(2*G+lambda)*(dx/dy);
			coefficientsMatrix[sigma_P][v_S]+=-(2*G+lambda)*(dx/dy);
			coefficientsMatrix[sigma_P][P_P]+=-alpha*dx;
		}
		else if(j==0 && i==0) // NW
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_E=getUDisplacementFVPosition(i,j+1);
			v_P=getVDisplacementFVPosition(i,j);
			v_S=getVDisplacementFVPosition(i+1,j);
			P_P=getPressureFVPosition(i,j);

			coefficientsMatrix[sigma_P][u_P]+=-lambda*0.5;
			coefficientsMatrix[sigma_P][u_E]+=lambda*0.5;
			coefficientsMatrix[sigma_P][v_P]+=(2*G+lambda)*(dx/dy)*0.5;
			coefficientsMatrix[sigma_P][v_S]+=-(2*G+lambda)*(dx/dy)*0.5;
			coefficientsMatrix[sigma_P][P_P]+=-alpha*dx*0.5;			
		}
	}

	coefficientsMatrix[sigma_P][sigma_P]+=dx*0.5;

	return;
}

void coefficientsAssembly::assemblyMacroPorosityMatrix(double dx, double dy, double dt, double G,
	double lambda, double alpha, double K, double mu_f, double Q, double phi, double phiM,
	double KM, double QM)
{
	double alpham=alpha*phi/(1-phi-phiM);
	double alphaM=alpha*phiM/(1-phi-phiM);
	double propCoef=3*0.4/min(phi*phi,phiM*phiM);
	NPM=NP;

	increaseMacroPorosityCoefficientsMatrixSize();
	addUDisplacementToXMomentum(dx,dy,G,lambda);
	addVDisplacementToXMomentum(dx,dy,G,lambda);
	addPressureToXMomentum(dy,alpham);
	addMacroPressureToXMomentum(dy,alphaM);
	addBCToXMomentum(dx,dy,G);

	addUDisplacementToYMomentum(dx,dy,G,lambda);
	addVDisplacementToYMomentum(dx,dy,G,lambda);
	addPressureToYMomentum(dx,alpham);
	addMacroPressureToYMomentum(dx,alphaM);
	addBCToYMomentum(dx,dy,G);

	addTransientToContinuity(dx,dy,dt,Q);
	addMacroTransientToContinuity(dx,dy,dt,QM);
	addFluidFlowToContinuity(dx,dy,dt,K,mu_f,alpham,G,lambda);
	addMacroFluidFlowToContinuity(dx,dy,dt,KM,mu_f,alpham,alphaM,G,lambda);
	addDisplacementToContinuity(dx,dy,dt,alpham,G,lambda);
	addMacroDisplacementToContinuity(dx,dy,dt,alphaM,G,lambda);
	addMacroPressuretoContinuity(propCoef,K,mu_f);
	addBCToContinuity();
	addMacroBCToContinuity();

	assemblySparseMatrix(coefficientsMatrix);

	return;
}

void coefficientsAssembly::increaseMacroPorosityCoefficientsMatrixSize()
{
	int rowNo=coefficientsMatrix.size()+NP;
	int colNo=coefficientsMatrix[0].size()+NP;

	coefficientsMatrix.resize(rowNo);
	for(int i=0; i<rowNo; i++)
		coefficientsMatrix[i].resize(colNo);

	return;
}

void coefficientsAssembly::addMacroPressureToXMomentum(double dy, double alpha)
{
	if(gridType=="staggered") addStaggeredMacroPressureToXMomentum(dy,alpha);
	else if(gridType=="collocated")
	{
		if(interpScheme=="CDS") addCDSMacroPressureToXMomentum(dy,alpha);
		if(interpScheme=="1DPIS") addCDSMacroPressureToXMomentum(dy,alpha);
		if(interpScheme=="I2DPIS") addCDSMacroPressureToXMomentum(dy,alpha);
		if(interpScheme=="C2DPIS") addCDSMacroPressureToXMomentum(dy,alpha);
	}

	return;
}

void coefficientsAssembly::addStaggeredMacroPressureToXMomentum(double dy, double alpha)
{
	int u_P;
	int P_P, P_W;
	int FVCounter;
	int i, j;

	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(j==0) // FV on the western border
		{
			P_P=getMacroPressureFVPosition(i,j);

			coefficientsMatrix[u_P][P_P]-=-alpha*dy;
		}
		else if(j==uDisplacementFVIndex[0].size()-1) // FV on the eastern border
		{
			P_W=getMacroPressureFVPosition(i,j-1);

			coefficientsMatrix[u_P][P_W]-=alpha*dy;
		}
		else // FV not on the western or eastern border
		{
			P_P=getMacroPressureFVPosition(i,j);
			P_W=getMacroPressureFVPosition(i,j-1);

			coefficientsMatrix[u_P][P_P]-=-alpha*dy;
			coefficientsMatrix[u_P][P_W]-=alpha*dy;
		}
	}

	return;
}

void coefficientsAssembly::addCDSMacroPressureToXMomentum(double dy, double alpha)
{
	int u_P;
	int P_P, P_E, P_W;
	int FVCounter;
	int i, j;
	double value=1;
	
	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{	
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(j==0) // FV on the western border
		{
			P_P=getMacroPressureFVPosition(i,j);
			P_E=getMacroPressureFVPosition(i,j+1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[u_P][P_P]-=-0.5*alpha*dy*value;
			coefficientsMatrix[u_P][P_E]-=-0.5*alpha*dy*value;

			value=1;
		}
		else if(j==uDisplacementFVIndex[0].size()-1) // FV on the eastern border
		{
			P_P=getMacroPressureFVPosition(i,j);
			P_W=getMacroPressureFVPosition(i,j-1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[u_P][P_P]-=0.5*alpha*dy*value;
			coefficientsMatrix[u_P][P_W]-=0.5*alpha*dy*value;

			value=1;
		}
		else // FV not on the western or eastern border
		{
			P_E=getMacroPressureFVPosition(i,j+1);
			P_W=getMacroPressureFVPosition(i,j-1);

			if(i==0 || i==uDisplacementFVIndex.size()-1) value=0.5;

			coefficientsMatrix[u_P][P_E]-=-0.5*alpha*dy*value;
			coefficientsMatrix[u_P][P_W]-=0.5*alpha*dy*value;

			value=1;
		}
	}
	
	return;
}

void coefficientsAssembly::addMacroPressureToYMomentum(double dx, double alpha)
{
	if(gridType=="staggered") addStaggeredMacroPressureToYMomentum(dx,alpha);
	else if(gridType=="collocated")
	{
		if(interpScheme=="CDS") addCDSMacroPressureToYMomentum(dx,alpha);
		if(interpScheme=="1DPIS") addCDSMacroPressureToYMomentum(dx,alpha);
		if(interpScheme=="I2DPIS") addCDSMacroPressureToYMomentum(dx,alpha);
		if(interpScheme=="C2DPIS") addCDSMacroPressureToYMomentum(dx,alpha);
	}

	return;
}

void coefficientsAssembly::addStaggeredMacroPressureToYMomentum(double dx, double alpha)
{
	int v_P;
	int P_P, P_N;
	int FVCounter;
	int i, j;
	
	for(FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(i==0) // FV on the northern border
		{
			P_P=getMacroPressureFVPosition(i,j);

			coefficientsMatrix[v_P][P_P]-=alpha*dx;
		}
		else if(i==vDisplacementFVIndex.size()-1) // FV on the southern border
		{
			P_N=getMacroPressureFVPosition(i-1,j);

			coefficientsMatrix[v_P][P_N]-=-alpha*dx;
		}
		else // FV not on the northern or southern border
		{
			P_P=getMacroPressureFVPosition(i,j);
			P_N=getMacroPressureFVPosition(i-1,j);

			coefficientsMatrix[v_P][P_P]-=alpha*dx;
			coefficientsMatrix[v_P][P_N]-=-alpha*dx;
		}
	}

	return;
}

void coefficientsAssembly::addCDSMacroPressureToYMomentum(double dx, double alpha)
{
	int v_P;
	int P_P, P_N, P_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(i==0) // FV on the northern border
		{
			P_P=getMacroPressureFVPosition(i,j);
			P_S=getMacroPressureFVPosition(i+1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[v_P][P_P]-=0.5*alpha*dx*value;
			coefficientsMatrix[v_P][P_S]-=0.5*alpha*dx*value;

			value=1;
		}
		else if(i==vDisplacementFVIndex.size()-1) // FV on the southern border
		{
			P_P=getMacroPressureFVPosition(i,j);
			P_N=getMacroPressureFVPosition(i-1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[v_P][P_P]-=-0.5*alpha*dx*value;
			coefficientsMatrix[v_P][P_N]-=-0.5*alpha*dx*value;

			value=1;
		}
		else // FV not on the northern or southern border
		{
			P_N=getMacroPressureFVPosition(i-1,j);
			P_S=getMacroPressureFVPosition(i+1,j);

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[v_P][P_N]-=-0.5*alpha*dx*value;
			coefficientsMatrix[v_P][P_S]-=0.5*alpha*dx*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addMacroTransientToContinuity(double dx, double dy, double dt, double Q)
{	
	int P_P;
	int FVCounter;
	int i, j;
	double Mp=(1/Q)*(dx*dy/dt);
	int borderCounter=0;

	for(FVCounter=0; FVCounter<NPM; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getMacroPressureFVPosition(i,j);

		if(gridType=="staggered")
		{
			coefficientsMatrix[P_P][P_P]+=Mp;
		}
		else if(gridType=="collocated")
		{
			if(i==0 || i==pressureFVIndex.size()-1) borderCounter++;
			if(j==0 || j==pressureFVIndex[0].size()-1) borderCounter++;

			coefficientsMatrix[P_P][P_P]+=Mp/pow(2,borderCounter);

			borderCounter=0;
		}
	}
}

void coefficientsAssembly::addMacroDisplacementToContinuity(double dx, double dy, double dt,
	double alpha, double G, double lambda)
{
	if(gridType=="staggered") addStaggeredMacroDisplacementToContinuity(dx,dy,dt,alpha);
	else if(gridType=="collocated")
	{
		if(interpScheme=="CDS") addCDSMacroDisplacementToContinuity(dx,dy,dt,alpha);
		if(interpScheme=="I2DPIS") addI2DPISMacroDisplacementToContinuity(dx,dy,dt,alpha);
	}
	
	return;
}

void coefficientsAssembly::addStaggeredMacroDisplacementToContinuity(double dx, double dy,
	double dt, double alpha)
{
	int u_P, u_E;
	int v_P, v_S;
	int P_P;
	int FVCounter;
	int i, j;

	for(FVCounter=0; FVCounter<NPM; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);
		u_E=getUDisplacementFVPosition(i,j+1);
		v_P=getVDisplacementFVPosition(i,j);
		v_S=getVDisplacementFVPosition(i+1,j);
		P_P=getMacroPressureFVPosition(i,j);

		coefficientsMatrix[P_P][u_P]-=alpha*(dy/dt);
		coefficientsMatrix[P_P][v_P]-=-alpha*(dx/dt);
		coefficientsMatrix[P_P][u_E]-=-alpha*(dy/dt);
		coefficientsMatrix[P_P][v_S]-=alpha*(dx/dt);
	}

	return;
}

void coefficientsAssembly::addCDSMacroDisplacementToContinuity(double dx, double dy, double dt,
	double alpha)
{
	int P_P;
	int u_P, u_E, u_W;
	int v_P, v_N, v_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getMacroPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			v_P=getVDisplacementFVPosition(i,j);
			v_S=getVDisplacementFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][v_P]-=-0.5*alpha*(dx/dt)*value;
			coefficientsMatrix[P_P][v_S]-=0.5*alpha*(dx/dt)*value;

			value=1;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			v_P=getVDisplacementFVPosition(i,j);
			v_N=getVDisplacementFVPosition(i-1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][v_P]-=0.5*alpha*(dx/dt)*value;
			coefficientsMatrix[P_P][v_N]-=-0.5*alpha*(dx/dt)*value;

			value=1;
		}
		else
		{
			v_N=getVDisplacementFVPosition(i-1,j);
			v_S=getVDisplacementFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][v_N]-=-0.5*alpha*(dx/dt)*value;
			coefficientsMatrix[P_P][v_S]-=0.5*alpha*(dx/dt)*value;

			value=1;
		}

		if(j==0) // Western border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_E=getUDisplacementFVPosition(i,j+1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][u_P]-=0.5*alpha*(dy/dt)*value;
			coefficientsMatrix[P_P][u_E]-=-0.5*alpha*(dy/dt)*value;

			value=1;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_W=getUDisplacementFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][u_P]-=-0.5*alpha*(dy/dt)*value;
			coefficientsMatrix[P_P][u_W]-=0.5*alpha*(dy/dt)*value;

			value=1;
		}
		else
		{
			u_E=getUDisplacementFVPosition(i,j+1);
			u_W=getUDisplacementFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][u_E]-=-0.5*alpha*(dy/dt)*value;
			coefficientsMatrix[P_P][u_W]-=0.5*alpha*(dy/dt)*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addI2DPISMacroDisplacementToContinuity(double dx, double dy, double dt,
	double alpha)
{
	int P_P;
	int u_P, u_E, u_W, u_N, u_S, u_NE, u_NW, u_SE, u_SW;
	int v_P, v_E, v_W, v_N, v_S, v_NE, v_NW, v_SE, v_SW;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NPM; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getMacroPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{	
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_E=getVDisplacementFVPosition(i,j+1);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_SE=getVDisplacementFVPosition(i+1,j+1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_P]-=-alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_E]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_W]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_SE]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_SW]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				coefficientsMatrix[P_P][u_P]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				coefficientsMatrix[P_P][u_P]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_E=getVDisplacementFVPosition(i,j+1);
				v_W=getVDisplacementFVPosition(i,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_NW=getVDisplacementFVPosition(i-1,j-1);

				coefficientsMatrix[P_P][u_E]-=-0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=0.25*alpha*(dy/dt);
				coefficientsMatrix[P_P][v_P]-=-(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_P]-=alpha*(dx/dt);
				coefficientsMatrix[P_P][v_N]-=-(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_E]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_W]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_NE]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_NW]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
			}
		}
		else
		{
			if(j==0) // Western border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=-(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_P]-=alpha*(dy/dt);
				coefficientsMatrix[P_P][u_E]-=-(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_N]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_NE]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_S]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_SE]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_N=getUDisplacementFVPosition(i-1,j);
				u_S=getUDisplacementFVPosition(i+1,j);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				coefficientsMatrix[P_P][u_P]-=(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_P]-=-alpha*(dy/dt);
				coefficientsMatrix[P_P][u_W]-=(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_N]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_S]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_NW]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_SW]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_N]-=-0.25*alpha*(dx/dt);
				coefficientsMatrix[P_P][v_S]-=0.25*alpha*(dx/dt);
			}
			else
			{
				u_E=getUDisplacementFVPosition(i,j+1);
				u_W=getUDisplacementFVPosition(i,j-1);
				u_NE=getUDisplacementFVPosition(i-1,j+1);
				u_NW=getUDisplacementFVPosition(i-1,j-1);
				u_SE=getUDisplacementFVPosition(i+1,j+1);
				u_SW=getUDisplacementFVPosition(i+1,j-1);
				v_N=getVDisplacementFVPosition(i-1,j);
				v_S=getVDisplacementFVPosition(i+1,j);
				v_NE=getVDisplacementFVPosition(i-1,j+1);
				v_NW=getVDisplacementFVPosition(i-1,j-1);
				v_SE=getVDisplacementFVPosition(i+1,j+1);
				v_SW=getVDisplacementFVPosition(i+1,j-1);

				coefficientsMatrix[P_P][u_E]-=-(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_W]-=(alpha*dy*dy*dy)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_NE]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_NW]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_SE]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][u_SW]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_N]-=-(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_S]-=(alpha*dx*dx*dx)/(2*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_NE]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_NW]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_SE]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][v_SW]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt);
			}
		}
	}

	return;
}

void coefficientsAssembly::addMacroFluidFlowToContinuity(double dx, double dy, double dt, double K,
	double mu_f, double alpham, double alphaM, double G, double lambda)
{
	if(gridType=="staggered") addStaggeredMacroFluidFlowToContinuity(dx,dy,K,mu_f);
	else if(gridType=="collocated") 
	{
		addCollocatedMacroFluidFlowToContinuity(dx,dy,K,mu_f);

		if(interpScheme=="I2DPIS")
		{
			addI2DPISMacroFluidFlowToMicroContinuity(dx,dy,dt,alphaM,G);
			// addI2DPISMacroFluidFlowToMacroContinuity(dx,dy,dt,alphaM,G);
			// addI2DPISMicroFluidFlowToMacroContinuity(dx,dy,dt,alpham,G);
		}
	}

	return;
}

void coefficientsAssembly::addStaggeredMacroFluidFlowToContinuity(double dx, double dy, double K,
	double mu_f)
{
	int P_P, P_E, P_W, P_N, P_S;
	int FVCounter;
	int i, j;
	int bcType;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getMacroPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			P_S=getMacroPressureFVPosition(i+1,j);

			coefficientsMatrix[P_P][P_S]-=(K/mu_f)*(dx/dy);
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dx/dy);

			bcType=boundaryConditionType[0][3];
			if(bcType==1) coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dx/dy);
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			P_N=getMacroPressureFVPosition(i-1,j);

			coefficientsMatrix[P_P][P_N]-=(K/mu_f)*(dx/dy);
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dx/dy);

			bcType=boundaryConditionType[2][3];
			if(bcType==1) coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dx/dy);
		}
		else
		{
			P_N=getMacroPressureFVPosition(i-1,j);
			P_S=getMacroPressureFVPosition(i+1,j);

			coefficientsMatrix[P_P][P_N]-=(K/mu_f)*(dx/dy);
			coefficientsMatrix[P_P][P_S]-=(K/mu_f)*(dx/dy);
			coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dx/dy);
		}

		if(j==0) // Western border
		{
			P_E=getMacroPressureFVPosition(i,j+1);

			coefficientsMatrix[P_P][P_E]-=(K/mu_f)*(dy/dx);
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dy/dx);

			bcType=boundaryConditionType[1][3];
			if(bcType==1) coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dy/dx);
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			P_W=getMacroPressureFVPosition(i,j-1);

			coefficientsMatrix[P_P][P_W]-=(K/mu_f)*(dy/dx);
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dy/dx);

			bcType=boundaryConditionType[3][3];
			if(bcType==1) coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dy/dx);
		}
		else
		{
			P_E=getMacroPressureFVPosition(i,j+1);
			P_W=getMacroPressureFVPosition(i,j-1);

			coefficientsMatrix[P_P][P_E]-=(K/mu_f)*(dy/dx);
			coefficientsMatrix[P_P][P_W]-=(K/mu_f)*(dy/dx);
			coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dy/dx);
		}
	}

	return;
}

void coefficientsAssembly::addCollocatedMacroFluidFlowToContinuity(double dx, double dy, double K,
	double mu_f)
{
	int P_P, P_E, P_W, P_N, P_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getMacroPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			P_S=getMacroPressureFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_S]-=(K/mu_f)*(dx/dy)*value;
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dx/dy)*value;

			value=1;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			P_N=getMacroPressureFVPosition(i-1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_N]-=(K/mu_f)*(dx/dy)*value;
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dx/dy)*value;

			value=1;
		}
		else
		{
			P_N=getMacroPressureFVPosition(i-1,j);
			P_S=getMacroPressureFVPosition(i+1,j);

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			coefficientsMatrix[P_P][P_N]-=(K/mu_f)*(dx/dy)*value;
			coefficientsMatrix[P_P][P_S]-=(K/mu_f)*(dx/dy)*value;
			coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dx/dy)*value;

			value=1;
		}

		if(j==0) // Western border
		{
			P_E=getMacroPressureFVPosition(i,j+1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_E]-=(K/mu_f)*(dy/dx)*value;
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dy/dx)*value;

			value=1;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			P_W=getMacroPressureFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_W]-=(K/mu_f)*(dy/dx)*value;
			coefficientsMatrix[P_P][P_P]+=(K/mu_f)*(dy/dx)*value;

			value=1;
		}
		else
		{
			P_E=getMacroPressureFVPosition(i,j+1);
			P_W=getMacroPressureFVPosition(i,j-1);

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			coefficientsMatrix[P_P][P_E]-=(K/mu_f)*(dy/dx)*value;
			coefficientsMatrix[P_P][P_W]-=(K/mu_f)*(dy/dx)*value;
			coefficientsMatrix[P_P][P_P]+=2*(K/mu_f)*(dy/dx)*value;

			value=1;
		}
	}

	return;
}

void coefficientsAssembly::addI2DPISMacroFluidFlowToMicroContinuity(double dx, double dy, double dt,
	double alpha, double G)
{
	int P_P, PM_E, PM_W, PM_N, PM_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_S=getMacroPressureFVPosition(i+1,j);
	
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_N=getMacroPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;		
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_N=getMacroPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_N=getMacroPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else
		{
			if(j==0) // Western border
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_N=getMacroPressureFVPosition(i-1,j);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_N=getMacroPressureFVPosition(i-1,j);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);				
			}
			else
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_N=getMacroPressureFVPosition(i-1,j);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
			}
		}
	}
	
	return;
}

void coefficientsAssembly::addI2DPISMacroFluidFlowToMacroContinuity(double dx, double dy, double dt,
	double alpha, double G)
{
	int P_P, PM_E, PM_W, PM_N, PM_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getMacroPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_S=getMacroPressureFVPosition(i+1,j);
	
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_N=getMacroPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;		
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_N=getMacroPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_N=getMacroPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else
		{
			if(j==0) // Western border
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_N=getMacroPressureFVPosition(i-1,j);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_N=getMacroPressureFVPosition(i-1,j);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);				
			}
			else
			{
				PM_E=getMacroPressureFVPosition(i,j+1);
				PM_W=getMacroPressureFVPosition(i,j-1);
				PM_N=getMacroPressureFVPosition(i-1,j);
				PM_S=getMacroPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
			}
		}
	}
	
	return;
}

void coefficientsAssembly::addI2DPISMicroFluidFlowToMacroContinuity(double dx, double dy, double dt,
	double alpha, double G)
{
	int P_P, PM_E, PM_W, PM_N, PM_S;
	int FVCounter;
	int i, j;
	double value=1;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getMacroPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				PM_E=getPressureFVPosition(i,j+1);
				PM_S=getPressureFVPosition(i+1,j);
	
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getPressureFVPosition(i,j-1);
				PM_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				PM_E=getPressureFVPosition(i,j+1);
				PM_W=getPressureFVPosition(i,j-1);
				PM_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				PM_E=getPressureFVPosition(i,j+1);
				PM_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;		
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getPressureFVPosition(i,j-1);
				PM_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx)/(16*G*dt)*dy;
			}
			else
			{
				PM_E=getPressureFVPosition(i,j+1);
				PM_W=getPressureFVPosition(i,j-1);
				PM_N=getPressureFVPosition(i-1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx)/(16*G*dt)*dy;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx)/(16*G*dt)*dy;
			}
		}
		else
		{
			if(j==0) // Western border
			{
				PM_E=getPressureFVPosition(i,j+1);
				PM_N=getPressureFVPosition(i-1,j);
				PM_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				PM_W=getPressureFVPosition(i,j-1);
				PM_N=getPressureFVPosition(i-1,j);
				PM_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy)/(16*G*dt)*dx;
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);				
			}
			else
			{
				PM_E=getPressureFVPosition(i,j+1);
				PM_W=getPressureFVPosition(i,j-1);
				PM_N=getPressureFVPosition(i-1,j);
				PM_S=getPressureFVPosition(i+1,j);

				coefficientsMatrix[P_P][PM_N]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_S]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_E]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][PM_W]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
				coefficientsMatrix[P_P][P_P]+=2*(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt);
			}
		}
	}
	
	return;
}

void coefficientsAssembly::addMacroPressuretoContinuity(double propCoef, double K, double mu_f)
{
	int i, j;
	int P_P, PM_P;

	for(int FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);
		PM_P=getMacroPressureFVPosition(i,j);

		coefficientsMatrix[P_P][P_P]-=-propCoef*K/mu_f;
		coefficientsMatrix[P_P][PM_P]-=propCoef*K/mu_f;
		coefficientsMatrix[PM_P][P_P]-=propCoef*K/mu_f;
		coefficientsMatrix[PM_P][PM_P]-=-propCoef*K/mu_f;
	}

	return;
}

void coefficientsAssembly::addMacroBCToContinuity()
{
	if(gridType=="collocated")
		for(int counter=0; counter<4; counter++)
			addMacroDirichletBCToContinuity(counter);

	return;
}

void coefficientsAssembly::addMacroDirichletBCToContinuity(int counter)
{
	int P_P;
	int bcType;
	int i, j;

	bcType=boundaryConditionType[counter][3];

	if(bcType==1)
	{
		switch(counter)
		{
			case 0:
				for(j=0; j<pressureFVIndex[0].size(); j++)
				{
					P_P=getMacroPressureFVPosition(0,j);

					coefficientsMatrix[P_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[P_P][P_P]+=1;
				}
				break;

			case 1:
				for(i=0; i<pressureFVIndex.size(); i++)
				{
					P_P=getMacroPressureFVPosition(i,0);

					coefficientsMatrix[P_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[P_P][P_P]+=1;
				}
				break;

			case 2:
				for(j=0; j<pressureFVIndex[0].size(); j++)
				{
					P_P=getMacroPressureFVPosition(pressureFVIndex.size()-1,j);

					coefficientsMatrix[P_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[P_P][P_P]+=1;
				}
				break;

			case 3:
				for(i=0; i<pressureFVIndex.size(); i++)
				{
					P_P=getMacroPressureFVPosition(i,pressureFVIndex[0].size()-1);

					coefficientsMatrix[P_P].assign(Nu+Nv+NP+NPM,0);
					coefficientsMatrix[P_P][P_P]+=1;
				}
				break;
		}
	}

	return;
}