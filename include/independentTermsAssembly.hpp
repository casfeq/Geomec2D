/*
	This header is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The class defined 
	here contains the functions for assembly of the independent terms of the linear system which
	represents the discretized problem of poroelasticity.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.
*/

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class independentTermsAssembly
{
public:
	// Class variables
	vector<double> independentTermsArray;
	vector<vector<int>> boundaryConditionType;
	vector<vector<double>> boundaryConditionValue;
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

	// Class functions
	void resizeIndependentTermsArray();
	int getUDisplacementFVPosition(int,int);
	int getVDisplacementFVPosition(int,int);
	int getPressureFVPosition(int,int);
	int getMacroPressureFVPosition(int,int);
	int getFVPosition(int,int,int);
	void zeroIndependentTermsArray();
	void assemblyIndependentTermsArray(double,double,double,double,double,double,double,double,
		double,double,double,vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,
		int);
	void addUDisplacement(double,double,double,double);
	void addVDisplacement(double,double,double,double,double,double);
	void addPressure(double,double,double,double,double,double,double,vector<vector<double>>,
		vector<vector<double>>,vector<vector<double>>,int,double,double);
	void addBC(int);
	void addStaggeredUDisplacement(double,double,double,double);
	void addStaggeredVDisplacement(double,double,double,double,double,double);
	void addStaggeredPressure(double,double,double,double,double,double,double,
		vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,int);
	void addCollocatedUDisplacement(double,double,double,double);
	void addCollocatedVDisplacement(double,double,double,double,double,double);
	void addCollocatedPressure(double,double,double,double,double,double,double,
		vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,int,double,double);
	void addCDSDisplacement(double,double,double,double,vector<vector<double>>,
		vector<vector<double>>,int);
	void add1DPISPressure(double,double,double,double,double,double,vector<vector<double>>,int);
	void addI2DPISDisplacement(double,double,double,double,vector<vector<double>>,
		vector<vector<double>>,int);
	void addI2DPISPressure(double,double,double,double,double,vector<vector<double>>,int);
	void addC2DPISDisplacement(double,double,double,double,double,double,vector<vector<double>>,
		vector<vector<double>>,int);
	void addC2DPISPressure(double,double,double,double,double,double,vector<vector<double>>,int);
	void addDirichletBC(int,int);
	void increaseMandelIndependentTermsArray();
	void assemblyMandelIndependentTermsArray(double,double);
	void addMandelRigidMotion();
	void addMandelForce(double,double);
	void assemblyMacroIndependentTermsArray(double,double,double,double,double,double,double,double,
		double,double,double,vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,
		vector<vector<double>>,int,double,double,double,double);
	void increaseMacroIndependentTermsArray();
	void addMacroPressure(double,double,double,double,double,double,double,vector<vector<double>>,
		vector<vector<double>>,vector<vector<double>>,int,double,double);
	void addStaggeredMacroPressure(double,double,double,double,double,double,double,
		vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,int);
	void addCollocatedMacroPressure(double,double,double,double,double,double,double,
		vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,int,double,double);
	void addCDSMacroDisplacement(double,double,double,double,vector<vector<double>>,
		vector<vector<double>>,int);
	void addI2DPISMacroDisplacement(double,double,double,double,vector<vector<double>>,
		vector<vector<double>>,int);
	void addI2DPISMacroPressure(double,double,double,double,double,vector<vector<double>>,int);

	// Constructor
	independentTermsAssembly(vector<vector<int>>,vector<vector<double>>,int,int,int,
		vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,
		vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,string,
		string);

	// Destructor
	~independentTermsAssembly();
};

independentTermsAssembly::independentTermsAssembly(vector<vector<int>> bcType,
	vector<vector<double>> bcValue, int numberOfActiveUDisplacementFV,
	int numberOfActiveVDisplacementFV, int numberOfActivePressureFV, vector<vector<int>> idU,
	vector<vector<int>> idV, vector<vector<int>> idP, vector<vector<int>> cooU,
	vector<vector<int>> cooV, vector<vector<int>> cooP, vector<vector<int>> horFStatus,
	vector<vector<int>> verFStatus, string myGridType, string myInterpScheme)
{
	boundaryConditionType=bcType;
	boundaryConditionValue=bcValue;
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

	resizeIndependentTermsArray();
}

independentTermsAssembly::~independentTermsAssembly(){}

void independentTermsAssembly::resizeIndependentTermsArray()
{
	int rowNo, colNo;

	rowNo=NP+Nv+Nu;
	independentTermsArray.resize(rowNo);

	return;
}

int independentTermsAssembly::getUDisplacementFVPosition(int x, int y)
{
	int uDisplacementFVPosition;

	uDisplacementFVPosition=uDisplacementFVIndex[x][y]-1;

	return uDisplacementFVPosition;
}

int independentTermsAssembly::getVDisplacementFVPosition(int x, int y)
{
	int vDisplacementFVPosition;

	vDisplacementFVPosition=vDisplacementFVIndex[x][y]+Nu-1;

	return vDisplacementFVPosition;
}

int independentTermsAssembly::getPressureFVPosition(int x, int y)
{
	int pressureFVPosition;

	pressureFVPosition=pressureFVIndex[x][y]+Nu+Nv-1;

	return pressureFVPosition;
}

int independentTermsAssembly::getMacroPressureFVPosition(int x, int y)
{
	int pressureFVPosition;

	pressureFVPosition=pressureFVIndex[x][y]+Nu+Nv+NP-1;

	return pressureFVPosition;
}

int independentTermsAssembly::getFVPosition(int variable, int i, int j)
{
	switch(variable)
	{
		case 0:
			return getUDisplacementFVPosition(i,j);
			break;

		case 1:
			return getVDisplacementFVPosition(i,j);
			break;

		case 2:
			return getPressureFVPosition(i,j);
			break;

		case 3:
			return getMacroPressureFVPosition(i,j);
			break;
	}

	return 0;
}

void independentTermsAssembly::zeroIndependentTermsArray()
{
	int rowNo=independentTermsArray.size();

	for(int i=0; i<rowNo; i++)
	{
		independentTermsArray[i]=0;
	}

	return;
}

void independentTermsAssembly::assemblyIndependentTermsArray(double dx, double dy, double dt,
	double G, double lambda, double alpha, double K, double mu_f, double Q, double rho, double g,
	vector<vector<double>> uField, vector<vector<double>> vField, vector<vector<double>> pField,
	int timeStep)
{
	zeroIndependentTermsArray();
	addUDisplacement(dx,dy,G,lambda);
	addVDisplacement(dx,dy,G,lambda,rho,g);
	addPressure(Q,dx,dy,dt,K,mu_f,alpha,uField,vField,pField,timeStep,G,lambda);
	addBC(3);

	return;
}

void independentTermsAssembly::addUDisplacement(double dx, double dy, double G, double lambda)
{
	if(gridType=="staggered") addStaggeredUDisplacement(dx,dy,G,lambda);
	else if(gridType=="collocated") addCollocatedUDisplacement(dx,dy,G,lambda);

	return;
}

void independentTermsAssembly::addVDisplacement(double dx, double dy, double G, double lambda,
	double rho, double g)
{
	if(gridType=="staggered") addStaggeredVDisplacement(dx,dy,G,lambda,rho,g);
	else if(gridType=="collocated") addCollocatedVDisplacement(dx,dy,G,lambda,rho,g);

	return;
}

void independentTermsAssembly::addPressure(double Q, double dx, double dy, double dt, double K,
	double mu_f, double alpha, vector<vector<double>> uField, vector<vector<double>> vField,
	vector<vector<double>> pField, int timeStep, double G, double lambda)
{
	if(gridType=="staggered") addStaggeredPressure(Q,dx,dy,dt,K,mu_f,alpha,uField,vField,pField,
		timeStep);
	else if(gridType=="collocated") addCollocatedPressure(Q,dx,dy,dt,K,mu_f,alpha,uField,vField,
		pField,timeStep,G,lambda);

	return;
}

void independentTermsAssembly::addBC(int var)
{
	if(gridType=="collocated")
	{
		for(int border=0; border<4; border++)
			for(int variable=0; variable<var; variable++)
				addDirichletBC(border,variable);
	}

	return;
}

void independentTermsAssembly::addStaggeredUDisplacement(double dx, double dy, double G,
	double lambda)
{
	int u_P;
	int FVCounter;
	int i, j;
	int bcType;
	double bcValue;

	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(i==0) // Northern border
		{	
			bcType=boundaryConditionType[0][0];
			bcValue=boundaryConditionValue[0][0];
			if(bcType==1) independentTermsArray[u_P]+=2*G*(dx/dy)*bcValue;
			else if(bcType==-1) independentTermsArray[u_P]+=bcValue*dx;
		}
		else if(i==uDisplacementFVIndex.size()-1) // Southern border
		{
			bcType=boundaryConditionType[2][0];
			bcValue=boundaryConditionValue[2][0];
			if(bcType==1) independentTermsArray[u_P]+=2*G*(dx/dy)*bcValue;
			else if(bcType==-1) independentTermsArray[u_P]+=bcValue*dx;
		}

		if(j==0) // Western border
		{
			bcType=boundaryConditionType[1][0];
			bcValue=boundaryConditionValue[1][0];
			if(bcType==1)
			{
				independentTermsArray[u_P]=bcValue;
				break;
			}
			else if(bcType==-1) independentTermsArray[u_P]+=bcValue*dy;
		}
		else if(j==uDisplacementFVIndex[0].size()-1) // Eastern border
		{
			bcType=boundaryConditionType[3][0];
			bcValue=boundaryConditionValue[3][0];
			if(bcType==1)
			{
				independentTermsArray[u_P]=bcValue;
				break;
			}
			else if(bcType==-1) independentTermsArray[u_P]+=bcValue*dy;
		}
	}

	return;
}

void independentTermsAssembly::addStaggeredVDisplacement(double dx, double dy, double G,
	double lambda, double rho, double g)
{
	int v_P;
	int FVCounter;
	int i, j;
	int borderCounter;
	int bcType;
	double bcValue;

	for(FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		borderCounter=1;

		if(j==0) // Western border
		{
			bcType=boundaryConditionType[1][1];
			bcValue=boundaryConditionValue[1][1];
			if(bcType==1) independentTermsArray[v_P]+=2*G*(dy/dx)*bcValue;
			else if(bcType==-1) independentTermsArray[v_P]+=bcValue*dy;
		}
		else if(j==vDisplacementFVIndex[0].size()-1) // Eastern border
		{
			bcType=boundaryConditionType[3][1];
			bcValue=boundaryConditionValue[3][1];
			if(bcType==1) independentTermsArray[v_P]+=2*G*(dy/dx)*bcValue;
			else if(bcType==-1) independentTermsArray[v_P]+=bcValue*dy;
		}

		if(i==0) // Northern border
		{	
			bcType=boundaryConditionType[0][1];
			bcValue=boundaryConditionValue[0][1];
			if(bcType==1)
			{
				independentTermsArray[v_P]=bcValue;
				break;
			}
			else if(bcType==-1) independentTermsArray[v_P]+=bcValue*dx;

			borderCounter++;
		}
		else if(i==vDisplacementFVIndex.size()-1) // Southern border
		{
			bcType=boundaryConditionType[2][1];
			bcValue=boundaryConditionValue[2][1];
			if(bcType==1)
			{
				independentTermsArray[v_P]=bcValue;
				break;
			}
			else if(bcType==-1) independentTermsArray[v_P]+=bcValue*dx;

			borderCounter++;
		}

		independentTermsArray[v_P]+=rho*g*dx*dy/borderCounter;
	}

	return;
}

void independentTermsAssembly::addStaggeredPressure(double Q, double dx, double dy, double dt,
	double K, double mu_f, double alpha, vector<vector<double>> uField,
	vector<vector<double>> vField, vector<vector<double>> pField, int timeStep)
{
	int u_P, u_E, v_P, v_S, P_P;
	double uP, uE, vP, vS, PP;
	int FVCounter;
	int i, j;
	int bcType;
	double bcValue;
	double MP=(1/Q)*(dx*dy)/dt;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);
		u_E=getUDisplacementFVPosition(i,j+1);
		v_P=getVDisplacementFVPosition(i,j);
		v_S=getVDisplacementFVPosition(i+1,j);
		P_P=getPressureFVPosition(i,j);

		uP=uField[u_P][timeStep];
		uE=uField[u_E][timeStep];
		vP=vField[v_P-Nu][timeStep];
		vS=vField[v_S-Nu][timeStep];
		PP=pField[P_P-Nu-Nv][timeStep];

		independentTermsArray[P_P]+=MP*PP;
		independentTermsArray[P_P]-=alpha*(dy/dt)*uP;
		independentTermsArray[P_P]-=-alpha*(dy/dt)*uE;
		independentTermsArray[P_P]-=-alpha*(dx/dt)*vP;
		independentTermsArray[P_P]-=alpha*(dx/dt)*vS;

		if(i==0) // Northern border
		{
			bcType=boundaryConditionType[0][2];
			bcValue=boundaryConditionValue[0][2];

			if(bcType==1) independentTermsArray[P_P]+=2*(K/mu_f)*(dx/dy)*bcValue;
			else if(bcType==0) independentTermsArray[P_P]+=(K/mu_f)*bcValue*dx;
			else if(bcType==-1) independentTermsArray[P_P]+=bcValue*dx;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			bcType=boundaryConditionType[2][2];
			bcValue=boundaryConditionValue[2][2];

			if(bcType==1) independentTermsArray[P_P]+=2*(K/mu_f)*(dx/dy)*bcValue;
			else if(bcType==0) independentTermsArray[P_P]-=(K/mu_f)*bcValue*dx;
			else if(bcType==-1) independentTermsArray[P_P]-=bcValue*dx;
		}

		if(j==0) // Western border
		{
			bcType=boundaryConditionType[1][2];
			bcValue=boundaryConditionValue[1][2];

			if(bcType==1) independentTermsArray[P_P]+=2*(K/mu_f)*(dy/dx)*bcValue;
			else if(bcType==0) independentTermsArray[P_P]-=(K/mu_f)*bcValue*dy;
			else if(bcType==-1) independentTermsArray[P_P]-=bcValue*dy;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			bcType=boundaryConditionType[3][2];
			bcValue=boundaryConditionValue[3][2];

			if(bcType==1) independentTermsArray[P_P]+=2*(K/mu_f)*(dy/dx)*bcValue;
			else if(bcType==0) independentTermsArray[P_P]+=(K/mu_f)*bcValue*dy;
			else if(bcType==-1) independentTermsArray[P_P]+=bcValue*dy;
		}
	}

	return;
}

void independentTermsAssembly::addCollocatedUDisplacement(double dx, double dy, double G,
	double lambda)
{
	int u_P;
	int FVCounter;
	int i, j;
	int bcType;
	double bcValue;
	double value=1;

	for(FVCounter=0; FVCounter<Nu; FVCounter++)
	{
		i=uDisplacementFVCoordinates[FVCounter][0]-1;
		j=uDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);

		if(i==0) // Northern border
		{
			bcValue=boundaryConditionValue[0][0];

			if(j==0 || j==uDisplacementFVIndex[0].size()-1) value=0.5;

			independentTermsArray[u_P]+=bcValue*dx*value;

			value=1;
		}

		if(i==uDisplacementFVIndex.size()-1) // Southern border
		{
			bcValue=boundaryConditionValue[2][0];

			if(j==0 || j==uDisplacementFVIndex[0].size()-1) value=0.5;

			independentTermsArray[u_P]+=bcValue*dx*value;

			value=1;
		}

		if(j==0) // Western border
		{
			bcValue=boundaryConditionValue[1][0];

			if(i==0 || i==uDisplacementFVIndex.size()-1) value=0.5;

			independentTermsArray[u_P]+=bcValue*dy*value;

			value=1;
		}
		
		if(j==uDisplacementFVIndex[0].size()-1) // Eastern border
		{
			bcValue=boundaryConditionValue[3][0];

			if(i==0 || i==uDisplacementFVIndex.size()-1) value=0.5;

			independentTermsArray[u_P]+=bcValue*dy*value;

			value=1;
		}
	}

	return;
}

void independentTermsAssembly::addCollocatedVDisplacement(double dx, double dy, double G,
	double lambda, double rho, double g)
{
	int v_P;
	int FVCounter;
	int i, j;
	int bcType;
	double bcValue;
	double value=1;
	int borderCounter=0;
	
	for(FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(i==0) // Northern border
		{
			bcValue=boundaryConditionValue[0][1];

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) value=0.5;

			independentTermsArray[v_P]+=bcValue*dx*value;

			value=1;
			borderCounter++;
		}

		if(i==vDisplacementFVIndex.size()-1) // Southern border
		{
			bcValue=boundaryConditionValue[2][1];

			if(j==0 || j==vDisplacementFVIndex[0].size()-1) value=0.5;

			independentTermsArray[v_P]+=bcValue*dx*value;

			value=1;
			borderCounter++;
		}

		if(j==0) // Western border
		{
			bcValue=boundaryConditionValue[1][1];

			if(i==0 || i==vDisplacementFVIndex.size()-1) value=0.5;

			independentTermsArray[v_P]+=bcValue*dy*value;

			value=1;
			borderCounter++;
		}
		
		if(j==vDisplacementFVIndex[0].size()-1) // Eastern border
		{
			bcValue=boundaryConditionValue[3][1];

			if(i==0 || i==vDisplacementFVIndex.size()-1) value=0.5;

			independentTermsArray[v_P]+=bcValue*dy*value;

			value=1;
			borderCounter++;
		}

		independentTermsArray[v_P]+=rho*g*dx*dy/pow(2,borderCounter);
		borderCounter=0;
	}

	return;
}

void independentTermsAssembly::addCollocatedPressure(double Q, double dx, double dy, double dt,
	double K, double mu_f, double alpha, vector<vector<double>> uField,
	vector<vector<double>> vField, vector<vector<double>> pField, int timeStep, double G,
	double lambda)
{
	int P_P;
	double PP;
	int FVCounter;
	int i, j;
	double MP=(1/Q)*(dx*dy)/dt;
	int borderCounter=0;
	double sizeFV=1;
	int bcType;
	double bcValue;
	
	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);
		PP=pField[P_P-Nu-Nv][timeStep];

		if(i==0) // Northern border
		{
			bcType=boundaryConditionType[0][2];
			bcValue=boundaryConditionValue[0][2];

			if(j==0 || j==pressureFVIndex[0].size()-1) sizeFV=0.5;
			if(bcType==0) independentTermsArray[P_P]+=(K/mu_f)*bcValue*dx*sizeFV;

			sizeFV=1;
			borderCounter++;
		}
		if(i==pressureFVIndex.size()-1) // Southern border
		{
			bcType=boundaryConditionType[2][2];
			bcValue=boundaryConditionValue[2][2];

			if(j==0 || j==pressureFVIndex[0].size()-1) sizeFV=0.5;
			if(bcType==0) independentTermsArray[P_P]-=(K/mu_f)*bcValue*dx*sizeFV;

			sizeFV=1;
			borderCounter++;
		}
		if(j==0) borderCounter++;
		if(j==pressureFVIndex[0].size()-1) borderCounter++;

		independentTermsArray[P_P]+=MP*PP/pow(2,borderCounter);
		borderCounter=0;
	}

	if(interpScheme=="CDS") addCDSDisplacement(dx,dy,dt,alpha,uField,vField,timeStep);
	else if(interpScheme=="1DPIS")
	{
		addCDSDisplacement(dx,dy,dt,alpha,uField,vField,timeStep);
		add1DPISPressure(dx,dy,dt,alpha,G,lambda,pField,timeStep);
	}
	else if(interpScheme=="I2DPIS")
	{
		addI2DPISDisplacement(dx,dy,dt,alpha,uField,vField,timeStep);
		addI2DPISPressure(dx,dy,dt,alpha,G,pField,timeStep);
	}
	else if(interpScheme=="C2DPIS")
	{
		addC2DPISDisplacement(dx,dy,dt,alpha,G,lambda,uField,vField,timeStep);
		addC2DPISPressure(dx,dy,dt,alpha,G,lambda,pField,timeStep);
	}

	return;
}

void independentTermsAssembly::addCDSDisplacement(double dx, double dy, double dt, double alpha,
	vector<vector<double>> uField, vector<vector<double>> vField, int timeStep)
{
	int P_P;
	int u_P, u_E, u_W;
	int v_P, v_N, v_S;
	double uP, uE, uW;
	double vP, vN, vS;
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

			vP=vField[v_P-Nu][timeStep];
			vS=vField[v_S-Nu][timeStep];

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=-0.5*alpha*(dx/dt)*vP*value;
			independentTermsArray[P_P]-=0.5*alpha*(dx/dt)*vS*value;

			value=1;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			v_P=getVDisplacementFVPosition(i,j);
			v_N=getVDisplacementFVPosition(i-1,j);

			vP=vField[v_P-Nu][timeStep];
			vN=vField[v_N-Nu][timeStep];

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=0.5*alpha*(dx/dt)*vP*value;
			independentTermsArray[P_P]-=-0.5*alpha*(dx/dt)*vN*value;

			value=1;
		}
		else
		{
			v_N=getVDisplacementFVPosition(i-1,j);
			v_S=getVDisplacementFVPosition(i+1,j);

			vN=vField[v_N-Nu][timeStep];
			vS=vField[v_S-Nu][timeStep];

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=-0.5*alpha*(dx/dt)*vN*value;
			independentTermsArray[P_P]-=0.5*alpha*(dx/dt)*vS*value;

			value=1;
		}
		
		if(j==0) // Western border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_E=getUDisplacementFVPosition(i,j+1);

			uP=uField[u_P][timeStep];
			uE=uField[u_E][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=0.5*alpha*(dy/dt)*uP*value;
			independentTermsArray[P_P]-=-0.5*alpha*(dy/dt)*uE*value;

			value=1;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_W=getUDisplacementFVPosition(i,j-1);

			uP=uField[u_P][timeStep];
			uW=uField[u_W][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=-0.5*alpha*(dy/dt)*uP*value;
			independentTermsArray[P_P]-=0.5*alpha*(dy/dt)*uW*value;

			value=1;
		}
		else
		{
			u_E=getUDisplacementFVPosition(i,j+1);
			u_W=getUDisplacementFVPosition(i,j-1);

			uE=uField[u_E][timeStep];
			uW=uField[u_W][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=-0.5*alpha*(dy/dt)*uE*value;
			independentTermsArray[P_P]-=0.5*alpha*(dy/dt)*uW*value;

			value=1;
		}
	}

	return;
}

void independentTermsAssembly::add1DPISPressure(double dx, double dy, double dt, double alpha,
	double G, double lambda, vector<vector<double>> pField, int timeStep)
{
	int P_P, P_E, P_W, P_N, P_S;
	double PP, PE, PW, PN, PS;
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
			P_P=getPressureFVPosition(i,j);
			P_S=getPressureFVPosition(i+1,j);

			PP=pField[P_P-Nu-Nv][timeStep];
			PS=pField[P_S-Nu-Nv][timeStep];

		
			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=-(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*PP*value;
			independentTermsArray[P_P]-=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*PS*value;

			value=1;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			P_P=getPressureFVPosition(i,j);
			P_N=getPressureFVPosition(i-1,j);

			PP=pField[P_P-Nu-Nv][timeStep];
			PN=pField[P_N-Nu-Nv][timeStep];

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=-(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*PP*value;
			independentTermsArray[P_P]-=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*PN*value;

			value=1;
		}
		else
		{
			P_P=getPressureFVPosition(i,j);
			P_N=getPressureFVPosition(i-1,j);
			P_S=getPressureFVPosition(i+1,j);

			PP=pField[P_P-Nu-Nv][timeStep];
			PN=pField[P_N-Nu-Nv][timeStep];
			PS=pField[P_S-Nu-Nv][timeStep];

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=-2*(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*PP*value;
			independentTermsArray[P_P]-=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*PN*value;
			independentTermsArray[P_P]-=(alpha*alpha*dy)/(8*(2*G+lambda)*dt)*dx*PS*value;

			value=1;
		}
		
		if(j==0) // Western border
		{
			P_P=getPressureFVPosition(i,j);
			P_E=getPressureFVPosition(i,j+1);

			PP=pField[P_P-Nu-Nv][timeStep];
			PE=pField[P_E-Nu-Nv][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=-(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*PP*value;
			independentTermsArray[P_P]-=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*PE*value;

			value=1;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			P_P=getPressureFVPosition(i,j);
			P_W=getPressureFVPosition(i,j-1);

			PP=pField[P_P-Nu-Nv][timeStep];
			PW=pField[P_W-Nu-Nv][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=-(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*PP*value;
			independentTermsArray[P_P]-=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*PW*value;

			value=1;
		}
		else
		{
			P_P=getPressureFVPosition(i,j);
			P_E=getPressureFVPosition(i,j+1);
			P_W=getPressureFVPosition(i,j-1);

			PP=pField[P_P-Nu-Nv][timeStep];
			PE=pField[P_E-Nu-Nv][timeStep];
			PW=pField[P_W-Nu-Nv][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=-2*(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*PP*value;
			independentTermsArray[P_P]-=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*PE*value;
			independentTermsArray[P_P]-=(alpha*alpha*dx)/(8*(2*G+lambda)*dt)*dy*PW*value;

			value=1;
		}
	}

	return;
}

void independentTermsAssembly::addI2DPISDisplacement(double dx, double dy, double dt,
	double alpha, vector<vector<double>> uField, vector<vector<double>> vField, int timeStep)
{
	int P_P;
	int u_P, u_E, u_W, u_N, u_S, u_NE, u_NW, u_SE, u_SW;
	int v_P, v_E, v_W, v_N, v_S, v_NE, v_NW, v_SE, v_SW;
	double uP, uE, uW, uN, uS, uNE, uNW, uSE, uSW;
	double vP, vE, vW, vN, vS, vNE, vNW, vSE, vSW;
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

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{	
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vE=vField[v_E-Nu][timeStep];
				vW=vField[v_W-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];
				vSE=vField[v_SE-Nu][timeStep];
				vSW=vField[v_SW-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vP;
				independentTermsArray[P_P]-=-alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vS;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vE;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vW;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vSE;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vSW;
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

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];

				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vE=vField[v_E-Nu][timeStep];
				vW=vField[v_W-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vNE=vField[v_NE-Nu][timeStep];
				vNW=vField[v_NW-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vP;
				independentTermsArray[P_P]-=-alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vN;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vE;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vW;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vNE;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vNW;
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

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				uN=uField[u_N][timeStep];
				uS=uField[u_S][timeStep];
				uNE=uField[u_NE][timeStep];
				uSE=uField[u_SE][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=-(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uP;
				independentTermsArray[P_P]-=-alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uE;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uN;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uNE;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uS;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uSE;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
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

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				uN=uField[u_N][timeStep];
				uS=uField[u_S][timeStep];
				uNW=uField[u_NW][timeStep];
				uSW=uField[u_SW][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uP;
				independentTermsArray[P_P]-=alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uW;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uN;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uS;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uNW;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uSW;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				uNE=uField[u_NE][timeStep];
				uNW=uField[u_NW][timeStep];
				uSE=uField[u_SE][timeStep];
				uSW=uField[u_SW][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];
				vNE=vField[v_NE-Nu][timeStep];
				vNW=vField[v_NW-Nu][timeStep];
				vSE=vField[v_SE-Nu][timeStep];
				vSW=vField[v_SW-Nu][timeStep];

				independentTermsArray[P_P]-=-(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uE;
				independentTermsArray[P_P]-=(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uW;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uNE;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uNW;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uSE;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uSW;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vN;
				independentTermsArray[P_P]-=(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vS;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vNE;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vNW;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vSE;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vSW;
			}
		}
	}

	return;
}

void independentTermsAssembly::addI2DPISPressure(double dx, double dy, double dt, double alpha,
	double G, vector<vector<double>> pField, int timeStep)
{
	int P_P, P_E, P_W, P_N, P_S;
	double PP, PE, PW, PN, PS;
	int FVCounter;
	int i, j;

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

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt)
					*(PS-PP);
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_N=getPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt)
					*(PN-PP);
			}
		}
		else
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt)
					*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt)
					*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt)
					*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt)
					*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt)
					*(PN-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt)
					*(PS-PP);
			}
		}
	}

	return;
}

void independentTermsAssembly::addC2DPISDisplacement(double dx, double dy, double dt,
	double alpha, double G, double lambda, vector<vector<double>> uField,
	vector<vector<double>> vField, int timeStep)
{
	int P_P;
	int u_P, u_E, u_W, u_N, u_S, u_NE, u_NW, u_SE, u_SW;
	int v_P, v_E, v_W, v_N, v_S, v_NE, v_NW, v_SE, v_SW;
	double uP, uE, uW, uN, uS, uNE, uNW, uSE, uSW;
	double vP, vE, vW, vN, vS, vNE, vNW, vSE, vSW;
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

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{	
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				uSE=uField[u_SE][timeStep];
				uSW=uField[u_SW][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vE=vField[v_E-Nu][timeStep];
				vW=vField[v_W-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];
				vSE=vField[v_SE-Nu][timeStep];
				vSW=vField[v_SW-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt)*uE;
				independentTermsArray[P_P]-=((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt)*uW;
				independentTermsArray[P_P]-=-((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt)*uSE;
				independentTermsArray[P_P]-=-((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt)*uSW;
				independentTermsArray[P_P]-=((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt)*vS;
				independentTermsArray[P_P]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vE;
				independentTermsArray[P_P]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vW;
				independentTermsArray[P_P]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vSE;
				independentTermsArray[P_P]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vSW;
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

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];

				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				uNE=uField[u_NE][timeStep];
				uNW=uField[u_NW][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vE=vField[v_E-Nu][timeStep];
				vW=vField[v_W-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vNE=vField[v_NE-Nu][timeStep];
				vNW=vField[v_NW-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt)*uE;
				independentTermsArray[P_P]-=((G+lambda)*dx*dy)/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)
					*alpha*(dx/dt)*uW;
				independentTermsArray[P_P]-=-((G+lambda)*dx*dy)
					/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)*alpha*(dx/dt)*uNE;
				independentTermsArray[P_P]-=-((G+lambda)*dx*dy)
					/(16*(2*G+lambda)*dx*dx+16*G*dy*dy)*alpha*(dx/dt)*uNW;
				independentTermsArray[P_P]-=-((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt)*vN;
				independentTermsArray[P_P]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vE;
				independentTermsArray[P_P]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vW;
				independentTermsArray[P_P]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vNE;
				independentTermsArray[P_P]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vNW;
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

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				uN=uField[u_N][timeStep];
				uS=uField[u_S][timeStep];
				uNE=uField[u_NE][timeStep];
				uSE=uField[u_SE][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];
				vNE=vField[v_NE-Nu][timeStep];
				vSE=vField[v_SE-Nu][timeStep];

				independentTermsArray[P_P]-=-((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uN;
				independentTermsArray[P_P]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uNE;
				independentTermsArray[P_P]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uS;
				independentTermsArray[P_P]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uSE;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
				independentTermsArray[P_P]-=((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*vN;
				independentTermsArray[P_P]-=((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*vS;
				independentTermsArray[P_P]-=-((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*vNE;
				independentTermsArray[P_P]-=-((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*vSE;
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

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				uN=uField[u_N][timeStep];
				uS=uField[u_S][timeStep];
				uNW=uField[u_NW][timeStep];
				uSW=uField[u_SW][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];
				vNW=vField[v_NW-Nu][timeStep];
				vSW=vField[v_SW-Nu][timeStep];

				independentTermsArray[P_P]-=((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uN;
				independentTermsArray[P_P]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uS;
				independentTermsArray[P_P]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uNW;
				independentTermsArray[P_P]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uSW;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
				independentTermsArray[P_P]-=((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*vN;
				independentTermsArray[P_P]-=((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*vS;
				independentTermsArray[P_P]-=-((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*vNW;
				independentTermsArray[P_P]-=-((G+lambda)*dx*dx)/(16*G*dx*dx+16*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*vSW;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				uNE=uField[u_NE][timeStep];
				uNW=uField[u_NW][timeStep];
				uSE=uField[u_SE][timeStep];
				uSW=uField[u_SW][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];
				vNE=vField[v_NE-Nu][timeStep];
				vNW=vField[v_NW-Nu][timeStep];
				vSE=vField[v_SE-Nu][timeStep];
				vSW=vField[v_SW-Nu][timeStep];

				independentTermsArray[P_P]-=-((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=((2*G+lambda)*dy*dy)
					/(2*G*dx*dx+2*(2*G+lambda)*dy*dy)*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uNE;
				independentTermsArray[P_P]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uNW;
				independentTermsArray[P_P]-=-(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uSE;
				independentTermsArray[P_P]-=(G*dx*dx)/(4*G*dx*dx+4*(2*G+lambda)*dy*dy)*alpha
					*(dy/dt)*uSW;
				independentTermsArray[P_P]-=-((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt)*vN;
				independentTermsArray[P_P]-=((2*G+lambda)*dx*dx)/(2*(2*G+lambda)*dx*dx+2*G*dy*dy)
					*alpha*(dx/dt)*vS;
				independentTermsArray[P_P]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vNE;
				independentTermsArray[P_P]-=-(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vNW;
				independentTermsArray[P_P]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vSE;
				independentTermsArray[P_P]-=(G*dy*dy)/(4*(2*G+lambda)*dx*dx+4*G*dx*dx)
					*alpha*(dx/dt)*vSW;
			}
		}
	}

	return;
}

void independentTermsAssembly::addC2DPISPressure(double dx, double dy, double dt, double alpha,
	double G, double lambda, vector<vector<double>> pField, int timeStep)
{
	int P_P, P_E, P_W, P_N, P_S;
	double PP, PE, PW, PN, PS;
	int FVCounter;
	int i, j;

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

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(8*(2*G+lambda)*dx*dx+8*G*dy*dy)
					*alpha*(dx/dt)*(PS-PP);
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_N=getPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(8*(2*G+lambda)*dx*dx+8*G*dy*dy)
					*alpha*(dx/dt)*(PN-PP);
			}
		}
		else
		{
			if(j==0) // Western border
			{
				P_E=getPressureFVPosition(i,j+1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(8*G*dx*dx+8*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(8*G*dx*dx+8*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else
			{
				P_E=getPressureFVPosition(i,j+1);
				P_W=getPressureFVPosition(i,j-1);
				P_N=getPressureFVPosition(i-1,j);
				P_S=getPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv][timeStep];
				PN=pField[P_N-Nu-Nv][timeStep];
				PS=pField[P_S-Nu-Nv][timeStep];

				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(8*G*dx*dx+8*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*(PE-PP);
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(8*G*dx*dx+8*(2*G+lambda)*dy*dy)
					*alpha*(dy/dt)*(PW-PP);
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(8*(2*G+lambda)*dx*dx+8*G*dy*dy)
					*alpha*(dx/dt)*(PN-PP);
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(8*(2*G+lambda)*dx*dx+8*G*dy*dy)
					*alpha*(dx/dt)*(PS-PP);
			}
		}
	}

	return;
}

void independentTermsAssembly::addDirichletBC(int border, int variable)
{	
	int bcType, var, i, j;
	double bcValue;
	vector<vector<int>> FVIndex;

	bcType=boundaryConditionType[border][variable];
	bcValue=boundaryConditionValue[border][variable];

	switch(variable)
	{
		case 0:
			FVIndex=uDisplacementFVIndex;
			break;

		case 1:
			FVIndex=vDisplacementFVIndex;
			break;

		case 2:
			FVIndex=pressureFVIndex;
			break;

		case 3:
			FVIndex=pressureFVIndex;
			break;
	}

	if(bcType==1)
	{
		switch(border)
		{
			case 0:
				for(j=0; j<FVIndex[0].size(); j++)
				{
					var=getFVPosition(variable,0,j);

					independentTermsArray[var]=bcValue;
				}
				break;

			case 1:
				for(i=0; i<FVIndex.size(); i++)
				{
					var=getFVPosition(variable,i,0);

					independentTermsArray[var]=bcValue;
				}
				break;

			case 2:
				for(j=0; j<FVIndex[0].size(); j++)
				{
					var=getFVPosition(variable,FVIndex.size()-1,j);

					independentTermsArray[var]=bcValue;
				}
				break;

			case 3:
				for(i=0; i<FVIndex.size(); i++)
				{
					var=getFVPosition(variable,i,FVIndex[0].size()-1);

					independentTermsArray[var]=bcValue;
				}
				break;
		}
	}

	return;
}

void independentTermsAssembly::increaseMandelIndependentTermsArray()
{
	int rowNo=independentTermsArray.size()+1;

	independentTermsArray.resize(rowNo);

	return;
}

void independentTermsAssembly::assemblyMandelIndependentTermsArray(double F, double L)
{
	addMandelRigidMotion();
	addMandelForce(F,L);

	return;	
}

void independentTermsAssembly::addMandelRigidMotion()
{
	int i, j;
	int v_P;

	for(int FVCounter=0; FVCounter<Nv; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		v_P=getVDisplacementFVPosition(i,j);

		if(j!=vDisplacementFVIndex[0].size()-1 && i==0) // North but not eastern border
			independentTermsArray[v_P]=0;
	}

	return;
}

void independentTermsAssembly::addMandelForce(double F, double L)
{
	independentTermsArray[independentTermsArray.size()-1]+=F*L;

	return;
}

void independentTermsAssembly::increaseMacroIndependentTermsArray()
{
	NPM=NP;

	int rowNo=independentTermsArray.size()+NPM;

	independentTermsArray.resize(rowNo);

	return;
}

void independentTermsAssembly::assemblyMacroIndependentTermsArray(double dx, double dy, double dt,
	double G, double lambda, double alpha, double K, double mu_f, double Q, double rho, double g,
	vector<vector<double>> uField, vector<vector<double>> vField, vector<vector<double>> pField,
	vector<vector<double>> pMField, int timeStep, double phi, double phiM, double KM, double QM)
{
	double alpham=alpha*phi/(1-phi-phiM);
	double alphaM=alpha*phiM/(1-phi-phiM);

	zeroIndependentTermsArray();
	addUDisplacement(dx,dy,G,lambda);
	addVDisplacement(dx,dy,G,lambda,rho,g);
	addPressure(Q,dx,dy,dt,K,mu_f,alpham,uField,vField,pField,timeStep,G,lambda);
	addMacroPressure(QM,dx,dy,dt,KM,mu_f,alphaM,uField,vField,pMField,timeStep,G,lambda);
	addBC(4);

	return;
}

void independentTermsAssembly::addMacroPressure(double Q, double dx, double dy, double dt, double K,
	double mu_f, double alpha, vector<vector<double>> uField, vector<vector<double>> vField,
	vector<vector<double>> pField, int timeStep, double G, double lambda)
{
	if(gridType=="staggered") addStaggeredMacroPressure(Q,dx,dy,dt,K,mu_f,alpha,uField,vField,
		pField,timeStep);
	else if(gridType=="collocated") addCollocatedMacroPressure(Q,dx,dy,dt,K,mu_f,alpha,uField,
		vField,pField,timeStep,G,lambda);

	return;
}

void independentTermsAssembly::addStaggeredMacroPressure(double Q, double dx, double dy, double dt,
	double K, double mu_f, double alpha, vector<vector<double>> uField,
	vector<vector<double>> vField, vector<vector<double>> pField, int timeStep)
{
	int u_P, u_E, v_P, v_S, P_P;
	double uP, uE, vP, vS, PP;
	int FVCounter;
	int i, j;
	int bcType;
	double bcValue;
	double MP=(1/Q)*(dx*dy)/dt;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);
		u_E=getUDisplacementFVPosition(i,j+1);
		v_P=getVDisplacementFVPosition(i,j);
		v_S=getVDisplacementFVPosition(i+1,j);
		P_P=getMacroPressureFVPosition(i,j);

		uP=uField[u_P][timeStep];
		uE=uField[u_E][timeStep];
		vP=vField[v_P-Nu][timeStep];
		vS=vField[v_S-Nu][timeStep];
		PP=pField[P_P-Nu-Nv-NP][timeStep];

		independentTermsArray[P_P]+=MP*PP;
		independentTermsArray[P_P]-=alpha*(dy/dt)*uP;
		independentTermsArray[P_P]-=-alpha*(dy/dt)*uE;
		independentTermsArray[P_P]-=-alpha*(dx/dt)*vP;
		independentTermsArray[P_P]-=alpha*(dx/dt)*vS;

		if(i==0) // Northern border
		{
			bcType=boundaryConditionType[0][3];
			bcValue=boundaryConditionValue[0][3];

			if(bcType==1) independentTermsArray[P_P]+=2*(K/mu_f)*(dx/dy)*bcValue;
			else if(bcType==0) independentTermsArray[P_P]+=(K/mu_f)*bcValue*dx;
			else if(bcType==-1) independentTermsArray[P_P]+=bcValue*dx;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			bcType=boundaryConditionType[2][3];
			bcValue=boundaryConditionValue[2][3];

			if(bcType==1) independentTermsArray[P_P]+=2*(K/mu_f)*(dx/dy)*bcValue;
			else if(bcType==0) independentTermsArray[P_P]-=(K/mu_f)*bcValue*dx;
			else if(bcType==-1) independentTermsArray[P_P]-=bcValue*dx;
		}

		if(j==0) // Western border
		{
			bcType=boundaryConditionType[1][3];
			bcValue=boundaryConditionValue[1][3];

			if(bcType==1) independentTermsArray[P_P]+=2*(K/mu_f)*(dy/dx)*bcValue;
			else if(bcType==0) independentTermsArray[P_P]-=(K/mu_f)*bcValue*dy;
			else if(bcType==-1) independentTermsArray[P_P]-=bcValue*dy;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			bcType=boundaryConditionType[3][3];
			bcValue=boundaryConditionValue[3][3];

			if(bcType==1) independentTermsArray[P_P]+=2*(K/mu_f)*(dy/dx)*bcValue;
			else if(bcType==0) independentTermsArray[P_P]+=(K/mu_f)*bcValue*dy;
			else if(bcType==-1) independentTermsArray[P_P]+=bcValue*dy;
		}
	}

	return;
}

void independentTermsAssembly::addCollocatedMacroPressure(double Q, double dx, double dy, double dt,
	double K, double mu_f, double alpha, vector<vector<double>> uField,
	vector<vector<double>> vField, vector<vector<double>> pField, int timeStep, double G,
	double lambda)
{
	int P_P;
	double PP;
	int FVCounter;
	int i, j;
	double MP=(1/Q)*(dx*dy)/dt;
	int borderCounter=0;
	double sizeFV=1;
	int bcType;
	double bcValue;
	
	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getMacroPressureFVPosition(i,j);
		PP=pField[P_P-Nu-Nv-NP][timeStep];

		if(i==0) // Northern border
		{
			bcType=boundaryConditionType[0][3];
			bcValue=boundaryConditionValue[0][3];

			if(j==0 || j==pressureFVIndex[0].size()-1) sizeFV=0.5;
			if(bcType==0) independentTermsArray[P_P]+=(K/mu_f)*bcValue*dx*sizeFV;

			sizeFV=1;
			borderCounter++;
		}
		if(i==pressureFVIndex.size()-1) // Southern border
		{
			bcType=boundaryConditionType[2][3];
			bcValue=boundaryConditionValue[2][3];

			if(j==0 || j==pressureFVIndex[0].size()-1) sizeFV=0.5;
			if(bcType==0) independentTermsArray[P_P]-=(K/mu_f)*bcValue*dx*sizeFV;

			sizeFV=1;
			borderCounter++;
		}
		if(j==0) borderCounter++;
		if(j==pressureFVIndex[0].size()-1) borderCounter++;

		independentTermsArray[P_P]+=MP*PP/pow(2,borderCounter);
		borderCounter=0;
	}

	if(interpScheme=="CDS") addCDSMacroDisplacement(dx,dy,dt,alpha,uField,vField,timeStep);
	else if(interpScheme=="I2DPIS")
	{
		addI2DPISMacroDisplacement(dx,dy,dt,alpha,uField,vField,timeStep);
		addI2DPISMacroPressure(dx,dy,dt,alpha,G,pField,timeStep);
	}

	return;
}

void independentTermsAssembly::addCDSMacroDisplacement(double dx, double dy, double dt,
	double alpha, vector<vector<double>> uField, vector<vector<double>> vField, int timeStep)
{
	int P_P;
	int u_P, u_E, u_W;
	int v_P, v_N, v_S;
	double uP, uE, uW;
	double vP, vN, vS;
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

			vP=vField[v_P-Nu][timeStep];
			vS=vField[v_S-Nu][timeStep];

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=-0.5*alpha*(dx/dt)*vP*value;
			independentTermsArray[P_P]-=0.5*alpha*(dx/dt)*vS*value;

			value=1;
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			v_P=getVDisplacementFVPosition(i,j);
			v_N=getVDisplacementFVPosition(i-1,j);

			vP=vField[v_P-Nu][timeStep];
			vN=vField[v_N-Nu][timeStep];

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=0.5*alpha*(dx/dt)*vP*value;
			independentTermsArray[P_P]-=-0.5*alpha*(dx/dt)*vN*value;

			value=1;
		}
		else
		{
			v_N=getVDisplacementFVPosition(i-1,j);
			v_S=getVDisplacementFVPosition(i+1,j);

			vN=vField[v_N-Nu][timeStep];
			vS=vField[v_S-Nu][timeStep];

			if(j==0 || j==pressureFVIndex[0].size()-1) value=0.5;

			independentTermsArray[P_P]-=-0.5*alpha*(dx/dt)*vN*value;
			independentTermsArray[P_P]-=0.5*alpha*(dx/dt)*vS*value;

			value=1;
		}
		
		if(j==0) // Western border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_E=getUDisplacementFVPosition(i,j+1);

			uP=uField[u_P][timeStep];
			uE=uField[u_E][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=0.5*alpha*(dy/dt)*uP*value;
			independentTermsArray[P_P]-=-0.5*alpha*(dy/dt)*uE*value;

			value=1;
		}
		else if(j==pressureFVIndex[0].size()-1) // Eastern border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_W=getUDisplacementFVPosition(i,j-1);

			uP=uField[u_P][timeStep];
			uW=uField[u_W][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=-0.5*alpha*(dy/dt)*uP*value;
			independentTermsArray[P_P]-=0.5*alpha*(dy/dt)*uW*value;

			value=1;
		}
		else
		{
			u_E=getUDisplacementFVPosition(i,j+1);
			u_W=getUDisplacementFVPosition(i,j-1);

			uE=uField[u_E][timeStep];
			uW=uField[u_W][timeStep];

			if(i==0 || i==pressureFVIndex.size()-1) value=0.5;

			independentTermsArray[P_P]-=-0.5*alpha*(dy/dt)*uE*value;
			independentTermsArray[P_P]-=0.5*alpha*(dy/dt)*uW*value;

			value=1;
		}
	}

	return;
}

void independentTermsAssembly::addI2DPISMacroDisplacement(double dx, double dy, double dt,
	double alpha, vector<vector<double>> uField, vector<vector<double>> vField, int timeStep)
{
	int P_P;
	int u_P, u_E, u_W, u_N, u_S, u_NE, u_NW, u_SE, u_SW;
	int v_P, v_E, v_W, v_N, v_S, v_NE, v_NW, v_SE, v_SW;
	double uP, uE, uW, uN, uS, uNE, uNW, uSE, uSW;
	double vP, vE, vW, vN, vS, vNE, vNW, vSE, vSW;
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
				u_P=getUDisplacementFVPosition(i,j);
				u_E=getUDisplacementFVPosition(i,j+1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{	
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_S=getVDisplacementFVPosition(i+1,j);

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vE=vField[v_E-Nu][timeStep];
				vW=vField[v_W-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];
				vSE=vField[v_SE-Nu][timeStep];
				vSW=vField[v_SW-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vP;
				independentTermsArray[P_P]-=-alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vS;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vE;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vW;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vSE;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vSW;
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

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				u_P=getUDisplacementFVPosition(i,j);
				u_W=getUDisplacementFVPosition(i,j-1);
				v_P=getVDisplacementFVPosition(i,j);
				v_N=getVDisplacementFVPosition(i-1,j);

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];

				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				vP=vField[v_P-Nu][timeStep];
				vE=vField[v_E-Nu][timeStep];
				vW=vField[v_W-Nu][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vNE=vField[v_NE-Nu][timeStep];
				vNW=vField[v_NW-Nu][timeStep];

				independentTermsArray[P_P]-=-0.25*alpha*(dy/dt)*uE;
				independentTermsArray[P_P]-=0.25*alpha*(dy/dt)*uW;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vP;
				independentTermsArray[P_P]-=-alpha*(dx/dt)*vP;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vN;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vE;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vW;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vNE;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vNW;
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

				uP=uField[u_P][timeStep];
				uE=uField[u_E][timeStep];
				uN=uField[u_N][timeStep];
				uS=uField[u_S][timeStep];
				uNE=uField[u_NE][timeStep];
				uSE=uField[u_SE][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=-(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uP;
				independentTermsArray[P_P]-=-alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=-(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uE;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uN;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uNE;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uS;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uSE;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
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

				uP=uField[u_P][timeStep];
				uW=uField[u_W][timeStep];
				uN=uField[u_N][timeStep];
				uS=uField[u_S][timeStep];
				uNW=uField[u_NW][timeStep];
				uSW=uField[u_SW][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];

				independentTermsArray[P_P]-=(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uP;
				independentTermsArray[P_P]-=alpha*(dy/dt)*uP;
				independentTermsArray[P_P]-=(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uW;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uN;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uS;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uNW;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uSW;
				independentTermsArray[P_P]-=-0.25*alpha*(dx/dt)*vN;
				independentTermsArray[P_P]-=0.25*alpha*(dx/dt)*vS;
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

				uE=uField[u_E][timeStep];
				uW=uField[u_W][timeStep];
				uNE=uField[u_NE][timeStep];
				uNW=uField[u_NW][timeStep];
				uSE=uField[u_SE][timeStep];
				uSW=uField[u_SW][timeStep];
				vN=vField[v_N-Nu][timeStep];
				vS=vField[v_S-Nu][timeStep];
				vNE=vField[v_NE-Nu][timeStep];
				vNW=vField[v_NW-Nu][timeStep];
				vSE=vField[v_SE-Nu][timeStep];
				vSW=vField[v_SW-Nu][timeStep];

				independentTermsArray[P_P]-=-(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uE;
				independentTermsArray[P_P]-=(alpha*dy*dy*dy)/(2*(dy*dy+dx*dx)*dt)*uW;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uNE;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uNW;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uSE;
				independentTermsArray[P_P]-=(alpha*dx*dx*dy)/(4*(dx*dx+dy*dy)*dt)*uSW;
				independentTermsArray[P_P]-=-(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vN;
				independentTermsArray[P_P]-=(alpha*dx*dx*dx)/(2*(dy*dy+dx*dx)*dt)*vS;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vNE;
				independentTermsArray[P_P]-=-(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vNW;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vSE;
				independentTermsArray[P_P]-=(alpha*dx*dy*dy)/(4*(dx*dx+dy*dy)*dt)*vSW;
			}
		}
	}

	return;
}

void independentTermsAssembly::addI2DPISMacroPressure(double dx, double dy, double dt, double alpha,
	double G, vector<vector<double>> pField, int timeStep)
{
	int P_P, P_E, P_W, P_N, P_S;
	double PP, PE, PW, PN, PS;
	int FVCounter;
	int i, j;

	for(FVCounter=0; FVCounter<NP; FVCounter++)
	{
		i=pressureFVCoordinates[FVCounter][0]-1;
		j=pressureFVCoordinates[FVCounter][1]-1;

		P_P=getPressureFVPosition(i,j);

		if(i==0) // Northern border
		{
			if(j==0) // Western border
			{
				P_E=getMacroPressureFVPosition(i,j+1);
				P_S=getMacroPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv-NP][timeStep];
				PS=pField[P_S-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getMacroPressureFVPosition(i,j-1);
				P_S=getMacroPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv-NP][timeStep];
				PS=pField[P_S-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else
			{
				P_E=getMacroPressureFVPosition(i,j+1);
				P_W=getMacroPressureFVPosition(i,j-1);
				P_S=getMacroPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv-NP][timeStep];
				PW=pField[P_W-Nu-Nv-NP][timeStep];
				PS=pField[P_S-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt)
					*(PS-PP);
			}
		}
		else if(i==pressureFVIndex.size()-1) // Southern border
		{
			if(j==0) // Western border
			{
				P_E=getMacroPressureFVPosition(i,j+1);
				P_N=getMacroPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv-NP][timeStep];
				PN=pField[P_N-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getMacroPressureFVPosition(i,j-1);
				P_N=getMacroPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv-NP][timeStep];
				PN=pField[P_N-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
			}
			else
			{
				P_E=getMacroPressureFVPosition(i,j+1);
				P_W=getMacroPressureFVPosition(i,j-1);
				P_N=getMacroPressureFVPosition(i-1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv-NP][timeStep];
				PW=pField[P_W-Nu-Nv-NP][timeStep];
				PN=pField[P_N-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dx)/(16*G*dt)*dy*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt)
					*(PN-PP);
			}
		}
		else
		{
			if(j==0) // Western border
			{
				P_E=getMacroPressureFVPosition(i,j+1);
				P_N=getMacroPressureFVPosition(i-1,j);
				P_S=getMacroPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv-NP][timeStep];
				PN=pField[P_N-Nu-Nv-NP][timeStep];
				PS=pField[P_S-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt)
					*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else if(j==pressureFVIndex[0].size()-1) // Eastern border
			{
				P_W=getMacroPressureFVPosition(i,j-1);
				P_N=getMacroPressureFVPosition(i-1,j);
				P_S=getMacroPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PW=pField[P_W-Nu-Nv-NP][timeStep];
				PN=pField[P_N-Nu-Nv-NP][timeStep];
				PS=pField[P_S-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt)
					*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PN-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy)/(16*G*dt)*dx*(PS-PP);
			}
			else
			{
				P_E=getMacroPressureFVPosition(i,j+1);
				P_W=getMacroPressureFVPosition(i,j-1);
				P_N=getMacroPressureFVPosition(i-1,j);
				P_S=getMacroPressureFVPosition(i+1,j);

				PP=pField[P_P-Nu-Nv][timeStep];
				PE=pField[P_E-Nu-Nv-NP][timeStep];
				PW=pField[P_W-Nu-Nv-NP][timeStep];
				PN=pField[P_N-Nu-Nv-NP][timeStep];
				PS=pField[P_S-Nu-Nv-NP][timeStep];

				independentTermsArray[P_P]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt)
					*(PE-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dx*dy*dy*dy)/(8*G*(dx*dx+dy*dy)*dt)
					*(PW-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt)
					*(PN-PP);
				independentTermsArray[P_P]-=(alpha*alpha*dy*dx*dx*dx)/(8*G*(dx*dx+dy*dy)*dt)
					*(PS-PP);
			}
		}
	}

	return;
}