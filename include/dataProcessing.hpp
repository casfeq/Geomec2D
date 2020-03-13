/*
	This header is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The class defined 
	here contains the functions for post-processing of the data obtained with the solution of the
	discretized problem of poroelasticity.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.
*/

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

class dataProcessing
{
public:
	// Class variables
	int aisleNo;
	vector<vector<vector<double>>> uDisplacement3DField;
	vector<vector<vector<double>>> vDisplacement3DField;
	vector<vector<vector<double>>> pressure3DField;
	vector<vector<vector<double>>> macroPressure3DField;
	vector<vector<vector<double>>> strain3DField;
	string gridType;
	vector<double> mandelRoots;
	double mandelTransCoef;

	// Class structures
	struct errorNorm
	{
		double p=0;
		double eps=0;
		double u=0;
		double v=0;
	} myErrorNorm;

	// Class functions
	void resize3DFields(vector<vector<int>>,vector<vector<int>>,vector<vector<int>>);
	void storeUDisplacem3DField(vector<vector<int>>,vector<vector<double>>);
	void storeVDisplacem3DField(vector<vector<int>>,vector<vector<double>>);
	void storePressure3DField(vector<vector<int>>,vector<vector<double>>);
	void storeStrain3DField(vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,double,
		double,vector<vector<double>>,vector<vector<double>>);
	void export1DFieldToTxt(vector<double>,string);
	void exportSealedColumnNumericalSolution(double,double,double,int,string);
	void exportSealedColumnAnalyticalSolution(double,double,double,double,double,double,double,
		double,double,int,double,string);
	void exportTerzaghiNumericalSolution(double,double,double,int,string);
	double terzaghiPAnalyticalSolution(double,double,double,double,double,double,double,double,
		double,double,double);
	double terzaghiVAnalyticalSolution(double,double,double,double,double,double,double,double,
		double,double,double);
	void exportTerzaghiAnalyticalSolution(double,double,double,double,double,double,double,double,
		double,int,double,string);
	void getTerzaghiErrorNorm(double,double,double,double,double,double,double,double,double,double,
		double,double,double);
	void saveTerzaghiErrorNorm(double,double,double,double);
	void exportMandelNumericalSolution(double,double,double,double,double,int,string);
	void findMandelRoots(double,double,double,double,double,double,double);
	double mandelPAnalyticalSolution(double,double,double,double,double);
	double mandelEAnalyticalSolution(double,double,double,double,double,double,double,double,
		double,double);
	double mandelUAnalyticalSolution(double,double,double,double,double,double,double,double,double,
		double,double);
	double mandelVAnalyticalSolution(double,double,double,double,double,double,double,double,double,
		double,double);
	void exportMandelAnalyticalSolution(double,double,double,double,double,double,double,double,
		double,double,double,double,int,string);
	void storeMacroPressure3DField(vector<vector<int>>,vector<vector<double>>);
	void exportMacroPressureHSolution(double,double,double,int,string);
	void exportMacroPressureTSolution(double,double,double,int,string);
	void exportStripfootSolution(double,double,double,double,int,string);

	// Constructor
	dataProcessing(vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,
		vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,string,string,double,
		double);

	// Destructor
	~dataProcessing();
};

dataProcessing::dataProcessing(vector<vector<int>> idU, vector<vector<int>> idV,
	vector<vector<int>> idP, vector<vector<double>> uField, vector<vector<double>> vField,
	vector<vector<double>> pField, string myGridType, string myInterpScheme, double dx, double dy)
{
	aisleNo=pField[0].size();
	if(myGridType=="staggered") gridType=myGridType;
	else if(myGridType=="collocated") gridType=myGridType+"+"+myInterpScheme;

	resize3DFields(idU,idV,idP);
	storeUDisplacem3DField(idU,uField);
	storeVDisplacem3DField(idV,vField);
	storePressure3DField(idP,pField);
	storeStrain3DField(idU,idV,idP,dx,dy,uField,vField);
}

dataProcessing::~dataProcessing(){}

void dataProcessing::resize3DFields(vector<vector<int>> idU, vector<vector<int>> idV,
	vector<vector<int>> idP)
{
	int rowNo, colNo;

	// Resize uDisplacement3DField
	rowNo=idU.size();
	colNo=idU[0].size();
	uDisplacement3DField.resize(rowNo);
	for(int i=0; i<rowNo; i++)
	{
		uDisplacement3DField[i].resize(colNo);

		for(int j=0; j<colNo; j++)
			uDisplacement3DField[i][j].resize(aisleNo);
	}

	// Resize vDisplacement3DField
	rowNo=idV.size();
	colNo=idV[0].size();
	vDisplacement3DField.resize(rowNo);
	for(int i=0; i<rowNo; i++)
	{
		vDisplacement3DField[i].resize(colNo);

		for(int j=0; j<colNo; j++)
			vDisplacement3DField[i][j].resize(aisleNo);
	}

	// Resize pressure3DField, strain3DField and macroPressure3DField
	rowNo=idP.size();
	colNo=idP[0].size();
	pressure3DField.resize(rowNo);
	strain3DField.resize(rowNo);
	macroPressure3DField.resize(rowNo);
	for(int i=0; i<rowNo; i++)
	{
		pressure3DField[i].resize(colNo);
		strain3DField[i].resize(colNo);
		macroPressure3DField[i].resize(colNo);

		for(int j=0; j<colNo; j++)
		{
			pressure3DField[i][j].resize(aisleNo);
			strain3DField[i][j].resize(aisleNo);
			macroPressure3DField[i][j].resize(aisleNo);
		}
	}

	return;
}

void dataProcessing::storeUDisplacem3DField(vector<vector<int>> idU, vector<vector<double>> uField)
{	
	int rowNo=idU.size();
	int colNo=idU[0].size();
	int position;
	double value;

	for(int i=0; i<rowNo; i++)
	{
		for(int j=0; j<colNo; j++)
		{
			position=idU[i][j]-1;

			for(int k=0; k<aisleNo; k++)
			{
				value=uField[position][k];
				uDisplacement3DField[i][j][k]=value;
			}
		}
	}

	return;
}

void dataProcessing::storeVDisplacem3DField(vector<vector<int>> idV, vector<vector<double>> vField)
{	
	int rowNo=idV.size();
	int colNo=idV[0].size();
	int position;
	double value;

	for(int i=0; i<rowNo; i++)
	{
		for(int j=0; j<colNo; j++)
		{
			position=idV[i][j]-1;

			for(int k=0; k<aisleNo; k++)
			{
				value=vField[position][k];
				vDisplacement3DField[i][j][k]=value;
			}
		}
	}

	return;
}

void dataProcessing::storePressure3DField(vector<vector<int>> idP, vector<vector<double>> pField)
{	
	int rowNo=idP.size();
	int colNo=idP[0].size();
	int position;
	double value;

	for(int i=0; i<rowNo; i++)
	{
		for(int j=0; j<colNo; j++)
		{
			position=idP[i][j]-1;

			for(int k=0; k<aisleNo; k++)
			{
				value=pField[position][k];
				pressure3DField[i][j][k]=value;
			}
		}
	}

	return;
}

void dataProcessing::storeStrain3DField(vector<vector<int>> idU, vector<vector<int>> idV,
	vector<vector<int>> idP, double dx, double dy, vector<vector<double>> uField,
	vector<vector<double>> vField)
{
	int rowNo=idP.size();
	int colNo=idP[0].size();
	int u_P, u_E, u_W, v_P, v_N, v_S;
	double eP, uP, uE, uW, vP, vN, vS;

	if(gridType=="staggered")
	{
		for(int i=0; i<rowNo; i++)
		{
			for(int j=0; j<colNo; j++)
			{
				u_P=idU[i][j]-1;
				u_E=idU[i][j+1]-1;
				v_P=idV[i][j]-1;
				v_S=idV[i+1][j]-1;

				for(int k=0; k<aisleNo; k++)
				{
					uP=uField[u_P][k];
					uE=uField[u_E][k];
					vP=vField[v_P][k];
					vS=vField[v_S][k];
					eP=(uE-uP)/dx+(vP-vS)/dy;
					strain3DField[i][j][k]=eP;
				}
			}
		}
	}
	else
	{
		for(int i=0; i<rowNo; i++)
		{
			for(int j=0; j<colNo; j++)
			{
				if(i==0 || i==rowNo-1) v_P=idV[i][j]-1;
				if(i>0) v_N=idV[i-1][j]-1;
				if(i<rowNo-1) v_S=idV[i+1][j]-1;

				if(j==0 || j==colNo-1) u_P=idU[i][j]-1;
				if(j>0) u_W=idU[i][j-1]-1;
				if(j<colNo-1) u_E=idU[i][j+1]-1;

				for(int k=0; k<aisleNo; k++)
				{
					eP=0;

					if(i==0 || i==rowNo-1) vP=vField[v_P][k];
					if(i>0) vN=vField[v_N][k];
					if(i<rowNo-1) vS=vField[v_S][k];

					if(j==0 || j==colNo-1) uP=uField[u_P][k];
					if(j>0) uW=uField[u_W][k];
					if(j<colNo-1) uE=uField[u_E][k];

					if(i==0) eP+=(vP-vS)/dy;
					else if(i==rowNo-1) eP+=(vN-vP)/dy;
					else eP+=(vN-vS)/(2*dy);

					if(j==0) eP+=(uE-uP)/dx;
					else if(j==colNo-1) eP+=(uP-uW)/dx;
					else eP+=(uE-uW)/(2*dx);
					
					strain3DField[i][j][k]=eP;
				}
			}
		}
	}
	return;
}

void dataProcessing::export1DFieldToTxt(vector<double> my1DField, string fieldName)
{
	int rowNo=my1DField.size();
	string fileName="../export/"+fieldName+"_"+gridType+"-grid.txt";

	ofstream myFile(fileName);
	if(myFile.is_open())
	{
		for(int i=0; i<rowNo; i++)
		{
			myFile << my1DField[i];
			myFile << "\n";
		}

		myFile.close();
	}

	return;
}

void dataProcessing::exportSealedColumnNumericalSolution(double dy, double dt, double Ly,
	int timeStep, string pairName)
{
	string fieldName;
	vector<double> yCoordP;yCoordP.resize(pressure3DField.size());
	vector<double> pField;pField.resize(pressure3DField.size());
	vector<double> yCoordV;yCoordV.resize(vDisplacement3DField.size());
	vector<double> vField;vField.resize(vDisplacement3DField.size());

	// Exports yCoordP
	if(gridType=="staggered") yCoordP[0]=Ly-dy/2;
	else yCoordP[0]=Ly;
	for(int i=1; i<yCoordP.size(); i++) yCoordP[i]=yCoordP[i-1]-dy;
	fieldName="sealedColumn_"+pairName+"_YPNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(yCoordP,fieldName);

	// Exports pField
	if(pressure3DField[0].size()%2==0)
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=pressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[i][midCols][timeStep]/2;
			pField[i]+=pressure3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=pressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[i][midCols][timeStep];
		}
	fieldName="sealedColumn_"+pairName+"_PNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(pField,fieldName);

	// Exports yCoordV
	for(int i=0; i<yCoordV.size(); i++) yCoordV[i]=Ly-i*dy;
	fieldName="sealedColumn_"+pairName+"_YVNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(yCoordV,fieldName);

	// Exports vField
	if(vDisplacement3DField[0].size()%2==2)
		for(int i=0; i<vField.size(); i++)
		{
			int midCols=vDisplacement3DField[0].size()/2;
			vField[i]=0;
			vField[i]+=vDisplacement3DField[i][midCols][timeStep]/2;
			vField[i]+=vDisplacement3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<vField.size(); i++)
		{
			int midCols=vDisplacement3DField[0].size()/2;
			vField[i]=0;
			vField[i]+=vDisplacement3DField[i][midCols][timeStep];
		}
	fieldName="sealedColumn_"+pairName+"_VNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(vField,fieldName);

	return;
}

void dataProcessing::exportSealedColumnAnalyticalSolution(double Ly, double alpha, double Q,
	double rho, double g, double rho_f, double M, double sigmab, double dt, int timeStep, double c,
	string pairName)
{
	int pointsExact=5000;
	double yValue;
	double dy=Ly/pointsExact;
	vector<double> yExact;
	vector<double> vExact;
	vector<double> pExact;
	string fieldName;
	double tValue=dt*timeStep;
	double sigma0;

	yExact.resize(pointsExact+1);
	vExact.resize(pointsExact+1);
	pExact.resize(pointsExact+1);

	yExact[0]=Ly;
	for(int i=1; i<pointsExact+1; i++) yExact[i]=yExact[i-1]-dy;
	
	for(int i=0; i<pointsExact+1; i++)
	{
		vExact[i]=-(rho-alpha*rho_f)*g/(2*M)*yExact[i]*yExact[i];
		vExact[i]+=(sigmab+0.5*rho*g*Ly)/(M+alpha*alpha*Q)*yExact[i];
		vExact[i]+=(rho-alpha*rho_f)*g*Ly/(2*M)*yExact[i];
		pExact[i]=-alpha*Q/(M+alpha*alpha*Q)*(sigmab+0.5*rho*g*Ly);
		pExact[i]+=rho_f*g*(yExact[i]-0.5*Ly);
	}

	fieldName="sealedColumn_"+pairName+"_YExact";
	export1DFieldToTxt(yExact,fieldName);

	fieldName="sealedColumn_"+pairName+"_VExact";
	export1DFieldToTxt(vExact,fieldName);

	fieldName="sealedColumn_"+pairName+"_PExact";
	export1DFieldToTxt(pExact,fieldName);

	ofstream myFile("../export/solveSealedColumnRunInfo.txt",fstream::app);
	if(myFile.is_open())
	{
		myFile << timeStep << "\n";

		myFile.close();
	}

	return;
}

void dataProcessing::exportTerzaghiNumericalSolution(double dy, double dt, double Ly, int timeStep,
	string pairName)
{
	string fieldName;
	vector<double> yCoordP;yCoordP.resize(pressure3DField.size());
	vector<double> pField;pField.resize(pressure3DField.size());
	vector<double> yCoordE;yCoordE.resize(strain3DField.size());
	vector<double> eField;eField.resize(strain3DField.size());
	vector<double> yCoordV;yCoordV.resize(vDisplacement3DField.size());
	vector<double> vField;vField.resize(vDisplacement3DField.size());

	// Exports yCoordP
	if(gridType=="staggered") yCoordP[0]=Ly-dy/2;
	else yCoordP[0]=Ly;
	for(int i=1; i<yCoordP.size(); i++) yCoordP[i]=yCoordP[i-1]-dy;
	fieldName="terzaghi_"+pairName+"_YPNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	if(gridType=="staggered") yCoordP.insert(yCoordP.begin(),Ly);
	export1DFieldToTxt(yCoordP,fieldName);

	// Exports pField
	if(pressure3DField[0].size()%2==0)
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=pressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[i][midCols][timeStep]/2;
			pField[i]+=pressure3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=pressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[i][midCols][timeStep];
		}
	fieldName="terzaghi_"+pairName+"_PNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	if(gridType=="staggered") pField.insert(pField.begin(),0);
	export1DFieldToTxt(pField,fieldName);

	// Exports yCoordE
	if(gridType=="staggered") yCoordE[0]=Ly-dy/2;
	else yCoordE[0]=Ly;
	for(int i=1; i<yCoordE.size(); i++) yCoordE[i]=yCoordE[i-1]-dy;
	fieldName="terzaghi_"+pairName+"_YENumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(yCoordE,fieldName);

	// Exports eField
	if(strain3DField[0].size()%2==0)
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=strain3DField[0].size()/2;
			eField[i]=0;
			eField[i]+=strain3DField[i][midCols][timeStep]/2;
			eField[i]+=strain3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<eField.size(); i++)
		{
			int midCols=strain3DField[0].size()/2;
			eField[i]=0;
			eField[i]+=strain3DField[i][midCols][timeStep];
		}
	fieldName="terzaghi_"+pairName+"_ENumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(eField,fieldName);

	// Exports yCoordV
	for(int i=0; i<yCoordV.size(); i++) yCoordV[i]=Ly-i*dy;
	fieldName="terzaghi_"+pairName+"_YVNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(yCoordV,fieldName);

	// Exports vField
	if(vDisplacement3DField[0].size()%2==2)
		for(int i=0; i<vField.size(); i++)
		{
			int midCols=vDisplacement3DField[0].size()/2;
			vField[i]=0;
			vField[i]+=vDisplacement3DField[i][midCols][timeStep]/2;
			vField[i]+=vDisplacement3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<vField.size(); i++)
		{
			int midCols=vDisplacement3DField[0].size()/2;
			vField[i]=0;
			vField[i]+=vDisplacement3DField[i][midCols][timeStep];
		}
	fieldName="terzaghi_"+pairName+"_VNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(vField,fieldName);

	return;
}

double dataProcessing::terzaghiPAnalyticalSolution(double yValue, double tValue, double L,
	double sigmab, double M, double alpha, double Q, double rho, double g, double rho_f, double c)
{
	double pValue;
	double As, P0;

	As=sigmab/(M+alpha*alpha*Q)+(rho*g*L)/(2*(M+alpha*alpha*Q))+0.5*(rho-alpha*rho_f)*(g*L/M);
	P0=0.5*rho_f*g*L+alpha*Q*(rho-alpha*rho_f)*g*L/(2*M)-alpha*Q*As;

	double pi=4*atan(1);
	double value, sum=0;

	for(int j=1; j<301; j++)
	{
		value=cos((2*j-1)*(pi*yValue)/(2*L));
		value=value*exp(-(2*j-1)*(2*j-1)*(c*pi*pi*tValue)/(4*L*L));
		value=value/(2*j-1);
		value=value*pow(-1,j+1);
		sum=sum+value;			
	}

	pValue=((4*P0)/pi)*sum;
	pValue-=rho_f*g*(L-yValue);

	return pValue;
}

double dataProcessing::terzaghiVAnalyticalSolution(double yValue, double tValue, double L,
	double sigmab, double M, double alpha, double Q, double rho, double g, double rho_f, double c)
{	
	double vValue;
	double pi=4*atan(1);
	double value, sum=0;
	double As, P0;

	As=sigmab/(M+alpha*alpha*Q)+(rho*g*L)/(2*(M+alpha*alpha*Q))+0.5*(rho-alpha*rho_f)*(g*L/M);
	P0=0.5*rho_f*g*L+alpha*Q*(rho-alpha*rho_f)*g*L/(2*M)-alpha*Q*As;

	for(int j=1; j<301; j++)
	{
		value=sin((2*j-1)*(pi*yValue)/(2*L));
		value=value*exp(-(2*j-1)*(2*j-1)*(c*pi*pi*tValue)/(4*L*L));
		value=value/((2*j-1)*(2*j-1));
		value=value*pow(-1,j+1);
		sum=sum+value;			
	}

	vValue=(8*alpha*L*P0)/(pi*pi*M)*sum;
	vValue+=sigmab*yValue/M;
	vValue+=(g/M)*(rho-alpha*rho_f)*(L*yValue-0.5*yValue*yValue);

	return vValue;
}

void dataProcessing::exportTerzaghiAnalyticalSolution(double Ly, double alpha, double Q, double rho,
	double g, double rho_f, double M, double sigmab, double dt, int timeStep, double c,
	string pairName)
{
	int pointsExact=5000;
	double yValue;
	double dy=Ly/pointsExact;
	vector<double> yExact;
	vector<double> vExact;
	vector<double> pExact;
	string fieldName;
	double tValue=dt*timeStep;

	yExact.resize(pointsExact+1);
	vExact.resize(pointsExact+1);
	pExact.resize(pointsExact+1);

	yExact[0]=Ly;
	for(int i=1; i<pointsExact+1; i++) yExact[i]=yExact[i-1]-dy;
	
	for(int i=0; i<pointsExact+1; i++)
	{
		vExact[i]=terzaghiVAnalyticalSolution(yExact[i],tValue,Ly,sigmab,M,alpha,Q,rho,g,rho_f,c);
		pExact[i]=terzaghiPAnalyticalSolution(yExact[i],tValue,Ly,sigmab,M,alpha,Q,rho,g,rho_f,c);
	}

	fieldName="terzaghi_"+pairName+"_YExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(yExact,fieldName);

	fieldName="terzaghi_"+pairName+"_VExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(vExact,fieldName);

	fieldName="terzaghi_"+pairName+"_PExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(pExact,fieldName);

	ofstream myFile("../export/solveTerzaghiRunInfo.txt",fstream::app);
	if(myFile.is_open())
	{
		myFile << timeStep << "\n";

		myFile.close();
	}

	return;
}

void dataProcessing::getTerzaghiErrorNorm(double dy, double dt, double h, double Ly, double P0,
	double c, double alpha, double M, double sigmab, double Q, double rho, double g, double rho_f)
{
	myErrorNorm.p=0;
	myErrorNorm.v=0;
	int rowNo=pressure3DField.size();
	int colNo=pressure3DField[0].size();
	double time=(pressure3DField[0][0].size()-1)*dt;
	double pExact;
	double pNumeric;
	double vExact;
	double vNumeric;
	double yValue;
	double dividend, divisor;

	// Calculates pressure error norm
	if(gridType=="staggered") yValue=Ly-dy/2;
	else yValue=Ly;

	dividend=0;
	divisor=0;
	for(int i=0; i<rowNo; i++)
	{
		for(int j=0; j<colNo; j++)
		{
			pExact=terzaghiPAnalyticalSolution(yValue,time,Ly,sigmab,M,alpha,Q,rho,g,rho_f,c);
			pNumeric=pressure3DField[i][j][aisleNo-1];
			dividend+=pow((pExact-pNumeric)*h,2);
			divisor+=pow(h,2);
		}

		yValue=yValue-dy;
	}

	myErrorNorm.p=sqrt(dividend/divisor);

	// Calculates displacement error norm
	yValue=Ly;
	if(gridType=="staggered") rowNo=rowNo+1;

	dividend=0;
	divisor=0;
	for(int i=0; i<rowNo; i++)
	{
		for(int j=0; j<colNo; j++)
		{
			vExact=terzaghiVAnalyticalSolution(yValue,time,Ly,sigmab,M,alpha,Q,rho,g,rho_f,c);
			vNumeric=vDisplacement3DField[i][j][aisleNo-1];
			dividend+=pow((vExact-vNumeric)*h,2);
			divisor+=pow(h,2);
		}

		yValue=yValue-dy;
	}

	myErrorNorm.v=sqrt(dividend/divisor);

	saveTerzaghiErrorNorm(dt,h,myErrorNorm.p,myErrorNorm.v);

	return;
}

void dataProcessing::saveTerzaghiErrorNorm(double dt, double h, double pError, double vError)
{
	string fileName="../export/terzaghiErrorNorm_dt="+to_string(dt)+"_"+gridType+"-grid.txt";

	ofstream myFile(fileName,fstream::app);
	if(myFile.is_open())
	{
		myFile << h << "\t" << pError << "\t" << vError << "\n";

		myFile.close();
	}

	fileName="../export/terzaghiErrorNorm_h="+to_string(h)+"_"+gridType+"-grid.txt";

	ofstream myFile2(fileName,fstream::app);
	if(myFile2.is_open())
	{
		myFile2 << dt << "\t" << pError << "\t" << vError << "\n";

		myFile2.close();
	}

	return;
}

void dataProcessing::exportMandelNumericalSolution(double dx, double dy, double dt, double Lx,
	double Ly, int timeStep, string pairName)
{
	string fieldName;
	vector<double> xCoordP;xCoordP.resize(pressure3DField[0].size());
	vector<double> pField;pField.resize(pressure3DField[0].size());
	vector<double> xCoordE;xCoordE.resize(strain3DField[0].size());
	vector<double> eField;eField.resize(strain3DField[0].size());
	vector<double> xCoordU;xCoordU.resize(uDisplacement3DField[0].size());
	vector<double> uField;uField.resize(uDisplacement3DField[0].size());
	vector<double> yCoordV;yCoordV.resize(vDisplacement3DField.size());
	vector<double> vField;vField.resize(vDisplacement3DField.size());
	
	// Exports xCoordP
	if(gridType=="staggered") xCoordP[0]=dx/2;
	else xCoordP[0]=0;
	for(int i=1; i<xCoordP.size(); i++) xCoordP[i]=xCoordP[i-1]+dx;
	fieldName="mandel_"+pairName+"_XPNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	if(gridType=="staggered") xCoordP.insert(xCoordP.end(),Lx);
	export1DFieldToTxt(xCoordP,fieldName);

	// Exports pField
	if(pressure3DField.size()%2==0)
		for(int i=0; i<pField.size(); i++)
		{
			int midRows=pressure3DField.size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[midRows][i][timeStep]/2;
			pField[i]+=pressure3DField[midRows-1][i][timeStep]/2;
		}
	else
		for(int i=0; i<pField.size(); i++)
		{
			int midRows=pressure3DField.size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[midRows][i][timeStep];
		}
	fieldName="mandel_"+pairName+"_PNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	if(gridType=="staggered") pField.insert(pField.end(),0);
	export1DFieldToTxt(pField,fieldName);

	// Exports xCoordE
	if(gridType=="staggered") xCoordE[0]=dx/2;
	else xCoordE[0]=0;
	for(int i=1; i<xCoordE.size(); i++) xCoordE[i]=xCoordE[i-1]+dx;
	fieldName="mandel_"+pairName+"_XENumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(xCoordE,fieldName);

	// Exports eField
	if(strain3DField.size()%2==0)
		for(int i=0; i<eField.size(); i++)
		{
			int midRows=strain3DField.size()/2;
			eField[i]=0;
			eField[i]+=strain3DField[midRows][i][timeStep]/2;
			eField[i]+=strain3DField[midRows-1][i][timeStep]/2;
		}
	else
		for(int i=0; i<eField.size(); i++)
		{
			int midRows=strain3DField.size()/2;
			eField[i]=0;
			eField[i]+=strain3DField[midRows][i][timeStep];
		}
	fieldName="mandel_"+pairName+"_ENumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(eField,fieldName);
	
	// Exports xCoordU
	for(int i=0; i<xCoordU.size(); i++) xCoordU[i]=i*dx;
	fieldName="mandel_"+pairName+"_XUNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(xCoordU,fieldName);
	
	// Exports uField
	if(uDisplacement3DField.size()%2==2)
		for(int i=0; i<=uField.size(); i++)
		{
			int midRows=uDisplacement3DField.size()/2;
			uField[i]=0;
			uField[i]+=uDisplacement3DField[midRows][i][timeStep]/2;
			uField[i]+=uDisplacement3DField[midRows-1][i][timeStep]/2;
		}
	else
		for(int i=0; i<uField.size(); i++)
		{
			int midRows=uDisplacement3DField.size()/2;
			uField[i]=0;
			uField[i]+=uDisplacement3DField[midRows][i][timeStep];
		}
	fieldName="mandel_"+pairName+"_UNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(uField,fieldName);

	// Exports yCoordV
	for(int i=0; i<yCoordV.size(); i++) yCoordV[i]=Ly-i*dy;
	fieldName="mandel_"+pairName+"_YVNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(yCoordV,fieldName);

	// Exports vField
	if(vDisplacement3DField[0].size()%2==2)
		for(int i=0; i<vField.size(); i++)
		{
			int midCols=vDisplacement3DField[0].size()/2;
			vField[i]=0;
			vField[i]+=vDisplacement3DField[i][midCols][timeStep]/2;
			vField[i]+=vDisplacement3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<vField.size(); i++)
		{
			int midCols=vDisplacement3DField[0].size()/2;
			vField[i]=0;
			vField[i]+=vDisplacement3DField[i][midCols][timeStep];
		}
	fieldName="mandel_"+pairName+"_VNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(vField,fieldName);

	return;
}

void dataProcessing::findMandelRoots(double P0, double F, double L, double alpha, double M,
	double lambda, double Q)
{
	int i, j;
	int mandelRootsNumber=500;
	int maxIterations=500;
	double tolerance=1e-12;
	double residual, lowerResidual, upperResidual;
	double pi=4*atan(1);
	double root, lowerBound, upperBound;

	mandelTransCoef=(F*M)/(P0*alpha*(lambda-M));

	for(i=0; i<mandelRootsNumber; i++)
	{
		lowerBound=i*pi+pi/4;
		upperBound=lowerBound+pi/2-1e-5;
		residual=1;
		j=0;

		while(fabs(residual)>tolerance && j<maxIterations)
		{
			lowerResidual=tan(lowerBound)-mandelTransCoef*lowerBound;
			upperResidual=tan(upperBound)-mandelTransCoef*upperBound;
			root=(lowerBound+upperBound)/2;
			residual=tan(root)-mandelTransCoef*root;
			if(residual*lowerResidual>0) lowerBound=root;
			else upperBound=root;
			j++;
		}

		mandelRoots.push_back(root);
	}

	return;
}

double dataProcessing::mandelPAnalyticalSolution(double xValue, double tValue, double c, double P0,
	double L)
{
	double pValue=0;
	double dividend, divisor;
	double xi;

	for(int i=0; i<mandelRoots.size(); i++)
	{
		xi=mandelRoots[i]/L;
		dividend=cos(xi*xValue)-cos(xi*L);
		dividend=dividend*sin(xi*L);
		dividend=dividend*exp(-xi*xi*c*tValue);
		divisor=xi*L;
		divisor=divisor-sin(xi*L)*cos(xi*L);
		pValue=pValue+dividend/divisor;
	}

	pValue=pValue*2*P0;

	return pValue;
}

double dataProcessing::mandelEAnalyticalSolution(double xValue, double tValue, double c, double F,
	double L, double alpha, double Q, double M, double lambda, double P0)
{
	double eValue=0;
	double dividend, divisor;
	double xi;

	for(int i=0; i<mandelRoots.size(); i++)
	{
		xi=mandelRoots[i]/L;
		dividend=(alpha/M)*cos(xi*xValue)+(1/(alpha*Q))*cos(xi*L);
		dividend=dividend*sin(xi*L);
		dividend=dividend*exp(-xi*xi*c*tValue);
		divisor=xi*L;
		divisor=divisor-sin(xi*L)*cos(xi*L);
		eValue=eValue+dividend/divisor;
	}

	eValue=eValue*2*P0+F/(M+lambda);

	return eValue;
}

double dataProcessing::mandelUAnalyticalSolution(double xValue, double tValue, double c, double F,
	double L, double alpha, double Q, double G, double M, double lambda, double P0)
{
	double uValue=0;
	double dividend, divisor;
	double xi;

	for(int i=0; i<mandelRoots.size(); i++)
	{
		xi=mandelRoots[i]/L;
		dividend=(alpha/(M*xi))*sin(xi*xValue);
		dividend=dividend-((lambda+alpha*alpha*Q)/(2*alpha*G*Q))*cos(xi*L)*xValue;
		dividend=dividend*sin(xi*L)*exp(-xi*xi*c*tValue);
		divisor=xi*L;
		divisor=divisor-sin(xi*L)*cos(xi*L);
		uValue=uValue+dividend/divisor;
	}

	uValue=uValue*2*P0;
	uValue=uValue-((F*lambda)/(2*G*(M+lambda)))*xValue;

	return uValue;
}

double dataProcessing::mandelVAnalyticalSolution(double yValue, double tValue, double F, double M,
	double lambda, double L, double K, double mu_f, double c, double alpha, double P0)
{
	double vValue=0;
	double dividend, divisor;
	double xi;
	double G=(M-lambda)/2;

	for(int i=0; i<mandelRoots.size(); i++)
	{
		xi=mandelRoots[i]/L;
		dividend=sin(xi*L)*cos(xi*L)*exp(-xi*xi*c*tValue);
		divisor=xi*L-sin(xi*L)*cos(xi*L);
		vValue=vValue+dividend/divisor;
	}

	vValue=vValue*((K*M*P0)/(alpha*c*G*mu_f));
	vValue=vValue+((F*M)/(2*G*(M+lambda)));
	vValue=vValue*yValue;

	return vValue;
}

void dataProcessing::exportMandelAnalyticalSolution(double L, double H, double c, double P0,
	double alpha, double Q, double M, double lambda, double F, double K, double mu_f, double dt,
	int timeStep, string pairName)
{
	int pointsExact=5000;
	double xValue;
	double tValue=timeStep*dt;
	double dx=L/pointsExact;
	double dy=H/pointsExact;
	vector<double> xExact;
	vector<double> yExact;
	vector<double> uExact;
	vector<double> vExact;
	vector<double> pExact;
	vector<double> eExact;
	string fieldName;
	double G=(M-lambda)/2;

	xExact.resize(pointsExact+1);
	yExact.resize(pointsExact+1);
	uExact.resize(pointsExact+1);
	vExact.resize(pointsExact+1);
	pExact.resize(pointsExact+1);
	eExact.resize(pointsExact+1);

	xExact[0]=0;
	for(int i=1; i<pointsExact+1; i++) xExact[i]=xExact[i-1]+dx;

	yExact[0]=0;
	for(int i=1; i<pointsExact+1; i++) yExact[i]=yExact[i-1]+dy;
	
	for(int i=0; i<pointsExact+1; i++)
	{
		pExact[i]=mandelPAnalyticalSolution(xExact[i],tValue,c,P0,L);
		eExact[i]=mandelEAnalyticalSolution(xExact[i],tValue,c,F,L,alpha,Q,M,lambda,P0);
		uExact[i]=mandelUAnalyticalSolution(xExact[i],tValue,c,F,L,alpha,Q,G,M,lambda,P0);
		vExact[i]=mandelVAnalyticalSolution(yExact[i],tValue,F,M,lambda,L,K,mu_f,c,alpha,P0);
	}

	fieldName="mandel_"+pairName+"_XExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(xExact,fieldName);

	fieldName="mandel_"+pairName+"_YExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(yExact,fieldName);

	fieldName="mandel_"+pairName+"_UExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(uExact,fieldName);

	fieldName="mandel_"+pairName+"_VExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(vExact,fieldName);

	fieldName="mandel_"+pairName+"_PExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(pExact,fieldName);

	fieldName="mandel_"+pairName+"_EExact_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	export1DFieldToTxt(eExact,fieldName);

	ofstream myFile("../export/solveMandelRunInfo.txt",fstream::app);
	if(myFile.is_open())
	{
		myFile << timeStep << "\n";

		myFile.close();
	}

	return;
}

void dataProcessing::storeMacroPressure3DField(vector<vector<int>> idP,
	vector<vector<double>> pField)
{	
	int rowNo=idP.size();
	int colNo=idP[0].size();
	int position;
	double value;

	for(int i=0; i<rowNo; i++)
	{
		for(int j=0; j<colNo; j++)
		{
			position=idP[i][j]-1;

			for(int k=0; k<aisleNo; k++)
			{
				value=pField[position][k];
				macroPressure3DField[i][j][k]=value;
			}
		}
	}

	return;
}

void dataProcessing::exportMacroPressureHSolution(double dy, double h, double Ly, int timeStep,
	string pairName)
{
	string fieldName;
	vector<double> yCoordP;yCoordP.resize(pressure3DField.size());
	vector<double> pField;pField.resize(pressure3DField.size());
	vector<double> pMField;pMField.resize(macroPressure3DField.size());

	// Exports yCoordP
	if(gridType=="staggered") yCoordP[0]=Ly-dy/2;
	else yCoordP[0]=Ly;
	for(int i=1; i<yCoordP.size(); i++) yCoordP[i]=yCoordP[i-1]-dy;
	fieldName="terzaghi_"+pairName+"_YPNumeric_h="+to_string(h)+"_timeStep="+to_string(timeStep);
	// if(gridType=="staggered") yCoordP.insert(yCoordP.begin(),Ly);
	export1DFieldToTxt(yCoordP,fieldName);

	// Exports pField
	if(pressure3DField[0].size()%2==0)
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=pressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[i][midCols][timeStep]/2;
			pField[i]+=pressure3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=pressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[i][midCols][timeStep];
		}
	fieldName="terzaghi_"+pairName+"_PNumeric_h="+to_string(h)+"_timeStep="+to_string(timeStep);
	// if(gridType=="staggered") pField.insert(pField.begin(),0);
	export1DFieldToTxt(pField,fieldName);

	// Exports pMField
	if(macroPressure3DField[0].size()%2==0)
		for(int i=0; i<pMField.size(); i++)
		{
			int midCols=macroPressure3DField[0].size()/2;
			pMField[i]=0;
			pMField[i]+=macroPressure3DField[i][midCols][timeStep]/2;
			pMField[i]+=macroPressure3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<pMField.size(); i++)
		{
			int midCols=macroPressure3DField[0].size()/2;
			pMField[i]=0;
			pMField[i]+=macroPressure3DField[i][midCols][timeStep];
		}
	fieldName="terzaghi_"+pairName+"_MacroPNumeric_h="+to_string(h)+"_timeStep="+
		to_string(timeStep);
	// if(gridType=="staggered") pMField.insert(pMField.begin(),0);
	export1DFieldToTxt(pMField,fieldName);

	return;
}

void dataProcessing::exportMacroPressureTSolution(double dy, double dt, double Ly, int timeStep,
	string pairName)
{
	string fieldName;
	vector<double> yCoordP;yCoordP.resize(pressure3DField.size());
	vector<double> pField;pField.resize(pressure3DField.size());
	vector<double> pMField;pMField.resize(macroPressure3DField.size());

	// Exports yCoordP
	if(gridType=="staggered") yCoordP[0]=Ly-dy/2;
	else yCoordP[0]=Ly;
	for(int i=1; i<yCoordP.size(); i++) yCoordP[i]=yCoordP[i-1]-dy;
	fieldName="terzaghi_"+pairName+"_YPNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	// if(gridType=="staggered") yCoordP.insert(yCoordP.begin(),Ly);
	export1DFieldToTxt(yCoordP,fieldName);

	// Exports pField
	if(pressure3DField[0].size()%2==0)
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=pressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[i][midCols][timeStep]/2;
			pField[i]+=pressure3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=pressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=pressure3DField[i][midCols][timeStep];
		}
	fieldName="terzaghi_"+pairName+"_PNumeric_dt="+to_string(dt)+"_timeStep="+to_string(timeStep);
	// if(gridType=="staggered") pField.insert(pField.begin(),0);
	export1DFieldToTxt(pField,fieldName);

	// Exports pMField
	if(macroPressure3DField[0].size()%2==0)
		for(int i=0; i<pMField.size(); i++)
		{
			int midCols=macroPressure3DField[0].size()/2;
			pMField[i]=0;
			pMField[i]+=macroPressure3DField[i][midCols][timeStep]/2;
			pMField[i]+=macroPressure3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<pMField.size(); i++)
		{
			int midCols=macroPressure3DField[0].size()/2;
			pMField[i]=0;
			pMField[i]+=macroPressure3DField[i][midCols][timeStep];
		}
	fieldName="terzaghi_"+pairName+"_MacroPNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	// if(gridType=="staggered") pMField.insert(pMField.begin(),0);
	export1DFieldToTxt(pMField,fieldName);

	return;
}

void dataProcessing::exportStripfootSolution(double dx, double dy, double dt, double Ly,
	int timeStep, string pairName)
{
	string fileName;
	vector<vector<double>> xCoord, yCoord;
	double position;
	int rowNo=pressure3DField.size();
	int colNo=pressure3DField[0].size();

	fileName="../export/stripfoot_"+pairName+"xCoord_"+gridType+"-grid.txt";
	ofstream xCoordFile(fileName);
	if(xCoordFile.is_open())
	{
		for(int i=0; i<rowNo; i++)
		{
			if(gridType=="staggered") position=dx/2;
			else position=0;

			for(int j=0; j<colNo; j++)
			{
				xCoordFile << position;
				xCoordFile << "\t";
				position+=dx;
			}

			xCoordFile << "\n";
		}

		xCoordFile.close();
	}

	fileName="../export/stripfoot_"+pairName+"yCoord_"+gridType+"-grid.txt";
	ofstream yCoordFile(fileName);
	if(yCoordFile.is_open())
	{
		if(gridType=="staggered") position=Ly-dy/2;
		else position=Ly;

		for(int i=0; i<rowNo; i++)
		{
			for(int j=0; j<colNo; j++)
			{
				yCoordFile << position;
				yCoordFile << "\t";
			}

			yCoordFile << "\n";
			position-=dy;
		}

		yCoordFile.close();
	}

	fileName="../export/stripfoot_"+pairName+"_PNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep)+"_"+gridType+"-grid.txt";
	ofstream pFile(fileName);
	fileName="../export/stripfoot_"+pairName+"_MacroPNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep)+"_"+gridType+"-grid.txt";
	ofstream pMFile(fileName);
	if(pFile.is_open() && pMFile.is_open())
	{
		for(int i=0; i<rowNo; i++)
		{
			for(int j=0; j<colNo; j++)
			{
				pFile << pressure3DField[i][j][timeStep];
				pFile << "\t";
				pMFile << macroPressure3DField[i][j][timeStep];
				pMFile << "\t";
			}

			pFile << "\n";
			pMFile << "\n";
		}

		pFile.close();
		pMFile.close();
	}

	return;
}