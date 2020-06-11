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

class doubleDataProcessing
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
	void storeMacroPressure3DField(vector<vector<int>>,vector<vector<double>>);
	void exportSealedDoubleNumericalSolution(double,double,double,int,string);
	void exportSealedDoubleAnalyticalSolution(double,double,double,double,double,double,double,
		double,double,int,string);
	void exportDrainedDoubleNumericalSolution(double,double,double,int,string);
	double storagePhi1AnalyticalSolution(double,double,double,double,double,double,double,double,
		double);
	double storagePhi2AnalyticalSolution(double,double,double,double,double,double,double,double,
		double);
	void exportStorageAnalyticalSolution(double,double,double,double,double,double,double,double,
		double,double,double,double,int,string);
	double leakingAnalyticalSolution(double,double,double,double,double,double,double,double,
		double,double);
	void exportLeakingAnalyticalSolution(double,double,double,double,double,double,double,double,
		double,double,double,double,double,int,string);

	// Constructor
	doubleDataProcessing(vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,
		vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,string,string,double,
		double);

	// Destructor
	~doubleDataProcessing();
};

doubleDataProcessing::doubleDataProcessing(vector<vector<int>> idU, vector<vector<int>> idV,
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

doubleDataProcessing::~doubleDataProcessing(){}

void doubleDataProcessing::resize3DFields(vector<vector<int>> idU, vector<vector<int>> idV,
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

void doubleDataProcessing::storeUDisplacem3DField(vector<vector<int>> idU,
	vector<vector<double>> uField)
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

void doubleDataProcessing::storeVDisplacem3DField(vector<vector<int>> idV,
	vector<vector<double>> vField)
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

void doubleDataProcessing::storePressure3DField(vector<vector<int>> idP, 
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
				pressure3DField[i][j][k]=value;
			}
		}
	}

	return;
}

void doubleDataProcessing::storeStrain3DField(vector<vector<int>> idU, vector<vector<int>> idV,
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

void doubleDataProcessing::export1DFieldToTxt(vector<double> my1DField, string fieldName)
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

void doubleDataProcessing::storeMacroPressure3DField(vector<vector<int>> idP,
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

void doubleDataProcessing::exportSealedDoubleNumericalSolution(double dy, double dt, double Ly,
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
	fieldName="sealedDouble_"+pairName+"_YPNumeric_dt="+to_string(dt)+"_timeStep="+
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
	fieldName="sealedDouble_"+pairName+"_PPoreNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(pField,fieldName);

	// Exports pField
	if(macroPressure3DField[0].size()%2==0)
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=macroPressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=macroPressure3DField[i][midCols][timeStep]/2;
			pField[i]+=macroPressure3DField[i][midCols-1][timeStep]/2;
		}
	else
		for(int i=0; i<pField.size(); i++)
		{
			int midCols=macroPressure3DField[0].size()/2;
			pField[i]=0;
			pField[i]+=macroPressure3DField[i][midCols][timeStep];
		}
	fieldName="sealedDouble_"+pairName+"_PFracNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(pField,fieldName);

	// Exports yCoordV
	for(int i=0; i<yCoordV.size(); i++) yCoordV[i]=Ly-i*dy;
	fieldName="sealedDouble_"+pairName+"_YVNumeric_dt="+to_string(dt)+"_timeStep="+
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
	fieldName="sealedDouble_"+pairName+"_VNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(vField,fieldName);

	return;
}

void doubleDataProcessing::exportSealedDoubleAnalyticalSolution(double Ly, double alpha1,
	double alpha2, double M, double S11, double S12, double S22, double sigmab, double dt,
	int timeStep, string pairName)
{
	int pointsExact=5000;
	double yValue;
	double dy=Ly/pointsExact;
	vector<double> yExact;
	vector<double> vExact;
	vector<double> pExact;
	string fieldName;

	yExact.resize(pointsExact+1);
	vExact.resize(pointsExact+1);
	pExact.resize(pointsExact+1);

	yExact[0]=Ly;
	for(int i=1; i<pointsExact+1; i++) yExact[i]=yExact[i-1]-dy;
	
	for(int i=0; i<pointsExact+1; i++)
	{
		pExact[i]=-((alpha1+alpha2)/(M*(S11+2*S12+S22)+(alpha1+alpha2)*(alpha1+alpha2)))*sigmab;
	}

	fieldName="sealedDouble_"+pairName+"_YExact";
	export1DFieldToTxt(yExact,fieldName);

	fieldName="sealedDouble_"+pairName+"_VExact";
	export1DFieldToTxt(vExact,fieldName);

	fieldName="sealedDouble_"+pairName+"_PExact";
	export1DFieldToTxt(pExact,fieldName);

	ofstream myFile("../export/solveSealedDoubleRunInfo.txt",fstream::app);
	if(myFile.is_open())
	{
		myFile << timeStep << "\n";

		myFile.close();
	}

	return;
}

void doubleDataProcessing::exportDrainedDoubleNumericalSolution(double dy, double dt, double Ly,
	int timeStep, string pairName)
{
	string fieldName;
	vector<double> yCoordP;yCoordP.resize(pressure3DField.size());
	vector<double> pField;pField.resize(pressure3DField.size());
	vector<double> pMField;pMField.resize(macroPressure3DField.size());

	// Exports yCoordP
	if(gridType=="staggered") yCoordP[0]=Ly-dy/2;
	else yCoordP[0]=Ly;
	for(int i=1; i<yCoordP.size(); i++) yCoordP[i]=yCoordP[i-1]-dy;
	fieldName="drainedDouble_"+pairName+"_YPNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
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
	fieldName="drainedDouble_"+pairName+"_PPoreNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	if(gridType=="staggered") pField.insert(pField.begin(),0);
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
	fieldName="drainedDouble_"+pairName+"_PFracNumeric_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	if(gridType=="staggered") pMField.insert(pMField.begin(),0);
	export1DFieldToTxt(pMField,fieldName);

	return;
}

double doubleDataProcessing::storagePhi1AnalyticalSolution(double yValue, double tValue, double L,
	double p0, double b1, double b2, double c1, double c2, double q1)
{
	double phiValue;
	double pi=4*atan(1);
	double value, sum=0;

	for(int j=1; j<301; j++)
	{
		value=cos((2*j-1)*(pi*yValue)/(2*L));
		value=value*exp(-(2*j-1)*(2*j-1)*(q1*pi*pi*tValue)/(4*L*L));
		value=value/(2*j-1);
		value=value*pow(-1,j+1);
		sum=sum+value;			
	}

	phiValue=((4*p0)/pi)*sum;
	phiValue=(c1/q1+b1*c2/(c2-q1))*phiValue;

	return phiValue;
}

double doubleDataProcessing::storagePhi2AnalyticalSolution(double yValue, double tValue, double L,
	double p0, double b1, double b2, double c1, double c2, double q2)
{
	double phiValue;
	double pi=4*atan(1);
	double value, sum=0;

	for(int j=1; j<301; j++)
	{
		value=cos((2*j-1)*(pi*yValue)/(2*L));
		value=value*exp(-(2*j-1)*(2*j-1)*(q2*pi*pi*tValue)/(4*L*L));
		value=value/(2*j-1);
		value=value*pow(-1,j+1);
		sum=sum+value;			
	}

	phiValue=((4*p0)/pi)*sum;
	phiValue=(b2*c1/(c1-q2)+c2/q2)*phiValue;

	return phiValue;
}

void doubleDataProcessing::exportStorageAnalyticalSolution(double Ly, double alpha1,
	double alpha2, double M, double S11, double S12, double S22, double KPore, double KFrac,
	double mu_f, double sigmab, double dt, int timeStep, string pairName)
{
	int pointsExact=5000;
	double yValue;
	double dy=Ly/pointsExact;
	vector<double> yExact;
	vector<double> vExact;
	vector<double> pPoreExact;
	vector<double> pFracExact;
	string fieldName;
	double tValue=dt*timeStep;

	double alpha=alpha1+alpha2;
	double p0=-(alpha/(M*(S11+2*S12+S22)+alpha*alpha))*sigmab;
	double c1=(KPore/mu_f)/(S11+alpha1*alpha1/M);
	double c2=(KFrac/mu_f)/(S22+alpha2*alpha2/M);
	double b1=(S12+alpha1*alpha2/M)/(S11+alpha1*alpha1/M);
	double b2=(S12+alpha1*alpha2/M)/(S22+alpha2*alpha2/M);
	double q1=(c1+c2+sqrt((c1-c2)*(c1-c2)+4*b1*b2*c1*c2))/(2*(1-b1*b2));
	double q2=(c1+c2-sqrt((c1-c2)*(c1-c2)+4*b1*b2*c1*c2))/(2*(1-b1*b2));
	double phi1, phi2;

	yExact.resize(pointsExact+1);
	vExact.resize(pointsExact+1);
	pPoreExact.resize(pointsExact+1);
	pFracExact.resize(pointsExact+1);

	yExact[0]=Ly;
	for(int i=1; i<pointsExact+1; i++) yExact[i]=yExact[i-1]-dy;
	
	for(int i=0; i<pointsExact+1; i++)
	{
		phi1=storagePhi1AnalyticalSolution(yExact[i],tValue,Ly,p0,b1,b2,c1,c2,q1);
		phi2=storagePhi2AnalyticalSolution(yExact[i],tValue,Ly,p0,b1,b2,c1,c2,q2);

		pPoreExact[i]=(c1-q2)*q1/(b1*b2*c1*q1*q2-c1*(c1-q2)*(c2-q1))*(-(c2-q1)*phi1+b1*q2*phi2);
		pFracExact[i]=(c2-q1)*q2/(b1*b2*c2*q1*q2-c2*(c1-q2)*(c2-q1))*(b2*q1*phi1-(c1-q2)*phi2);
	}

	fieldName="storageDouble_"+pairName+"_YExact_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(yExact,fieldName);

	fieldName="storageDouble_"+pairName+"_VExact_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(vExact,fieldName);

	fieldName="storageDouble_"+pairName+"_PPoreExact_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(pPoreExact,fieldName);

	fieldName="storageDouble_"+pairName+"_PFracExact_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(pFracExact,fieldName);

	ofstream myFile("../export/solveStorageDoubleRunInfo.txt",fstream::app);
	if(myFile.is_open())
	{
		myFile << timeStep << "\n";

		myFile.close();
	}

	return;
}

double doubleDataProcessing::leakingAnalyticalSolution(double yValue, double tValue, double L,
	double p0, double b1, double b2, double c1, double c2, double q1, double q2)
{
	double solution;
	double pi=4*atan(1);
	double value1, sum1=0;
	double value2, sum2=0;
	double A1, A2;

	A1=-((c1-q2)*(c2-q1)*c2*q1+b2*c1*(c2-q1)*q1*q2)/(b1*b2*c1*c2*q1*q2-c1*c2*(c1-q2)*(c2-q1));
	A2=(b1*(c1-q2)*c2*q1*q2+c1*(c1-q2)*(c2-q1)*q2)/(b1*b2*c1*c2*q1*q2-c1*c2*(c1-q2)*(c2-q1));

	for(int j=1; j<301; j++)
	{
		value1=cos((2*j-1)*(pi*yValue)/(2*L));
		value1=value1*(exp(-(2*j-1)*(2*j-1)*(q1*pi*pi*tValue)/(4*L*L))-1);
		value1=value1/((2*j-1)*(2*j-1)*(2*j-1));
		value1=value1*pow(-1,j);
		sum1=sum1+value1;

		value2=cos((2*j-1)*(pi*yValue)/(2*L));
		value2=value2*(exp(-(2*j-1)*(2*j-1)*(q2*pi*pi*tValue)/(4*L*L))-1);
		value2=value2/((2*j-1)*(2*j-1)*(2*j-1));
		value2=value2*pow(-1,j);
		sum2=sum2+value2;
	}


	solution=sum1*L*L*p0/q1*16*A1/(pi*pi*pi)*(c1/q1+b1*c2/(c2-q1));
	solution+=sum2*L*L*p0/q1*16*A2/(pi*pi*pi)*(c2/q2+b2*c1/(c1-q2));

	return solution;
}

void doubleDataProcessing::exportLeakingAnalyticalSolution(double Ly, double alpha1,
	double alpha2, double M, double S11, double S12, double S22, double KPore, double KFrac,
	double mu_f, double sigmab, double leak, double dt, int timeStep, string pairName)
{
	int pointsExact=5000;
	double yValue;
	double dy=Ly/pointsExact;
	vector<double> yExact;
	vector<double> vExact;
	vector<double> pPoreExact;
	vector<double> pFracExact;
	string fieldName;
	double tValue=dt*timeStep;

	yExact.resize(pointsExact+1);
	vExact.resize(pointsExact+1);
	pPoreExact.resize(pointsExact+1);
	pFracExact.resize(pointsExact+1);

	double alpha=alpha1+alpha2;
	double p0=-(alpha/(M*(S11+2*S12+S22)+alpha*alpha))*sigmab;
	double c1=(KPore/mu_f)/(S11+alpha1*alpha1/M);
	double c2=(KFrac/mu_f)/(S22+alpha2*alpha2/M);
	double b1=(S12+alpha1*alpha2/M)/(S11+alpha1*alpha1/M);
	double b2=(S12+alpha1*alpha2/M)/(S22+alpha2*alpha2/M);
	b1=0;
	b2=0;
	double q1=(c1+c2+sqrt((c1-c2)*(c1-c2)+4*b1*b2*c1*c2))/(2*(1-b1*b2));
	double q2=(c1+c2-sqrt((c1-c2)*(c1-c2)+4*b1*b2*c1*c2))/(2*(1-b1*b2));
	q1=c1;
	q2=c2;
	double leak1=leak/(S11+alpha1*alpha1/M);
	double leak2=leak/(S22+alpha2*alpha2/M);

	yExact[0]=Ly;
	for(int i=1; i<pointsExact+1; i++) yExact[i]=yExact[i-1]-dy;
	
	for(int i=0; i<pointsExact+1; i++)
	{
		pPoreExact[i]=leak2/b2*leakingAnalyticalSolution(yExact[i],tValue,Ly,p0,b1,b2,c1,c2,q1,q2);
		pFracExact[i]=leak1/b1*leakingAnalyticalSolution(yExact[i],tValue,Ly,p0,b1,b2,c1,c2,q1,q2);

	}

	fieldName="leakingDouble_"+pairName+"_YExact_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(yExact,fieldName);

	fieldName="leakingDouble_"+pairName+"_VExact_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(vExact,fieldName);

	fieldName="leakingDouble_"+pairName+"_PPoreExact_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(pPoreExact,fieldName);

	fieldName="leakingDouble_"+pairName+"_PFracExact_dt="+to_string(dt)+"_timeStep="+
		to_string(timeStep);
	export1DFieldToTxt(pFracExact,fieldName);

	ofstream myFile("../export/solveLeakingDoubleRunInfo.txt",fstream::app);
	if(myFile.is_open())
	{
		myFile << timeStep << "\n";

		myFile.close();
	}

	return;
}