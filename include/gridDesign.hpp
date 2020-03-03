/*
	This header is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The class defined 
	here contains the functions for creation of the grid for discretizing the poroelasticity 
	problem.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.
*/

#include <iostream>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

class gridDesign
{
public:
	// Class variables
	int sPointsNo; // No of points of the surface s
	int Nx; // No of FV in x
	int Ny; // No of FV in y
	string gridType;
	double Lx; // Grid size in x [m]
	double Ly; // Grid size in y [m]
	vector<vector<double>> sCoordinates; // Coordinates of the surface s
	double dx, dy; // FV sizes in x and y [m]
	double dt; // Time-step [s]
	double h; // Characteristic length [m]
	double Lt; // Total simulation time [s]
	int Nt; // No of time-steps + 1 (initial condition)
	vector<vector<double>> centroidsXCoordinates;
	vector<vector<double>> centroidsYCoordinates;
	vector<vector<int>> generalFVIndex;
	int numberOfActiveGeneralFV;
	vector<vector<int>> uDisplacementFVIndex;
	int numberOfActiveUDisplacementFV;
	vector<vector<int>> vDisplacementFVIndex;
	int numberOfActiveVDisplacementFV;
	vector<vector<int>> generalFVCoordinates;
	vector<vector<int>> uDisplacementFVCoordinates;
	vector<vector<int>> vDisplacementFVCoordinates;
	vector<vector<int>> horizontalFacesStatus;
	vector<vector<int>> verticalFacesStatus;
	vector<vector<double>> uDisplacementField;
	vector<vector<double>> vDisplacementField;
	vector<vector<double>> pressureField;

	// Class functions
	void buildGrid(string);
	void computeFVMesh();
	void resizeVectors(string);
	void computeCentroidCoordinates();
	bool inPolygon(double,double);
	void buildGeneralFVIndexes();
	void buildFaces();
	void buildUDisplacementFVIndexes(string);
	void buildVDisplacementFVIndexes(string);
	void buildFieldVectors();

	// Constructor
	gridDesign(int,int,int,double,double,double,string,vector<vector<double>>);

	// Destructor
	~gridDesign();
};

gridDesign::gridDesign(int numberOfXFV, int numberOfYFV, int numberOfTimeSteps, double gridSizeX, 
	double gridSizeY, double totalSimulationTime, string myGridType,
	vector<vector<double>> coordinatesOfSurfacePoints)
{
	Nx=numberOfXFV;
	Ny=numberOfYFV;
	Nt=numberOfTimeSteps;
	Lx=gridSizeX;
	Ly=gridSizeY;
	Lt=totalSimulationTime;
	gridType=myGridType;
	sCoordinates=coordinatesOfSurfacePoints;
	buildGrid(gridType);
}

gridDesign::~gridDesign(){}

void gridDesign::buildGrid(string gridType)
{
	computeFVMesh();
	resizeVectors(gridType);
	computeCentroidCoordinates();
	buildGeneralFVIndexes();
	buildFaces();
	buildUDisplacementFVIndexes(gridType);
	buildVDisplacementFVIndexes(gridType);
	buildFieldVectors();

	return;
}

void gridDesign::computeFVMesh()
{
	dx=Lx/Nx;
	dy=Ly/Ny;
	dt=Lt/(Nt-1);
	h=sqrt(dx*dy);
	sPointsNo=sCoordinates.size();

	return;
}

void gridDesign::resizeVectors(string gridType)
{
	int rowNo; // No of rows
	int colNo; // No of columns
	
	// Resize verticalFacesStatus
	rowNo=Ny;
	colNo=Nx+1;
	verticalFacesStatus.resize(rowNo);
	for(int i=0; i<rowNo; ++i)
	{
		verticalFacesStatus[i].resize(colNo);
	}

	// Resize uDisplacementFVIndex
	if(gridType=="staggered") rowNo=Ny;
	else if(gridType=="collocated") rowNo=Ny+1;
	if(gridType=="staggered") colNo=Nx+1;
	else if(gridType=="collocated") colNo=Nx+1;
	uDisplacementFVIndex.resize(rowNo);
	for(int i=0; i<rowNo; ++i)
	{
		uDisplacementFVIndex[i].resize(colNo);
	}

	// Resize horizontalFacesStatus
	rowNo=Ny+1;
	colNo=Nx;
	horizontalFacesStatus.resize(rowNo);
	for(int i=0; i<rowNo; ++i)
		horizontalFacesStatus[i].resize(colNo);

	// Resize vDisplacementFVIndex
	if(gridType=="staggered") rowNo=Ny+1;
	else if(gridType=="collocated") rowNo=Ny+1;
	if(gridType=="staggered") colNo=Nx;
	else if(gridType=="collocated") colNo=Nx+1;
	vDisplacementFVIndex.resize(rowNo);
	for(int i=0; i<rowNo; ++i)
		vDisplacementFVIndex[i].resize(colNo);

	// Resize centroidsXCoordinates, centroidsYCoordinates and generalFVIndex
	if(gridType=="staggered") rowNo=Ny;
	else if(gridType=="collocated") rowNo=Ny+1;
	if(gridType=="staggered") colNo=Nx;
	else if(gridType=="collocated") colNo=Nx+1;
	centroidsXCoordinates.resize(rowNo);
	centroidsYCoordinates.resize(rowNo);
	generalFVIndex.resize(rowNo);
	for(int i=0; i<rowNo; ++i)
	{
		centroidsXCoordinates[i].resize(colNo);
		centroidsYCoordinates[i].resize(colNo);
		generalFVIndex[i].resize(colNo);
	}

	// Resize uDisplacementFVCoordinates
	if(gridType=="staggered") rowNo=(Nx+1)*Ny;
	else if(gridType=="collocated") rowNo=(Nx+1)*(Ny+1);
	colNo=2;
	uDisplacementFVCoordinates.resize(rowNo);
	for(int i=0; i<rowNo; ++i)
		uDisplacementFVCoordinates[i].resize(colNo);

	// Resize vDisplacementFVCoordinates
	if(gridType=="staggered") rowNo=Nx*(Ny+1);
	else if(gridType=="collocated") rowNo=(Nx+1)*(Ny+1);
	colNo=2;
	vDisplacementFVCoordinates.resize(rowNo);
	for(int i=0; i<rowNo; ++i)
		vDisplacementFVCoordinates[i].resize(colNo);

	// Resize generalFVCoordinates
	if(gridType=="staggered") rowNo=Nx*Ny;
	else if(gridType=="collocated") rowNo=(Nx+1)*(Ny+1);
	colNo=2;
	generalFVCoordinates.resize(rowNo);
	for(int i=0; i<rowNo; ++i)
		generalFVCoordinates[i].resize(colNo);

	return;
}

void gridDesign::computeCentroidCoordinates()
{
	if(gridType=="staggered")
	{
		for(int i=0; i<Ny; i++)
		{
			for(int j=0; j<Nx; j++)
			{
				centroidsXCoordinates[i][j]=(j+0.5)*dx;
				centroidsYCoordinates[i][j]=Ly-(i+0.5)*dy;
			}
		}
	}
	else if(gridType=="collocated")
	{
		for(int i=0; i<Ny+1; i++)
		{
			for(int j=0; j<Nx+1; j++)
			{
				centroidsXCoordinates[i][j]=j*dx;
				centroidsYCoordinates[i][j]=Ly-i*dy;
			}
		}
	}

	return;
}

bool gridDesign::inPolygon(double x, double y)
{
	bool isInPolygon=0;
	int i, j;

	for(i=0, j=sPointsNo-1; i<sPointsNo; j=i++)
	{
		if(((sCoordinates[i][1]>=y)!=(sCoordinates[j][1]>y)) && (x<=(sCoordinates[j][0]
			-sCoordinates[i][0])*(y-sCoordinates[i][1])/(sCoordinates[j][1]-sCoordinates[i][1])
			+sCoordinates[i][0]))
    	isInPolygon=!isInPolygon;
	}

	return isInPolygon;
}

void gridDesign::buildGeneralFVIndexes()
{
	double x;
	double y;
	int ct=1;
	int i, j;
	int rowNo, colNo;
	bool isInPolygon;

	if(gridType=="staggered")
	{
		rowNo=Ny;
		colNo=Nx;
	}
	else if(gridType=="collocated")
	{
		rowNo=Ny+1;
		colNo=Nx+1;
	}

	for(i=0; i<rowNo; i++)
	{
		for(j=0; j<colNo; j++)
		{
			x=centroidsXCoordinates[i][j];
			y=centroidsYCoordinates[i][j];
			isInPolygon=inPolygon(x,y);

			if(isInPolygon) // If inside the polygon, receives a counter number
			{
				if(j==1 && generalFVIndex[i][0]==0)
				{
					generalFVIndex[i][0]=ct;
					generalFVCoordinates[ct-1][0]=i+1;
					generalFVCoordinates[ct-1][1]=1;
					ct=ct+1;
				}
				generalFVIndex[i][j]=ct;
				generalFVCoordinates[ct-1][0]=i+1;
				generalFVCoordinates[ct-1][1]=j+1;
				ct=ct+1;
			}
			else
			{
				generalFVIndex[i][j]=0; // If outside the polygon, receives a zero
			}
		}
	}

	numberOfActiveGeneralFV=ct-1;

	return;
}

void gridDesign::buildFaces()
{
	for(int i=0; i<Ny; i++)
	{
		for(int j=0; j<Nx; j++)
		{
			if(generalFVIndex[i][j]!=0) // If not zero, than its an active FV
			{
				verticalFacesStatus[i][j]=verticalFacesStatus[i][j]+1;
					// Right face status added by 1
				verticalFacesStatus[i][j+1]=verticalFacesStatus[i][j+1]+1;
					// Left face status added by 1
				horizontalFacesStatus[i][j]=horizontalFacesStatus[i][j]+1;
					// Upper face added by 1
				horizontalFacesStatus[i+1][j]=horizontalFacesStatus[i+1][j]+1;
					// Lower face status added by 1
			}
		}
	}

	return;
}

void gridDesign::buildUDisplacementFVIndexes(string gridType)
{
	int ct=0;
	int faceStatus;

	if(gridType=="collocated")
	{
		uDisplacementFVIndex=generalFVIndex;
		uDisplacementFVCoordinates=generalFVCoordinates;
		numberOfActiveUDisplacementFV=numberOfActiveGeneralFV;

		return;
	}

	for(int i=0; i<Ny; i++)
	{
		for(int j=0; j<Nx+1; j++)
		{
			faceStatus=verticalFacesStatus[i][j];
			if(faceStatus!=0)
			{
				ct++;
				uDisplacementFVIndex[i][j]=ct;
				uDisplacementFVCoordinates[ct-1][0]=i+1;
				uDisplacementFVCoordinates[ct-1][1]=j+1;
			}
			else
			{
				uDisplacementFVIndex[i][j]=0;
			}
		}
	}

	numberOfActiveUDisplacementFV=ct;

	return;
}

void gridDesign::buildVDisplacementFVIndexes(string gridType)
{
	int ct=0;
	int faceStatus;

	if(gridType=="collocated")
	{
		vDisplacementFVIndex=generalFVIndex;
		vDisplacementFVCoordinates=generalFVCoordinates;
		numberOfActiveVDisplacementFV=numberOfActiveGeneralFV;

		return;
	}

	for(int i=0; i<Ny+1; i++)
	{
		for(int j=0; j<Nx; j++)
		{
			faceStatus=horizontalFacesStatus[i][j];
			if(faceStatus!=0)
			{	
				ct++;
				vDisplacementFVIndex[i][j]=ct;
				vDisplacementFVCoordinates[ct-1][0]=i+1;
				vDisplacementFVCoordinates[ct-1][1]=j+1;
			}
			else
			{
				vDisplacementFVIndex[i][j]=0;
			}
		}
	}
	
	numberOfActiveVDisplacementFV=ct;

	return;
}

void gridDesign::buildFieldVectors()
{
	int rowNo, colNo;

	// Resize uDisplacementField
	rowNo=numberOfActiveUDisplacementFV;
	colNo=Nt;
	uDisplacementField.resize(rowNo);
	for(int i=0; i<rowNo; i++)
	{
		uDisplacementField[i].resize(colNo);
	}

	// Resize vDisplacementField
	rowNo=numberOfActiveVDisplacementFV;
	colNo=Nt;
	vDisplacementField.resize(rowNo);
	for(int i=0; i<rowNo; i++)
	{
		vDisplacementField[i].resize(colNo);
	}

	// Resize pressureField
	rowNo=numberOfActiveGeneralFV;
	colNo=Nt;
	pressureField.resize(rowNo);
	for(int i=0; i<rowNo; i++)
	{
		pressureField[i].resize(colNo);
	}

	return;
}