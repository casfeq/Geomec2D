/*
	This source code implements a Finite Volume Method for discretization and solution of the 
	consolidation problem as part of a master's thesis entitled "Analysis Numerical Schemes in
	Collocated and Staggered Grids for Poroelasticity Problems". This source code uses the functions
	predefined for the solution of the problem presented and solved by Terzaghi [2]. The governing
	equations are discretized within the FVM and the resulting linear system of equations is solved
	with a LU Factorization found in PETSc [1]. The numerical error in comparison to the analytical
	solution is plotted against the meshsize to test the convergence of the method.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.

 	[1] BALAY et al. PETSc User Manual. Technical Report, Argonne National Laboratory, 2017.
 	[2] TERZAGHI, K. Erdbaumechanik auf Bodenphysikalischer Grundlage. Franz Deuticke, Leipzig,
 	1925.
*/

#include "customPrinter.hpp"
#include "exportRunInfo.hpp"
#include "benchmarking.hpp"

int main(int argc, char** args)
{	
	string myGridType=args[1];
	string myInterpScheme=args[2];
	string myMedium=args[3];

/*		PROPERTIES IMPORT
	----------------------------------------------------------------*/	

	poroelasticProperties myProperties;
	ifstream inFile;
	inFile.open("../input/"+myMedium+".txt");
	if(!inFile)
	{
		cout << "Unable to open properties file.";
		exit(1);
	}
	getline(inFile,myProperties.pairName);
	myProperties.pairName=myMedium;
	inFile >> myProperties.shearModulus;
	inFile >> myProperties.bulkModulus;
	inFile >> myProperties.solidBulkModulus;
	inFile >> myProperties.solidDensity;
	inFile >> myProperties.fluidBulkModulus;
	inFile >> myProperties.porosity;
	inFile >> myProperties.permeability;
	inFile >> myProperties.fluidViscosity;
	inFile >> myProperties.fluidDensity;
	inFile.close();
	
/*		GRID DEFINITION
	----------------------------------------------------------------*/
	// Consolidation coefficient
	double storativity,porosity,fluidViscosity,permeability,fluidCompressibility,
		solidCompressibility,bulkCompressibility,longitudinalModulus,alpha;
	porosity=myProperties.phi;
	fluidViscosity=myProperties.mu_f;
	permeability=myProperties.K;
	fluidCompressibility=myProperties.c_f;
	solidCompressibility=myProperties.c_s;
	bulkCompressibility=1/(myProperties.lambda+2*myProperties.G/3);
	longitudinalModulus=2*myProperties.G+myProperties.lambda;
	alpha=1-solidCompressibility/bulkCompressibility;
	storativity=porosity*fluidCompressibility+(alpha-porosity)*solidCompressibility;
	double c=(permeability/fluidViscosity)/(storativity+alpha*alpha/longitudinalModulus);

	double Lt=5e5; // [s]

	vector<int> timeSteps=
	{
		{101},
		{201},
		{401},
		{501}
	};

	vector<int> meshSize=
	{
		{3},
		{4},
		{5},
		{6},
		{7},
		{8},
		{9},
		{10},
		{11},
		{12},
		{13},
		{14},
		{15}	

	};
	
/*		OTHER PARAMETERS
	----------------------------------------------------------------*/

	double g=0; // m/s^2
	double sigmab=-10e3; // Pa
	
/*		PETSC INITIALIZE
	----------------------------------------------------------------*/

	PetscErrorCode ierr;

	ierr=PetscInitialize(&argc,&args,(char*)0,NULL);CHKERRQ(ierr);

/*		EXECUTE CONVERGENCE TEST
	----------------------------------------------------------------*/

	int timeStepsNo=timeSteps.size();
	int meshSizeNo=meshSize.size();
	double dt;
	double h;

	cout << "Grid type: " << myGridType << "\n";
	cout << "Interpolation scheme: " << myInterpScheme << "\n";
	cout << "Tested Terzaghi for: \n";

	createConvergenceRunInfo(myGridType,myInterpScheme,"Terzaghi");
	for(int i=0; i<timeStepsNo; i++)
	{
		dt=Lt/(timeSteps[i]-1);
		exportConvergenceRunInfo(dt,"Terzaghi");

		for(int j=0; j<meshSizeNo; j++)
		{
			ierr=convergence(myGridType,myInterpScheme,timeSteps[i],meshSize[j],Lt,g,sigmab,
				myProperties);CHKERRQ(ierr);
		}
	}

/*		PETSC FINALIZE
	----------------------------------------------------------------*/
	
	ierr=PetscFinalize();CHKERRQ(ierr);

	return ierr;
};