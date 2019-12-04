/*
	This source code implements a Finite Volume Method for discretization and solution of the 
	consolidation problem as part of a master's thesis entitled "Analysis Numerical Schemes in
	Collocated and Staggered Grids for Poroelasticity Problems". This source code uses the functions
	predefined for the solution of the problems presented and solved by Terzaghi [3] and Mandel [2].
	The governing equations are discretized within the FVM and the resulting linear system of
	equations is solved with a LU Factorization found in PETSc [1].

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.

 	[1] BALAY et al. PETSc User Manual. Technical Report, Argonne National Laboratory, 2017.
	[2] MANDEL, J. Consolidation Des Sols (Étude Mathématique). Géotechnique, v. 3, n. 7, pp. 287-
	299, 1953.
 	[3] TERZAGHI, K. Erdbaumechanik auf Bodenphysikalischer Grundlage. Franz Deuticke, Leipzig,
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
	porosity=myProperties.porosity;
	fluidViscosity=myProperties.fluidViscosity;
	permeability=myProperties.permeability;
	fluidCompressibility=1/myProperties.fluidBulkModulus;
	solidCompressibility=1/myProperties.solidBulkModulus;
	bulkCompressibility=1/myProperties.bulkModulus;
	longitudinalModulus=myProperties.bulkModulus+4*myProperties.shearModulus/3;
	alpha=1-solidCompressibility/bulkCompressibility;
	storativity=porosity*fluidCompressibility+(alpha-porosity)*solidCompressibility;
	double consolidationCoefficient=(permeability/fluidViscosity)/(storativity+
		alpha*alpha/longitudinalModulus);

	int Nt=501;
	int mesh=5;
	double h=1./mesh;
	double consolidationTime=h*h/consolidationCoefficient;
	double Lt=(Nt-1)*consolidationTime; // [s]
	double dt=Lt/(Nt-1);
/*		OTHER PARAMETERS
	----------------------------------------------------------------*/

	double g=-9.8; // m/s^2
	double columnLoad=-10e3; // Pa
	double mandelLoad=-10e4; // N/m
	double stripLoad=-10e3; // Pa
	
/*		PETSC INITIALIZE
	----------------------------------------------------------------*/

	PetscErrorCode ierr;
	ierr=PetscInitialize(&argc,&args,(char*)0,NULL);CHKERRQ(ierr);

/*		SOLVE BENCHMARKING PROBLEMS
	----------------------------------------------------------------*/

	cout << "Grid type: " << myGridType << "\n";
	cout << "Interpolation scheme: " << myInterpScheme << "\n";
	cout << "Minimum time-step: " << consolidationTime/6 << "\n";
	cout << "Solved sealed column for: \n";
	createSolveRunInfo(myGridType,myInterpScheme,"SealedColumn");
	exportSolveRunInfo(dt,"SealedColumn_"+myMedium);
	ierr=sealedColumn(myGridType,myInterpScheme,Nt,mesh,Lt,g,columnLoad,myProperties);
		CHKERRQ(ierr);
	cout << "Solved Terzaghi for: \n";
	createSolveRunInfo(myGridType,myInterpScheme,"Terzaghi");
	exportSolveRunInfo(dt,"Terzaghi_"+myMedium);
	ierr=terzaghi(myGridType,myInterpScheme,Nt,mesh,Lt,g,columnLoad,myProperties);CHKERRQ(ierr);
	cout << "Solved Mandel for: \n";
	createSolveRunInfo(myGridType,myInterpScheme,"Mandel");
	exportSolveRunInfo(dt,"Mandel_"+myMedium);
	ierr=mandel(myGridType,myInterpScheme,Nt,mesh,Lt,0,mandelLoad,myProperties);CHKERRQ(ierr);
	
/*		PETSC FINALIZE
	----------------------------------------------------------------*/
	
	ierr=PetscFinalize();CHKERRQ(ierr);

	return ierr;
};