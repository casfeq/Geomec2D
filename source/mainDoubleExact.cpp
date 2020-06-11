/*
	This source code implements a Finite Volume Method for discretization and solution of the consolidation problems with double porosity. The governing equations are discretized within the FVM and the resulting linear system of equations is solved with a LU Factorization found in PETSc [1].

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2020.
	
 	[1] BALAY et al. PETSc User Manual. Technical Report, Argonne National Laboratory, 2017.
*/

#include "customPrinter.hpp"
#include "exportRunInfo.hpp"
#include "benchmarking.hpp"

int main(int argc, char** args)
{	
	string myGridType=args[1];
	string myInterpScheme=args[2];
	string myMedium=args[3];
	vector<int> problemsSolved;
	problemsSolved.push_back(atoi(args[4])/100%10);
	problemsSolved.push_back(atoi(args[4])/10%10);
	problemsSolved.push_back(atoi(args[4])%10);
	double pore=stod(args[5]);
	double frac=stod(args[6]);

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
	double porosity;
	inFile >> porosity;
	inFile >> myProperties.permeability;
	inFile >> myProperties.fluidViscosity;
	inFile >> myProperties.fluidDensity;
	inFile.close();
	myProperties.macroPorosity=frac*porosity;
	myProperties.porosity=pore*porosity;
	myProperties.macroPermeability=myProperties.permeability*1e3;
	
/*		GRID DEFINITION
	----------------------------------------------------------------*/

	// Consolidation coefficient
	double storativity,fluidViscosity,permeability,fluidCompressibility,solidCompressibility,
		bulkCompressibility,longitudinalModulus,alpha;
	porosity=myProperties.porosity+myProperties.macroPorosity;
	fluidViscosity=myProperties.fluidViscosity;
	permeability=(myProperties.macroPermeability+myProperties.permeability)/2;
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
	double dt=consolidationTime/5;
	double Lt=(Nt-1)*dt;
	
/*		OTHER PARAMETERS
	----------------------------------------------------------------*/

	double columnLoad=-10e3; // Pa
	
/*		PETSC INITIALIZE
	----------------------------------------------------------------*/

	PetscErrorCode ierr;
	ierr=PetscInitialize(&argc,&args,(char*)0,NULL);CHKERRQ(ierr);

/*		SOLVE BENCHMARKING PROBLEMS
	----------------------------------------------------------------*/

	cout << "Grid type: " << myGridType << "\n";
	cout << "Interpolation scheme: " << myInterpScheme << "\n";
	cout << "Minimum time-step: " << consolidationTime/6 << "\n";
	cout << "Medium:" << myProperties.pairName << "\n";
	for(int i=0; i<3; i++)
	{
		if(problemsSolved[i]==2)
		{
			cout << "Solved Sealed Column (Double Porosity) for: \n";
			createSolveRunInfo(myGridType,myInterpScheme,"SealedDouble");
			exportSolveRunInfo(dt,"SealedDouble_"+myMedium);
			ierr=sealedDouble(myGridType,myInterpScheme,Nt,mesh,Lt,0,columnLoad,myProperties);
				CHKERRQ(ierr);
		}
		else if(problemsSolved[i]==4)
		{
			cout << "Solved Drained Column (Double Porosity / Storage) for: \n";
			createSolveRunInfo(myGridType,myInterpScheme,"StorageDouble");
			exportSolveRunInfo(dt,"StorageDouble_"+myMedium);
			ierr=storageDouble(myGridType,myInterpScheme,Nt,mesh,Lt,0,columnLoad,myProperties);
				CHKERRQ(ierr);
		}
		else if(problemsSolved[i]==8)
		{
			cout << "Solved Drained Column (Double Porosity / Leaking) for: \n";
			createSolveRunInfo(myGridType,myInterpScheme,"LeakingDouble");
			exportSolveRunInfo(dt,"LeakingDouble_"+myMedium);
			ierr=leakingDouble(myGridType,myInterpScheme,Nt,mesh,Lt,0,columnLoad,myProperties);
				CHKERRQ(ierr);
		}
	}
	
/*		PETSC FINALIZE
	----------------------------------------------------------------*/
	
	ierr=PetscFinalize();CHKERRQ(ierr);

	return ierr;
};