/*
	This header is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The class defined 
	here contains the functions for the solution of the linear system which represents the 
	discretized problem of poroelasticity. The linear system of equations is solved with LU 
	Factorization found in PETSc [1].
	
 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.

 	[1] BALAY et al. PETSc User Manual. Technical Report, Argonne National Laboratory, 2017.
*/

#include <iostream>
#include <petscksp.h>
#include <string>
#include <vector>

using namespace std;

class linearSystemSolver
{
public:
	// Class variables
	vector<vector<double>> coefficientsMatrix;
	vector<double> sparseCoefficientsRow;
	vector<double> sparseCoefficientsColumn;
	vector<double> sparseCoefficientsValue;
	vector<double> independentTermsArray;
	vector<vector<double>> uField;
	vector<vector<double>> vField;
	vector<vector<double>> pField;
	vector<vector<double>> pMField;
	int Nu;
	int Nv;
	int NP;
	int Nt;
	vector<vector<int>> uDisplacementFVIndex;
	vector<vector<int>> vDisplacementFVIndex;
	vector<vector<int>> pressureFVIndex;
	vector<vector<int>> uDisplacementFVCoordinates;
	vector<vector<int>> vDisplacementFVCoordinates;
	vector<vector<int>> pressureFVCoordinates;
	PetscErrorCode ierr;
	Mat coefficientsMatrixPETSc;
	Vec independentTermsArrayPETSc;
	Vec linearSystemSolutionPETSc;
	IS perm, iperm;
	MatFactorInfo info;

	// Class functions
	int getUDisplacementFVPosition(int,int);
	int getVDisplacementFVPosition(int,int);
	int getPressureFVPosition(int,int);
	int coefficientsMatrixLUFactorization();
	int createPETScArrays();
	int zeroPETScArrays();
	int setRHSValue(vector<double>);
	int solveLinearSystem();
	int setFieldValue(int);
	double mandelErrorCalculation(string,double,double,int,double,double,double,double);
	double mandelStaggeredErrorCalculation(double,double,int,double,double,double,double);
	double mandelCollocatedErrorCalculation(double,double,int,double,double,double,double);
	int createMacroPressureField(vector<vector<double>>);
	int setMacroFieldValue(int);

	// Constructor
	linearSystemSolver(vector<vector<double>>,vector<double>,vector<double>,vector<double>,
		vector<vector<double>>,vector<vector<double>>,vector<vector<double>>,int,int,int,int,
		vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,vector<vector<int>>,
		vector<vector<int>>,vector<vector<int>>);

	// Destructor
	~linearSystemSolver();
};

linearSystemSolver::linearSystemSolver(vector<vector<double>> myCoefficientsMatrix, 
	vector<double> mySparseCoefficientsRow, vector<double> mySparseCoefficientsColumn, 
	vector<double> mySparseCoefficientsValue, vector<vector<double>> uDisplacementField,
	vector<vector<double>> vDisplacementField, vector<vector<double>> pressureField, int myNu,
	int myNv, int myNP, int myNt, vector<vector<int>> idU, vector<vector<int>> idV,
	vector<vector<int>> idP, vector<vector<int>> cooU, vector<vector<int>> cooV,
	vector<vector<int>> cooP)
{
	coefficientsMatrix=myCoefficientsMatrix;
	sparseCoefficientsRow=mySparseCoefficientsRow;
	sparseCoefficientsColumn=mySparseCoefficientsColumn;
	sparseCoefficientsValue=mySparseCoefficientsValue;
	uField=uDisplacementField;
	vField=vDisplacementField;
	pField=pressureField;
	Nu=myNu;
	Nv=myNv;
	NP=myNP;
	Nt=myNt;
	uDisplacementFVIndex=idU;
	vDisplacementFVIndex=idV;
	pressureFVIndex=idP;
	uDisplacementFVCoordinates=cooU;
	vDisplacementFVCoordinates=cooV;
	pressureFVCoordinates=cooP;

	return;
}

linearSystemSolver::~linearSystemSolver()
{
	MatDestroy(&coefficientsMatrixPETSc);
	VecDestroy(&independentTermsArrayPETSc);
	VecDestroy(&linearSystemSolutionPETSc);
}

int linearSystemSolver::getUDisplacementFVPosition(int x, int y)
{
	int uDisplacementFVPosition;

	uDisplacementFVPosition=uDisplacementFVIndex[x][y]-1;

	return uDisplacementFVPosition;
}

int linearSystemSolver::getVDisplacementFVPosition(int x, int y)
{
	int vDisplacementFVPosition;

	vDisplacementFVPosition=vDisplacementFVIndex[x][y]+Nu-1;

	return vDisplacementFVPosition;
}

int linearSystemSolver::getPressureFVPosition(int x, int y)
{
	int pressureFVPosition;

	pressureFVPosition=pressureFVIndex[x][y]+Nu+Nv-1;

	return pressureFVPosition;
}

int linearSystemSolver::coefficientsMatrixLUFactorization()
{
	PetscInt n=coefficientsMatrix.size();
	PetscInt nonZeroEntries=sparseCoefficientsValue.size();
	PetscInt rowNo, colNo;
	PetscScalar value;

	ierr=MatCreate(PETSC_COMM_WORLD,&coefficientsMatrixPETSc);CHKERRQ(ierr);
	ierr=MatSetSizes(coefficientsMatrixPETSc,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
	ierr=MatSetFromOptions(coefficientsMatrixPETSc);CHKERRQ(ierr);
	ierr=MatSetUp(coefficientsMatrixPETSc);CHKERRQ(ierr);
	ierr=MatZeroEntries(coefficientsMatrixPETSc);CHKERRQ(ierr);

	for(int i=0; i<nonZeroEntries; i++)
	{	
		rowNo=sparseCoefficientsRow[i];
		colNo=sparseCoefficientsColumn[i];
		value=sparseCoefficientsValue[i];
		ierr=MatSetValue(coefficientsMatrixPETSc,rowNo,colNo,value,ADD_VALUES);CHKERRQ(ierr);
	}

	ierr=MatAssemblyBegin(coefficientsMatrixPETSc,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr=MatAssemblyEnd(coefficientsMatrixPETSc,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr=MatGetOrdering(coefficientsMatrixPETSc,MATORDERINGRCM,&perm,&iperm);CHKERRQ(ierr);

	ierr=MatFactorInfoInitialize(&info);CHKERRQ(ierr);
	info.fill=1.0;
	info.dt=0;
	info.dtcol=0;
	info.zeropivot=0;
	info.pivotinblocks=0;

	ierr=MatLUFactor(coefficientsMatrixPETSc,perm,iperm,&info);CHKERRQ(ierr);

	return ierr;
}

int linearSystemSolver::createPETScArrays()
{
	PetscInt n=coefficientsMatrix.size();

	ierr=VecCreateSeq(PETSC_COMM_SELF,n,&independentTermsArrayPETSc);CHKERRQ(ierr);
	ierr=VecSet(independentTermsArrayPETSc,0.0);CHKERRQ(ierr);
	ierr=VecAssemblyBegin(independentTermsArrayPETSc);CHKERRQ(ierr);
	ierr=VecAssemblyEnd(independentTermsArrayPETSc);CHKERRQ(ierr);
	ierr=VecDuplicate(independentTermsArrayPETSc,&linearSystemSolutionPETSc);CHKERRQ(ierr);

	return ierr;
}

int linearSystemSolver::zeroPETScArrays()
{
	ierr=VecZeroEntries(independentTermsArrayPETSc);CHKERRQ(ierr);
	ierr=VecAssemblyBegin(independentTermsArrayPETSc);CHKERRQ(ierr);
	ierr=VecAssemblyEnd(independentTermsArrayPETSc);CHKERRQ(ierr);
	ierr=VecZeroEntries(linearSystemSolutionPETSc);CHKERRQ(ierr);
	ierr=VecAssemblyBegin(linearSystemSolutionPETSc);CHKERRQ(ierr);
	ierr=VecAssemblyEnd(linearSystemSolutionPETSc);CHKERRQ(ierr);

	return ierr;
}

int linearSystemSolver::setRHSValue(vector<double> independentTermsArray)
{
	PetscInt n=independentTermsArray.size();
	PetscScalar value;

	for(int i=0; i<n; i++)
	{
		value=independentTermsArray[i];
		ierr=VecSetValue(independentTermsArrayPETSc,i,value,ADD_VALUES);CHKERRQ(ierr);
	}

	ierr=VecAssemblyBegin(independentTermsArrayPETSc);CHKERRQ(ierr);
	ierr=VecAssemblyEnd(independentTermsArrayPETSc);CHKERRQ(ierr);

	return ierr;
}

int linearSystemSolver::solveLinearSystem()
{	
	ierr=MatSolve(coefficientsMatrixPETSc,independentTermsArrayPETSc,linearSystemSolutionPETSc);
		CHKERRQ(ierr);
	ierr=VecAssemblyBegin(linearSystemSolutionPETSc);CHKERRQ(ierr);
	ierr=VecAssemblyEnd(linearSystemSolutionPETSc);CHKERRQ(ierr);

	return ierr;
}

int linearSystemSolver::setFieldValue(int timeStep)
{
	PetscScalar value;

	for(int i=0; i<Nu; i++)
	{
		ierr=VecGetValues(linearSystemSolutionPETSc,1,&i,&value);CHKERRQ(ierr);
		uField[i][timeStep]=value;
	}

	for(int i=Nu; i<Nu+Nv; i++)
	{
		ierr=VecGetValues(linearSystemSolutionPETSc,1,&i,&value);CHKERRQ(ierr);
		vField[i-Nu][timeStep]=value;
	}

	for(int i=Nu+Nv; i<Nu+Nv+NP; i++)
	{
		ierr=VecGetValues(linearSystemSolutionPETSc,1,&i,&value);CHKERRQ(ierr);
		pField[i-Nu-Nv][timeStep]=value;
	}

	return ierr;
}

double linearSystemSolver::mandelErrorCalculation(string gridType, double dx, double dy,
	int timeStep, double M, double lambda, double alpha, double F)
{
	double error;

	if(gridType=="staggered") error=mandelStaggeredErrorCalculation(dx,dy,timeStep,M,lambda,alpha,
		F);
	else if(gridType=="collocated") error=mandelCollocatedErrorCalculation(dx,dy,timeStep,M,lambda,
		alpha,F);

	return error;
}

double linearSystemSolver::mandelStaggeredErrorCalculation(double dx, double dy, int timeStep,
	double M, double lambda, double alpha, double F)
{
	int i, j;
	int u_P, u_E, v_P, v_S, P_P;
	double currStrain, currPressure;
	int Nx=vDisplacementFVIndex[0].size();
	double stressTrial=0;

	for(int FVCounter=0; FVCounter<Nx; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		u_P=getUDisplacementFVPosition(i,j);
		u_E=getUDisplacementFVPosition(i,j+1);
		v_P=getVDisplacementFVPosition(i,j)-Nu;
		v_S=getVDisplacementFVPosition(i+1,j)-Nu;
		P_P=getPressureFVPosition(i,j)-(Nu+Nv);

		currStrain=(uField[u_E][timeStep]-uField[u_P][timeStep])/dx+(vField[v_P][timeStep]-
			vField[v_S][timeStep])/dy;
		currPressure=pField[P_P][timeStep];

		stressTrial=stressTrial+(M+lambda)*currStrain-2*alpha*currPressure;
	}

	double error=stressTrial*dx-F;

	return error;
}

double linearSystemSolver::mandelCollocatedErrorCalculation(double dx, double dy, int timeStep,
	double M, double lambda, double alpha, double F)
{
	int i, j;
	int u_P, u_E, u_W, v_P, v_S, P_P;
	double currStrain, currPressure;
	int Nx=vDisplacementFVIndex[0].size();
	double stressTrial=0;

	for(int FVCounter=0; FVCounter<Nx; FVCounter++)
	{
		i=vDisplacementFVCoordinates[FVCounter][0]-1;
		j=vDisplacementFVCoordinates[FVCounter][1]-1;

		if(i==0) // Western border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_E=getUDisplacementFVPosition(i,j+1);
			v_P=getVDisplacementFVPosition(i,j)-Nu;
			v_S=getVDisplacementFVPosition(i+1,j)-Nu;
			P_P=getPressureFVPosition(i,j)-(Nu+Nv);

			currStrain=(uField[u_E][timeStep]-uField[u_P][timeStep])/dx+(vField[v_P][timeStep]-
				vField[v_S][timeStep])/dy;
			currPressure=pField[P_P][timeStep];

			stressTrial=stressTrial+(M+lambda)*currStrain/2-alpha*currPressure;
		}
		else if(i==Nx-1) // Eastern border
		{
			u_P=getUDisplacementFVPosition(i,j);
			u_W=getUDisplacementFVPosition(i,j-1);
			v_P=getVDisplacementFVPosition(i,j)-Nu;
			v_S=getVDisplacementFVPosition(i+1,j)-Nu;
			P_P=getPressureFVPosition(i,j)-(Nu+Nv);

			currStrain=(uField[u_P][timeStep]-uField[u_W][timeStep])/dx+(vField[v_P][timeStep]-
				vField[v_S][timeStep])/dy;
			currPressure=pField[P_P][timeStep];

			stressTrial=stressTrial+(M+lambda)*currStrain/2-alpha*currPressure;
		}
		else
		{
			u_W=getUDisplacementFVPosition(i,j-1);
			u_E=getUDisplacementFVPosition(i,j+1);
			v_P=getVDisplacementFVPosition(i,j)-Nu;
			v_S=getVDisplacementFVPosition(i+1,j)-Nu;
			P_P=getPressureFVPosition(i,j)-(Nu+Nv);

			currStrain=(uField[u_E][timeStep]-uField[u_W][timeStep])/dx+(vField[v_P][timeStep]-
				vField[v_S][timeStep])/dy;
			currPressure=pField[P_P][timeStep];

			stressTrial=stressTrial+(M+lambda)*currStrain-2*alpha*currPressure;	
		}
	}

	double error=stressTrial*dx-F;

	return error;
}

int linearSystemSolver::createMacroPressureField(vector<vector<double>> pressureField)
{
	pMField=pressureField;

	return ierr;
}

int linearSystemSolver::setMacroFieldValue(int timeStep)
{
	PetscScalar value;

	for(int i=Nu+Nv+NP; i<Nu+Nv+NP+NP; i++)
	{
		ierr=VecGetValues(linearSystemSolutionPETSc,1,&i,&value);CHKERRQ(ierr);
		pMField[i-Nu-Nv-NP][timeStep]=value;
	}

	return ierr;
}