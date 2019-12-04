/*
	This header is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Poroelasticity Problems". The functions defined
	here are responsible for creating reports which allows the classification of the results.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.
*/

#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

void createConvergenceRunInfo(string gridType, string interpScheme, string problemName)
{
	string fileName="../export/convergence"+problemName+"RunInfo.txt";

	ofstream myFile(fileName);
	if(myFile.is_open())
	{
		myFile << gridType << "\n" << interpScheme << "\n";

		myFile.close();
	}

	return;
}

void exportConvergenceRunInfo(double dt, string problemName)
{
	string fileName="../export/convergence"+problemName+"RunInfo.txt";

	ofstream myFile(fileName,fstream::app);
	if(myFile.is_open())
	{
		myFile << fixed << setprecision(6) << dt << "\n";

		myFile.close();
	}

	return;
}

void createSolveRunInfo(string gridType, string interpScheme, string problemName)
{
	string fileName="../export/solve"+problemName+"RunInfo.txt";

	ofstream myFile(fileName);
	if(myFile.is_open())
	{
		myFile << gridType << "\n" << interpScheme << "\n";

		myFile.close();
	}

	return;
}

void exportSolveRunInfo(double dt, string problemName)
{
	string fileName="../export/solve"+problemName+"RunInfo.txt";

	ofstream myFile(fileName,fstream::app);
	if(myFile.is_open())
	{
		myFile << fixed << setprecision(6) << dt << "\n";

		myFile.close();
	}

	return;
}