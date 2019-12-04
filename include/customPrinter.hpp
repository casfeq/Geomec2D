/*
	This header is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Poroelasticity Problems". The template functions
	defined here are used for printing on the terminal scalars, vectors and matrices.

 	Written by FERREIRA, C. A. S.

 	Florianópolis, 2019.
*/

#include <iostream>
#include <string>
#include <vector>

#define printscalar(scalar) PRINTSCALAR(#scalar,(scalar))
#define printvector(vector) PRINTVECTOR(#vector,(vector))
#define printmatrix(matrix) PRINTMATRIX(#matrix,(matrix))
#define printsparse(sparse) PRINTSPARSE(#sparse,(sparse))

using namespace std;

template<typename scalarType>
void PRINTSCALAR(string name, scalarType value)
{
	cout << name << "=" << value;

    return;
}

template<typename vecType>
void PRINTVECTOR(string name, vecType value)
{
	int rowNo=value.size();

	cout << name << "=\n";

	for(int i=0; i<rowNo; i++)
	{
		cout << "Row "<< i+1 << ":\t" << value[i] << "\n";
		// cout << value[i] << "\n";
	}

	return;
}

template<typename matType>
void PRINTMATRIX(string name, matType value)
{
	int rowNo=value.size();
	int colNo=value[0].size();

	cout << name << "=\n";

	for(int i=0; i<rowNo; i++)
	{
		for(int j=0; j<colNo; j++)
		{
			cout << value[i][j] << "\t";
		}

		cout << "\n";
	}
}

template<typename sparseType>
void PRINTSPARSE(string name, sparseType value)
{
	int rowNo=value.size();
	int colNo=value[0].size();

	cout << name << "=\n";

	for(int i=0; i<rowNo; i++)
	{
		cout << "Row " << i+1 << ":\t";

		for(int j=0; j<colNo; j++)
		{
			if(value[i][j]!=0)
			{
				cout << "(" << j+1 << "," << value[i][j] << "), ";
			}
		}

		cout << "\n";
	}
}

void newline()
{
	cout << "\n";
}

void tab()
{
	cout << "\t";
}