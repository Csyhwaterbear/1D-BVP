#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include "Process.h"
using namespace std;

void Process::Input(const string& file_path)
{
	this->file_path = file_path;
	ParseFile();
}

void Process::ParseFile()
{
	ifstream inputFile(file_path);
	if (!inputFile.is_open())
	{
		cerr << "Unable to open input file: " << file_path << endl;
		return;
	}

	string line;
	string currentSection = "";
	while (getline(inputFile, line))
	{
		if (line[0] == '*')
		{
			if (line[1] != '*')
			{
				currentSection = line.substr(1);
			}
			continue;
		}
		// cout << currentSection << endl;
		if (currentSection == "alpha")
		{
			stringstream ss(line);
			int number;
			double value;
			char comma;
			ss >> number >> comma >> value;
			alpha.push_back({number,value});
		}
		if (currentSection == "beta")
		{
			stringstream ss(line);
			int number;
			double value;
			char comma;
			ss >> number >> comma >> value;
			beta.push_back({number,value});
		}
		if (currentSection == "force")
		{
			stringstream ss(line);
			int number;
			double value;
			char comma;
			ss >> number >> comma >> value;
			force.push_back({number,value});
		}
		if (currentSection == "nodal coordinate")
		{
			stringstream ss(line);
			int number;
			double value;
			char comma;
			ss >> number >> comma >> value;
			nodal_cord.push_back({number,value});
		}
		if (currentSection == "nodal flux")
		{
			stringstream ss(line);
			int number;
			double value;
			char comma;
			ss >> number >> comma >> value;
			nodal_flux.push_back({number,value});
		}
		if (currentSection == "element data")
		{
			stringstream ss(line);
			Element_Data data;
			char comma;
			ss >> data.e >> comma >> data.a >> comma >> data.b >> comma >> data.f >> comma >> data.type;
			if (data.type == "1DC0L,")
			{
				int nb, ne;
				ss >> nb >> comma >> ne;
				data.n = {nb, ne};
			} 
			else if (data.type == "1DC0Q,")
			{
				int nb, nm, ne;
				ss >> nb >> comma >> nm >> comma >> ne;
				data.n = {nb, nm, ne};
			}
			Element.push_back(data);
		}
		if (currentSection == "left end BC")
		{
			stringstream ss(line);
			Boundary_Cond data;
			char comma;
			getline(ss, data.type, ',');
			ss >> data.value_1 >> comma >> data.value_2;
			LBC = data;
		}
		if (currentSection == "right end BC")
		{
			stringstream ss(line);
			Boundary_Cond data;
			char comma;
			getline(ss, data.type, ',');
			ss >> data.value_1 >> comma >> data.value_2;
			RBC = data;
		}		
	}
	inputFile.close();
}

void Process::PrintData()
{
	for(const auto& a : alpha)
	{
		cout << "alpha index " << a.index << ", alpha value " << a.value << endl;
	}
	cout << endl;
	for(const auto& b : beta)
	{
		cout << "beta index " << b.index << ", beta value " << b.value << endl;
	}
	cout << endl;
	for(const auto& f : force)
	{
		cout << "force index " << f.index << ", force value " << f.value << endl;
	}
	cout << endl;
	for(const auto& nc : nodal_cord)
	{
		cout << "nodal_cord index " << nc.index << ", nodal_cord value " << nc.value << endl;
	}
	cout << endl;
	for(const auto& nf : nodal_flux)
	{
		cout << "nodal_flux index " << nf.index << ", nodal_flux value " << nf.value << endl;
	}
	cout << endl;
	for (const auto& element : Element) 
	{
		cout << "Element: " << element.e << ", "
		<< "Alpha: " << element.a << ", "
		<< "Beta: " << element.b << ", "
		<< "Force: " << element.f << ", "
		<< "Type: " << element.type << ", "
		<< "Node: ";
		for (const auto& node : element.n)
		{
			cout << node << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "left BC type: " << LBC.type << " and values " << LBC.value_1 <<", " << LBC.value_2 << endl;
	cout << endl;
	cout << "right BC type: " << RBC.type << " and values " << RBC.value_1 <<", " << RBC.value_2 << endl;
	cout << endl;
}

void Process::Print_Matrix(vector<vector<double>> A)
{
	for (const auto& row : A)
	{
		for (const auto& element : row)
		{
			cout << setw(10) << element << " ";
		}
		cout << endl;
	}
}

void Process::Print_Vector(vector<double> A) {
	for (const auto& element : A)
	{
		cout << element << endl;
	}
	cout << endl;
}

vector<double> Process::Solve(vector<vector<double>> & A, vector<double> & B)
{
	int n = A.size();
	if ( n == 0 || A[0].size() != n || B.size() != n )
	{
		cerr << "Error: Invalid input dimensions!" << endl;
		return {};
	}
	vector<vector<double>> Matrix = A;
	for (int i = 0; i < n; i++)
	{
		Matrix[i].push_back(B[i]);
	}
	Forward(Matrix);
	return Backward(Matrix);
}

tuple< vector<double>, vector<int>, vector<int> > Process::CSR(vector<vector<double>> A)
{
	if ( Symmetric_Cond( A ) )
	{
		int n = A.size();
		vector<double> K;
		vector<int> C, R;
		for (int i = 0; i < n; i++)
		{
			for (int j = i; j < n; j++)
			{
				if ( abs(A[i][j]) >= 1e-9 )
				{
					K.push_back(A[i][j]);
					C.push_back(j);
					if( j == i )
					{
						R.push_back(K.size()-1);
					}
				}
			}
		}
		return {K, C, R};
	}
	else
	{
		cerr << "Error: Invalid input Matrix dimensions or not a symmetric matrix!" << endl;
		return {{}, {}, {}};
	}
}

tuple<vector<vector<double>>, vector<double>> Process::Build()
{
	int n = nodal_cord.size();
	vector<vector<double>> K(n, vector<double>(n));
	vector<double> D(n), F(n);
	for ( int i = 0; i < Element.size(); i++ )
	{
		double L = abs( nodal_cord[ Element[i].n.back()-1 ].value - nodal_cord[ Element[i].n.front()-1 ].value );
		double A = alpha[(Element[i].a-1)].value;
		double B = beta[(Element[i].b-1)].value;
		double f = force[(Element[i].f-1)].value;
		if (Element[i].type == "1DC0L,")
		{
			int nb = Element[i].n[0] - 1, ne = Element[i].n[1] - 1;
			
			K[nb][nb] += A / L + B * L / 3;
			K[ne][ne] += A / L + B * L / 3;
			K[nb][ne] += -A / L + B * L / 6;
			K[ne][nb] += -A / L + B * L / 6;
			
			F[nb] += f * L / 2;
			F[ne] += f * L / 2;
		}
		else if (Element[i].type == "1DC0Q,")
		{
			int nb = Element[i].n[0] - 1, nm = Element[i].n[1] - 1, ne = Element[i].n[2] - 1;

			K[nb][nb] += (7 * A) / (3 * L) + (4* B * L) / 30;
			K[ne][ne] += (16 * A) / (3 * L) + (16 * B * L) / 30;
			K[nm][nm] += (7 * A) / (3 * L) + (4 * B * L) / 30;
			
			K[nb][ne] += (A) / (3 * L) - (B * L) / 30;
			K[ne][nb] += (A) / (3 * L) - (B * L) / 30;
			
			K[nb][nm] += (-8 * A) / (3 * L) + (2 * B * L) / 30;
			K[nm][nb] += (-8 * A) / (3 * L) + (2 * B * L) / 30;
			
			K[ne][nm] += (-8 * A) / (3 * L) + (2 * B * L) / 30;
			K[nm][ne] += (-8 * A) / (3 * L) + (2 * B * L) / 30;
			
			F[nb] += (f * L) / 6;
			F[nm] += (2 * f * L) / 3;
			F[ne] += (f * L) / 6;
		}
	}
	for ( int i = 0; i < nodal_flux.size(); i++ )
	{
		F[ nodal_flux[i].index-1 ] += nodal_flux[i].value;
	}
	return {K, F};
}

bool Process::Symmetric_Cond(vector<vector<double>> A)
{
	// check shape sizes
	int n = A.size(), row, col;
	for (const auto& row : A)
	{
		if (row.size() != n)
		{
			return false;
		}
	}
	
	// checking with transposed
	for (int i = 0; i < n * n; i++)
	{
		row = i / n;
		col = i % n;
		if (col > row && A[row][col] != A[col][row])
		{
			return false;
		}
	}
	return true;	
}

void Process::Forward(vector<vector<double>> & Matrix)
{
	int n = Matrix.size();
	for ( int i = 0; i < n; i++ )
	{
		int maxRow = i;
		double factor; 
		for (int k = i + 1; k < n; k++)
		{
			if ( abs( Matrix[k][i] ) > abs( Matrix[maxRow][i] ) )
			{
				maxRow = k;
			}
		}
		swap( Matrix[i],  Matrix[maxRow] );
		if ( abs( Matrix[i][i] ) < 1e-9 )
		{
			cerr << "Matrix is singular or nearly singular!" << endl;
			return;
		}
		for (int k = i + 1; k < n; k++)
		{
			factor = Matrix[k][i] / Matrix[i][i];
			for (int j = i; j <= n; j++)
			{
				Matrix[k][j] -= factor * Matrix[i][j];
			}
		}
	}
}

vector<double> Process::Backward(vector<vector<double>> & Matrix)
{
	vector<vector<double>> temp = Matrix;
	int n = Matrix.size();
	vector<double> x(n);

	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = temp[i][n] / temp[i][i];
		for (int k = i - 1; k >= 0; k--)
		{
			temp[k][n] -= temp[k][i] * x[i];
		}
	}
	return x;
}
