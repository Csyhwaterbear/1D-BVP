#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ctime>
#include <limits>
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
			ss >> data.c >> comma >> data.d;
			LBC = data;
		}
		if (currentSection == "right end BC")
		{
			stringstream ss(line);
			Boundary_Cond data;
			char comma;
			getline(ss, data.type, ',');
			ss >> data.c >> comma >> data.d;
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
		<< "Type: " << element.type << " Node: ";
		for (const auto& node : element.n)
		{
			cout << node << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "left BC type: " << LBC.type << " and values " << LBC.c <<", " << LBC.d << endl;
	cout << endl;
	cout << "right BC type: " << RBC.type << " and values " << RBC.c <<", " << RBC.d << endl;
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
	cout << endl;
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
			//L /= 2.0;
			int nb = Element[i].n[0] - 1, nm = Element[i].n[1] - 1, ne = Element[i].n[2] - 1;
			// textbook (T4L4-7)
			K[nb][nb] += (7 * A) / (3 * L) + (4* B * L) / 30;
			K[nm][nm] += (16 * A) / (3 * L) + (16 * B * L) / 30;
			K[ne][ne] += (7 * A) / (3 * L) + (4 * B * L) / 30;
			
			K[nb][ne] += (A) / (3 * L) - (B * L) / 30;
			K[ne][nb] += (A) / (3 * L) - (B * L) / 30;
			
			K[nb][nm] += (-8 * A) / (3 * L) + (2 * B * L) / 30;
			K[nm][nb] += (-8 * A) / (3 * L) + (2 * B * L) / 30;
			
			K[ne][nm] += (-8 * A) / (3 * L) + (2 * B * L) / 30;
			K[nm][ne] += (-8 * A) / (3 * L) + (2 * B * L) / 30;
			
			F[nb] += (1 * f * L) / 6;
			F[nm] += (4 * f * L) / 6;
			F[ne] += (1 * f * L) / 6;
		}
	}
	for ( int i = 0; i < nodal_flux.size(); i++ )
	{
		F[ nodal_flux[i].index-1 ] += nodal_flux[i].value;
	}
	Print_Matrix(K);
	Print_Vector(F);
	ABC(K, F, LBC, RBC);
	Print_Matrix(K);
	Print_Vector(F);
	return {K, F};
}

vector<double> Process::Solution()
{
	auto [K, F] = Build();
	vector<double> solution = Solve(K, F);
	if (!solution.empty())
	{
		if (LBC.type == "EBC")
        {
            solution.insert(solution.begin(), LBC.c);
        }
        if (RBC.type == "EBC")
        {
            solution.push_back(RBC.c);
        }
	}
	else
	{
		cerr << "Error: Solution vector is empty!" << endl;
	}
	return solution;
}

void Process::PrintFormattedOutput(ofstream& outFile)
{
	outFile << "		           1D Boundary Value Program         " << endl;
	outFile << "		         Finite Elements for Engineers       " << endl;
	outFile << "		          (c) 2025, S. CSYH                  " << endl;
	outFile << "		              Version: 1.0. Rev: 1           " << endl;
	outFile << "		---------------------------------------------" << endl;
	outFile << endl;
	
	// Print details of the finite element model
	outFile << "          DETAILS OF FINITE ELEMENT MODEL" << endl;
 	outFile << "          -------------------------------" << endl;
	outFile << "               Number of Nodes: " << nodal_cord.size() << endl;
	outFile << "            Number of Elements: " << Element.size() << endl;
	outFile << endl;
	
	// Print nodal data
	outFile << "                                  Nodal Data" << endl;
	outFile << "                           -------------------------" << endl;
	outFile << "                           Node      X Coor         " << endl;
	outFile << "                           -------------------------" << endl;
	for (const auto& nc : nodal_cord)
	{
		outFile << std::left << "                           " << std::setw(20) << nc.index <<  nc.value << endl;
	}
	outFile << endl;

	// Print nodal boundary conditions
	outFile << "                           Nodal Boundary Conditions" << endl;
	outFile << "               --------------------------------------------------" << endl;
	outFile << "               Node      BC        Value 1        Value 2        " << endl;
	outFile << "               --------------------------------------------------" << endl;
	if (LBC.type == "MBC")
	{
        outFile << std::left << "                " << std::setw(10) << 1 << std::setw(14) << LBC.type << std::setw(15) << LBC.c << std::setw(15) << LBC.d << std::endl;
    }
    else if (LBC.type == "EBC")
    {
        outFile << std::left << "                " << std::setw(10) << 1 << std::setw(14) << LBC.type << std::setw(15) << LBC.c << std::setw(15) << " " << std::endl;
    }
    else if (LBC.type == "NBC")
    {
        outFile << std::left << "                " << std::setw(10) << 1 << std::setw(14) << LBC.type << std::setw(15) << " " << std::setw(15) << LBC.d << std::endl;
    }
	if (RBC.type == "MBC")
	{
        outFile << std::left << "                " << std::setw(10) << 1 << std::setw(14) << RBC.type << std::setw(15) << RBC.c << std::setw(15) << RBC.d << std::endl;
    }
    else if (RBC.type == "EBC")
    {
        outFile << std::left << "                " << std::setw(10) << 1 << std::setw(14) << RBC.type << std::setw(15) << RBC.c << std::setw(15) << " " << std::endl;
    }
    else if (RBC.type == "NBC")
    {
        outFile << std::left << "                " << std::setw(10) << 1 << std::setw(14) << RBC.type << std::setw(15) << " " << std::setw(15) << RBC.d << std::endl;
    }
	outFile << endl;
	
	// Print applied nodal flux
	outFile << "                              Applied Nodal Flux" << endl;
	outFile << "                           -------------------------" << endl;
	outFile << "                           Node      Flux           " << endl;
	outFile << "                           -------------------------" << endl;
	for (const auto& nf : nodal_flux)
	{
		outFile << "                           " << nf.index << "         " << nf.value << "              " << endl;
	}
	outFile << endl;
	
	// Print element data
	outFile << "                                 Element Data" << endl;
	outFile << "       -----------------------------------------------------------------" << endl;
	outFile << "       Element   Type           Int. Order     List of Nodes ...        " << endl;
	outFile << "       -----------------------------------------------------------------" << endl;
	for (const auto& element : Element)
	{
		outFile << "       " << element.e << "         ";
		if ( element.type == "1DC0L,")
		{
			outFile << "C0-Linear" << "      " << 2 << "              ";
		}
		else if ( element.type == "1DC0Q,")
		{
			outFile << "C0-Quadratic" << "   " << 3 << "              ";
		}
		for (const auto& node : element.n)
		{
			outFile << node << "    ";
		}
		outFile << endl;
	}
	outFile << endl;

	// Print element properties
	outFile << "                               Element Properties" << endl;
	outFile << "     ----------------------------------------------------------------------" << endl;
	outFile << "     Element   Qty            x**2           x              constant       " << endl;
	outFile << "     ----------------------------------------------------------------------" << endl;
	for (const auto& element : Element)
	{
		outFile << "    " << element.e << "          Alpha          0              0              " << alpha[element.a - 1].value << "          " << endl;
		outFile << "               Beta             0              0              " << beta[element.b - 1].value << "          " << endl;
		outFile << "               f                    0              0              " << force[element.f - 1].value << "          " << endl;
	}
	outFile << endl;

	// Print nodal values
	vector<double> solution = Solution();
	outFile << "                                  Nodal Values" << endl;
	outFile << "                    ----------------------------------------" << endl;
	outFile << "                    Node      Location       Value          " << endl;
	outFile << "                    ----------------------------------------" << endl;
	for (size_t i = 0; i < solution.size(); ++i)
	{
		outFile << "                    " << i + 1 << "             " << nodal_cord[i].value << "                " << solution[i] << "          " << endl;
	}
	outFile << endl;
	
	// Print element flux
	vector<double> flux = Flux(solution);
	outFile << "                                  Element Flux" << endl;
	outFile << "                    ----------------------------------------" << endl;
	outFile << "                    Element   Location       Flux           " << endl;
	outFile << "                    ----------------------------------------" << endl;
	for (size_t i = 0; i < flux.size(); ++i)
	{
		outFile << "                    " << i + 1 << "               " << nodal_cord[i].value << "              " << flux[i] << "          " << endl;
	}
	outFile << endl;
}

vector<double> Process::Flux(const vector<double> &Sol)
{
	int n = Sol.size();
	vector<double> flux(Element.size(), 0.0);
	for (int i = 0; i < Element.size(); i++)
	{
		double L = abs( nodal_cord[ Element[i].n.back()-1 ].value - nodal_cord[ Element[i].n.front()-1 ].value );
		double A = alpha[(Element[i].a-1)].value;
		flux[i] = -(Sol[i+1] - Sol[i]) / L * A;
	}
	return flux;
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

void Process::ABC(vector<vector<double>>& K, vector<double>& F, const Boundary_Cond& LBC, const Boundary_Cond& RBC)
// Apply Boundary Conditions
{
	int n = K.size();
	// left Boundary Condition
	if (LBC.type == "EBC")
	{
        	for (int i = 0; i < n; i++)
        	{
			F[i] -= K[i][0] * LBC.c;
		}
		// Remove row and column for node 1
		K.erase(K.begin());
		for (auto& row : K) row.erase(row.begin());
		F.erase(F.begin());
	}
	else if (LBC.type == "NBC")
	{
		// Natural BC: Add LBC.c to F_1
		F[0] += LBC.d;
	}
	else if (LBC.type == "MBC")
	{
		// Mixed BC: Modify K and F
		K[0][0] -= LBC.c; // c_a
		F[0] += LBC.d;    // d_a * c
	}
	// Right Boundary Condition
	if (RBC.type == "EBC")
	{
	 	for (int i = 0; i < n-1; i++)
		{
			F[i] -= K[i][n-1] * RBC.c;
		}
		// Remove row and column for node n
		K.pop_back();
		for (auto& row : K) row.pop_back();
		F.pop_back();
	}
	else if (RBC.type == "NBC")
	{
		// Natural BC: Add RBC.c to F_n
		F[n-1] -= RBC.d;
	}
	else if (RBC.type == "MBC")
	{
		// Mixed BC: Modify K and F
		K[n-1][n-1] += RBC.c; // c_b
		F[n-1] -= RBC.d;      // d_b * d
	}
}
