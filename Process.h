#ifndef PROCESS_H
#define PROCESS_H
#include <string>
#include <vector>
#include <map>

using namespace std;

class Process
{
	public:
		void Get_Debug_Level(int num);
		void Input(const std::string& file_path);		// Get input file path
		void ParseFile();	// Reads the file once and stores data
		void PrintData();	// Prints all input values
		void Print_Matrix(vector<vector<double>> A);	// Print Matrix
		void Print_Vector(vector<double> A);			// Print Vector
		vector<double> Solve(vector<vector<double>> & A, vector<double> & B);
/*		tuple<vector<double>, vector<int>, vector<int>> CSR(vector<vector<double>> A);	// CSR store*/
		tuple<vector<vector<double>>, vector<double>> Build();	// KDF build
		void ABC(vector<vector<double>>& K, vector<double>& F);	// Apply Boundary Conditions
		vector<double> Solution();
		void PrintFormattedOutput(ofstream& outFile, double absError, double relError);
		vector<double> Flux(const vector<double> &Sol);
	private:
		
		struct Element_Value
		{
			// alpha, beta, force, nodal_coordinate, nodel_flux
			int index;
			double value;
		};
		
		struct Element_Data
		{
			int e, a, b, f;
			vector<int> n;
			// element, alpha, beta, force, node
			string type;
			// element type
		};

		struct Boundary_Cond
		{
			double c, d;
			// c, d
			string type;
			// element type
		};
		
		int Debug_Level;
		string file_path;	// Input file path
		
		vector<Element_Value>	alpha;
		vector<Element_Value>	beta;
		vector<Element_Value>	force;
		vector<Element_Value>	nodal_cord;
		vector<Element_Value>	nodal_flux;
		vector<Element_Data>	Element;
		Boundary_Cond LBC, RBC;

		static void Forward(vector<vector<double>> & Matrix);
		static vector<double> Backward(vector<vector<double>> & Matrix);
};

#endif // PROCESS_H
