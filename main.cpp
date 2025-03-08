# include <iostream>
# include <fstream>
# include <limits>
# include <algorithm>
# include <cctype>
# include <string>
# include <cmath>
# include "Process.h"
using namespace std;

std::string trim(const std::string& str) {
	auto start = str.begin();
	while (start != str.end() && std::isspace(*start))
	{
		start++;
	}
	auto end = str.end();
	while (end != start && std::isspace(*(end - 1)))
	{
		end--;
	}
	return std::string(start, end);
}

double computeL2Norm(const vector<double>& vec) {
	double norm = 0.0;
	for (const auto& value : vec)
	{
		norm += value * value;
	}
	return sqrt(norm);
}

pair<double, double> computeErrors(const vector<double>& R, const vector<double>& F)
{
	double absError = computeL2Norm(R); // Absolute error: ||R||
	double relError = absError / computeL2Norm(F); // Relative error: ||R|| / ||F||
	return {absError, relError};
}

int main()
{
	// input for Debug Level
	Process test;
	int Debug_Level = 0;
	string inputFilePath, outputFilePath;
	
	// cout << "Enter debug level" << endl;
	// cout << "Enter 0: Regular output" << endl;
	// cout << "Enter 1: Regular output with System stiffness matrix and the system load vector." << endl;
	// cout << "Enter 2: Regular output Element stiffness matrix, element load vector, system  stiffness matrix and the system load vector." << endl;
	// cin >> Debug_Level;
	test.Get_Debug_Level(Debug_Level);
	
	//cin.ignore(numeric_limits<streamsize>::max(), '\n');
	
	// Input for input file 
	cout << "Enter input File name, like Input_Name.txt" << endl;
	getline(cin, inputFilePath);
	inputFilePath = "./" + trim(inputFilePath);
	
	ifstream inputFile(inputFilePath);
	while (!inputFile.is_open())
	{
		cerr << "Unable to open input file: " << inputFilePath << endl;
		cout << "Please enter a valid input file name: ";
		getline(cin, inputFilePath);
		inputFilePath = "./" + trim(inputFilePath);
		inputFile.open(inputFilePath);
	}
	inputFile.close();
	
	// Input for output file
	cout << "Enter output File name, like Output_Name.out" << endl;
	getline(cin, outputFilePath);
	outputFilePath = "./" + trim(outputFilePath);
	
	ofstream outFile(outputFilePath);
	if (!outFile.is_open()) {
		cerr << "Error: Unable to open output file: " << outputFilePath << endl;
		return 1;
	}
	
	// Computation
	test.Input( inputFilePath );
//	test.PrintData();
	auto [K, F] = test.Build();
	vector<double> D;
	D = test.Solution();
	vector<FluxResult> flux;
	flux = test.Flux(D);
	
	vector<double> R(F.size(), 0.0);
	for (size_t i = 0; i < K.size(); ++i)
	{
		for (size_t j = 0; j < K[i].size(); ++j)
		{
			R[i] += K[i][j] * D[j];
		}
		R[i] -= F[i]; // R = KD - F
	}
	auto [absError, relError] = computeErrors(R, F);

	// Print formatted output to the file
	test.PrintFormattedOutput(outFile, absError, relError);
	
	// Close the file
	outFile.close();
	
	// Print 
	ifstream inFile(outputFilePath);
	if (!inFile.is_open()) {
		cerr << "Error: Unable to open output file for reading: " << outputFilePath << endl;
		return 1; // Exit the program with an error code
	}

	cout << "\nContents of the output file:\n";
	cout << "----------------------------------------\n";
	string line;
	while (getline(inFile, line))
	{
	cout << line << endl;
	}
	cout << "----------------------------------------\n";

	// Close the file
	inFile.close();
	cout << "Output written to: " << outputFilePath << endl;

	return 0;
}
