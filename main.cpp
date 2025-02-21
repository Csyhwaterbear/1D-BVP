# include <iostream>
# include <fstream>
# include <limits>
# include <algorithm>
# include <cctype>
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

int main()
{
	// input 
	int Debug_Level;
	string inputFilePath, outputFilePath;
	
	cout << "Enter debug level" << endl;
	cout << "Enter 0: Regular output" << endl;
	cout << "Enter 1: Regular output with System stiffness matrix and the system load vector." << endl;
	cout << "Enter 2: Regular output Element stiffness matrix, element load vector, system  stiffness matrix and the system load vector." << endl;
	cin >> Debug_Level;
	
	cin.ignore(numeric_limits<streamsize>::max(), '\n');
	
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

	cout << "Enter output File name, like Output_Name.out" << endl;
	getline(cin, outputFilePath);
	outputFilePath = "./" + trim(outputFilePath);
	
	
	Process test1;
	
	test1.Input( inputFilePath );
	test1.PrintData();

	vector<vector<double>> K;
	vector<double> F, D, flux;
	tie(K, F) = test1.Build();
	D = test1.Solution();
	flux = test1.Flux(D);

	// Create output file name
	size_t lastDot = inputFilePath.find_last_of(".");
	// string outputFilePath = inputFilePath.substr(0, lastDot) + ".out";
	
	// Open the output file
	ofstream outFile(outputFilePath);
	if (!outFile.is_open())
	{
		cerr << "Error: Unable to open output file: " << outputFilePath << endl;
		return 1;
	}

	// Print formatted output to the file
	test1.PrintFormattedOutput(outFile);
	
	// Close the file
	outFile.close();
	
	cout << "Output written to: " << outputFilePath << endl;

	return 0;
}
