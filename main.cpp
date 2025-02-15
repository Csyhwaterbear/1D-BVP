# include <iostream>
# include "Process.h"
using namespace std;

int main()
{
	Process test1;
	test1.Input("./Test Model 2.txt");
	test1.PrintData();

	vector<vector<double>> A = {
		{2.2, 1.1, -19.1},
		{-3.2 , -1.0, 2.3},
		{-2.3, 1.1, 2.4}
	};
	vector<double> B = {8.8, -11.7, -3.5};
	test1.Print_Matrix(A);
	vector<double> X = test1.Solve(A, B);
	if ( !X.empty() )
	{
		cout << "Solution vector X:" << endl;
		for (const auto& x : X)
		{
			cout << x << endl;
		}
	}
	
	vector<vector<double>> C = {
		{1.0, 2.0, 3.0},
		{2.0, 4.0, 5.0},
		{3.0, 5.0, 6.0}
	};
	cout << "check symmetric with tolence less than 1e-9 " << test1.Symmetric_Cond(C) << endl << endl;
	/*
	vector<double> K;
	vector<int> Q,P;
	tie(K, Q, P) = test1.CSR(C);
	// OR auto [K, Q, P] = Sol.CSR(C);
	cout << "K: ";
	for (double val : K) cout << val << " ";
	cout << endl;

	cout << "Q: ";
	for (int val : Q) cout << val << " ";
	cout << endl;

	cout << "P: ";
	for (int val : P) cout << val << " ";
	cout << endl;

	cout << endl;
	*/
	vector<vector<double>> K;
	vector<double> F;
	tie(K, F) = test1.Build();
	test1.Print_Matrix(K);
	test1.Print_Vector(F);
	test1.Print_Vector(test1.Solve(K, F));


	return 0;
}
