

#include "matrix.h"
#include "power.h"
#include <cassert>
int main(){
	Matrix<double> A(3,3);
	printmatrix(A);
	for(int i = 1 ; i <= 3; i++){
		vector<double> a;
		for(int j = 1 ; j <=3; j++){
			a.push_back(i* j);
			A.setRow(i-1, a);
		}
	}
	printmatrix(A);
	Matrix<double> B(3,3);
	for(int i = 1 ; i <= 3; i++){
		vector<double> a;
		for(int j = 1 ; j <=3; j++){
			a.push_back(1);
			B.setRow(i-1, a);
		}
	}
	Matrix<double> C = A *B;
	printmatrix(C);
	Matrix<double> D = A + B;
	printmatrix(D);
	Matrix<double> E = A - B;
	printmatrix(E);
	Matrix<double> F = A * 2;
	printmatrix(F);
	Matrix<double> G = A / 2.0;
	printmatrix(G);
	Matrix<double> H = 2.0/A;
	printmatrix(H);

	vector<double> a(3);
	a[0] =1 ; a[1] = 1; a[2] = 1;

	vector<double> I = A * a;
	
	cout << I[0] << "  " <<I[1] << "  " << I[2] << endl;
	
	vector<double> J = a *A;
	
	cout << J[0] << "  " <<J[1] << "  " << J[2] << endl;

	cout <<"max  " << power(A) << endl;
	Matrix<double> AA(2,2);
	vector<double> aa(2);
	aa[0] = 3; aa[1] = 1;
	AA.setRow(0,aa);
	vector<double> b(2);
	b[0] = 3; b[1] = 1;
	AA.setRow(1, b);
	cout <<"max  " << power(AA) << endl;
	return 0;
}
