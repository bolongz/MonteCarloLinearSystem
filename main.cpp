


#include "matrix.h"
#include "power.h"
//#include "neumann.h"
#include "pdf.h"
#include <cassert>
#include<stdio.h>
#include "stdlib.h"
#include "neumann.h"
int main(int argc, char *argv[]){
	freopen("../in.txt", "r", stdin);
	
	int size = 3;
	Matrix<double> A(size,size);
	for(int i = 0; i < size; i++){
		vector<double> a;
		for(int j = 0 ; j < size; j++){
			double input;
			scanf("%lf",  &input);
			a.push_back(input);
		}
		A.setRow(i, a);
	}
	vector<double> b;
	for(int i = 0 ;i < size; i++){
		double input;
		scanf("%lf",  &input);
		b.push_back(input);
	}
	printmatrix(A);

	for(int i = 0 ; i < size; i++){

		cout << b[i] << "  "  << endl;
	}
	Matrix<double> I(size,size);
	I.unit();
	A = I - A;
	cout << " H " << endl;
	printmatrix(A);
	double lamda = power(A);
	cout << lamda << endl;
	if(abs(lamda) >= 1){
		cout << "Can not convergence !" << endl;
		exit(1);
	}
	Neumann neumann1;
	std::vector<double> res = neumann1.absorbing(A,b);
	cout << " final----------------" << endl;
	for(int i = 0 ; i  <res.size(); i++){

		cout << res[i] << endl;
	}
	

	Neumann neumann;

	std::vector<double> res1 = neumann.nonabsorbing(A,b);
	cout << " final----------------" << endl;
	for(int i = 0 ; i  <res1.size(); i++){

		cout << res1[i] << endl;
	}
	return 0;
}
