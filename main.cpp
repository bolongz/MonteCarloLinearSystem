/* Due to the generation problems. Here we only test strictly diagonal dominant problems */

#include "matrix.h"
#include "power.h"
#include "pdf.h"
#include <cassert>
#include<stdio.h>
#include "stdlib.h"
#include <string.h>
#include "neumann.h"
#include "time.h"

int main(int argc, char *argv[]){

	FILE *equ = fopen(argv[2], "r");
	if(equ == NULL){
		printf("errori\n");
		exit(1);
	}
	int r, c; //define the row and col;
	fscanf(equ, "%d %d", &r, &c);
	Matrix<double> A(r, c);
	for(int i = 0 ; i < r ;i++){
		for(int j = 0 ; j < c; j++){
			fscanf(equ, "%lf", &A[i][j]);
		}
	}
	vector<double> b(r, 0);
	for(int i = 0 ; i < r ;i++){
		fscanf(equ, "%lf", &b[i]);
	}
	bool _sdd = true;
	for(int i = 0; i< r; i++){
		double sum = 0.0;
		for(int j = 0; j< r; j++){
			if(i != j){
				sum += abs(A[i][j]);
			}
		}
		if(abs(A[i][i]) < sum){
			cout << "not sdd " << endl;
			_sdd = false;
		}

	}

	if(!_sdd){
		Matrix<double> I(r,c); //define the unit matrix
		I.unit();
		A = I - A; // calculate the H;
		double lamda = power(A); // check the largest eigenvalue for the matrix by power method
		if(abs(lamda) >= 1){
			cout << "Can not convergence !" << endl;
			exit(1);
		}
	}else{
		for(int i = 0 ; i<r;i++){
			for(int j = 0; j < r; j++){
				if(i != j)
					A[i][j] = -A[i][j] / A[i][i];
			}
			b[i] = b[i]/A[i][i];
			A[i][i] = 0.0;
		}
	}
	
	//printmatrix(A);
	double p;
/*	if(argc>= 4){
		cout << "wrong parameters: You should input the expectation precision fistly" << endl;
		exit(1);
	}else */
	int nwalks = 0;
	if(argc < 4){
		p = 0.01;
	}else{
		p = atof(argv[4]);
		nwalks = atoi(argv[4]);
	}
	
	Timer tmr;
	Neumann<double> neumann;
	std::vector<double> res;
	if(argc >= 6 && argv[5][1] == 'm'){
		if(argv[6][0] == 'n'){
			cout << "Method with nonAbsorbing Matrix:" << endl;
			res = neumann.nonabsorbing(A,b, p);
		}else if(argv[6][0] == 'a'){
			cout << "Method with Absorbing Matrix:" << endl;
			res = neumann.absorbing(A,b, p);
		}else if (argv[6][0] == 'B') {
			cout << "Method with Absorbing Matrix using threads:" << endl;
			res = neumann.absorbing_UsingThreads(A, b, p);
		}else if (argv[6][0] == 'c') {
			cout << "Method with nonAbsorbing Matrix using threads:" << endl;
			res = neumann.nonabsorbing_UsingThreads(A, b, p);
		}else if(argv[6][0] == 'b'){
			res = neumann.backwards(A,b,nwalks);
		}else{
			cout << "Wrong input type paramenters" << endl;
		}
	}else{
			cout << "Method with Absorbing Matrix:" << endl;
			res = neumann.absorbing(A,b, p);
	}

	cout << "Solution for the Matrix: " << endl;

	for(int i = 0 ; i  < (int)res.size(); i++){
			//cout << res[i] << endl;
			printf(" %+.4e\n",  res[i]);;
	}

	cout <<"Time Cost: " << tmr.elapsed() <<"s" << endl;
	
	return 0;

}
