#include<iostream>
#include<string.h>
#include<cmath>
#include<stdlib.h>
#include<stdio.h>
#include "../matrix.h"
#include<fstream>
using namespace std;
int main(int argc, char *argv[]){
	ofstream myfile;
	myfile.open("testm.m", ios::out);

	int _r = atoi(argv[1]);
	int c = atoi(argv[2]);
	Matrix<double> A(_r, c);
	cout << _r << " " << c << endl;
	srand((unsigned)time(NULL));
	for(int i = 0 ; i< _r; i++){
		for(int j = 0 ;j  <c ;j++){
			double r = double(rand()%10);
			A[i][j] = r;
		}
	}
	
	for(int i = 0 ; i< _r; i++){
		double max = 0.0;
		for(int j = 0 ;j  <c ;j++){
			if(max < A[i][j]){
				max = A[i][j];
			}
		}
		A[i][i] = max * _r + 1.0;
	}
	
	for(int i = 0 ; i< _r; i++){
		for(int j = 0 ;j  <c ;j++){
			cout << A[i][j] << " " ;
		}
		cout << endl;
	}
	
	cout << endl;
	vector<double> b(_r, 0);
	for(int i = 0 ; i< _r; i++){
			double r = double(rand());
			b[i] = r;
			cout << b[i] << " ";
	}
	
	myfile << "A = [";
	for(int i = 0 ; i< _r; i++){

		for(int j = 0 ;j  <c ;j++){
			if(j != c-1){
				myfile << A[i][j] <<",";
			}else{
				myfile << A[i][j];
			}
		}
		myfile << ";" << endl;
	}
	myfile <<"];" << endl;
	myfile << "b = [";	
	for(int i = 0 ; i< _r; i++){
		
		if( i != _r -1){
			myfile << b[i] <<",";
		}else{
			myfile << b[i];
		}
	}
	myfile<<"];" << endl;
	myfile << "A\\(b')" << endl;
	cout << endl;
	
	return 0;
}

