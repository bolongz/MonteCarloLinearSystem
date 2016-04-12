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
	
//	printmatrix(A);
	
	Matrix<double> I(r,c); //define the unit matrix
	I.unit();
	
	A = I - A; // calculate the H;
	
	double lamda = power(A); // check the largest eigenvalue for the matrix by power method
	
	if(abs(lamda) >= 1){
		cout << "Can not convergence !" << endl;
		exit(1);
	}


	double p;
/*	if(argc>= 4){

		cout << "wrong parameters: You should input the expectation precision fistly" << endl;
		exit(1);
	}else */
	if(argc < 4){
		p = 0.01;
	}else{
		p = atof(argv[4]);
	}
	
	Timer tmr;
	Neumann neumann;
	std::vector<double> res;
	if(argc >= 6 && argv[5][1] == 'm'){
		if(argv[6][0] == 'n'){
			cout << "Method with nonAbsorbing Matrix:" << endl;
			res = neumann.nonabsorbing(A,b, p);
		}else if(argv[6][0] == 'a'){
			cout << "Method with Absorbing Matrix:" << endl;

			res = neumann.absorbing(A,b, p);
		}else{
			cout << "Wrong input type paramenters" << endl;
		}
	}else if(argv[4][0] == 'b'){
		
			res = neumann.backwards(A,b, p);
	}else{
			cout << "Method with Absorbing Matrix:" << endl;
			res = neumann.absorbing(A,b, p);
	}

	cout << "Solution for the Matrix: " << endl;

	for(int i = 0 ; i  <res.size(); i++){
			cout << res[i] << endl;
	}

	cout <<"Time Cost: " << tmr.elapsed() <<"s" << endl;
	
	return 0;

}

	/*string filename(argv[1]);
	ifstream fin(filename.c_str(), ios::in);
	string line;
	while(getline(fin,line)){
		size_t pos1 = 0,pos2 = 0;
		while(true){
			pos1 = line.find_first_of("[,;", pos2);
			if(pos1 == string::npos) break;
			pos2 = line.find_first_of(",;]", pos1 + 1);
			if(pos2 == string::npos) break;
			if(line[pos2]
			string p = line.substr(pos1 + 1, pos2 - pos1 -1);
			points.push_back(atof(p.c_str()));

		}
	}

	*/


	
	
	/*	freopen("in.txt", "r", stdin);
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
*/

/*
#define line_size 512
int main(int argc, char *argv[])
{
	FILE *input;
	
	char *token;

	if (argc < 2)
	{
		Matrix<double> A(3, 3);
		printmatrix(A);
	}
	else
	{
		if ((inFPtr = fopen(argv[1], "r")) == NULL)
		{
			cout << "Error: unable to open the inputed file.\n";
			exit(0);
		}

		Matrix<double> A(0, 0);
		while (fgets(line, line_size, inFPtr) != NULL)
		{
			if (line[0] != '#')
			{
				vector<double> vect;
				token = strtok(line, ", ");
				while (token != NULL)
				{
					vect.push_back(strtod(token, NULL));
					token = strtok(NULL, ", ");
				}
				A.addRow(vect);
			}
		}
		printmatrix(A);
	}
	return 0;
}*/
