#ifndef _PDF_H_
#define _PDF_H_
#include "matrix.h"
#include "math.h"
//transition probility with non absorbing 
template<typename Type>
const Matrix<Type> nonAbsorbing(const Matrix<Type> &A, Matrix<Type> &P){
	int row = A.rows();
	int col = A.cols();
	Matrix<Type> T(row, col);
	Type sum = 0.0;
	for(int i = 0 ; i < row ; i++){
		sum = 0.0;
		for(int j = 0 ;j < col; j++){
			sum += abs(A[i][j]);
		}
		for(int j = 0 ; j < col; j++){

			P[i][j] = abs(A[i][j]) / sum;
		}
	}
	for(int i = 0 ; i < row; i++){
		T[i][0] = P[i][0];
		for(int j = 1 ;j <col; j++){
			T[i][j] = T[i][j-1] + P[i][j];
		}
		T[i][col-1] = 1.0;
	}
	return T;
}
template<typename Type>
const Matrix<Type> Absorbing(const Matrix<Type> &A, double scale, Matrix<Type> &P){
	int row = A.rows();
	int col = A.cols();
	Matrix<Type> T(row + 1, col +1);
	Type sum = 0.0;
	for(int i = 0 ; i < row ; i++){
		sum = 0.0;
		double sum1 = 0.0;
		for(int j = 0 ;j < col; j++){
			sum += abs(A[i][j]);
		}
		sum = double((double(col) + scale)) / col * sum;
		for(int j = 0 ; j < col; j++){
			P[i][j] = abs(A[i][j]) / sum;
			sum1 += P[i][j];
		}
		P[i][col] = 1 - sum1;
	}
	for(int i = 0 ; i < row; i++){
		T[i][0] = P[i][0];
		for(int j = 1 ;j <col; j++){
			T[i][j] = T[i][j-1] + P[i][j];
		}
		T[i][col] = 1.0;
	}
	return T;
}

template<typename Type>
const Matrix<Type> backwards_p(const Matrix<Type> &A, double scale, Matrix<Type> &P){
	int row = A.rows();
	int col = A.cols();
	Matrix<Type> B(row, col);
	for(int i = 0 ; i < row; i++){
		for(int j = 0 ; j < col; j++){
			B[i][j] = A[j][i];
		}
	}
	Matrix<Type> T(row + 1, col +1);
	Type sum = 0.0;
	for(int i = 0 ; i < row ; i++){
		sum = 0.0;
		double sum1 = 0.0;
		for(int j = 0 ;j < col; j++){
			sum += abs(B[i][j]);
		}
		sum = double((double(col) + scale)) / col * sum;
		for(int j = 0 ; j < col; j++){
			P[i][j] = abs(B[i][j]) / sum;
			sum1 += P[i][j];
		}
		P[i][col] = 1 - sum1;
	}
	for(int i = 0 ; i < row; i++){
		T[i][0] = P[i][0];
		for(int j = 1 ;j <col; j++){
			T[i][j] = T[i][j-1] + P[i][j];
		}
		T[i][col] = 1.0;
	}
	return T;
}
#endif
