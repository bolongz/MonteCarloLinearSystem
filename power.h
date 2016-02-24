
#include "matrix.h"
#include <cmath>

template<typename Type>

Type power(const Matrix<Type> &A){

	int size = A.cols();
	vector<Type> v(size,1), u(size,1);
	Type lamda, lamda2 = 1000, err = 1000;
	while( err > 1e-5){
		v = A * u;
		lamda = v[0];
		for(int i = 1 ; i < size; i++){
			if( abs(lamda)  < abs(v[i])){
				lamda = v[i];
			}
		}
		for(int i = 0 ; i < size; i++){
			u[i] = v[i] / lamda;
		}
		err = lamda2 - lamda;
		lamda2 = lamda;
	}
	return lamda;
}

