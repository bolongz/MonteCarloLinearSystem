#ifndef _NEUMANN_H_
#define _NEUEMAN_H_

#include "matrix.h"
#include "pdf.h"
#include <time.h>
#include<stdlib.h>
#include<algorithm>

class Neumann{
private:
	size_t row;
	size_t col;
	double err, sum1, sum2, x,err_w;
	int next, step, times;
public:
	void init(){
		err = 10000.0;
		sum1 = 0.0; sum2 = 0.0 ;x = 0.0;
		next  = 0; step = 100; times = 1;
		err_w = 1e-6;
	}
	template <typename Type>
	std::vector<Type> absorbing(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1){
		row = A.rows();col = A.cols();
		Matrix<Type> P(row +1, col +1);
		Matrix<Type> t = Absorbing(A, 0.2, P);
		size_t size = b.size();
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));	
		for(int i = 0 ; i < size; i++){
			init();
			/*double err = 10000.0;
			double sum1 = 0.0, sum2 = 0.0, x = 0.0;
			int next = 0, step = 100, times = 1;
			cout << "step i" << endl;
			*/while( err > _err){

				while(step--){
					double v = 1.0;
					int index = i, next = 0;
					while(next != col){
						double r = double(rand())/RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) -t[index].begin();
					/*	cout << " ---------------------" << endl;
						cout << r << "  " <<" " << next << " "  <<  t[index][next] << "  " << t[index][next-1] << endl;
						cout << " ---------------------" << endl;
					*/	if(next == col) continue;	
						v = v * A[index][next] /  P[index][next];
						index = next;
					}
					v = v * b[index] /P[index][col];
					sum1 += v;
					sum2 += v * v ;
				}

				step = 100;
				int total = step * times;
				x = sum1 / total;
				double _err = (sum2 - sum1 / total) / total / total ;
				times++;
				err = sqrt(_err) / x; 
			//	cout <<total << "  "  <<  err <<  "    "   << x << " " << sum1  << endl;
			}
			res[i] = x;
		}
		return res;
	}
	template <typename Type>
	std::vector<Type> nonabsorbing(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1){
		row = A.rows(); col = A.cols();
		Matrix<Type> P(row, col);
		Matrix<Type> t = nonAbsorbing(A, P);
		size_t size = b.size();
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));	
		for(int i = 0 ; i < size; i++){
			init();
			/*double err = 10000.0;
			double sum1 = 0.0, sum2 = 0.0, x = 0.0, ee = 0.000001;
			int next = 0, step = 100 ,times = 1;
			cout << "step i" << endl;
			*/while( err > _err){
				while(step--){
					double v = 0.0, w = 1.0;
					int index = i, next = 0;
					//cout << "start" << endl;
					while(abs(w) > err_w){
						double r = double(rand())/RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) -t[index].begin();
						w = w * A[index][next] / P[index][next];
	/*					cout << " ---------------------" << endl;
						cout << r << "  " <<" " << next << " "  <<  t[index][next] << "  " << t[index][next-1] << endl;
						cout << " ---------------------" << endl;
				
	*/					v = v + w * b[next];
						//cout << v << endl;
						index = next;
					}
					v = v + b[i];
					//cout << "break" << endl;
					sum1 += v;
					sum2 += v * v ;
				}
				step = 100;
				int total = step * times;
				x = sum1 / total;
				double _err = (sum2 - sum1 / total) / total / total ;
				times++;
				err = sqrt(_err) / x; 
				//cout <<total << "  "  <<  err <<  "    "   << x << " " << sum1  << endl;
			}
			res[i] = x;
		}
		return res;
	}
};
#endif
