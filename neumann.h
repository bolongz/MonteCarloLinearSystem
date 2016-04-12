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
		next  = 0; step = 20000; times = 1;
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
		int total;
		for(int i = 0 ; i < size; i++){
			init();
			while( err > _err){
				while(step--){
					double v = 1.0;
					int index = i, next = 0;
					while(next != col){
						double r = double(rand())/RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) -t[index].begin();
						if(next == col) continue;	
						v = v * A[index][next] /  P[index][next];
						index = next;
					}
					v = v * b[index] /P[index][col];
					sum1 += v;
					sum2 += v * v ;
				}

				step = 100;
				total = step * times;
				if(total % 200000 == 0){
					cout << "Calculating x[" << i<< "]: " << total  << " Random walks generated" << endl;
				}
				x = sum1 / total;
				double __err = (sum2 - sum1 / total) / total / total ;
				times++;
				err = sqrt(__err) / x; 
			}
			cout << "Calculating x[" << i<< "]: " << total  << " Random walks generated" << endl;
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
		int total;
		for(int i = 0 ; i < size; i++){
			init();
			while( err > _err){
				while(step--){
					double v = 0.0, w = 1.0;
					int index = i, next = 0;
					while(abs(w) > err_w){
						double r = double(rand())/RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) -t[index].begin();
						w = w * A[index][next] / P[index][next];
						v = v + w * b[next];
						index = next;
					}
					v = v + b[i];
					sum1 += v;
					sum2 += v * v ;
				}
				step = 100;
				total = step * times;
				if(total % 200000 == 0){
					cout << "Calculating x[" << i<< "]: " << total  << " Random walks generated" << endl;
				}
				x = sum1 / total;
				double __err = (sum2 - sum1 / total) / total / total ;
				times++;
				err = sqrt(__err) / x; 
			}
			cout << "Calculating x[" << i<< "]: " << total  << " Random walks generated" << endl;
			res[i] = x;
		}
		return res;
	}
	
	template <typename Type>
	std::vector<Type> backwards(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1){
		row = A.rows();col = A.cols();
		Matrix<Type> P(row +1, col +1);
		vector<double> pi(row, 0); //define the initial probability.
		for(int i = 0 ; i< row; i++){
			pi[i] = double(i + 1) / double(row);
		}
		pi[row - 1] = 1.0;
		double _sum1[row], _sum2[row];	
		int count[row];
		memset(count, 0, sizeof(count));
		memset(_sum1, 0, sizeof(_sum1));
		Matrix<Type> t = backwards_p(A, 0.4, P);
		size_t size = b.size();
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));	
		int total, _i, _times = 1;
		int en = 100;
		init();
		while(en--){
				while(step--){
					double _r = double(rand())/ RAND_MAX;
					_i= upper_bound(pi.begin(), pi.end(), _r) -pi.begin();
					double v = 1.0;
					int index = _i, next = 0;
					while(next != col){
						double r = double(rand())/RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) -t[index].begin();
						if(next == col) continue;	
						v = v * A[next][index] /  P[index][next];
						index = next;
					}
					for(int i = 0 ; i < row; i++){
						double _v = v;
						if(i == index){
							if(_i == 0){
								_v = _v / P[index][col] * (b[_i]/pi[_i]);
							}else{
								_v = _v /P[index][col] * (b[_i]/(pi[_i] - pi[_i-1]));
							}
							_sum1[i] += _v;
						}
					}
				}
				init();
				//step = 20000;
				total = step * _times;
				_times++;
				cout << "Backwards: " << total  << " Random walks generated" << endl;
		}
		cout <<total << endl;
//		total = en * step;
		for(int i = 0 ; i < row ;i++){
			res[i] = _sum1[i] /total;
		}
		return res;
	}
};
#endif
