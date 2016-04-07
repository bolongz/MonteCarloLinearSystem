#ifndef _NEUMANN_H_
#define _NEUEMAN_H_

#include "matrix.h"
#include "pdf.h"
#include <time.h>
#include<stdlib.h>
#include<algorithm>
#include <pthread.h>
#include <unistd.h>

template <typename Type>
class Neumann{
private:
	size_t row;
	size_t col;
	double err, sum1, sum2, x, err_w;
	int next, step, times;
	
	int threadNum[100];
	pthread_t threadID[100];
	const Matrix<Type> A_ForThreading;
	const std::vector<Type> b_ForThreading;
	Matrix<Type> P_ForThreading, t_ForThreading;
	std::vector<Type> res_ForThreading;
public:
	void init() {
		err = 10000.0;
		sum1 = 0.0; sum2 = 0.0; x = 0.0;
		next = 0; step = 100; times = 1;
		err_w = 1e-6;
	}
	//template <typename Type>
	std::vector<Type> absorbing(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1) {
		row = A.rows(); col = A.cols();
		Matrix<Type> P(row + 1, col + 1);
		Matrix<Type> t = Absorbing(A, 0.2, P);
		size_t size = b.size();
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));
		int total;
		for (int i = 0; i < size; i++) {
			init();
			while (err > _err) {
				while (step--) {
					double v = 1.0;
					int index = i, next = 0;
					while (next != col) {
						double r = double(rand()) / RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) - t[index].begin();
						if (next == col) continue;
						v = v * A[index][next] / P[index][next];
						index = next;
					}
					v = v * b[index] / P[index][col];
					sum1 += v;
					sum2 += v * v;
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
	//template <typename Type>
	std::vector<Type> nonabsorbing(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1) {
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
					sum2 += v * v;
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

	//Threaded Absorbing
	void *abs_thread(void *arg)
	{
		double err_trd = 10000.0, _err_trd = 0.1, v = 1.0, r;
		double sum1_trd = 0.0, sum2_trd = 0.0, x_trd = 0.0;
		int step_trd = 100, times_trd = 1, index = *(int *)arg, next_trd, total;

		while (err_trd > _err_trd)
		{
			while (step_trd--)
			{
				v = 1.0;
				index = *(int *)arg;
				next_trd = 0;
				while (next_trd != A_ForThreading.cols)
				{
					r = double(rand()) / RAND_MAX;
					next_trd = upper_bound(t_ForThreading[index].begin(), t_ForThreading[index].end(), r) - t_ForThreading[index].begin();
					if (next_trd == A_ForThreading.cols) continue;
					v = v * A_ForThreading[index][next_trd] / P_ForThreading[index][next_trd];
					index = next_trd;
				}
				v = v * b_ForThreading[index] / P_ForThreading[index][A_ForThreading.cols];
				sum1_trd += v;
				sum2_trd += v * v;
			}
			step_trd = 100;
			total = step_trd* times_trd;
			if (total % 200000 == 0) {
				cout << "Calculating x[" << *(int *)arg << "]: " << total << " Random walks generated" << endl;
			}
			x_trd = sum1_trd / total;
			double __err_trd = (sum2_trd - sum1_trd / total) / total / total;
			times_trd++;
			err_trd = sqrt(__err_trd) / x_trd;
		}
		cout << "Calculating x[" << *(int *)arg << "]: " << total << " Random walks generated" << endl;
		res_ForThreading[*(int *)arg] = x_trd;
	}

	std::vector<Type> absorbing_UsingThreads(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1) {
		row = A.rows(); col = A.cols();
		A_ForThreading = A;
		b_ForThreading = b;
		P_ForThreading = new Matrix<Type>(row + 1, col + 1);
		t_ForThreading = Absorbing(A, 0.2, P_ForThreading);
		size_t size = b.size();
		res_ForThreading = new std::vector<Type>(size);
		srand((unsigned)time(NULL));
		for (int i = 0; i < size; i++) {
			threadNum[i] = i;
			pthread_create(&threadID[i], NULL, &abs_thread, &threadNum[i]);
		}
		for (int i = 0; i < size; i++) {
			pthread_join(threadID[i], NULL);
		}
		return res_ForThreading;
	}

	//Threaded NonAbsorbing
	void *nonabs_thread(void *arg)
	{
		double err_trd = 10000.0, _err_trd = 0.1, v = 0.0, r, w = 1.0, err_w_trd = 1e-6;;
		double sum1_trd = 0.0, sum2_trd = 0.0, x_trd = 0.0;
		int step_trd = 100, times_trd = 1, index = *(int *)arg, next_trd, total;

		while (err_trd > _err_trd)
		{
			while (step_trd--)
			{
				v = 0.0;
				w = 1.0;
				index = *(int *)arg;
				next_trd = 0;
				while (abs(w) > err_w_trd) {
					r = double(rand()) / RAND_MAX;
					next_trd = upper_bound(t_ForThreading[index].begin(), t_ForThreading[index].end(), r) - t_ForThreading[index].begin();
					w = w * A_ForThreading[index][next_trd] / P_ForThreading[index][next_trd];
					v = v + w * b_ForThreading[next_trd];
					index = next_trd;
				}
				v = v + b_ForThreading[*(int *)arg];
				sum1_trd += v;
				sum2_trd += v * v;
			}
			step_trd = 100;
			total = step_trd* times_trd;
			if (total % 200000 == 0) {
				cout << "Calculating x[" << *(int *)arg << "]: " << total << " Random walks generated" << endl;
			}
			x_trd = sum1_trd / total;
			double __err_trd = (sum2_trd - sum1_trd / total) / total / total;
			times_trd++;
			err_trd = sqrt(_err_trd) / x_trd;
		}
		res_ForThreading[*(int *)arg] = x_trd;
	}

	std::vector<Type> nonabsorbing_UsingThreads(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1) {
		row = A.rows(); col = A.cols();
		A_ForThreading = A;
		b_ForThreading = b;
		P_ForThreading = new Matrix<Type>(row, col);
		t_ForThreading = nonAbsorbing(A, P_ForThreading);
		size_t size = b.size();
		res_ForThreading = new std::vector<Type>(size);
		srand((unsigned)time(NULL));
		for (int i = 0; i < size; i++) {
			threadNum[i] = i;
			pthread_create(&threadID[i], NULL, &nonabs_thread, &threadNum[i]);
		}
		for (int i = 0; i < size; i++) {
			pthread_join(threadID[i], NULL);
		}
		return res_ForThreading;
	}
};
#endif
