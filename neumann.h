#ifndef _NEUMANN_H_
#define _NEUEMAN_H_

#include "matrix.h"
#include "pdf.h"
#include <time.h>
#include<stdlib.h>
#include<algorithm>
#include <thread>

template <typename Type>
class Neumann{
private:
	size_t row;
	size_t col;
	double err, sum1, sum2, x, err_w;
	int next, step, times;
	int hops;	
	Matrix<Type> A_ForThreading;
	std::vector<Type> b_ForThreading;
	Matrix<Type> P_ForThreading, t_ForThreading;
	std::vector<Type> res_ForThreading;
	std::vector<int> Thr_boun;
	double _err_Glob;
public:
	void init() {
		err = 10000.0;
		sum1 = 0.0; sum2 = 0.0; x = 0.0;
		next = 0; step = 20000; times = 1;
		err_w = 1e-6;
		hops = 0;
	}
	//template <typename Type>
	std::vector<Type> absorbing(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1) {
		row = A.rows(); col = A.cols();
		Matrix<Type> P(row + 1, col + 1);
		Matrix<Type> t = Absorbing(A, 0.2, P);
		//printmatrix(A);
//		printmatrix(P);
//		printmatrix(t);
		size_t size = b.size();
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));
		int total = 0, _total = 0, _hops = 0;
		for (int i = 0; i < size; i++) {
			init();
			total = 0;
			cout << "Calculating x[" << i << "]... "  << endl;
			while (err > _err) {
				int cc = step;
				while (step--) {
					double v = 1.0;
					int index = i, next = 0;
					while (next != col) {
						double r = double(rand()) / RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) - t[index].begin();
						if (next == col) continue;
						if(abs(P[index][next]) > 1e-6 ){
							v = v * A[index][next] / P[index][next];
						}else{
							v = 0;
						}
						hops++;_hops++;
						index = next;
					}
					v = v * b[index] / P[index][col];
					sum1 += v;
					sum2 += v * v;
				}

				step = 1;
				total = total + cc;
				if (total % 200000 == 0) {
					cout  << total << " Random walks generated" << endl;
				}
				x = sum1 / total;
				double __err = (sum2 - sum1 / total) / total / total;
				times++;
				err = sqrt(__err) / x;
			}
			cout << endl;
			cout << "Total random walks: " << total <<endl;
			cout <<  "Average hops: " << hops/total <<  endl;
			res[i] = x;
			_total += total;
		}
		
		cout << endl;
		cout << "Avegage total random walks: " << _total <<endl;
		cout <<  "Average hops: " << _hops/_total <<  endl;
		
		return res;
	}
	//template <typename Type>
	std::vector<Type> nonabsorbing(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1) {
		row = A.rows(); col = A.cols();
		Matrix<Type> P(row, col);
		Matrix<Type> t = nonAbsorbing(A, P);
		size_t size = b.size();
		//printmatrix(P);
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));
		int total = 0, _total = 0, _hops = 0;;
		for (int i = 0; i < size; i++) {
			init();
			total = 0;
			cout << "Calculating x[" << i << "]... " << endl;
			while (err > _err) {
				int cc = step;
				while (step--) {
					double v = 0.0, w = 1.0;
					int index = i, next = 0;
					while (abs(w) > err_w) {
						double r = double(rand()) / RAND_MAX;
						next = upper_bound(t[index].begin(), t[index].end(), r) - t[index].begin();
						if(P[index][next] > 1e-6){
							w = w * A[index][next] / P[index][next];
						}else{
							w = 0;
						}
						hops++;_hops++;
						v = v + w * b[next];
						index = next;
					}

					v = v + b[i];
					sum1 += v;
					sum2 += v * v;
				}
				step = 1;;
				total = total + cc;
				if (total % 200000 == 0) {
					//cout << "Calculating x[" << i << "]: " << total << " Random walks generated" << endl;
					cout  << total << " Random walks generated" << endl;
				}
				x = sum1 / total;
				double __err = (sum2 - sum1 / total) / total / total;
				times++;
				err = sqrt(__err) / x;
				if(total >= 5000000) break;
			}
			cout << endl;
			cout << "Total random walks: " << total <<endl;
			cout <<  "Average hops: " << hops/total <<  endl;
			_total += total;
			res[i] = x;
		}
		cout << "Avegage total random walks: " << _total <<endl;
		cout <<  "Average hops: " << _hops/_total <<  endl;
		return res;
	}
	
	std::vector<Type> backwards(const Matrix<Type> &A, const std::vector<Type> &b, int nwalks){
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
		Matrix<Type> t = backwards_p(A, 0.2, P);
		size_t size = b.size();
		std::vector<Type> res(size);
		srand((unsigned)time(NULL));	
		int total, _i, _times = 1;
		int en = 100;
		init();
		step = 20000;
		en = nwalks;
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
						if(P[index][next] > 1e-6){
							v = v * A[next][index] /  P[index][next];
						}else{
							v = 0;
						}
						hops++;
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
				int cc = hops;
				init();
				hops = cc;
				step = 20000;
				total = step * _times;
				_times++;
				cout << "Backwards :" << total << " Random number generated" << endl;
		}
//		total = en * step;
		for(int i = 0 ; i < row ;i++){
			res[i] = _sum1[i] /total;
		}
				cout << endl;
				cout << "Total random walks: " << total <<endl;
				cout <<  "Average hops: " << hops/total <<  endl;
		return res;
	}
	//Threaded Absorbing
	void abs_thread(int arg)
	{
		int startingBoun = 0;
		double err_trd, _err_trd = _err_Glob, v, r, __err_trd;
		double sum1_trd, sum2_trd, x_trd;
		int step_trd, times_trd, index, next_trd, total;
		if (arg > 0)
			startingBoun = Thr_boun[arg - 1];
		for (int i = startingBoun; i < Thr_boun[arg]; i++)
		{
			err_trd = 10000.0;
			sum1_trd = 0.0; sum2_trd = 0.0; x_trd = 0.0;
			step_trd = 100; times_trd = 1;

			while (err_trd > _err_trd)
			{
				while (step_trd--)
				{
					v = 1.0;
					index = i;
					next_trd = 0;
					while (next_trd != col)
					{
						r = double(rand()) / RAND_MAX;
						next_trd = upper_bound(t_ForThreading[index].begin(), t_ForThreading[index].end(), r) - t_ForThreading[index].begin();
						if (next_trd == col) continue;
						if (P_ForThreading[index][next_trd] != 0)
							v = v * A_ForThreading[index][next_trd] / P_ForThreading[index][next_trd];
						else
							v = 0;
						index = next_trd;
					}
					v = v * b_ForThreading[index] / P_ForThreading[index][col];
					sum1_trd += v;
					sum2_trd += v * v;
				}
				step_trd = 100;
				total = step_trd* times_trd;
				if (total % 200000 == 0) {
					cout << "Calculating x[" << i << "]: " << total << " Random walks generated" << endl;
				}
				x_trd = sum1_trd / total;
				__err_trd = (sum2_trd - sum1_trd / total) / total / total;
				times_trd++;
				err_trd = sqrt(__err_trd) / x_trd;
			}
			cout << "Calculating x[" << i << "]: " << total << " Random walks generated" << endl;
			res_ForThreading[i] = x_trd;
		}
	}

	//template <typename Type>
	std::vector<Type> absorbing_UsingThreads(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1, int thrds = 3) {
		std::thread threads[100];
		_err_Glob = _err;
		row = A.rows(); col = A.cols();
		A_ForThreading = A;
		b_ForThreading = b;
		P_ForThreading = Matrix<Type>(row + 1, col + 1);
		t_ForThreading = Absorbing(A, 0.2, P_ForThreading);
		size_t size = b.size();
		res_ForThreading = std::vector<Type>(size);
		srand((unsigned)time(NULL));
		Thr_boun = std::vector<int>();
		int takenJobs = 0;
		int thrdsNum = thrds;
		if (thrdsNum > size)
			thrdsNum = size;
		for (int i = 0; i < thrdsNum; i++) {
			takenJobs = (int)ceil((size - takenJobs) / (thrdsNum - i)) + takenJobs;
			Thr_boun.push_back(takenJobs);
			threads[i] = std::thread(&Neumann<Type>::abs_thread, this, i);
		}
		for (int i = 0; i < thrdsNum; i++) {
			threads[i].join();
		}
		return res_ForThreading;
	}

	//Threaded NonAbsorbing
	void nonabs_thread(int arg)
	{
		int startingBoun = 0;
		double err_trd, _err_trd = _err_Glob, v, r, w, err_w_trd, __err_trd;
		double sum1_trd, sum2_trd, x_trd;
		int step_trd, times_trd, index, next_trd, total;
		if (arg > 0)
			startingBoun = Thr_boun[arg - 1];
		for (int i = startingBoun; i < Thr_boun[arg]; i++)
		{
			err_trd = 10000.0; err_w_trd = 1e-6;
			sum1_trd = 0.0; sum2_trd = 0.0; x_trd = 0.0;
			step_trd = 100; times_trd = 1;

			while (err_trd > _err_trd)
			{
				while (step_trd--)
				{
					v = 0.0;
					w = 1.0;
					index = i;
					next_trd = 0;
					while (abs(w) > err_w_trd) {
						r = double(rand()) / RAND_MAX;
						next_trd = upper_bound(t_ForThreading[index].begin(), t_ForThreading[index].end(), r) - t_ForThreading[index].begin();
						if (P_ForThreading[index][next_trd] != 0)
							w = w * A_ForThreading[index][next_trd] / P_ForThreading[index][next_trd];
						else
							w = 0;
						v = v + w * b_ForThreading[next_trd];
						index = next_trd;
					}
					v = v + b_ForThreading[i];
					sum1_trd += v;
					sum2_trd += v * v;
				}
				step_trd = 100;
				total = step_trd* times_trd;
				if (total % 200000 == 0) {
					cout << "Calculating x[" << i << "]: " << total << " Random walks generated" << endl;
				}
				x_trd = sum1_trd / total;
				__err_trd = (sum2_trd - sum1_trd / total) / total / total;
				times_trd++;
				err_trd = sqrt(__err_trd) / x_trd;
			}
			cout << "Calculating x[" << i << "]: " << total << " Random walks generated" << endl;
			res_ForThreading[i] = x_trd;
		}
	}

	//template <typename Type>
	std::vector<Type> nonabsorbing_UsingThreads(const Matrix<Type> &A, const std::vector<Type> &b, double _err = 0.1, int thrds = 3) {
		std::thread threads[100];
		_err_Glob = _err;
		row = A.rows(); col = A.cols();
		A_ForThreading = A;
		b_ForThreading = b;
		P_ForThreading = Matrix<Type>(row, col);
		t_ForThreading = nonAbsorbing(A, P_ForThreading);
		size_t size = b.size();
		res_ForThreading = std::vector<Type>(size);
		srand((unsigned)time(NULL));
		Thr_boun = std::vector<int>();
		int takenJobs = 0;
		int thrdsNum = thrds;
		if (thrdsNum > size)
			thrdsNum = size;
		for (int i = 0; i < thrdsNum; i++) {
			takenJobs = (int)ceil((size - takenJobs) / (thrdsNum - i)) + takenJobs;
			Thr_boun.push_back(takenJobs);
			threads[i] = std::thread(&Neumann<Type>::nonabs_thread, this, i);
		}
		for (int i = 0; i < thrdsNum; i++) {
			threads[i].join();
		}
		return res_ForThreading;
	}
};
#endif
