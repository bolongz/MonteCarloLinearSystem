#ifndef _MARKOV_CHAIN_H_
/* Author:Omar Azhar
*/

/* *************************************************************
*                                 markovChine.h
*
* Class template for implementing Markov Chine
*
**************************************************************/
#define _MARKOV_CHAIN_H_

#include "matrix.h"
#include <stdlib.h>
#include <time.h>

/**
* Provide a radom state 
*/
template<typename Type>
int init_MC(const Matrix<Type> A, int row_index)
{
	int i;
	double source = A(0, 0);
	double randm = ((double)rand()) / RAND_MAX;
	for (i = 0; randm > source & (i < A.cols()); i++)
	{
		source = source + A(row_index, i+1);
	}
	return i;
}


/**
* Provide a sequence
*/
template<typename Type>
vector<int> provideSequence(const Matrix<Type> Init_prob_vector, const Matrix<Type> trans_matrix, int sequSize)
{
	int i, index;
	vector<int> result;
	result.resize(sequSize);
	int init = init_MC(Init_prob_vector, 0);
	result[0] = index = init;
	for (i = 1; i<sequSize; i++)
	{
		index = init_MC(trans_matrix, index);
		result[i] = index;
	}
	return result;
}

#endif

