/* Author:Bolong Zhang
*/

/* *************************************************************
 *                                 matrix.h
 *
 * Class template of matrix which is designed for basic linear algebra
 * operations such as:
 *              A1 + A2   A1 += A2
 *              A1 - A2   A1 -= A2
 *				A1 *A2
 *              A * x    x * A    A *= x
 *				
 **************************************************************/


#ifndef MATRIX_H
#define MATRIX_H


#include <vector>
#include<cassert>
#include <iostream>
using namespace std;

namespace matrix
{

    template <typename Type>
    class Matrix
    {
    public:
        // constructors and destructor
        Matrix();
		Matrix(int, int);
		Matrix( const Matrix<Type> &A );
		void unit();
        // accessors
        const vector<Type> &operator[]( int i ) const;
        vector<Type> &operator[]( int i );
        Type& operator()( int row, int column );
        const Type& operator()( int row, int column ) const;
        int rows() const;
        int cols() const ;
		void setScalar(const Type &x);
		void setRow(const size_t &, const vector<Type> &a);
		void addRow(const vector<Type> &a);
	private:
		typedef vector<Type> vec;
		vector<vec> _matrix;
		// row number, column number
        int	 row;
        int	 col;
        void init( int rows, int columns );
    };
	
	template<typename Type>
	void printmatrix(const Matrix<Type>& );

	template<typename Type, typename Class>
	Matrix<Type> operator+( const Matrix<Type>&, const Matrix<Class>& );
	template<typename Type, typename Class>
	Matrix<Type> operator-( const Matrix<Type>&, const Matrix<Class>& );
	template<typename Type, typename Class>
	vector<Type> operator*( const Matrix<Type>&, const vector<Class>& );
	template<typename Type, typename Class>
	vector<Type> operator*( const vector<Class> &, const Matrix<Type>& );
	template<typename Type, typename Class>
	Matrix<Type> operator*( const Matrix<Type>&, const Matrix<Class>& );
	template<typename Type, typename Class>
	Matrix<Type> operator*( const Matrix<Type>&, const Class&);
	template<typename Type, typename Class>
	Matrix<Type> operator*( const Class &, const Matrix<Type>&);
	template<typename Type, typename Class>
	Matrix<Type> operator/( const Class &, const Matrix<Type>&);
	template<typename Type, typename Class>
	Matrix<Type> operator/( const Matrix<Type>&, const Class&);
	
	#include "imatrix.h"
}

using namespace matrix;

#endif
// MATRIX_H
