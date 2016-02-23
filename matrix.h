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

        // accessors
        const vector<Type> &operator[]( int i );
        const Type& operator()( int row, int column );
        const int rows();
        const int cols();
		void setScalar(const Type &x);

		
		typedef vector<Type> vec;
		vector<vec> _matrix;
		// row number, column number and total number
        int	 row;
        int	 col;
        void init( int rows, int columns );
		void setRow(const size_t &, const vector<Type> &a);
    };
    // class Matrix
    template<typename Type>
    ostream& operator<<( ostream&, const Matrix<Type>& );
    template<typename Type>
    istream& operator>>( istream&, Matrix<Type>& );

    template<typename Type>
    Matrix<Type> operator+( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator-( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator*( const Matrix<Type>&, const vector<Type>& );
    template<typename Type>
    Matrix<Type> operator*( const vector<Type> &, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator*( const Matrix<Type>&, const Matrix<Type>& );
	template<typename Type>
	Matrix<Type> operator*( const Matrix<Type>&, const Type&);
	template<typename Type>
	Matrix<Type> operator*( const Type &, const Matrix<Type>&);
	template<typename Type>
	Matrix<Type> operator/( const Type &, const Matrix<Type>&);
	template<typename Type>
	Matrix<Type> operator*( const Matrix<Type>&, const Type&);

}

using namespace matrix;

#endif
// MATRIX_H
