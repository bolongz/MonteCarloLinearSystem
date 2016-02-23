

#include "matrix.h"
/*****************************************************************************
 *
 * Implementation for Matrix class.
 *
 *****************************************************************************/


/**
 * initialize
 */
template <typename Type>
void Matrix<Type>::init( int rows, int columns )
{
	row = rows;
	col = columns;
	_matrix.resize(row);
	for(int i = 0; i < row; i++){
		_matrix[i].resize(col);
		for(int j = 0; j < col; j++){
			_matrix[i][j] = 0;
		}
	}
}


/**
 * Insert the ith row
 */
template <typename Type>
void Matrix<Type>::setRow(const size_t &i,  const vector<Type> &a )
{ // if the length of a is large than col of the matrix
// we take the first col elements of a. Otherwise, the rest element 
// in the i row is zero.
	for(int j = 0 ; j < col;j++){	
		_matrix[i][j] = a[j];
	}

}


/**
 * set matrix by a scalar
 */
template <typename Type>
void Matrix<Type>::setScalar( const Type &x )
{
	for(int i = 0; i < row; i++){
		for(int j = 0; j < col; j++){
			_matrix[i][j] = x * _matrix[i][j];
		}
	}
}


/**
 * constructors and destructor
 */
template <typename Type>
Matrix<Type>::Matrix()
:row(0), col(0)
{
	_matrix.resize(0);
}

template <typename Type>
Matrix<Type>::Matrix( const Matrix<Type> &A )
{
	row = A.rows();
	col = A.cols();
	init(row, col);
	for(int i = 0 ; i < row ; i++){
		for(int j = 0 ; j < col ;j++){

			_matrix[i][j] = A[i][j];
		}
	}
}

/**
 * overload operator [] for 0-offset access
 */

template <typename Type>
inline const vector<Type> &Matrix<Type>::operator[]( int i )
{
	return _matrix[i];
}


/**
 * overload operator () for 1-offset access
 */
template <typename Type>
inline const Type& Matrix<Type>::operator()( int row, int column )
{
	return  _matrix[row][column];
}


template <typename Type>
const int Matrix<Type>::rows()
{
    return row;
}

template <typename Type>
const int Matrix<Type>::cols()
{
    return col;
}


/**
 * Overload the output stream function.
 */
template <typename Type>
ostream& operator<<( ostream &out, const Matrix<Type> &A )
{
	int rows = A.rows();
	int columns = A.cols();

	out << "size: " << rows << " by " << columns << "\n";
	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			out << A[i][j] << "\t";
		out << "\n";
	}

	return out;
}


/**
 * Overload the intput stream function.
 */
template <typename Type>
istream& operator>>( istream &in, Matrix<Type> &A )
{
	int rows, columns;
	in >> rows >> columns;

	if( !( rows == A.rows() && columns == A.cols() ) )
		A.resize( rows, columns );

	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			in >> A[i][j];

	return in;
}

/**
 * matrix-matrix addition
 */
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp(A1);		
	for(int i = 0 ; i <  A1.rows(); i++){
		for(int j = 0; j < A1.cols(); j++){
			tmp[i][j] = A1[i][j] + A2[i][j];
		}
	}
	return tmp;
}

/**
 * matrix-matrix subtraction
 */
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A1, const Matrix<Type> &A2 ){
	Matrix<Type> tmp(A1);		
	for(int i = 0 ; i <  A1.rows(); i++){
		for(int j = 0; j < A1.cols(); j++){
			tmp[i][j] = A1[i][j] + A2[i][j];
		}
	}
	return tmp;
}


/**
 * matrix-scaling multiplication
 */
template <typename Type>
inline Matrix<Type> operator*( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp = A.setScaler(x);
	
	return tmp;
}

template <typename Type>
inline Matrix<Type> operator*( const Type &x, const Matrix<Type> &A )
{
	return A * x;
}


/**
 * matrix-matrix multiplication
 */
template <typename Type>
Matrix<Type> operator*( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	assert( A1.cols() == A2.rows() );

	int rows = A1.rows();
	int columns = A2.cols();
	int K = A1.cols();
	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
		{
            tmp[i][j] = 0;
			for( int k=0; k<K; ++k )
			    tmp[i][j] += A1[i][k] * A2[k][j];
		}

	return tmp;
}


/**
 * matrix-vector multiplication
 */
template <typename Type>
vector<Type> operator*( const Matrix<Type> &A, const vector<Type> &b )
{
	assert( A.cols() == b.size() );

	int rows = A.rows();
	int columns = A.cols();
	vector<Type> tmp(rows);
	for( int i=0; i<rows; ++i )
	{
		Type sum = 0;
		for( int j=0; j<columns; ++j )
			sum += A[i][j] * b[j];
		tmp[i] = sum;
	}

	return tmp;
}

template <typename Type>
vector<Type> operator*( const vector<Type> &b, const Matrix<Type> &A )
{
	assert( A.rows() == b.size() );

	int cols = A.cols();
	int size = A.rows();
	vector<Type> tmp(cols);
	for( int i=0; i<cols; ++i )
	{
		Type sum = 0;
		for( int j=0; j<size; ++j )
			sum += A[i][j] * b[j];
		tmp[i] = sum;
	}

	return tmp;
}


/**
 * matrix-scalar division
 */
template <typename Type>
inline Matrix<Type> operator/( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp /= x;
}

template <typename Type>
Matrix<Type> operator/( const Type &x, const Matrix<Type> &A )
{
	int rows = A.rows();
	int clumns = A.cols();

	Matrix<Type> tmp( rows,clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = x / A[i][j];

	return tmp;
}
int main(){

	return 0;
}
