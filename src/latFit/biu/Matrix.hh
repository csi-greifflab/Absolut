// $Id: Matrix.hh,v 1.2 2016/08/08 12:41:59 mmann Exp $
#ifndef BIU_MATRIX_HH
#define BIU_MATRIX_HH

#include <iostream>
#include <vector>

namespace biu {


		/** A generic matrix that allows matrix operations.
		 */
	template <class T> class Matrix {
	protected:
			//! matrix dimensions = number of rows
		size_t rows;
			//! matrix dimensions = number of columns
		size_t cols;
			//! the data store
		T **v;
	public:
			
			//! creates an empty matrix of dimensions (0,0) 
		Matrix();
			//! creates an unitialised matrix of given dimensions 
			//! @param r the row number of the matrix
			//! @param c the column number of the matrix
		Matrix(const size_t r, const size_t c);
			//! creates a matrix and fills all elements
			//! @param r the row number of the matrix
			//! @param c the column number of the matrix
			//! @param val the default value for the matrix entries
		Matrix(const size_t r, const size_t c, const T &val);
			//! creates a matrix that contains a copy of the given matrix 
			//! @param r the row number of the matrix
			//! @param c the column number of the matrix
			//! @param mat an array containing the matrix in the form 
			//!          element(row,col) = mat[row*c+col] 
			//!          with row = (0,r] && col = (0,c]
		Matrix(const size_t r, const size_t c, const T mat[]);
			//! copy constructor including matrix resizing
			//! @param mat the matrix this matrix should be a copy of
		Matrix(const Matrix<T> &mat);
			//! assignment including resizing
			//! @param mat the matrix this matrix should be a copy of
			//! @return *this 
		Matrix& operator = (const Matrix<T> &mat);
			//! all elements = val
			//! @param val the new value for all matrix entries
			//! @return *this 
		Matrix& operator = (const T &val);
			//! matrix + matrix
			//! row and column number have to be equal
			//! @param mat the matrix to add
			//! @return a new matrix representing (this + mat)
		Matrix<T> operator + (const Matrix<T> &mat);
			//! matrix + constant
			//! @param c the constant to add to all elements
			//! @return a new matrix representing (this + c)
		Matrix<T> operator + (const T &c);
			//! matrix - matrix
			//! row and column number have to be equal
			//! @param mat the matrix to substract
			//! @return a new matrix representing (this - mat)
		Matrix<T> operator - (const Matrix<T> &mat);
			//! matrix - constant
			//! @param c the constant to substract from all elements
			//! @return a new matrix representing (this - c)
		Matrix<T> operator - (const T &c);
			//! matrix * matrix
			//! row number have to equal column number and vice versa
			//! @param mat the matrix to multiply with
			//! @return a new matrix representing (this * mat)
		Matrix<T> operator * (const Matrix<T> &mat) const;
			//! matrix * constant
			//! @param c the constant to multiply all elements with
			//! @return a new matrix representing (this * c)
		Matrix<T> operator * (const T &c) const;
			//! special case of *
			//! column number has to equal vec size
			//! @param vec the column vector to multiply with
			//! @return a new column vector representing (this * vec)
		std::vector<T> operator * (const std::vector<T> &vec) const;
			//! matrix *= matrix
			//! row number have to equal column number and vice versa
			//! @param mat the matrix to multiply with
			//! @return *this
		Matrix<T>& operator *= (const Matrix<T> &mat);
			//! matrix *= constant
			//! @param c the constant to multiply all elements with
			//! @return *this
		Matrix<T>& operator *= (const T& c);
			//! test for equality
			//! @param mat the matrix to compare with
			//! @return true if dimensions and all elements are equal, 
			//!         false otherwise 
		inline bool operator == (const Matrix<T> &mat) const;
			//! test if all elements == val
			//! @param val the value to compare with
			//! @return true if all elements are equal, 
			//!         false otherwise 
		inline bool operator == (const T &val) const;
			//! test for inequality
			//! @param mat the matrix to compare with
			//! @return true if dimensions or at least one element are is not equal, 
			//!         false otherwise 
		inline bool operator != (const Matrix<T> &mat) const;
			//! test if at least one element doesnt equal val
			//! @param val the value to compare with
			//! @return true if at least one element are is not equal, 
			//!         false otherwise 
		inline bool operator != (const T &val) const;
			//! access to row 
			//! @param row the row to access
			//! @return array access to the row
		inline T* const operator[](const size_t row);
			//! constant access to row i
			//! @param row the row to access
			//! @return constant array access to the row
		inline const T* const operator[](const size_t row) const;
			//! constant access to element (r,c)
			//! @param r the row to access
			//! @param c the column to access
			//! @return the matrix element at position (r,c) 
		inline T at(const size_t r, const size_t c) const;
			//! direct access to element (r,c)
			//! @param r the row to access
			//! @param c the column to access
			//! @return the matrix element at position (r,c) 
		inline T& at(const size_t r, const size_t c);
			//! number of rows
			//! @return row number
		inline size_t numRows() const;
			//! number of columns
			//! @return column number
		inline size_t numColumns() const;
			//! access to a matrix column
			//! @param col the column to access
			//! @return the column vector
		std::vector<T> columnVec(const size_t col) const;
			//! Resizes the matrix. If the the dimensions are increased, the 
			//! new fields are not initialised.
			//! @param row the new row number
			//! @param col the new column number
		void resize(const size_t row, const size_t col);
			//! Resizes the matrix. If the the dimensions are increased, the 
			//! new fields are filled with the given default value
			//! @param row the new row number
			//! @param col the new column number
			//! @param defVal the default value for new matrix elements
		void resize(const size_t row, const size_t col, const T& defVal);
			//! destruction and freeing memory 
		~Matrix();
		
	};

} // namespace biu

  //! Reads (m.rows*m.cols) matrix elements from stream assuming that they are 
  //! whitespace separated.
  //! @param in the output stream to write to
  //! @param m the matrix to print
  //! @return the changed output stream out
template <class T> inline
std::istream&
operator>>(std::istream& in,  biu::Matrix<T>& m);

  //! Prints the matrix elements blank (' ') separated and linewise to stream.
  //! @param out the output stream to write to
  //! @param m the matrix to print
  //! @return the changed output stream out
template <class T> inline
std::ostream&
operator << (std::ostream& out, const biu::Matrix<T>& m);

#include "biu/Matrix.icc"
	
#endif // define MATRIX_HH
