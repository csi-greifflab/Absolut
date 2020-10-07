// $Id: SquareMatrix.hh,v 1.2 2016/08/08 12:41:57 mmann Exp $
#ifndef BIU_SQUAREMATRIX_HH_
#define BIU_SQUAREMATRIX_HH_


#include "biu/Matrix.hh"
#include <biu/assertbiu.hh>

namespace biu
{

		//! Generic square matrix with several matrix operations.
		//!
		//! @author Martin Mann <mmann@@informatik.uni-freiburg.de>
	template <class T, size_t matrixDim>
	class SquareMatrix : public Matrix<T>
	{
	public:
			//! constructs an empty matrix of size 0
		SquareMatrix();
			//! copy construction including resizing
		SquareMatrix(const SquareMatrix& m);
			//! copy construction including resizing
		SquareMatrix(const Matrix<T>& m);
	
			//! construction and initialisation via twodimensional array
		SquareMatrix(const T m[matrixDim][matrixDim]);
	
	
		~SquareMatrix();
		
		static SquareMatrix<T,matrixDim> createSquareMatrix(const T  data[matrixDim][matrixDim]);
		
	};
	
	
	/* IMPLEMENTIERUNG */
	
	//construction
	template <class T, size_t matrixDim> inline 
	SquareMatrix<T,matrixDim>::SquareMatrix() :
		Matrix<T>(matrixDim, matrixDim)
	{
	}
	
	//! construction from array
	template <class T, size_t matrixDim> inline 
	SquareMatrix<T,matrixDim>::SquareMatrix(const T m[matrixDim][matrixDim])
	{
		*this = createSquareMatrix(m);
	}
	
	//! copy constructor
	template <class T, size_t matrixDim> inline
	SquareMatrix<T,matrixDim>::SquareMatrix(const SquareMatrix& m) :
		Matrix<T>(m)
	{
	}
	
	template <class T, size_t matrixDim> inline
	SquareMatrix<T,matrixDim>::SquareMatrix(const Matrix<T>& m) :
		Matrix<T>(m)
	{
		// m has to be a square matrix of the correct size
		assertbiu(m.numRows() == matrixDim && m.numColumns() == matrixDim, "matrices differ in dimension");	
	}
	
	
	
	// destruction
	template <class T, size_t matrixDim> inline
	SquareMatrix<T,matrixDim>::~SquareMatrix ()
	{
	}
	
	// construct and fill
	#include <iostream>
	template <class T, size_t matrixDim> inline 
	SquareMatrix<T,matrixDim>
	SquareMatrix<T,matrixDim>::createSquareMatrix(const T  data[matrixDim][matrixDim]) 
	{
		SquareMatrix<T,matrixDim> matrix;
		for (size_t row=0; row < matrixDim; row++)
		{
			for (size_t col=0; col<matrixDim; col++) {
				matrix[row][col] = data[row][col];
			}
		}
		return matrix;
	}
	
	

} // namespace biu

#endif /*SQUAREMATRIX_HH_*/
