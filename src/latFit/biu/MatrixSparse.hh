// $Id: MatrixSparse.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_MATRIX_SPARSE_HH
#define BIU_MATRIX_SPARSE_HH

#include "biu/HashMap.hh"
// include best available hash or map
#if HAVE_UNORDERED_MAP == 1
	#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP == 1
	#include <tr1/unordered_map>
#elif HAVE_GNU_HASH_MAP == 1
	#include <ext/hash_map>
#else
	#include <map>
#endif

#include <vector>
#include <iostream>


namespace biu {


		/** 
		 * A generic sparse matrix that allows matrix operations.
		 * 
		 * The data is stored in a map of maps, i.e. for each row a map of
		 * non-empty entries is maintained.
		 */
	template <class T> class MatrixSparseR {
	public:
		// set typedef for best available hash or map
		#if HAVE_UNORDERED_MAP == 1
			typedef typename std::unordered_map< size_t, T > EntryMap;
			typedef typename std::unordered_map< size_t, EntryMap > RowMap;
		#elif HAVE_TR1_UNORDERED_MAP == 1
			typedef typename std::tr1::unordered_map< size_t, T > EntryMap;
			typedef typename std::tr1::unordered_map< size_t, EntryMap > RowMap;
		#elif HAVE_GNU_HASH_MAP == 1
			typedef typename __gnu_cxx::hash_map< size_t, T > EntryMap;
			typedef typename __gnu_cxx::hash_map< size_t, EntryMap > RowMap;
		#else
			typedef typename std::map< size_t, T > EntryMap;
			typedef typename std::map< size_t, EntryMap > RowMap;
		#endif
		
	protected:
			//! matrix dimensions = number of rows
		size_t rows;
			//! matrix dimensions = number of columns
		size_t cols;
			//! the data store
		RowMap row2col2entry;
			//! the default values used to create new entries
		T defaultVal;
		
	public:
			
			//! creates an empty matrix of dimensions (0,0)
			//! @param defVal the default value to use for new entries
		explicit 
		MatrixSparseR(const T &defVal);
			//! creates a matrix of given dimension
			//! @param r the row number of the matrix
			//! @param c the column number of the matrix
			//! @param val the default value for new matrix entries
		MatrixSparseR(const size_t r, const size_t c, const T &val);
			//! copy constructor including matrix resizing
			//! @param mat the matrix this matrix should be a copy of
		MatrixSparseR(const MatrixSparseR<T> &mat);
			//! assignment including resizing
			//! @param mat the matrix this matrix should be a copy of
			//! @return *this 
		MatrixSparseR& operator = (const MatrixSparseR<T> &mat);
			//! special case of *
			//! column number has to equal vec size
			//! @param vec the column vector to multiply with
			//! @return a new column vector representing (this * vec)
		std::vector<T> operator * (const std::vector<T> &vec) const;
			//! test for equality
			//! @param mat the matrix to compare with
			//! @return true if dimensions and all elements are equal, 
			//!         false otherwise 
		inline bool operator == (const MatrixSparseR<T> &mat) const;
			//! test for inequality
			//! @param mat the matrix to compare with
			//! @return true if dimensions or at least one element are is not equal, 
			//!         false otherwise 
		inline bool operator != (const MatrixSparseR<T> &mat) const;
			//! constant access to element (r,c)
			//! @param r the row to access
			//! @param c the column to access
			//! @return the matrix element at position (r,c) 
		inline T at(const size_t r, const size_t c) const;
			//! explicit constant access to element (r,c)
			//! @param r the row to access
			//! @param c the column to access
			//! @return the matrix element at position (r,c) 
		inline T atConst(const size_t r, const size_t c) const;
			//! direct access to element (r,c)
			//! @param r the row to access
			//! @param c the column to access
			//! @return the matrix element at position (r,c) 
		inline T& at(const size_t r, const size_t c);
			//! checks if element (r,c) is non-empty
			//! @param r the row to access
			//! @param c the column to access
			//! @return wether or not (r,c) exists 
		inline bool exists(const size_t r, const size_t c) const;
			//! number of rows
			//! @return row number
		inline size_t numRows() const;
			//! number of columns
			//! @return column number
		inline size_t numColumns() const;
			//! access to a matrix row
			//! @param row the row to access
			//! @return the row vector
		std::vector<T> rowVec(const size_t row) const;
			//! access to a matrix column
			//! @param col the column to access
			//! @return the column vector
		std::vector<T> columnVec(const size_t col) const;
			//! access to the non-empty matrix row entries
			//! @param row the row to access
			//! @return the non-empty row entries
		EntryMap rowValues(const size_t row) const;
			//! access to the non-empty matrix column entries
			//! @param col the column to access
			//! @return the non-empty column entries
		EntryMap columnValues(const size_t col) const;
			//! Resizes the matrix. If the the dimensions are increased, the 
			//! new fields will be initialised with current default value
			//! @param row the new row number
			//! @param col the new column number
		void resize(const size_t row, const size_t col);
			//! Resizes the matrix. If the the dimensions are increased, 
			//! sets the default value of new fields to the given value
			//! @param row the new row number
			//! @param col the new column number
			//! @param defVal the default value for new matrix elements
		void resize(const size_t row, const size_t col, const T& defVal);
			//! Sets the default value for new fields
			//! @param defVal the default value for new matrix elements
		void setDefaultValue( const T& defVal );
			//! Access to the default value for new fields
			//! @return the default value for new matrix elements
		const T& getDefaultValue( void ) const;
			//! destruction and freeing memory 
		~MatrixSparseR();
	};

} // namespace biu

namespace biu {


		/** 
		 * A generic sparse matrix that allows matrix operations.
		 * 
		 * The data is stored in a map of maps, i.e. for each column a map of
		 * non-empty entries is maintained.
		 */
	template <class T> class MatrixSparseC {
	public:
		// set typedef for best available hash or map
		#if HAVE_UNORDERED_MAP == 1
			typedef typename std::unordered_map< size_t, T > EntryMap;
			typedef typename std::unordered_map< size_t, EntryMap > ColMap;
		#elif HAVE_TR1_UNORDERED_MAP == 1
			typedef typename std::tr1::unordered_map< size_t, T > EntryMap;
			typedef typename std::tr1::unordered_map< size_t, EntryMap > ColMap;
		#elif HAVE_GNU_HASH_MAP == 1
			typedef typename __gnu_cxx::hash_map< size_t, T > EntryMap;
			typedef typename __gnu_cxx::hash_map< size_t, EntryMap > ColMap;
		#else
			typedef typename std::map< size_t, T > EntryMap;
			typedef typename std::map< size_t, EntryMap > ColMap;
		#endif
		
	protected:
			//! matrix dimensions = number of rows
		size_t rows;
			//! matrix dimensions = number of columns
		size_t cols;
			//! the data store
		ColMap col2row2entry;
			//! the default values used to create new entries
		T defaultVal;
		
	public:
			
			//! creates an empty matrix of dimensions (0,0)
			//! @param defVal the default value to use for new entries
		explicit 
		MatrixSparseC(const T &defVal);
			//! creates a matrix of given dimension
			//! @param r the row number of the matrix
			//! @param c the column number of the matrix
			//! @param val the default value for new matrix entries
		MatrixSparseC(const size_t r, const size_t c, const T &val);
			//! copy constructor including matrix resizing
			//! @param mat the matrix this matrix should be a copy of
		MatrixSparseC(const MatrixSparseC<T> &mat);
			//! assignment including resizing
			//! @param mat the matrix this matrix should be a copy of
			//! @return *this 
		MatrixSparseC& operator = (const MatrixSparseC<T> &mat);
			//! special case of *
			//! column number has to equal vec size
			//! @param vec the column vector to multiply with
			//! @return a new column vector representing (this * vec)
		std::vector<T> operator * (const std::vector<T> &vec) const;
			//! test for equality
			//! @param mat the matrix to compare with
			//! @return true if dimensions and all elements are equal, 
			//!         false otherwise 
		inline bool operator == (const MatrixSparseC<T> &mat) const;
			//! test for inequality
			//! @param mat the matrix to compare with
			//! @return true if dimensions or at least one element are is not equal, 
			//!         false otherwise 
		inline bool operator != (const MatrixSparseC<T> &mat) const;
			//! constant access to element (r,c)
			//! @param r the row to access
			//! @param c the column to access
			//! @return the matrix element at position (r,c) 
		inline T at(const size_t r, const size_t c) const;
			//! explicit constant access to element (r,c)
			//! @param r the row to access
			//! @param c the column to access
			//! @return the matrix element at position (r,c) 
		inline T atConst(const size_t r, const size_t c) const;
			//! direct access to element (r,c)
			//! @param r the row to access
			//! @param c the column to access
			//! @return the matrix element at position (r,c) 
		inline T& at(const size_t r, const size_t c);
			//! checks if element (r,c) is non-empty
			//! @param r the row to access
			//! @param c the column to access
			//! @return wether or not (r,c) exists 
		inline bool exists(const size_t r, const size_t c) const;
			//! number of rows
			//! @return row number
		inline size_t numRows() const;
			//! number of columns
			//! @return column number
		inline size_t numColumns() const;
			//! access to a matrix row
			//! @param row the row to access
			//! @return the row vector
		std::vector<T> rowVec(const size_t row) const;
			//! access to a matrix column
			//! @param col the column to access
			//! @return the column vector
		std::vector<T> columnVec(const size_t col) const;
			//! access to the non-empty matrix row entries
			//! @param row the row to access
			//! @return the non-empty row entries
		EntryMap rowValues(const size_t row) const;
			//! access to the non-empty matrix column entries
			//! @param col the column to access
			//! @return the non-empty column entries
		EntryMap columnValues(const size_t col) const;
			//! Resizes the matrix. If the the dimensions are increased, the 
			//! new fields will be initialised with current default value
			//! @param row the new row number
			//! @param col the new column number
		void resize(const size_t row, const size_t col);
			//! Resizes the matrix. If the the dimensions are increased, 
			//! sets the default value of new fields to the given value
			//! @param row the new row number
			//! @param col the new column number
			//! @param defVal the default value for new matrix elements
		void resize(const size_t row, const size_t col, const T& defVal);
			//! Sets the default value for new fields
			//! @param defVal the default value for new matrix elements
		void setDefaultValue( const T& defVal );
			//! Access to the default value for new fields
			//! @return the default value for new matrix elements
		const T& getDefaultValue( void ) const;
			//! destruction and freeing memory 
		~MatrixSparseC();
	};

} // namespace biu

//! Prints the matrix elements blank (' ') separated and linewise to stream.
//! @param out the output stream to write to
//! @return the changed output stream out
template< class T >
std::ostream&
operator << (std::ostream& out, const biu::MatrixSparseR<T> & m);

//! Prints the matrix elements blank (' ') separated and linewise to stream.
//! @param out the output stream to write to
//! @return the changed output stream out
template< class T >
std::ostream&
operator << (std::ostream& out, const biu::MatrixSparseC<T> & m);


#include "biu/MatrixSparse.icc"
	
#endif // define BIU_MATRIX_SPARSE_HH
