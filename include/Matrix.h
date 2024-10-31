// Matrix.h
//

#ifndef MATRIX_H
#define MATRIX_H

#include "MatExpr.h"
#include <stdexcept>
#include <iostream>
#include <cmath>

namespace soundmath {
	template <class T>
	//! A matrix based on expression metatemplates
	class Matrix : public MatrixBase<T, Matrix<T> > {
	private:
		int m_rows, m_columns, m_totalSize;
		T* m_data;
		T** m_m;
	public:
		typedef T value_type;

		explicit Matrix (int rows, int columns) : MatrixBase<T, Matrix<T> >() {
			#ifdef RANGECHECK
				if (0>= rows || 0 >= columns) {
					throw std::range_error ("invalid numbers of rows/columns specified");
				}
			#endif
			m_rows = rows;
			m_columns = columns;
			m_totalSize  = (m_rows * m_columns);
			m_data  = new T[m_totalSize];
			m_m = new T*[m_rows];
			for (int i = 0; i < m_rows; ++i) m_m[i] = m_data +  i * m_columns;
			null ();
		}
		Matrix () : MatrixBase<T, Matrix<T> >() {
			m_rows = 0;
			m_columns = 0;
			m_totalSize = 0; 
			m_data = 0;
			m_m = 0;
		}		
		Matrix (const Matrix<T>& x) : MatrixBase<T, Matrix<T> >(),
			m_rows (x.m_rows), m_columns (x.m_columns), 
			m_totalSize (m_rows * m_columns), m_data (new T[m_totalSize]),
			m_m (new T*[m_rows]) {
			for (int i = 0; i < m_rows; ++i) m_m[i] = m_data + i * m_columns;
			assignFrom (x);
		}
		~Matrix() {
			delete [] m_m;
			delete [] m_data;
		}
		// assignments
		Matrix<T>& operator= (const Matrix<T>& rhs) {
			return assignFrom (rhs);
		}
		template <class X> Matrix<T>& operator= (const Expr2<T, X>& rhs) {
			return assignFrom (rhs);
		}
		template <class M> Matrix<T>& operator= (const MatrixBase<T, M>& rhs) {
			return assignFrom (rhs);
		}
		Matrix<T>& operator= (T rhs) {
			return assignFrom (rhs);
		}
		template <class Closure> Matrix<T>& operator= (const Closure& rhs) {
			rhs.assignTo (m_m);
			return *this;
		}
		// getters/setters
		T& operator() (int i, int j) {
		#ifdef RANGECHECK
			if (i < 0 || m_rows <= i) throw std::range_error ("out of row range");
			if (j < 0 || m_columns <= j) throw std::range_error ("out of column range");
		#endif
			return m_m[i][j];
		}
		T operator() (int i, int j) const {
		#ifdef RANGECHECK
			if (i < 0 || m_rows <= i) {
				throw std::range_error ("out of row range");
			}
			if (j < 0 || m_columns <= j) {
				throw std::range_error ("out of column range");
			}
		#endif
			return m_m[i][j];
		}
		T* operator[] (int i) {
		#ifdef RANGECHECK
			if (i < 0 || m_rows <= i) throw std::range_error ("out of row range");
		#endif
			return m_m[i];
		}
		const T* operator[] (int i) const {
		#ifdef RANGECHECK
			if (i < 0 || m_rows <= i) {
				throw std::range_error ("out of row range");
			}
		#endif
			return m_m[i];
		}	
		// others
		int size () const {
			return m_totalSize;
		}
		void clear () {
			delete [] m_m;
			delete [] m_data;
			
			m_rows = 0;
			m_columns = 0;
			m_totalSize = 0; 
			m_data = 0;
			m_m = 0;		
		}
		void resize (int i, int j) {
		#ifdef RANGECHECK
			if (0 >= i || 0 >= j) {
				throw std::runtime_error ("invalid rows/columns specified");
			} 
		#endif
			Matrix<T> temp = (*this);
			
			delete [] m_m;
			delete [] m_data;
			m_rows = i;
			m_columns = j;
			m_totalSize = (i * j);
			
			m_data = new T[m_totalSize];
			m_m  = new T*[m_rows];
			for (int ii = 0; ii < m_rows; ++ii) m_m[ii] = m_data +  ii * m_columns;
			null ();			

			for (int ii = 0; ii < std::min (temp.rows (), i); ++ii) {
				for (int jj = 0; jj < std::min (temp.cols (), j); ++jj) {
					m_m[ii][jj] = temp[ii][jj];
				}
			}
		}		
		int rows () const {
			return m_rows;
		}
		int cols () const {
			return m_columns;
		}		
		// set null to all elements of this matrix
		void null () {
			for (int i = 0; i < m_rows; i++) {
				for (int j = 0; j < m_columns; j++) m_m[i][j] = 0;
			}
		}		
		// set this matrix to unity
		inline void unit () {
			int row = std::min (m_rows, m_columns);
			m_rows = m_columns = row;
		
			for (int i = 0; i < m_rows; ++i) {
				for (int j = 0; j < m_columns; ++j) {
					m_m[i][j] = i == j ? T(1) : T(0);
				}
			}
			return;
		}
		T** data () {
			return m_m;
		}
		// algebra
		// unary inversion operator
		inline Matrix<T> inv () {
			T a1, a2, *rowptr;
		
			if (m_rows != m_columns) {
				throw std::runtime_error ("cannot invert a non-square matrix");
			}
			Matrix<T> temp (m_rows, m_columns);
			temp.unit ();
			for (int k = 0; k < m_rows; ++k) {
				int indx = pivot (k);
				if (indx == -1) {
					throw std::runtime_error ("cannot invert a singular matrix");
				}
		
				if (indx != 0) {
					rowptr = temp.m_m[k];
					temp.m_m[k] = temp.m_m[indx];
					temp.m_m[indx] = rowptr;
				}
				a1 = m_m[k][k];
				for (int j = 0; j < m_rows; ++j) {
					m_m[k][j] /= a1;
					temp.m_m[k][j] /= a1;
				}
				for (int i = 0; i < m_rows; ++i) {
					if (i != k) {
						a2 = m_m[i][k];
						for (int j = 0; j < m_rows; ++j) {
							m_m[i][j] -= a2 * m_m[k][j];
							temp.m_m[i][j] -= a2 * temp.m_m[k][j];
						}
					}
				}
			}
			return temp;
		}		
		// calculate the determinant of a matrix
		inline T det () const {
			T piv,detVal = T(1);
		
			if (m_rows != m_columns) {
				throw std::runtime_error ("cannot compute determinant a non-square matrix");
			}
		
			Matrix<T> temp (*this);
			for (int k = 0; k < m_rows; ++k) {
				int indx = temp.pivot (k);
				if (indx == -1) return 0;
				if (indx != 0) detVal = - detVal;
				detVal = detVal * temp.m_m[k][k];
				for (int i = k + 1; i < m_rows; ++i) {
					piv = temp.m_m[i][k] / temp.m_m[k][k];
					for (int j = k + 1; j < m_rows; ++j)
						temp.m_m[i][j] -= piv * temp.m_m[k][j];
				}
			}
		return detVal;
		}
		// calculate the norm of a matrix (Frobenius)
		inline T norm () {
			T retVal = T(0);
		
			for (int i = 0; i < m_rows; ++i) {
				for (int j = 0; j < m_columns; ++j) {
					retVal += m_m[i][j] * m_m[i][j];
				}
			}
			retVal = sqrt (retVal);
			return retVal;
		}		
	private:
		// private partial pivoting method
		int pivot (int row) {
			int k = row;
			double amax, temp;
		
			amax = -1;
			for (int i = row; i < m_rows; i++) {
				if ((temp = fabs ((T) m_m[i][row])) > amax && temp != 0.0) {
					amax = temp;
					k = i;
				}
			}
			if (m_m[k][row] == T(0)) return -1;
			if (k != row) {
				T* rowptr = m_m[k];
				m_m[k] = m_m[row];
				m_m[row] = rowptr;
				return k;
			}
			return 0;
		}
	};
	// overloaded operators
	// unary transpose operator
	template <typename T>
	inline Matrix<T> operator~ (const Matrix<T>& m) {
		Matrix<T> temp (m.cols (), m.rows ());
	
		for (int i = 0; i < m.rows (); ++i)
			for (int j=0; j < m.cols (); ++j) {
				T x = m (i, j);
				temp (j, i) = x;
			}
		return temp;
	}
	// unary inversion operator
	template <typename T>
	inline Matrix<T> operator! (const Matrix<T>& m) {
		Matrix<T> temp = m;
		return temp.inv ();
	}	
	template <typename T>	
	// logical equal-to operator
	inline bool operator== (const Matrix<T>& m1, const Matrix<T>& m2) {
		if (m1.rows () != m2.rows () || m1.cols () != m2.cols ()) {
			return false;
		}
	
		for (int i = 0; i < m1.rows (); ++i) {
			for (int j = 0; j < m1.cols (); ++j) {
				if (m1 (i,j) != m2 (i,j)) {
					return false;
				}
			}
		}
	
		return true;
	}
	// logical no-equal-to operator
	template <typename T>
	inline bool	operator != (const Matrix<T>& m1, const Matrix<T>& m2) {
		return (m1 == m2) ? false : true;
	}	
}

#endif	// MATRIX_H

// EOF
