// Vector.h
//

#ifndef VECTOR_H
#define VECTOR_H

#include "VecExpr.h"
#include <stdexcept>
#include <cmath>

namespace soundmath {
	
	template <class T>
	//! A vector based on expression metatemplates
	class Vector : public VectorBase<T, Vector<T> > {
	private:
		int m_size;
		T* m_data;
	public:
		typedef T value_type;
		typedef T* iterator;
		typedef const T* const_iterator;
		typedef std::reverse_iterator <const_iterator> const_reverse_iterator;
		typedef std::reverse_iterator <iterator> reverse_iterator;	

		explicit Vector (int size) : VectorBase<T, Vector<T> >() {
		#ifdef SIZECHECK
			if (0 >= size) {
				throw std::runtime_error ("invalid size specified");
			}
		#endif
			m_size = size;
			m_data = new T[m_size]; 		
			null ();
		}
		Vector () : VectorBase<T, Vector<T> >() {
			m_size = 0;
			m_data = 0; 
		}		
		Vector (const Vector<T>& x) : VectorBase<T, Vector<T> >(), 
			m_size (x.m_size), m_data (new T[m_size]) {
			assignFrom (x);
		}
		Vector (T* x, int size) : VectorBase<T, Vector<T> >() {
			if (0 >= size) {
				throw std::runtime_error ("invalid size specified");
			}
			m_size = size;
			m_data = new T[m_size];		
			assignFrom (x);
		}		
		~Vector() {
			delete [] m_data;
		}
		// assignments
		template <class X> Vector<T>& operator= (const Expr1<T, X>& rhs) {
			return assignFrom (rhs);
		}
		template <class V> Vector<T>& operator= (const VectorBase<T, V>& rhs) {
			return assignFrom (rhs);
		}
		Vector<T>& operator= (const Vector<T>& rhs) {
			return assignFrom (rhs);
		}
		Vector<T>& operator= (T rhs) {
			return assignFrom (rhs);
		}
		template <class Closure> Vector<T>& operator= (const Closure& rhs) {
			rhs.assignTo (m_data);
			return *this;
		}
		// getters/setters
		T& operator() (int i) {
			return m_data[i];
		}
		T operator() (int i) const {
			return m_data[i];
		}
		T& operator[] (int i) {
			return m_data[i];
		}
		T operator[] (int i) const {
			return m_data[i];
		}
		const T* data () {
			return m_data;
		}
		// iterators
		iterator begin () { 
			return iterator (m_data); 
		}
		const_iterator begin () const { 
			return const_iterator (m_data); 
		}
		iterator end ()	{ 
			return iterator (m_data + m_size); 
		}
		const_iterator end () const	{ 
			return const_iterator (m_data + m_size); 
		}
		reverse_iterator rbegin () { 
			return reverse_iterator (end ());
		}
		const_reverse_iterator rbegin () const { 
			return const_reverse_iterator (end ()); 
		}
		reverse_iterator rend () { 
			return reverse_iterator (begin ()); 
		}
		const_reverse_iterator rend () const { 
			return const_reverse_iterator (begin ());
		}	
		// others
		int size () const {
			return m_size;
		}
		void clear () {
			delete [] m_data;
			m_data = 0;
			m_size = 0;
		}
		void resize (int n) {
		#ifdef SIZECHECK
			if (0 >= n) {
				throw std::runtime_error ("invalid size specified");
			} 
		#endif
			Vector<T> temp = (*this);
			
			delete [] m_data;
			m_data = new T[m_size = n];
			null ();
			
			for (int i = 0; i < std::min (temp.size (), n); ++i) {
				m_data[i] = temp[i];
			}
		}		
		void push_back (const T& v) {
			resize (m_size + 1);
			m_data[m_size - 1] = v;
		}
		void pop_back () {
			if (1 >= m_size) {
				throw std::runtime_error ("cannot pop_back last element");
			}
			resize (m_size - 1);
		}
		void null () {
			for (int i = 0; i < m_size; ++i) {
				m_data[i] = (T) 0;
			}
		}	
		// algebra
		// calculate the norm of a vector (Frobenius)
		inline T norm () {
			T retVal = T(0);
		
			for (int i = 0; i < m_size; ++i) {
				retVal += m_data[i] * m_data[i];
			}

			retVal = sqrt (retVal);
			return retVal;
		}		
	};
	// overloaded operators
	template <typename T>	
	// logical equal-to operator
	inline bool operator== (const Vector<T>& v1, const Vector<T>& v2) {
		if (v1.size () != v2.size ()) {
			return false;
		}
	
		for (int i = 0; i < v1.size (); ++i) {
			if (v1 (i) != v2 (i)) {
				return false;
			}
		}
	
		return true;
	}
	// logical no-equal-to operator
	template <typename T>
	inline bool	operator != (const Vector<T>& v1, const Vector<T>& v2) {
		return (v1 == v2) ? false : true;
	}		
}

#endif // VECTOR_H

// EOF


