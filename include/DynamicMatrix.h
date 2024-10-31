// DynamicMatrix.h
// 

#ifndef DYNAMICMATRIX_H
#define DYNAMICMATRIX_H 

//#include "Vector.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <vector>

#include <iostream>

namespace soundmath {
	//! A dynamic sized matrix (vector of vectors) -------------------------------- //
	template <typename T>
	class DynamicMatrix {
        std::vector<std::vector<T> > m_data;
	public:
		std::vector<T>& operator[] (unsigned i) {
			return m_data[i];
		}
		unsigned int size () const {
			return m_data.size ();
		}
		void resize (int v) {
			m_data.resize (v);
		}
		void clear () {
			for (int i = 0; i < m_data.size (); ++i) {
				m_data[i].resize (0);
			}
		}
		void push_back (std::vector<T>& v) {
			m_data.push_back (v);
		}
		int rows () const { return m_data.size (); }
	};
}

#endif	// MATRIX_H 

// EOF
