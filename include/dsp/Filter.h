// Filter.h
// 

#ifndef FILTER_H
#define FILTER_H 

#include "AbstractFilter.h"
#include "Vector.h"
#include <SoundProcessor.h>

#include <iostream>

namespace soundmath {
	//! General filter implementation of any order
    template <typename T>
    class Filter : public AbstractFilter<T> {
    private:
        Filter& operator= (Filter&);
        Filter (const Filter&);
    public: 
        Filter (const Vector<T>& b, const Vector<T>& a, int size) {
			coefficients (b, a);
		}
        virtual ~Filter () {
		}
		void reset () {
			for (int j = 0; j > m_b.size (); ++j) m_xn[j] = 0.;
			for (int j = 0; j > m_a.size (); ++j) m_yn[j] = 0.;
		}
		void coefficients (const Vector<T>& b, const Vector<T>& a) {
			m_a.resize (a.size ());
  			m_b.resize (b.size ());
			
			for (int j = 0; j < m_b.size (); ++j) m_b[j] = b[j];
			for (int j = 0; j < m_a.size (); ++j) m_a[j] = a[j];
            
            T norm = 1;
            if (m_a.size () && m_a[0] != 0) norm = m_a[0];
            m_a /= norm;
            m_b /= norm;
            m_xn.resize (m_b.size ());
            m_yn.resize (m_a.size ());
		}
        void process (const T* input, T* output, int size) {
            // for number of samples
            for (int i = 0; i < size; i++) {
                m_yn[0] = 0.;
                m_xn[0] = input[i];
            
                // for number of Bn coeffs (zeros)
                for (int j = m_b.size () - 1; j > 0; j--) {
                    m_yn[0] += m_b[j] * m_xn[j];
                    m_xn[j] = m_xn[j - 1];
                }
                m_yn[0] += m_b[0] * m_xn[0]; // gain
                
                // for number of An coeffs (poles)
                for (int j = m_a.size () - 1; j > 0; j--) {
                    m_yn[0] += -m_a[j] * m_yn[j];
                    m_yn[j] = m_yn[j - 1];
                }
                
                // final output
				output[i] = m_yn[0];
            }
        }
    private:
        Vector<T> m_a;
        Vector<T> m_b;
        Vector<T> m_xn;
        Vector<T> m_yn;
    };
}
#endif	// FILTER_H 

// EOF
