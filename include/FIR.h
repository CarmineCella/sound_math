// FIR.h
// 


#ifndef FIR_H
#define FIR_H

#include "AbstractFilter.h"

namespace soundmath {
	//! A FIR filter of any order (see also Filter)
	template <typename T>
	class FIR : public AbstractFilter<T> {
	public:
		FIR (const T* coefs, int ncoeff, float gain = 1.) {
			m_length = ncoeff;
			m_impulseResponse = new T[m_length];
			m_delayLine = new T[m_length];
			impulse (coefs, ncoeff);        
			m_gain = gain;
		}
		virtual ~FIR () {
			delete [] m_impulseResponse;
			delete [] m_delayLine;
		}
		virtual void process (const T* input, T* output, int size) {
			for (int j = 0; j < size; ++j) {
				output[j] = 0.;
				m_delayLine[0] = input[j] * m_gain;
				for (unsigned int i = m_length - 1; i > 0; i--) {
					output[j] += m_impulseResponse[i] * m_delayLine[i];
					m_delayLine[i] = m_delayLine[i - 1];
				}
				output[j] += m_impulseResponse[0] * m_delayLine[0];
			}
		}
		void impulse (const T* coefs, int ncoeff) {
			for (int i = 0; i < ncoeff; ++i) {
				m_impulseResponse[i] = coefs[i];
				m_delayLine[i] = 0.;
			}
		}
	private:
		T m_gain;
		int m_length;
		T* m_delayLine;
		T* m_impulseResponse;
	};
}
#endif	// FIR_H 

// EOF
