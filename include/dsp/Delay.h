// Delay.h
// 

#ifndef DELAY_H
#define DELAY_H 

#include "AbstractFilter.h"

namespace soundmath {
	template <typename T>
	//! Single-tapped delay line: for each buffer returns a buffer of same size
	class Delay : public AbstractFilter<T> {
	public:
	    Delay (double sr, double dur, double feedback) : 
			m_feedback (feedback), m_ptr (0) {
	        m_samples = (long) (sr * dur);
			allocate ();
	    }
	    Delay (double sr, int samples, double feedback) : 
			m_feedback (feedback), m_ptr (0) {
	        m_samples = samples;
			allocate ();
	    }	
	    virtual ~Delay () {
	        delete [] m_delay;
	    }
//		const T& tap (long pos) const {
//			return m_delay[pos];
//		}
	    long length () const { return m_samples; }
	    void reset () {
	        for (long i = 0; i < m_samples; ++i) {
	            m_delay[i] = (T) 0.;
	        }
	    }
		void feedback (double f) {
			m_feedback = f;
		}	
	    virtual void process (const T* input, T* output, int size) {
	        for (long i = 0; i < size; ++i) {
				output[i] = m_delay[m_ptr];
				//undenormalise (output[i]);
				m_delay[m_ptr++] = input[i] +  (m_feedback * output[i]);
				m_ptr %= m_samples;		
	        }
	    }
	protected:
	    T* m_delay;
	    double m_feedback;
		long m_samples;
	    long m_ptr;
	private:
		void allocate () {
			m_delay = new T[m_samples];
	        reset ();
		}
	};
}
#endif	// DELAY_H 

// EOF
