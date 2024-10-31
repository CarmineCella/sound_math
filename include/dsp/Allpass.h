// Allpass.h
// 

#ifndef ALLPASS_H
#define ALLPASS_H 

#include "Delay.h"

namespace soundmath {
	template <typename T>
	//! All-pass filter: for each buffer returns a buffer of same size
	class Allpass : public Delay <T> {
	public:
		typedef Delay <T> Base;
	    Allpass (double sr, double dur, double feedback) : 
			Delay <T> (sr, dur, feedback) {}
	    Allpass (double sr, int samples, double feedback) : 
			Delay <T> (sr, samples, feedback) {}		
	    virtual void process (const T* input, T* output, int size) {
	        for (long i = 0; i < size; ++i) {
				double buff = Base::m_delay[Base::m_ptr];
				//undenormalise (buff);
				output[i] = -(input[i]) + buff;
				Base::m_delay[Base::m_ptr++] = input[i] + (Base::m_feedback * buff);
				Base::m_ptr %= Base::m_samples;		
	
	        }
	    }
	};
}

#endif	// ALLPASS_H 

// EOF
