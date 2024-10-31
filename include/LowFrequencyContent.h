// LowFrequencyContent.h
// 

#ifndef LOWFREQUENCYCONTENT_H
#define LOWFREQUENCYCONTENT_H 

#include "VariState.h"
#include "EnvelopeFollower.h"

namespace soundmath {
	template <typename T>
	//! Computes low frequency content by filtering at the given frequency with 36 dB rolloff and computing envelope
	class LowFrequencyContent {
	private:
		LowFrequencyContent& operator= (LowFrequencyContent&);
		LowFrequencyContent (const LowFrequencyContent&);
	public: 
		LowFrequencyContent (double sr, double maxFreq) {
			m_max = maxFreq;
			m_ft1 = new VariState<T> (sr, maxFreq);
			m_ft2 = new VariState<T> (sr, maxFreq);
			m_ft3 = new VariState<T> (sr, maxFreq);
			m_env = new EnvelopeFollower<T> ();
		}
		virtual ~LowFrequencyContent () {
			delete m_ft1;
			delete m_ft2;
			delete m_env;
		}
		T process (const T* in, T* out, int len) {
			m_ft1->process (in, out, len);
			m_ft2->process (out, out, len);
			m_ft3->process (out, out, len);
			
			return m_env->process (out, len);
		};
	private:
		double m_max;
		VariState <T>* m_ft1;
		VariState <T>* m_ft2;
		VariState <T>* m_ft3;
		EnvelopeFollower <T>* m_env;
	};
}

#endif	// LOWFREQUENCYCONTENT_H 

// EOF
