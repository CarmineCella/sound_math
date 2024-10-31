// FreqShift.h
// 

#ifndef FREQSHIFT_H
#define FREQSHIFT_H 

#include "AbstractFilter.h"
#include "Biquad.h"
#include "OscillatorF.h"
#include <cstring>

#define TABLENFS 1024
#define TABLEN 4096

namespace soundmath {
	//! A frequency shifter based on single-band modulation	
	template <typename T>
	class FreqShift : public AbstractFilter<T> {
	private:
		FreqShift& operator= (FreqShift&);
		FreqShift (const FreqShift&);
	public: 
		FreqShift (T sr, int maxBufferSize) {
			
			//memory allocation
			m_tmp1 = new T[maxBufferSize];
			m_tmp2 = new T[maxBufferSize];
			m_tmp3 = new T[maxBufferSize];
			m_tmp4 = new T[maxBufferSize];
			
			m_b1.reset (sr, ALLPASS, 2749, 0.24, 0, false);
			m_b2.reset (sr, ALLPASS, 600, 0.24, 0, false);
					
			m_tab = new T[TABLENFS + 1];
			std::memset (m_tab, 0, sizeof (T) * TABLENFS + 1);
			T a = 1;
			OscillatorF<T>::gen (m_tab, TABLENFS, &a, 1, 0);
			m_tab[TABLENFS] = m_tab[0];
			
			m_osc1 = new OscillatorF<T> (sr, m_tab, TABLENFS);
			m_osc2 = new OscillatorF<T> (sr, m_tab, TABLENFS);
			
			m_osc1->phase (0);
			m_osc2->phase (0.25);
			
			//shiftAmount (1);
			
			//amp
			m_ampBuffer	= new T[maxBufferSize];
			m_tabamp = new T[TABLEN + 1]; // GUARD POINT
			std::memset (m_tabamp, 0, sizeof (T) * TABLEN + 1);
			
			T a_amp = 1;
			OscillatorF<T>::gen (m_tabamp, TABLEN, &a_amp, 1, 0);
			m_tabamp[TABLEN] = m_tabamp[0];
			m_amp = new OscillatorF<T> (sr, m_tabamp, TABLEN);
			
			frequency(1);
			amount(0);
			
		}
		
		virtual ~FreqShift () {
			delete m_osc1;
			delete m_osc2;
					
			delete [] m_tab;
			
			if (m_tmp1) delete[] m_tmp1;
			if (m_tmp2) delete[] m_tmp2;
			if (m_tmp3) delete[] m_tmp3;
			if (m_tmp4) delete[] m_tmp4;
			
			
			if (m_amp)	delete m_amp;
			if (m_tabamp) delete [] m_tabamp;
			if (m_ampBuffer) delete [] m_ampBuffer;
			
		}
		
		void amount (T amp) {
			m_vamp = amp;
		}
		
		T amount () {
			return m_vamp;
		}
		
		void frequency (T freq) {
			m_amp->frequency (freq);
		}
		
		T frequency () {
			return m_amp->frequency();
		}
		
		T* process (const T* in, T* out, int len) {
			
			//amp
			m_amp->process (m_ampBuffer, len);
			T f = m_ampBuffer[0];
			m_osc1->frequency (f * m_vamp);
			m_osc2->frequency (f * m_vamp);
			
			m_b1.process (in, m_tmp1, len); // rephase 1
			m_b2.process (in, m_tmp2, len); // rephase 2
			
			m_osc1->process (m_tmp3, len); // sin w
			m_osc2->process (m_tmp4, len); // cos w
			
			for (int i = 0; i < len; ++i) {
				out[i] = (m_tmp1[i] * m_tmp3[i] + m_tmp2[i] * m_tmp4[i]);
			}
					
			return out;
		}
	
	private:
		
		
		Biquad<T> m_b1;
		Biquad<T> m_b2;
		OscillatorF<T>* m_osc1;
		OscillatorF<T>* m_osc2;
		T* m_tmp1;
		T* m_tmp2;
		T* m_tmp3;
		T* m_tmp4;
		T* m_tab;
		int m_len;
		
		//osc for frequency
		OscillatorF<T>* m_amp;
		T* m_tabamp;
		T m_vamp;
		T* m_ampBuffer;
	};
}

#endif	// FREQSHIFT_H 

// EOF
