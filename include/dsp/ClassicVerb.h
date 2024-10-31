// ClassicVerb.h
// 

#ifndef CLASSICVERB_H
#define CLASSICVERB_H 

#include "FIR.h"
#include "Delay.h"
#include "Comb.h"
#include "Allpass.h"

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace soundmath {
//#define undenormalise(sample) if(((*(unsigned int*)&sample)&0x7f800000)==0) sample=0.0f

	static const int COMBS_LEN[] = {
		1116,
		1188,
		1277,
		1356,
		1422,
		1491,
		1557,
		1617,
		1663,
		1699,
		1741,
		1801,
		1867,
		1889, 
		1931, 
		1973, 
		1993, 
		2017,
		2113,
		2293,
		2467,
		2647
	};
	static const int MAXCOMBS = sizeof (COMBS_LEN) / sizeof (int);

	static const int ALLPASS_LEN[] = {
		556,
		441,
		341,
		225,
		181,
		163,
		151,
		137,
		113,
		91,
		68,
		53,
		37
	};
	static const int MAXALLPASS = sizeof (ALLPASS_LEN) / sizeof (int);

	//! Classic reverberation tool based on comb and allpass filters
	template <typename T>
	class ClassicVerb {
	private:
		ClassicVerb& operator= (ClassicVerb&);
		ClassicVerb (const ClassicVerb&);
	public: 
		ClassicVerb (float sr, int bsize, int density, int diffusion, int offset, float ratio = 1.) :
			m_t60 (2.3), m_minDamp (.2), m_maxDamp (.3), m_early (.35), m_direct (.7), m_tail (.3) {		
			m_sr = sr;
			if (m_sr <= 0) throw std::runtime_error ("invalid sr specified");
			if (density <= 0 || density > MAXCOMBS) {
				std::stringstream err;
				err << "too many density units; max allowed is " << MAXCOMBS;
				throw std::runtime_error (err.str ());
			}
			if (diffusion <= 0 || diffusion > MAXALLPASS) {
				std::stringstream err;
				err << "too many diffusion units; max allowed is " << MAXCOMBS;
				throw std::runtime_error (err.str ());			
			}
			m_ratio = ratio;
			if (m_ratio <= 0) throw std::runtime_error ("invalid ratio specified");
			m_offset = offset;
			if (offset < 0) throw std::runtime_error ("invalid offset specified");
			m_size = bsize;
			if (m_size <= 0) throw std::runtime_error ("invalid size specified");
		
			alloc (density, diffusion, m_size);
			m_earlyRef = 0;
			earlyRef (0.05, 18, .9999);
			m_predelay = 0;
			predelay (0.0001);
		
			m_rescale = 0.001 + (1. / (diffusion * diffusion + density)) / diffusion;
		}
		virtual ~ClassicVerb () {
			for (unsigned int i = 0; i < m_combs.size (); ++i) {
				delete m_combs[i];
				delete [] m_cbuff[i];
			}
		
			for (unsigned int i = 0; i < m_allpass.size (); ++i) {
				delete m_allpass[i];
				delete [] m_abuff[i];
			}
			delete m_earlyRef;
			delete m_predelay;
	
			delete [] m_delay;
			delete [] m_tmp0;
			delete [] m_tmp1;
			delete [] m_tmp2;
		}
		void earlyRef (float dur, int density, float diffusion) {
			if (m_earlyRef) {
				delete m_earlyRef;
			}
			if (dur < 0) throw std::runtime_error ("invalid total duration for early reflections");

			int samples = (int) (dur * m_sr);
			if (samples < density) samples = density;
			if (samples <= 0) samples = 1;
				
			T* c = new T[samples];
			T* d = new T[samples];
			for (int i = 0; i < samples; ++i) c[i] = 0.;
			c[samples - 1] = 1. / exp (1.);

			int sign = 1;
			for (int i = 0; i < density; ++i) {
				int r = (int) frand (0, samples);
				c[r] = (float) sign / exp ((T) r / samples);
				sign -= sign;
			}		

			Allpass <T> alp (m_sr, dur, diffusion);
			alp.process (c, d, samples);
			//for (int i = 0; i < samples; ++i) std::cout << c[i] << std::endl;
			for (int i = 0; i < samples; ++i) {
				c[i] += d[i];
			}
			m_earlyRef = new FIR<T> (c, samples);
			delete [] c;
			delete [] d;
		}	
		void earlyRef (float dur, const T* coefs) {
			if (m_earlyRef) {
				delete m_earlyRef;
			}
			if (dur < 0) throw std::runtime_error ("invalid total duration for early reflections");

			int samples = dur * m_sr;
			if (samples <= 0) samples = 1;
				
			float* c = new float[samples];

			for (int i = 0; i < samples; ++i) {
				c[i] = coefs[i];
			}		

			m_earlyRef = new FIR<T> (c, samples);
			delete [] c;
		}	
		void predelay (float dur) {
			if (m_predelay) {
				delete m_predelay;
			}
			if (dur < 0) throw std::runtime_error ("invalid predelay specified");
			m_predelay = new Delay <T> (m_sr, dur, 0);
		}
		void t60 (float t) {
			m_t60 = t;
			if (m_t60 <= 0) throw std::runtime_error ("invalid t60 specified");
		
			for (unsigned int i = 0; i < m_combs.size (); ++i) {
				float s = (float) ((COMBS_LEN[i] + m_offset)) * m_ratio;
				float fb = pow (10.0, (-3.0 * s / (m_t60 * m_sr)));
				m_combs[i]->feedback (fb);
			}
		}
		void damp (float min, float max) {
			if (min < 0 || max < 0 || min > 1 || max > 1 || min > max) {
				throw std::runtime_error ("invalid damping specified");
			}
			for (unsigned int i = 0; i < m_combs.size (); ++i) {
				m_combs[i]->damp (frand (min, max));
			}
		}
		void gains (float direct, float early, float tail) {
			m_direct = direct;
			if (m_direct < 0 || m_direct > 1) throw std::runtime_error ("invalid gain for direct path");
			m_early = early;
			if (m_early < 0 || m_early > 1) throw std::runtime_error ("invalid gain for early reflections");
			m_tail = tail;
			if (m_tail < 0 || m_tail > 1) throw std::runtime_error ("invalid gain for tail reverberation");		
		}
		virtual T* process (const T* input, T* output) {
			for (int j = 0; j < m_size; ++j) {
				m_delay[j] = input[j] * m_rescale;
				output[j] = input[j] * m_direct;
			}
		
			m_predelay->process (m_delay, m_tmp0, m_size);
			//m_earlyRef->process (m_tmp0, m_tmp1, m_size);
						
			for (unsigned int i = 0; i < m_combs.size (); ++i) {
				m_combs[i]->process (m_tmp0, m_cbuff[i], m_size);
			}
		
			memset (m_tmp2, 0, sizeof (T) * m_size);
			for (unsigned int i = 0; i < m_combs.size (); ++i) {
				for (int j = 0; j < m_size; ++j) {
					m_tmp2[j] += m_cbuff[i][j];
				}        
			}
	
			m_allpass[0]->process (m_tmp2, m_abuff[0], m_size);
			for (unsigned int i = 1; i < m_allpass.size (); ++i) {
				m_allpass[i]->process (m_abuff[i - 1], m_abuff[i], m_size);
			}
		
			for (int j = 0; j < m_size; ++j) {
				output[j] += (m_abuff[m_allpass.size () - 1][j] * m_tail + (m_tmp1[j] * m_early));
				//output[j] += m_tmp1[j];
			}   

			return output;  	
		}
	private:
		void alloc (int density, int diffusion, int bsize) {
			for (int i = 0; i < density; ++i) {
				float s = (float) ((COMBS_LEN[i] + m_offset)) * m_ratio * m_sr / 44100.;
				float fb = pow (10.0, (-3.0 * s / (m_t60 * m_sr)));
				Comb<T>* c = new Comb<T> (m_sr, (int) s, fb);
				c->damp (frand (m_minDamp, m_maxDamp));
				m_combs.push_back (c);
			
				T* t = new T[bsize];
                memset (t, 0, sizeof (T) * bsize);
				m_cbuff.push_back (t);
			}
			for (int i = 0; i < diffusion; ++i) {
				m_allpass.push_back (new Allpass<T> (m_sr, ALLPASS_LEN[i], .7));		
				T* t = new T[bsize];
                memset (t, 0, sizeof (T) * bsize);
				m_abuff.push_back (t);
			}
		
			m_delay = new T[bsize];
            memset (m_delay, 0, sizeof (T) * bsize);
			m_tmp0 = new T[bsize];
            memset (m_tmp0, 0, sizeof (T) * bsize);
			m_tmp1 = new T[bsize];
            memset (m_tmp1, 0, sizeof (T) * bsize);
			m_tmp2 = new T[bsize];
            memset (m_tmp2, 0, sizeof (T) * bsize);
		}
	
		float frand (float min, float max) {
			float f = 0;
			short r = abs (rand ());
			f = fabs ((float) r / 32768);
			f *= (max - min);
			f += min;
			return f;
		}

		float m_sr;
		float m_t60;
		float m_minDamp;
		float m_maxDamp;
		float m_early;
		float m_direct;
		float m_tail;
		float m_rescale;
		float m_ratio;
		int m_offset;
		int m_size;
		std::vector <Comb <T>* > m_combs;
		std::vector <Allpass <T>* > m_allpass;
		std::vector <T*> m_cbuff;
		std::vector <T*> m_abuff;
		T* m_delay;
		T* m_tmp0;
		T* m_tmp1;
		T* m_tmp2;
	
		Delay<T>* m_predelay;
		FIR<T>* m_earlyRef; 
	};
}
#endif	// CLASSICVERB_H 

// EOF
