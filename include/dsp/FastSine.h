// FastSine.h
// 


#ifndef FASTSINE_H
#define FASTSINE_H 

#include <cmath>

namespace soundmath {
	//! A fast generator of sinusoids
	template <typename T>
	class FastSine {
	public:	
		FastSine (T sr, T freq, T amp) {
			m_sr = sr;
			frequency (freq);
			m_z1 = .5;
			m_z2 = 0.;	
			amplitude (amp);
		}
		virtual ~FastSine () {}
		void frequency (T f) {
			m_a = 2.f * (T) sin (M_PI * f / m_sr);
	
		}
		void amplitude (T a) {
			m_amp = a;
		}
		T* process (T* out, int len) {
			for (int i = 0; i < N; ++i) {
				m_z1 -= m_a * m_z2;
				m_z2 += m_a * m_z1;				
				out[i] = m_z1 * m_amp;
			}

			return out;
		}
	private:
		T m_sr;
		T m_a;
		T m_z1;
		T m_z2;
		T m_amp;
	};
}

#endif	// FASTSINE_H 

// EOF
