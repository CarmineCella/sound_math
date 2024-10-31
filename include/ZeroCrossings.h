// ZeroCrossings.h
// 


#ifndef ZEROCROSSINGS_H
#define ZEROCROSSINGS_H 

namespace soundmath {
	template <typename T>
	//! Zero-crossings detector: for each buffer returns the mean value
	class ZeroCrossings {
	public: 
		ZeroCrossings (T sr = 44100, T feedback = .8) {
			m_rate = 1. / sr;
			m_feedback = feedback;
			m_gain = 1. - m_feedback;
	
			m_val = 0;
			m_z = 0;
			m_posit = 0;
			m_negat = 0;		
		}
		virtual ~ZeroCrossings () {}
		T process (const T* in, int len) {
			long count = 0;
			T sum = 0;
			for (int i = 0; i < len; ++i) {
				m_val = -1;
				if (in[i] >= 0.f) {
					++m_posit;
		
					if (m_negat != 0) {
						m_z = m_z * m_feedback + (T) (m_posit + m_negat) * m_gain;
						m_val = (T) m_z; //1. / (m_z * m_rate);
						m_posit = 0;
						m_negat = 0;
					}
				}
			
				if (in[i] < 0.f) {
					++m_negat;
				}
				
				if (m_val > 0) {
					sum += m_val;
					++count;
				}
			}
			
			return (sum == 0 || count == 0) ? 0 : 1. / ((sum / count) * m_rate);
		}
	private:
		T m_feedback;
		T m_gain;
				
		long m_posit;
		long m_negat;
		T m_rate;
		T m_z;
		T m_val;
	};

}
#endif	// ZEROCROSSINGS_H 

// EOF
