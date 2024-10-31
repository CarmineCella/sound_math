// EnvelopeFollowerBuffer.h
// 

#ifndef ENVELOPEFOLLOWERBUFFER_H
#define ENVELOPEFOLLOWERBUFFER_H 

#include <cmath>

namespace soundmath {
	template <typename T>
	//! Follows amplitude envelope of a signal: for each buffer returns a whole buffer
	class EnvelopeFollowerBuffer {
	public:
		EnvelopeFollowerBuffer (T feedback = .999) {
			m_val = 0;
			m_feedback = feedback;
			m_gain = 1. - m_feedback;
		}
		virtual ~EnvelopeFollowerBuffer () {}
		void reset () { m_val = 0; }
		void process (const T* in, T* out, int len, T& max, T& min) {
			max = in[0];
			min = in[0];
			for (int i = 0; i < len; ++i) {
				T abs = fabs (in[i]);
				m_val = out[i] = abs > m_val ? abs : m_val * m_feedback + abs * m_gain;
				if (max < m_val) max = m_val;
				if (min > m_val) min = m_val;
			}
		}
		void setMax (T m) { if (m > m_val) m_val = m; }
	private:
		T m_val;
		T m_feedback;
		T m_gain;
	};
}

#endif	// ENVELOPEFOLLOWERBUFFER_H 

// EOF
