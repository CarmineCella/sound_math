// EnvelopeFollower.h
// 

#ifndef ENVELOPEFOLLOWER_H
#define ENVELOPEFOLLOWER_H 

#include <cmath>

namespace soundmath {
	template <typename T>
	//! Follows amplitude envelope of a signal: for each buffer returns the max value
	class EnvelopeFollower {
	public:
		EnvelopeFollower (double feedback = .999) {
			m_val = 0;
			m_feedback = feedback;
			m_gain = 1. - m_feedback;
		}
		virtual ~EnvelopeFollower () {}
		void reset () { m_val = 0; }
		T process (const T* in, int len) {
			T max = in[0];
			for (int i = 0; i < len; ++i) {
				double abs = fabs (in[i]);
				m_val = abs > m_val ? abs : m_val * m_feedback + abs * m_gain;
				if (max < m_val) max = m_val;
			}
			return max;
		}
		void setMax (T m) { if (m > m_val) m_val = m; }
	private:
		T m_val;
		double m_feedback;
		double m_gain;
	};
}
#endif	// ENVELOPEFOLLOWER_H 

// EOF
