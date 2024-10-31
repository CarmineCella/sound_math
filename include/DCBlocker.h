// DCBlocker.h
// 

#ifndef DCBLOCKER_H
#define DCBLOCKER_H 

#include "AbstractFilter.h"
#include <cmath>

namespace soundmath {
	//! Basic DC blocker.
	template <typename T>
	class DCBlocker : public AbstractFilter {
	private:
		DCBlocker& operator= (DCBlocker&);
		DCBlocker (const DCBlocker&);
	public: 
		DCBlocker (T sr, T gain = .99) 
			: m_sr (sr), m_gain (gain), m_cutoff (10) {
			reset ();
		}

		virtual ~DCBlocker () {}

		void cutoff (T c) {
			m_cutoff = c;
			T theta = (2 * M_PI * m_cutoff) / m_sr;
			m_gain = 1 - theta;
		}
		void gain (T g) { m_gain = g; }
		void reset () {
			m_x1 = 0;
			m_y1 = 0;
		}

		void process (const T* in, T* out, int N) {
			for (int i = 0; i < N; ++i) {
				out[i] = in[i] - m_x1 + (m_gain * m_y1);
				m_y1 = out[i];
				m_x1 = in[i];
			}
		}
	
	private:
		T m_sr;
		T m_gain;
		T m_cutoff;
		T m_x1;
		T m_y1;
	};
}
#endif	// DCBLOCKER_H 

// EOF

