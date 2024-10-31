// Compressor.h
//

#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include "AbstractFilter.h"
#include <cmath>

namespace soundmath {
	//! A simple dynamic compressor
	template <typename T>
	class Compressor : public AbstractFilter {
	private:
		Compressor& operator= (Compressor&);
		Compressor (const Compressor&);
	public:
		Compressor () {
			m_threshold = 1.f;
			m_attack = m_release = m_envelope_decay = 0.f;
			m_outputF = 1.f;
	
			m_transfer_A = 0.f;
			m_transfer_B = 1.f;
	
			m_env = 0.f;
			m_gain = 1.f;
			m_ratio = 1;
		}
	
		void threshold (T value) {
			m_threshold = value;
			m_transfer_B = m_outputF * pow (m_threshold, -m_transfer_A);
		}
		
		T threshold () {
			return m_threshold;
		}
	
	
		void ratio (T value) {
			m_transfer_A = value-1.f;
			m_transfer_B = m_outputF * pow (m_threshold, -m_transfer_A);
			m_ratio = value;
		}
	
		T ratio () {
			return m_ratio;
		}
	
		void attack (T value) {
			m_attack = exp (-1.f / value);
		}
	
	
		void release (T value) {
			m_release = exp (-1.f / value);
			m_envelope_decay = exp (-4.f / value); /* = exp(-1/(0.25*value)) */
		}
	
	
		void output (T value) {
			m_outputF = value;
			m_transfer_B = m_outputF * pow (m_threshold, -m_transfer_A);
		}
	
	
		void reset () {
			m_env = 0.f;
			m_gain = 1.f;
		}
	
		void process (const T *input, T *output, int frames) {
			T det, transfer_gain;
			for (int i = 0; i < frames; i++) {
				det = fabs (input[i]);
				det += 10e-30f; // anti-denormals
	
				m_env = det >= m_env ? det : det + m_envelope_decay * (m_env - det);
	
				transfer_gain = m_env > m_threshold ? pow ((T) m_env, (T) m_transfer_A) * m_transfer_B : m_outputF;
	
				m_gain = transfer_gain < m_gain ?
					   transfer_gain + m_attack  * (m_gain - transfer_gain):
					   transfer_gain + m_release * (m_gain - transfer_gain);
	
				output[i] = input[i] * m_gain;
			}
		}
	private:
		T   m_threshold;
		T   m_attack, m_release, m_envelope_decay;
		T   m_outputF;
		T   m_transfer_A, m_transfer_B;
		T   m_env, m_gain;
		T	m_ratio;
	};
}

#endif	// COMPRESSOR_H

// EOF
