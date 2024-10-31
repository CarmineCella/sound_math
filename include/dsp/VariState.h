// VariState.h// #ifndef VARISTATE_H#define VARISTATE_H #include "AbstractFilter.h"#include "filterType.h"#include <cmath>namespace soundmath {	template <typename T>	//! Variable-state filter (Chamberlin version): for each buffer returns a buffer of same size	class VariState : public AbstractFilter {	private:		VariState& operator= (VariState&);		VariState (const VariState&);	public:		VariState (T sr, T cutoff = 440,			T quality = .707, filterType t = LOPASS) {			type (t);						if (cutoff > 0) {				m_cutoff = cutoff;				m_f1 = 2 * sin (M_PI * m_cutoff / sr);			}			else {				m_f1 = 2 * sin (M_PI * 440 / sr);			}			if (quality > 0) {				m_q1 = 1 / quality;			}			else {				m_q1 = 1 / .707;			}			reset ();			m_sicvt = M_PI / sr;			m_sr  = sr;		}		virtual ~VariState () {}		void cutoff (T c) {			if (c > 0) {				m_cutoff = c;				m_f1 = 2 * sin (M_PI * m_cutoff / m_sr);			}		}		T cutoff () const {			return m_cutoff;		}		void quality (T q) {			if (q > 0) {				m_q1 = 1 / q;			}		}		void reset () {			m_data[0] = m_data[1] = m_data[2] = m_data[3] = 0;		}		void type (filterType t) {			switch (t) {				default:				case LOPASS:					m_type = 0;				break;				case HIPASS:					m_type = 1;				break;				case BANDPASS1:				case BANDPASS2:					m_type = 2;				break;				case NOTCH:					m_type = 3;				break;					}					}		T* process (const T* in, T* out, int len) {			for (int i = 0; i < len; ++i) {				m_data[0] = m_data[0] +	m_f1 * m_data[2];				m_data[1] = in[i] -	m_data[0] - m_q1 * m_data[2];				m_data[2] = m_f1 * m_data[1] + m_data[2];				m_data[3] = m_data[1] + m_data[0];				out[i] = m_data[m_type];			}			return out;		}	private:		T m_sr;		T m_cutoff;		int m_type;		T m_data[4];		T m_q1;		T m_f1;		T m_sicvt;	};}#endif	// VARISTATE_H // EOF