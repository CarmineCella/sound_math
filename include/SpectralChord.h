// SpectralChord.h
// 

#ifndef SPECTRALCHORD_H
#define SPECTRALCHORD_H 

#include "HarmonicDesigner.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>

namespace soundmath {
	const char* methods[] = {"geometric", "companded harmonic", "random", "model based"};
	
	//! Basic struct to describe an oscillator for SpectralChord
	template <typename T>
	struct OscilDescriptor {
		T amp;
		T freq;
		T phi;
		T panL;
		T panR;
	};
	
	//! Basic struct to describe spectral generation for SpectralChord
	template <typename T>
	struct Parameters {
		T start;
		T dur;
		int partials;
		int thickness;
		T spread;
		T fund;
		T coeff;
		T slope;
	};
	
	template <typename T>
	class SpectralChord {
	private:
		SpectralChord& operator= (SpectralChord&);
		SpectralChord (const SpectralChord&);
	
	public:
		SpectralChord (T srate, T* tab, int tablen) {
			m_tab = tab;
			m_tablen = tablen;
			
			m_stereoSpread = 23;
			m_oscbank = 0;
			m_sr = srate;
			
			Parameters<T> params;
			params.start = 0; 
			params.dur = 10;
			params.fund = 110;
			params.partials = 5;
			params.thickness = 3;
			params.spread = .001;
			params.coeff = 1.32;
			params.slope = .5;
			
			HarmonicDesigner<T> designer;		
			design (&params, &designer);
			
			m_dline = new T[m_stereoSpread];
			memset (m_dline , 0, sizeof (T) * m_stereoSpread);
			m_ptr = 0;
		}
		virtual ~SpectralChord () {
			delete [] m_dline;
			if (m_oscbank) delete [] m_oscbank;
		}
		
		void design (Parameters<T>* params, AbstractDesigner<T>* designer) {
			if (m_oscbank) delete [] m_oscbank;
			m_oscbank = new OscilDescriptor<T>[params->partials * params->thickness];
			
			int rpartials = params->partials * params->thickness;
			T tamp = 0;
			T fcurr = params->fund;
			for (int i = 0; i < params->partials && fcurr < m_sr / 2; ++i) {
				T acurr = 1.;
				fcurr = designer->design (params->fund, params->coeff, i + 1, acurr);
							
				for (int k = 0; k < params->thickness; ++k)  {
					m_oscbank[i * params->thickness + k].freq = fcurr + k * (float) rand () / RAND_MAX * fcurr * params->spread;
					m_oscbank[i * params->thickness + k].amp = acurr;
					
					tamp += m_oscbank[i * params->thickness + k].amp;
					
					m_oscbank[i * params->thickness + k].phi = (float) rand () / RAND_MAX;
					m_oscbank[i * params->thickness + k].panL = (float) rand () / RAND_MAX;
					m_oscbank[i * params->thickness + k].panR = 1. - m_oscbank[i].panL;
				}
				rpartials = (i + 1) * params->thickness;
			}
			
			for (int i = 0; i < rpartials; ++i) {
				m_oscbank[i].amp /= tamp; 
			}
			m_partials = rpartials;
			
			m_startSample = (int) (params->start * m_sr);
			m_samples = (int) (params->dur * m_sr);
			m_eincr1 = 1. / ((T) m_samples * params->slope);
			m_eincr2 = 1. / (m_samples - ((T) m_samples * params->slope));
			m_envLimit = (m_samples * params->slope);
			m_time = 0;
			m_env = 0;
		}
		
		int totalPartials () const { return m_partials; }
		int start () const { return m_startSample; }
		int length () const { return m_samples; }
		
		void process (T* output, int n) {
			T sicvt = m_tablen / m_sr;
			T* l = output;
			T* r = output + 1;
			while (n--) {
				T tmp = 0;
				OscilDescriptor<T>* osc = m_oscbank;
				int pcurr = m_partials;
			
				//m_env = (m_time < m_envLimit ? m_env += m_eincr1 : m_env -= m_eincr2);
				//m_env = hanning (m_time, m_samples);
				m_env = blackmann (m_time, m_samples);
				
				while (pcurr--) {
					T incr = osc->freq * sicvt;
					int pos = (int) osc->phi;
					T frac = osc->phi - pos;
			
					T v = m_tab[pos] * (1. - frac) + m_tab[pos + 1] * frac;
					osc->phi += incr;
					if (osc->phi >= m_tablen) osc->phi -= m_tablen;
					v *= osc->amp;
						
					*l += v * osc->panL * m_env;
					tmp += v * osc->panR * m_env;
					++osc;
				};
				l += 2;
				
				*r += m_dline[m_ptr];
				m_dline[m_ptr++] = tmp;
				m_ptr %= m_stereoSpread;
				r += 2;
				
				++m_time;
				if (m_time > m_samples) m_time -= m_samples;
			}
		}
	private:
		T hanning (long pos, long samples) {
			return .5 - .5 * (cos ((2. * M_PI * (T) pos) / samples));
		}
	
		T blackmann (long pos, long samples) {
			return .42 - .5 * cos ((2. * M_PI * (T) pos) / (samples - 1)) + .08 *
				cos ((4. * M_PI * (T) pos) / (samples - 1));
		}
		T m_sr;
		int m_partials;
		
		OscilDescriptor<T>* m_oscbank;
		T* m_dline;
		int m_ptr;
		int m_stereoSpread;
		int m_startSample;
		int m_samples;
		
		T m_env;
		T m_eincr1;
		T m_eincr2;
		int m_envLimit;
		int m_time;
		
		T* m_tab;
		int m_tablen;
	};
}

#endif	// SPECTRALCHORD_H 

// EOF

