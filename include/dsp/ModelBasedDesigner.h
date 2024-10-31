// ModelBasedDesigner.h
// 

#ifndef MODELBASEDDESIGNER_H
#define MODELBASEDDESIGNER_H 

#include "AbstractDesigner.h"
#include "FFT.h"

namespace soundmath {
	//! A spectrum designer based on model analysis
	template <typename T>
	class ModelBasedDesigner : public AbstractDesigner<T> {
	public:
		ModelBasedDesigner (int partials = 512) {
			m_peaks = new Peak<T>[partials];
			m_partials = partials;
			m_sumAmp = 0;
			m_freqNorm = 20;
		}
		virtual ~ModelBasedDesigner () {
			//if (m_peaks) delete [] m_peaks;
		}
		void analyse (T sr, T* data, int len, int msecPos, int N) {
			int startPos = (int) ((T) msecPos / 1000. * sr); // position in samples
			
			if (startPos < 1 || startPos > (len - N)) return;
			
			//std::cout << "start sample = " << startPos << std::endl;
			
			T* cbuf = new T[2 * N];
			memset (cbuf, 0, sizeof (T) * N * 2);
			
			T* win = new T[N];
			makeWindow<T> (win, N, .5, .5, 0); // hanning
			
			T* amps = new T[N];
			T* freq = new T[N];
			
			for (int i = startPos; i < startPos + N; ++i) {
				if (i >= len) break;
				cbuf[2 * (i - startPos)] = data[i] * win[i - startPos];
				cbuf[2 * (i - startPos) + 1] = 0;
			}
			
			fft (cbuf, N,  -1);
			ampFreqParabolic (cbuf, amps, freq, N, sr);
			
			std::vector<int> maxima;
			locmax2AmpFreq<T> (amps, freq, N / 2, maxima, 20.);
	
			//std::cout << "peaks = " << maxima.size () << std::endl;
			
			m_partials = maxima.size ();
			if (m_partials <= 0) return;
			
			if (m_peaks) delete [] m_peaks;
			m_peaks = new Peak<T>[m_partials];				
			
			T sumAmp = 0;
			for (int i = 0; i < m_partials; ++i) {
				m_peaks[i].freq = freq[maxima[i]];
				m_peaks[i].amp = amps[maxima[i]];
				sumAmp += amps[maxima[i]];
			}
			
			// FIXME: la normalizzazione non e' del tutto corretta poiche' assume
			//        che la nuova fondamentale sia il primo picco....
			m_freqNorm = m_peaks[0].freq; // before reordering
			if (m_freqNorm == 0) m_freqNorm += 20;
			
			sortSpectrum (m_peaks, m_partials);
			//
			//for (int i = 0; i < m_partials; ++i) {
			//	std::cout << m_peaks[i].freq << " Hz, " << m_peaks[i].amp << std::endl;
			//}
			//sumAmp *= 3; // FIXME
			
			delete [] cbuf;
			delete [] win;
			delete [] amps;
			delete [] freq;
		}
		virtual T design (T f0, T dev, int part, T& amp) {
			if (part >= m_partials || part < 1 || m_partials == 0) return 0;
			
			amp = m_peaks[part - 1].amp;
			return (m_peaks[part - 1].freq / m_freqNorm) * f0 * dev;
		}
		const char* type () const {
			return "ModelBasedDesigner";
		}
		
	private:
		T m_sumAmp;
		T m_freqNorm;
		Peak<T>* m_peaks;
		int m_partials;
	};
}

#endif	// MODELBASEDDESIGNER_H 

// EOF
