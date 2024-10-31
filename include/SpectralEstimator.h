// SpectralEstimator.h
// 

#ifndef SPECTRALESTIMATOR_H
#define SPECTRALESTIMATOR_H 

#include "AbstractEstimator.h"
#include "algorithms.h"
#include "FFT.h"
#include <fstream>
#include <iostream>
#include <cstring>

#define DELTA_THRESHOLD .2
#define AMP_RATIO .25
#define MAX_HARMONIC 6

namespace soundmath {
	//! Fundamental frequency estimator based on spectrum analysis
	template <typename T>
	class SpectralEstimator : public AbstractEstimator {
	private:
		SpectralEstimator& operator= (SpectralEstimator&);
		SpectralEstimator (const SpectralEstimator&);
	public: 
		SpectralEstimator (T sr, int winSize, int fftPad, T threshold) :
			m_sr (sr), m_winSize (winSize), m_threshold (threshold) {
			m_follower = new EnvelopeFollower<T> ();
			int WnPow2 = (int) pow (2., ceil (std::log (winSize) / std::log (2.)));
			m_N  = WnPow2 << fftPad;
			m_NN = m_N << 1;
			m_N2 = m_N >> 1;
			m_fft = createFFT<T> (m_N);

			m_wksp = new T[m_NN];
			m_window = new T[winSize];
			makeWindow<T> (m_window, m_winSize, .5, .5, 0);
			m_amp = new T[m_N];
			m_freq = new T[m_N];
			m_phi = new T[m_N];
			memset (m_phi, 0, sizeof (T) * m_N);		
			m_peaks = new int[m_N];
			m_peakAmp = new T[m_N];
			m_peakFreq = new T[m_N];
		}
		virtual ~SpectralEstimator () {
			delete [] m_window;
			delete [] m_peakAmp;
			delete [] m_peakFreq;
			delete [] m_peaks;
			delete [] m_phi;
			delete [] m_freq;
			delete [] m_amp;
			delete [] m_wksp;
			delete m_fft;
			delete m_follower;
		}
		T process (const T* in, T& amp) {
			static int cc = 0;
			T f0 = -1; // default to no f0
		
			amp = m_follower->process (in, m_winSize);			
			// std::cout << "amp: " << amp << std::endl;
			if (amp > m_threshold) {
				memset (m_wksp, 0, sizeof (T) * m_NN);
				for (int i = 0; i < m_winSize; ++i) {
					m_wksp[2 * i] = in[i] * m_window[i];
				}
				m_fft->forward (m_wksp);
	
				ampFreqPhaseDiff (m_wksp, m_amp, m_freq, m_phi, m_N2, m_winSize, m_sr);
				int num = locmax (m_amp, m_N2, m_peaks);
				//std::cout << "\tpeaks: " << num << std::endl;
				if (num != 0) { // peaks found
					memset (m_peakAmp, 0, sizeof (T) * m_N);
					memset (m_peakFreq, 0, sizeof (T) * m_N);
					for (int i = 0; i < num; ++i) {
						m_peakAmp[i] = m_amp[m_peaks[i]];
						m_peakFreq[i] = m_freq[m_peaks[i]];
					}
					int maxPos = 0;
					maximum <T> (m_peakAmp, num, maxPos);
					T meanAmp = mean<T> (m_peakAmp, num);
					T centroid = speccentr<T> (m_peakAmp, m_peakFreq, num);

					//std::cout << "\t\tmax peak is: " << maxPos << " = " 
							  //<< m_freq[m_peaks[maxPos]] << " " << m_peakFreq[maxPos] 
							  //<< " centroid: " << centroid << std::endl;

					f0 = m_freq[m_peaks[maxPos]]; // highest peak is also f0

					// if (cc == 75) {

					// 	double tmp[m_N2];
					// 	memset (tmp, 0, sizeof (double) * m_N2);
					// 	for (int i = 0; i < num; ++i) {
					// 		tmp[m_peaks[i]] = m_amp[m_peaks[i]];
					// 	}
					// 	std::ofstream outamp ("amp.raw", std::ios::binary);
					// 	std::ofstream outfreq ("freq.raw", std::ios::binary);
					// 	std::ofstream outpeaks ("peaks.raw", std::ios::binary);
						
					// 	for (int i = 0; i < m_N2; ++i) {
					// 		outamp.write ((char*) &m_amp[i], sizeof (T));
					// 		outfreq.write ((char*) &m_freq[i], sizeof (T));
					// 		outpeaks.write ((char*) &tmp[i], sizeof (T));
					// 	}

					// }

					for (int i = 0; i < maxPos; ++i) {
						T ratio = m_freq[m_peaks[maxPos]] / m_freq[m_peaks[i]];
						T delta = ratio - ((int) ratio);
							
						//std::cout << "\t\t\t" << i << " -> ratio: " << ratio
								  ////<< ", amp: " << m_amp[m_peaks[i]] << ", freq: " 
								  //<<  m_freq[m_peaks[i]] << " - mean amp: " << meanAmp << std::endl;
						//std::cout << "\t\t\tdelta: " << delta << " hpos: " 
								  //<< (T) (maxPos - i - 1) + DELTA_THRESHOLD << std::endl << std::endl;
						if (
							// ratio < (T) i + DELTA_THRESHOLD && ratio > (T) i < DELTA_THRESHOLD
							delta < DELTA_THRESHOLD
							&& m_amp[m_peaks[i]] > m_peakAmp[maxPos] * AMP_RATIO
							&& i < MAX_HARMONIC
							&& m_freq[m_peaks[i]] < centroid
							) {
							f0 = m_freq[m_peaks[i]]; // got f0
							//std::cout << "\t\t\tf0 is peak " << i << " = " << f0 << std::endl;
							break;
						}								 	   	
					}

				}
				
			}	
			++cc;
			return f0;
		}
	private:
		EnvelopeFollower<T>* m_follower;
		AbstractFFT<T>* m_fft;
		T* m_wksp;
		T* m_amp;
		T* m_freq;
		T* m_phi;	
		int* m_peaks;
		T* m_peakAmp;
		T* m_peakFreq;
		T* m_window;
		T m_sr;
		int m_winSize;
		int m_N;
		int m_NN;
		int m_N2;
		T m_threshold;
	};
}


#endif	// SPECTRALESTIMATOR_H 

// EOF
