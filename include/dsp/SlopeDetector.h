// SlopeDetector.h
// 

#ifndef SLOPEDETECTOR_H
#define SLOPEDETECTOR_H 

#include "OnePole.h"
//#include "EnvelopeFollowerBuffer.h"

//#define USE_LINEAR_REGRESSION

namespace soundmath {
	//! An onset detector based on slope calculation
	template <typename T>
	class SlopeDetector {
	private:
		SlopeDetector& operator= (SlopeDetector&);
		SlopeDetector (const SlopeDetector&);
	public: 
		SlopeDetector (int sr, int maxBufferSize){
			m_sr = sr;
			m_threshold = 0.01;
			m_timeGateInSamples = 50 * (((T) m_sr) / 1000);
			m_timeCount = 0;
			m_rms = new T[maxBufferSize];
			m_opol =  new OnePole<T> (m_sr, 80, 1);
			//m_ef = new EnvelopeFollowerBuffer<T> (.9999);
		}
		
		virtual ~SlopeDetector () {
			delete[] m_rms;
		}
		
		void setTimeGateMs (T ms) {
			m_timeGateInSamples = ms * (((T) m_sr) / 1000);
		}
		
		void setThreshold (T t) {
			m_threshold = t;
		}
		void setSmoothingHz (T hz) {
			m_opol->setup (hz);
		}
		T process (T* buf, int size, T& rms) {
			//static WavOutFile out ("debug.wav", 44100, 16, 1);
			//T min = 0, max = 0;
			//m_ef->process (buf, m_rms, size, max, min);
			//out.write (m_rms, size);
			
			//rms
			rms = 0;
			//memcpy (m_rms, buf, sizeof (T) * size);
	
			// onepole
			m_opol->process (buf, m_rms, size);
				
	#ifdef USE_LINEAR_REGRESSION
			//LINEAR REGRESSION
			T sx = 0;
			T sy = 0;
			T stt = 0;
			T sts = 0;
			
			for (int i = 0; i < size; ++i){
				sx += i;
				sy += m_rms[i];
			}
			
			for (int i = 0; i < size; ++i){
				T t = i - (sx / size);
				stt += (t * t);
				sts += t * m_rms[i];
			}
			
			T slope = (sts / stt) ;
			T intercept = (sy - sx * slope) / size;
	#else
			// DIFF IN RMS
			T* prms = &m_rms[0];
			int n = size;
			T min = m_rms[0];
			T max = m_rms[0];        
			while (n--) {
				*prms = sqrt (*prms * *prms);
				
				if (*prms > rms) {
					rms = *prms;
				}
				if (*prms < min) {
					min = *prms;
				}
				if (*prms > max) {
					max = *prms;
				}        
				prms++;
			}
			T slope = max - min;        
	#endif
	
			// timegate
			m_timeCount += size;
			if (m_timeCount < m_timeGateInSamples) {
				return 0;
			}
			
			// threshold
			if (slope < m_threshold) { // no onset detected
				slope = 0;
			} else	{	// onset detected
				m_timeCount = 0;   		
			}
	
			return slope;
		}
		
	
	private:
		T* m_rms;
		T m_threshold;    
		int m_timeCount; //timegate
		T m_sr;
		int m_timeGateInSamples;
		OnePole <T>* m_opol;
	//    EnvelopeFollowerBuffer<T>* m_ef;
	};
}

#endif	// SLOPEDETECTOR_H 

// EOF
