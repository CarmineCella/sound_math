// OnsetDetector.h
// 

#ifndef ONSETDETECTOR_H
#define ONSETDETECTOR_H 

#include "OnePole.h"

namespace soundmath {
	//! An onset detector based on energy
	template <typename T>
	class OnsetDetector {
	private:
		OnsetDetector& operator= (OnsetDetector&);
		OnsetDetector (const OnsetDetector&);
	public: 
		OnsetDetector (int sr){
			m_sr = sr;
			m_threshold = 0.01;
			m_timeGateInSamples = 50 * (((T) m_sr) / 1000);
			m_timeCount = 0;
			m_opol =  new OnePole<T> (m_sr, 80, 1);
			m_count = 0;
			m_ringBuffer = NULL;
			setMedianSize (3);
		}
		
		virtual ~OnsetDetector () {
			if (m_ringBuffer)	delete[] m_ringBuffer;
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
		T process (T* buf, T* tmp, int size) {
			// copy buffer
			memcpy (tmp, buf, size * sizeof (T));
	
	
			// square
			for (int i = 0; i < size; ++i){
				tmp[i] *= tmp[i]; 
			}
			
			// onepole
			m_opol->process (tmp, tmp, size);
			
			// snapshot
			T v = tmp[0];
			//for (int i = 0; i < size; ++i) {
			//    v += tmp[i];
			//}
			//v /= size;
			
			// sqrroot
			v = sqrtf (v);
			
			// median
			m_ringBuffer[m_count] = v;
			m_count = (m_count + 1) % m_medianSize;
			T mean = 0;
			for (int i = 0; i < m_medianSize; ++i){
				mean += m_ringBuffer[i];
			}
			mean /= m_medianSize;
			
			// timegate
			m_timeCount += size;
			if (m_timeCount < m_timeGateInSamples) return 0; // doing timegate after ring buffer to avoid discontinuities
			
			// diff median
			T diff = v - mean;
			
			// threshold
			if (diff < m_threshold) { // no onset detected
				diff = 0;
			} else	{	// onset detected
				m_timeCount = 0;
			}
			
			return diff;
		}
		
		void setMedianSize (int size) {
			if (m_ringBuffer) {
				delete m_ringBuffer;
			}
			
			m_ringBuffer = new T[size];
			
			for(int i=0; i<size; ++i){
				m_ringBuffer[i] = 0;
			}
			m_medianSize = size;
		}
	
	private:
		T m_threshold;    
		OnePole <T>* m_opol;
		T* m_ringBuffer;
		int m_count;
		int m_medianSize;
		int m_timeCount; //timegate
		T m_sr;
		int m_timeGateInSamples;
	};
}

#endif	// ONSETDETECTOR_H 

// EOF
