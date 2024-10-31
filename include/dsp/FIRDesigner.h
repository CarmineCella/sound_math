// FIRDesigner.h
// 

#ifndef FIRDESIGNER_H
#define FIRDESIGNER_H

#include <cmath>

namespace soundmath {
	template <typename T>
	class FIRDesigner {
	public:
		void designLowPass (T* h, int N, T sr, T cutoff) {
			T s = 0;
			for (int i = 0; i < N; ++i) {
				T radians = 2. * M_PI * cutoff / sr;
				T t = ((T) i - (N / 2));
				
				if (t == 0) {
					h[i] = radians;
					continue;
				}
				
				h[i] = sin (radians * t) / (t);
				h[i] *= (0.54 - 0.46 * cos (2. * M_PI * (T) i / N));
				s += h[i];
			}
			
			for (int i = 0; i < N; ++i) {
				h[i] /= s;
			}   
		}
		void designHighPass (T* h, int N, T sr, T cutoff) {
			designLowPass  (h, N, sr, cutoff);
			invertSpectrum (h, N);
		}
	
		void designBandReject (T* h, int N, T sr, T cutoff, T bw) {
			designLowPass  (h, N, sr, cutoff);
			T* tmp = new T[N];
			designHighPass  (tmp, N, sr, cutoff + bw);
			
			for (int i = 0; i < N; ++i) {
				h[i] += tmp[i];
			}
			
			delete [] tmp;
		}
		void designBandPass (T* h, int N, T sr, T cutoff, T bw) {
			designBandReject (h, N, sr, cutoff, bw);
			invertSpectrum (h, N);
		}
	
	private:
		void invertSpectrum (T* h, int N) {
			for (int i = 0; i < N; ++i) {
				 h[i] = -h[i];
			}
			h[(N - 1) / 2] += 1.;
		}
		 
    };
}

#endif	// FIRDESIGNER_H 

// EOF
