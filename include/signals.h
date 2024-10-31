// signals.h
// 


#ifndef SIGNALS_H
#define SIGNALS_H

#include "constants.h"
#include <cmath>

namespace soundmath {
	// ------------------------------------------------------------------//
	
	template <typename T>
	inline T logAmplitude (const T& input) {
		// sound pressure level (minimum 32 bits)
		if (input < P0) return 0;
		T f = 20. * std::log10 ((T) (input / P0));
		return f;
	}
	
	template <typename T>
		void conv (T* X, T* Y, T* Z, int lenx, int leny, T scale) {
		T *zptr, s, *xp, *yp;
		int i, n, n_lo, n_hi;
		int lenz = lenx + leny - 1;
		
		zptr = Z;
		for (i = 0; i < lenz; i++) {
			s=0.0;
			n_lo= 0 > (i - leny + 1) ? 0 : i - leny + 1;
			n_hi = lenx - 1 < i ? lenx - 1 : i;
			xp = X + n_lo;
			yp = Y + i-n_lo;
			for (n = n_lo; n <= n_hi; n++) {
				s += *xp * *yp;
				xp++;
				yp--;
			}
			*zptr = s * scale;
			zptr++;
		}
	}
	template <typename T>
		T acf (const T* in, T* out, int len) {	
		int len2 = len >> 1;
		T scale = 1. / len;
		
		for (int i = 0; i < len2; ++i) {
			T acfVal = 0.;
	
			for (int j = 0; j < len2; ++j) {
				acfVal += in[i + j] * in[j];
			}
			out[i] = acfVal * scale;
		}
	}
		
	template <typename T>
		T amdf (const T* in, T* out, int len) {	
		int len2 = len >> 1;
		T scale = 1. / len;
		
		for (int i = 0; i < len2; ++i) {
			T amdfVal = 0.;
	
			for (int j = 0; j < len2; ++j) {
				amdfVal += fabs (in[j] - in[i + j]);
			}
			out[i] = amdfVal * scale;
		}
	}	
		
	template <typename T>
		void gdelay (const T* f, T* phase, T* t, int size) {
		/// unwrap phase data 
		for (int n = 1; n < size; ++n) { 
			while (phase[n] - phase[n - 1] <= -180.) {
				phase[n] = phase[n] + (T) 360.; 
			}
			while (phase[n] - phase[n - 1] >= 180.) {
				phase[n] = phase[n] - (T) 360.; 
			}
		}
		 
		// calculate group delay 
		for (int n = 1; n < size - 1; ++n) {
			t[n] = A1 * (((phase[n] - phase[n - 1]) / (f[n] - f[n - 1]))
						 + ((phase[n + 1] - phase[n]) / (f[n + 1] - f[n]))); 
		}
		t[0] =  A2 * (((phase[1] - phase[0]) / (f[1] - f[0]))); 
		t[size - 1] = A2 * (((phase[size - 1] - phase[size - 2]) / (f[size - 1] - f[size - 2])));
	}
		
}
#endif	// SIGNALS_H 

// EOF
