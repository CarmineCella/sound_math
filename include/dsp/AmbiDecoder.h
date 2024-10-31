// AmbiDecoder.h
// 

#ifndef AMBIDECODER_H
#define AMBIDECODER_H

#include "constants.h"

namespace soundmath {
	struct AmbiDecoder {
		virtual ~AmbiDecoder () {}
		virtual void op (int N, const float* ww, const float* xx, const float* yy, float* lr) = 0;
	};
	
	struct MidSide : public AmbiDecoder {
		virtual ~MidSide () {}
		void op (int N, const float* ww, const float* xx, const float* yy, float* lr)  {
			for	(int i = 0; i < N; ++i) {
				float mid = (ww[i] * SQRT2) + xx[i];
				lr[2 * i] = (mid + yy[i]) / SQRT2;
				lr[2 * i + 1] = (mid - yy[i]) / SQRT2;
			}
		}
	};
	
	struct Blumlein : public AmbiDecoder {
		virtual ~Blumlein () {}
		void op (int N, const float* ww, const float* xx, const float* yy, float* lr)  {
			for	(int i = 0; i < N; ++i) {
				lr[2 * i] = (xx[i] + yy[i]) / SQRT2;
				lr[2 * i + 1] = (xx[i] - yy[i]) / SQRT2;
			}
		}
	};
	
	struct UHJ : public AmbiDecoder {
		virtual ~UHJ () {}
		void op (int N, const float* ww, const float* xx, const float* yy, float* lr)  {
			for	(int i = 0; i < N; ++i) {
				float s = 0.9396926 * ww[i] + 0.1855740 * xx[i];
				float d = -1. * (-0.3420201 * ww[i] + 0.5098604 * xx[i]) + 0.6554516 * yy[i];
			
				lr[2 * i] = (d + d) / 2.;
				lr[2 * i + 1] = (s - d) / 2.;
			}
		}
	};
}

#endif	// AMBIDECODER_H 

// EOF
