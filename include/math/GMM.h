// GMM.h
// 

#ifndef GMM_H
#define GMM_H 

#include "Matrix.h"
#include "Choelsky.h"
#include <cmath>

namespace soundmath {
    
	template <typename T>
    //! Basic structure for GMM analysis
	struct PreGMM {
		static int mmstat;
		struct Mat_mm : Matrix <T> {
			Mat_mm () : Matrix <T> (mmstat, mmstat) {}
		};
		PreGMM (int mm) {
			mmstat = mm;
		}
	};
	template <typename T>
	int PreGMM<T>::mmstat = -1;
	
    
	template <typename T>
        //! Data clustering by means of Gaussian mixture models
		struct GMM : PreGMM <T> {
		int nn, kk, mm;
		Matrix<T> data, means, resp;
		std::vector<T> frac, lndets;
		std::vector<typename PreGMM <T>::Mat_mm> sig;
		T loglike;
		GMM (Matrix<T>& ddata, Matrix<T>& mmeans) : PreGMM <T> (ddata.cols ()),
				nn (ddata.rows ()), kk (mmeans.rows ()),
				mm (PreGMM <T>::mmstat), data (ddata), means (mmeans),
				resp (nn, kk), frac (kk), lndets (kk), sig (kk) {
			int i, j, k;
			for (k  =0; k < kk; k++) {
				frac[k] = 1. / kk;
				for (i = 0; i < mm; i++) {
					for (j = 0;j < mm; j++) sig[k][i][j] = 0.;
					sig[k][i][i] = 1.0e-10;
				}
			}
			estep ();
			mstep ();
		}
		T estep () {
			int k, m, n;
			T tmp, sum, max, oldloglike;
			std::vector<T> u (mm), v (mm);
			oldloglike = loglike;
			for (k = 0; k < kk; k++) {
				Cholesky<T> choltmp (sig[k]);
				lndets[k] = choltmp.logdet ();
				for (n = 0; n < nn; n++) {
					for (m = 0; m < mm; m++) u[m] = data[n][m] - means[k][m];
					choltmp.elsolve (u, v);
					for (sum = 0., m = 0; m < mm; m++) sum += squared (v[m]);
					resp[n][k] = -0.5 * (sum + lndets[k]) + log (frac[k]);
				}
			}
			loglike = 0;
			for (n = 0; n < nn; n++) {
				max = -99.9e99;
				for (k = 0; k < kk; k++) if (resp[n][k] > max) max = resp[n][k];
				for (sum = 0., k = 0; k < kk; k++) sum += std::exp (resp[n][k] - max);
				tmp = max + log(sum);
				for (k = 0; k < kk; k++) resp[n][k] = std::exp (resp[n][k] - tmp);
				loglike +=tmp;
			}
			return loglike - oldloglike;
		}
		void mstep () {
			int j, n, k, m;
			T wgt, sum;
			for (k = 0; k < kk; k++) {
				wgt = 0.;
				for (n = 0; n < nn; n++) wgt += resp[n][k];
				frac[k] = wgt / nn;
				for (m = 0; m < mm; m++) {
					for (sum = 0., n = 0; n < nn; n++) sum += resp[n][k] * data[n][m];
					means[k][m] = sum / wgt;
					for (j = 0;j < mm; j++) {
						for (sum = 0., n = 0; n < nn; n++) {
							sum += resp[n][k]*
								   (data[n][m] - means[k][m]) * (data[n][j] - means[k][j]);
						}
						sig[k][m][j] = sum / wgt;
					}
				}
			}
		}
	};
}

#endif	// GMM_H 

// EOF
