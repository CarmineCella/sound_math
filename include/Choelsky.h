// Choelsky.h
// 

#ifndef CHOELSKY_H
#define CHOELSKY_H

#include <vector>
#include <stdexcept>
#include <cmath>

namespace soundmath {
    
	template <typename T>
	//! Choelsky decomposition
    struct Cholesky {
		int n;
		Matrix<T> el;
		Cholesky (Matrix<T> &a) : n (a.rows ()), el (a) {
			int i, j, k;
			std::vector<T> tmp;
			T sum;
			if ((int) el.cols () != n) throw std::runtime_error ("need square matrix");
			for (i=0;i<n;i++) {
				for (j = i;j < n; j++) {
					for (sum = el[i][j], k = i - 1; k >= 0; k--) sum -= el[i][k] * el[j][k];
					if (i == j) {
						if (sum <= 0.0)
							throw std::runtime_error ("Cholesky failed");
						el[i][i] = sqrt(sum);
					} else el[j][i] = sum / el[i][i];
				}
			}
			for (i = 0; i < n; i++) for (j = 0; j < i; j++) el[j][i] = 0.;
		}
		void solve (const std::vector<T> &b, std::vector<T> &x) {
			int i, k;
			T sum;
			if (b.size () != n || x.size () != n) throw std::runtime_error ("bad lengths in Cholesky");
			for (i=0;i<n;i++) {
				for (sum = b[i], k = i - 1; k >= 0; k--) sum -= el[i][k] * x[k];
				x[i] = sum / el[i][i];
			}
			for (i = n - 1; i >= 0; i--) {
				for (sum = x[i], k = i + 1; k < n; k++) sum -= el[k][i] *x[k];
				x[i] = sum / el[i][i];
			}
		}
		void elmult (const std::vector<T> &y, std::vector<T> &b) {
			int i, j;
			if (b.size () != n || y.size () != n) throw std::runtime_error ("bad lengths");
			for (i = 0; i < n; i++) {
				b[i] = 0.;
				for (j = 0;j <= i; j++) b[i] += el[i][j] * y[j];
			}
		}
		void elsolve (const std::vector<T> &b, std::vector<T> &y) {
			int i, j;
			T sum;
			if ((int) b.size () != n || (int) y.size () != n) throw std::runtime_error ("bad lengths");
			for (i = 0; i < n; i++) {
				for (sum = b[i], j = 0; j < i; j++) sum -= el[i][j] * y[j];
				y[i] = sum / el[i][i];
			}
		}
		void inverse (Matrix<T>& ainv) {
			int i, j, k;
			T sum;
			ainv.resize (n, n);
			for (i = 0; i < n; i++) for (j = 0; j <= i; j++) {
					sum = (i == j ? 1. : 0.);
					for (k = i - 1; k >= j; k--) sum -= el[i][k] * ainv[j][k];
					ainv[j][i]= sum / el[i][i];
				}
			for (i = n - 1; i >= 0; i--) for (j = 0; j <= i;j++) {
					sum = (i < j ? 0. : ainv[j][i]);
					for (k = i + 1; k < n; k++) sum -= el[k][i] * ainv[j][k];
					ainv[i][j] = ainv[j][i] = sum / el[i][i];
				}
		}
		T logdet () {
			T sum = 0.;
			for (int i = 0; i < n; i++) sum += log (el[i][i]);
			return 2. * sum;
		}
	};
}

#endif	// CHOELSKY_H 

// EOF
