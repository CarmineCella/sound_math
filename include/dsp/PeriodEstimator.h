// PeriodEstimator.h
//

#ifndef PERIODESTIMATOR_H
#define PERIODESTIMATOR_H

#include "AbstractEstimator.h"
#include "algorithms.h"
#include "VariState.h"
#include "ZeroCrossings.h"
#include "EnvelopeFollower.h"
#include "DCBlocker.h"
#include <vector>
#include <algorithm>

#define PE_MIN_LEVELS 1
#define PE_MIN_SPLIT 2

namespace soundmath {
	//! Fundamental frequency estimator based on period estimation by means of subband filtering
	template <typename T>
	class PeriodEstimator : public AbstractEstimator {
	private:
		PeriodEstimator& operator= (PeriodEstimator&);
		PeriodEstimator (const PeriodEstimator&);
	public:
		PeriodEstimator (T sr, int N, T min = 77.782,
						 T max = 1396.913, int levels = 6, int split = 10,
						 T threshold = .01, T smoothingRate = 2000) :
			m_sr (sr), m_N (N),
			m_min (min), m_max (max), m_levels (levels), m_split (split), 
			m_threshold (threshold), m_num (0) {
			if (m_min > m_max) {
				m_min = m_max;
			}
			if (m_levels < PE_MIN_LEVELS) {
				m_levels = PE_MIN_LEVELS;
			}
			if (m_split < PE_MIN_SPLIT) {
				m_split = PE_MIN_SPLIT;
			}
	
			T current = min;
			while (current < max) {
				m_bins.push_back (current *= 1.4142); // half oct.
				m_num++;
			}
	   
			m_filters = new VariState<T>*[m_num * m_levels];
			m_zeroes = new ZeroCrossings<T>*[m_num];
			m_envFollowers = new EnvelopeFollower<T>*[m_num];
			m_filtered = new T*[m_num];
			for (int i = 0; i < m_num; i++) {
				m_filtered[i] = new T[N];
				memset (m_filtered[i], 0, sizeof (T) * m_N);

				for (int k = 0; k < m_levels; k++) {
					m_filters[(i * m_levels) + k] = new VariState<T> (sr, m_bins[i], LOPASS);
				}    
				m_zeroes[i] = new ZeroCrossings<T> (sr);
				m_envFollowers[i] = new EnvelopeFollower<T> ();
			}

			m_block = new DCBlocker<T> (sr);
			m_final = new VariState<T> (sr, smoothingRate, LOPASS);
	
			m_periods = new T[m_num];
			memset (m_periods, 0, sizeof (T) * m_num);
			m_followed = new T[m_num];
			memset (m_followed, 0, sizeof (T) * m_num);
	
			m_tmp = new T[m_N];
			memset (m_tmp, 0, sizeof (T) * m_N);

			m_counter = 0;
			m_period = 0;	
		}
		virtual ~PeriodEstimator () {
			delete [] m_filters;
			delete [] m_zeroes;
			delete m_final;
			delete [] m_envFollowers;
			delete m_block;
			for (int i = 0; i < m_num; i++) {
				delete [] m_filtered[i];
			}
			delete [] m_filtered;

			delete [] m_periods;
			delete [] m_followed;
		}
		void threshold (T t) { 
			m_threshold = t; 
		}

		void reset () {
			m_block->reset ();
			for (int i = 0; i < m_num; i++) {
				for (int k = 0; k < m_levels; k++) {
					m_filters[(i * m_levels) + k]->reset ();
				}  
				m_zeroes[i]->reset ();
				m_envFollowers[i]->reset (); 
			}        	    
			m_final->reset ();
			memset (m_filtered, 0, sizeof (T) * m_num);
			memset (m_periods, 0, sizeof (T) * m_num);
			memset (m_followed, 0, sizeof (T) * m_num);
		}

		void split (T s) { 
			if (s < PE_MIN_SPLIT) {
				s = PE_MIN_SPLIT; 
			}
			m_split = s; 
		}

		void filterQuality (T q) {
			for (int i = 0; i < m_num; i++) {
				for (int k = 0; k < m_levels; k++) {
					m_filters[(i * m_levels) + k]->quality (q);
				}
			} 
		}

		void scaleBinsTo (T m) {
			for (int i = 0; i < m_num; i++) {
				m_bins[i] *= m;
				for (int k = 0; k < m_levels; k++) {
					m_filters[(i * m_levels) + k]->cutoff (m_bins[i]);
				}   
			}   
		}

		void feedback (T f) {
			for (int i = 0; i < m_num; i++) {
				m_zeroes[i]->feedback (f);
			}  	
		}

		void smoothingRate (T f) {
			m_final->cutoff (f);
		}

		void smoothingQuality (T f) {
			m_final->quality (f);
		}

		T min () const { return m_min; }
		T max () const { return m_max; }

		T process (const T* in, T& max) {
			// m_block->process (in, m_tmp, m_N);
			// std::copy (in, in + m_N, m_tmp);
			for (int i = 0; i < m_num; i++) {
				std::copy (in, in + m_N, m_filtered[i]); // dc
				for (int k = 0; k < m_levels; k++) {
					m_filters[(i * m_levels) + k]->process (m_filtered[i], m_tmp, m_N);
					std::copy (m_tmp, m_tmp + m_N, m_filtered[i]); // dc;
				}
				m_followed[i] = m_envFollowers[i]->process (m_filtered[i], m_N);
				m_periods[i] = m_zeroes[i]->process (m_filtered[i], m_N);
			}
			int pos = 0;
			max = maximum<T> (m_followed, m_num, pos);
			//std::cout << "max : " << max << std::endl;
			int filt = 0;
			if (max < m_threshold) {
				m_val = -1;
			}
			else {
				T maxSplit = max / m_split;
				m_val = -1;
				while (filt < m_num) {
					if (m_followed[filt] >= maxSplit) {
						if (m_periods[filt] <= m_bins[filt] && m_periods[filt] >= m_min) {
							m_val = m_periods[filt]; //m_final->process (m_periods[filt]);
							break;
						}
						else if (filt < m_num - 1 && m_periods[filt] <= m_bins[filt + 1]
								 && m_periods[filt] >= m_min) {
							m_val = m_periods[filt]; // m_final->process (m_periods[filt]);
							break;
						}
						else {
							break;
						}
					}
					filt++;
				}
			}

			return m_val;
		}

		virtual T get () const { return m_val; }
	private:
		T m_sr;
		int m_N;
		T m_min, m_max;
		int m_levels;
		T m_split;
		T m_threshold;
		long m_num;
		VariState<T>** m_filters;
		ZeroCrossings<T>** m_zeroes;
		EnvelopeFollower<T>** m_envFollowers;
		VariState<T>* m_final;
		DCBlocker<T>* m_block;
		T** m_filtered;
		T* m_followed;
		T* m_periods;
		std::vector<T> m_bins;
		long m_counter, m_period;
		T m_val;
		T* m_tmp;
	};
}
#endif	// PERIODESTIMATOR_H

// EOF


