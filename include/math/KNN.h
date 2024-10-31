// KNN.h
// 

#ifndef KNN_H
#define KNN_H 

#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>

#include <iostream>

namespace soundmath {
	template <typename T>
	//! Basic structure to describe data observations
	struct Observation {
		Observation () : classlabel (0) {}
		
		std::vector<T> attributes;
		T classlabel;
	};
	
	template <typename T>
	struct Item {
		Item () : distance (0), classlabel (0) {}
		T distance;
		T classlabel;
	};
	
	template <typename T> //, template <typename X> class DistanceType>
	//! Data clustering by means of K-nearest neighbours
	class KNN {
	private:
		KNN& operator= (KNN&);
		KNN (const KNN&);
	public: 
		KNN (int K, int features) {
			m_K = K;
			m_features = features;
			
			m_freq = new int[m_K];
			memset (m_freq, 0, sizeof (int) * m_K);
			
			m_knn = new Item<T>[m_K];
			memset (m_knn, 0, sizeof (Item<T>) * m_K);
		}
		virtual ~KNN () {
			delete [] m_freq;
			delete [] m_knn;
		}
		
		unsigned int addObservation (Observation<T>* v) {
			if (v->attributes.size () != m_features) {
				std::cout << v->attributes.size () << " " << m_features << std::endl;
				throw std::runtime_error ("invalid number of features for the given observation");
				
			}
			m_trainingSet.push_back (v);
			return m_trainingSet.size ();
		}
		unsigned int samples () const {
			return m_trainingSet.size ();
		}
		T classify (const Observation<T>& v) {
			T dd = 0;
			int maxn = 0;
			double mfreqC = 0;
			memset (m_freq, 0, sizeof (int) * m_K);
			
			for (int i = 0; i < m_K; ++i) m_knn[i].distance = std::numeric_limits<T>::max ();
			for (unsigned int i = 0; i < m_trainingSet.size (); ++i) {
				dd = distance (&m_trainingSet[i]->attributes[0], &v.attributes[0], m_features);
				maxn = max (m_knn);
				if (dd < m_knn[maxn].distance) {
					m_knn[maxn].distance = dd;
					m_knn[maxn].classlabel = m_trainingSet[i]->classlabel;
				}
			}
			for (int i = 0; i < m_K; ++i) m_freq[i] = 1;
			
			for (int i = 0; i < m_K; ++i) {
				for (int j = 0; j < m_K; ++j) {
					if ((i != j) && (m_knn[i].classlabel == m_knn[j].classlabel)) {
						m_freq[i] += 1;
					}
				}
			}
			
			int mfreq = 1;
			mfreqC = m_knn[0].classlabel;
			
			for (int i = 0; i < m_K; ++i) {
				if (m_freq[i] > mfreq)  {
					mfreq = m_freq[i];
					mfreqC = m_knn[i].classlabel;
				}
			}
			return mfreqC;
		}
	
		Item<T>* m_knn;
	protected:
		std::vector<Observation<T>*> m_trainingSet;
		int m_K;
		int* m_freq;
		int m_features;
		
	private:
		T distance (const T* v1, const T* v2, int size) {
			T d = 0.0;
			T tem = 0.0;
			for (int i = 0; i < size; ++i) {
				tem += (v1[i] - v2[i]) * (v1[i] - v2[i]);
			}
			d = std::sqrt (tem);
			return d;
		}
		int max (Item<T>* knn)	{
			int maxNo = 0;
			if (m_K > 1){
				for (int i = 1; i < m_K; ++i) {
					if (knn[i].distance > knn[maxNo].distance) maxNo = i;
				}
			}
			return maxNo;
		}    
	};
}

#endif	// KNN_H 

// EOF
