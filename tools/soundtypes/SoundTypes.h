// SoundTypes.h
//

#ifndef SOUNDTYPES_H
#define SOUNDTYPES_H

#include "WavFile.h"
#include "Matrix.h"
#include "Vector.h"
#include "DynamicMatrix.h"
#include "algorithms.h"
#include "signals.h"
#include "FFT.h"
#include "MFCC.h"
#include "features.h"
#include "ModelingSpace.h"
#include "GMM.h"
#include "BlockVocoder.h"

#include <vector>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <map>

// NB: differences with sparkle:
//	   output hop is variable thus corrupting audio quality during stretches
//	   there are no phase transformations
//     thres is no traditional xsynth
//     there is no fft padding

// BRUNO - REMATCH
//1 - calcoli media Mi e stddev Si diì ogni descrittore i nel tuo database
//2 - registri un pezzettino del segnale che ti arriva in input, in modo che sia il piu significativo possibile
//3 - di questo segnale, calcoli media M e stddev S
//4 - fissi un coefficiente di mapping C (che puoi cambiare in tempo reale) che vari ad esempio tra 0.1 e 2 (tanto per fare un esempio)
//5 - al momento del match, scali ogni descrittore di input Xi  tra (Mi - C*Si) e (Mi + C*Si) in modo che piu o meno il segnale sia scalato in [-1,1]
//6 - poi fai il match con i descrittori, anch'essi scalati in maniera analoga
//
//

// PROB MERGE MATLAB
// I re-arrange the target and source transition matrices so that the same rows/columns
// correspond to the closest types/centroids. If there are no matches for a type,
// a new row/column of zeros is inserted to the opposite matrix. Once both matrix
// have the same size and have been re-arranged, the weighted sum is taken (M=a*T+(1-a)*S)


extern const int MAXBANDS = 48;
extern const int MAX_REPETITIONS = 3;
extern const double DENORM = .0000001;

namespace soundmath {
struct SoundPacket  {
	SoundPacket () {
		typeNumber = 0;
		instances = 0;
		compliant = false;
	}

	int typeNumber;
	int instances;
	bool compliant;
};

template <typename T>
class SoundTypes {
private:
	SoundTypes& operator= (SoundTypes&);
	SoundTypes (const SoundTypes&);
public:
	SoundTypes (const ModelingSpace<T>* space) {

		p = (ModelingSpace<T>*) space;

		workspace = new T[2 * p->N];
		workspace2 = new T[2 * p->N];
		freq = new T[p->N];
		amp = new T[p->N];
		oldAmp = new T[p->N];
		memset (oldAmp, 0, sizeof (T) * p->N);
		phi = new T[p->N];
		memset (phi, 0, sizeof (T) * p->N);

		buffer = new T[p->N];
		window = new T[p->N];
		createWindow (window);

		result = new T[p->N];

		ft = soundmath::createFFT<T> (p->N);
		mfcc = new soundmath::MFCC<T> (p->sr, MAXBANDS, p->N);

	}
	virtual ~SoundTypes () {
		delete [] result;
		delete [] window;
		delete [] workspace;
		delete [] workspace2;
		delete [] buffer;
		delete [] amp;
		delete [] oldAmp;
		delete [] freq;
		delete [] phi;

		delete ft;
	}
    ModelingSpace<T>* modelingSpace () const {
        return p;
    }
	void analyse (const T* data, int maxBin, int totSamp, int nchnls) {
		int pointer = 0;
		T tmp = 0, ctrd = 0, sprd = 0, f0 = 0;
		T winEnergy = 0;
		for (int i = 0; i < p->N; ++i) {
			winEnergy += (window[i] * window[i]);
		}

		do {
			//memset(workspace, 0, sizeof (T) * 2 * p->N);
			//memset(buffer, 0, sizeof (T) * p->N);

			//int r = p->N > totSamp - pointer ? totSamp - pointer : p->N;
			int r = pointer + p->N > (totSamp) ? pointer + p->N - (totSamp) : 0;
			pointer -= r;
			for (int i = 0; i < p->N; ++i) {
				workspace[2 * i] = 0.;
				workspace[2 * i + 1] = 0.;
				for (int j = 0; j < nchnls; ++j) {
					int pp = nchnls * (i + pointer) + j;
					// int pp = (nchnls * pointer) + (nchnls * i + j);

					workspace[2 * i] += (T)  data[pp];;

				}

				workspace[2 * i] *= window[i];
				workspace[2 * i] /= nchnls;
				buffer[i] = workspace[2 * i];
			}

	#ifdef USE_SLOW_FFT
			soundmath::fft<T> (workspace, p->N, -1);
	#else
			ft->forward(workspace);
	#endif

			soundmath::ampFreqBins<T> (workspace, amp, freq, p->N, p->sr);

			//getchar ();
			switch (p->scale) {
				case Features::LINEAR:
					// nothing to do for linear
					break;
				case Features::LOG:
					for (int i = 0; i < p->N; ++i) {
						amp[i] = soundmath::logAmplitude(amp[i]);
					}
					break;
				case Features::POWER:
					for (int i = 0; i < p->N; ++i) {
						T a = amp[i];
						amp[i] = a * a;
					}
					break;
			}

			if (p->switches[Features::SPECTRAL_CENTROID]) {
				ctrd = soundmath::speccentr<T> (amp, freq, maxBin);
				p->matrix[Features::SPECTRAL_CENTROID].push_back(ctrd);
			}
			if (p->switches[Features::SPECTRAL_SPREAD]) {
				sprd = soundmath::specspread<T> (amp, freq, maxBin, ctrd);
				p->matrix[Features::SPECTRAL_SPREAD].push_back(sprd);
			}
			if (p->switches[Features::SPECTRAL_SKEWNESS]) {
				tmp = soundmath::specskew<T> (amp, freq, maxBin, ctrd, sprd);
				p->matrix[Features::SPECTRAL_SKEWNESS].push_back(tmp);
			}
			if (p->switches[Features::SPECTRAL_KURTOSIS]) {
				tmp = soundmath::speckurt<T> (amp, freq, maxBin, ctrd, sprd);
				p->matrix[Features::SPECTRAL_KURTOSIS].push_back(tmp);
			}
			if (p->switches[Features::SPECTRAL_FLUX]) {
				tmp = soundmath::specflux<T> (amp, oldAmp, maxBin);
				p->matrix[Features::SPECTRAL_FLUX].push_back(tmp);
			}
			if (p->switches[Features::SPECTRAL_IRREGULARITY]) {
				tmp = soundmath::specirr<T> (amp, maxBin);
				p->matrix[Features::SPECTRAL_IRREGULARITY].push_back(tmp);
			}
			if (p->switches[Features::SPECTRAL_DECREASE]) {
				tmp = soundmath::specdecr<T> (amp, maxBin);
				p->matrix[Features::SPECTRAL_DECREASE].push_back(tmp);
			}
			if (p->switches[Features::SPECTRAL_SLOPE]) {
				tmp = soundmath::specslope<T> (amp, freq, maxBin);
				p->matrix[Features::SPECTRAL_SLOPE].push_back(tmp);
			}
			if (p->switches[Features::SPECTRAL_FLATNESS]) {
				tmp = soundmath::specflat<T> (amp, maxBin);
				p->matrix[Features::SPECTRAL_FLATNESS].push_back(tmp);
			}
			if (p->switches[Features::SPECTRAL_CREST]) {
				tmp = soundmath::speccrest<T> (amp, maxBin);
				p->matrix[Features::SPECTRAL_CREST].push_back(tmp);
			}
			if (p->switches[Features::HIGH_FREQUENCY_CONTENT]) {
				tmp = soundmath::hfc<T> (amp, maxBin);
				p->matrix[Features::HIGH_FREQUENCY_CONTENT].push_back(tmp);
			}
			if (p->switches[Features::MFCC1]) {
				tmp = mfcc->getCoeff (amp, 0);
				p->matrix[Features::MFCC1].push_back(tmp);
			}
			if (p->switches[Features::MFCC2]) {
				tmp = mfcc->getCoeff (amp, 1);
				p->matrix[Features::MFCC2].push_back(tmp);
			}
			if (p->switches[Features::MFCC3]) {
				tmp = mfcc->getCoeff (amp, 2);
				p->matrix[Features::MFCC3].push_back(tmp);
			}
			if (p->switches[Features::MFCC4]) {
				tmp = mfcc->getCoeff (amp, 3);
				p->matrix[Features::MFCC4].push_back(tmp);
			}
			if (p->switches[Features::MFCC5]) {
				tmp = mfcc->getCoeff (amp, 4);
				p->matrix[Features::MFCC5].push_back(tmp);
			}
			if (p->switches[Features::MFCC6]) {
				tmp = mfcc->getCoeff (amp, 5);
				p->matrix[Features::MFCC6].push_back(tmp);
			}
			if (p->switches[Features::MFCC7]) {
				tmp = mfcc->getCoeff (amp, 6);
				p->matrix[Features::MFCC7].push_back(tmp);
			}
			if (p->switches[Features::MFCC8]) {
				tmp = mfcc->getCoeff (amp, 7);
				p->matrix[Features::MFCC8].push_back(tmp);
			}
			if (p->switches[Features::MFCC9]) {
				tmp = mfcc->getCoeff (amp, 8);
				p->matrix[Features::MFCC9].push_back(tmp);
			}
			if (p->switches[Features::MFCC10]) {
				tmp = mfcc->getCoeff (amp, 9);
				p->matrix[Features::MFCC10].push_back(tmp);
			}
			if (p->switches[Features::MFCC11]) {
				tmp = mfcc->getCoeff (amp, 10);
				p->matrix[Features::MFCC11].push_back(tmp);
			}
			if (p->switches[Features::MFCC12]) {
				tmp = mfcc->getCoeff (amp, 11);
				p->matrix[Features::MFCC12].push_back(tmp);
			}
			if (p->switches[Features::MFCC13]) {
				tmp = mfcc->getCoeff (amp, 12);
				p->matrix[Features::MFCC13].push_back(tmp);
			}
			if (p->switches[Features::MFCC14]) {
				tmp = mfcc->getCoeff (amp, 13);
				p->matrix[Features::MFCC14].push_back(tmp);
			}
			if (p->switches[Features::MFCC15]) {
				tmp = mfcc->getCoeff (amp, 14);
				p->matrix[Features::MFCC15].push_back(tmp);
			}
			if (p->switches[Features::F0]) {
				// f0 = acfF0Estimate<T> (sr, buffer, result, p->N);
				f0 = soundmath::fftF0Estimate<T> (amp, freq, maxBin);
				p->matrix[Features::F0].push_back(f0);
			}
			if (p->switches[Features::INHARMONICITY]) {
				T maxAmp = 0;
				tmp = soundmath::inharmonicity<T> (amp, freq, maxBin, f0, p->sr, maxAmp);
				p->matrix[Features::INHARMONICITY].push_back(tmp);
			}
			if (p->switches[Features::TOTAL_ENERGY]) {
				tmp = soundmath::energy<T> (buffer, p->N, winEnergy);
				p->matrix[Features::TOTAL_ENERGY].push_back(tmp);
			}
			if (p->switches[Features::ZERO_CROSSINGS]) {
				tmp = soundmath::zcr<T> (buffer, p->N);
				p->matrix[Features::ZERO_CROSSINGS].push_back(tmp);
			}

			pointer += p->hop;
			++p->frames;
			if (r != 0) break;
		} while (true); //pointer <= (totSamp));
	}

	void finaliseAnalysis () {
		// scale on number of frames
		p->clusters = round (p->clusters * p->frames);

    	// backup energy
		p->env.resize (p->matrix[Features::TOTAL_ENERGY].size ());
		for (unsigned int i = 0; i < p->env.size (); ++i) p->env[i] = p->matrix[Features::TOTAL_ENERGY][i];

		// normalize features
		if (p->normalizeFeatures) {
			p->normalize (p->normalizeFeatures);
		}
	}
	void restoreClusters () {
		p->clusters /= p->frames;
	}

	void clusterAnalysis (T* data, int totSamp, int nchnls) {
					// clustering: create a matrix with only the COMPUTED features
			int n = p->frames;
			int m = p->descriptors.size ();
			soundmath::Matrix<T> matrix (n, m);
			for (unsigned int i = 0; i < p->descriptors.size (); ++i) {
				int d = p->descriptors[i];
				for (int j = 0; j < p->frames; ++j) {
					matrix[j][i] =  p->matrix[d][j]; // matrix is rows by cols
				}
			}
			p->centroids.resize (p->clusters, matrix.cols ());

			// ----------------------------------------------------------- TYPES

			getLabels (matrix);

			p->types.resize (p->clusters, p->N);
			makeDictionary (data, totSamp, nchnls);

	}

	void probabilityAnalysis () {
		int order = p->ngrams;
		int length = p->frames;
		p->hiMarkov.resize (order);

		const int* le = &p->labels[0] + length;
		std::vector<int> tmp_seq;
		for (int len = 1; len <= order; ++len) {
			typename ModelingSpace<T>::Occurences& matrix (p->hiMarkov[len - 1]);
			tmp_seq.resize (len, 0);
			for (const int* sb = &p->labels[0], *se = sb + len; se != le; ++sb, ++se) {
				std::copy (sb, se, tmp_seq.begin());
				std::pair<typename ModelingSpace<T>::Occurences::iterator, bool> ir =
					matrix.insert (std::make_pair(tmp_seq,  std::vector<T>()));
				if (ir.second) ir.first->second.resize (p->classes.size(), 0);
				++ir.first->second[*se];
			}
		}

		smoothProbabilities ();
	}
	void rematchTypes (soundmath::Matrix<T>& rtypes, soundmath::Matrix<T>& rcentroids) {
		// apply expansion
		if (p->rematchExpansionsApply) {
			for (int i = 0; i < p->centroids.rows (); ++i) {
				for (int j = 0; j < p->centroids.cols (); ++j) {
					p->centroids[i][j] *= p->rematchExpansions[j];
				}
			}
		}

		std::vector<int> rematch (p->centroids.rows ());
		//for (unsigned int i = 0; i < rematch.size (); ++i) rematch[i] = i; // identity permutation
		for (int i = 0; i < p->centroids.rows (); ++i) {
			//std::cout << "comparing " << i << " with " << std::endl;
			std::vector<T> distances;
			for (int j = 0; j < rcentroids.rows (); ++j) {

				T d = p->typesFun(&p->centroids[i][0], &rcentroids[j][0], p->centroids.cols ());
				distances.push_back (d);
				//std::cout << "\t" << j << " dist = " << d << std::endl;
			}
			int minpos = 0;
			soundmath::minimum (&distances[0], distances.size (), minpos);
			rematch[i] = minpos; // < rcentroids.rows () ? minpos : rcentroids.rows () - 1;
			//std::cout << "rematch " << i << " to " << minpos << std::endl << std::endl;
		}

		// rematch fix to avoid too many duplicates
		int oldpos = -1;
		int count = 0;
		for (unsigned int i = 0; i < rematch.size (); ++i) {
			if (rematch[i] == oldpos) {
				++count;
			}
			else {
				for (int j = 0; j < count; ++j) {
					rematch[i - count + j] = (rematch[i - count + j] + (j + 1)) % rcentroids.rows ();
				}

				count = 0;
			}
			oldpos = rematch[i];
		}

	//	for (int j = 0; j < count; ++j) {
	//		rematch[rematch.size () - 1 - count + j] += (j + 1) % rcentroids.rows ();
	//	}

		// FIXME: this is the good one!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		soundmath::BlockVocoder<T> voc (p->N);
		voc.crossMode (p->rematchMode);

		T amountIncr = (p->rematchAmountE - p->rematchAmountI) / p->types.rows();
		T amount = p->rematchAmountI;

		int tcount = 0;
		while (tcount < p->types.rows()) {
			memset (workspace, 0, sizeof (T)* 2 * p->N);
			memset (workspace2, 0, sizeof (T)* 2 * p->N);

			for (int i = 0; i < p->types.cols (); ++i) {
				workspace[2 * i] = rtypes[rematch[tcount]][i]; // rematch
				workspace2[2 * i] = p->types[tcount][i];
			}

			voc.process (workspace, workspace2, 1., 1., p->N, p->N, p->rematchC, 0, amount);

			for (int i = 0; i < p->types.cols (); ++i) {
				p->types[tcount][i] = workspace2[2 * i] / p->N;
			}

			++tcount;
			amount += amountIncr;
		}

//		int tcount = 0;
//		while (tcount < p->types.rows()) {
//
//			for (int i = 0; i < p->types.cols (); ++i) {
//				p->types[tcount][i] = (1. - amount) * p->types[tcount][i]  + amount * rtypes[rematch[tcount]][i]; // rematch
//			}
//
//			++tcount;
//			amount += amountIncr;
//		}
//
		// remove expansion
		if (p->rematchExpansionsApply) {
			for (int i = 0; i < p->centroids.rows (); ++i) {
				for (int j = 0; j < p->centroids.cols (); ++j) {
					p->centroids[i][j] /= p->rematchExpansions[j];
				}
			}
		}
	}

	void mergeProbabilities (typename ModelingSpace<T>::Occurences& rMarkov) {
		typename ModelingSpace<T>::Occurences& occ = p->hiMarkov[p->order - 1];
		typename ModelingSpace<T>::Occurences& occ1 = rMarkov;
		typename  ModelingSpace<T>::Occurences::iterator it = occ.begin ();
		typename  ModelingSpace<T>::Occurences::iterator it1 = occ1.begin ();

		while (it != occ.end()) {
			for (unsigned int i = 0; i < it->second.size(); ++i) {
				it->second[i] = (it->second[i] * (1. - p->mergeAmount)) +  (it1->second[i] * p->mergeAmount);
			}
			++it;
			++it1;
		}
	}
	void synthesize (T* synth, int totSamp, std::vector<T>& en) {
		// probabilities generation
		typename ModelingSpace<T>::Occurences::iterator probiter = p->hiMarkov[p->order - 1].begin();
		int cols = probiter->second.size ();
		int rows = p->hiMarkov[p->order - 1].size ();
		soundmath::Matrix<T> dist (rows, cols);

		if (p->decompose == 2) {
			int k = 0;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					for (dist[i][j] = k = 0; k <= j; k++) {
						T v = probiter->second[k];
						dist[i][j] += v;
					}
				}
				++probiter;
			}
		}

		soundmath::BlockVocoder<T> pv (p->N);
		int NN = p->N << 1;

		T norm = 0;
		for (int i = 0; i < p->N; ++i) {
			norm += (window[i]);
		}
		norm = 1. / ((norm * (p->N / p->hop)));

	//	std::vector<T> typePhasors ((int) p->clusters, 0);

	//	int typeFrame1 = 0;
	//	int rp1 = 0;
	//	int typeFrame2 = 0;
	//	int rp2 = 0;

		int state = 0;
		T I = p->hop / p->stretchI;
		T Ie = (T) p->hop / p->stretchE;

		T pincr = (p->pitchE - p->pitchI) / p->frames;
		T iincr = (Ie - I) / p->frames;
		T tincr = (p->thresholdE - p->thresholdI) / p->frames;
		T kincr = (p->formantsE - p->formantsI) / p->frames;

		T pp = p->pitchI;
		T ii = I;
		T tt = p->thresholdI;
		T kk = p->formantsI;

		if (p->decompose == Features::GENERATE) {
			int currentType = 0;
			std::vector<int> memory (p->order);
			typename ModelingSpace<T>::Occurences::iterator it = p->hiMarkov[p->order - 1].begin ();
			std::copy (it->first.begin (), it->first.end (), memory.begin ());

			int prevState = 0;
			int delooper = 0;

			int t = 0;
			T probWeightIncr = (p->probWeightE - p->probWeightI) / p->frames;
			T probWeight = p->probWeightI;
			int pointer = 0;
			while (t < p->frames) {
	//			int c = 0;
	//			do {
					T c0 =  soundmath::wchoice (dist[state], cols);
					T c1 = rand () % cols;
					T pdf = (T) rand () / RAND_MAX;
					int c = (int) (pdf > probWeight ? c1 : c0);
	//				std::cout << "choice = " << c << ", " << c0 << ", " << c1 << ", " << pdf << std::endl;
	//			} while (!checkCompliance (c, p, centroids));

				currentType = c;
	//			std::cout << "ctype = " << currentType << std::endl;

	//			currentType = (p->decompose == 2 ? rand () % cols : labels[t]); // TEST - random genration

				int jitter = p->jitter == 0 ? p->jitter : rand () % p->jitter;

				memset (workspace, 0, sizeof (T) * NN);
	//			int vv = (int) typePhasors[currentType];
	//
	//			typeFrame1 = classes[currentType][vv];
	//			rp1 = (int) ((T) p->hop * typeFrame1);
	//			typePhasors[currentType] += (p->phaseincr);
	//			if (typePhasors[currentType] >= classes[currentType].size ()) {
	//				typePhasors[currentType] -= classes[currentType].size ();
	//			}
	//
				for (int j = 0; j < (int) p->N; ++j) {
	//				T sample = 0;
	//				for (int k = 0; k < nchnls; ++k) {
	//					sample += (T) data[(nchnls * rp1) + (nchnls * j + k)];
	//				}
	//
	//				if (p->algorithm == Features::SEQUENCE) {
	//					workspace[2 * j] = sample;
	//				} else {
						workspace[2 * j] = p->types[currentType][j];
	//				}

					workspace[2 * j + 1] = 0;
				}

				// transformation: in -> workspace, out -> workspace2
				//pv.process (workspace, workspace2, pp, kk, p->hop, ii, p->C, tt, 0);
				//if (checkCompliance (currentType))
				pv.process (workspace, workspace2, pp, kk, p->hop, ii, p->C, tt, 0);
				//else pv.process (workspace, workspace2, 1, 1, p->hop, p->hop, 0, 0, 0);

				// overlapp-add
				for (int j = 0; j < p->N; ++j) {
					synth[j + pointer + jitter] += (workspace2[2 * j] * window[j] * norm) * p->rescaling;
				}

				// time-varying parameters
				pointer += ii;

				probWeight += probWeightIncr;
				pp += pincr;
				kk += kincr;
				ii += iincr;
				tt += tincr;

				slideIn (memory, currentType);
				int i = 0;
				for (i = 0, it = p->hiMarkov[p->order - 1].begin(); it != p->hiMarkov[p->order - 1].end (); ++it, ++i) {
					if (it->first == memory) {
						state = i;
						break;
					}
				}
				if (it == p->hiMarkov[p->order - 1].end () || state >= rows) {
					state = (int) soundmath::frand (0, rows - 1);
					it = p->hiMarkov[p->order - 1].begin ();
					std::advance (it, state);
					std::copy (it->first.begin(), it->first.end (), memory.begin ());
				}
				if (state == prevState) ++delooper;

				if (delooper == MAX_REPETITIONS) {
					state = (int) soundmath::frand (0, rows - 1);
					it = p->hiMarkov[p->order - 1].begin ();
					std::advance (it, state);
					std::copy (it->first.begin(), it->first.end (), memory.begin ());
					delooper = 0;
				}
				prevState = state;

				++t; // frame counting
			}
		} else if (p->decompose == Features::REBUILD) { // always interpolated
			std::vector<SoundPacket> packets;
			int lastType = -1;

			for (int i = 0; i < p->frames; ++i) {
				if (packets.size () == 0 || p->labels[i] != lastType) {
					state = p->labels[i];

					SoundPacket pk;
					pk.typeNumber = state;
					pk.instances = 1;
					pk.compliant = checkCompliance (state);
					packets.push_back (pk);
				} else {
					++packets[packets.size () - 1].instances;
				}
				lastType = state;
			}
			// last packet guard point for interpolation
			packets.push_back (packets[packets.size () - 1]);

			int pointer = 0;
			int cc = 0;
			for (unsigned int t = 0; t < packets.size () - 1; ++t) {
				SoundPacket pk = packets[t];
	//			typeFrame1 = classes[pk.typeNumber][(int) typePhasors[pk.typeNumber]];
	//			typePhasors[pk.typeNumber] += p->phaseincr;
	//			rp1 = (int) ((T) p->hop * typeFrame1);
	//			if (typePhasors[pk.typeNumber] >= classes[pk.typeNumber].size ()) {
	//				typePhasors[pk.typeNumber] -= classes[pk.typeNumber].size ();
	//			};

				SoundPacket pk_next = packets[t + 1];
	//			typeFrame2 = classes[pk_next.typeNumber][(int) typePhasors[pk_next.typeNumber]];
	//			typePhasors[pk_next.typeNumber] += p->phaseincr;
	//			rp1 = (int) ((T) p->hop * typeFrame1);
	//			if (typePhasors[pk_next.typeNumber] >= classes[pk_next.typeNumber].size ()) {
	//				typePhasors[pk_next.typeNumber] -= classes[pk_next.typeNumber].size ();
	//			};

				T interp = 1. / pk.instances;
				for (int i = 0; i < pk.instances; ++i)	{
					int jitter = p->jitter == 0 ? p->jitter : rand () % p->jitter;

					T alpha = (T) i * interp;
					memset (workspace, 0, sizeof (T) * NN);
					for (int j = 0; j < (int) p->N; ++j) {
	//					T sample1 = 0;
	//					T sample2 = 0;
	//					for (int k = 0; k < nchnls; ++k) {
	//						sample1 += (T) data[(nchnls * rp1) + (nchnls * j + k)];
	//						sample2 += (T) data[(nchnls * rp2) + (nchnls * j + k)];
	//					}
	//
	//					if (p->algorithm == Features::SEQUENCE) {
	//						workspace[2 * j] = (1. - alpha) * sample1 +
	//							alpha * sample2;
	//					} else {
							workspace[2 * j] = (1. - alpha) * p->types[pk.typeNumber][j] +
								alpha * p->types[pk_next.typeNumber][j];
	//					}
						workspace[2 * j + 1] = 0;
					}

					// transformation: in -> workspace, out -> workspace2
					//if (pk.compliant)
					pv.process (workspace, workspace2, pp, kk, p->hop, ii, p->C, tt, 0);
					//else pv.process (workspace, workspace2, 1, 1, p->hop, p->hop, 0, 0, 0);

					if (p->reshape) {
						T max = 0;
						for (int j = 0; j < p->N; ++j) {
							T v = fabs (workspace2[2 * j]);
							if (max <= v) max = v;
						}
						max = (max == 0 ? 1 : max);

						// overlapp-add
						for (int j = 0; j < p->N; ++j) {
							synth[j + pointer + jitter] += (workspace2[2 * j] * window[j] / max * en[cc] * p->rescaling);
						}
					} else {
						for (int j = 0; j < p->N; ++j) {
							synth[j + pointer + jitter] += (workspace2[2 * j] * window[j] * norm * p->rescaling);
						}
					}


					// time-varying parameters
					pointer += ii;

					pp += pincr;
					kk += kincr;
					ii += iincr;
					tt += tincr;

					++cc;
				}
			}

		}
	}
	void stereofy (T* in, T* out, int monoSamples, T width, int delaySamples) {
		T* delay = new T[delaySamples];
		memset (delay, 0, sizeof (T) * delaySamples);
		int ptr = 0;
		T sqrt2 = sqrt(2);

		for (int i = 0; i < monoSamples; ++i) {
			T output = delay[ptr];
			delay[ptr++] = in[i];
			ptr %= delaySamples;

			out[2 * i] = (in[i] + width * output) / sqrt2;
			out[2 * i + 1] = (in[i] - width * output) / sqrt2;
		}
	}

private:
	void createWindow (T* window) {
		if (p->wintype == "hamming") {
			soundmath::makeWindow<T> (window, p->N, .46, .54, 0);
		} else if (p->wintype == "hanning") {
			soundmath::makeWindow<T> (window, p->N, .5, .5, 0);
		} else if (p->wintype == "blackman") {
			soundmath::makeWindow<T> (window, p->N, .42, .5, .08);
		} else if (p->wintype == "bartlett") {
			T amp = 0.;
			T incr = (1. / (2 * p->N));
			for (int i = 0; i < p->N; ++i) {
				window[i] = amp;
				if (i < (p->N / 2)) {
					amp += incr;
				} else {
					amp -= incr;
				}
			}
		} else {
			throw std::runtime_error ("invalid window type specified");
		}
	}

	T clusterDispersion(soundmath::Matrix<T>& data, int m) {
		T W = 0;
		for (int i = 0; i < p->clusters; ++i) {
			T sumDist = 0;
			for (int j = 0; j < p->classes[i].size(); ++j) {
				for (int k = j + 1; k < p->classes[i].size(); ++k) {
					int p1 = p->classes[i][j];
					int p2 = p->classes[i][k];
					T d = p->distFun(data[p1], data[p2], m);

					sumDist += (d * d);
				}
			}

			T cd = (sumDist / (2. * p->classes[i].size()));
			p->dispersions.push_back(cd);
			W += cd;
		}
		// check
		if (std::isnan(W) || std::isinf(W)) return -1;
		else return W;
	}


	void makeDictionary (T* data, int totSamp, int nchnls) {

		// types building
		for (int t = 0; t < (int) p->clusters; ++t) {
			// compute distances from the center of the cluster
			std::vector <T> distances (p->classes[t].size ());
			T totDistance = 0;

			for (int j = 0 ; j < p->classes[t].size (); ++j) {
				std::vector <T> features;
				if (p->dimensions == 0) { // pca not applied
					for (unsigned int i = 0 ; i < p->descriptors.size (); ++i) {
						if (p->switches[i] == true) {
							features.push_back (p->matrix[i][p->classes[t][j]]);
						}
					}
				} else { // pca applied
					for (unsigned int i = 0; i < p->pca.size (); ++i) {
						int c = (int) p->classes[t][j];
						T f = p->pca[i][c];
						features.push_back (f);
					}
				}

				T d = p->distFun (&features[0], p->centroids[t], features.size ());
				distances[j] = d;
				totDistance += d;
			}

			if (p->algorithm == Features::RANDOM) {
				if (p->classes[t].size ()) { // skip empty centriods
					int rv = (int) rand ();
					T qr = (T) rv / RAND_MAX;
					qr *= p->classes[t].size ();
					int q = (int) qr;
					int pointer = (int) ((T) p->hop *  p->classes[t][q]);
					int r = pointer + p->N > (totSamp) ? pointer + p->N - (totSamp) : 0;
					pointer -= r;
					for (int i = 0; i < p->N; ++i) {
						T sample = 0;
						for (int j = 0; j < nchnls; ++j) {
							int pos = (nchnls * pointer) + (nchnls * i + j);
							sample += (T) data[pos];
						}
						p->types[t][i] = (sample / nchnls) * window[i];
					}
				}
			} else if (p->algorithm == Features::AVERAGE || p->algorithm == Features::WEIGHTED) {
				for (int q = 0; q < p->classes[t].size (); ++q) {
					T cnorm = 1. / p->classes[t].size ();

					T d = fabs (distances[q] / (totDistance > soundmath::EPS ? totDistance : 1));
					if (p->algorithm == Features::WEIGHTED) cnorm = 1. / (d > soundmath::EPS ? d : 1);

					int frame = p->classes[t][q];
					int pointer = (int) ((T) p->hop * frame);
					//int r = p->N > totSamp - pointer ? totSamp - pointer : p->N;
					int r = pointer + p->N > (totSamp) ? pointer + p->N - (totSamp) : 0;
					pointer -= r;
					for (int i = 0; i < p->N; ++i) {
						T sample = 0;
						for (int j = 0; j < nchnls; ++j) {
							sample += (T) data[(nchnls * pointer) + (nchnls * i + j)];
						}
						p->types[t][i] += (sample / nchnls * cnorm) * window[i];
					}
				}
			} else if (p->algorithm == Features::CLOSEST) {
				if (p->classes[t].size ()) { // skip empty centriods
					int q = 0;
					soundmath::minimum (&distances[0], distances.size (), q);

					int pointer = (int) ((T) p->hop *  p->classes[t][q]);
					//int r = p->N > totSamp - pointer ? totSamp - pointer : p->N;
					int r = pointer + p->N > (totSamp) ? pointer + p->N - (totSamp) : 0;
					pointer -= r;
					for (int i = 0; i < p->N; ++i) {
						T sample = 0;
						for (int j = 0; j < nchnls; ++j) {
							sample += (T) data[(nchnls * pointer) + (nchnls * i + j)];
						}
						p->types[t][i] = (sample / nchnls) * window[i];
					}
				}
			}
			// NB: Features::SEQUENCE is handled on generation, and is NOT save to dictionary
		}
	}
	void getLabels (soundmath::Matrix<T>& data1) {

		//FILE* stream = fopen ("work/data.txt","r");

		p->labels.resize (p->frames);

		int n = data1.rows();
		int m = data1.cols();

		// PCA reduction
		int dimsAfterReduction = p->dimensions ? p->dimensions : m;
		if (p->dimensions) {
			soundmath::Matrix<T> symmat(m, m);

			soundmath::covmat<T> (data1.data(), n, m, symmat.data());
			T* evals = new T[m];
			T* interm = new T[m];

			soundmath::tred2<T> (symmat.data(), m, evals - 1, interm - 1);
			soundmath::tqli(evals - 1, interm - 1, m, symmat.data());

			// projection to eigenvector
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					interm[j] = data1[i][j];
				}
				for (int k = 0; k < m; k++) {
					data1[i][k] = 0.;
					for (int k2 = 0; k2 < p->dimensions; k2++) {
						data1[i][k] += interm[k2] * symmat[k2][m - k]; // + 1];
					}
				}
			}
			delete [] evals;
			delete [] interm;

			data1.resize(data1.rows(), dimsAfterReduction);
			p->pca.clear();
			p->pca.resize (dimsAfterReduction);

			for (int i = 0; i < dimsAfterReduction; ++i) {
				for (int j = 0; j < n; ++j) {
					p->pca[i].push_back(data1[j][i]);
				}
			}
		}

		int h = (int) p->clusters;
		if (p->clustAlgo == Features::KMEANS) {
			soundmath::kmeans <T> (data1.data(), n, dimsAfterReduction, h,
					0.00001, &p->labels[0], p->centroids.data());
		} else if (p->clustAlgo == Features::GMM) {
			// for (int h = 1 ; h < clusters ; ++h) {
			T ratio = (T) data1.rows() / h;
			// initial guess for gmm
			soundmath::Matrix<T> means(h, dimsAfterReduction);
			for (int i = 0; i < means.rows(); ++i) {
				for (int j = 0; j < means.cols(); ++j) {
					int pos = i * ratio;
					// int pos = ((rand () % data1.rows ()) + (rand () % 10)) % data1.rows ();
					means[i][j] = data1[pos][j];
				}
			}
			soundmath::GMM<T> gmm(data1, means);
			for (int i = 0; i < (int) p->clusters; ++i) {
				for (int j = 0; j < dimsAfterReduction; ++j) {
					p->centroids[i][j] = gmm.means[i][j];
				}
			}

			for (int i = 0; i < n; ++i) {
				int max = 0;
				soundmath::maximum(gmm.resp[i], p->clusters, max);
				p->labels[i] = max;
			}
			p->loglike = gmm.loglike;
		}
		//cout << h << " " << gap << endl;
		//}
		p->classes.clear();
		p->classes.resize(h);
		for (int j = 0; j < n; ++j) {
			p->classes[p->labels[j]].push_back(j);
		}

		p->W = clusterDispersion (data1, dimsAfterReduction);
	}

	bool checkCompliance (int typeNumber) {
		bool compliant = true;
		if (p->decompositionFiltersApply) {
			for (int i = 0; i < p->centroids.cols(); ++i) { // assumes there is at least one features FIXME

				if (p->centroids[typeNumber][i] < p->decompositionFiltersMin[i] ||
					p->centroids[typeNumber][i] > p->decompositionFiltersMax[i]) {
					compliant = false;
					break;

				}
			}
		}
		return compliant;
	}

	void slideIn (std::vector<int>& v, int s) {
		for (unsigned int i = 0; i < v.size () - 1; ++i) {
			v[i] = v[i + 1];
		}
		v[v.size () - 1] = s;
	}

	void smoothProbabilities () {
		for (unsigned int i = 0; i < p->hiMarkov.size (); ++i) {
			typename ModelingSpace<T>::Occurences& level = p->hiMarkov[i];
			for (typename ModelingSpace<T>::Occurences::iterator it = level.begin (); it != level.end (); ++it) {

				T* win = new T[p->smoothWin];
				T* temp0 = new T[p->classes.size ()];
				T* temp1 = new T[p->smoothWin + p->classes.size ()];
				memset (temp1, 0, sizeof (T) * p->smoothWin + p->classes.size ());

				if (p->smoothWin == 1) {
					win[0] = 1;
				} else {
					soundmath::makeWindow<T> (win, p->smoothWin, .5, .5, 0.);
				}

				for (unsigned int j = 0; j < p->classes.size (); ++j) {
					temp0[j] = it->second[j];
				}

				soundmath::conv <T> (temp0, win, temp1, p->classes.size (), p->smoothWin, 1.);

				// normalization
				T s = soundmath::sum (temp1, p->classes.size ());

				for (unsigned int j = 0; j < p->classes.size (); ++j) {
					it->second[j] = temp1[j] / (s > DENORM ? s : 1);
				}

			}
		}
	}

	ModelingSpace<T>* p;

	soundmath::AbstractFFT<T>* ft;
	soundmath::MFCC<T>* mfcc;

	T* workspace;
	T* workspace2;
	T* freq;
	T* amp;
	T* oldAmp;
	T* phi;
	T* buffer;
	T* window;
	T* result;
};
}
#endif	// SOUNDTYPES_H

// EOF
