// tools.h
//

#ifndef TOOLS_H
#define TOOLS_H

#include "SoundTypes.h"
#include "ModelingSpace.h"
#include "WavFile.h"
#include "utilities.h"

#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>

extern const int SPACING = 14;

template <typename T>
void saveAndPrintLog (const std::string& outputPosition, ModelingSpace<T>& p,
	int nchnls) {
	std::stringstream tmp;
	p.printLog (nchnls, std::cout);
	tmp.str ("");
	tmp << outputPosition << ".log.txt";
	std::ofstream logfile (tmp.str ().c_str ());
	p.printLog (nchnls, logfile);
	logfile.close();
}

template <typename T>
void saveAndPrintModelings (const std::string& outputPosition, ModelingSpace<T>& p) {
	std::stringstream tmp;
	tmp.str ("");
	tmp << outputPosition << ".modelings.txt";

	std::ofstream modelings (tmp.str ().c_str ());

	std::cout << std::endl << "[modelings]" << std::endl;
	std::cout << "means     : ";


	int loop = 0;
	for (unsigned int i = 0; i < p.means.size (); ++i) {
		std::cout << std::setw (SPACING) << p.means[i];
		modelings <<  std::setw (SPACING) << p.means[i];
		if ((++loop % 4) == 0) std::cout << std::endl
												   << "            ";
	}
	modelings << std::endl;
	std::cout << std::endl << "stddevs   : ";
	loop = 0;
	for (unsigned int i = 0; i < p.stddevs.size (); ++i) {
		std::cout << std::setw (SPACING) << p.stddevs[i];
		modelings << std::setw (SPACING) << p.stddevs[i];
		if ((++loop % 4) == 0) std::cout << std::endl
												   << "            ";
	}
	modelings << std::endl;
	std::cout << std::endl << std::endl;
}

template <typename T>
void loadTypes (soundmath::Matrix<T>& rtypes, soundmath::Matrix<T>& rcentroids,
		ModelingSpace<T>& p) {
	std::stringstream centrs;
	centrs << p.rematchTypes;

	std::ifstream ctrin (centrs.str ().c_str ());
	if (!ctrin.good ()) {
		throw std::runtime_error ("cannot open file for centroids");
	}
	int tcount = 0;

	while (!ctrin.eof () && tcount < rcentroids.rows ()) {

		for (int i = 0; i < rcentroids.cols (); ++i) {
			float v = 0;
			ctrin >> v;
			rcentroids[tcount][i]  = v ; // * (max - min); // rescale target
		}
		++tcount;
	}

	std::stringstream dict;
	dict << soundmath::removeExtension (p.rematchTypes) << ".dictionary.wav";

	soundmath::WavInFile in (dict.str ().c_str ());

	if (in.getSampleRate () != p.sr || in.getNumChannels () != 1) {
		throw std::runtime_error ("invalid type for dictionary");
	}

	tcount = 0;
	while (!in.eof () && tcount < rtypes.rows()) { // should be the same
		in.read (rtypes[tcount], rtypes.cols ());
		++tcount;
	}
}

template <typename T>
void loadProbabilities (typename ModelingSpace<T>::Occurences& rMarkov,
	ModelingSpace<T>& p) {
	typename  ModelingSpace<T>::Occurences::iterator it = rMarkov.begin ();

	std::ifstream in (p.mergeProb.c_str());
	if (in.good ()) {
		while (!in.eof () && it != rMarkov.end()) {
			std::string inp;
			std::string opcode;

			std::getline (in, inp, '\n');

			if (inp.size () == 0)  {
				continue;
				++it;
			}

			std::istringstream istr (inp, std::ios_base::out);

			std::vector <T> tokens (it->second.size ());
			unsigned int row = 0;
			while (!istr.eof () && row <= tokens.size ()) {
				istr >> opcode;
				tokens[row] = (atof (opcode.c_str ()));

				++row;
			}
			std::cout << std::endl;
			// nb: this will create normalized rows always (inserting zeros if needed)
			for (unsigned int i = 0; i < it->second.size(); ++i) {
				it->second[i] = tokens[i];
			}
			++it;
		}

		// nb: missing rows will be fixed during merge
	}
	else throw std::runtime_error ("cannot open file for probabilities merging");
}

template <typename T>
void synthesizeAndSave (const std::string& outputPosition,
	soundmath::SoundTypes<T>& st, ModelingSpace<T>& p, int totSamp) {
	// allocate space for synthesized samples
	T mstretch = 1. / ((p.stretchI + p.stretchE) / 2);
	int synsam = (int) ((T) (totSamp + p.N + p.jitter) * mstretch) + (p.N);
	T* synth = new T[synsam];
	memset (synth, 0, sizeof (T) * synsam);

	T* stereoOut = new T[synsam * 2];
	memset (stereoOut, 0, sizeof (T) * 2 * synsam);

	st.synthesize (synth, totSamp, p.env);

	std::stringstream tmp;
	tmp.str ("");
	tmp << outputPosition << ".output.wav";

	if (p.stereoWidth) {
		st.stereofy (synth, stereoOut, synsam, p.stereoWidth, p.stereoDelay);
		soundmath::WavOutFile synthesized (tmp.str ().c_str (), p.sr, 16, 2);
		synthesized.write (stereoOut, synsam * 2);
	}
	else {
		soundmath::WavOutFile synthesized (tmp.str ().c_str (), p.sr, 16, 1);
		synthesized.write (synth, synsam);
	}

	delete [] synth;
	delete [] stereoOut;
}

template <typename T>
void printAnalysis (ModelingSpace<T>& p, soundmath::SoundTypes<T>& st) {

	int* counters = new int[(int) p.clusters];
	memset (counters, 0, sizeof (int) * (int) p.clusters);

	for (int i = 0; i < p.frames; ++i) {
		++counters[p.labels[i]];
	}

	T avg = (T) p.frames / (int) p.clusters;

	for (int i = 0; i < (int) p.clusters; ++i) {
		std::cout << "[cluster " << i + 1 << "]" << std::endl;
		std::cout << "element(s): " << std::setw (SPACING) << counters[i];
		if (counters[i] < (avg / 4)) std::cout << "*" << std::endl;
		else std::cout << std::endl;

		int loop = 0;
		std::cout << "centroid  : ";
		for (int j = 0; j < p.dimensionsAfterReduction (); ++j) {
			std::cout << std::setw (SPACING) << p.centroids[i][j];
			if ((++loop % 4) == 0) std::cout << std::endl
													   << "            ";
		}
		std::cout << std::endl << std::endl;
	}

	if (p.showLabels) {
		std::cout << std::endl << "[labels]" << std::endl;
		int loop = 0;
		for (int j = 0; j < p.frames; ++j) {
			std::cout << std::setw (SPACING) << j << "::" << p.labels[j];
			if ((++loop % 4) == 0) std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	delete [] counters;
}

template <typename T>
void printSummary (ModelingSpace<T>& p, soundmath::SoundTypes<T>& st) {
	std::cout << "within-cluster dispersion: ";
	if (p.W < 0) std::cout << "undefined" << std::endl;
	else std::cout << p.W << std::endl;
	if (p.clustAlgo == Features::GMM) {
		std::cout << "bayesian criterion (BIC) : " << -2 * p.loglike
			+ p.dimensionsAfterReduction () * log (p.frames) << std::endl;
	}
	std::cout << std::endl << p.frames << " element(s), " << (int) p.clusters
		 << " cluster(s) ["
		 << (float) p.clusters / p.frames * 100 << "%]" << std::endl;

}

#endif	// TOOLS_H

// EOF
