// ModelingSpace.h
// 

#ifndef MODELINGSPACE_H
#define MODELINGSPACE_H 

#include "Matrix.h"
#include "DynamicMatrix.h"
#include "algorithms.h"
#include "utilities.h"

#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>

struct Features {
	enum Algorithm {
		RANDOM = 0,
		AVERAGE = 1,
		WEIGHTED = 2,
		CLOSEST = 3,
		SEQUENCE = 4
	};

	enum Performance {
		NONE = 0,
		REBUILD = 1,
		GENERATE = 2
	};
	
	enum XMode {
		MULTIPLICATIVE = 1,
		FLATTENING = 2,
		MORPHING = 3
	};
	enum Clustering {
		KMEANS = 0, 
		GMM = 1
	};

	enum AmpliutdeScale {
		LINEAR,
		LOG,
		POWER
	};

	enum DescriptorID {
		SPECTRAL_CENTROID = 0,
		SPECTRAL_SPREAD = 1,
		SPECTRAL_SKEWNESS = 2,
		SPECTRAL_KURTOSIS = 3,
		SPECTRAL_FLUX = 4,
		SPECTRAL_IRREGULARITY = 5,	
		SPECTRAL_DECREASE = 6,
		SPECTRAL_SLOPE = 7,	
		SPECTRAL_FLATNESS = 8, 
		SPECTRAL_CREST = 9, 
		HIGH_FREQUENCY_CONTENT = 10,
		MFCC1 = 11,
		MFCC2 = 12,
		MFCC3 = 13,
		MFCC4 = 14,
		MFCC5 = 15,
		MFCC6 = 16,
		MFCC7 = 17,
		MFCC8 = 18,
		MFCC9 = 19,
		MFCC10 = 20,
		MFCC11 = 21,
		MFCC12 = 22,
		MFCC13 = 23,
		MFCC14 = 24,		
		MFCC15 = 25,
		F0 = 26,
		INHARMONICITY = 27,
		ZERO_CROSSINGS = 28,		
		TOTAL_ENERGY = 29,
		DESCRIPTORS_NUMBER = TOTAL_ENERGY + 1
	};

	const std::string& id2str (DescriptorID n) {
		// NB: must be in the same order of the ENUM!
		static const std::string NAMES[2 * DESCRIPTORS_NUMBER] = {
			"centroid",		 "spectral",
			"spread",		 "spectral",
			"skewness",		 "spectral",
			"kurtosis",		 "spectral",
			"flux",			 "spectral",
			"irregularity",	 "spectral",
			"decrease",		 "spectral",
			"slope",		 "spectral",
			"flatness",	     "spectral",
			"crest",         "spectral",
			"hfc",			 "spectral",
			"mfcc1",		 "perceptual",
			"mfcc2",		 "perceptual",
			"mfcc3",		 "perceptual",
			"mfcc4",		 "perceptual",
			"mfcc5",		 "perceptual",
			"mfcc6",		 "perceptual",
			"mfcc7",		 "perceptual",
			"mfcc8",		 "perceptual",
			"mfcc9",		 "perceptual",
			"mfcc10",		 "perceptual",
			"mfcc11",		 "perceptual",
			"mfcc12",		 "perceptual",
			"mfcc13",		 "perceptual",
			"mfcc14",		 "perceptual",				
			"mfcc15",		 "perceptual",				
			"f0",		     "harmonic",
			"inharmonicity", "harmonic",
			"zcr",			 "temporal",			
			"energy",		 "temporal"
		};	
		
#ifdef PLOT_FEATURES
			std::ofstream dot ("features_graph.txt");
			dot << "digraph G {" << std::endl;
			dot << "\trankdir = LR" << std::endl;
			dot << "\tnode [shape = box]" << std::endl << std::endl;
			for (int i = 0; i < DESCRIPTORS_NUMBER; ++i) {
				dot << "\t\"" << NAMES[2 * i + 1] << "\" -> \"" << NAMES[2 * i] << "\"" << std::endl;
			}
			dot << "}" << std::endl;
			dot.close ();
		}
#endif
		return NAMES[(int) 2 * n];
	}

	DescriptorID str2id (const std::string& name) {
		DescriptorID pos = (DescriptorID) -1;
		for (int i = 0; i < DESCRIPTORS_NUMBER; ++i) {
			if (name == id2str ((DescriptorID) (i))) {
				pos = (DescriptorID) i;
				break;
			}
		}
		
		if (pos == -1) throw std::runtime_error ("invalid descriptor specified");
		return pos;
	}	
};

template <typename T>
class ModelingSpace : public Features {
public:
	typedef std::map <std::vector<int>, std::vector<T> > Occurences;	
	typedef typename Occurences::iterator OccurencesIt;	
	typedef typename std::vector<Occurences>::iterator VOccurencesIt;
	
	ModelingSpace () {
      	matrix.resize (DESCRIPTORS_NUMBER);
		switches.resize (DESCRIPTORS_NUMBER);
		for (int i = 0; i < DESCRIPTORS_NUMBER; ++i) switches[i] = false;
		
		analysisPerformance = "compute";
		sr = 44100.;
		wintype = "blackman";
		N = 4096;
		hop = 512;
		limit = 5512.5;
		scale = POWER;
		modeling = false;
		normalizeFeatures = false;
		
		typesPerformance = "compute";
		showLabels = false;
		dimensions = 0;
		clustAlgo = GMM;
		clusters = .1;
		rematchTypes = "none";
		rematchAmountI = 0;	
		rematchAmountE = 0;	
		rematchC = 0;
		rematchExpansionsApply = false;
		rematchExpansions.clear ();
		rematchMode = MORPHING;
		typesFun = &soundmath::taxicab;
		
		probabilitiesPerformance = "compute";
		ngrams = 1;
		probWeightI = 1;
		probWeightE = 1;
		smoothWin = 1;
		mergeProb = "none";
		mergeAmount = 0;
		
		decompose = 1;
		reshape = true;
		order = 1;
		decompositionFiltersApply = false;
		//phaseincr = 1;
		distFun = &soundmath::edistance;
		stretchI = 1.;
		stretchE = 1.;
		pitchI = 1.;
		pitchE = 1.;
		formantsI = 1.;
		formantsE = 1.;
		C = 0;
		thresholdI = 0;
		thresholdE = 0;
		rescaling = .25;
		jitter = 0;
		algorithm = RANDOM;		
		stereoWidth = 0;
		stereoDelay = 23;

		verbose = false;
		inputFile = "UNDEFINED";
		outputPath = ".";
        
   		clear ();
	}
	virtual ~ModelingSpace () {
		// delete m_matrix;
	}
	void normalize (int method) {
		for (int i = 0; i < matrix.rows(); ++i) {
			if (matrix[i].size ()) {
				switch (method) {
					default:
					case 1:
						soundmath::normalize(&matrix[i][0], &matrix[i][0],  matrix[i].size ());
					break;
						
					case 2:
						soundmath::normalize2(&matrix[i][0], &matrix[i][0],  matrix[i].size ());
					break;
				}
			}
		}
	}
	
	void computeModelings (int frames) {		
		for (unsigned int i = 0; i < descriptors.size (); ++i) {
			float m = soundmath::centroid (&matrix[Features::TOTAL_ENERGY][0], 
				&matrix[descriptors[i]][0], frames);
			means.push_back (m);
			stddevs.push_back (soundmath::stddev (&matrix[Features::TOTAL_ENERGY][0],
				&matrix[descriptors[i]][0], m, frames));
		}			
				
	}
	
	int dimensionsAfterReduction () const {
		return dimensions ? dimensions : centroids.cols ();
	}

	void configure (std::ifstream& config) {
		int line = 0;
		while (!config.eof ()) {
			std::string inp;
			std::string opcode;

			++line;
			std::getline (config, inp, '\n');

			if (inp.size () == 0) continue;
			inp = soundmath::trim (inp);
			std::istringstream istr (inp, std::ios_base::out);

			std::vector <std::string> tokens;
			while (!istr.eof ()) {
				istr >> opcode;		
				if (opcode.size ())	tokens.push_back (opcode);	
			}   									

			if (tokens[0][0] == ';') continue;
			if (tokens.size () < 2) {
				std::stringstream err;
				err << "invalid syntax at line " << line;
				throw std::runtime_error (err.str ());
			}


			// ANALYSIS ------------------------------------------------------------
			
			if (tokens[0] == "analysis.performance") {
				if (tokens[1] == "compute" || tokens[1] == "load") analysisPerformance = tokens[1];
				else throw std::runtime_error ("invalid performance specified for analysis");
			} 			else if (tokens[0] == "analysis.sr") {
				 sr = atof (tokens[1].c_str ());
			} else if (tokens[0] == "analysis.win.type") {
				 wintype = tokens[1];
			} else if (tokens[0] == "analysis.win.size") {
				 N = atoi (tokens[1].c_str ());
			} else if (tokens[0] == "analysis.win.hop") {
				 hop = atoi (tokens[1].c_str ());
			} else if (tokens[0] == "analysis.maxfreq") {
				 limit = atof (tokens[1].c_str ());
			} else if (tokens[0] == "analysis.ampscale") {
				if (tokens[1] == "lin") {		
					 scale = Features::LINEAR;
				} else if (tokens[1] == "log") {		
					 scale = Features::LOG;
				} else if (tokens[1] == "pow") {
					 scale = Features::POWER;
				} else {
					throw std::runtime_error ("invalid amplitude scale specified");
				}
			} else if (tokens[0] == "analysis.normalize") {
				 normalizeFeatures = (int) atoi (tokens[1].c_str ());
			} else if (tokens[0] == "analysis.modeling") {
				if (tokens[1] == "yes") modeling = true;
				else modeling = false;
			} else if (tokens[0] == "analysis.descriptors") {
				for (unsigned int i = 1; i < tokens.size (); ++i) {
					Features::DescriptorID d =  str2id (tokens[i]);
					if (! find (d)) {
						 descriptors.push_back (d);
						 switches[(int) d] = true;
					}
				}
			}

	//		else if (tokens[0] == "analysis.mfccs") {
	//			 mfccNum = (bool) atoi (tokens[1].c_str ());
	//		}

			// TYPES ---------------------------------------------------------------
			else if (tokens[0] == "types.performance") {
				if (tokens[1] == "compute" || tokens[1] == "load") typesPerformance = tokens[1];
				else throw std::runtime_error ("invalid performance specified for types");
			} 
			else if (tokens[0] == "types.showlabels") {
				if (tokens[1] == "yes") showLabels = true;
				else showLabels = false;
			} else if (tokens[0] == "types.algorithm") {
				if (tokens[1] == "gmm") {
					 clustAlgo = Features::GMM;
				} else if (tokens[1] == "kmeans") {
					 clustAlgo = Features::KMEANS;
				} else {
					throw std::runtime_error ("invalid clustering algorithm pecified");
				}
			} else if (tokens[0] == "types.clusters") {
				 clusters = atof (tokens[1].c_str ());
			} else if (tokens[0] == "types.dimensions") {
				 dimensions = atoi (tokens[1].c_str ());
			} else if (tokens[0] == "types.dictionary") {
				if (tokens[1] == "random") {
					 algorithm = Features::RANDOM;
				} else if (tokens[1] == "average") {
					 algorithm = Features::AVERAGE;
				} else if (tokens[1] == "weighted") {
					 algorithm = Features::WEIGHTED;
				}  else if (tokens[1] == "closest") {
					 algorithm = Features::CLOSEST;
				}  
	//			else if (tokens[1] == "sequence") {
	//				 algorithm = Features::SEQUENCE;
	//			}			
				else {
					throw std::runtime_error ("invalid dictionary specified");
				}
			} 
			
			else if (tokens[0] == "types.rematch.file") {
				 rematchTypes = tokens[1];
			} else if (tokens[0] == "types.rematch.amount") {
				size_t sep = tokens[1].find (":");
				if (sep == std::string::npos) {
					 rematchAmountI = atof (tokens[1].c_str ());
					 rematchAmountE =  rematchAmountI;
				} else {
					 rematchAmountI = atof (tokens[1].substr (0, sep).c_str ());
					 rematchAmountE  = atof (tokens[1].substr (sep + 1, tokens[1].size () - 1).c_str ());
				}
			}  else if (tokens[0] == "types.rematch.envelope") {
				 rematchC = atof (tokens[1].c_str ());
			} else if (tokens[0] == "types.rematch.weights.apply") {
				if (tokens[1] == "yes") rematchExpansionsApply = true;
				else rematchExpansionsApply = false;
			}
			else if (tokens[0] == "types.rematch.weights") {
				for (unsigned int i = 1; i < tokens.size (); ++i) {
					 rematchExpansions.push_back (atof (tokens[i].c_str ()));
				}	
			} else if (tokens[0] == "types.rematch.xmode") {
				if (tokens[1] == "mult") {
					 rematchMode = Features::MULTIPLICATIVE;
				} else if (tokens[1] == "flatten") {
					 rematchMode = Features::FLATTENING;
				} else if (tokens[1] == "morph") {
					 rematchMode = Features::MORPHING;
				} else {
					throw std::runtime_error ("invalid performance specified");
				}
			}  else if (tokens[0] == "types.rematch.distance") {
				if (tokens[1] == "euclid") {
					 typesFun = &soundmath::edistance;
				} else if (tokens[1] == "mahalanobis") {
					 typesFun = &soundmath::mahalanobis;
				} else if (tokens[1] == "taxicab") {
					 typesFun = &soundmath::taxicab;
				}  else if (tokens[1] == "cosim") {
					 typesFun = &soundmath::cosineSimilarity;
				} else if (tokens[1] == "kullback") {
					 typesFun = &soundmath::kullbackLeibler;
				} else {
					throw std::runtime_error ("invalid distance specified");
				}
			}

			// PROBABILITIES -------------------------------------------------------
			else if (tokens[0] == "probabilities.performance") {
				if (tokens[1] == "compute" || tokens[1] == "load") probabilitiesPerformance = tokens[1];
				else throw std::runtime_error ("invalid performance specified for probabilities");
			} 
			else if (tokens[0] == "probabilities.ngrams") {
				 ngrams = atof (tokens[1].c_str ());
			} else if (tokens[0] == "probabilities.merge.file") {
				 mergeProb = tokens[1];
			} else if (tokens[0] == "probabilities.smoothwin") {
				 smoothWin = atoi (tokens[1].c_str ());
			}  
			else if (tokens[0] == "probabilities.merge.amount") {
				 mergeAmount = atof (tokens[1].c_str ());
			} 

			// DECOMPOSITION -------------------------------------------------------
			else if (tokens[0] == "synthesis.performance") {
				if (tokens[1] == "none") {
					 decompose = Features::NONE;
				} else if (tokens[1] == "rebuild") {
					 decompose = Features::REBUILD;
				} else if (tokens[1] == "generate") {
					 decompose = Features::GENERATE;
				} else {
					throw std::runtime_error ("invalid performance specified for synthesis");
				}
			} else if (tokens[0] == "synthesis.reshape") {
				if (tokens[1] == "yes") reshape = true;
				else reshape = false;
			}
	//		else if (tokens[0] == "synthesis.phasor") {
	//			 phaseincr = atof (tokens[1].c_str ());
	//		} 
			else if (tokens[0] == "synthesis.gen.order") {
				 order = atoi (tokens[1].c_str ());
			}  else if (tokens[0] == "synthesis.gen.weight") {
				size_t sep = tokens[1].find (":");
				if (sep == std::string::npos) {
					 probWeightI = atof (tokens[1].c_str ());
					 probWeightE =  probWeightI;
				} else {
					 probWeightI = atof (tokens[1].substr (0, sep).c_str ());
					 probWeightE  = atof (tokens[1].substr (sep + 1, tokens[1].size () - 1).c_str ());
				}
			} else if (tokens[0] == "synthesis.distance") {
				if (tokens[1] == "euclid") {
					 distFun = &soundmath::edistance;
				} else if (tokens[1] == "mahalanobis") {
					 distFun = &soundmath::mahalanobis;
				} else if (tokens[1] == "taxicab") {
					 distFun = &soundmath::taxicab;
				}  else if (tokens[1] == "cosim") {
					 distFun = &soundmath::cosineSimilarity;
				} else if (tokens[1] == "kullback") {
					 distFun = &soundmath::kullbackLeibler;
				} else {
					throw std::runtime_error ("invalid distance specified");
				}
			} else if (tokens[0] == "synthesis.stretch") {
				size_t sep = tokens[1].find (":");
				if (sep == std::string::npos) {
					 stretchI = atof (tokens[1].c_str ());
					 stretchE =  stretchI;
				} else {
					 stretchI = atof (tokens[1].substr (0, sep).c_str ());
					 stretchE  = atof (tokens[1].substr (sep + 1, tokens[1].size () - 1).c_str ());
				}
			} else if (tokens[0] == "synthesis.pitch") {
				size_t sep = tokens[1].find (":");
				if (sep == std::string::npos) {
					 pitchI = atof (tokens[1].c_str ());
					 pitchE =  pitchI;
				} else {
					 pitchI = atof (tokens[1].substr (0, sep).c_str ());
					 pitchE  = atof (tokens[1].substr (sep + 1, tokens[1].size () - 1).c_str ());
				}
			}  else if (tokens[0] == "synthesis.formants") {
				size_t sep = tokens[1].find (":");
				if (sep == std::string::npos) {
					 formantsI = atof (tokens[1].c_str ());
					 formantsE =  formantsI;
				} else {
					 formantsI = atof (tokens[1].substr (0, sep).c_str ());
					 formantsE  = atof (tokens[1].substr (sep + 1, tokens[1].size () - 1).c_str ());
				}
			} else if (tokens[0] == "synthesis.envelope") {
				 C = atoi (tokens[1].c_str ());
			} else if (tokens[0] == "synthesis.threshold") {
				size_t sep = tokens[1].find (":");
				if (sep == std::string::npos) {
					 thresholdI = atof (tokens[1].c_str ());
					 thresholdE =  thresholdI;
				} else {
					 thresholdI = atof (tokens[1].substr (0, sep).c_str ());
					 thresholdE  = atof (tokens[1].substr (sep + 1, tokens[1].size () - 1).c_str ());
				}
			} else if (tokens[0] == "synthesis.rescaling") {
				 rescaling = atof (tokens[1].c_str ());
			} else if (tokens[0] == "synthesis.jitter") {
				 jitter = atoi (tokens[1].c_str ());
			}else if (tokens[0] == "synthesis.stereo.width") {
				 stereoWidth = atof (tokens[1].c_str ());
			}else if (tokens[0] == "synthesis.stereo.delay") {
				 stereoDelay = atoi (tokens[1].c_str ());
			} 
			else if (tokens[0] == "synthesis.filters.apply") {
				if (tokens[1] == "yes") decompositionFiltersApply = true;
				else decompositionFiltersApply = false;
			} 
			else if (tokens[0] == "synthesis.filters.min") {
				for (unsigned int i = 1; i < tokens.size (); ++i) {
					 decompositionFiltersMin.push_back (atof (tokens[i].c_str ()));
				}	
			} else if (tokens[0] == "synthesis.filters.max") {
				for (unsigned int i = 1; i < tokens.size (); ++i) {
					 decompositionFiltersMax.push_back (atof (tokens[i].c_str ()));
				}	
			}

			// GLOBAL --------------------------------------------------------------
			else if (tokens[0] == "global.verbose") {
				 verbose = (bool) atoi (tokens[1].c_str ());
			} else if (tokens[0] == "global.input") {
				 inputFile = tokens[1];
			} else if (tokens[0] == "global.outputpath") {
				 outputPath = tokens[1];
			} else {
				std::stringstream err;
				err << "invalid token in configuration file at line " << line;
				throw std::runtime_error (err.str ());
			} 
		}
	}

	int checkParameters (int totSamp) {
		int maxBin = 0;
		// ERRORS
		if (sr < 8000) {
			throw std::runtime_error ("unsopported sampling rate");
		}
		if ( descriptors.size () < 1) {
			throw std::runtime_error ("invalid number of descriptors");
		}
		if ((((~ N +  N) &  N) !=  N) ||  N < 2) {
			throw std::runtime_error ("window length must be a power of two");
		}

		if ( hop < 1 ||  hop >  N) {
			throw std::runtime_error ("invalid hop size");
		}		

		float step = sr /  N;
		if ( limit < step ||  limit > (sr * .5)) {
			throw std::runtime_error ("invalid limit speficied");
		}
		maxBin = (int) ( limit / step);

		if ( clusters < 0 ||  clusters > 1) {
			throw std::runtime_error ("invalid number of centroids");
		}
		if ( dimensions < 0 ||  dimensions > (int)  descriptors.size ()) {
			throw std::runtime_error ("invalid dimensions for pca reduction");
		}
		if ( ngrams < 0 || ngrams > totSamp - 1) {
			std::cout << ngrams << " " << std::cout << totSamp << std::endl;
			throw std::runtime_error ("invalid ngrams level");
		}
		if ( order >  ngrams ||  order < 0) {
			throw std::runtime_error ("invalid order for decomposition");
		}

	//	if ( phaseincr <= 0) {	
	//		throw std::runtime_error ("invalid value for phasor");
	//	}

		if ( probWeightI < 0 ||  probWeightI > 1) {
			throw std::runtime_error ("invalid amout for probabilities initial weight");
		}	

		if ( probWeightE < 0 ||  probWeightE > 1) {
			throw std::runtime_error ("invalid amout for probabilities final weight");
		}	

		if ( smoothWin <= 0) {
			throw std::runtime_error ("invalid size for probabilities smoothing window");
		}	

		if ( mergeAmount < 0 ||  mergeAmount > 1) {
			throw std::runtime_error ("invalid amout for merge");
		}	

		if ( rematchAmountI < 0 ||  rematchAmountI > 1) {
			throw std::runtime_error ("invalid amount for initial types rematch");
		}	

		if ( rematchAmountE < 0 ||  rematchAmountE > 1) {
			throw std::runtime_error ("invalid amount for final types rematch");
		}		

		if ((rematchExpansionsApply == true && rematchExpansions.size () != descriptors.size ()) 
			&&  rematchTypes != "none") {
			throw std::runtime_error ("inconsistent number of weights for types rematch");
		}

		if ( modeling == true &&  switches[Features::TOTAL_ENERGY] != true) {
			throw std::runtime_error ("energy needed to compute temporal modelings");
		}
		if ( switches[Features::SPECTRAL_SPREAD] == true && 
			 switches[Features::SPECTRAL_CENTROID] != true) {
			throw std::runtime_error ("centroid needed to compute spread");

		}
		if (reshape == true && switches[Features::TOTAL_ENERGY] != true) {
			throw std::runtime_error ("energy needed to perform rebuild");
		}

		if ((switches[Features::SPECTRAL_SKEWNESS] == true || switches[Features::SPECTRAL_KURTOSIS] == true) &&
			 switches[Features::SPECTRAL_SPREAD] != true) {
			throw std::runtime_error ("centroid and spread needed to compute skewness or kurtosis");
		}
		
		if ( switches[Features::INHARMONICITY] == true &&
			! switches[Features::F0] == true) {
			throw std::runtime_error ("f0 needed to compute inharmonicity");
		}

//		if ( (decompositionFiltersApply == true && decompose != Features::NONE) && 
//				( decompositionFiltersMin.size () !=  decompositionFiltersMax.size () ||
//				  decompositionFiltersMin.size () !=  descriptors.size () ||
//				  decompositionFiltersMax.size () !=  descriptors.size ())) {
//				throw std::runtime_error ("inconsistent number of filters for decomposition");
//		}

		// WARNINGS
//		if ( rematchTypes != "none" &&  algorithm == Features::SEQUENCE) {
//			std::cout << "warning: rematch incompatible with selected decomposition" << std::endl;
//		}
		if (( rematchMode == Features::MULTIPLICATIVE ||  rematchMode == Features::MORPHING) 
				&&  rematchC != 0) {
			std::cout << "warning: incorrect number of coefficients for selected xmode" << std::endl;
		}

		if (( pitchI != 1 ||  pitchE != 1 ||  formantsI != 1 ||  formantsE != 1) &&  C == 0) {
			std::cout << "warning: pitch/formants move without envelope preservation" << std::endl;
		}
		if ( decompose != Features::NONE &&  clusters == 0) {
			std::cout << "warning: decomposition will not be performed; not enough clusters" << std::endl;
		}
		
		return maxBin;
	}

	void printLog (int nchnls, std::ostream& out) {
		out << "input file    = [" <<  inputFile << "]" << std::endl;
		out << "output path   = [" <<  outputPath << "]" << std::endl << std::endl;

		out << "sr            = " << sr << " Hz" << std::endl;
		out << "channels      = " << nchnls << std::endl;
		out << "window type   = " <<  wintype << std::endl;		
		out << "window size   = " <<  N << " samples" << std::endl;
		out << "hop size      = " <<  hop << " samples" << std::endl;
		out << "freq. limit   = " <<  limit << std::endl;
		out << "amp  scale    = ";

		switch ( scale) {
			case Features::LINEAR: out << "linear" << std::endl; break;
			case Features::LOG: out << "log" << std::endl; break;
			case Features::POWER: out << "power" << std::endl; break;
		}

		out << "pitch shift   = " <<  pitchI << " -> " <<  pitchE << std::endl;
		out << "formants move = " <<  formantsI << " -> " <<  formantsE << std::endl;
		out << "time stretch  = " <<  stretchI << " -> " <<  stretchE << std::endl;
		out << "denoising     = " <<  thresholdI << " -> " <<  thresholdE << std::endl;

		if ( clusters) {
			out << "clustering    = " << 
				( clustAlgo == Features::KMEANS ? "kmeans" : "gmm") << std::endl;

			if ( dimensions) {
				out << "pca dimens.   = " <<  dimensions << std::endl;
			}
			out << "ngrams        = " <<  ngrams << std::endl;

			if ( rematchTypes != "none") {
				out << "types rematch = " <<  rematchAmountI * 100  << "% -> " <<  rematchAmountE * 100  
					<< "%, mode " <<  rematchMode << " (weights ";
					int loop = 0;
					for (unsigned int i = 0; i <  descriptors.size (); ++i) {
						out <<  rematchExpansions[i] << " ";
						if ((++loop % 4) == 0) out << std::endl << "              ";
					}
					out << ")" << std::endl;
			}		
			if ( mergeProb != "none") {
				out << "prob. merging = " <<  mergeAmount * 100 << "%" << std::endl;
			}
			if ( decompose) {
				out << "decomposing   = ";
				switch ( algorithm) {
					case Features::RANDOM: out << "random, "; break;
					case Features::AVERAGE: out << "average, "; break;
					case Features::WEIGHTED: out << "weighted, "; break;
					case Features::CLOSEST: out << "closest, "; break;								 
					case Features::SEQUENCE: /*out << "sequence (phasor = " <<  phaseincr << "), "; */ break;								 
				}

				switch ( decompose) {
					case Features::REBUILD: out << "rebuild"; 
						if (reshape) out << " (reshaped)";
					break;
					case Features::GENERATE: out << "generate (order = " <<  order <<
							", weight = " <<  probWeightI << " -> " <<  probWeightE << ")"; break;
				}

				out << std::endl;
			}

		}

		out << "feature(s)    = ";

		int loop = 0;
		for (unsigned int i = 0; i <  descriptors.size (); ++i) {
			out <<  id2str ( descriptors[i]) << " ";
			if ((++loop % 4) == 0) out << std::endl << "              ";
		}
		if ( normalizeFeatures) out << "(normalized - " <<  normalizeFeatures << ")";
		out << std::endl;
	}
	void clear () {
		//for (int i = 0; matrix.rows (); ++i) {
        //    matrix[i].clear ();
        //}
	
        matrix.clear ();
		pca.clear ();
		means.clear ();
		stddevs.clear ();
		labels.clear ();
		centroids.clear ();
		classes.clear ();
		dispersions.clear ();
		types.clear ();
		hiMarkov.clear ();
		env.clear ();
	
		W = 0;
		loglike = 0;
        
        frames = 0;
	}
	
	void saveTypes (const std::string& outputPosition) {	
		std::stringstream tmp;
		tmp.str ("");
		tmp << outputPosition << ".centroids.dictionary.wav";

		soundmath::WavOutFile audioDictionary (tmp.str ().c_str (), sr, 16, 1);
		for (int i = 0; i < (int) clusters; ++i) {
			audioDictionary.write (types[i], N);
		}	

		tmp.str ("");
		tmp << outputPosition << ".labels.raw";			
		soundmath::serialize <int> (tmp.str(), &labels[0], labels.size());

		tmp.str ("");
		tmp << outputPosition << ".dispersions.raw";
		soundmath::serialize <T> (tmp.str(), &dispersions[0], dispersions.size ());


		tmp.str ("");
		tmp << outputPosition << ".centroids.txt";
		std::ofstream ctrds (tmp.str().c_str());


		ctrds << (int) clusters << " " << dimensionsAfterReduction () << std::endl;
		for (int i = 0; i < (int) clusters; ++i) {
			for (int j = 0; j < dimensionsAfterReduction (); ++j) {
				ctrds << std::setw (15) << centroids[i][j];
			}
			ctrds << std::endl;
		}
		
		ctrds.close ();
	}	

	void saveProbabilities (const std::string& outputPosition) {
		std::stringstream tmp;

		std::ofstream plot;

		tmp.str ("");
		tmp << outputPosition << ".probabilityGraph.txt";


		plot.open (tmp.str ().c_str());
		plot << "digraph G {" << std::endl;
		plot << "\trankdir = LR" << std::endl;
		plot << "node [shape = circle, color=blue]" << std::endl << std::endl;


		for (unsigned int i = 0; i < hiMarkov.size (); ++i) {
			typename ModelingSpace<T>::Occurences& level = hiMarkov[i];

			std::ofstream probs;
		
			std::stringstream fname;
			fname << outputPosition << ".probabilities" << i + 1 << ".txt";

			probs.open (fname.str().c_str());
			//probs << level.size() << " " << level.begin ()->second.size () << std::endl;
			
			for (typename ModelingSpace<T>::Occurences::iterator it = level.begin (); it != level.end (); ++it) {
				int cols = it->second.size ();
				for (unsigned int j = 0; j < it->first.size(); ++j) {
					probs << it->first[j];
					if (j != it->first.size () - 1) probs << " ";
				}					
				probs << "\t\t";
				for (int j = 0; j < cols; ++j) {
					probs << it->second[j] << "\t";
				}
				probs << std::endl;
				
				std::stringstream tmp;
				tmp << "\"";
				for (unsigned int j = 0; j < it->first.size(); ++j) {
					tmp << it->first[j];
					if (j != it->first.size () - 1) tmp << " ";
				}
				tmp << "\"";

				for (unsigned int j = 0; j < it->second.size(); ++j) {
					if (it->second[j] != 0) {
						plot << "\t" << tmp.str () << " -> " << j << " [label = " 
							 << it->second[j] << "]" << std::endl;

						if (it->first.size () == 1)  {
							plot << tmp.str () << "[style=filled, fillcolor=red]" << std::endl;
						}
					}					
				}
			}				

			probs.close ();

		}

		
		plot << "}" << std::endl;
		plot.close ();
	}

	void saveFeatures (const std::string& outputPosition) {
		for (unsigned int i = 0; i < descriptors.size (); ++i) {
			std::string n = (std::string) id2str ((Features::DescriptorID) 
										  descriptors[i]) + (std::string) ".raw";

			std::stringstream tmp;
			tmp.str ("");
			tmp << outputPosition << "." << n;
			soundmath::serialize <T> (tmp.str(), &matrix[descriptors[i]][0], matrix[descriptors[i]].size());
		}

		if (clusters) {
			for (int i = 0; i < dimensions; ++i) {
				std::stringstream n0;
				n0 << outputPosition << ".pca_dim" << i + 1 << ".raw";
				soundmath::serialize <T> (n0.str (), &pca[i][0], pca[i].size());
			}
		}
	}

	void loadTypes (const std::string& inputPosition) {
		std::stringstream centrs;
		centrs << inputPosition << ".centroids.txt";

		//std::cout << centrs.str () << std::endl;
		std::ifstream ctrin (centrs.str ().c_str ());
		if (!ctrin.good ()) {
			throw std::runtime_error ("cannot open file for centroids");
		}	

		int rows = 0, cols = 0;
		ctrin >> rows;
		ctrin >> cols;
		centroids.resize (rows, cols);
		
		int rcount = 0;
		while (!ctrin.eof () && rcount < rows) {
			int ccount = 0;
			while (!ctrin.eof () && ccount < cols) {
				float v = 0;
				ctrin >> v; 
				centroids[rcount][ccount]  = v ; // * (max - min); // rescale target
				//std::cout << v <<  "\t";
				++ccount;
			}
			++rcount;
			//std::cout << std::endl;
		}
		//getchar ();

		std::stringstream dict; 
		dict << soundmath::removeExtension (centrs.str ()) << ".dictionary.wav";

		soundmath::WavInFile in (dict.str ().c_str ());

		if (in.getSampleRate () != sr || in.getNumChannels () != 1) {
			throw std::runtime_error ("invalid type for dictionary");
		}

		int tcount = 0;
		types.resize(rows, N);
		while (!in.eof () && tcount < rows) { // should be the same
			in.read (types[tcount], types.cols ());
			++tcount;
		}

		clusters = rows;
	}

	void loadProbabilities (const std::string& inputPosition) {
		hiMarkov.resize (order);
		for (int i = 0; i < order; ++i) {
			typename  ModelingSpace<T>::Occurences& it = hiMarkov[i];
			
			std::stringstream prbsfile;
			prbsfile << inputPosition << ".probabilities" << i  + 1<< ".txt"; // 1-based files
			std::cout << prbsfile.str () << std::endl;
			
			std::ifstream in (prbsfile.str ().c_str ());
			if (in.good ()) {

				while (!in.eof ()) {
					std::string inp;
					std::string opcode;

					std::getline (in, inp, '\n');

					inp = soundmath::trim (inp);
					if (inp.size () == 0)  {
						continue;
					}

					std::stringstream linp;
					linp << inp;
					
					std::vector <int> indexes;
					std::vector <T> tokens;

					int symNum = 0;
					while (!linp.eof ()) {
						linp >> opcode;			
						if (symNum <= i) {
							indexes.push_back (atoi (opcode.c_str ())); // symbol
							//std::cout << "symbol: " << atoi (opcode.c_str ()) << std::endl;
						}
						else {
							tokens.push_back (atof (opcode.c_str ()));	
							//std::cout << "prob: " << atof (opcode.c_str ()) << std::endl;
						}
						++symNum;
					}   									
					//std::cout << std::endl;
					//it[indexes] = tokens;
					it.insert(std::make_pair<std::vector<int>, std::vector<T> > (indexes, tokens));
				}
				//getchar ();
			}	
			else throw std::runtime_error ("cannot open file for probabilities");
		}
	}
	
	void loadFeatures (const std::string& inputPosition) {
		for (unsigned int i = 0; i < descriptors.size (); ++i) {
			std::string n = (std::string) id2str ((Features::DescriptorID) 
										  descriptors[i]) + (std::string) ".raw";
 
			std::stringstream tmp;
			tmp.str ("");
			tmp << inputPosition << "." << n;
			soundmath::deserialize (tmp.str(), &matrix[descriptors[i]][0], matrix[descriptors[i]].size());
		}

		if (clusters) {
			for (int i = 0; i < dimensions; ++i) {
				std::stringstream n0;
				n0 << inputPosition << ".pca_dim" << i + 1 << ".raw";
				soundmath::deserialize <T> (n0.str (), &pca[i][0], pca[i].size());
			}
		}
		
		
		// scale on number of st.nframes ()
		clusters = round (clusters * matrix[descriptors[0]].size ());		
		
    	// backup energy
		if (normalizeFeatures == 1)	soundmath::unnormalize(&matrix[Features::TOTAL_ENERGY][0], &env[0], matrix[Features::TOTAL_ENERGY].size());
		else if (normalizeFeatures == 2) soundmath::unnormalize2(&matrix[Features::TOTAL_ENERGY][0], &env[0], matrix[Features::TOTAL_ENERGY].size());		
	}
	
	// analysis-dependant data
	soundmath::DynamicMatrix<T> matrix;
	soundmath::DynamicMatrix<T> pca;
	std::vector <T> means;
	std::vector <T> stddevs;
	std::vector<int> labels;
	soundmath::Matrix<float> centroids;
	soundmath::DynamicMatrix<float> classes;
	std::vector <T> dispersions;
	soundmath::Matrix<T> types;	
	std::vector<Occurences> hiMarkov;
	std::vector<T> env;
	
	T W;
	T loglike;
    int frames;
			
	// configuration data
	std::string analysisPerformance;
	T sr;
	std::string wintype;
	int N;
	int hop;
	T limit;
	bool showLabels;
	AmpliutdeScale scale;
	bool modeling;
	std::vector<Features::DescriptorID> descriptors;
	std::vector <bool> switches;
	
	std::string typesPerformance;
	int dimensions;
	Clustering clustAlgo;
	T clusters;
	std::string rematchTypes;
	T rematchAmountI;				
	T rematchAmountE;
	T rematchC;
	std::vector<T> rematchExpansions;
	bool rematchExpansionsApply;
	XMode rematchMode;
	T (*typesFun)(const T*, const T*, int);
	
	std::string probabilitiesPerformance;
	int ngrams;
	T probWeightI;
	T probWeightE;
	int smoothWin;
	std::string mergeProb;
	T mergeAmount;
	
	int decompose;
	bool reshape;
	//T phaseincr;
	int order;
	T (*distFun)(const T*, const T*, int);
	T stretchI;
	T stretchE;
	T pitchI;
	T pitchE;
	T formantsI;
	T formantsE;
	T C;
	T thresholdI;
	T thresholdE;
	T rescaling;
	int jitter;
	Algorithm algorithm;
	int normalizeFeatures;
	T stereoWidth;
	int stereoDelay;
	bool decompositionFiltersApply;
	std::vector<T> decompositionFiltersMin;
	std::vector<T> decompositionFiltersMax;
	
	bool verbose;
	std::string inputFile;
	std::string outputPath;
	
	bool find (Features::DescriptorID d) {
		bool found = false;
		for (unsigned int i = 0; i < descriptors.size (); ++i) {
			if (descriptors[i] == d) {
				found = true;
				break;
			}
			
		}
		return found;
	}
};

#endif	// MODELINGSPACE_H 

// EOF
