// sounder.cpp
//


#include "Tokenizer.h"
#include "SpectralChord.h"
#include "GeometricDesigner.h"
#include "HarmonicDesigner.h"
#include "RandomDesigner.h"
#include "ModelBasedDesigner.h"

#include "WavFile.h"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>

using namespace std;
using namespace soundmath;

const int BSIZE = 1024;
const int TABLEN = 4096;
const float SR = 44100;

template <typename T>
void mul (T* buff, T v, int N) {
	for (int i = 0; i < N; ++i) {
		buff[i] *= v;
	}
}

template <typename T>
void checkParams (Parameters<T>* params, AbstractDesigner<T>* designer, int line) {
	if (params->start < 0 || params->dur <= 0) {
		stringstream tmp;
		tmp << "invalid start/duration at line " << line;
		throw runtime_error (tmp.str ());
	}
	if (params->fund <= 0) {
		stringstream tmp;
		tmp << "invalid fundamental at line " << line;
		throw runtime_error (tmp.str ());
	}
	if (params->partials < 2 || params->thickness < 1) {
		stringstream tmp;
		tmp << "invalid partials/thickness at line " << line;
		throw runtime_error (tmp.str ());
	}
	
	if (params->spread < 0) {
		stringstream tmp;
		tmp << "invalid spread at line " << line;
		throw runtime_error (tmp.str ());
	}

	string t = designer->type ();
	if (t == "GeometricDesigner" && params->coeff <= 1) {
		std::cout << params->coeff << std::endl;
		stringstream tmp;
		tmp << "invalid coefficients at line " << line;
		throw runtime_error (tmp.str ());
	}
	if (t == "HarmonicDesigner" && params->coeff <= 0) {
		std::cout << params->coeff << std::endl;
		stringstream tmp;
		tmp << "invalid coefficients at line " << line;
		throw runtime_error (tmp.str ());
	}
	if (t == "RandomDesigner" && (params->coeff <= 0 || params->coeff > 1)) {
		std::cout << params->coeff << std::endl;
		stringstream tmp;
		tmp << "invalid coefficient at line " << line;
		throw runtime_error (tmp.str ());
	}

	if (params->slope < 0 || params->slope > 1) {
		stringstream tmp;
		tmp << "invalid slope at line " << line;
		throw runtime_error (tmp.str ());
	}	
}

template <typename T>
int configure (T sr, std::ifstream& config, std::vector<SpectralChord<T>*>& events,
			   T* tab, int tablen) {
	int line = 0;
	int tsamps = 0;
	while (!config.eof ()) {
		std::string inp;
		std::string opcode;

		++line;
		std::getline (config, inp, '\n');
		
		if (inp.size () == 0) continue;
		
		std::istringstream istr (inp, std::ios_base::out);

		std::vector <std::string> tokens;
		while (!istr.eof ()) {
			istr >> opcode;			
			tokens.push_back (opcode);	
		}   									
		
		if (tokens[0][0] == ';') continue;
		if (tokens.size () != 9) {
			std::stringstream err;
			err << "invalid syntax at line " << line;
			throw std::runtime_error (err.str ());
		}
		
		SpectralChord<T>* gs = new SpectralChord<T> (sr, tab, tablen);
		
		AbstractDesigner<T>* designer = 0;
		if (tokens[0] == "geom") {
			designer = new GeometricDesigner<T> ();
		} else if (tokens[0] == "harm") {
			designer = new HarmonicDesigner<T> ();
		} else if (tokens[0] == "rand") {
			designer = new RandomDesigner<T> ();
		} else if (tokens[0].find ("model") != string::npos) {
            Tokenizer<char> tok (tokens[0], ":");
            if (tok.size () == 3) {
                string name = tok[1].c_str ();
                T time = atof (tok[2].c_str ());
				
				WavInFile in (name.c_str ());
				if (in.getNumChannels () != 1) {
					throw runtime_error ("invalid number of channels in the model");
				}
				ModelBasedDesigner<T>* mb = new ModelBasedDesigner<T> ();
				designer = mb;
				T* data = new T[in.getNumSamples ()];
				in.read (data, in.getNumSamples ());
				int msecPos = time * 1000;
				int N = 4096;
				
				mb->analyse (in.getSampleRate (), data, in.getNumSamples (), msecPos, N);
				delete [] data;
			} else {
                throw runtime_error ("invalid syntax for model based designer");
            }

		} else {
			std::stringstream err;
			err << "invalid design at line " << line;
			throw std::runtime_error (err.str ());
		}
		
		Parameters<T> params;
		params.start = atof (tokens[1].c_str ());
		params.dur = atof (tokens[2].c_str ());
		params.fund = atof (tokens[3].c_str ());
		params.partials = atoi (tokens[4].c_str ());
		params.thickness = atoi (tokens[5].c_str ());
		params.spread = atof (tokens[6].c_str ());
		params.coeff = atof (tokens[7].c_str ());
		params.slope = atof (tokens[8].c_str ());
				
		checkParams (&params, designer, line);

		cout << "[new event found]" << std::endl;	
		cout << "\ttype       : " << designer->type () << std::endl;
		cout << "\tcoefficient: " << params.coeff << std::endl;
		cout << "\tstart      : " << params.start << " sec." << std::endl;
		cout << "\tduration   : " << params.dur << " sec." << std::endl;
		cout << "\tfundamental: " << params.fund << " Hz" << std::endl;
		cout << "\tpartials   : " << params.partials << std::endl;
		cout << "\tthickness  : " << params.thickness << std::endl;
		cout << "\tspread     : " << params.spread << std::endl << std::endl;
		
		gs->design (&params, designer);
		
		delete designer;
		
		events.push_back (gs);
		
		int s = (int) ((params.start + params.dur) * sr);
		if (s > tsamps) tsamps = s;
 	}
	
	return tsamps;
}


int main (int argc, char* argv[]) {
	cout << "[sounder, ver. 1.00]" << endl << endl;
	cout << "spectral sound synthesis" << endl;
	cout << "(c) 2011-2024 by Carmine E. Cella" << endl << endl;
	
	srand (time (NULL));
	
	try {
		
		if (argc < 3 || argc > 4) {
			throw runtime_error ("syntax is 'sounder source.txt out.wav [scale]'");
		}
		
		ifstream cfg (argv[1]);
		if (!cfg.good ()) {
			throw runtime_error  ("cannot open source file");
		}
		
		cout << "parsing source..." << endl << endl; cout.flush ();

 		float* buff = new float[2 * BSIZE];
		float* tab = new float[TABLEN + 1];
		memset (tab, 0, sizeof (float) * (TABLEN + 1));
		for (int i = 0; i < TABLEN; ++i) {
			tab[i] = sin (2. * M_PI * i / TABLEN);
		}
		tab[TABLEN] = tab[0];
		
		vector<SpectralChord<float>*> events;
		int tsamps = configure<float> (SR, cfg, events, tab, TABLEN);

		int buffers = tsamps / BSIZE;
		
		WavOutFile out (argv[2], SR, 16, 2);
		
		float scale = 1;
		if (argc == 4) {
			scale = atof (argv[3]);
		}
		
		cout << "computing....."; cout.flush ();
		for (int j = 0; j < buffers; ++j) {
			memset (buff, 0, sizeof (float) * 2 * BSIZE);
			for (unsigned int i = 0; i < events.size (); ++i) {
				SpectralChord<float>* c = events[i];
				int sb = c->start () / BSIZE;
				int eb = (c->length () / BSIZE) + sb;
				if (sb <= j && j <= eb) c->process (buff, BSIZE);
			}
			mul (buff, scale, 2 * BSIZE);
			out.write (buff, 2 * BSIZE);
		}
		cout << "done" << endl;
		
		delete [] buff;
		delete [] tab;
	}
	catch (exception& e) {
		cout << "Error: " << e.what () << endl;
	}
	catch (...) {
		cout << "Fatal error: unknown exception" << endl;
	}
	return 0;
}

// EOF

