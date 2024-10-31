// combsyn.cpp
// 

#include "OscillatorF.h"
#include "Comb.h"
#include "Delay.h"
#include "WavFile.h"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <cstring>
#include <ctime>
#include <memory>
#include <cstdlib>

using namespace std;
using namespace soundmath;

void scale (float* inout, float factor, int size) {
	for (int i = 0; i < size; ++i)  {
		inout[i] *= factor;
	}
}

float bartlett (float in, int pos, int samples, int slope) {
	if (pos <= slope) {
		return ((float) pos / slope) * in;
	} else if (pos >= (samples - slope)) {
		return (float) (samples - pos) / slope * in;
	} else return in;
}

float frand (float min, float max) {
	float f = 0;
	short r = (short) (rand ());
	f = fabs ((float) r / 32768);
	f *= (max - min);
	f += min;
	return f;
}
		
// -------------------------------------------------------------------------- //

const int BSIZE = 1024;
const int ARGMIN = 15;

int main (int argc, char* argv[]) {
	try {
		srand (time (NULL));
		
		cout << "[combsyn, ver. 1.00]" << endl << endl;
		cout << "minimalistic comb-based synthesizer" << endl;
		cout << "(c) 2011-2024 by Carmine E. Cella" << endl << endl;

		if (argc < ARGMIN) {
			throw runtime_error ("combSyn in.wav out.wav dur tmin tmax tmine tmaxe smin smax dmin dmax norm sl f1...fn\n"\
				"where:\n"\
				"\tdur            is total time in seconds\n" \
				"\ttmin - tmax    are min and max for time reduction coefficients at the beginning\n" \
				"\ttmine - tmaxe  are min and max for time reduction coefficients at the end\n" \
				"\tsmin - smax    are min and max for 'soundness' coefficients\n" \
				"\tdmin - dmax    are min and max for damping\n" \
				"\tnorm           is the normalization coefficient\n" \
				"\tsl             is the slope for the envelope shape\n" \
				"\tfreq1....freqn are the frequencies to be synthesized\n"
			);
		}
		
		WavInFile input (argv[1]);
		//if (input.getNumChannels () != 1) throw runtime_error ("invalid format: only mono 16 bit files allowed");
		int ch = input.getNumChannels ();
		float chnorm = 1. / ch;
		float sr = input.getSampleRate ();		
		int samples = input.getNumSamples ();
		float nF = sr / (float) samples; // natural frequency
		
		WavOutFile output (argv[2], sr, 16, 2);
		
		float dur = atof (argv[3]);
		if (dur <= 0) throw runtime_error ("invalid duration");
		float tmin = atof (argv[4]);
		float tmax = atof (argv[5]);
		if (tmax < tmin || tmin <= 0 || tmax <= 0) throw runtime_error ("invalid tmin/tmax parameters [0.eps : 1]");
		float tmine = atof (argv[6]);
		float tmaxe = atof (argv[7]);
		if (tmaxe < tmine || tmine <= 0 || tmaxe <= 0) throw runtime_error ("invalid tmine/tmaxe parameters [0.eps : 1]");		
		float smin = atof (argv[8]);
		float smax = atof (argv[9]);
		if (smax < smin || smin < 0 || smax < 0 || smin >= 1 || smax >= 1) throw runtime_error ("invalid smin/smax parameters [0 : .99999]");
		float dmin = atof (argv[10]);
		float dmax = atof (argv[11]);
		if (dmax < dmin || dmin < 0 || dmax < 0 || dmin >= 1 || dmax >= 1) throw runtime_error ("invalid dmin/dmax parameters [0 : .99999]");
		
		float norm = atof (argv[12]);
		if (norm < 0) throw runtime_error ("invalid normalization [>0]");
		
		float sl = atof (argv[13]);
		if (sl < 0 || sl > .5) throw runtime_error ("invalid slope [0 : .5]");

		int outsamp = (int) (dur * sr);
		int buffers = (int) (((float) outsamp / (float) BSIZE));
			
		vector <float> frequencies;
		vector <float> pannings;
		vector <float> pitches;
		vector <float> pincrs;
		
		for (int i = ARGMIN - 1; i < argc; ++i) {
			float f = atof (argv[i]);
			frequencies.push_back (f);
			pannings.push_back (frand (0, 1));
			float pin = nF * frand (tmin, tmax);
			float pend = nF * frand (tmine, tmaxe);
			pincrs.push_back ((pend - pin) / buffers);			
			pitches.push_back (pin);
		}
		
		float* buff = new float[BSIZE];
		float* obuffS = new float[BSIZE];
		float* obuff = new float[BSIZE];

		float* leftright = new float[2 * BSIZE];

		vector <float> sound ((int) (2. * outsamp), 0);
		
		cout << "reading data...";
		cout.flush ();
		float* data2 = new float[ch * samples];
		input.read (data2, ch * samples);
		float* data = new float[samples + 1]; // guardpoint

		for (int i = 0; i < samples; ++i) {
			data[i] = 0;
			for (int j = 0; j < ch; ++j) {
				data[i] += data2[i * ch + j];
			}
			data[i] *= chnorm;
		}
		data[samples] = data[samples - 1];
		
		cout << "done" << endl << endl;
		cout.flush ();
		
		OscillatorF<float>** osc = new OscillatorF<float>*[frequencies.size ()];
		Comb<float>** cb = new Comb<float>*[frequencies.size ()];
		Delay<float>** dl = new Delay<float>*[frequencies.size ()];
		
		float stot = 0.;
		for (unsigned int i = 0; i < frequencies.size (); ++i) {
			osc[i] = new OscillatorF<float> (sr, data, samples);
			osc[i]->frequency (nF * pitches[i]);
			
			float s = frand (smin, smax);
			stot += s;
			cb[i] = new Comb<float> (sr, (float) (1. / frequencies[i]), s);
			cb[i]->damp (dmin, dmax);
			int d = (int) frand (11, 73);
			dl[i] = new Delay<float> (sr, d, 0); // stereo decorrelation			
		}
		float cnorm = 1. - (stot / frequencies.size ());
		if (cnorm < .1) cnorm = .1;
		
		float anorm = 1. / frequencies.size ();
		
		cout << "total streams    : " << frequencies.size () << endl;
		cout << "out buffers      : " << buffers << endl;	
		cout << "out samples      : " << samples << endl;	
		cout << "natural frequency: " << nF << endl;		
		cout << "normalizations   : " << cnorm << ", " << anorm << endl << endl;
	
		cout << "computing... ";
		cout.flush ();
		int pos = 0;
		
		for (int i = 0; i < buffers; ++i) {
			memset (leftright, 0, sizeof (float) * 2 * BSIZE);
						
			for (unsigned int j = 0; j < frequencies.size (); ++j) {
				osc[j]->process (buff, BSIZE);
				osc[j]->frequency (pitches[j]);
				pitches[j] += pincrs[j];
				scale (buff, cnorm, BSIZE); // comb normalization
				cb[j]->process (buff, obuff, BSIZE);
				scale (obuff, anorm, BSIZE); // additive normalization
				
				dl[j]->process (obuff, obuffS, BSIZE);
				
				for (int k = 0; k < BSIZE; ++k) {
					leftright[k * 2] += obuff[k] * pannings[j];
					leftright[k * 2 + 1] += obuffS[k] * (1. - pannings[j]);
				}
			}
				
			std::copy (leftright, leftright + (2 * BSIZE), &sound[pos]);
			pos += 2 * BSIZE;
		}
		
		pos = 0;
		float slope = (int) (sl * outsamp);
		for (int i = 0; i < outsamp; ++i) {
			sound[i * 2] = bartlett (sound[i * 2], pos, outsamp, slope) * norm;
			sound[i * 2 + 1] = bartlett (sound[i * 2 + 1], pos, outsamp, slope) * norm;
			++pos;
		}
		
		cout << "done\nwriting...   ";
		cout.flush ();
		output.write ((float*) &sound[0], sound.size ());
		cout << "done\n" << endl;
	    
	    delete [] data;
		delete [] data2;
	    delete [] buff;
	    delete [] obuff;
	    delete [] leftright;
	    
		for (unsigned int i = 0; i < frequencies.size (); ++i) {
			delete osc[i];
			delete cb[i];
			delete dl[i];
		}	    
		
		delete [] osc;
		delete [] cb;
		delete [] dl;
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

