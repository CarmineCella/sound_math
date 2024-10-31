// sparkle.cpp
//

#include "BlockVocoder.h"
#include "WavFile.h"
#include "GetOpt.h"
#include "Tokenizer.h"

#include <iostream>
#include <stdexcept>
#include <ctime>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace soundmath;

// TODO: whitening, convolution, filtering, transients detection
//

int main (int argc, char* argv[])  {
	srand (time (NULL));

	cout << "[sparkle, ver. 1.00]" << endl << endl;
	cout << "advanced phase vocoder" << endl;
	cout << "(c) 2011-2024 by Carmine E. Cella" << endl << endl;

	// preparing buffers for allocation...
	float* input = 0;
    float* input2 = 0;
	float* output = 0;
	float* wksp = 0;
	float* wksp2 = 0;
	float* window = 0;
	float* outdata = 0;
	float* indata = 0;

	try  {
		GetOpt opt ("s:p:k:w:f:o:c:t:z:d:m:h:v");
		
		float stretch = 1;
		float P = 1;
		float K = 1;
		int Wn = 2048;
		int F = 1;
		float olap = 8;
		int C = 0;
		float threshold = 0;
		float Te = 0;
		float Se = stretch;
		float Pe = P;
		float Ke = K;
		int crossMod = 0;
        float cross = 0;
		float Ce = 0;
        string crossFile;
        int phiMod = 0;
        bool verbose = false;
        
		int c = 0;
		while ((c = opt.parse (argc, argv)) != -1) {
			switch (c) {
                case 's': {
                    string v = opt.arg ();
                    size_t sep = v.find (":");
                    if (sep == string::npos) {
                        stretch = atof (opt.arg ());
                        Se = stretch;
                    } else {
                        stretch = atof (v.substr (0, sep).c_str ());
                        Se = atof (v.substr (sep + 1, v.size () - 1).c_str ());
                    }
                    break;
                }
                case 'p': {
                    string v = opt.arg ();
                    size_t sep = v.find (":");
                    if (sep == string::npos) {
                        P = atof (opt.arg ());
                        Pe = P;
                    } else {
                        P = atof (v.substr (0, sep).c_str ());
                        Pe = atof (v.substr (sep + 1, v.size () - 1).c_str ());
                    }
                    break;
                }
                case 'k': {
                    string v = opt.arg ();
                    size_t sep = v.find (":");
                    if (sep == string::npos) {
                        K = atof (opt.arg ());
                        Ke = K;
                    } else {
                        K = atof (v.substr (0, sep).c_str ());
                        Ke = atof (v.substr (sep + 1, v.size () - 1).c_str ());
                    }
                    break;
                }
                case 'w':
                    Wn = atoi (opt.arg ());
                    break;
                case 'f':
                    F = atoi (opt.arg ());
                    break;		
                case 'o':
                    olap = atof (opt.arg ());
                    break;
                case 'c':
                    C = atoi (opt.arg ());
                    break;
                 case 't': {
                    string v = opt.arg ();
                    size_t sep = v.find (":");
                    if (sep == string::npos) {
                        threshold = atof (opt.arg ());
                        Te = threshold;
                    } else {
                        threshold = atof (v.substr (0, sep).c_str ());
                        Te = atof (v.substr (sep + 1, v.size () - 1).c_str ());
                    }
                    break;
                }
                case 'm': {
                    string v = opt.arg ();
					Tokenizer<char> tok (v.c_str (), ":");
					if (tok.size () == 3) {
						crossMod = atoi (tok[0].c_str ());
						cross = atof (tok[1].c_str ());
						Ce = cross;
						crossFile = tok[2];
					} else if (tok.size () == 4) {
						crossMod = atoi (tok[0].c_str ());
						cross = atof (tok[1].c_str ());
						Ce = atof (tok[2].c_str ());
						crossFile = tok[3];
					} else {
						throw runtime_error ("invalid syntax for cross-synthesis/morphing");
                    }

                  if (crossMod <= 0 || crossMod > 3) {
                        throw runtime_error ("invalid mode for cross-synthesis/morphing");
                    }
                    break;
                }
                case 'h':
                    phiMod = atoi (opt.arg ());
                    if (phiMod <= 0 || phiMod > 2) {
                        throw runtime_error ("invalid mode for phase modification");
                    }
                    break;
                case 'v':
                    verbose = true;
                    break;
            }
		}
		
		if (argc - opt.index () != 2) {
			throw runtime_error (
                "syntax is 'sparkle [options] in.wav out.wav'\n\n" \

                "-s<x> <x:y>               time stretch ratio\n"	\
                "-p<x> <x:y>               pitch shift ratio\n"		\
                "-k<x> <x:y>               formant move ratio\n\n"	\

                "-w<x                      window size (not needed to be power of 2)\n" \
                "-f<x>                     fft padding (shifting factor)\n" \
                "-o<x>                     synthesis overlap factor\n" \
                "-c<x>                     envelope coefficients (0 = no preservation)\n\n" \

                "-t<x> <x:y>               denoising threshold (0 for no denoising)\n" \
                "-m<t:x:file> <t:x:y:file> t = 1 cross-synthesis (multiplicative)\n" \
				"                          t = 2 cross-synthesis (spectral flattener)\n"
				"                          t = 3 morphing\n" \
                "-h<t>                     phase modification:\n" \
				"                          t = 1 robotization (requires high overlap)\n" \
				"                          t = 2 whispering (requires small windows)\n\n"					  \

                "-v                        enable verbose mode\n\n" \
    
                "in and out hop sizes are computed depending on -o and -s;\n" \
                "-s and -p arguments can be also in the form <x:y> where x and y are\n" \
                "the values of the parameter at the beginning and at the end of file;\n" \
                "-f equal to 0 sets fft size to the next power of 2 of window size,\n" \
                "-f equal to 1 sets fft size to the second power of 2 and so on;\n"\
                "default values are: w = 2048, f = 1, o = 8, c = 0, t = 0\n\n" \

                "examples:\n\n" \

                "sparkle -s2 in.wav out.wav               [stretch in.wav two times]\n" \
                "sparkle -p1.5 -c70 in.wav out.wav        [transpose in.wav preserving formants]\n" \
                "sparkle -k1:.5 -c70 in.wav out.wav       [move formants of in.wav]\n" \
                "sparkle -s1:3 in.wav out.wav             [stretch in.wav from 1 to 3 times slower]\n" \
                "sparkle -t1 in.wav out.wav               [apply a denoise on in.wav]\n"
                "sparkle -m2:.5:noise.wav in.wav out.wav  [50% cross-synth (type 2) with noise.wav]\n"
            );
		}

		WavInFile in (argv[opt.index ()]);
		WavOutFile out (argv[opt.index () + 1], in.getSampleRate (), 
            in.getNumBits (), in.getNumChannels ());
		int sr = in.getSampleRate ();
		int ch = in.getNumChannels ();
		int tot = in.getNumSamples ();

        WavInFile* in2 = 0;
        if (crossMod) {
                in2 = new WavInFile (crossFile.c_str ());
                if ((int) in2->getNumChannels  () != ch) {
                    throw runtime_error ("invalid number of channels for secondary input");
                }
        }   
        
		int WnPow2 = (int) pow (2., ceil (std::log (Wn) / std::log (2.)));
		int N  = F ? WnPow2 << F : WnPow2;
		int NN = N << 1;
		int offset = (N - Wn) / 2;

		float I = ((float) Wn / olap);
		float D = I / stretch; 
		float De = (float) I / Se;
		float Dmean = (D + De) / 2;

		int frames = (int) ceil ((float) tot / Dmean);
		int outSize = (int) (I * frames);

        if (!(Wn > 0 && N > 0 && D > 0 && I > 0 && P > 0 && Pe > 0 
			  && stretch > 0 && Se > 0 && D <= N && I <= N && C >= 0) ||
			((((~N + 1) & N) != N) || N < 4)) {
            throw runtime_error ("invalid parameters specified; check the values");
        }

        if (verbose) {
            cout << "sr            : " << sr << " Hz" << endl;
            cout << "channels      : " << ch << endl;
            cout << "input samples : " << ch * tot << endl;
            cout << "output samples: " << ch * outSize << endl;
            cout << "total frames  : " << frames << endl;
            cout << "stretch       : " << stretch << "x";
            if (Se != stretch) cout << " -> " << Se << "x" << endl;
            else cout << endl;
            cout << "pitch         : " << P << "x";
            if (Pe != P) cout << " -> " << Pe << "x" << endl;
            else cout << endl;
            cout << "win size      : " << Wn << " samples" << endl;
            cout << "fft size      : " << N << " samples" << endl;
            cout << "input hop     : " << D;
            if (Se != stretch) cout << " -> " << De << " samples" << endl;
            else cout << " samples" << endl; 
            cout << "output hop    : " << I << " samples" << endl;
            cout << "envelope      : " << C;
            if (!C) cout << " (no preservation)" << endl;
            else cout << endl;
            cout << "threshold     : " << threshold;
            if (threshold != Te) cout << " -> " << Te << "x" << endl;
            else cout << endl;
			stringstream tmp;
			tmp << crossFile.c_str () << " (type " << crossMod << ")";
            cout << "cross/morph   : " << (crossMod ?  tmp.str (): "none") << endl;
            cout << "phase modifier: ";
            if (phiMod == 1) cout << "robot";
            else if (phiMod == 2) cout << "whisper";
            else cout << "none";
            cout << endl << endl;
        }

        // warnings
        if (C && (P == 1 && Pe == 1 && K == 1 && Ke == 1 && !crossMod)) cout << "[warning: envelope detection without transposition, formant move or x-synth]" << endl << endl;
		if (crossMod == 2 && !C) cout << "[warning: x-synth (type 2) without envelope detection]" << endl << endl;
		if ((P != 1 || Pe != 1) && !C) cout << "[warning: pitch-shift without envelope detection]" << endl << endl;    
		if ((K != 1 || Ke != 1) && !C) cout << "[warning: formant move without envelope detection]" << endl << endl;    
		if (phiMod == 1 && olap < 16) cout << "[warning: robotization works better with high overlap]" << endl  << endl;
        if (phiMod == 2 && N > 256) cout << "[warning: whisperization works better with small windows]" << endl << endl;
        
        // channels
		cout << "allocating........."; cout.flush ();

       	input = new float[ch * tot + N];
		memset (input, 0, sizeof (float) * (ch * tot + N));
        
        if (crossMod) {
            input2 = new float[ch * tot + N]; // same size of input 
            memset (input2, 0, sizeof (float) * (ch * tot + N));
        }
        
		// output buffers
		output = new float[ch * outSize + NN]; // is + NN correct??
		memset (output, 0, sizeof (float) * (ch * outSize + NN));

		// working buffers
		wksp = new float[NN];
		memset (wksp, 0, sizeof (float) * NN);
		wksp2 = new float[NN];
		memset (wksp2, 0, sizeof (float) * NN);

		outdata = new float[Wn];
		memset (outdata, 0, sizeof (float) * Wn);
		indata = new float[Wn];	
		memset (indata, 0, sizeof (float) * Wn);
		
		// window
		window = new float[Wn];
		hanningz (window, Wn);

		// normalization
		float norm = 0;
		for (int i = 0; i < Wn; ++i) {
			norm += (window[i]);
		}
		//float pamp = P < Pe ? P : Pe;
	    norm = 1. / ((norm * (N / I))); // * (pamp < 1 ? pamp : 1));
		
		float pincr = (Pe - P) / frames;
		float kincr = (Ke - K) / frames;
		float dincr = (De - D) / frames;
		float cincr = (Ce - cross) / frames;
		float tincr = (Te - threshold) / frames;

		// the phase vocoder object
		BlockVocoder<float> pv (N);

		cout << "done\nreading............"; cout.flush ();
	    in.read (input, ch * tot);
        if (crossMod) {
			int r = in2->read (input2, ch * tot);
			for (int i = r; i < ch * tot; ++i) {
				input2[i] = input2[i - r]; // loop cross-data if shorter than carrier
			}
		}

		pv.phaseMode (phiMod);
		pv.crossMode (crossMod);

		cout << "done\nprocessing........."; cout.flush ();
		clock_t tic = clock ();
		for (int j = 0; j < ch; ++j) {
			float read = 0;
			float write = 0;
			int f = 0;
			float pp = P;
			float dd = D;
			float cc = cross;
			float tt = threshold;
			float kk = K;

			pv.reset ();
			while (read < tot) {		
				// slide in
				 memset (wksp, 0, sizeof (float) * NN);
                 if (crossMod) memset (wksp2, 0, sizeof (float) * NN);
				 for (int i = (int) read ; i < (int) (read + Wn); ++i) {
				 	wksp[2 * (i + offset - (int) read)] = input[ch * i + j] * window[i - (int) read];
                     if (crossMod) wksp2[2 * (i + offset - (int) read)] = input2[ch * i + j] * window[i - (int) read];
				 }

				 // transformation: in -> wksp, out -> wksp2
				 pv.process (wksp, wksp2, pp, kk, dd, I, C, tt, cc);

				 // overlapp-add
				 for (int i = (int) write; i < (int) (write + Wn); ++i) {
				 	output[ch * i + j] += (wksp2[2 * (i + offset - (int) write)]
				 						   * window[i - (int) write] * norm);
				 }
					
				read += dd;
				write += I;

				pp += pincr;
				kk += kincr;
				dd += dincr;
				cc += cincr;
				tt += tincr;

				++f;
			}
		}
		clock_t toc = clock ();

		cout << "done\nwriting............"; cout.flush ();
		out.write (output, ch * outSize);
		cout << "done" << endl << endl;

		cout << "processing time: " << (float) (toc - tic) / CLOCKS_PER_SEC << " sec." << endl;
	} catch (bad_alloc& memmoryAllocationException) {
		cout << "memory allocation exception: "
			 << memmoryAllocationException.what ()
			 << endl;
	} catch (exception& e) {
		cout << "error: " << e.what () << endl;
	} catch (...)  {
		cout << "fatal error; bailing out.." << endl;
	}
	
	// clean up
	delete [] input;
    delete [] input2;
	delete [] output;

	delete [] wksp;
	delete [] wksp2;
	delete [] outdata;
	delete [] indata;
	
	delete [] window;

	return 0;
}
	
// eof
