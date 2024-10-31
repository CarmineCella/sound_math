// pconv.cpp
//

#include "WavFile.h"
#include "BlockConv.h"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace soundmath;

const int N = 4096; //feed buffer size: FFT size will be 2 * N

// -------------------------------------------------------------------------- //

int main (int argc, char* argv[]) {
    cout << "[pconv, ver. 1.00]" << endl << endl;
    
    cout << "fast-convolution engine" << endl << endl;
    cout << "written by carmine e. cella" << endl;
    cout << "(c) 2011 www.sonic-pad.com" << endl << endl;
    
	try {
		if (argc < 5 || argc > 6) {
			throw runtime_error ("syntax is 'pconv in.wav ir.wav out.wav scale [mix]'");
		}

		WavInFile input (argv[1]);
		if (input.getNumBits () != 16) 
            throw runtime_error ("invalid format for input file (only 16 bit allowed)");
		unsigned int sr = input.getSampleRate ();
        int inch = input.getNumChannels ();
        
		WavInFile irinput (argv[2]);
		if (irinput.getNumBits () != 16 || irinput.getSampleRate () != sr) 
            throw runtime_error ("invalid format for impulse response (only 16 bit allowed)");
		int irsamples = irinput.getNumSamples ();
        int irch = irinput.getNumChannels ();
        
        if (inch > 2 || irch > 2) throw runtime_error ("unsupported combination of channels");
        
		float *imp0 = new float[irsamples];
        float *imp1 = new float[irsamples];
		float *impL = new float[irsamples / 2];
		float *impR = new float[irsamples / 2];
                
		float *out0 	= new float[N];
		float *out1 	= new float[N];
        float *out2 	= new float[N / 2];
		float *out3 	= new float[N / 2];
		float *in		= new float[N];
        float *inL		= new float[N / 2];
        float *inR		= new float[N / 2];        

        float mix = 0;
        if (argc == 6) mix = atof (argv[5]);
    
        int pattern = inch * 10 + irch;
        cout << "convolving (" << inch << "->" << irch << ")..."; cout.flush ();
        
        clock_t tic = clock ();
        switch (pattern) {
            case 11: {
                //read response
                irinput.read (imp0, irsamples);
                memset (imp1, 0, irsamples * sizeof (float));
                
                WavOutFile output (argv[3], sr, 16, 1);

                float scale = atof (argv[4]);
                BlockConv<float> convolver (
                    imp0, imp1, irsamples, N, scale);

                int blocks = convolver.blocks ();

                while (!input.eof ()) {
                    input.read (in, N);
                    convolver.process (in, out0, out1);
                    if (mix) for (int i = 0; i < N; ++i) out0[i] += in[i] * mix;
                    output.write (out0, N);
                }
                
                // sound tail
                memset (in, 0, sizeof (float) * N);
                while (--blocks) {
                    convolver.process (in, out0, out1);
                    output.write (out0, N);
                }
            }
            break;
            case 12: {
                //read response
                irinput.read (imp0, irsamples);
                for (int i = 0; i < irsamples / 2; ++i) {
                    impL[i] = imp0[i * 2];
                    impR[i] = imp0[i * 2 + 1];
                }

                WavOutFile output (argv[3], sr, 16, 2);

                float scale = atof (argv[4]);
                BlockConv<float> convolver (impL, impR, irsamples / 2, N, scale);

                int blocks = convolver.blocks ();

                while (!input.eof ()) {
                    input.read (in, N);

                    convolver.process (in, out0, out1);

                    for (int i = 0; i < N; ++i) {
                        if (mix) { out0[i] += in[i] * mix; out1[i] += in[i] * mix; }
                        //write output
                        output.write (&out0[i], 1);
                        output.write (&out1[i], 1);
                    }
                }
                
                // sound tail
                memset (in, 0, sizeof (float) * N);
                while (--blocks) {
                    convolver.process (in, out0, out1);
                    for (int i = 0; i < N; ++i) {
                        //write output
                        output.write (&out0[i], 1);
                        output.write (&out1[i], 1);
                    }
                }
            }
            break;
            case 21: {
                //read response
                irinput.read (imp0, irsamples);
                memset (imp1, 0, irsamples * sizeof (float));
                
                WavOutFile output (argv[3], sr, 16, 2);

                float scale = atof (argv[4]);
                BlockConv<float> convolverL (imp0, imp1, irsamples, N / 2, scale);
                BlockConv<float> convolverR (imp0, imp1, irsamples, N / 2, scale);
                
                int blocksL = convolverL.blocks ();

                while (!input.eof ()) {
                    input.read (in, N);
                    for (int i = 0; i < N / 2; ++i) {
                        inL[i] = in[i * 2];
                        inR[i] = in[i * 2 + 1];
                    }    
                    convolverL.process (inL, out0, out1);
                    convolverR.process (inR, out2, out3);
                    
                    for (int i = 0; i < N / 2; ++i) {
                        if (mix) { out0[i] += inL[i] * mix; out2[i] += inR[i] * mix; }
                        //write output
                        output.write (&out0[i], 1);
                        output.write (&out2[i], 1);
                    }
                }
                
                // sound tail
                memset (inL, 0, sizeof (float) * N / 2);
                memset (inR, 0, sizeof (float) * N / 2);
                while (--blocksL) { // left and right channels length is the same
                    convolverL.process (inL, out0, out1);
                    convolverR.process (inR, out2, out3);
                    for (int i = 0; i < N / 2; ++i) {
                        //write output
                        output.write (&out0[i], 1);
                        output.write (&out2[i], 1);
                    }
                }
            }
            break;        
            case 22: {
                //read response
                irinput.read (imp0, irsamples);
                for (int i = 0; i < irsamples / 2; ++i) {
                    impL[i] = imp0[i * 2];
                    impR[i] = imp0[i * 2 + 1];
                }

                memset (imp1, 0, irsamples * sizeof (float));
                WavOutFile output (argv[3], sr, 16, 2);

                float scale = atof (argv[4]);
                BlockConv<float> convolverL (impL, imp1, irsamples / 2, N / 2, scale);
                BlockConv<float> convolverR (impR, imp1, irsamples / 2, N / 2, scale);
                
                int blocksL = convolverL.blocks ();

                while (!input.eof ()) {
                    input.read (in, N);
                    for (int i = 0; i < N / 2; ++i) {
                        inL[i] = in[i * 2];
                        inR[i] = in[i * 2 + 1];
                    }    
                    convolverL.process (inL, out0, out1);
                    convolverR.process (inR, out2, out3);
                    
                    for (int i = 0; i < N / 2; ++i) {
                        if (mix) { out0[i] += inL[i] * mix; out2[i] += inR[i] * mix; }
                        //write output
                        output.write (&out0[i], 1);
                        output.write (&out2[i], 1);
                    }
                }
                
                // sound tail
                memset (inL, 0, sizeof (float) * N / 2);
                memset (inR, 0, sizeof (float) * N / 2);
                while (--blocksL) { // left and right channels length is the same
                    convolverL.process (inL, out0, out1);
                    convolverR.process (inR, out2, out3);
                    for (int i = 0; i < N / 2; ++i) {
                        //write output
                        output.write (&out0[i], 1);
                        output.write (&out2[i], 1);
                    }
                }
            }
            break;
            default:
                throw runtime_error ("unsupported combination of channels");
		}
        clock_t toc = clock ();
        
        cout << "done" << endl << endl;
        cout << "performance time: " << 
            (float) (toc - tic) / CLOCKS_PER_SEC << " sec." << endl;
            
		delete [] imp0;
        delete [] imp1;
		delete [] impL;
		delete [] impR;
		delete [] out0;
		delete [] out1;
		delete [] out2;
		delete [] out3;        
		delete [] in;
        delete [] inL;
        delete [] inR;
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

