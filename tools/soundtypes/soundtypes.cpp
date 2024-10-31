// soundtypes.cpp
// 

#include "utilities.h"
#include "tools.h"

#include "WavFile.h"
#include "SoundTypes.h"

#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <ctime>

//#include "memCheck.h"

using namespace std;
using namespace soundmath;

int main (int argc, char* argv[]) {
	clock_t tic = clock ();
	
	srand (time (NULL));

	float* data = 0;
	
	cout << "[soundtypes, ver. 1.00]" << endl << endl
		 << "analysis and synthesis framework" << endl
		 << "written by carmine e. cella" << endl
		 << "(c) 2009-2013 www.soundtypes.com" << endl << endl;
	try {
		if (argc != 2) {
			throw runtime_error (
			    "syntax is 'soundtypes config_file.txt'");
		}

		// allocate and configure basic structures -----------------------------
		ifstream config (argv[1]);
		if (!config.good ()) {
			throw runtime_error ("cannot open configuration file");
		}

		ModelingSpace<float> p;
		p.configure (config);

		SoundTypes<float> st (&p);		
		string outputPosition = p.outputPath + "/" + removePath (p.inputFile);
		
		const char* nn = p.inputFile.c_str ();
		soundmath::WavInFile in (nn);

		float sr = in.getSampleRate ();
		
		if (sr != p.sr) {
			throw runtime_error ("invalid sr for input file");
		}
		
		int nchnls = in.getNumChannels ();
		int totSamp = in.getNumSamples ();

		// check parameters, read data and compute low-level features ----------
		int maxBin = p.checkParameters (totSamp);		
		saveAndPrintLog (outputPosition, p, nchnls);

		cout << "\nreading data.............."; cout.flush ();
		data = new float[totSamp * nchnls + p.N];
		memset (data, 0, sizeof (float) * (totSamp * nchnls + p.N));
		in.read (data, totSamp * nchnls); // data reading	
		cout << "done" << endl;
		
		if (p.analysisPerformance == "compute"){
			cout << "computing features........"; cout.flush ();
			st.analyse (data, maxBin, totSamp, nchnls);
			p.saveFeatures (outputPosition);
			st.finaliseAnalysis ();
			cout << "done" << endl;
		} else {
			cout << "loading features.............."; cout.flush ();
			p.loadFeatures (outputPosition);
			cout << "done" << endl;
		}
		
		// compute temporal modelings ------------------------------------------
		if (p.modeling) {
			p.computeModelings (p.frames);
			saveAndPrintModelings (outputPosition, p);
		} 

		// compute sound-types (types and probabilities) -----------------------
		if ((int) p.clusters) {
			if (p.typesPerformance == "compute") {
				cout << "computing types..........."; cout.flush ();			
				st.clusterAnalysis (data, totSamp, nchnls);
				cout << "done" << endl;
			} else {
                cout << "loading types............."; cout.flush ();			
                p.loadTypes(outputPosition);
                cout << "done" << endl;
            }

            if (p.probabilitiesPerformance == "compute") {
                cout << "computing probabilities..."; cout.flush ();			
                st.probabilityAnalysis ();
                cout << "done" << endl;
            } else {
                cout << "loading probabilities....."; cout.flush ();			
                p.loadProbabilities (outputPosition);
                cout << "done" << endl;
            }
            
			if (p.rematchTypes == "none") { // save only if not rematching
				p.saveTypes (outputPosition);
			} else {
				cout << "rematching types.........."; cout.flush ();
				ModelingSpace<float> p1;
				p1.loadTypes (p.rematchTypes);
				st.rematchTypes (p1.types, p1.centroids);
				cout << "done" << endl;
			}		
			
			if (p.mergeProb == "none") { // save only if not merging
				p.saveProbabilities  (outputPosition) ;
				
			} else {
				cout << "merging probabilities....."; cout.flush ();
				// copy to get same geometry (will be overwritten)
				//ModelingSpace<float>::Occurences rMarkov = p.hiMarkov[p.order - 1];

				//loadProbabilities (rMarkov, p);
				ModelingSpace<float> p1;
				p1.loadProbabilities (p.mergeProb);
				st.mergeProbabilities (p1.hiMarkov[p1.order - 1]); // fixme: check
				cout << "done" << endl;
			}	
			
			
			// synthesize sounds -----------------------------------------------
			if (p.decompose != 0) {
				cout << "synthesizing/writing......"; cout.flush ();
				synthesizeAndSave (outputPosition, st, p, totSamp);
				cout << "done" << endl;
			}		
						
			cout << endl;

			if (p.verbose) printAnalysis (p, st);
			printSummary (p, st);

		} 
		else {
			if (p.verbose) cout << endl << "cluster analsysis not performed" 
								<< endl;
		}
	}
	catch (exception& e) {
		cout << endl << endl << "Error: " << e.what () << endl;
	}
	catch (...) {
		cout << endl << endl << "Fatal error: unknown exception" << endl;
	}


	clock_t toc = clock ();	

	cout << (float) (toc - tic) / CLOCKS_PER_SEC << " sec(s) performance time" << endl;
	
	delete [] data;
	return 0;
}

// EOF
