// scoregen.cpp
// 

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const char* methods[] = {"geom", "harm", "rand"};

double frand (double min, double max) {
	double f = 0;
	short r = (short) (rand ());
	f = fabs ((double) r / 32768);
	f *= (max - min);
	f += min;
	return f;
}

int main (int argc, char* argv[]) {
	srand (time (NULL));
	try {
		if (argc != 7) {
			throw runtime_error ("syntax is scoregen score.txt num f0min f0max pmin pmax");
		}
		
		ofstream out (argv[1]);
		if (!out.good ()) {
			throw runtime_error ("cannot create score file");
		}
		int events = atoi (argv[2]);
		if (events < 1) {
			throw runtime_error ("invalid number of events");
		}
		
		double maxstart = 30;
		double mindur = 40;
		double maxdur = 90;
		
		double fmin = atof (argv[3]);
		double fmax = atof (argv[4]);
		int pmin = atof (argv[5]);
		int pmax = atof (argv[6]);
		int thickmin = 2;
		int thickmax = 5;
		double spreadmin = .001;
		double spreadmax = .05;
		double cmin1 = .6;
		double cmax1 = 1.8;
		
		double cmin[3];
		double cmax[3];
		
		cmin[0]= 1;
		cmax[0] = 1.8;

		cmin[1]= .5;
		cmax[1] = 1.7;
		
		cmin[2]= .1;
		cmax[2] = .99;
		
		double start = 0;
		for (int i = 0; i < events; ++i) {
			int m = (int) frand (0, 3);
			out << methods[m] << " ";
			out << start << " ";
			start = frand (start, start + maxstart);
			out << frand (mindur, maxdur) << " ";
			out << frand (fmin, fmax) << " ";
			out << (int) frand (pmin, pmax) << " ";
			out << (int) frand (thickmin, thickmax) << " ";
			out << frand (spreadmin, spreadmax) << " ";
			double c = frand (cmin[m], cmax[m]);
			out << c << " ";
			double slo = frand (.1, .9);
			out << slo << endl;
		}
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

