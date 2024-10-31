// constants.h
// 

#ifndef CONSTANTS_H
#define CONSTANTS_H 

#include <cmath>

namespace soundmath {
    extern const double EPS = .0000000000001;
    extern const double EPS2 = 2.22044604925031e-16;
	extern const double MINIMUM_REASONABLE_THRESHOLD = 0.005;
	
    extern const double PI = 4. * atan ((double) 1.);
    extern const double TWOPI = 8. * atan ((double) 1.);
	
    extern const double LOG2 = (double) log ((double) 2.);
	extern const double LOG2OF10 = 3.32192809488736234787;
    extern const double SQRT2 = sqrt (2.);
 
    extern const double P0 = 1 / (std::pow (2.0, 32.) / 2 - 1);
	extern const double A1 = ((double) -1. / 720.);
	extern const double A2 = ((double) -1. / 360.);	
	
	extern const double RMAX = 0x7fffffff;
}
#endif	// CONSTANTS_H 

// EOF
