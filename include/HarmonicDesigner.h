// HarmonicDesigner.h
// 

#ifndef HARMONICDESIGNER_H
#define HARMONICDESIGNER_H 

#include "AbstractDesigner.h"

namespace soundmath {
	//! A spectrum designer based on compressed/expanded harmonic series
	template <typename T>
	class HarmonicDesigner : public AbstractDesigner <T> {
	public:
		virtual T design (T f0, T dev, int part, T& amp) {
			//T f =  f0 * part * sqrt (1. + dev * (part * part));
			T f =  f0 * pow ((T) part, dev);
			amp = 1. / f;
			return f;
		}
		const char* type () const {
			return "HarmonicDesigner";
		}
		
	};
}

#endif	// HARMONICDESIGNER_H 

// EOF
