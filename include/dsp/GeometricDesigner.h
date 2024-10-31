// GeometricDesigner.h
// 

#ifndef GEOMETRICDESIGNER_H
#define GEOMETRICDESIGNER_H 

#include "AbstractDesigner.h"

namespace soundmath {
	//! A spectrum designer based on geometric laws
	template <typename T>
	class GeometricDesigner : public AbstractDesigner <T> {
	public:
		virtual T design (T f0, T dev, int part, T& amp) {
			if (dev <= 1) return 0;
			
			
			
			T freq = f0;
			while (part--) freq *= dev;
			amp = 1. / freq;	
			
			return freq;
		}
		const char* type () const {
			return "GeometricDesigner";
		}		
	};
}

#endif	// GEOMETRICDESIGNER_H 

// EOF
