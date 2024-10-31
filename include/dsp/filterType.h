// filterType.h
// 

#ifndef FILTERTYPE_H
#define FILTERTYPE_H 

namespace soundmath {
	//! Types of supported filters
	enum filterType {
		LOPASS = 0,
		HIPASS = 1,
		BANDPASS1 = 2,
		BANDPASS2 = 3,
		NOTCH = 4,
		ALLPASS = 5,
		PEAKING = 6,
		LOSHELF = 7,
		HISHELF = 8
	};
}

#endif	// FILTERTYPE_H 

// EOF
