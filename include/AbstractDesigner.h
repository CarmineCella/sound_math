// AbstractDesigner.h
// 

#ifndef ABSTRACTDESIGNER_H
#define ABSTRACTDESIGNER_H 

namespace soundmath {
    //! Base class for all spectral designers
	template <typename T>
	class AbstractDesigner {
	public:
		virtual ~AbstractDesigner () {
		}
		virtual T design (T f0, T dev, int part, T& amp) = 0;
		virtual const char* type () const  = 0;
	};
}
#endif	// ABSTRACTDESIGNER_H 

// EOF
