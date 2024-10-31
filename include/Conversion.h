// Conversion.h
//

#ifndef CONVERSION_H
#define CONVERSION_H
	
namespace soundmath {
	
	template <class T, class U>
	//! Basic metaprogramming tool to check conversions between types.
	class Conversion {
		typedef char Small;
		class Big { char dummy[2]; };
		static Small Test (U);
		static Big Test (...);
		static T MakeT (); // not defined, only declared!
	public:
		enum { exists = sizeof (Test (MakeT ())) == sizeof (Small) };
		enum { exists2Way = exists &&  Conversion <U, T>::exists };
		enum { sameType = false };
	};
	template <class T>
	class Conversion <T, T> {
	public:
		enum { exists = 1, exists2Way = 1, sameType = 1 };
	};
	
	//! Metatemplate representation of compile-time errors
	template<int> struct CompileTimeError;
	template<> struct CompileTimeError<true> {};

}

#define SUPERSUBCLASS(T, U) \
	(metatools::Conversion <const U*, const T*>::exists && \
	!metatools::Conversion <const T*, const void*>::sameType)
	
#define SUPERSUBCLASS_STRICT(T, U) \
	(SUPERSUBCLASS(T, U) && \
	!metatools::Conversion <const T, const U>::sameType)

#define STATIC_CHECK(expr, msg) \
	{ metatools::CompileTimeError<((expr) != 0)> ERROR_##msg; (void)ERROR_##msg; } 

#endif	// CONVERSION_H

// EOF

