// Typelist.h
//


#ifndef TYPELIST_H
#define TYPELIST_H

namespace soundmath {
	class NullType;
	

	template <class T, class U>
	//! Basic metaprogramming tool to create lists of types.
	struct Typelist {
		typedef T Head;
		typedef U Tail;
	};
	
	
	template <class TList>
	//! Compute the lenght of a Typelist
	struct Length;
	template <> struct Length <NullType> {
		enum { value = 0 };
	};
	template <class T, class U>
	struct Length <Typelist <T, U> > {
		enum { value = 1 + Length <U>::value };
	};
	
	template <class TList, unsigned int index>
	//! Return the type at the requested position
	struct TypeAt;
	template <class Head, class Tail, unsigned int i>
	struct TypeAt <Typelist <Head, Tail>, i> {
		typedef typename TypeAt <Tail, i - 1>::Result Result;
	};
	template <class Head, class Tail>
	struct TypeAt <Typelist <Head, Tail>, 0> {
		typedef Head Result;
	};
	
	
	template <class TList, class T>
	//! Return the index position of the specified type
	struct IndexOf;
	
	template <class T>
	struct IndexOf <NullType, T> {
		enum { value = -1 };
	};
	
	template <class T, class Tail>
	struct IndexOf <Typelist<T, Tail>, T> {
		enum { value = 0 };
	};
	
	template <class Head, class Tail, class T>
	struct IndexOf <Typelist <Head, Tail>, T> {
	private:
		enum { temp = IndexOf <Tail, T>::value };
	public:
		enum { value = (temp == -1 ? -1 : 1 + temp) };
	};

	
	template <typename Iterator, typename TList, int N>
	//! Copy a typelist to a generic container using an output iterator
	struct Copy {
		static Iterator iterate (Iterator a) {
			a = Copy <Iterator, TList, N - 1>::iterate (a);
			*a = typeid (typename metatools::TypeAt <TList, N - 1>::Result);
			return ++a;
		}
	}; 

	template <typename Iterator, typename TList> 
	struct Copy <Iterator, TList, 0> {
		static Iterator iterate (Iterator a) {
			return a;
		}
	};	
			
	
	
	template <
		typename T1	 = NullType, typename T2  = NullType, typename T3  = NullType,
		typename T4	 = NullType, typename T5  = NullType, typename T6  = NullType,
		typename T7	 = NullType, typename T8  = NullType, typename T9  = NullType,
		typename T10 = NullType, typename T11 = NullType, typename T12 = NullType,
		typename T13 = NullType, typename T14 = NullType, typename T15 = NullType,
		typename T16 = NullType, typename T17 = NullType, typename T18 = NullType>
	//! Basic metaprogramming tool to create Typelists
	//! with different number of arguments
	struct MakeTypelist {
	private:
		typedef typename MakeTypelist <
			T2 , T3 , T4 , 
			T5 , T6 , T7 , 
			T8 , T9 , T10, 
			T11, T12, T13,
			T14, T15, T16, 
			T17, T18>::Result TailResult;
	
	public:
		typedef Typelist <T1, TailResult> Result;
	};
	
	template <>
	struct MakeTypelist <> {
		typedef NullType Result;
	};
}

#define TYPELIST_0 NullType
#define TYPELIST_1(T1) metatools::Typelist<T1, metatools::NullType>
#define TYPELIST_2(T1, T2) metatools::Typelist<T1, TYPELIST_1(T2) >
#define TYPELIST_3(T1, T2, T3) metatools::Typelist<T1, TYPELIST_2(T2, T3) >
#define TYPELIST_4(T1, T2, T3, T4) metatools::Typelist<T1, TYPELIST_3(T2, T3, T4) >
#define TYPELIST_5(T1, T2, T3, T4, T5) metatools::Typelist<T1, TYPELIST_4(T2, T3, T4, T5) >

#endif	// TYPELIST_H

// EOF

