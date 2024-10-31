// metaTags.h
// 


#ifndef METATAGS_H
#define METATAGS_H 

namespace soundmath {
	typedef struct { char a[1]; }  meta_no_tag;
	typedef struct { char a[2]; }  meta_yes_tag;
	
	struct default_tag {};

    template <typename T>
    meta_no_tag has_tag_typedef_function (...);
    
    template <typename T>
    meta_yes_tag has_tag_typedef_function (typename T::tag_type const volatile *);
    
	template <typename T>	
	struct has_tag_typedef {
		static const bool result = sizeof (
			has_tag_typedef_function<T>(0)) == sizeof (meta_yes_tag);
    };   
    
	template<>
	struct has_tag_typedef<void> {
		static const bool result = false;
	};
    
	template<bool condition, typename T1, typename T2>
	struct if_b {
		typedef T1 result;
	};
	
	template<typename T1, typename T2>
	struct if_b<false, T1, T2>    {
		typedef T2 result;
	};
	
	//! Metatemplate tags
	template <typename T, typename Default = default_tag>
	struct tag_traits {
	private:
		// default to vectors
		typedef typename if_b<has_tag_typedef<T>::result, T, Default>::result temp;		
	public:
		typedef typename temp::tag_type tag_type;
	};
}

#endif	// METATAGS_H 

// EOF
