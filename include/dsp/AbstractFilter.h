// AbstractFilter.h
// 

#ifndef ABSTRACTFILTER_H
#define ABSTRACTFILTER_H 

namespace soundmath {
    //! Base class for all processing units (filters)
	template <typename T>
	class AbstractFilter {
	private:
		AbstractFilter& operator= (AbstractFilter&);
		AbstractFilter (const AbstractFilter&);
	public: 
		AbstractFilter () {
		}
		virtual ~AbstractFilter () {
		}
		virtual void process (const T* input, T* output, int size) = 0;
	};
}
#endif	// ABSTRACTFILTER_H 

// EOF
