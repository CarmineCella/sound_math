// AbstractEstimator.h
// 

#ifndef ABSTRACTESTIMATOR_H
#define ABSTRACTESTIMATOR_H 

namespace soundmath {
    //! Base class for all pitch detectors
	template <typename T>
	class AbstractEstimator {
	private:
		AbstractEstimator& operator= (AbstractEstimator&);
		AbstractEstimator (const AbstractEstimator&);
	public: 
		AbstractEstimator () {
		}
		virtual ~AbstractEstimator () {
		}
		
		virtual T process (const T* in, T& amp) = 0;
	};
}
#endif	// ABSTRACTESTIMATOR_H 

// EOF
