// ParamSmooth.h
// 


#ifndef PARAMSMOOTH_H
#define PARAMSMOOTH_H 

namespace soundmath {
	template <typename T>
	//! Parameter smoothing filter: sample-by-sample version
	class ParamSmooth {
	public:
	    ParamSmooth (double smooth = .5) {
	    	smoothing (smooth);
	    	reset ();
	    };
	    virtual ~ParamSmooth () {}
	    inline T process (T in) { m_z = (in * m_b) + (m_z * m_a); return m_z; }
	    void smoothing (double s) {
	    	m_a = s;
	    	m_b = 1.f - m_a;
	    }
	    void reset () {
	    	m_z = 0;
	    }
	private:
	    T m_a, m_b, m_z;
	};
}
#endif	// PARAMSMOOTH_H 

// EOF
