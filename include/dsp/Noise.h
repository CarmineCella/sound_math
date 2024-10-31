// Noise.h
// 

#ifndef NOISE_H
#define NOISE_H 

namespace soundmath {
	class Noise {
	private:
		Noise& operator= (Noise&);
		Noise (const Noise&);
	public: 
		Noise (int seed = 1) : m_seed (seed) {}
		
		virtual ~Noise () {
		}
		
		virtual ~Noise () {}
		void process (T* output, int N) {
			for (int i = 0; i < N; ++i) {
				m_seed *= 16807;
				output[i] = (T) m_seed * 4.6566129e-010f;
			}
		}
	private:
		int m_seed;
	};
}

#endif	// NOISE_H 

// EOF
