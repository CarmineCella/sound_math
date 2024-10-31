// RandomDesigner.h
// 

#ifndef RANDOMDESIGNER_H
#define RANDOMDESIGNER_H 


namespace soundmath {
	//! A spectrum designer based on random selection
	template <typename T>
	class RandomDesigner : public AbstractDesigner <T> {
	public:
		virtual T design (T f0, T dev, int part, T& amp) {
			T freq = f0;
			while (part--) freq += ((T) rand () / RAND_MAX * freq * dev);
			amp = 1. / freq;	
			
			return freq;
		}
		const char* type () const {
			return "RandomDesigner";
		}
		
	};
}

#endif	// RANDOMDESIGNER_H 

// EOF
