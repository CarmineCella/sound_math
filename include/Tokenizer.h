// Tokenizer.h
//

#ifndef TOKENIZER_H
#define TOKENIZER_H

#include <string>
#include <vector>
#include <iterator>

namespace soundmath {
	template <typename Ch,
		typename Traits = std::char_traits<Ch>,
		typename Alloc = std::allocator<Ch> >
	//! Templated string tokenizer; internal usage.
	class Tokenizer {
	public:
		typedef std::basic_string<Ch, Traits, Alloc> T;
	
		// Generic types
		typedef T& reference;
		typedef const T& const_reference;
		typedef typename std::vector<T>::iterator iterator;
		typedef typename std::vector<T>::const_iterator const_iterator;
		typedef std::reverse_iterator<iterator> reverse_iterator;
		typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
		typedef typename std::vector<T>::size_type size_type;
		typedef T value_type;
	
		Tokenizer () {}
		//! Initialize a tokenizer delimited by one character
		explicit Tokenizer (const T& init, Ch delim) { split(init, delim); }
		//! Initialize a tokenizer delimited by character selection
		explicit Tokenizer (const T& init, const T& delim) { split(init, delim); }
	
		// Assign a new base string
		void assign (const T& init, Ch delim) { tokens.clear(); split(init, delim); }
		void assign (const T& init, const T& delim) { tokens.clear(); split(init, delim); }
	
		// Iteration
		iterator begin () { return tokens.begin (); }
		const_iterator begin() const { return tokens.begin (); }
		iterator end () { return tokens.end (); }
		const_iterator end () const { return tokens.end (); }
		reverse_iterator rbegin () { return tokens.rbegin (); }
		const_reverse_iterator rbegin () const { return tokens.rbegin (); }
		reverse_iterator rend () { return tokens.rend (); }
		const_reverse_iterator rend () const { return tokens.rend (); }
	
		// Capacity
		size_type size () const { return tokens.size (); }
		size_type empty () const { return tokens.empty (); }
	
		// Element access
		reference operator[](size_type index) { return tokens[index]; }
		const_reference operator[](size_type index) const { return tokens[index]; }
		reference at (size_type index) { return tokens.at(index); }
		const_reference at (size_type index) const { return tokens.at (index); }
		private:
		std::vector<T> tokens;
	
		void split (T base, Ch delim) { split(base, T(1, delim)); }
		void split (T base, const T& delim);
	};
	
	template <typename Ch, typename Traits, typename Alloc>
	void Tokenizer<Ch, Traits, Alloc>::split (T base, const T& delim) {
		typename T::size_type first = 0, last;
	
		while (first != T::npos) {
			last = base.find_first_of (delim, first);
			tokens.push_back (base.substr (first, last - first));
			first = base.find_first_not_of (delim, last);
		}
	}
	
	typedef Tokenizer<char> tokenizer;
	typedef Tokenizer<wchar_t> wtokenizer;
}

#endif // TOKENIZER_H

// EOF

