// TypeInfo.h
//

#ifndef TYPEINFO_H
#define TYPEINFO_H

#include <typeinfo>
#include <cassert>

namespace soundmath {
	//! This class is a std::type_info wrapper to let the user to create
	//! collections.
	class TypeInfo {
	public:
		TypeInfo ();
		TypeInfo (const std::type_info&);
		const std::type_info& get () const;
		bool before (const TypeInfo& rhs) const;
		const char* name () const;
	private:
		const std::type_info* m_info;
	};
	
	inline TypeInfo::TypeInfo () {
		class Nil {};
		m_info = &typeid (Nil);
		assert (m_info);
	}
	
	inline TypeInfo::TypeInfo (const std::type_info& ti)
		: m_info(&ti) {
		assert (m_info); 
	}
	
	inline bool TypeInfo::before (const TypeInfo& rhs) const {
		assert (m_info);
		return m_info->before (*rhs.m_info) != 0;
	}
	
	inline const std::type_info& TypeInfo::get () const {
		assert (m_info);
		return *m_info;
	}
	
	inline const char* TypeInfo::name () const {
		assert (m_info);
		return m_info->name ();
	}
	
	inline bool operator== (const TypeInfo& lhs, const TypeInfo& rhs) { 
		return (lhs.get () == rhs.get ()) != 0;
	}
	
	inline bool operator< (const TypeInfo& lhs, const TypeInfo& rhs) { 
		return lhs.before(rhs);
	}
	
	inline bool operator!= (const TypeInfo& lhs, const TypeInfo& rhs) {
		return !(lhs == rhs);
	}    
	
	inline bool operator> (const TypeInfo& lhs, const TypeInfo& rhs) {
		return rhs < lhs;
	}
	
	inline bool operator<= (const TypeInfo& lhs, const TypeInfo& rhs) {
		return !(lhs > rhs);
	}
	 
	inline bool operator>= (const TypeInfo& lhs, const TypeInfo& rhs) { 
		return !(lhs < rhs);
	}
}

#endif	// TYPEINFO_H

// EOF

