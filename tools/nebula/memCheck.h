// memCheck.h
//

#ifndef MEMCHECK_H
#define MEMCHECK_H

#include <cstddef>  // for size_t

void* operator new (std::size_t, const char*, long);
void* operator new[] (std::size_t, const char*, long);
#define new new (__FILE__, __LINE__)

#endif	// MEMCHECK_H

// EOF

