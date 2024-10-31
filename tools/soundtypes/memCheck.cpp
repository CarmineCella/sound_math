// memCheck.cpp
//

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

#undef new

namespace {
	struct Info {
		void* ptr;
		long amount;
		const char* file;
		long line;
	};

	const size_t MAXPTRS = 10000u;
	Info memMap[MAXPTRS];
	size_t nptrs = 0;
	FILE* logger = 0;
	long totalMemory = 0;
	char tempBuffer[100];

	int findPtr (void* p) {
		for (unsigned int i = 0; i < nptrs; ++i)
			if (memMap[i].ptr == p)
				return i;
		return -1;
	}

	void delPtr (void* p) {
		int pos = findPtr (p);
		assert (p >= 0);
	
		for (size_t i = pos; i < nptrs - 1; ++i)
			memMap[i] = memMap[i + 1];
		
		--nptrs;
	}

	struct Sentinel {	// must be the last memory object!!
		Sentinel () {
			logger = fopen ("memCheck.log", "w+b");
			assert (logger);
			
			sprintf (tempBuffer, "[memCheck session - %s, %s]\n\n", __DATE__, __TIME__);
			fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
			printf ("%s", tempBuffer);
			
			#ifdef _WIN32
			OutputDebugString (tempBuffer);
			#endif
		}
		~Sentinel () {
			if (nptrs > 0) {
				sprintf (tempBuffer, "\nLeaked memory:\n");
				fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
				printf ("%s", tempBuffer);
				
				#ifdef _WIN32
				OutputDebugString (tempBuffer);
				#endif

				
				for (size_t i = 0; i < nptrs; ++i) {
					sprintf (tempBuffer, "\t%p %ld bytes (file: %s, line %ld)\n",
						memMap[i].ptr, memMap[i].amount, memMap[i].file, memMap[i].line);
					fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
					printf ("%s", tempBuffer);

					#ifdef _WIN32
					OutputDebugString (tempBuffer);
					#endif

					totalMemory += memMap[i].amount;
				}
				
				sprintf (tempBuffer, "\t%ld bytes totally\n", totalMemory);
				fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
				printf ("%s", tempBuffer);

				#ifdef _WIN32
				OutputDebugString (tempBuffer);
				#endif
			}
			else {
				sprintf (tempBuffer, "\nNo user memory leaks!\n");
				fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
				printf ("%s", tempBuffer);
				
				#ifdef _WIN32
				OutputDebugString (tempBuffer);
				#endif				
			}
		}
	};

	Sentinel s;
}

void* operator new (size_t siz, const char* file, long line) {
	void* p = malloc (siz);
	if (nptrs == MAXPTRS) {
		sprintf (tempBuffer, "\nError: memory map too small (increase MAXPTRS)\n");
		fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
		printf ("%s", tempBuffer);
				
		#ifdef _WIN32
		OutputDebugString (tempBuffer);
		#endif
			
		exit (1);
	}

	memMap[nptrs].ptr = p;
	memMap[nptrs].amount = siz;
	memMap[nptrs].file = file;
	memMap[nptrs].line = line;
	++nptrs;

	sprintf (tempBuffer, "Allocated %u bytes at address %p (file: %s, line: %ld)\n", 
		siz, p, file, line);
	fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
	printf ("%s", tempBuffer);
	
	#ifdef _WIN32
	OutputDebugString (tempBuffer);
	#endif

	return p;
}

void* operator new[] (size_t siz, const char* file, long line) {
	return operator new (siz, file, line);	
}

void operator delete (void* p) {
	if (findPtr (p) >= 0) {
		free (p);
		assert (nptrs > 0);
		delPtr (p);

		sprintf (tempBuffer, "Deleted memory at address %p\n", p);
		fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
		printf ("%s", tempBuffer);
				
		#ifdef _WIN32
		OutputDebugString (tempBuffer);
		#endif
	}

	else if (!p) {
		sprintf (tempBuffer, "Attempt to delete unknown pointer: %p\n", p);
		fwrite (tempBuffer, strlen (tempBuffer), 1, logger);
		printf ("%s", tempBuffer);
				
		#ifdef _WIN32
		OutputDebugString (tempBuffer);
		#endif
	}
}

void operator delete[](void* p) {
	operator delete(p);
}

// EOF

