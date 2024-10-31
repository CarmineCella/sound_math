#sounder makefile

SOURCES=sounder.cpp 

HEADERS=

EXECUTABLE=sounder

CC=g++

CXXFLAGS=-I../../include/dsp -I../../include/inout -I../../include/helpers -I../../include/math -I../../include/metatools -c -Wall -Wno-write-strings -arch i386 -arch x86_64 -mmacosx-version-min=10.7 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/

LIBDIRNAME=

LIBS=

LDFLAGS= -arch i386 -arch x86_64 -mmacosx-version-min=10.7 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/

OBJECTS= $(SOURCES:.cpp=.o)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LIBDIRNAME) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CXXFLAGS) $< -o $@

release: CXXFLAGS+= -O3
release: all

debug: CXXFLAGS+= -g3 -DDEBUG
debug: all

all: $(SOURCES) $(EXECUTABLE)

clean:
	rm -f -r *.o $(EXECUTABLE).app *~ $(EXECUTABLE)

