#modied from htslib makefile
#FLAGS=-O3 -std=c++11
GITVERSION:= $(shell git log -1 --pretty='%h')

FLAGS=-ggdb -std=c++11 -DGITVERSION="\"${GITVERSION}\""

CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

all: ngsngs #qualconvert

.PHONY: all clean test

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif


-include $(OBJ:.o=.d)

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp >$*.d

ngsngs: $(OBJ)
	$(CXX) $(FLAGS)  -o ngsngs *.o $(HTS_LIBDIR) -lz -llzma -lbz2 -lpthread -lcurl

else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

ngsngs: $(OBJ)
	$(CXX) $(FLAGS)  -o ngsngs *.o -lz -llzma -lbz2 -lpthread -lcurl -lhts

#qualconvert: $(CXX) $(FLAGS) -o qualconvert ReadQualConverter.cpp -lz -llzma -lbz2 -lpthread -lcurl -lhts -lgsl -lgslcblas

endif

clean:	
	rm  -f ngsngs *.o *.d

test:
	echo "Subprograms is being tested";
	cd test; sh testAll.sh;
