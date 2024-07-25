##made by rasmus H 18juli 2022. Its summer!!
CC  ?= gcc
CXX ?= g++

ALWAYS_FLAGS=-std=c++11

FLAGS=-O3

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl

## add isins crypto trick, copied from ANGSD makefile, which is based on something from samtools
CRYPTO_TRY=$(shell echo 'int main(){}'|g++ -x c++ - -lcrypto 1>/dev/null 2>/dev/null; echo $$?)
ifeq "$(CRYPTO_TRY)" "0"
$(info Crypto library is available to link; adding -lcrypto to LIBS)
LIBS += -lcrypto
else
$(info Crypto library is not available to link; will not use -lcrypto)
endif


CFLAGS += $(FLAGS) $(ALWAYS_FLAGS) 
CXXFLAGS += $(FLAGS) $(ALWAYS_FLAGS) 

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

all: ngsngs

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always)
PACKAGE_RELEASE := $(shell git describe --tags --abbrev=0)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define NGSNGS_RELEASE "$(PACKAGE_RELEASE)"' > $@
	echo '#define NGSNGS_VERSION "$(PACKAGE_VERSION)"' >> $@

.PHONY: all clean test

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
LIBS += -lhts
endif

-include $(OBJ:.o=.d)

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp >$*.d

ngsngs: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o ngsngs *.o $(HTS_LIBDIR) $(LIBS)

amplicon:	
	$(CXX) $(FLAGS) Amplicon.cpp -o amplicon mrand.o Briggs.o NtSubModels.o add_indels.o sample_qscores.o HelpPage.o Amplicon_cli.o NGSNGS_misc.o RandSampling.o $(HTS_LIBDIR) $(LIBS) -D __WITH_MAIN_AMPLICON__ 

else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

ngsngs: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o ngsngs *.o $(LIBS)

amplicon:	
	$(CXX) $(FLAGS) Amplicon.cpp -o amplicon mrand.o Briggs.o NtSubModels.o add_indels.o sample_qscores.o HelpPage.o Amplicon_cli.o NGSNGS_misc.o RandSampling.o $(LIBS) -D __WITH_MAIN_AMPLICON__ 

#$(CXX) Amplicon.cpp $(FLAGS) -o amplicon $(LIBS) mrand.o Briggs.o NtSubModels.o add_indels.o HelpPage.o Amplicon_cli.o NGSNGS_misc.o -D __WITH_MAIN__ 

endif

clean:	
	rm -f ngsngs amplicon *.o *.d version.h a.out
	rm -f test/*.fa test/*.fq test/*.sam test/*.txt

test:
	echo "Subprograms is being tested"
	cd test; ./testAll.sh

force:
