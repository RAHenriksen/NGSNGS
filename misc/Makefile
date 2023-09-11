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


# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
LIBS += $(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
LIBS += -lhts
endif

PROGRAMS = QualConvert LenConvert MisConvert MetaMisConvert RandRef RandVar VCFtoBED

default: all

all: $(PROGRAMS)

LenConvert: MapDMGlen.cpp
	$(CXX) $< -D __WITH_MAIN__ $(CXXFLAGS) -o $@ $(LIBS)
MisConvert: MapDMGmis.cpp
	$(CXX) $< -D __WITH_MAIN__ $(CXXFLAGS) -o $@ $(LIBS)
MetaMisConvert: MetaDMGdeam.cpp
	$(CXX) $< -D __WITH_MAIN__ $(CXXFLAGS) -o $@ $(LIBS)
RandRef: RandRef.cpp
	$(CXX) $< -D __WITH_MAIN__ ../mrand.o ../fasta_sampler.o ../RandSampling.o $(CXXFLAGS) -o $@ $(LIBS)
RandVar: RandVariant.cpp
	$(CXX) $< -D __WITH_MAIN__ ../mrand.o ../fasta_sampler.o ../RandSampling.o $(CXXFLAGS) -o $@ $(LIBS)
QualConvert: ReadQualConverter.cpp
	$(CXX) $< -D __WITH_MAIN__ $(CXXFLAGS) -o $@ $(LIBS)
VCFtoBED: VCFtoBED.cpp
	$(CXX) $< -D __WITH_MAIN__ $(CXXFLAGS) -o $@ $(LIBS)
clean:	
	rm -f $(PROGRAMS) *~ 

.PHONY: all clean