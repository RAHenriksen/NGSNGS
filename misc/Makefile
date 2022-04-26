#modied from htslib makefile

CC  ?= gcc -std=c++11
CXX ?= g++

PROGRAMS = QualConvert
#PROGRAMS_MISC = XXX
all: $(PROGRAMS) #$(PROGRAMS_MISC)

QualConvert: ReadQualConverter.cpp
	$(CXX) ReadQualConverter.cpp -D __WITH_MAIN__ -lz -o QualConvert

clean:	
	rm  -f *.o $(PROGRAMS) *~ 

#rm  -f *.o $(PROGRAMS) $(PROGRAMS_MISC) *~
