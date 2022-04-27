#modied from htslib makefile

CC  ?= gcc -std=c++11
CXX ?= g++

PROGRAMS = QualConvert LenConvert MisConvert
#PROGRAMS_MISC = XXX
all: $(PROGRAMS) #$(PROGRAMS_MISC)

QualConvert: ReadQualConverter.cpp
	$(CXX) ReadQualConverter.cpp -D __WITH_MAIN__ -lz -o QualConvert

LenConvert: MapDMGlen.cpp
	$(CXX) MapDMGlen.cpp -D __WITH_MAIN__ -lz -o LenConvert

MisConvert: MapDMGlen.cpp
	$(CXX) MapDMGmis.cpp -D __WITH_MAIN__ -lz -o MisConvert

clean:	
	rm  -f *.o $(PROGRAMS) *~ 

#rm  -f *.o $(PROGRAMS) $(PROGRAMS_MISC) *~
