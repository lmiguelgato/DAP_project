CPP		:= g++
CPPFLAGS:= -g -O2 -Wall

BIN		:= bin
SRC		:= src
INCLUDE	:= include
LIB		:= lib

LIBRARIES	:= -ljack -lfftw3 -lm -lsndfile

EXECUTABLE1	:= gcc_beamformer
EXECUTABLE2	:= gcc_beamformer_offline
EXECUTABLE3	:= ReadMicWavs

all: $(BIN)/$(EXECUTABLE1) $(BIN)/$(EXECUTABLE2) $(BIN)/$(EXECUTABLE3)

run: clean all
	./$(BIN)/$(EXECUTABLE1)
	./$(BIN)/$(EXECUTABLE2)
	./$(BIN)/$(EXECUTABLE3)

$(BIN)/$(EXECUTABLE1): $(SRC)/gcc_beamformer.cpp $(SRC)/tools/*.cpp
	$(CPP) $(CPPFLAGS) -I$(INCLUDE) -I$(INCLUDE)/Eigen -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(EXECUTABLE2): $(SRC)/gcc_beamformer_offline.cpp $(SRC)/tools/*.cpp
	$(CPP) $(CPPFLAGS) -I$(INCLUDE) -I$(INCLUDE)/Eigen -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(EXECUTABLE3): $(SRC)/jack_interface/ReadMicWavs.cpp
	$(CPP) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

clean:
	-rm -f $(BIN)/$(EXECUTABLE1)
	-rm -f $(BIN)/$(EXECUTABLE2)
	-rm -f $(BIN)/$(EXECUTABLE3)