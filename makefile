all:
	g++ gcc_beamformer.cpp -ljack -lfftw3 -lm -lsndfile -I/usr/include/eigen3 tools/*.cpp -o gcc_beamformer
	g++ beamformer.cpp -ljack -lfftw3 -lm -lsndfile -I/usr/include/eigen3 tools/*.cpp -o beamformer
	g++ -g -O2 -I/usr/include/ ReadMicWavs.cpp -o ReadMicWavs -ljack -lsndfile