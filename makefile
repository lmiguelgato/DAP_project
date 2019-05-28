all:
	g++ gcc_beamformer.cpp -ljack -lfftw3 -lm -lsndfile -I/usr/include/eigen3 tools/*.cpp -o dap_project
	g++ -g -O2 -I/usr/include/ ReadMicWavs.cpp -o ReadMicWavs -ljack -lsndfile