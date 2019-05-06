all:
	g++ main.cpp -ljack -lfftw3 -lm -lsndfile -lsamplerate -I/usr/include/eigen3 -o pda_project
	g++ -g -O2 -I/usr/include/ ReadMicWavs.cpp -o ReadMicWavs -ljack -lsndfile