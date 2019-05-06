all:
	g++ main.cpp -ljack -lfftw3 -lm -lsndfile -lsamplerate -I/usr/include/eigen3 -o pda_project