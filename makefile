all:
	g++ main.cpp -ljack -lfftw3 -lm -lsndfile -lsamplerate -o pda_project