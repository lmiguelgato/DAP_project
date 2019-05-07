all:
	g++ main.cpp -ljack -lfftw3 -lm -I/usr/include/eigen3 tools/*.cpp -o dap_project
	g++ -g -O2 -I/usr/include/ ReadMicWavs.cpp -o ReadMicWavs -ljack -lsndfile