all:
	gcc -std=gnu99 -o pda_project main.c -ljack -lfftw3 -lm