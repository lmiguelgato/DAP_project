#ifndef ANGLE_MANIPULATION
#define ANGLE_MANIPULATION
#include <math.h>

double state2angle (double* state);
void angle2state (double angle, double* state);
double angleRedundancy (double*, double*, double);
void angleTranslation (double *);

#endif