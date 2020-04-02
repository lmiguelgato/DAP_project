#ifndef ANGLE_MANIPULATION
#define ANGLE_MANIPULATION
#include <math.h>
#define DEG2RAD 0.017453292519943		// useful to convert from degrees to radians
#define RAD2DEG 57.295779513082323f		// useful to convert from radians to degrees

double state2angle (double* state);
void angle2state (double angle, double* state);
double angleRedundancy (double*, double*, double);
void angleTranslation (double *, float *);
void angleTranslation (double *);

#endif