#ifndef ANGLE_MANIPULATION
#define ANGLE_MANIPULATION
#ifndef DEG2RAD
#define DEG2RAD 0.017453292519943		// useful to convert from degrees to radians
#endif
#ifndef RAD2DEG
#define RAD2DEG 57.295779513082323f		// useful to convert from radians to degrees
#endif

#include <math.h>

double state2angle (double* state);
void angle2state (double angle, double* state);
double angleRedundancy (double*, double);
void angleTranslation (double *);

#endif