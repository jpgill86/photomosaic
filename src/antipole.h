#ifndef ANTIPOLE_H
#define ANTIPOLE_H 

#include <stdint.h>

#define DIM 3

struct point {
   uint8_t vec[DIM];
   int id;
};

double dist( struct point p1, struct point p2 );

#endif /* ANTIPOLE_H */

