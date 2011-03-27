#include <math.h>       /* sqrt, pow */
#include <stdint.h>     /* uint8_t */
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* rand, time */
#include "antipole.h"


const int DIM = 3;         /* dimensionality of the mean RGB data */
typedef uint8_t VEC_TYPE;  /* data type of the mean RGB data */

// Calculate the Euclidian distance between two points
double
dist( struct ap_Point *p1, struct ap_Point *p2 ) {
   int i;
   double sum = 0;
   for( i = 0; i < DIM; i++ )
      sum += pow( ((VEC_TYPE*)p1->vec)[i] - ((VEC_TYPE*)p2->vec)[i] , 2);
   return sqrt( sum );
}


int main() {
   printf(" ----- PHOTOMOSAIC ----- \n");

   int i;
   struct ap_Point x, y;
   srand(time(NULL));
   x.vec = malloc( sizeof( VEC_TYPE[DIM] ) );
   y.vec = malloc( sizeof( VEC_TYPE[DIM] ) );
   for( i = 0; i < DIM; i++ ) {
      ((VEC_TYPE*)x.vec)[i] = rand() % 256;
      ((VEC_TYPE*)y.vec)[i] = rand() % 256;
   }

   for( i = 0; i < DIM; i++ )
      printf("%d ", ((VEC_TYPE*)x.vec)[i]);
   printf("\n");
   for( i = 0; i < DIM; i++ )
      printf("%d ", ((VEC_TYPE*)y.vec)[i]);
   printf("\n");
   printf("%f\n", dist(&x,&y));

   return 0;
}

