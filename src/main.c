#include <stdio.h>      /* printf */
#include <stdlib.h>     /* rand, time */
#include "antipole.h"

int main() {
   printf(" ----- PHOTOMOSAIC ----- \n");

   int i;
   struct ap_Point x, y;
   srand(time(NULL));
   for( i = 0; i < DIM; i++ ) {
      x.vector[i] = rand() % 256;
      y.vector[i] = rand() % 256;
   }

   for( i = 0; i < DIM; i++ )
      printf("%d ", x.vector[i]);
   printf("\n");
   for( i = 0; i < DIM; i++ )
      printf("%d ", y.vector[i]);
   printf("\n");
   printf("%f\n", dist(&x,&y));

   return 0;
}

