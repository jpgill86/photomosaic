#include <stdio.h>
#include <stdlib.h>
#include "antipole.h"

int main() {
   int i;
   struct point x, y;
   srand(time(NULL));
   for( i = 0; i < DIM; i++ ) {
      x.vec[i] = rand() % 256;
      y.vec[i] = rand() % 256;
   }

   for( i = 0; i < DIM; i++ )
      printf("%d ", x.vec[i]);
   printf("\n");
   for( i = 0; i < DIM; i++ )
      printf("%d ", y.vec[i]);
   printf("\n");
   printf("%f\n", dist(x,y));

   return 0;
}

