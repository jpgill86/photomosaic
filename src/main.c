#include <math.h>       /* sqrt, pow */
#include <stdint.h>     /* uint8_t */
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* rand, time */
#include "antipole.h"


const int DIM = 2;         /* dimensionality of the mean RGB data */
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

   // Create and initialize an array of ap_Points
   // with random vector data
   int i, j, n = 5;
   struct ap_Point arr[n];
   srand(time(NULL));
   for( i = 0; i < n; i++ ) {
      arr[i].id = i;
      arr[i].vec = calloc( DIM, sizeof( VEC_TYPE ) );
      for( j = 0; j < DIM; j++ ) {
         ((VEC_TYPE*)arr[i].vec)[j] = rand() % 256;
      }
   }

   // Print the vector for each ap_Point
   for( i = 0; i < n; i++ ) {
      printf("arr[%d].vec = { ", i);
      for( j = 0; j < DIM; j++ ) {
         printf("%d ", ((VEC_TYPE*)arr[i].vec)[j]);
      }
      printf("}\n");
   }

   // Print the distances between each ap_Point
   for( i = 0; i < n; i++ )
      for( j = i+1; j < n; j++ )
            printf("dist(arr[%d], arr[%d]) = %f\n", i, j, dist(arr+i, arr+j));

   // Place the ap_Points in an ap_List and find
   // the 1-median
   struct ap_List *s = NULL;
   for( i = 0; i < n; i++ )
      add_point( &s, &arr[i], 0 );
   printf("1-median is id=%d\n", find_1_median( s, dist )->id);

   // Test for mem leaks
   for( i = 0; i < 10000000; i++ ) {
      free_list(s);
      s = NULL;
      for( j = 0; j < n; j++ )
         add_point( &s, &arr[j], 0 );
      find_1_median( s, dist );
   }

   printf(" ----------------------- \n");

   return 0;
}

