#include <assert.h>     /* assert */
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


int
main() {

   printf(" ----- PHOTOMOSAIC ----- \n");

   // Create and initialize an array of ap_Points
   // with random vector data
   int i, j, n = 20;
   struct ap_Point arr[n];
   int seed = time(NULL);
   srand(seed);
   printf("DIM = %d\n", DIM);
   printf("n = %d\n", n);
   printf("seed = %d\n", seed);
   for( i = 0; i < n; i++ ) {
      arr[i].id = i;
      arr[i].vec = calloc( DIM, sizeof( VEC_TYPE ) );
      arr[i].ancestors = NULL;
      for( j = 0; j < DIM; j++ ) {
         ((VEC_TYPE*)arr[i].vec)[j] = rand() % 256;
      }
   }

   /*
   // Print the vector for each ap_Point
   for( i = 0; i < n; i++ ) {
      printf("arr[%d].vec = { ", i);
      for( j = 0; j < DIM; j++ ) {
         printf("%d ", ((VEC_TYPE*)arr[i].vec)[j]);
      }
      printf("}\n");
   }
   */

   // Dump the vectors for Mathematica
   printf("data = {");
   for( i = 0; i < n; i++ ) {
      printf("{%d,%d}", ((VEC_TYPE*)arr[i].vec)[0], ((VEC_TYPE*)arr[i].vec)[1]);
      if( i < n-1 )
         printf(",");
   }
   printf("};\n");

   /*
   // Print the distances between each ap_Point
   for( i = 0; i < n; i++ )
      for( j = i+1; j < n; j++ )
            printf("dist(arr[%d], arr[%d]) = %f\n", i, j, dist(&arr[i], &arr[j]));
   */

   // Place the ap_Points in an ap_List and find
   // the 1-median
   struct ap_List *s = NULL;
   for( i = 0; i < n; i++ )
      add_point( &s, &arr[i], 0 );
   printf("exact 1-median    id=%d\n", exact_1_median( s, dist )->id);
   printf("approx 1-median   id=%d\n", approx_1_median( s, dist )->id);

   // Find the exact antipole pair
   int ap1, ap2, id1, id2;
   struct ap_List *index;
   exact_antipoles( s, &ap1, &ap2, dist );
   index = s;
   for( i = 0; i < ap1; i++ )
      index = index->next;
   id1 = index->p->id;
   index = s;
   for( i = 0; i < ap2; i++ )
      index = index->next;
   id2 = index->p->id;
   printf("exact antipoles   id=%d and id=%d\n", id1, id2);

   // Find the approx antipole pair
   approx_antipoles( s, &ap1, &ap2, dist );
   index = s;
   for( i = 0; i < ap1; i++ )
      index = index->next;
   id1 = index->p->id;
   index = s;
   for( i = 0; i < ap2; i++ )
      index = index->next;
   id2 = index->p->id;
   printf("approx antipoles  id=%d and id=%d\n", id1, id2);


   /*
   // Copy s into t
   struct ap_List *t = copy_list( s );
   struct ap_List *i0 = s, *j0 = t;

   // Move members of t into u
   printf("members of t:\n");
   i0 = t;
   while( i0 != NULL ) {
      printf("i0 = %p\n", i0);
      i0 = i0->next;
   }
   printf("moving members of t into u\n", i);
   struct ap_List *u = NULL;
   while( set_size(t) > 0 ) {
      move_list( 0, &t, &u );
   }
   printf("members of u:\n");
   i0 = u;
   while( i0 != NULL ) {
      printf("i0 = %p\n", i0);
      i0 = i0->next;
   }
   */


   /*
   // Test for mem leaks in exact_1_median
   for( i = 0; i< 1e7; i++ )
      exact_1_median( s, dist );
   */

   /*
   // Test for mem leaks in approx_1_median
   for( i = 0; i < 1e7; i++ )
      approx_1_median( s, dist );
   */

   /*
   // Test for mem leaks in exact_antipoles
   for( i = 0; i < 1e7; i++ )
      exact_antipoles( s, &ap1, &ap2, dist );
   */

   /*
   // Test for mem leaks in approx_antipoles
   for( i = 0; i < 1e7; i++ )
      approx_antipoles( s, &ap1, &ap2, dist );
   */

   /*
   // Test for mem leaks in free_lists and add_point
   for( i = 0; i < 1e7; i++ ) {
      free_list(s);
      s = NULL;
      for( j = 0; j < n; j++ )
         add_point( &s, &arr[j], 0 );
   }
   */

   /*
   // Test for mem leaks in set_size and move_list
   for( i = 0; i < 1e7; i++ ) {
      while( set_size(s) > 0 )
         move_list( 0, &s, &t );
      while( set_size(t) > 0 )
         move_list( 0, &t, &s );
   }
   */

   printf(" ----------------------- \n");

   return 0;
}

