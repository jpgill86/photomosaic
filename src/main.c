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
dist( ap_Point *p1, ap_Point *p2 ) {

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
   ap_Point arr[n];
   int seed = time(NULL);
   //seed = 1301506514;
   srand(seed);
   printf("DIM = %d\n", DIM);
   printf("n = %d\n", n);
   printf("seed = %d\n", seed);
   for( i = 0; i < n; i++ ) {
      arr[i].id = i;
      arr[i].vec = calloc( DIM, sizeof( VEC_TYPE ) );
      arr[i].ancestors = NULL;
      for( j = 0; j < DIM; j++ )
         ((VEC_TYPE*)arr[i].vec)[j] = rand() % 256;
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

   /*
   // Dump the vectors for Mathematica
   printf("data = {");
   for( i = 0; i < n; i++ ) {
      printf("{%d,%d}", ((VEC_TYPE*)arr[i].vec)[0], ((VEC_TYPE*)arr[i].vec)[1]);
      if( i < n-1 )
         printf(",");
   }
   printf("};\n");
   */

   /*
   // Print the distances between each ap_Point
   for( i = 0; i < n; i++ )
      for( j = i+1; j < n; j++ )
            printf("dist(arr[%d], arr[%d]) = %f\n", i, j, dist(&arr[i], &arr[j]));
   */

   // Place the ap_Points in an ap_List
   ap_List *s = NULL;
   for( i = 0; i < n; i++ )
      add_point( &s, &arr[i], 0 );

   /*
   // Find the 1-median
   ap_Point *median;
   exact_1_median( s, &median, dist );
   printf("exact 1-median    id=%d\n", median->id);
   approx_1_median( s, &median, dist );
   printf("approx 1-median   id=%d\n", median->id);

   // Find the antipole pair
   ap_Point *antipole_a, *antipole_b;
   exact_antipoles( s, &antipole_a, &antipole_b, dist );
   printf("exact antipoles   id=%d and id=%d\n", antipole_a->id, antipole_b->id);
   approx_antipoles( s, &antipole_a, &antipole_b, dist );
   printf("approx antipoles  id=%d and id=%d\n", antipole_a->id, antipole_b->id);
   */

   // Construct a tree and search it
   printf("build tree\n");
   ap_Tree *tree = build_tree( 0, s, 256*0.05*sqrt(DIM), NULL, NULL, DIM, dist );
   printf("range search\n");
   ap_List *range = NULL;
   ap_Point query_p;
   ap_Point *query = &query_p;
   query->id = -1;
   query->vec = calloc( DIM, sizeof( VEC_TYPE ) );
   query->ancestors = NULL;
   for( i = 0; i < DIM; i++ )
      ((VEC_TYPE*)query->vec)[i] = rand() % 256;
   //printf("query = (%d,%d)\n", ((VEC_TYPE*)query->vec)[0], ((VEC_TYPE*)query->vec)[1]);
   range_search( tree, query, 50, &range, dist );
   ap_List *index = range;
   printf("query results\n");
   printf("size=%d\n", list_size(range));
   while( index != NULL ) {
      printf("match has id=%d\n", index->p->id);
      index = index->next;
   }


   /*
   // Copy s into t
   ap_List *t = copy_list( s );
   ap_List *i0 = s, *j0 = t;

   // Move members of t into u
   printf("members of t:\n");
   i0 = t;
   while( i0 != NULL ) {
      printf("i0 = %p\n", i0);
      i0 = i0->next;
   }
   printf("moving members of t into u\n", i);
   ap_List *u = NULL;
   while( list_size(t) > 0 ) {
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
      exact_1_median( s, &median, dist );
   */

   /*
   // Test for mem leaks in approx_1_median
   for( i = 0; i < 1e7; i++ )
      approx_1_median( s, &median, dist );
   */

   /*
   // Test for mem leaks in exact_antipoles
   for( i = 0; i < 1e7; i++ )
      exact_antipoles( s, &antipole_a, &antipole_b, dist );
   */

   /*
   // Test for mem leaks in approx_antipoles
   for( i = 0; i < 1e7; i++ )
      approx_antipoles( s, &antipole_a, &antipole_b, dist );
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
   // Test for mem leaks in list_size and move_list
   for( i = 0; i < 1e7; i++ ) {
      while( list_size(s) > 0 )
         move_list( 0, &s, &t );
      while( list_size(t) > 0 )
         move_list( 0, &t, &s );
   }
   */

   printf(" ----------------------- \n");

   return 0;
}

