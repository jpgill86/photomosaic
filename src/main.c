/* main.c
 *
 * Copyright (c) 2011, Jeffrey P. Gill
 *
 * This file is part of photomosaic.
 *
 * photomosaic is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * photomosaic is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with photomosaic.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>     /* assert */
#include <math.h>       /* sqrt, pow */
#include <stdint.h>     /* uint8_t */
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* rand */
#include <time.h>       /* time */
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

   printf("(* ----- PHOTOMOSAIC ----- *)\n");

   // Create and initialize an array of ap_Points
   // with random vector data
   int i, j, n = 20;
   ap_Point data[n];
   int seed = time(NULL);
   srand(seed);
   printf("dim = %d;\n", DIM);
   printf("n = %d;\n", n);
   printf("seed = %d;\n", seed);
   for( i = 0; i < n; i++ ) {
      data[i].id = i + 1;
      data[i].vec = calloc( DIM, sizeof( VEC_TYPE ) );
      data[i].ancestors = NULL;
      for( j = 0; j < DIM; j++ )
         ((VEC_TYPE*)data[i].vec)[j] = rand() % 256;
   }

#ifdef DEBUG
   // Dump the data vectors for Mathematica
   printf("data = {");
   for( i = 0; i < n; i++ ) {
      printf("{");
      for( j = 0; j < DIM; j++ ) {
         printf("%d", ((VEC_TYPE*)data[i].vec)[j]);
         if( j < DIM-1 )
            printf(",");
      }
      if( i < n-1 )
         printf("},");
      else
         printf("}");
   }
   printf("};\n");
#endif

   // Place the ap_Points in an ap_PointList
   ap_PointList *s = NULL;
   for( i = 0; i < n; i++ )
      add_point( &s, &data[i], 0 );

   /*
   // Find the 1-median
   ap_Point *median;
   exact_1_median( s, &median, dist );
   printf("exactMedian = %d;\n", median->id);
   approx_1_median( s, &median, DIM, dist );
   printf("approxMedian = %d;\n", median->id);

   // Find the antipole pair
   ap_Point *antipole_a, *antipole_b;
   exact_antipoles( s, &antipole_a, &antipole_b, dist );
   printf("exactAntipoles = {%d,%d};\n", antipole_a->id, antipole_b->id);
   approx_antipoles( s, &antipole_a, &antipole_b, DIM, dist );
   printf("approxAntipoles = {%d,%d};\n", antipole_a->id, antipole_b->id);
   */

   // Construct a tree
   printf("(* build tree *)\n");
   ap_Tree *tree = build_tree( s, 256*0.05*sqrt(DIM), NULL, NULL, DIM, dist );

   // Construct a set of query points
   int n_query = 5;
   ap_Point query[n_query];
   for( i = 0; i < n_query; i++ ) {
      query[i].id = i + 1;
      query[i].vec = calloc( DIM, sizeof( VEC_TYPE ) );
      query[i].ancestors = NULL;
      for( j = 0; j < DIM; j++ )
         ((VEC_TYPE*)query[i].vec)[j] = rand() % 256;
   }

#ifdef DEBUG
   // Dump the query vectors for Mathematica
   printf("query = {");
   for( i = 0; i < n_query; i++ ) {
      printf("{");
      for( j = 0; j < DIM; j++ ) {
         printf("%d", ((VEC_TYPE*)query[i].vec)[j]);
         if( j < DIM-1 )
            printf(",");
      }
      if( i < n_query-1 )
         printf("},");
      else
         printf("}");
   }
   printf("};\n");
#endif

   // Perform a range search on the query
   printf("(* range search *)\n");
   double range = 40;
   printf("range = %f;\n", range);
   ap_PointList *results[n_query];
   for( i = 0; i < n_query; i++ ) {
      results[i] = NULL;
      range_search( tree, &query[i], range, &results[i], dist );
   }

#ifdef DEBUG
   // Dump the search results for Mathematica
   printf("(* query results *)\n");
   ap_PointList *index;
   printf("results = {");
   for( i = 0; i < n_query; i++ ) {
      printf("{");
      for( index = results[i]; index != NULL; index = index->next ) {
         printf("%d", index->p->id);
         if( index->next != NULL )
            printf(",");
      }
      if( i < n_query-1 )
         printf("},");
      else
         printf("}");
   }
   printf("};\n");
#endif


   /*
   // Copy s into t
   ap_PointList *t = copy_list( s );
   ap_PointList *i0 = s, *j0 = t;

   // Move members of t into u
   printf("members of t:\n");
   i0 = t;
   while( i0 != NULL ) {
      printf("(* i0 = %p *)\n", i0);
      i0 = i0->next;
   }
   printf("moving members of t into u\n", i);
   ap_PointList *u = NULL;
   while( list_size(t) > 0 ) {
      move_point( 0, &t, &u );
   }
   printf("members of u:\n");
   i0 = u;
   while( i0 != NULL ) {
      printf("(* i0 = %p *)\n", i0);
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
   // Test for mem leaks in build_tree, make_cluster,
   // free_tree, and free_cluster
   for( i = 0; i < 1e6; i++ ) {
      free_tree( tree );
      tree = build_tree( s, 256*0.05*sqrt(DIM), NULL, NULL, DIM, dist );
   }
   */

   /*
   // Test for mem leaks in free_list and add_point
   for( i = 0; i < 1e7; i++ ) {
      free_list(s);
      s = NULL;
      for( j = 0; j < n; j++ )
         add_point( &s, &data[j], 0 );
   }
   */

   /*
   // Test for mem leaks in list_size and move_point
   for( i = 0; i < 1e7; i++ ) {
      while( list_size(s) > 0 )
         move_point( 0, &s, &t );
      while( list_size(t) > 0 )
         move_point( 0, &t, &s );
   }
   */

   return 0;
}

