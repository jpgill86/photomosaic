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
#define PCS "%d"           /* the data type printf conversion specifier */
#define VEC_DOMAIN 256     /* the range over which the data can fall */
#define RAND_DATA rand()%VEC_DOMAIN  /* macro for generating random numbers of the right type */
//typedef double VEC_TYPE;
//#define PCS "%f"
//#define VEC_DOMAIN 1.0
//#define RAND_DATA VEC_DOMAIN*(double)rand()/(double)RAND_MAX

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

   int i, j;
   int n_data = 20, n_query = 10, n_neighbor = 5;
   double bounded_radius = VEC_DOMAIN * 0.05 * sqrt(DIM);
   double range = VEC_DOMAIN * 0.1;
   int seed = time(NULL);
   srand(seed);

   ap_Point *data[n_data];
   ap_Point *query[n_query];
   ap_PointList *results[n_query];
   ap_PointList *s;
   ap_Tree *tree;

   printf("(* parameters *)\n");
   printf("dim = %d;\n", DIM);
   printf("nData = %d;\n", n_data);
   printf("nQuery = %d;\n", n_query);
   printf("nNeighbor = %d;\n", n_neighbor);
   printf("domain = %f;\n", (double)VEC_DOMAIN);
   printf("range = %f;\n", range);
   printf("seed = %d;\n", seed);

   // Create a random data array
   printf("(* creating data points... ");
   for( i = 0; i < n_data; i++ ) {
      data[i] = malloc( sizeof( ap_Point ) );
      data[i]->id = i + 1;
      data[i]->vec = calloc( DIM, sizeof( VEC_TYPE ) );
      data[i]->ancestors = NULL;
      for( j = 0; j < DIM; j++ )
         ((VEC_TYPE*)data[i]->vec)[j] = RAND_DATA;
   }
   printf("done *)\n");

#ifdef DEBUG
   // Dump the data vectors for Mathematica
   printf("data = {");
   for( i = 0; i < n_data; i++ ) {
      printf("{");
      for( j = 0; j < DIM; j++ ) {
         printf(PCS, ((VEC_TYPE*)data[i]->vec)[j]);
         if( j < DIM-1 )
            printf(",");
      }
      if( i < n_data-1 )
         printf("},");
      else
         printf("}");
   }
   printf("};\n");
#endif

   // Place the ap_Points in an ap_PointList
   s = NULL;
   printf("(* creating data set from points... ");
   for( i = 0; i < n_data; i++ )
      add_point( &s, data[i], 0 );
   printf("done *)\n");

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
   printf("(* building tree... *)\n");
   tree = build_tree( s, bounded_radius, NULL, NULL, DIM, dist );
   printf("(* ... done *)\n");

   // Construct a set of query points
   printf("(* creating query points... ");
   for( i = 0; i < n_query; i++ ) {
      query[i] = malloc( sizeof( ap_Point ) );
      query[i]->id = i + 1;
      query[i]->vec = calloc( DIM, sizeof( VEC_TYPE ) );
      query[i]->ancestors = NULL;
      for( j = 0; j < DIM; j++ )
         ((VEC_TYPE*)query[i]->vec)[j] = RAND_DATA;
   }
   printf("done *)\n");

#ifdef DEBUG
   // Dump the query vectors for Mathematica
   printf("query = {");
   for( i = 0; i < n_query; i++ ) {
      printf("{");
      for( j = 0; j < DIM; j++ ) {
         printf(PCS, ((VEC_TYPE*)query[i]->vec)[j]);
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
   printf("(* performing range search... ");
   for( i = 0; i < n_query; i++ ) {
      results[i] = NULL;
      range_search( tree, query[i], range, &results[i], dist );
   }
   printf("done *)\n");

#ifdef DEBUG
   // Dump the range search results for Mathematica
   printf("rangeResults = {");
   ap_PointList *index;
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

   // Perform a nearest neighbor search on the query
   printf("(* performing nearest neighbor search... ");
   for( i = 0; i < n_query; i++ ) {
      free_list( results[i] );
      results[i] = NULL;
      nearest_neighbor_search( tree, query[i], n_neighbor, &results[i], dist );
   }
   printf("done *)\n");

#ifdef DEBUG
   // Dump the nearest neighbor search results for Mathematica
   printf("nearestNeighborResults = {");
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
   // Check for sane heap behavior
   ap_Heap *heap = create_heap( true, -1 );
   for( i = 0; i < n_data; i++ )
      heap_insert( heap, data[i], dist( query[0], data[i] ) );
   printf("\n");
   for( i = 0; i < heap->size; i++ )
      printf("(* h id=%d\tdist=%f *)\n", ((ap_Point*)heap->items[i])->id, heap->dists[i]);
   printf("\n");
   while( heap->size > n_neighbor ) {
      printf("(* r id=%d\tdist=%f *)\n", ((ap_Point*)heap->items[0])->id, heap->dists[0]);
      heap_pop( heap );
   }
   printf("\n");
   for( i = 0; i < heap->size; i++ )
      printf("(* h id=%d\tdist=%f *)\n", ((ap_Point*)heap->items[i])->id, heap->dists[i]);
   printf("\n");
   ap_PointList *v = heap_to_list( heap );
   while( v != NULL ) {
      printf("(* l id=%d\tdist=%f *)\n", v->p->id, v->dist);
      v = v->next;
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
      approx_1_median( s, &median, DIM, dist );
   */

   /*
   // Test for mem leaks in exact_antipoles
   for( i = 0; i < 1e7; i++ )
      exact_antipoles( s, &antipole_a, &antipole_b, dist );
   */

   /*
   // Test for mem leaks in approx_antipoles
   for( i = 0; i < 1e7; i++ )
      approx_antipoles( s, &antipole_a, &antipole_b, DIM, dist );
   */

   /*
   // Test for mem leaks in adapted_approx_antipoles
   for( i = 0; i < 2e8; i++ )
      adapted_approx_antipoles( s, &antipole_a, &antipole_b, bounded_radius, dist );
   */

   /*
   // Test for mem leaks in build_tree, make_cluster,
   // free_tree, and free_cluster
   for( i = 0; i < 1e6; i++ ) {
      free_tree( tree );
      tree = build_tree( s, bounded_radius, NULL, NULL, DIM, dist );
   }
   */

   /*
   // Test for mem leaks in free_list and add_point
   for( i = 0; i < 1e7; i++ ) {
      free_list(s);
      s = NULL;
      for( j = 0; j < n_data; j++ )
         add_point( &s, data[j], 0 );
   }
   */

   /*
   // Test for mem leaks in list_size, move_point, and
   // move_nth_point
   ap_PointList *t = NULL;
   for( i = 0; i < 1e7; i++ ) {
      while( s != NULL )
         move_nth_point( 0, &s, &t );
      while( list_size(t) > 0 )
         move_point( t->p, &t, &s );
   }
   */

   /*
   // Test for mem leaks in create_heap, heap_insert,
   // heap_remove, and free_heap
   ap_Heap *heap = NULL;
   for( i = 0; i < 2e7; i++ ) {
      free_heap( heap );
      heap = create_heap( false, -1 );
      for( j = 0; j < n_data; j++ )
         heap_insert( heap, data[j], dist( query[0], data[j] ) );
      while( heap->size > n_neighbor )
         heap_pop( heap );
   }
   */

   /*
   // Test for mem leaks in heap_to_list
   ap_Heap *heap = create_heap( false, -1 );
   for( i = 0; i < n_data; i++ )
      heap_insert( heap, data[i], dist( query[0], data[i] ) );
   for( i = 0; i < 2e6; i++ ) {
      free_list( heap_to_list( heap ) );
   }
   */

   /*
   // Test for mem leaks in range_search
   for( i = 0; i < 1e6; i++ ) {
      for( j = 0; j < n_query; j++ ) {
         free_list( results[j] );
         results[j] = NULL;
         range_search( tree, query[j], range, &results[j], dist );
      }
   }
   */

   /*
   // Test for mem leaks in nearest_neighbor_search
   for( i = 0; i < 3e5; i++ ) {
      for( j = 0; j < n_query; j++ ) {
         free_list( results[j] );
         results[j] = NULL;
         nearest_neighbor_search( tree, query[j], n_neighbor, &results[j], dist );
      }
   }
   */

   /*
   // Naive range search
   printf("(* performing naive range search... ");
   int k;
   double d;
   for( i = 0; i < 1e6; i++ ) {
      for( j = 0; j < n_query; j++ ) {
         free_list( results[j] );
         results[j] = NULL;
         for( k = 0; k < n_data; k++ ) {
            d = dist( query[j], data[k] );
            if( d <= range )
               add_point( &results[j], data[k], d );
         }
      }
   }
   printf("done *)\n");

#ifdef DEBUG
   // Dump the naive range search results for Mathematica
   printf("naiveRangeResults = {");
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
   */

   /*
   // Naive nearest neighbor search
   printf("(* performing naive nearest neighbor search... ");
   int k;
   ap_Heap *point_pq = NULL;
   for( i = 0; i < 1e6; i++ ) {
      for( j = 0; j < n_query; j++ ) {
         free_list( results[j] );
         results[j] = NULL;
         free_heap( point_pq );
         point_pq = create_heap( true, n_neighbor );;
         for( k = 0; k < n_data; k++ )
            nearest_neighbor_search_try_point( point_pq, data[k], dist( query[j], data[k] ) );
         results[j] = heap_to_list( point_pq );
      }
   }
   printf("done *)\n");

#ifdef DEBUG
   // Dump the naive nearest neighbor search results for Mathematica
   printf("naiveNearestNeighborResults = {");
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
   */

   printf("(* ----------------------- *)\n");
   return 0;
}

