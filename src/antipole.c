/* antipole.c
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

#include <assert.h>  /* assert */
#include <math.h>    /* fmax */
#include <stdio.h>   /* printf */
#include <stdlib.h>  /* NULL, rand */
#include "antipole.h"

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

// Create an ap_Tree that serves as the root, an internal
// node, or a leaf for the tree data structure. Non-leaves
// contain the identities of two antipole points, a left
// subtree and a right subtree which each contain the
// subset of points that is nearest its respective antipole
// point, and the radii of the subsets. Leaves contain a
// cluster of points.
ap_Tree*
build_tree( ap_PointList *set, double target_radius, ap_Point *antipole_a, ap_Point *antipole_b, int dimensionality, DIST_FUNC ) {

#ifdef DEBUG
   static int depth = -1;
   depth++;
#endif

   // Create the new ap_Tree
   ap_Tree *new_tree = malloc( sizeof( ap_Tree ) );
   assert( new_tree );

#ifdef DEBUG
   if( depth == 0 )
      printf("tree = {");
   printf("%ld->%d,", (long)new_tree, list_size( set ));
#endif

   // Determine if this tree is an internal node or a leaf
   if( antipole_a == NULL || antipole_b == NULL ) {
      adapted_approx_antipoles( set, &antipole_a, &antipole_b, target_radius, dist );
      if( antipole_a == NULL || antipole_b == NULL ) {
         // If it is a leaf, create a cluster from the set and return
         // the leaf
         new_tree->is_leaf = 1;
         new_tree->cluster = make_cluster( set, dimensionality, dist );
#ifdef DEBUG
         depth--;
#endif
         return new_tree;
      }
   }

   // If this tree is an internal node, initialize it
   new_tree->is_leaf = 0;
   new_tree->a = antipole_a;
   new_tree->b = antipole_b;
   new_tree->radius_a = 0;
   new_tree->radius_b = 0;

   // For each point in the set, find the distance to each
   // antipole, store the distances in the point's ancestor
   // list, add the point to the subset belonging to the
   // nearest antipole, and update the radius of the subset if
   // necessary
   double dist_a, dist_b;
   ap_PointList *index, *set_a = NULL, *set_b = NULL;
   for( index = set; index != NULL; index = index->next ) {
      dist_a = dist( new_tree->a, index->p );
      dist_b = dist( new_tree->b, index->p );
      add_point( &(index->p->ancestors), new_tree->a, dist_a );
      add_point( &(index->p->ancestors), new_tree->b, dist_b );
      if( dist_a < dist_b ) {
         add_point( &set_a, index->p, dist_a );
         new_tree->radius_a = fmax( dist_a, new_tree->radius_a );
      } else {
         add_point( &set_b, index->p, dist_b );
         new_tree->radius_b = fmax( dist_b, new_tree->radius_b );
      }
   }

   // Build subtrees as children for this node using the two
   // point subsets
   check_for_antipoles( set_a, target_radius, new_tree->a, &antipole_a, &antipole_b );
   new_tree->left = build_tree( set_a, target_radius, antipole_a, antipole_b, dimensionality, dist );
   check_for_antipoles( set_b, target_radius, new_tree->b, &antipole_a, &antipole_b );
   new_tree->right = build_tree( set_b, target_radius, antipole_a, antipole_b, dimensionality, dist );

#ifdef DEBUG
   printf("{%ld->%ld,%d},", (long)new_tree, (long)new_tree->left, new_tree->a->id);
   if( depth == 0 )
      printf("{%ld->%ld,%d}};\n", (long)new_tree, (long)new_tree->right, new_tree->b->id);
   else
      printf("{%ld->%ld,%d},", (long)new_tree, (long)new_tree->right, new_tree->b->id);
   depth--;
#endif

   free_list( set_a );
   free_list( set_b );

   return new_tree;
}


// Create an ap_Cluster owned by a leaf of the tree data
// structure containing a list of the points in the cluster
// (already determined to be sufficiently close to one
// another to group together), the identity of the geometric
// median of the cluster, and the cluster radius.
ap_Cluster*
make_cluster( ap_PointList *set, int dimensionality, DIST_FUNC ) {

   ap_PointList *index;
   double dist_centroid;

   // Create the new ap_Cluster and initialize it
   ap_Cluster *new_cluster = malloc( sizeof( ap_Cluster ) );
   assert( new_cluster );
   approx_1_median( set, &(new_cluster->centroid), dimensionality, dist );
   new_cluster->radius = 0;
   new_cluster->members = NULL;

   // For every point in the set (besides the centroid), find
   // the distance to the centroid, add the point to the list
   // of points in the cluster, and update the radius of the
   // cluster if necessary
   for( index = set; index != NULL; index = index->next ) {
      if( index->p != new_cluster->centroid ) {
         dist_centroid = dist( new_cluster->centroid, index->p );
         add_point( &(new_cluster->members), index->p, dist_centroid );
         new_cluster->radius = fmax( new_cluster->radius, dist_centroid );
      }
   }

   return new_cluster;
}


// Prepend to an ap_PointList an ap_Point with a distance
// value (to an ancestor, cluster centroid, or query,
// depending on the use of the ap_PointList). Requires the
// address of an ap_PointList pointer so that the list can
// be given a new first member. The ap_PointList pointer
// should be set to NULL before calling this function if the
// list is empty; otherwise the list may not terminate
// properly. Returns 1 if the point was added or 0 if it was
// not due to lack of uniqueness.
int
add_point( ap_PointList **set, ap_Point *p, double dist ) {

   // Check for uniqueness
   ap_PointList *index = *set;
   while( index != NULL ) {
      if( index->p == p )
         return 0;
      index = index->next;
   }

   // Add the point if it is unique to the list
   ap_PointList *new_list_member = malloc( sizeof( ap_PointList ) );
   assert( new_list_member );
   new_list_member->p = p;
   new_list_member->dist = dist;
   new_list_member->next = *set;
   *set = new_list_member;

   return 1;
}


// Move the first instance of point p found in the
// ap_PointList *from into the ap_PointList *to. Requires
// the addresses of the ap_PointList pointers so that the
// moved member can be prepended to *to, and so that if the
// first member of ap_PointList *from is moved, the second
// member can become the new first member. The ap_PointList
// pointer *to should be set to NULL before calling this
// function if the list is empty; otherwise the list may not
// terminate properly. Returns 1 if the point was moved or 0
// if it was not found
int
move_point( ap_Point *p, ap_PointList **from, ap_PointList **to ) {

   int i = 0;
   ap_PointList *before = NULL, *index = *from;

   // Find the first ap_PointList containing p, and store the
   // link preceding it for later user
   while( index != NULL && index->p != p ) {
      before = index;
      index = index->next;
      i++;
   }

   // Return if p was not found
   if( index == NULL )
      return 0;

   // Remove index from *from
   if( i == 0 )
      *from = index->next;
   else
      before->next = index->next;

   // Prepend index to *to
   index->next = *to;
   *to = index;

   return 1;
}


// Move the n-th member (start counting at zero) of the
// ap_PointList *from into the ap_PointList *to. Requires
// the addresses of the ap_PointList pointers so that the
// moved member can be prepended to *to, and so that if the
// first member of ap_PointList *from is moved, the second
// member can become the new first member. The ap_PointList
// pointer *to should be set to NULL before calling this
// function if the list is empty; otherwise the list may not
// terminate properly. Returns 1 if the point was moved or 0
// if n was too large.
int
move_nth_point( int n, ap_PointList **from, ap_PointList **to ) {

   int i;
   ap_PointList *before = NULL, *index = *from;

   // Find the n-th ap_PointList, and store the link preceding
   // it for later user
   for( i = 0; i < n && index != NULL; i++ ) {
      before = index;
      index = index->next;
   }

   // Return if n was too large
   if( index == NULL )
      return 0;

   // Remove index from *from
   if( n == 0 )
      *from = index->next;
   else
      before->next = index->next;

   // Prepend index to *to
   index->next = *to;
   *to = index;

   return 1;
}


// Create a new ap_PointList with the same contents as set.
ap_PointList*
copy_list( ap_PointList *set ) {

   ap_PointList *index = set, *new_list = NULL;
   while( index != NULL ) {
      add_point( &new_list, index->p, index->dist );
      index = index->next;
   }
   return new_list;
}


// Find the size of a set of points
int
list_size( ap_PointList *set ) {

   int size = 0;
   while( set != NULL ) {
      size++;
      set = set->next;
   }
   return size;
}


// Create a new ap_Heap object with a minimum capacity
ap_Heap*
create_heap( int capacity ) {

   ap_Heap *heap = malloc( sizeof( ap_Heap ) );
   assert( heap );
   heap->capacity = max( capacity, 10 );
   heap->size = 0;
   heap->items = calloc( heap->capacity, sizeof( void* ) );
   heap->dists = calloc( heap->capacity, sizeof( double ) );
   assert( heap->items && heap->dists );
   heap->max_index = -1;
   heap->max_dist  = -1;

   return heap;
}


// Increase the capacity of heap
void
heap_grow( ap_Heap *heap ) {

   int i;
   heap->capacity *= 2;
   void **temp_items  = calloc( heap->capacity, sizeof( void* ) );
   double *temp_dists = calloc( heap->capacity, sizeof( double ) );
   assert( temp_items && temp_dists );

   for( i = 0; i < heap->size; i++ ) {
      temp_items[i] = heap->items[i];
      temp_dists[i] = heap->dists[i];
   }
   free( heap->items );
   free( heap->dists );
   heap->items = temp_items;
   heap->dists = temp_dists;
}


// Swap items at positions i and j in heap. Returns 1 if a
// swap occurred or 0 if the indices were inapproriate
int
heap_swap( ap_Heap *heap, int i, int j ) {

   // Return if the indices match or are out of bounds
   if( i == j || i >= heap->size || j >= heap->size || i < 0 || j < 0 )
      return 0;

   void *temp_item  = heap->items[i];
   double temp_dist = heap->dists[i];
   heap->items[i] = heap->items[j];
   heap->dists[i] = heap->dists[j];
   heap->items[j] = temp_item;
   heap->dists[j] = temp_dist;

   return 1;
}


// Sift upward the item at position index in heap until it
// is properly positioned (the heap is a min-heap).
void
heap_sift_up( ap_Heap* heap, int index ) {

   // Return is index is the root
   if( index == 0 )
      return;

   int parent = ( index - 1 ) / 2;

   // If index has a smaller dist than its parent, swap them
   // and continue sifting upward
   if( heap->dists[index] < heap->dists[parent] ) {
      heap_swap( heap, index, parent );
      heap_sift_up( heap, parent );
   }
}


// Sift downward the item at position index in heap until it
// is properly positioned (the heap is a min-heap).
void
heap_sift_down( ap_Heap *heap, int index ) {

   // Return if index is a leaf
   if( index >= heap->size / 2 )
      return;

   int smallestChild;
   int left = 2 * index + 1;
   int right = 2 * index + 2;

   // Find the smallest child of index
   if( right >= heap->size ) {
      smallestChild = left;
   } else {
      if( heap->dists[left] < heap->dists[right] )
         smallestChild = left;
      else
         smallestChild = right;
   }

   // If index has a greater dist than one of its children,
   // swap them and continue sifting downward
   if( heap->dists[index] > heap->dists[smallestChild] ) {
      heap_swap( heap, index, smallestChild );
      heap_sift_down( heap, smallestChild );
   }
}


// Attempt to add an item to an ap_Heap with a distance
// value used for sorting (the heap is a min-heap). The item
// will be inserted only if the heap has fewer items than
// the specified limit, if a limit is not specified (use
// limit < 1). The item will replace the max item in the
// heap if the limit is reached and dist is smaller than
// max_dist. Returns 1 if the item was inserted or 0 if the
// heap was full or dist was too large.
int
heap_try_insert( ap_Heap *heap, int limit, void *item, double dist ) {

   int i;

   // Add the item to the heap if a limit is not given or the
   // current heap size is smaller than the limit
   if( limit < 1 || heap->size < limit ) {
      
      // Return if the item is already in the list
      for( i = 0; i < heap->size; i++ )
         if( item == heap->items[i] )
            return 0;

      // Grow the data arrays if necessary
      if( heap->size == heap->capacity )
         heap_grow( heap );

      // Insert the new item at the end of the array
      i = heap->size;
      heap->items[i] = item;
      heap->dists[i] = dist;
      heap->size++;

      // Update max_index and max_dist if the item is the new
      // max item or sift it upward
      if( dist > heap->max_dist ) {
         heap->max_index = i;
         heap->max_dist = dist;
      } else {
         heap_sift_up( heap, i );
      }

      return 1;
   } else {

      // If the current heap size is not smaller than the limit,
      // add the item only if dist is smaller than max_dist
      if( dist < heap->max_dist ) {

         // Return if the item is already in the list
         for( i = 0; i < heap->size; i++ )
            if( item == heap->items[i] )
               return 0;

         // Replace the max item with the new item
         i = heap->max_index;
         heap->items[i] = item;
         heap->dists[i] = dist;

         // Sift the item upward
         heap_sift_up( heap, i );

         // Update max_index and max_dist
         heap->max_index = -1;
         heap->max_dist  = -1;
         for( i = heap->size / 2; i < heap->size; i++ ) {
            if( heap->dists[i] > heap->max_dist ) {
               heap->max_index = i;
               heap->max_dist = heap->dists[i];
            }
         }

         return 1;
      } else {

         // Return if dist is not smaller than max_dist
         return 0;
      }
   }
}


// Remove the first item in the heap and return a pointer to
// it. Returns NULL if the heap is empty.
void*
heap_pop( ap_Heap *heap ) {

   // Return if the heap is empty
   if( heap->size == 0 )
      return NULL;

   void* min_item = heap->items[0];

   // Replace the root with the last item in the list
   heap->items[0] = heap->items[heap->size-1];
   heap->dists[0] = heap->dists[heap->size-1];
   heap->size--;

   // Sift the new root downward
   heap_sift_down( heap, 0 );

   return min_item;
}


// Create an ap_PointList from an ap_Heap. The points in the
// list will be sorted by dist in ascending order. This
// function assumes the items in the heap are ap_Points.
ap_PointList*
heap_to_list( ap_Heap *heap ) {

   int i, counter = 0;
   ap_PointList *new_list = NULL;

   // Create a copy of the heap (because it will be modified
   // below)
   ap_Heap *new_heap = malloc( sizeof( ap_Heap ) );
   assert( new_heap );
   new_heap->capacity = heap->capacity;
   new_heap->size = heap->size;
   new_heap->items = calloc( new_heap->capacity, sizeof( void* ) );
   new_heap->dists = calloc( new_heap->capacity, sizeof( double ) );
   assert( new_heap->items && new_heap->dists );
   for( i = 0; i < heap->size; i++ ) {
      new_heap->items[i] = heap->items[i];
      new_heap->dists[i] = heap->dists[i];
   }
   new_heap->max_index = heap->max_index;
   new_heap->max_dist = heap->max_dist;

   // Sort the items in heap by dist in descending order
   for( i = new_heap->size - 1; i > 0; i-- ) {
      heap_swap( new_heap, 0, i );
      new_heap->size--;
      counter++;
      heap_sift_down( new_heap, 0 );
   }
   new_heap->size += counter;

   // Copy the points from the heap into a list
   for( i = 0; i < new_heap->size; i++ ) {
      add_point( &new_list, (ap_Point*)new_heap->items[i], new_heap->dists[i] );
   }
   free_heap( new_heap );

   return new_list;
}


// Recursively free up memory used by an ap_Tree.
void
free_tree( ap_Tree *tree ) {

   if( tree != NULL ) {
      if( tree->is_leaf ) {
         free_cluster( tree->cluster );
      } else {
         free_tree( tree->left );
         free_tree( tree->right );
      }
      free( tree );
   }
}


// Free up memory used by an ap_Cluster.
void
free_cluster( ap_Cluster *cluster ) {

   if( cluster != NULL ) {
      free_list( cluster->members );
      free( cluster );
   }
}


// Recursively free up memory used by an ap_PointList.
void
free_list( ap_PointList *set ) {

   if( set != NULL ) {
      free_list( set->next );
      free( set );
   }
}


// Free up memory used by an ap_Heap.
void free_heap( ap_Heap *heap ) {

   if( heap != NULL ) {
      free( heap->items );
      free( heap->dists );
      free( heap );
   }
}


// Find the exact geometric median of a set of points and
// store it in median.
void
exact_1_median( ap_PointList *set, ap_Point **median, DIST_FUNC ) {

   *median = NULL;

   int i, j, d, size = list_size( set );
   ap_PointList *i_list, *j_list;

   // Initialize the array of distance sums
   double sums[size];
   for( i = 0; i < size; i++ )
      sums[i] = 0;

   // Calculate the distance between each pair of points and
   // add each distance to the distance sums for each point in
   // the pair
   for( i = 0, i_list = set; i < size; i++, i_list = i_list->next ) {
      for( j = i + 1, j_list = i_list->next; j < size; j++, j_list = j_list->next ) {
         d = dist( i_list->p, j_list->p );
         sums[i] += d;
         sums[j] += d;
      }
   }

   // Identify the point with the minimum distance sum and
   // store it in median
   double min_sum = -1;
   for( i = 0, i_list = set; i < size; i++, i_list = i_list->next ) {
      if( sums[i] < min_sum || min_sum < 0 ) {
         min_sum = sums[i];
         *median = i_list->p;
      }
   }
}


// Find an approximation for the geometric median of a set
// of points and store it in median. The user should
// initialize the random number generator using srand.
void
approx_1_median( ap_PointList *set, ap_Point **median, int dimensionality, DIST_FUNC ) {

   *median = NULL;

   ap_PointList *contestants = copy_list( set ), *tournament, *winners;
   int i, contestants_size = list_size( contestants ), tournament_size = dimensionality + 1, winners_size;
   int final_round_size = max( pow( tournament_size, 2 ) - 1, round( sqrt( list_size( set ) ) ) );

   // Hold a series of rounds of tournaments
   while( contestants_size > final_round_size ) {
      // Find the winners that will continue to the next round
      winners = NULL;
      winners_size = 0;
      while( contestants_size >= 2 * tournament_size ) {
         tournament = NULL;
         // Move tournament_size random members of contestants into
         // tournament
         for( i = 0; i < tournament_size; i++ ) {
            move_nth_point( rand() % contestants_size, &contestants, &tournament );
            contestants_size--;
         }
         // Find the winner of this tournament and discard the losers
         exact_1_median( tournament, median, dist );
         move_point( *median, &tournament, &winners );
         winners_size++;
         free_list( tournament );
      }
      // Find the winner among the remaining contestants and
      // discard the losers
      exact_1_median( contestants, median, dist );
      move_point( *median, &contestants, &winners );
      winners_size++;
      free_list( contestants );

      // Fill the pool of contestants with all the winners in
      // preparation for the next round
      contestants = winners;
      contestants_size = winners_size;
   }
   
   // Find the overall winner and discard the losers
   exact_1_median( contestants, median, dist );
   free_list( contestants );
}


// Find the two points in the set that are farthest from one
// another and store them in antipole_a and antipole_b.
void
exact_antipoles( ap_PointList *set, ap_Point **antipole_a, ap_Point **antipole_b, DIST_FUNC ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   ap_PointList *i_list, *j_list;
   double d, max_dist = -1;

   // Calculate the distance between each pair of points and
   // if a distance is greater than any found yet, make the
   // pair of points the new antipole pair
   for( i_list = set; i_list != NULL; i_list = i_list->next ) {
      for( j_list = i_list->next; j_list != NULL; j_list = j_list->next ) {
         d = dist( i_list->p, j_list->p );
         if( d > max_dist ) {
            *antipole_a = i_list->p;
            *antipole_b = j_list->p;
            max_dist = d;
         }
      }
   }
}


// Find an approximation for the antipole pair of a set
// of points and store them in antipole_a and antipole_b.
// The user should initialize the random number generator
// using srand.
void
approx_antipoles( ap_PointList *set, ap_Point **antipole_a, ap_Point **antipole_b, int dimensionality, DIST_FUNC ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   ap_PointList *contestants = copy_list( set ), *tournament, *winners;
   int i, contestants_size = list_size( contestants ), tournament_size = dimensionality + 1, winners_size;
   int final_round_size = max( pow( tournament_size, 2 ) - 1, round( sqrt( list_size( set ) ) ) );

   // Hold a series of rounds of tournaments
   while( contestants_size > final_round_size ) {
      // Find the winners that will continue to the next round
      winners = NULL;
      winners_size = 0;
      while( contestants_size >= 2 * tournament_size ) {
         tournament = NULL;
         // Move tournament_size random members of contestants into
         // tournament
         for( i = 0; i < tournament_size; i++ ) {
            move_nth_point( rand() % contestants_size, &contestants, &tournament );
            contestants_size--;
         }
         // Find the winners of this tournament and discard the losers
         exact_antipoles( tournament, antipole_a, antipole_b, dist );
         move_point( *antipole_a, &tournament, &winners );
         move_point( *antipole_b, &tournament, &winners );
         winners_size += 2;
         free_list( tournament );
      }
      // Find the winners among the remaining contestants and
      // discard the losers
      exact_antipoles( contestants, antipole_a, antipole_b, dist );
      move_point( *antipole_a, &contestants, &winners );
      move_point( *antipole_b, &contestants, &winners );
      winners_size += 2;
      free_list( contestants );

      // Fill the pool of contestants with all the winners in
      // preparation for the next round
      contestants = winners;
      contestants_size = winners_size;
   }
   
   // Find the overall winners and discard the losers
   exact_antipoles( contestants, antipole_a, antipole_b, dist );
   free_list( contestants );
}


// Search for any two points whose distance from
// one another is greater than the target cluster
// diameter
void
adapted_approx_antipoles( ap_PointList *set, ap_Point **antipole_a, ap_Point **antipole_b, double target_radius, DIST_FUNC ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   ap_PointList *i_list, *j_list;

   // Calculate the distance between each pair of points and
   // if a distance is greater than the target cluster diameter
   // make the pair of points the new antipole pair and stop
   // searching
   for( i_list = set; i_list != NULL; i_list = i_list->next ) {
      for( j_list = i_list->next; j_list != NULL; j_list = j_list->next ) {
         if( dist( i_list->p, j_list->p ) > 2 * target_radius ) {
            *antipole_a = i_list->p;
            *antipole_b = j_list->p;
            return;
         }
      }
   }
}


// Search the set for a point that when paired with ancestor
// could serve as an antipole pair.
void
check_for_antipoles( ap_PointList *set, double target_radius, ap_Point *ancestor, ap_Point **antipole_a, ap_Point **antipole_b ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   // Search the set for a point whose distance to the ancestor
   // is greater than 2*target_radius and save it and the
   // ancestor as antipoles
   ap_PointList *index_i, *index_j;
   for( index_i = set; index_i != NULL; index_i = index_i->next ) {
      for( index_j = index_i->p->ancestors; index_j != NULL; index_j = index_j->next ) {
         if( index_j->p == ancestor ) {
            if( index_j->dist > 2 * target_radius ) {
               *antipole_a = ancestor;
               *antipole_b = index_i->p;
               return;
            }
         }
      }
   }
}


// Search the tree recursively to find all points within
// range of query and place them in out.
void
range_search( ap_Tree *tree, ap_Point *query, double range, ap_PointList **out, DIST_FUNC ) {

   if( !tree->is_leaf ) {
      // Calculate the distance between query and the antipoles
      // and store these values in the ancestor list for query
      double dist_a = dist( tree->a, query );
      double dist_b = dist( tree->b, query );
      /*
      int a_added = add_point( &(query->ancestors), tree->a, dist_a );
      int b_added = add_point( &(query->ancestors), tree->b, dist_b );
      */

      // If either antipole is within range, add it to out
      if( dist_a <= range )
         add_point( out, tree->a, dist_a );
      if( dist_b <= range )
         add_point( out, tree->b, dist_b );

      // Use the triangle inequality to determine if each subtree
      // is within range of the query, and descend those subtrees
      // that are
      if( dist_a <= range + tree->radius_a )
         range_search( tree->left, query, range, out, dist );
      if( dist_b <= range + tree->radius_b )
         range_search( tree->right, query, range, out, dist );

      /*
      // Once the search has returned from both subtrees, remove
      // the antipoles of tree from the ancestor list for query if
      // they were not ancestors in subtrees above this one, since
      // they will no longer be relevant to further searches
      ap_PointList *trash = NULL;
      if( b_added )
         move_nth_point( 0, &(query->ancestors), &trash );
      if( a_added )
         move_nth_point( 0, &(query->ancestors), &trash );
      free_list( trash );
      */
   } else {
      // If tree is a leaf, search its cluster for points within
      // range of query
      range_search_cluster( tree->cluster, query, range, out, dist );
   }
}


// Find all members of the cluster that are within range of
// query and place them in out.
void
range_search_cluster( ap_Cluster *cluster, ap_Point *query, double range, ap_PointList **out, DIST_FUNC ) {

   // Calculate the distance between the query and the centroid
   // and add it to out if it is within range
   double d, dist_centroid = dist( cluster->centroid, query );
   if( dist_centroid <= range )
      add_point( out, cluster->centroid, dist_centroid );

   // Use the triangle inequality with the cluster radius to
   // determine if the entire cluster can be excluded as a
   // group
   if( dist_centroid > range + cluster->radius )
      return;

   // Use the triangle inequality with the cluster radius to
   // determine if the entire cluster can be included as a
   // group
   if( dist_centroid <= range - cluster->radius ) {
      ap_PointList *index = cluster->members;
      while( index != NULL ) {
         add_point( out, index->p, -1 );
         index = index->next;
      }
      return;
   }

   // Check each member of the cluster
   ap_PointList *index = cluster->members;
   /*
   ap_PointList *query_ancestors, *cluster_ancestors;
   */
   while( index != NULL ) {
      // Use the triangle inequality with the cluster member's
      // distance to centroid to determine if the point is
      // definitely out of range
      if( dist_centroid > range + index->dist )
         goto next_cluster_member;

      // Use the triangle inequality with the cluster member's
      // distance to centroid to determine if the point is
      // definitely within range
      if( dist_centroid <= range - index->dist ) {
         add_point( out, index->p, -1 );
         goto next_cluster_member;
      }

      /*
      // Check the ancestors of the query and the member of the
      // cluster
      query_ancestors = query->ancestors;
      cluster_ancestors = index->p->ancestors;
      while( query_ancestors != NULL ) {
         assert( query_ancestors->p == cluster_ancestors->p );

         // Use the triangle inequality with the cluster member's
         // distance to ancestor to determine if the point is
         // definitely out of range
         if( query_ancestors->dist > range + cluster_ancestors->dist )
            goto next_cluster_member;

         // Use the triangle inequality with the cluster member's
         // distance to ancestor to determine if the point is
         // definitely within range
         if( query_ancestors->dist <= range - cluster_ancestors->dist ) {
            add_point( out, index->p, -1 );
            goto next_cluster_member;
         }

         query_ancestors = query_ancestors->next;
         cluster_ancestors = cluster_ancestors->next;
      }
      */

      // Finally, if all methods of using precalculated distances
      // to rule-out or rule-in the cluster member have failed,
      // calculate the distance between the query and the cluster
      // member and add it to out if it is within range
      d = dist( index->p, query );
      if( d <= range )
         add_point( out, index->p, d );

next_cluster_member:
      index = index->next;
   }
}


// Search the tree using a priority queue for the subtrees
// to find the k points nearest the query and place them in
// out.
void
nearest_neighbor_search( ap_Tree *tree, ap_Point *query, int k, ap_PointList **out, DIST_FUNC ) {

   double dist_a, dist_b;
   ap_Tree *index;
   ap_Heap *tree_pq = create_heap( 0 ), *point_pq = create_heap( k );

   // Initialize the tree priority queue with the root of the
   // tree
   heap_try_insert( tree_pq, 0, tree, -1 );

   // Search through the subtrees in order of proximity to the
   // query until there are no more subtrees to search or the
   // remaining subtrees are all farther away to the query than
   // the k points already found
   while( tree_pq->size > 0 ) {

      // If the next nearest subtree is farther away than the
      // farthest member of point_pq, then quit
      if( point_pq->size == k && tree_pq->dists[0] > point_pq->max_dist )
         break;

      // Pop off the next subtree in the tree priority queue
      index = (ap_Tree*)heap_pop( tree_pq );

      if( !index->is_leaf ) {
         // Calculate the distance between query and the antipoles,
         // store these values in the ancestor list for query, and
         // add each antipole to the point priority queue if it is
         // nearer than than the queue's farthest member
         dist_a = dist( index->a, query );
         dist_b = dist( index->b, query );
         heap_try_insert( point_pq, k, index->a, dist_a );
         heap_try_insert( point_pq, k, index->b, dist_b );
         /*
         add_point( &(query->ancestors), index->a, dist_a );
         add_point( &(query->ancestors), index->b, dist_b );
         */

         // Add the subtree's children to the tree priority queue
         heap_try_insert( tree_pq, 0, index->left,  dist_a - index->radius_a );
         heap_try_insert( tree_pq, 0, index->right, dist_b - index->radius_b );
      } else {
         // If tree is a leaf, search its cluster for points that
         // should be added to the point priority queue
         nearest_neighbor_search_cluster( index->cluster, query, k, point_pq, dist );
      }
   }

   // Convert the point priority queue into an ap_PointList
   // and store it in out
   *out = heap_to_list( point_pq );

   // Free up the memory used by the tree and point priority
   // queues
   free_heap( tree_pq );
   free_heap( point_pq );
}


// Find any members of the cluster that are nearer to the
// query than any of the k points already found in the point
// priority queue and place them in point_pq.
void
nearest_neighbor_search_cluster( ap_Cluster *cluster, ap_Point *query, int k, ap_Heap *point_pq, DIST_FUNC ) {

   // Calculate the distance between the query and the centroid
   // and add it to point_pq if it is nearer than the queue's
   // farthest member
   double d, dist_centroid = dist( cluster->centroid, query );
   heap_try_insert( point_pq, k, cluster->centroid, dist_centroid );

   // Use the triangle inequality with the cluster radius to
   // determine if the entire cluster can be excluded as a
   // group
   if( point_pq->size == k && dist_centroid > point_pq->max_dist + cluster->radius )
      return;

   // Check each member of the cluster
   ap_PointList *index = cluster->members;
   /*
   ap_PointList *query_ancestors, *cluster_ancestors;
   */
   while( index != NULL ) {
      // Use the triangle inequality with the cluster member's
      // distance to centroid to determine if the point is
      // definitely farther away than the farthest member of
      // point_pq
      if( point_pq->size == k && dist_centroid > point_pq->max_dist + index->dist )
         goto next_cluster_member;

      // Use the triangle inequality with the cluster member's
      // distance to centroid to determine if the point is
      // definitely nearer than the farthest member of point_pq
      if( dist_centroid <= point_pq->max_dist - index->dist ) {
         d = dist( index->p, query );
         heap_try_insert( point_pq, k, index->p, d );
         goto next_cluster_member;
      }

      /*
      // Check the ancestors of the query and the member of the
      // cluster
      for( query_ancestors = query->ancestors; query_ancestors != NULL; query_ancestors = query_ancestors->next ) {
         for( cluster_ancestors = index->p->ancestors; cluster_ancestors != NULL; cluster_ancestors = cluster_ancestors->next ) {
            if( query_ancestors->p == cluster_ancestors->p ) {

               // Use the triangle inequality with the cluster member's
               // distance to ancestor to determine if the point is
               // definitely farther away than the farthest member of
               // point_pq
               if( point_pq->size == k && query_ancestors->dist > point_pq->max_dist + cluster_ancestors->dist )
                  goto next_cluster_member;

               // Use the triangle inequality with the cluster member's
               // distance to ancestor to determine if the point is
               // definitely nearer than the farthest member of point_pq
               if( query_ancestors->dist <= point_pq->max_dist - cluster_ancestors->dist ) {
                  d = dist( index->p, query );
                  heap_try_insert( point_pq, k, index->p, d );
                  goto next_cluster_member;
               }
            }
         }
      }
      */

      // Finally, if all methods of using precalculated distances
      // to rule-out or rule-in the cluster member have failed,
      // calculate the distance between the query and the cluster
      // member and add it to point_pq if it is nearer than the
      // queue's farthest member
      d = dist( index->p, query );
      heap_try_insert( point_pq, k, index->p, d );

next_cluster_member:
      index = index->next;
   }
}


