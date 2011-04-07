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


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
               TREE CONSTRUCTION FUNCTIONS
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


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
      first_approx_antipoles( set, &antipole_a, &antipole_b, target_radius, dist );
      if( antipole_a == NULL || antipole_b == NULL ) {
         // If it is a leaf, create a cluster from the set and return
         // the leaf
         new_tree->is_leaf = true;
         new_tree->cluster = build_cluster( set, dimensionality, dist );
#ifdef DEBUG
         depth--;
#endif
         return new_tree;
      }
   }

   // If this tree is an internal node, initialize it
   new_tree->is_leaf = false;
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
   ap_PointList *set_a = NULL, *set_b = NULL;
   while( set != NULL ) {
      dist_a = dist( new_tree->a, set->p );
      dist_b = dist( new_tree->b, set->p );
      add_point( &(set->p->ancestors), new_tree->a, dist_a );
      add_point( &(set->p->ancestors), new_tree->b, dist_b );
      if( dist_a < dist_b ) {
         add_point( &set_a, set->p, dist_a );
         new_tree->radius_a = fmax( dist_a, new_tree->radius_a );
      } else {
         add_point( &set_b, set->p, dist_b );
         new_tree->radius_b = fmax( dist_b, new_tree->radius_b );
      }
      set = set->next;
   }

   // Build subtrees as children for this node using the two
   // point subsets
   check_ancestors_for_antipoles( set_a, target_radius, new_tree->a, &antipole_a, &antipole_b );
   new_tree->left = build_tree( set_a, target_radius, antipole_a, antipole_b, dimensionality, dist );
   check_ancestors_for_antipoles( set_b, target_radius, new_tree->b, &antipole_a, &antipole_b );
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
build_cluster( ap_PointList *set, int dimensionality, DIST_FUNC ) {

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
   while( set != NULL ) {
      if( set->p != new_cluster->centroid ) {
         dist_centroid = dist( new_cluster->centroid, set->p );
         add_point( &(new_cluster->members), set->p, dist_centroid );
         new_cluster->radius = fmax( new_cluster->radius, dist_centroid );
      }
      set = set->next;
   }

   return new_cluster;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                     SEARCH FUNCTIONS
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


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

   // Create the tree priority queue as a min-heap with no
   // maximum size
   ap_Heap *tree_pq = create_heap( false, 0 );

   // Create the point priority queue as a max-heap with a
   // maximum size k
   ap_Heap *point_pq = create_heap( true, k );

   // Initialize the tree priority queue with the root of the
   // tree
   heap_insert( tree_pq, tree, -1 );

   // Search through the subtrees in order of proximity to the
   // query until there are no more subtrees to search or the
   // remaining subtrees are all farther away to the query than
   // the k points already found
   while( tree_pq->size > 0 ) {

      // If point_pq already has k points and the next nearest
      // subtree is not nearer than the farthest member of
      // point_pq, then stop searching
      if( heap_is_full( point_pq ) && tree_pq->dists[0] >= point_pq->dists[0] )
         break;

      // Get the next subtree in the tree priority queue
      index = (ap_Tree*)heap_pop( tree_pq );

      if( !index->is_leaf ) {
         // Calculate the distance between query and the antipoles
         // and store these values in the ancestor list for query
         dist_a = dist( index->a, query );
         dist_b = dist( index->b, query );
         /*
         add_point( &(query->ancestors), index->a, dist_a );
         add_point( &(query->ancestors), index->b, dist_b );
         */

         // If either antipole is nearer to the query than the point
         // priority queue's farthest member, add it to point_pq
         nearest_neighbor_search_try_point( point_pq, index->a, dist_a );
         nearest_neighbor_search_try_point( point_pq, index->b, dist_b );

         // Add the subtree's children to the tree priority queue
         heap_insert( tree_pq, index->left,  dist_a - index->radius_a );
         heap_insert( tree_pq, index->right, dist_b - index->radius_b );
      } else {

         // If tree is a leaf, search its cluster for points that
         // should be added to the point priority queue
         nearest_neighbor_search_cluster( index->cluster, query, point_pq, dist );
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
nearest_neighbor_search_cluster( ap_Cluster *cluster, ap_Point *query, ap_Heap *point_pq, DIST_FUNC ) {

   // Calculate the distance between the query and the centroid
   // and add it to point_pq if it is nearer than the queue's
   // farthest member
   double d, dist_centroid = dist( cluster->centroid, query );
   nearest_neighbor_search_try_point( point_pq, cluster->centroid, dist_centroid );

   // Use the triangle inequality with the cluster radius to
   // determine if the entire cluster can be excluded as a
   // group
   if( heap_is_full( point_pq ) && dist_centroid >= point_pq->dists[0] + cluster->radius )
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
      if( heap_is_full( point_pq ) && dist_centroid > point_pq->dists[0] + index->dist )
         goto next_cluster_member;

      // Use the triangle inequality with the cluster member's
      // distance to centroid to determine if the point is
      // definitely nearer than the farthest member of point_pq
      if( dist_centroid <= point_pq->dists[0] - index->dist ) {
         d = dist( index->p, query );
         nearest_neighbor_search_try_point( point_pq, index->p, d );
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
               if( heap_is_full( point_pq ) && query_ancestors->dist > point_pq->dists[0] + cluster_ancestors->dist )
                  goto next_cluster_member;

               // Use the triangle inequality with the cluster member's
               // distance to ancestor to determine if the point is
               // definitely nearer than the farthest member of point_pq
               if( query_ancestors->dist <= point_pq->dists[0] - cluster_ancestors->dist ) {
                  d = dist( index->p, query );
                  nearest_neighbor_search_try_point( point_pq, index->p, d );
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
      nearest_neighbor_search_try_point( point_pq, index->p, d );

next_cluster_member:
      index = index->next;
   }
}


// Attempt to insert p into the point priority queue. Point
// p will be inserted if point_pq is not yet full. If
// point_pq is full, p will only be inserted if it is nearer
// than the farthest member of point_pq. In this case, the
// farthest member is removed to make room for p. Point p is
// never inserted if it already exists in point_pq. This
// function assumes that the items in point_pq are
// ap_Points. Returns true if p was inserted, or false
// otherwise.
bool
nearest_neighbor_search_try_point( ap_Heap *point_pq, ap_Point *p, double dist ) {

   // Return if p is already in the point priority queue
   int i;
   for( i = 0; i < point_pq->size; i++ )
      if( p == (ap_Point*)point_pq->items[i] )
         return false;

   // If the point priority queue is not full, insert p
   if( !heap_is_full( point_pq ) ) {
      return heap_insert( point_pq, p, dist );
   } else {

      // If the point priority queue is full and p is nearer than
      // the farthest member of point_pq, remove the farthest
      // member of point_pq and insert p
      if( dist < point_pq->dists[0] ) {
         heap_pop( point_pq );
         return heap_insert( point_pq, p, dist );
      } else {

         // If the point priority queue is full and p is not nearer
         // than the farthest member of point_pq, do not insert p
         return false;
      }
   }
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          GEOMETRIC MEDIAN AND ANTIPOLE FUNCTIONS
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

   ap_PointList *i, *j;
   double d, max_dist = -1;

   // Calculate the distance between each pair of points and
   // if a distance is greater than any found yet, make the
   // pair of points the new antipole pair
   for( i = set; i != NULL; i = i->next ) {
      for( j = i->next; j != NULL; j = j->next ) {
         d = dist( i->p, j->p );
         if( d > max_dist ) {
            *antipole_a = i->p;
            *antipole_b = j->p;
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
first_approx_antipoles( ap_PointList *set, ap_Point **antipole_a, ap_Point **antipole_b, double target_radius, DIST_FUNC ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   ap_PointList *i, *j;

   // Calculate the distance between each pair of points and
   // if a distance is greater than the target cluster diameter
   // make the pair of points the new antipole pair and stop
   // searching
   for( i = set; i != NULL; i = i->next ) {
      for( j = i->next; j != NULL; j = j->next ) {
         if( dist( i->p, j->p ) > 2 * target_radius ) {
            *antipole_a = i->p;
            *antipole_b = j->p;
            return;
         }
      }
   }
}


// Search the set for a point that when paired with ancestor
// could serve as an antipole pair.
void
check_ancestors_for_antipoles( ap_PointList *set, double target_radius, ap_Point *ancestor, ap_Point **antipole_a, ap_Point **antipole_b ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   ap_PointList *i, *j;

   // Search the set for a point whose distance to the ancestor
   // is greater than 2*target_radius and save it and the
   // ancestor as antipoles
   for( i = set; i != NULL; i = i->next ) {
      for( j = i->p->ancestors; j != NULL; j = j->next ) {
         if( j->p == ancestor ) {
            if( j->dist > 2 * target_radius ) {
               *antipole_a = ancestor;
               *antipole_b = i->p;
               return;
            }
         }
      }
   }
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                  LINKED LIST OPERATIONS
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


// Prepend to an ap_PointList an ap_Point with a distance
// value (to an ancestor, cluster centroid, or query,
// depending on the use of the ap_PointList). Requires the
// address of an ap_PointList pointer so that the list can
// be given a new first member. The ap_PointList pointer
// should be set to NULL before calling this function if the
// list is empty; otherwise the list may not terminate
// properly. Returns true if the point was added or false if
// it was not because the point already exists in the list.
bool
add_point( ap_PointList **list, ap_Point *p, double dist ) {

   // Return if the point is already in the list
   ap_PointList *index = *list;
   while( index != NULL ) {
      if( index->p == p )
         return false;
      index = index->next;
   }

   ap_PointList *new_list_member = malloc( sizeof( ap_PointList ) );
   assert( new_list_member );
   new_list_member->p = p;
   new_list_member->dist = dist;
   new_list_member->next = *list;
   *list = new_list_member;

   return true;
}


// Move the first instance of point p found in the
// ap_PointList *from into the ap_PointList *to. Requires
// the addresses of the ap_PointList pointers so that the
// moved member can be prepended to *to, and so that if the
// first member of ap_PointList *from is moved, the second
// member can become the new first member. The ap_PointList
// pointer *to should be set to NULL before calling this
// function if the list is empty; otherwise the list may not
// terminate properly. Returns true if the point was moved
// or false if it was not found.
bool
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
      return false;

   // Remove index from *from
   if( i == 0 )
      *from = index->next;
   else
      before->next = index->next;

   // Prepend index to *to
   index->next = *to;
   *to = index;

   return true;
}


// Move the n-th member (start counting at zero) of the
// ap_PointList *from into the ap_PointList *to. Requires
// the addresses of the ap_PointList pointers so that the
// moved member can be prepended to *to, and so that if the
// first member of ap_PointList *from is moved, the second
// member can become the new first member. The ap_PointList
// pointer *to should be set to NULL before calling this
// function if the list is empty; otherwise the list may not
// terminate properly. Returns true if the point was moved
// or false if n was too large.
bool
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
      return false;

   // Remove index from *from
   if( n == 0 )
      *from = index->next;
   else
      before->next = index->next;

   // Prepend index to *to
   index->next = *to;
   *to = index;

   return true;
}


// Create a new ap_PointList with the same contents as list
// (the order will be reversed).
ap_PointList*
copy_list( ap_PointList *list ) {

   ap_PointList *new_list = NULL;
   while( list != NULL ) {
      add_point( &new_list, list->p, list->dist );
      list = list->next;
   }
   return new_list;
}


// Find the size of the list of points
int
list_size( ap_PointList *list ) {

   int size = 0;
   while( list != NULL ) {
      size++;
      list = list->next;
   }
   return size;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                     HEAP OPERATIONS
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


// Create a new ap_Heap object. The heap is a min-heap if
// is_max_heap is false or a max-heap otherwise. If
// max_size < 1, the heap will have no maximum size.
ap_Heap*
create_heap( bool is_max_heap, int max_size ) {

   ap_Heap *heap = malloc( sizeof( ap_Heap ) );
   assert( heap );
   heap->is_max_heap = is_max_heap;
   heap->max_size = max_size;
   heap->size = 0;
   if( max_size < 1 )
      heap->capacity = 10;
   else
      heap->capacity = max_size;
   heap->items = calloc( heap->capacity, sizeof( void* ) );
   heap->dists = calloc( heap->capacity, sizeof( double ) );
   assert( heap->items && heap->dists );

   return heap;
}


// Return false if the heap has no maximum size or if the
// current heap size is smaller than the maximum size, or
// return true otherwise.
bool
heap_is_full( ap_Heap *heap ) {

   if( heap->max_size < 1 || heap->size < heap->max_size )
      return false;
   else
      return true;
}


// Increase the capacity of heap
void
heap_grow( ap_Heap *heap ) {

   heap->capacity *= 2;
   void  **temp_items = calloc( heap->capacity, sizeof( void* ) );
   double *temp_dists = calloc( heap->capacity, sizeof( double ) );
   assert( temp_items && temp_dists );

   int i;
   for( i = 0; i < heap->size; i++ ) {
      temp_items[i] = heap->items[i];
      temp_dists[i] = heap->dists[i];
   }

   free( heap->items );
   free( heap->dists );
   heap->items = temp_items;
   heap->dists = temp_dists;
}


// Swap items at positions i and j in heap. Returns true if
// a swap occurred or false if the indices were
// inapproriate.
bool
heap_swap( ap_Heap *heap, int i, int j ) {

   // Return if the indices match or are out of bounds
   if( i == j || i >= heap->size || j >= heap->size || i < 0 || j < 0 )
      return false;

   void *temp_item  = heap->items[i];
   double temp_dist = heap->dists[i];
   heap->items[i] = heap->items[j];
   heap->dists[i] = heap->dists[j];
   heap->items[j] = temp_item;
   heap->dists[j] = temp_dist;

   return true;
}


// Sift upward the item at position index in heap until it
// is properly positioned.
void
heap_sift_up( ap_Heap* heap, int index ) {

   int parent;

   // Repeatedly sift the item upward, stopping if it
   // becomes the root
   while( index > 0 ) {

      parent = ( index - 1 ) / 2;

      // If index has a greater/lesser dist than its parent and
      // the heap is a max-/min-heap, swap index and its parent
      // and continue sifting upward
      if( ( heap->is_max_heap && heap->dists[index] > heap->dists[parent] ) ||
         ( !heap->is_max_heap && heap->dists[index] < heap->dists[parent] ) ) {
         heap_swap( heap, index, parent );
         index = parent;
      } else {

         // Stop if index is properly positioned
         break;
      }
   }
}


// Sift downward the item at position index in heap until it
// is properly positioned.
void
heap_sift_down( ap_Heap *heap, int index ) {

   int left, right, child;

   // Repeatedly sift the item downward, stopping if it becomes
   // a leaf
   while( index < heap->size / 2 ) {

      left  = 2 * index + 1;
      right = 2 * index + 2;

      // If index only has one child, find its index
      if( right >= heap->size ) {
         child = left;
      } else {

         // Find the child with the greater/lesser dist if the heap
         // is a max-/min-heap
         if( ( heap->is_max_heap && heap->dists[left] > heap->dists[right] ) ||
            ( !heap->is_max_heap && heap->dists[left] < heap->dists[right] ) )
            child = left;
         else
            child = right;
      }

      // If index has a lesser/greater dist than one of its
      // children and the heap is a max-/min-heap, swap index and
      // the greater/lesser child and continue sifting upward
      if( ( heap->is_max_heap && heap->dists[index] < heap->dists[child] ) ||
         ( !heap->is_max_heap && heap->dists[index] > heap->dists[child] ) ) {
         heap_swap( heap, index, child );
         index = child;
      } else {

         // Stop if index is properly positioned
         break;
      }
   }
}


// Add an item to an ap_Heap with a distance value used for
// sorting. The item will be inserted only if the heap does
// not already contain the maximum number of items. Returns
// true if the item was inserted or false if the heap was
// full.
bool
heap_insert( ap_Heap *heap, void *item, double dist ) {

   int i;

   // Add the item to the heap only if it is not full
   if( !heap_is_full( heap ) ) {
      
      // Grow the data arrays if necessary
      if( heap->size == heap->capacity )
         heap_grow( heap );

      // Insert the new item at the end of the array
      i = heap->size;
      heap->items[i] = item;
      heap->dists[i] = dist;
      heap->size++;

      // Sift the new item upward
      heap_sift_up( heap, i );

      return true;
   } else {

      // Return if the heap has reached its maximum size
      return false;
   }
}


// Remove the first item in the heap and return a pointer to
// it. Returns NULL if the heap is empty.
void*
heap_pop( ap_Heap *heap ) {

   // Return if the heap is empty
   if( heap->size == 0 )
      return NULL;

   void* first_item = heap->items[0];

   // Replace the root with the last item in the list
   heap->items[0] = heap->items[heap->size-1];
   heap->dists[0] = heap->dists[heap->size-1];
   heap->size--;

   // Sift the moved item downward
   heap_sift_down( heap, 0 );

   return first_item;
}


// Create an ap_PointList from an ap_Heap. The points in the
// list will be sorted by distance in ascending/descending
// order if the heap is a max-/min-heap. This function
// assumes the items in the heap are ap_Points.
ap_PointList*
heap_to_list( ap_Heap *heap ) {

   int i;
   double dist;
   ap_PointList *new_list = NULL;

   // Create a copy of the heap (because it will be destroyed
   // below)
   ap_Heap *new_heap = malloc( sizeof( ap_Heap ) );
   assert( new_heap );
   new_heap->is_max_heap = heap->is_max_heap;
   new_heap->max_size = heap->max_size;
   new_heap->size = heap->size;
   new_heap->capacity = heap->capacity;
   new_heap->items = calloc( new_heap->capacity, sizeof( void* ) );
   new_heap->dists = calloc( new_heap->capacity, sizeof( double ) );
   assert( new_heap->items && new_heap->dists );
   for( i = 0; i < heap->size; i++ ) {
      new_heap->items[i] = heap->items[i];
      new_heap->dists[i] = heap->dists[i];
   }

   // Copy the points from the heap into a list (the sequence
   // of calls to add_point reverses the order)
   while( new_heap->size > 0 ) {
      dist = new_heap->dists[0];
      add_point( &new_list, (ap_Point*)heap_pop( new_heap ), dist );
   }
   free_heap( new_heap );

   return new_list;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
               MEMORY MANAGEMENT FUNCTIONS
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


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
free_list( ap_PointList *list ) {

   if( list != NULL ) {
      free_list( list->next );
      free( list );
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


