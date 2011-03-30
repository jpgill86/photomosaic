#include <assert.h>  /* assert */
#include <math.h>    /* fmax */
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
struct ap_Tree*
build_tree( struct ap_List *set, double target_radius, struct ap_List *antipoles, DIST_FUNC ) {

   // Create the new ap_Tree
   struct ap_Tree *new_tree = malloc( sizeof( struct ap_Tree ) );
   assert( new_tree );
      
   // Determine if this tree is an internal node or a leaf
   if( antipoles == NULL ) {
      antipoles = adapted_approx_antipoles( set, target_radius );
      if( antipoles == NULL ) {
         // If it is a leaf, create a cluster from the set and return
         // the leaf
         new_tree->is_leaf = 1;
         new_tree->cluster = make_cluster( set, dist );
         return new_tree;
      }
   }

   // If this tree is an internal node, initialize it
   new_tree->is_leaf = 0;
   new_tree->a = antipoles->p;
   new_tree->b = antipoles->next->p;
   new_tree->radius_a = 0;
   new_tree->radius_b = 0;

   // For each point in the set, find the distance to each
   // antipole, store the distances in the point's ancestor
   // list, add the point to the subset belonging to the
   // nearest antipole, and update the radius of the subset if
   // necessary
   double dist_a, dist_b;
   struct ap_List *index, *set_a = NULL, *set_b = NULL;
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
   free( set );

   // Build subtrees as children for this node using the two
   // point subsets
   new_tree->left  = build_tree( set_a, target_radius, check( set_a, target_radius, new_tree->a ), dist );
   new_tree->right = build_tree( set_b, target_radius, check( set_b, target_radius, new_tree->b ), dist );

   return new_tree;
}


// Create an ap_Cluster owned by a leaf of the tree data
// structure containing a list of the points in the cluster
// (already determined to be sufficiently close to one
// another to group together), the identity of the geometric
// median of the cluster, and the cluster radius.
struct ap_Cluster*
make_cluster( struct ap_List *set, DIST_FUNC ) {

   struct ap_List *index;
   double dist_centroid;

   // Create the new ap_Cluster and initialize it
   struct ap_Cluster *new_cluster = malloc( sizeof( struct ap_Cluster ) );
   assert( new_cluster );
   new_cluster->centroid = approx_1_median( set, dist );
   //new_cluster->size = 0;
   new_cluster->radius = 0;
   new_cluster->members = NULL;

   // For every point in the set (besides the centroid), find
   // the distance to the centroid, add the point to the list
   // of points in the cluster, and update the radius of the
   // cluster if necessary
   for( index = set; index != NULL; index = index->next ) {
      if( index->p != new_cluster->centroid ) {
         dist_centroid = dist( new_cluster->centroid, index->p );
         //index->p->dist_centroid = dist_centroid;
         add_point( &(new_cluster->members), index->p, dist_centroid );
         //new_cluster->size++;
         new_cluster->radius = fmax( new_cluster->radius, dist_centroid );
      }
   }

   free( set );
   return new_cluster;
}


// Prepend to an ap_List an ap_Point with a distance value
// (to an ancestor, cluster centroid, or query, depending on
// the use of the ap_List). Requires the address of an
// ap_List pointer so that the list can be given a new first
// member. The ap_List pointer should be set to NULL before
// calling this function if the list is empty; otherwise the
// list may not terminate properly.
void
add_point( struct ap_List **set, struct ap_Point *p, double dist ) {

   struct ap_List *new_list_member = malloc( sizeof( struct ap_List ) );
   assert( new_list_member );
   new_list_member->p = p;
   new_list_member->dist = dist;
   new_list_member->next = *set;
   *set = new_list_member;
}


// Create a new ap_List with the same contents as set.
struct ap_List*
copy_list( struct ap_List *set ) {

   struct ap_List *index = set, *new_list = NULL;
   while( index != NULL ) {
      add_point( &new_list, index->p, index->dist );
      index = index->next;
   }
   return new_list;
}


// Move the n-th member (start counting at zero) of the
// ap_List *from into the ap_List *to. Requires the
// addresses of the ap_List pointers so that the moved
// member can be prepended to *to, and so that if the first
// member of ap_List *from is moved, the second member can
// become the new first member. The ap_List pointer *to
// should be set to NULL before calling this function if the
// list is empty; otherwise the list may not terminate
// properly.
void
move_list( int n, struct ap_List **from, struct ap_List **to ) {

   int i;
   struct ap_List *before = NULL, *index = *from;

   // Find the ap_List to be moved in match, and store the
   // links to and from it for later user
   for( i = 0; i < n; i++ ) {
      before = index;
      index = index->next;
   }

   // Remove index from *from
   if( n == 0 )
      *from = index->next;
   else
      before->next = index->next;

   // Prepend index to *to
   index->next = *to;
   *to = index;
}


// Recursively free up memory used by an ap_List linked
// list.
void
free_list( struct ap_List *set ) {

   if( set != NULL ) {
      free_list( set->next );
      free( set );
   }
}


// Find the size of a set of points
int
set_size( struct ap_List *set ) {

   int size = 0;
   while( set != NULL ) {
      size++;
      set = set->next;
   }
   return size;
}


// Find the exact geometric median of a set of points
struct ap_Point*
exact_1_median( struct ap_List *set, DIST_FUNC ) {

   int i, j, d, size = set_size( set );
   struct ap_List *i_list, *j_list;

   // Initialize the array of distance sums
   double sums[ size ];
   for( i = 0; i < size; i++ )
      sums[ i ] = 0;

   // Calculate the distance between each pair of points and
   // add each distance to the distance sums for each point in
   // the pair
   for( i = 0, i_list = set; i < size; i++, i_list = i_list->next ) {
      for( j = i + 1, j_list = i_list->next; j_list != NULL; j++, j_list = j_list->next ) {
         d = dist( i_list->p, j_list->p );
         sums[ i ] += d;
         sums[ j ] += d;
      }
   }

   // Identify the point with the minimum distance sum
   double min_sum = -1;
   struct ap_Point *median = NULL;
   for( i = 0, i_list = set; i < size; i++, i_list = i_list->next ) {
      if( sums[ i ] < min_sum || min_sum < 0 ) {
         min_sum = sums[ i ];
         median = i_list->p;
      }
   }

   return median;
}


// Find an approximation for the geometric median of a set
// of points. The user should initialize the random number
// generator using srand.
struct ap_Point*
approx_1_median( struct ap_List *set, DIST_FUNC ) {

   struct ap_Point *median;
   struct ap_List *index, *contestants = copy_list( set ), *tournament, *winners;
   int i, contestants_size = set_size( contestants ), tournament_size = 3, winners_size;
   int final_round_size = min( pow( tournament_size, 2 ) - 1, round( sqrt( set_size( set ) ) ) );

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
            move_list( rand() % contestants_size, &contestants, &tournament );
            contestants_size--;
         }
         // Find the winner of this tournament and discard the losers
         median = exact_1_median( tournament, dist );
         for( i = 0, index = tournament; i < tournament_size; i++, index = index->next ) {
            if( index->p == median )
               break;
         }
         move_list( i, &tournament, &winners );
         winners_size++;
         free_list( tournament );
      }
      // Find the winner among the remaining contestants and
      // discard the losers
      median = exact_1_median( contestants, dist );
      for( i = 0, index = contestants; i < contestants_size; i++, index = index->next ) {
         if( index->p == median )
            break;
      }
      move_list( i, &contestants, &winners );
      winners_size++;
      free_list( contestants );

      // Fill the pool of contestants with all the winners in
      // preparation for the next round
      contestants = winners;
      contestants_size = winners_size;
   }
   
   // Find the overall winner and discard the losers
   median = exact_1_median( contestants, dist );
   free_list( contestants );

   return median;
}


// Find the two points in the set that are furthest from one
// another and place the indices in ap1 and ap2.
void
exact_antipoles( struct ap_List *set, int *ap1, int *ap2, DIST_FUNC ) {

   int i, j;
   struct ap_List *i_list, *j_list;
   double d, max_dist = 0;

   // Calculate the distance between each pair of points and
   // if a distance is greater than any found yet, make the
   // pair of points the new antipole pair
   for( i = 0, i_list = set; i_list != NULL; i++, i_list = i_list->next ) {
      for( j = i + 1, j_list = i_list->next; j_list != NULL; j++, j_list = j_list->next ) {
         d = dist( i_list->p, j_list->p );
         if( d > max_dist ) {
            *ap1 = i;
            *ap2 = j;
            max_dist = d;
         }
      }
   }
}


// Find an approximation for the antipole pair of a set
// of points and place the indices in ap1 and ap2. The user
// should initialize the random number generator using
// srand.
void
approx_antipoles( struct ap_List *set, int *ap1, int *ap2, DIST_FUNC ) {

   struct ap_Point *p1, *p2;
   struct ap_List *index, *contestants = copy_list( set ), *tournament, *winners;
   int i, contestants_size = set_size( contestants ), tournament_size = 3, winners_size;
   int final_round_size = min( pow( tournament_size, 2 ) - 1, round( sqrt( set_size( set ) ) ) );

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
            move_list( rand() % contestants_size, &contestants, &tournament );
            contestants_size--;
         }
         // Find the winners of this tournament and discard the losers
         exact_antipoles( tournament, ap1, ap2, dist );
         move_list( *ap2, &tournament, &winners );
         move_list( *ap1, &tournament, &winners );
         winners_size += 2;
         free_list( tournament );
      }
      // Find the winners among the remaining contestants and
      // discard the losers
      exact_antipoles( contestants, ap1, ap2, dist );
      move_list( *ap2, &contestants, &winners );
      move_list( *ap1, &contestants, &winners );
      winners_size += 2;
      free_list( contestants );

      // Fill the pool of contestants with all the winners in
      // preparation for the next round
      contestants = winners;
      contestants_size = winners_size;
   }
   
   // Find the overall winners and discard the losers
   exact_antipoles( contestants, ap1, ap2, dist );
   for( i = 0, index = contestants; i < *ap1; i++ )
      index = index->next;
   p1 = index->p;
   for( i = 0, index = contestants; i < *ap2; i++ )
      index = index->next;
   p2 = index->p;
   for( i = 0, index = set; index != NULL; i++, index = index->next )
      if( index->p == p1 )
         *ap1 = i;
   for( i = 0, index = set; index != NULL; i++, index = index->next )
      if( index->p == p2 )
         *ap2 = i;
   free_list( contestants );
}


// ...
struct ap_List*
adapted_approx_antipoles( struct ap_List *set, double target_radius ) {
}


// ...
struct ap_List*
check( struct ap_List *set, double target_radius, struct ap_Point *antipole ) {
}


