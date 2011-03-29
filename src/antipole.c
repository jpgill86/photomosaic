#include <assert.h>  /* assert */
#include <math.h>    /* fmax */
#include <stdlib.h>  /* NULL, rand */
#include "antipole.h"


// Create an ap_Tree that serves as the root, an internal
// node, or a leaf for the tree data structure. Non-leaves
// contain the identities of two antipole points, a left
// subtree and a right subtree which each contain the
// subset of points that is nearest its respective antipole
// point, and the radii of the subsets. Leaves contain a
// cluster of points.
struct ap_Tree*
build_tree( struct ap_List *set, double target_radius, struct ap_List *antipoles, DIST_FUNC ) {

   // Allocate space for the new ap_Tree
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

   new_tree->is_leaf = 0;

   // Store pointers to the antipoles in this set
   new_tree->a = antipoles->p;
   new_tree->b = antipoles->next->p;

   // Process the members of set
   struct ap_List *index = set, *set_a = NULL, *set_b = NULL;
   double dist_a, dist_b;
   new_tree->radius_a = 0;
   new_tree->radius_b = 0;
   while( index != NULL ) {
      // For each point in the set, find the distance to each
      // antipole, store the distances in the point's ancestor
      // list, add the point to the subset belonging to the
      // nearest antipole, and update the radius of the subset
      // if necessary
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

      index = index->next;
   }

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

   // Allocate space for the new ap_Cluster
   struct ap_Cluster *new_cluster = malloc( sizeof( struct ap_Cluster ) );
   assert( new_cluster );

   // Find the geometric median of the cluster
   new_cluster->centroid = approx_1_median( set, dist );

   // Process the members of the set
   struct ap_List *index = set;
   new_cluster->members = NULL;
   double dist_centroid;
   //new_cluster->size = 0;
   new_cluster->radius = 0;
   while( index != NULL ) {
      // For every point in the set (besides the centroid), find
      // the distance to the centroid, add the point to the list
      // of points in the cluster, and update the radius of the
      // cluster if necessary
      if( index->p != new_cluster->centroid ) {
         dist_centroid = dist( new_cluster->centroid, index->p );
         //index->p->dist_centroid = dist_centroid;
         add_point( &(new_cluster->members), index->p, dist_centroid );
         //new_cluster->size++;
         new_cluster->radius = fmax( new_cluster->radius, dist_centroid );
      }

      index = index->next;
   }

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

   int count = 0;
   struct ap_List *index = *from;
   struct ap_List *before = NULL, *match = index, *after = index->next;

   // Find the ap_List to be moved in match, and store the
   // links to and from it for later user
   while( count < n ) {
      before = index;
      index = index->next;
      match = index;
      after = index->next;
      count++;
   }

   // Remove match from *from
   if( *from == NULL || n == 0 )
      *from = after;
   else
      before->next = after;

   // Prepend match to *to
   match->next = *to;
   *to = match;
}


// Recursively free up memory used by an ap_List linked
// list.
void
free_list( struct ap_List *set ) {

   if( set != NULL ) {
      if( set->next != NULL )
         free_list( set->next );
      free( set );
   }
}


// Find the size of a set of points
int
set_size( struct ap_List *set ) {

   int size = 0;
   struct ap_List *index = set;
   while( index != NULL ) {
      size++;
      index = index->next;
   }
   return size;
}


// Find the exact geometric median of a set of points
struct ap_Point*
exact_1_median( struct ap_List *set, DIST_FUNC ) {

   struct ap_Point *median = NULL;
   double min_sum_dist = -1;
   double sum_dist;

   // Process the members of the set
   struct ap_List *i = set;
   while( i != NULL ) {
      // For each point in the set, sum its distances from every
      // other point
      sum_dist = 0;
      struct ap_List *j = set;
      while( j != NULL ) {
         if( i->p != j->p )
            sum_dist += dist( i->p, j->p );
         j = j->next;
      }

      // If the sum for a point is found to be the new minimum,
      // store the value and the identity of the candidate for the
      // geometric median
      if( min_sum_dist < 0 || sum_dist < min_sum_dist ) {
         min_sum_dist = sum_dist;
         median = i->p;
      }

      i = i->next;
   }

   return median;
}


// Find an approximation for the geometric median of a set
// of points. The user should initialize the random number
// generator using srand.
struct ap_Point*
approx_1_median( struct ap_List *set, DIST_FUNC ) {

   int i;
   int final_round_size = 3;
   int tournament_size = 3;
   struct ap_List *contestants = copy_list( set );

   // Hold a series of rounds of tournaments
   while( set_size( contestants ) > final_round_size ) {
      // Find the winners that will continue to the next round
      struct ap_List *winners = NULL;
      while( set_size( contestants ) >= 2 * tournament_size ) {
         struct ap_List *tournament = NULL;
         // Move tournament_size random members of contestants into
         // tournament
         for( i = 0; i < tournament_size; i++ )
            move_list( rand() % set_size( contestants ), &contestants, &tournament );
         // Find the winner of this tournament and discard the losers
         add_point( &winners, exact_1_median( tournament, dist ), 0 );
         free_list( tournament );
      }
      // Find the winner among the remaining contestants and
      // discard the losers
      add_point( &winners, exact_1_median( contestants, dist ), 0 );
      free_list( contestants );

      // Fill the pool of contestants with all the winners in
      // preparation for the next round
      contestants = winners;
   }
   
   // Find the overall winner and discard the losers
   struct ap_Point *median = exact_1_median( contestants, dist );
   free_list( contestants );

   return median;
}


// Find the two points in the set that are furthest from one
// another and return these as an antipole pair
struct ap_List*
exact_antipoles( struct ap_List *set, DIST_FUNC ) {

   struct ap_List *antipoles = NULL;
   double max_dist = 0;
   double d;

   // Create the antipole list if there are at least two
   // members in the set
   if( set != NULL ) {
      if( set->next != NULL ) {
         antipoles = malloc( sizeof( struct ap_List ) );
         antipoles->p = NULL;
         antipoles->dist = 0;
         antipoles->next = malloc( sizeof( struct ap_List ) );
         antipoles->next->p = NULL;
         antipoles->next->dist = 0;
         antipoles->next->next = NULL;
      }
   }

   // Process the members of the set
   struct ap_List *i = set, *j;
   while( i != NULL ) {
      j = i->next;
      while( j != NULL ) {
         // Calulate the distance between the two points, and if
         // their distance is the greatest found so far, make them
         // the new antipole pair
         d = dist( i->p, j->p );
         if( d > max_dist ) {
            max_dist = d;
            antipoles->p = i->p;
            antipoles->dist = max_dist;
            antipoles->next->p = j->p;
            antipoles->next->dist = max_dist;
         }
         j = j->next;
      }
      i = i->next;
   }

   return antipoles;
}


// ...
struct ap_List*
approx_antipoles( struct ap_List *set, DIST_FUNC ) {
}


// ...
struct ap_List*
adapted_approx_antipoles( struct ap_List *set, double target_radius ) {
}


// ...
struct ap_List*
check( struct ap_List *set, double target_radius, struct ap_Point *antipole ) {
}


