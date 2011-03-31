#include <assert.h>  /* assert */
#include <math.h>    /* fmax */
#include <stdlib.h>  /* NULL, rand */
#include "antipole.h"



#include <stdio.h>

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
build_tree( int level, ap_List *set, double target_radius, ap_Point *antipole_a, ap_Point *antipole_b, int dimensionality, DIST_FUNC ) {
   //printf("BUILD TREE --- %d\n", level);
   //printf(" |- size = %d\n", list_size(set));

   // Create the new ap_Tree
   ap_Tree *new_tree = malloc( sizeof( ap_Tree ) );
   assert( new_tree );

   // Determine if this tree is an internal node or a leaf
   if( antipole_a == NULL || antipole_b == NULL ) {
      adapted_approx_antipoles( set, &antipole_a, &antipole_b, target_radius, dist );
      if( antipole_a == NULL || antipole_b == NULL ) {
         //printf(" |- splitting condition not satisfied\n");
         // If it is a leaf, create a cluster from the set and return
         // the leaf
         new_tree->is_leaf = 1;
         new_tree->cluster = make_cluster( set, dimensionality, dist );
         return new_tree;
      }
   }

   //printf(" |- splitting condition satisfied: d = %f\n", dist(antipole_a,antipole_b));
   // If this tree is an internal node, initialize it
   new_tree->is_leaf = 0;
   new_tree->a = antipole_a;
   new_tree->b = antipole_b;
   new_tree->radius_a = 0;
   new_tree->radius_b = 0;
   //printf(" |- antipoles have id=%d and id=%d\n", antipole_a->id, antipole_b->id);

   // For each point in the set, find the distance to each
   // antipole, store the distances in the point's ancestor
   // list, add the point to the subset belonging to the
   // nearest antipole, and update the radius of the subset if
   // necessary
   double dist_a, dist_b;
   ap_List *index, *set_a = NULL, *set_b = NULL;
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
   //printf(" |- a_set has %d\n", list_size(set_a));
   //printf(" |- b_set has %d\n", list_size(set_b));

   // Build subtrees as children for this node using the two
   // point subsets
   level++;
   check_for_antipoles( set_a, target_radius, new_tree->a, &antipole_a, &antipole_b );
   //printf(" |- building child a\n");
   new_tree->left = build_tree( level, set_a, target_radius, antipole_a, antipole_b, dimensionality, dist );
   check_for_antipoles( set_b, target_radius, new_tree->b, &antipole_a, &antipole_b );
   //printf(" |- building child b\n");
   new_tree->right = build_tree( level, set_b, target_radius, antipole_a, antipole_b, dimensionality, dist );

   return new_tree;
}


// Create an ap_Cluster owned by a leaf of the tree data
// structure containing a list of the points in the cluster
// (already determined to be sufficiently close to one
// another to group together), the identity of the geometric
// median of the cluster, and the cluster radius.
ap_Cluster*
make_cluster( ap_List *set, int dimensionality, DIST_FUNC ) {
   //printf(" |- BUILD CLUSTER\n");
   //printf("    |- size = %d\n", list_size(set));

   ap_List *index;
   double dist_centroid;

   // Create the new ap_Cluster and initialize it
   ap_Cluster *new_cluster = malloc( sizeof( ap_Cluster ) );
   assert( new_cluster );
   approx_1_median( set, &(new_cluster->centroid), dimensionality, dist );
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
add_point( ap_List **set, ap_Point *p, double dist ) {

   // Check for uniqueness
   ap_List *index = *set;
   while( index != NULL ) {
      if( index->p == p )
         return;
      index = index->next;
   }

   // Add the point if it is unique to the list
   ap_List *new_list_member = malloc( sizeof( ap_List ) );
   assert( new_list_member );
   new_list_member->p = p;
   new_list_member->dist = dist;
   new_list_member->next = *set;
   *set = new_list_member;
}


// Move the first instance of point p found in the ap_List
// *from into the ap_List *to. Requires the addresses of the
// ap_List pointers so that the moved member can be
// prepended to *to, and so that if the first member of
// ap_List *from is moved, the second member can become the
// new first member. The ap_List pointer *to should be set
// to NULL before calling this function if the list is
// empty; otherwise the list may not terminate properly.
void
move_point( ap_Point *p, ap_List **from, ap_List **to ) {

   int i = 0;
   ap_List *before = NULL, *index = *from;

   // Find the first ap_List containing p, and store the link
   // preceding it for later user
   while( index != NULL && index->p != p ) {
      before = index;
      index = index->next;
      i++;
   }

   // Continue if p was found
   if( index != NULL ) {
      // Remove index from *from
      if( i == 0 )
         *from = index->next;
      else
         before->next = index->next;

      // Prepend index to *to
      index->next = *to;
      *to = index;
   }
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
move_nth_point( int n, ap_List **from, ap_List **to ) {

   int i;
   ap_List *before = NULL, *index = *from;

   // Find the n-th ap_List, and store the link preceding it
   // for later user
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


// Create a new ap_List with the same contents as set.
ap_List*
copy_list( ap_List *set ) {

   ap_List *index = set, *new_list = NULL;
   while( index != NULL ) {
      add_point( &new_list, index->p, index->dist );
      index = index->next;
   }
   return new_list;
}


// Recursively free up memory used by an ap_List linked
// list.
void
free_list( ap_List *set ) {

   if( set != NULL ) {
      free_list( set->next );
      free( set );
   }
}


// Find the size of a set of points
int
list_size( ap_List *set ) {

   int size = 0;
   while( set != NULL ) {
      size++;
      set = set->next;
   }
   return size;
}


// Find the exact geometric median of a set of points and
// store it in median
void
exact_1_median( ap_List *set, ap_Point **median, DIST_FUNC ) {

   *median = NULL;

   int i, j, d, size = list_size( set );
   ap_List *i_list, *j_list;

   // Initialize the array of distance sums
   double sums[ size ];
   for( i = 0; i < size; i++ )
      sums[ i ] = 0;

   // Calculate the distance between each pair of points and
   // add each distance to the distance sums for each point in
   // the pair
   for( i = 0, i_list = set; i < size; i++, i_list = i_list->next ) {
      for( j = i + 1, j_list = i_list->next; j < size; j++, j_list = j_list->next ) {
         d = dist( i_list->p, j_list->p );
         sums[ i ] += d;
         sums[ j ] += d;
      }
   }

   // Identify the point with the minimum distance sum and
   // store it in median
   double min_sum = -1;
   for( i = 0, i_list = set; i < size; i++, i_list = i_list->next ) {
      if( sums[ i ] < min_sum || min_sum < 0 ) {
         min_sum = sums[ i ];
         *median = i_list->p;
      }
   }
}


// Find an approximation for the geometric median of a set
// of points and store it in median. The user should
// initialize the random number generator using srand.
void
approx_1_median( ap_List *set, ap_Point **median, int dimensionality, DIST_FUNC ) {

   *median = NULL;

   ap_List *index, *contestants = copy_list( set ), *tournament, *winners;
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


// Find the two points in the set that are furthest from one
// another and store them in antipole_a and antipole_b.
void
exact_antipoles( ap_List *set, ap_Point **antipole_a, ap_Point **antipole_b, DIST_FUNC ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   ap_List *i_list, *j_list;
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
approx_antipoles( ap_List *set, ap_Point **antipole_a, ap_Point **antipole_b, int dimensionality, DIST_FUNC ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   ap_List *index, *contestants = copy_list( set ), *tournament, *winners;
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
adapted_approx_antipoles( ap_List *set, ap_Point **antipole_a, ap_Point **antipole_b, double target_radius, DIST_FUNC ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   ap_List *i_list, *j_list;

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
check_for_antipoles( ap_List *set, double target_radius, ap_Point *ancestor, ap_Point **antipole_a, ap_Point **antipole_b ) {

   *antipole_a = NULL;
   *antipole_b = NULL;

   // Search the set for a point whose distance to the ancestor
   // is greater than 2*target_radius and save it and the
   // ancestor as antipoles
   ap_List *index_i, *index_j;
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


// Search the tree recursively for all points within range
// of query and place them in out.
void
range_search( ap_Tree *tree, ap_Point *query, double range, ap_List **out, DIST_FUNC ) {
   //printf("ENTERING TREE\n");

   if( !tree->is_leaf ) {
      // Calculate the distance between query and the antipoles
      // and store these values in the ancestor list for query
      double dist_a = dist( tree->a, query );
      double dist_b = dist( tree->b, query );
      add_point( &(query->ancestors), tree->a, dist_a );
      add_point( &(query->ancestors), tree->b, dist_b );

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

      // Once the search has returned from both subtrees, remove
      // the antipoles of tree from the ancestor list for query,
      // since they will no longer be relevant to further searches
      ap_List *trash = NULL;
      move_point( tree->a, &(query->ancestors), &trash );
      move_point( tree->b, &(query->ancestors), &trash );
      free_list( trash );
   } else {
      // If tree is a leaf, search its cluster for points within
      // range of query
      range_visit_cluster( tree->cluster, query, range, out, dist );
   }
}


// Find the members of the cluster that are within range of
// query and add them to the list out.
void
range_visit_cluster( ap_Cluster *cluster, ap_Point *query, double range, ap_List **out, DIST_FUNC ) {
   //printf("VISIT CLUSTER: %p\n", cluster);

   // Calculate the distance between the query and the centroid
   // and add it to out if it is within range
   double d, dist_centroid = dist( cluster->centroid, query );
   if( dist_centroid <= range ) {
      //printf(" |- centroid id=%d within range (%f)\n", cluster->centroid->id, dist_centroid);
      add_point( out, cluster->centroid, dist_centroid );
   } else {
      //printf(" |- centroid out of range (%f)\n", dist_centroid);
   }

   //printf(" |- list_size members = %d\n", list_size(cluster->members));
   //printf(" |- cluster radius = %f\n", cluster->radius);

   // Use the triangle inequality with the cluster radius to
   // determine if the entire cluster can be excluded as a
   // group
   if( dist_centroid >= range + cluster->radius ) {
      //printf(" |- ALL EXCLUDED\n");
      return;
   }

   // Use the triangle inequality with the cluster radius to
   // determine if the entire cluster can be included as a
   // group
   if( dist_centroid <= range - cluster->radius ) {
      //printf(" |- ALL INCLUDED\n");
      ap_List *index = cluster->members;
      while( index != NULL ) {
         //printf(" |- id=%d\n", index->p->id);
         add_point( out, index->p, -1 );
         index = index->next;
      }
      return;
   }

   // Check each member of the cluster
   ap_List *index = cluster->members;
   ap_List *query_ancestors, *cluster_ancestors;
   //printf(" |- MEMBERS\n");
   while( index != NULL ) {
      //printf("    |- index = %p\n", index);
      // Use the triangle inequality with the cluster member's
      // distance to centroid to determine if the point is
      // definitely out of range
      if( dist_centroid >= range + index->dist )
         goto next_cluster_member;

      // Use the triangle inequality with the cluster member's
      // distance to centroid to determine if the point is
      // definitely within range
      if( dist_centroid <= range - index->dist ) {
         //printf("    |- tri_c id=%d\n", index->p->id);
         add_point( out, index->p, -1 );
         goto next_cluster_member;
      }

      // Check the ancestors of the query and the member of the
      // cluster
      //printf("    |- MEMBER ANCESTORS\n");
      for( query_ancestors = query->ancestors; query_ancestors != NULL; query_ancestors = query_ancestors->next ) {
         for( cluster_ancestors = index->p->ancestors; cluster_ancestors != NULL; cluster_ancestors = cluster_ancestors->next ) {
            //printf("       |- q_a id=%d  c_a id=%d\n", query_ancestors->p->id, cluster_ancestors->p->id);

            // Use the triangle inequality with the cluster member's
            // distance to ancestor to determine if the point is
            // definitely out of range
            if( query_ancestors->dist >= range + cluster_ancestors->dist )
               goto next_cluster_member;

            // Use the triangle inequality with the cluster member's
            // distance to ancestor to determine if the point is
            // definitely within range
            if( query_ancestors->dist <= range - cluster_ancestors->dist ) {
               //printf("    |- tri_anc id=%d\n", index->p->id);
               add_point( out, index->p, -1 );
               goto next_cluster_member;
            }
         }
      }

      // Finally, if all methods of using precalculated distances
      // to rule-out or rule-in the cluster member have failed,
      // calculate the distance between the query and the cluster
      // member and add it to out if it is within range
      d = dist( index->p, query );
      if( d <= range ) {
         //printf("    |- calc id=%d\n", index->p->id);
         add_point( out, index->p, d );
         goto next_cluster_member;
      }

next_cluster_member:
      index = index->next;
   }

   //printf(" |- EXIT CLUSTER\n");
}


