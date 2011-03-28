#include <assert.h>  /* assert */
#include <math.h>    /* fmax */
#include <stdlib.h>  /* NULL */
#include "antipole.h"


// Create an ap_Tree that serves as the root, an internal
// node, or a leaf for the tree data structure. Non-leaves
// contain the identities of two antipole points, a left
// subtree and a right subtree which each contain the
// sub-set of points that is nearest its respective antipole
// point, and the radii of the sub-sets. Leaves contain a
// cluster of points.
struct ap_Tree*
build_tree( struct ap_List *set, double radius, struct ap_List *antipoles, double (*dist)( struct ap_Point *p1, struct ap_Point *p2 ) ) {

   // Allocate space for the new ap_Tree
   struct ap_Tree *new_tree = malloc( sizeof( struct ap_Tree ) );
   assert( new_tree );
      
   // Determine if this tree is an internal node or a leaf
   if( antipoles == NULL ) {
      antipoles = adapted_approx_antipoles( set, radius );
      if( antipoles == NULL ) {
         // ...
         new_tree->is_leaf = 1;
         new_tree->cluster = make_cluster( set, dist );
         return new_tree;
      }
   }

   new_tree->is_leaf = 0;

   // ...
   new_tree->a = antipoles->p;
   new_tree->b = antipoles->next->p;

   // Process the members of set
   struct ap_List *index, *set_a, *set_b;
   double dist_a, dist_b;
   new_tree->radius_a = 0;
   new_tree->radius_b = 0;
   index->next = set;
   while( index->next != NULL ) {
      index = index->next;
      // Determine to which antipole the index is nearest
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

   // Build children
   new_tree->left  = build_tree( set_a, radius, check( set_a, radius, new_tree->a ), dist );
   new_tree->right = build_tree( set_b, radius, check( set_b, radius, new_tree->b ), dist );

   return new_tree;
}


// Create an ap_Cluster owned by a leaf of the tree data
// structure containing a list of the points in the cluster
// (already determined to be sufficiently close to one
// another to group together), the identity of the geometric
// median of the cluster, and the cluster radius.
struct ap_Cluster*
make_cluster( struct ap_List *set, double (*dist)( struct ap_Point *p1, struct ap_Point *p2 ) ) {

   // Allocate space for the new ap_Cluster
   struct ap_Cluster *new_cluster = malloc( sizeof( struct ap_Cluster ) );
   assert( new_cluster );

   // Find the geometric median of the cluster
   new_cluster->centroid = approx_1_median( set );

   // Process the members of the set
   struct ap_List *index;
   double dist_centroid;
   //new_cluster->size = 0;
   new_cluster->radius = 0;
   index->next = set;
   while( index->next != NULL ) {
      index = index->next;
      // For every point in the set (besides the centroid), find
      // the distance to the centroid, add the point to the list
      // of points in the cluster, and update the radius of the
      // cluster if necessary
      if( index->p != new_cluster->centroid ) {
         dist_centroid = dist( new_cluster->centroid, index->p );
         //index->p->dist_centroid = dist_centroid;
         add_point( &(new_cluster->neighbors), index->p, dist_centroid );
         //new_cluster->size++;
         new_cluster->radius = fmax( new_cluster->radius, dist_centroid );
      }
   }

   return new_cluster;
}


// Append to an ap_List an ap_Point with a distance value
// (to an ancestor or cluster centroid, depending on the use
// of the ap_List). Requires the address of an ap_List
// pointer so that if the ap_List pointer is not pointing to
// anything (i.e., the list is currently empty), its value
// can be changed to the memory address for a newly
// allocated ap_List (i.e., the list can be given a first
// member).
void
add_point( struct ap_List **set, struct ap_Point *p, double dist ) {

   // Create the new ap_List to be appended to the linked list
   struct ap_List *new_list_member = malloc( sizeof( struct ap_List ) );
   assert( new_list_member );
   new_list_member->p = p;
   new_list_member->dist = dist;
   new_list_member->next = NULL;

   // Determine whether the linked list is empty
   if( *set == NULL ) {
      // If empty, make the new ap_List pointer the head of the
      // list
      *set = new_list_member;
      return;
   } else {
      // Otherwise, append the new ap_List pointer to the end of
      // the list
      assert( (*set)->p != NULL );  // make sure the user did not call malloc on
                                    // the ap_List pointer (i.e., the first
                                    // member of the list should be initialized)
      struct ap_List *index = *set;
      while( index->next != NULL )
         index = index->next;
      index->next = new_list_member;
      return;
   }
}


// ...
struct ap_List*
adapted_approx_antipoles( struct ap_List *set, double radius ) {
}


// ...
struct ap_Point*
approx_1_median( struct ap_List *set ) {
}


// ...
struct ap_List*
check( struct ap_List *set, double radius, struct ap_Point *antipole ) {
}


