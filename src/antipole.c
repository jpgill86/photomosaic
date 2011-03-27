#include <math.h>    /* sqrt, pow */
#include <stdlib.h>  /* NULL */
#include "antipole.h"


// Calculate the Euclidian distance between two points
double
dist( struct ap_Point *p1, struct ap_Point *p2 ) {
   int i;
   double sum = 0;
   for( i = 0; i < DIM; i++ )
      sum += pow( p1->vector[i] - p2->vector[i] , 2);
   return sqrt( sum );
}


// ...
void
add_point( struct ap_List *set, struct ap_Point *p, double dist ) {
}


// ...
struct ap_List*
check( struct ap_List *set, double radius, struct ap_Point *antipole ) {
}


// ...
struct ap_Cluster*
make_cluster( struct ap_List *set ) {
}


// ...
struct ap_List*
adapted_approx_antipoles( struct ap_List *set, double radius ) {
}


// ...
struct ap_Tree*
build_tree( struct ap_List *set, double radius, struct ap_List *antipoles ) {
   // Tree to be returned
   struct ap_Tree *new_tree = malloc( sizeof( struct ap_Tree ) );
      
   // Determine if this tree is an internal node or a leaf
   if( antipoles == NULL ) {
      antipoles = adapted_approx_antipoles( set, radius );
      if( antipoles == NULL ) {
         // ...
         new_tree->is_leaf = 1;
         new_tree->cluster = make_cluster( set );
         return new_tree;
      }
   }

   new_tree->is_leaf = 0;

   // ...
   new_tree->a = antipoles->p;
   new_tree->b = antipoles->next->p;

   // Process the members of set
   struct ap_List *index;
   struct ap_List *set_a = malloc( sizeof( struct ap_List ) );
   struct ap_List *set_b = malloc( sizeof( struct ap_List ) );
   new_tree->radius_a = 0;
   new_tree->radius_b = 0;
   index->next = set;
   while( index->next != NULL ) {
      index = index->next;

      // Determine to which antipole the index is nearest
      double dist_a, dist_b;
      dist_a = dist( new_tree->a, index->p );
      dist_b = dist( new_tree->b, index->p );
      add_point( index->p->ancestors, new_tree->a, dist_a );
      add_point( index->p->ancestors, new_tree->b, dist_b );
      if( dist_a < dist_b ) {
         add_point( set_a, index->p, dist_a );
         new_tree->radius_a = fmax( dist_a, new_tree->radius_a );
      } else {
         add_point( set_b, index->p, dist_b );
         new_tree->radius_b = fmax( dist_b, new_tree->radius_b );
      }
   }

   // Build children
   new_tree->left  = build_tree( set_a, radius, check( set_a, radius, new_tree->a ) );
   new_tree->right = build_tree( set_b, radius, check( set_b, radius, new_tree->b ) );

   return new_tree;
}

