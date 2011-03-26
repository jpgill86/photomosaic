#ifndef ANTIPOLE_H
#define ANTIPOLE_H 

#include <stdint.h>

#define DIM 3     /* dimensionality of the mean RGB data */

struct point {
   int id;                       /* image id */
   uint8_t vector[DIM];          /* mean RGB data */
   //double dist_centroid;       /* distance to centroid of cluster */
   struct point_list* ancestors; /* list of all ancestors in tree */
};

struct point_list {
   struct point* p;              /* point in list */
   double dist;                  /* distance to ancestor or centroid */
   struct point_list* next;      /* pointer to next member of list */
};

struct cluster {
   struct point* centroid;       /* geometric median of cluster */
   double radius;                /* distance from centroid to farthest point in cluster */
   int size;                     /* number of points in cluster */
   struct point_list* list;      /* list of points in cluster */
};

struct tree {
   int is_leaf;                  /* can be a leaf or an internal node */
   struct point* a, b;           /* if internal node, pointers to antipoles */
   double radius_a, radius_b;    /* if internal node, distances from antipoles to their farthest point in cluster */
   struct tree* left;            /* if internal node, left branch */
   struct tree* right;           /* if internal node, right branch */
   struct cluster* c;            /* if leaf, pointer to cluster of points */
};

double dist( struct point p1, struct point p2 );

#endif /* ANTIPOLE_H */

