#ifndef ANTIPOLE_H
#define ANTIPOLE_H 

#define DIST_FUNC double (*dist)( struct ap_Point *p1, struct ap_Point *p2 )

struct ap_Point {
   int id;                       /* point id */
   void *vec;                    /* position vector */
   //double dist_centroid;       /* distance to centroid of cluster */
   struct ap_List *ancestors;    /* list of all ancestors in tree */
};

struct ap_List {
   struct ap_Point *p;           /* point in list */
   double dist;                  /* distance to ancestor, centroid, or query */
   struct ap_List *next;         /* pointer to next member of list */
};

struct ap_Cluster {
   struct ap_Point *centroid;    /* geometric median of cluster */
   //int size;                   /* number of points in cluster */
   double radius;                /* distance from centroid to farthest point in cluster */
   struct ap_List *members;      /* list of points in cluster */
};

struct ap_Tree {
   int is_leaf;                  /* can be a leaf or an internal node */
   struct ap_Point *a, *b;       /* if internal node, pointers to antipoles */
   double radius_a, radius_b;    /* if internal node, distances from antipoles to their farthest point in cluster */
   struct ap_Tree *left, *right; /* if internal node, left and right branches */
   struct ap_Cluster *cluster;   /* if leaf, pointer to cluster */
};

struct ap_Tree* build_tree( struct ap_List *set, double target_radius, struct ap_List *antipoles, DIST_FUNC );
struct ap_Cluster* make_cluster( struct ap_List *set, DIST_FUNC );

void add_point( struct ap_List **set, struct ap_Point *p, double dist );
void move_point( struct ap_Point *p, struct ap_List **from, struct ap_List **to );
void move_nth_point( int n, struct ap_List **from, struct ap_List **to );
struct ap_List* copy_list( struct ap_List *set );
void free_list( struct ap_List *set );
int set_size( struct ap_List *set );

struct ap_Point* exact_1_median( struct ap_List *set, DIST_FUNC );
struct ap_Point* approx_1_median( struct ap_List *set, DIST_FUNC );

void exact_antipoles( struct ap_List *set, struct ap_Point **ap1, struct ap_Point **ap2, DIST_FUNC );
void approx_antipoles( struct ap_List *set, struct ap_Point **ap1, struct ap_Point **ap2, DIST_FUNC );
struct ap_List* adapted_approx_antipoles( struct ap_List *set, double target_radius );
struct ap_List* check( struct ap_List *set, double target_radius, struct ap_Point *antipole );

#endif /* ANTIPOLE_H */

