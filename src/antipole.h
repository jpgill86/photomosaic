#ifndef ANTIPOLE_H
#define ANTIPOLE_H 

#define DIST_FUNC double (*dist)( ap_Point *p1, ap_Point *p2 )

typedef struct ap_Point ap_Point;
typedef struct ap_List ap_List;
typedef struct ap_Cluster ap_Cluster;
typedef struct ap_Tree ap_Tree;

struct ap_Point {
   int id;                    /* point id */
   void *vec;                 /* position vector */
   //double dist_centroid;    /* distance to centroid of cluster */
   ap_List *ancestors;        /* list of all ancestors in tree */
};

struct ap_List {
   ap_Point *p;               /* point in list */
   double dist;               /* distance to ancestor, centroid, or query */
   ap_List *next;             /* pointer to next member of list */
};

struct ap_Cluster {
   ap_Point *centroid;        /* geometric median of cluster */
   //int size;                /* number of points in cluster */
   double radius;             /* distance from centroid to farthest point in cluster */
   ap_List *members;          /* list of points in cluster */
};

struct ap_Tree {
   int is_leaf;               /* can be a leaf or an internal node */
   ap_Point *a, *b;           /* if internal node, pointers to antipoles */
   double radius_a, radius_b; /* if internal node, distances from antipoles to their farthest point in cluster */
   ap_Tree *left, *right;     /* if internal node, left and right branches */
   ap_Cluster *cluster;       /* if leaf, pointer to cluster */
};

ap_Tree* build_tree( int level, ap_List *set, double target_radius, ap_Point *antipole_a, ap_Point *antipole_b, int dimensionality, DIST_FUNC );
ap_Cluster* make_cluster( ap_List *set, int dimensionality, DIST_FUNC );

void add_point( ap_List **set, ap_Point *p, double dist );
void move_point( ap_Point *p, ap_List **from, ap_List **to );
void move_nth_point( int n, ap_List **from, ap_List **to );
ap_List* copy_list( ap_List *set );
void free_list( ap_List *set );
int list_size( ap_List *set );

void exact_1_median( ap_List *set, ap_Point **median, DIST_FUNC );
void approx_1_median( ap_List *set, ap_Point **median, int dimensionality, DIST_FUNC );

void exact_antipoles( ap_List *set, ap_Point **antipole_a, ap_Point **antipole_b, DIST_FUNC );
void approx_antipoles( ap_List *set, ap_Point **antipole_a, ap_Point **antipole_b, int dimensionality, DIST_FUNC );
void adapted_approx_antipoles( ap_List *set, ap_Point **antipole_a, ap_Point **antipole_b, double target_radius, DIST_FUNC );
void check_for_antipoles( ap_List *set, double target_radius, ap_Point *ancestor, ap_Point **antipole_a, ap_Point **antipole_b );

#endif /* ANTIPOLE_H */

