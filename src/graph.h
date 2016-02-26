/* graph type definitions and general purpose functions
 *
 * common graph types for HO, NI, K (KS also uses this at startup) the
 * data structures here are fundamentally the same, but little
 * differences in the algorithms make for little changes, and the
 * result is that this file is a mess of conditional compilation
 * directives. The good thing is that most such dependencies are
 * handled here, so not everything looks this bad.
 * $Id: graph.h,v 1.5 1997/05/16 03:41:52 mslevine Exp $ */

#ifndef GRAPH_H
#define GRAPH_H

#include "memassert.h"
#include <limits.h>

 /* Pick weight type. The options are DOUBLEWEIGHTS, INTWEIGHTS, and 
  * LONGLONGINTWEIGHTS (requires GNU C compiler).
  * This gives:
  * MAXWEIGHT = max number of type weight_t
  * WT_RD_FORMAT = scanf string for reading weight_t
  * fprintf_wt(fp, format, value) = function to print one value of type
  *    weight_t. The "format" string here can only be for placement purposes, 
  *    eg, fprintf_wt(stdout, "%20", minCap) to print minCap right justified 
  *    in a field of 20 spaces.
  *
  * We default to long long int if we're using gcc and long int otherwise.
  */    
#ifdef DOUBLEWEIGHTS
#include <float.h>
typedef double weight_t;
#define fprintf_wt(fp, format, val) fprintf(fp, format ".4f", val)
#define MAXWEIGHT MAXDOUBLE
#define WT_RD_FORMAT "lf"
#else /* DOUBLEWEIGHTS */
#define INTWEIGHTS
#define WT_RD_FORMAT "ld"
#if defined(__GNUC__) && !defined(__CHECKER__)
#define LONGLONGINTWEIGHTS
#include "fprintfll.h"
/* I shouldn't do this, but adams refuses to acknowledge that it has this 
   definition in its own header files. */
#ifndef LONG_LONG_MAX
#define LONG_LONG_MAX 9223372036854775807LL
#endif 
typedef long long int weight_t;
#define MAXWEIGHT LONG_LONG_MAX
#define fprintf_wt(fp, format, val) fprintfll(fp, format, val)
#else /* GNUC */
typedef long int weight_t;
#define MAXWEIGHT LONG_MAX
#define fprintf_wt(fp, format, val) fprintf(fp, format "ld", val)
#endif /* GNUC */
#endif /* DOUBLEWEIGHTS */

/* type of key in heaps. K with PST needs double, NI needs weight_t, rest don't
 * really care */
#if defined(K) && defined(PST)
typedef double heap_weight;
#else
typedef weight_t heap_weight;
#endif

/* resolve conflicts in options */
#ifdef NO_PR
#undef PR_34      /* NO_PR means !PR_34 too*/
#else
#if !defined(NO_PR34) && !defined(PR_34)
#define PR_34     /* !NO_PR but no indication of PR34 or not defaults to 34 */
#endif
#endif /* NO_PR */

#if defined(HO_INTERNAL) && !defined(HO)
#define HO
#endif

/* HO uses a different name for the capacity field */
#ifdef HO
#define CAPFIELD resCap     
#else
#define CAPFIELD cap
#endif

#ifdef K
#define MATULA
#endif

#ifdef MATULA
#define NI
#define ENABLE_GRAPH_BACKUPS
#endif

/* CONTRACT is defined when we intend to do contractions */
#if !defined(NO_PR) || defined(NI) || defined(KS) || defined(K)
#define CONTRACT
#endif


/************************************************ data types *******/

/*************** arc *****************/
typedef struct arcSt
{
  weight_t        CAPFIELD;          /* capacity */
  struct nodeSt   *head;           /* head node */
  struct arcSt    *rev;            /* reverse arc */
  struct arcSt    *next;           /* next in the arc list */
#ifdef CONTRACT 
  struct arcSt    *prev;           /* prev in the arc list */
#ifdef ENABLE_GRAPH_BACKUPS
  int              index;          /* index in currentArcs */
#endif
#endif
#ifdef K
  struct nodeSt   *lca;            /* LCA in tree */
#ifdef PST
  double           usage;
  double           cost;
#endif /* PST */
#endif /* K */
}
arc;

/*************** node *****************/
typedef struct nodeSt
{
  arc             *first;           /* first outgoing arc */
  arc             *last;            /* last outgoing arc */

#ifdef HO
  arc             *current;         /* current incident arc */
  weight_t         excess;          /* excess of the node */
  int              d;               /* distance from the sink */
  struct nodeSt   *bNext;           /* next node in bucket */
  struct nodeSt   *bPrev;           /* previous node in bucket */
  struct nodeSt   *nextA;           /* next active node in the bucket */
  unsigned int     status:2;        /* 1 regular, 0 frozen, 2 source */
  unsigned int     toContract:1;
#endif  /* HO */

#ifdef CONTRACT
#ifdef PR_34
  unsigned int     pr34scanned:1;   /* to keep linear time on 34 tests */
  unsigned int     age:28;          /* contract counts */
#endif /* PR_34 */
  struct nodeSt   *leader;          /* contracted component leader */
  weight_t         cap;             /* total adjacent arc capacity */
  arc             *auxArc;          /* for detecting parallel arcs */
  heap_weight      key;             /* heap key */
  int              heap_pos;        /* heap pos */
  int              index;           /* index as node (in nodes) */
#endif /* CONTRACT */

#ifdef K
  /* for k only */
#ifndef CONTRACT    /* k always needs cap */
  weight_t         cap;             /* weighted degree of vertex */
#endif
  weight_t         capdown;         /* sum of wght'd degrees of descendants */
  weight_t         rho;             /* sum of edges with this lca */
  weight_t         rhodown;         /* sum of edges crossing below */
  weight_t         cutdown;         /* cut value w/ descendants as one side */
  struct nodeSt   *parent;          /* union-find parent for LCA */
  void            *back;            /* pointer to dynamic tree node 
				       (see dyn_tree/dyn_tree.h) */
#ifdef GABOW
  arc             *used_first;      /* first in list of used edges */ 
  struct nodeSt   *next_ftree;      /* next active ftree */
  struct nodeSt   *prev_ftree;      /* prev active ftree */
  struct nodeSt   *ftree_label;     /* index of root of ftree vertex is in */
  unsigned int     has_no_labeled_edges:1;/* true iff !unused edges w/labels */
  unsigned int     active:1;        /* true iff active ftree root */
#endif /* GABOW */
#if defined(PST) && !defined(CONTRACT)
  arc *auxArc;
  heap_weight key;
  int heap_pos;
#endif
#endif  /* K */

#ifdef SAVECUT
  unsigned int     in_cut:1;        /* to save the cut */
#endif
} node;

#ifdef HO
typedef /* bucket */
 struct bucketSt
 {
   node             *first;          /* first regular node */
   node             *firstActive;    /* first active node */
 } bucket;

#include"qs.h"
#endif

/**************** graph ******************/

#define CONTEXTS 7
enum contexts {PRPRE12, PRPRE34, PR12, PR34, MAIN, GLUPDATE, PACKING}; 

typedef struct graphSt
{
  int      n,                  /* number of vertices */ 
           m,                  /* number of arcs */
           currentN,           /* current number of nodes */
           currentM;           /* curent number of arcs */

  node     *vertexalloc;       /* memory allocation for vertices; not to 
				  be used by anything but dimacsParse 
				  and freeGraph */
           
  node     *vertices,          /* array of vertices */
           *sentinelVertex;    /* next after last */
  arc      *arcs,              /* array of arcs */
           *sentinelArc;      /* next after last */

#ifdef CONTRACT
  node    **nodes;             /* current nodes (identified by some vertex) */
#ifdef ENABLE_GRAPH_BACKUPS
  arc     **currentArcs;       /* current arcs */
#endif
#endif

  weight_t  minCap;            /* minimum cut capacity seen */

  /* statistics */
  int       cutCount;          /* count of new cuts */
  float     dtime;             /* cut discovery time */
  int       PR1Cnt, PR2Cnt, PR3Cnt, PR4Cnt; /* counts of PR contractions */

  int       context;           /* decides who to charge edge scans to */
  long      edgeScanCnt[CONTEXTS];   /* number of edges scanned (by context)*/

#ifdef NI
  int       numScans;          /* node scans */
#endif

  /* HO needs much extra stuff */
#ifdef HO
  int   regN;                 /* number of regular nodes */
  int   sourceD;              /* source label */
  int   sinkD;                /* sink label */
  int   infD;                 /* label "infinity" */
  int   dMax;                 /* maximum label */
  int   aMax;                 /* maximum actie node label */
  bucket *buckets;             /* array of buckets */
  node   *source;              /* new source node */
  node   *sink;                /* sink node */
  
  queue  *BFSqueue;           /* queue for BFS */
  
  stack  *freezer;             /* stack of node layers */
  stack  *regNodes;            /* stack of regular nodes in globalUpdate */
  stack  *exStack;             /* stack of nodes to be excess contracted */
  
  /* more stats */

  int STCutCnt;                /* number of max-flow computations */
  int oneCnt;                  /* number of one node layer cntracts */
  int pushCnt;                 /* number of pushes */
  int relabelCnt;              /* number of relabels */
  
  int updateCnt;               /* number of global updates */
  int excessCnt;               /* number of contractions due to extra excess */
  int gapCnt;                  /* number of gaps */
  int gNodeCnt;                /* number of nodes after gap */  
  int relsSinceUpdate;         /* the number of relabels since last update */
  int totalSize;               /* total size of all s-t problems solved */
#endif 

  /* K needs much extra stuff too */
#ifdef K 
  char *is_descendant;    /* nxn array of descendant relation for a tree 
			     computed by dfs1_resp_desc, used in 
			     dfs2_resp[AB] */
  weight_t *subtreecuts;  /* cut values between all pairs of subtrees 
			     computed by dfs2_resp[AB] */
  struct tree *currentTree;  /* current tree for dfs computations */
  node *root;             /* root of trees */
  
#ifdef SAVECUT
  node *save1, *save2;    /* vertices that define cut in */
  struct tree *savet;            /* tree */
  char saveDescendant;    /* are the vertices comparable or not ? */
#endif /* SAVECUT */
#endif /* K */

}
graph;

/*********************************************** variables *****/

#ifdef CONTRACT
node **__tmp_node;     /* for use in ForAllNodes */
arc **__tmp_arc;       /* for use in ForAllArcs */
#endif

/************************************************ macros *******/

#define MAX( x, y ) ( ( (x) > (y) ) ?  x : y )
#define MIN( x, y ) ( ( (x) < (y) ) ? x : y )
#define ABS( x ) ( (x) >= 0 ) ? (x) : -(x)

#define Reverse( a ) ((a)->rev)                      /* reverse arc */
#define arcCap( a ) ((a)->CAPFIELD)                  /* arc capacity */
#define VERTEX_NAME(g, i) (1+(int)((i)-g->vertices))   /* vertex's name */

#define ForAllVertices(g, v) for (v=g->vertices; v!=g->sentinelVertex; v++)
#define ForAllIncidentArcs(v,a) for ( a = v->first; a != NULL; a = a->next )

#ifdef CONTRACT
#define ForAllNodes(g, v) \
  for ( __tmp_node=g->nodes; (v=*__tmp_node) != NULL; __tmp_node++)
#define ForAllArcs(g, a) \
  for ( __tmp_arc=g->currentArcs; (a=*__tmp_arc) != NULL; __tmp_arc++)
#define NODE(g, i) (g->nodes[(i)])   /* the ith node */
#define INDEX(i) ((i)->index)         /* node's index */
#define isLeader(v) ((v)==(v)->leader)
#else  /* CONTRACT */
#define ForAllNodes(g, v) for (v=g->vertices; v!=g->sentinelVertex; v++)
#define ForAllArcs(g, a)  for (a=g->arcs; a!=g->sentinelArc; a++)
#define NODE(g, i) (g->vertices+(i))
#define isLeader(v) (1)
#endif /* CONTRACT */

/************************************************ prototypes *******/

/* make a new graph */
graph *makeGraph(int n, int m);

/* return sum of all edge scans */
long totalEdgeScans(graph *g);

/* print current (contracted) graph */
void printGraph(graph *g);    

#ifdef SAVECUT
/* for every node, print which side of the cut it is on */
void printCut(graph *g);
#endif

/* save the cut that separates v and the rest of the graph */
/* we pass val in case there is no cap field */
void saveTrivialCut(graph *g, node *v, weight_t val);

/* clean up */
void freeGraph(graph *g);

/* read in a graph in DIMACS format */
graph *dimacsParse(FILE *fp);

#ifndef CONTRACT
/* this is just to clean up the input when we don't have contraction code */
void compact(graph *g);
#endif

#endif /* ! GRAPH_H */

