/* 
 * Gabow's round robin tree packing algorithm.
 * 
 * Description: Given a directed (multi)graph and a given vertex A,
 * this algorithm finds a maximal set of edge-disjoint spanning trees
 * subject to the constraint that A has in-degree 0, and all the other
 * vertices have the same in-degree as the number of trees. (The
 * in-degree constraints guarantee that this tree packing can be
 * modified to an A-arborescence, which I haven't done (See the
 * paper))
 *
 * Input: The graph is specified by a pointer to the node A, the
 * number of vertices, N, an array of length N, Nodes, that has the
 * node structures with edge_lists.
 *
 * Output: The resulting trees.
 *
 * Known Bugs: The final (failing) attempt to add another tree can
 * screw up the degree constraints. (It doesn't screw up the trees.)
 * This can easily be fixed by keeping a copy of the tree labels after
 * each successful tree addition and reverting to the last when adding
 * another tree fails; however, since the application I wrote this for
 * doesn't care about the degree constraints, and fixing the problem
 * requires a definite amount of added time and space, I decided to
 * leave the bug. 
 *
 * Reference: Harold N. Gabow, "A Matroid Approach to Finding Edge
 * Connectivity and Packing Arborescences", Journal of Computer and
 * System Sciences 50, 259-273 (1995)
 *
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 * $Id: gabow_pack.h,v 1.8 1997/05/16 15:30:31 mslevine Exp $ */

#ifndef GABOW_PACK_H
#define GABOW_PACK_H

#define K
#ifndef GABOW
#define GABOW
#endif
#include "graph.h"

#define first(t, v) ((t)->info[INDEX(v)].first)
#define savefirst(t, v) ((t)->info[INDEX(v)].savefirst)
#define parent_edge(t, v) ((t)->info[INDEX(v)].parent_edge)
#define ForAllTrees(p, t)  \
  for((t)=(p)->first_tree; (t)!=(p)->last_tree; (t)=(t)->next_tree)

#define ForAllTreeEdges(t, v, e) \
  for((e)=(t)->info[INDEX(v)].first; (e)!=NULL; (e)=(e)->next)
#define IsParentEdge(t, v, e) ((e) == ((t)->info[INDEX(v)].parent_edge))
#define HEAD(e) ((e)->head)

/* these two assume parent edges have been deleted */
#define HasChildren(t, v) ((t)->info[INDEX(v)].first != NULL)
#define HasOneChild(t, v) ((t)->info[INDEX(v)].first->next == NULL)

#define DeleteTreeEdge(t, v, e) \
{ \
  if ((e)->prev != NULL)  (e)->prev->next = (e)->next; \
  else (t)->info[INDEX(v)].first = (e)->next; \
  if ((e)->next != NULL)  (e)->next->prev = (e)->prev; \
}


/* one weighted input edge can become many (unweighted) tree edges, so
   we need another data structure */
typedef struct tedge {
  arc *e;                    /* the arc this derived from */
  node *head;                 
  struct tedge *rev;
  struct tree *tree;         /* spanning tree this edge belongs to */
  struct tedge *label;       /* labelling algorithm label */

  struct tedge *next;        /* next,prev pointers for doubly linked edge */
  struct tedge *prev;        /*   lists of trees */

  struct tedge *auxnext;     /* next for some extra lists */

  unsigned int startlabel:1; /* true iff NULL label means start */
} tedge;

/* the packing is basically an array of trees, where a tree is an array 
   indexed over the vertices of "treeinfo */

struct treeinfo
{
  int              depth;           /* depth of vertex in tree */
  struct tedge    *parent_edge;     /* edge to parent in tree */
  struct tedge    *first;           /* list of all incident edges in tree */
#ifdef SAVECUT
  struct tedge    *savefirst;       /* backup for tree destruction that
				       takes place when we use dyn trees */
#endif
  unsigned int     is_labeled:1;    /* true iff vertex labeled in tree */
};

typedef struct tree
{
  int name;                       /* number of this tree */
  struct treeinfo *info;          /* array of info for vertices */
  node *root_of_labeled_subtree;  /* index of vertex that is ... */
  unsigned int was_changed:1;     /* true if tree touched in last updates*/
  struct tree *next_tree;
} tree;

typedef struct packing
{
  int k;                   /* the number of trees */
  tree *first_tree;        /* first tree */
  tree *last_tree;         /* last tree */
  node *onesideroot;       /* root of subtree that defines one side of cut */
  void *tedge_allocator;   /* internal */
  tree *spare_trees;       /* internal */
} *packing;

packing gabow_pack(graph *g, int c, packing p);
void free_packing(packing p);
#ifndef NDEBUG
void printtree(node *v, tree *t, node *parent);
#endif

#endif
