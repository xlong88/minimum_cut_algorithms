/*----------------------------------------------------------------------
$Id: ks.h,v 1.4 1997/04/11 16:17:09 mslevine Exp $

TYPES for ks.c

Variable naming convention may seem quite odd.  it's a scheme called
hungarian; full text available on request, but to summarize what's
used here: 
base typenames:
  integer i, double dbl, count c, function fn, boolean f
  metavertex mv, metaedge me, metagraph mg
  weight wt
prefixes:
  p pointer, h handle (=pp), rg array
sufixes:
  Max: upper limit on value.  not a legal value
  Out: output 
  Tmp: temporary storage
Functions get capitalized.
----------------------------------------------------------------------*/

#include "graph.h"  /* for weight_t and other such fun things */

typedef int bool;
#define fFalse 0
#define fTrue 1

/*-------------------------------------------------------------------
struct metavertex

It was necessary for metavertices to support the disjoint set union
operations; thus the c field to store ranks.  similarly the pmv field,
which is not a next pointer but rather a parent pointer for the set
structure.  It may appear that the use of this disjoint set-union
makes pgphContract destructive since it changes parent pointers in the
vertices, but the fact to note is that one can back out of the changes
caused by a given call by clearing the parent pointers of all the
items which were roots at the time of the call.  Since the graph
structure contains precisely such a list, this is easy.

WARNING: this will not work if you use path compression to construct 
the list of original nodes on one side of the minimum cut.  This is 
the reason for the "safe" non-path-compressing version of FIND.

The i field is used in countingSort.

----------------------------------------------------------------------*/

typedef struct metavertex 
    {
    struct metavertex *pmv;  /*parent for set union*/
    int name;                /* name of vertex (number) */
    int c;                   /*rank for set union*/
    int i;                   /*index for counting sort of edges*/
    weight_t wt;               /*for accumulting degree*/
#ifdef ALLCUTS
    int iKey;			/* Key for identifying identical cuts */
#endif
#if (PR2==TIGHT)
    int iTime;                /*the last phase during which this
				vertex was responsible for a tight PR2
				test causing a contraction.  Used to
				ensure only one such contraction per
				vertex per phase.*/
#endif
#ifdef PR34
      struct metavertex *pmv34;
      weight_t wt34;
      weight_t wt34Tot;
      struct metaedge *pmeFirst;
      struct metaedge *pmeMax;
      bool fScanned;
#endif
    } metavertex;
/*type label mv*/

/*----------------------------------------------------------------------
struct metaedge

Fields mainly obvious.  The unexpected i field holds a random ineger
for each edge.  This allows each small cut we costruct to be
identified by a (probably) unique integer, namely the XOR of all the
tags of edges in the cut.  This in turn allows fast determination of
whether a cut has been see before (see fConditionalCallback).
pmgCompact() maintains the tags in the contracted edges.

----------------------------------------------------------------------*/

typedef struct metaedge {
    metavertex *rgpmv[2];
    weight_t wt;			/* weight */
    } metaedge;
/* type label me*/

typedef struct metagraph 
    {
    metaedge *rgme;
    metavertex **rgpmv;     
    int cme;
    int cmv;
    } metagraph;
/* type label mg*/

