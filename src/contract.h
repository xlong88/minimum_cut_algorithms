/* contract.h: routines for graph contraction 
* $Id: contract.h,v 1.4 1997/04/15 15:18:16 mslevine Exp $ */

#include "graph.h"

#if defined(CONTRACT) && !defined(CONTRACT_H)
#define CONTRACT_H

/* set auxArc fields of v's neighbors to point to v (used for compaction)
 * This function tries to be lazy by remembering the last vertex it
 * processed. Call with v=NULL to clear memory.  */
void setAuxArcs(graph *g, node *v);

/* delete incident arc a from v in g */
void DeleteArc(graph *g, node *v, arc *a);

/******* set-union contraction functions *******/

/* set-union contraction */
void setUnionContract(graph *g, node *v, node *w);

/* find the node this vertex is part of */
node *findLeader(node *v);

/* compactify set-union representation */
void compact(graph *g);

/******* compact contraction functions *********/

/* compact contraction */
void compactContract(graph *g, node *v, node *w, int type);

/* the type option above says what do do when v and w have a common
 * neighbor.  we either merge v's edge into w's or vice versa. this
 * makes a difference if compact_contract is called while scanning v's
 * edges. type should be one of: */
enum {MERGE_INTO_V,     /* append w's arcs to v's */
      MERGE_INTO_W};     /* prepend w's arcs to v's */

#ifdef K
/* special form of compactContract for use during K's 2-respect
   computations.  does NOT delete self-loops. does NOT update vertex
   capacities or leaders.  does NOT support being called inside a
   ForAllIncidentArcs() loop, as does comapctContract. */
void compactContractLite(graph *g, node *v, node *w);
#endif

/******* other functions ***********************/

/* print current (contracted) graph */
void printGraph(graph *g); 

#ifdef ENABLE_GRAPH_BACKUPS
/* save enough information so that contractions can be undone.
   DOES NOT SAVE EVERYTHING (use copyGraph for that) 
   in particular, will not save all flow information */ 
graph *makeBackupGraph(graph *g);

/* restore g from backup*/
void restoreBackupGraph(graph *g, graph *backup);

/* backups are not well formed, so need a special cleanup */
void freeBackupGraph(graph *backup);
#endif /* ENABLE_GRAPH_BACKUPS */

#endif
