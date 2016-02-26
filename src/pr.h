/* code for the Padberg-Rinaldi heuristics
 *
 * all codes use PRpreprocess() to handle preprocessing
 * NI makes use of PRPass() and sourcePR() for internal tests
 * KS has its own code for internal tests
 * HO makes use of sourcePR12() and sourcePR34()
 * K barely has internal tests; may run PRPass() again
 *
 * need a heap h to be defined and initialized somewhere outside!! 
 *
 * $Id: pr.h,v 1.5 1997/04/13 05:54:56 mslevine Exp $ */

#ifndef NO_PR

#include "graph.h"
#include "contract.h"

#ifndef AGE_INCR
#define AGE_INCR 2   /* increment for magic age field when node skipped */
#endif

#ifndef PR_H
#define PR_H

/* apply 1-4 as much as possible in linear time, iterating if reduce to 
   less than bailout fraction of nodes */
void PRpass(graph *g, double bailout, double bailout12, double bailout34,
	    char countAsPreprocessing);    

void sourcePR(graph *g, node *source, double bailout); 

/* PRpass and print stats */
void PRpreprocess(graph *g, double bailout, double bailout12, 
		  double bailout34); 

#endif /* ! PR_H */


/* during flow computations the arc "capacities" are actually residual
   capcities, so we need to take some special action to get the
   original capacities, which teh PR tests use. This requires some
   special functions. */
#if defined(HO_INTERNAL) && !defined(PR_HO_INTERNAL)
#define PR_HO_INTERNAL
void sourcePR12(graph *g, node *source);
void sourcePR34(graph *g, node *source);
#endif /* HO_INTERNAL */

#endif  /* ! NO_PR */
