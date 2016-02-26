/* Karger's sampling+tree-packing minimum cut algorithm 
 * (Gabow tree packing + O(n) treefixs)
 *
 * Description: Given an undirected graph with edge weights, this
 * algorithm finds a partition of the set of vertices into two
 * non-empty subsets such that the total weight of the edges crossing
 * between the two sets is minimized.
 *
 * Input: A global variable currentN that has the current number of
 * nodes and a global variable nodes that has an array of nodes with
 * edge lists is expected. This not the cleanest way to get input, but
 * it is consistent with other mincut implementations with whcih this
 * shares code, and it is not hard to fix either.
 *
 * Output: The value of the minimum cut. (Total weight of edges
 * crossing partition.) The vertex partition will be marked in the
 * in_cut field of the nodes if compiled with SAVECUT defined.
 *
 * Reference: David R. Karger, "Minimum Cuts in Near Linear Time", in 
 * proceedings of the ACM Symposium on the Theory of Computation 1996.
 *
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 * $Id: karger_mincut.h,v 1.6 1997/04/13 20:55:25 mslevine Exp $ */

#ifndef KARGER_MINCUT_H
#define KARGER_MINCUT_H

#define K

#include "graph.h"

/************************************ Constants ***************/

/* The following defines set the defaults for several tuneable
 * parameters.  If you want to compare different values, it may make
 * sense to edit the Makefile to produce several versions with
 * different setting of these parameters */
#ifdef GABOW
#define CHECK_GABOW_CUT
#endif

/* This algorithm is Monte Carlo, so we only guarantee that the answer
   is right with probability at least SUCCESSPROB. */
#ifndef SUCCESSPROB
#define SUCCESSPROB 0.95
#endif

/* If we get a cut estimate less than <value>*(previous estimate) the
   PR tests are run again. */
#ifndef PRCONT_MIN
#define PRCONT_MIN 0.95
#endif

/* This is a magic parameter. The sampling probability in karger_mincut 
 * is PROBCONST*log(n)/mincut. By the analysis in the paper, this needs 
 * to be 100-150 for high probability of success, but 6 seems to work
 * well... */
#ifndef PROBCONST
#define PROBCONST 150.
#endif

/* maximum density graph to use dynamic trees on */
#ifndef MAXDYNDENSITY
#define MAXDYNDENSITY 0.5
#endif

/* If using integer weights, bail to Gabow's mincut algorithm when
   sampling prob is > GABOWBAILOUT */
#ifndef GABOWBAILOUT
#define GABOWBAILOUT 0.4 
#endif

/* Number of trees to check for 2-respecting cuts. It seems that a
 * substantial fraction of the trees usually 2-respect, so we can just
 * check a few...  This can be a function of k, the number of trees.*/
#ifndef SAMPLE_2_RESP 
#define SAMPLE_2_RESP k 
#endif

/************************************ Prototypes ***************/

weight_t karger_mincut(graph *g, double successprob, double disrespect,
		       double gabowbailout, double maxdyndensity);

#endif
