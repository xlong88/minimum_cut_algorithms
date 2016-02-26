/* Routines related to the sparse certificate algorithm given by
 * Nagamochi and Ibaraki
 *
 * $Id: ni_sparse_cert.h,v 1.2 1997/04/11 16:17:10 mslevine Exp $ */

#include "graph.h"

#ifndef NI_SPARSE_CERT_H
#define NI_SPARSE_CERT_H

/********************************************* prototypes */

/* This function is really only an odd mutation of a sparse
 * certificate computation. We compute the certificate, but then
 * rather than just deleting extra edges, so that the graph is the
 * certificate, we contract all the extra edges.  For the purposes of
 * finding one minimum cut this action makes sense, because as long as
 * limit >= mincut-value, the extra edges cannot be mincut edges.  We
 * also do a little heuristic hunting for a better upper bound on the
 * minimum cut than limit, and if we find one use it instead of
 * limit. This includes looking at the last arc scanned, which is
 * always contracted. */
node *sparseCertContract(graph *g, weight_t limit);

/* run Matula's (2+epsilon)-approximation algorithm (minCap updated) 
 * non-destructive */
void matulaApprox(graph *g, double epsilon);

#endif
