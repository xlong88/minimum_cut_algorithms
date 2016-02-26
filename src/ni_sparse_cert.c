/* NI sparse certificate routines 
 * $Id: ni_sparse_cert.c,v 1.6 1997/05/16 03:41:54 mslevine Exp $ */

#include <stdlib.h>
#include <stdio.h>
#ifndef __CHECKER__
#include <assert.h>
#endif

#include "ni_sparse_cert.h"
#include "heap.h"
#include "timer.h"
#include "pr.h"
#include "contract.h"

/************************************************ shared macros *******/

#define InTree -1

/************************************* global variables **************/
int numScans=0;   /* operation counts */

extern heap   h;

/************************************* prototypes **************/
/* save cuts discovered by the alpha heuristic */
static void saveAlphaCut (graph *g, weight_t alphaP);

/************************************* definitions **************/

/* useful for cuts updated by alpha test */
static void saveAlphaCut (graph *g, weight_t alphaP )
{
#ifdef SAVECUT
  node *w;
#endif

  g->minCap = alphaP;
  g->dtime = timer ();
  g->cutCount++;

#ifdef SAVECUT
  ForAllVertices(g, w)
    if ( findLeader ( w ) -> key == InTree )
      w -> in_cut = 1;
    else
      w -> in_cut = 0;
#endif
}

/* NI sparse certificate computation (scan-first search) */
node *sparseCertContract(graph *g, weight_t limit)
{
  arc *a;
  node *w, *v, *last_node=NULL, *vl=g->nodes[0], *wl;
  weight_t newKey;
  weight_t alphaP = 0;
  int phaseScans = 0;
  int phaseContracts = 0;
  int saveN = g->currentN;

  hInsert(h, g->nodes[0]);

  while ( nonEmptyH(h) && (++phaseScans < saveN))
    {
      extractMax ( h, v );
      v -> key = InTree;
      
      /* NOI alpha heuristic */
      
      alphaP += v -> cap -  2.0 * v -> key;
      if ( alphaP < g->minCap && nonEmptyH( h )) 
	{
	  saveAlphaCut(g, alphaP);
	  if (g->minCap < limit) limit = g->minCap;
	}

      /* scan v */      
      vl = findLeader(v);
      ForAllIncidentArcs ( v, a ) 
	{
	  g->edgeScanCnt[g->context]++;
	  w = a->head;
	  if ( w->key > InTree ) {
	    newKey = w->key + a->cap;
	    if (( newKey >= limit ) && (vl != (wl=findLeader(w))))
	      {
		setUnionContract(g, vl, wl);
		last_node = vl;
		phaseContracts++;
	      }
	    
	    if ( w->key == 0 ) {
	      w->key = newKey;
	      hInsert ( h, w );
	    }
	    else {
	      w->key = newKey;
	      increaseKey ( h, w, newKey );
	    }
	  }
	}
    } 
  
  g->numScans += phaseScans;
  if ( phaseScans != saveN && g->minCap > 0) 
    {
      saveAlphaCut(g, 0);
      return NULL;
    }
    
  /* contract the last arc scanned if not contracted */
  assert(h.size == 1);
  extractMax(h, w);
  if (vl != (wl=findLeader(w)))
    {
      setUnionContract(g, vl, wl);
      last_node = vl;
      phaseContracts++;
    }

#ifdef VERBOSE
  printf("Phase: %d nodes contracted (%d left)\n", 
	 phaseContracts, currentN );
#endif

  return last_node;
}

void matulaApprox(graph *g, double epsilon)
{
  node *v;

  epsilon += 2;
  
  /* main loop */
  if (g->currentN > 2)  
    do{
      ForAllNodes(g, v) 
	{
	  v->key = 0;  /* remove junk left behind by PR */
#ifdef PR_34
	  v->age = 0;  /* clear contract counts */
#endif      
	}

      g->context = MAIN;
      v = sparseCertContract(g, g->minCap/epsilon);
      if (g->currentN <= 1) break;
      compact(g);

#ifndef NO_PR
#if PR_FREQ > 1
      if (!(numPhases % PR_FREQ))
#endif
	{
	  sourcePR(g, v, 0.9);
	  PRpass(g, 0.75, 0.5, 0.5, 0);
	}
#endif
    } while ( g->currentN > 2 );
}
