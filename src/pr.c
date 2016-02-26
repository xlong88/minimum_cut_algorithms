/* common code for the Padberg-Rinaldi tests */
/* $Id: pr.c,v 1.6 1997/04/13 05:54:56 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#ifdef __CHECKER__
#include "memassert.h"
#else
#include <assert.h>
#endif

#include "pr.h"
#include "timer.h"
#include "contract.h"

#ifdef PR_34
#include "heap.h"
extern heap h;
#endif

#ifdef HO
#include "ho.h"
#if !defined(NO_EXCESS_DETECTION) && !defined(NO_PR)
#define EXCESS_DETECTION
#endif
#endif

#ifdef HO_INTERNAL
#define arcCap2( a ) (((a -> resCap) + ( Reverse( a ) -> resCap )))
#define PRTest1(a) ((arcCap2(a) >= 2*g->minCap) ? 1 : 0)
#define PRTest2(a, vC, wC) ((arcCap2(a) >= MIN(vC, wC)) ? 1 : 0)
#else
#define PRTest1(a) ((arcCap(a) >= g->minCap) ? 1 : 0)
#define PRTest2(a, vC, wC) ((2*arcCap(a) >= MIN(vC, wC)) ? 1 : 0)

/* apply 1,2 eveywhere */
static void PR12Pass(graph *g);
#ifdef PR_34
/* apply 3,4 out of node */
static void PRTest34(graph *g, node *v); 
/* apply 3,4 as much as possible in linear time */
static void PR34Pass(graph *g);    
#endif /* PR_34 */

#endif /* HO_INTERNAL */

/* ================================================================ 

   Node status for tests 3&4:
   Unseen nodes have pr34scanned == 0

   NO_PR cannot be defined if PR34 is defined 

 ================================================================ */

#ifdef PR_34

/* do PR tests 3 and 4 out of v */
#ifdef HO_INTERNAL
void PRTest34_hoi(graph *g, node *v)
#else
void PRTest34(graph *g, node *v)
#endif
{
  node *w, *u;
  arc *vw, *wu, *a;
  weight_t vCap, sum4;
  int wHere;

  assert(v->leader == v);

  if (g->currentN <= 2) return;

  /* 
     mark neighbours of the v for tests 3,4.
     For the latter, if a neighbour is already marked,
     do not touch it or time may not be linear
  */

  setAuxArcs(g, v);
#ifndef PR_ONLY 
  v->pr34scanned = 1;
#endif

  ForAllIncidentArcs(v,vw) 
    {
      g->edgeScanCnt[g->context]++;
      w = vw->head;
      wHere = 1;

#ifdef HO_INTERNAL
      /* some nodes may get contracted by excess detection */
#ifdef EXCESS_DETECTION
      if (w == NULL || w->toContract) continue;
#endif
      if (((w->status != REGULAR) && (v->status != SOURCE)) ||
	   (w->status == SOURCE) || (w == g->sink))
	continue;
      sum4 = arcCap2(vw);
#else
      sum4 = arcCap(vw);
#endif  /* HO_INTERNAL */
      
      vCap = v->cap;

      if (!w->pr34scanned) 
	{
#ifndef PR_ONLY
	  w->pr34scanned = 1;
#endif
	  ForAllIncidentArcs(w,wu) 
	    {
	      u = wu->head;
    /* some nodes may get contracted by excess detection */
#ifdef EXCESS_DETECTION
	      if (u == NULL) continue;
#endif
	      g->edgeScanCnt[g->context]++;
	      a = u->auxArc;
	      
	      if ((a != NULL) && (a->head == v)) 
		{
		  /* accounting for test4*/
#ifdef HO_INTERNAL
		  sum4 += MIN(arcCap2(wu),arcCap2(a)); 
#else
		  sum4 += MIN(arcCap(wu), arcCap(a)); 
#endif

		  /* test 3 */
#ifdef HO_INTERNAL
		  if ((vCap <= (arcCap2(vw) + arcCap2(a))) &&
		      (w->cap <= (arcCap2(vw) + arcCap2(wu)))) 
#else
		  if ((vCap <= 2*(arcCap(vw) + arcCap(a))) &&
		      (w->cap <= 2*(arcCap(vw) + arcCap(wu)))) 
#endif
		    {
		      compactContract(g, v, w, MERGE_INTO_V);
		      wHere = 0; /* w gone */
		      vCap = v->cap; 
		      g->PR3Cnt++;
		      if (g->currentN <= 2) 
			return ;
		      
		      break;
		    }  
		}
	    }

#ifdef HO_INTERNAL
	  if (wHere && (sum4 >= 2*g->minCap)) 
#else
	  if (wHere && (sum4 >= g->minCap)) 
#endif
	    {
	      compactContract(g, v, w, MERGE_INTO_V);
	      vCap = v->cap; 
	      g->PR4Cnt++; 
	      if (g->currentN <= 2) 
		return;
	    }
	}
    }
}

#endif /* PR_34 */

#ifdef HO_INTERNAL
/* ================================================================ */
/* do PR tests 1,2 near source */  
void sourcePR12(graph *g, node *source)
{
  arc *a;
  node *w;
  weight_t vCap;
#ifdef VERBOSE
  int oldn = g->currentN;
  float tt = timer();
#endif

  g->context = PR12;
  ForAllIncidentArcs(source,a)
    {	  
      if ( g->currentN <= 2 ) break;
      g->edgeScanCnt[g->context]++; 

      w = a->head;
    /* some nodes may get contracted by excess detection */
#ifdef EXCESS_DETECTION
      if (w == NULL) continue;
#endif

      vCap = source->cap;
      assert(w == w->leader); 
      
      if (w != g->sink) { /* OK to contract frozen into source */
	if ((PRTest1(a) && ++g->PR1Cnt) || 
	    (PRTest2(a,vCap,w->cap) && ++g->PR2Cnt)) 
	  {
	    compactContract ( g, source, w, MERGE_INTO_V);  
	    vCap = source->cap;	  
	  }
      }
    }

#ifdef VERBOSE
  printf("c 12 at source did %d PR+excess contracts in %14.2f time \n",
	 oldn -g->currentN, timer()-tt);
#endif
}

#ifdef PR_34
void sourcePR34(graph *g, node *source)
{
  arc *a;
#ifdef VERBOSE
  int oldn = g->currentN;
  float tt = timer();
#endif
  
  g->context = PR34;

  /* unmark neighbours */
  /* should PrWork34 be charged here? I say no because we'll scan for 
     real and charge in the call to PRTest34() */
  ForAllIncidentArcs(source,a)       
    a->head->pr34scanned = 0;

  PRTest34_hoi(g, source);

#ifdef VERBOSE
  printf("c 34 at source did %d PR+excess contracts in %14.2f time \n",
	 oldn -g->currentN, timer()-tt);
#endif

}
#endif /* PR_34 */

#else  /* HO_INTERNAL */

/* ================================================================ */
/* apply PR tests 1 and 2 everywhere */
void PR12Pass (graph *g)
{
  node *v, *w;
  arc *a;
  weight_t vCap;
#ifdef VERBOSE
  int oldcnt = g->PR1Cnt+g->Pr2Cnt;
  float tt = timer();
#endif

  ForAllNodes (g, v) 
    {   
      vCap = v -> cap;
      
      /* store pointers back to v for merging */
      
      ForAllIncidentArcs ( v, a ) 
	{
	  if ( g->currentN <= 2 ) break;
	  w = a->head;  
	  g->edgeScanCnt[g->context]++;
	  assert(w == findLeader ( w ));  
#ifndef PR_ONLY
	  if ( v < w ) /* to insure linear time */
#endif
	    {
	      /* these counters work because C guarantees that a logic is 
		 evaluated left to right, and no more work than necessary 
		 is done */
	      if ((PRTest1(a) && ++g->PR1Cnt) || 
		  (PRTest2(a,vCap,w->cap) && ++g->PR2Cnt)) 
		{
		  compactContract(g, v, w, MERGE_INTO_W);
		  vCap = v->cap;
		}
	    }
	}
    }
  
#ifdef VERBOSE
  printf("c Phase of 12 -- did %d prcontracts in %14.2f time \n",
	 g->PR1Cnt + g->PR2Cnt - oldcnt, timer() - tt);
#endif
}

#ifdef PR_34
  /* apply PR tests 3 and 4 as much as possible in linear time.  the
   * strategy for keeping things linear is to mark nodes when we see
   * them and then not touch nodes that are marked. This is clearly
   * linear, and we can miss tests because an important edge is
   * attached to a node that got markes earlier. We try to avoid
   * getting stuck in a rut of always starting from a useless node by
   * trying the nodes in an order related to how long ago they were
   * last tried and how many nodes they got contarcted into them. */
void PR34Pass(graph *g)
{
#ifdef VERBOSE
  int old3cnt = g->PR3Cnt, old4cnt = g->PR4Cnt;
  float tt = timer();
#endif
  node *v;
  
  clearHeap(h);

  ForAllNodes(g, v) 
    {
      v -> pr34scanned = 0; 
      v->key = v->age;
      hInsert(h,v);
    }
      
  while ( nonEmptyH(h) && g->currentN > 2)
    {
      extractMax(h,v);
      if (v->leader == v)
	if (!v->pr34scanned)
	  {
	    v->age = 0;
	    PRTest34(g, v);
	  }
	else
	  v->age += AGE_INCR;
    }
  clearHeap(h);
        
#ifdef VERBOSE
  printf("c Phase of 34 -- did %d PR3 and %d PR4 in %14.2f time \n",
	 g->PR3Cnt-old3cnt, g->PR4Cnt-old4cnt, timer() - tt);
#endif

  return;
}
#endif /* PR_34 */


/* apply pr tests 1-4 as much as possible in linear time */
void PRpass(graph *g, double bailout, double bailout12, double bailout34,
	    char countAsPreprocessing)
{
  int oldN, oldN1;
  node *v;

  /* do the preprocessing */
  setAuxArcs(g, NULL);
  ForAllNodes(g, v)
    v->auxArc=NULL;
  do {
    oldN = g->currentN;
    do {
      oldN1 = g->currentN;
      g->context = countAsPreprocessing ? PRPRE12 : PR12;
      PR12Pass(g);
    } while ((g->currentN < oldN1 * bailout12) && (g->currentN >= 3));
   
#ifdef PR_34
    do { 
      oldN1 = g->currentN;
      g->context = countAsPreprocessing ? PRPRE34 : PR34;
      PR34Pass(g);
    } while ((g->currentN < oldN1 * bailout34) && (g->currentN >= 3)); 
#endif  /* PR_34 */

  } while ((g->currentN < oldN * bailout) && (g->currentN >= 3));
}

/* source pr test */
void sourcePR(graph *g, node *v, double bailout)
{
  node *w, *u;
  arc *vw, *wu, *a;
  weight_t vCap, sum4;
  int wHere, oldN;

  if (g->currentN <= 2) return;

  setAuxArcs(g, NULL);
  ForAllNodes(g, w) 
    {
      w->pr34scanned = 0;
      w->auxArc = NULL;
    }
  setAuxArcs(g, v);

  do
    {
      oldN = g->currentN;
      vCap = v->cap;
      ForAllIncidentArcs(v,vw) 
	{
	  w = vw->head;
	  g->edgeScanCnt[g->context=PR12]++;

	  /* do tests 1 and 2 */
	  if ((PRTest1(vw) && ++g->PR1Cnt) || 
	      (PRTest2(vw,vCap,w->cap) && ++g->PR2Cnt)) 
	    {
	      /* we also do 34, so we need to be conservative */
	      compactContract(g, v, w, MERGE_INTO_V);
	      if (g->currentN <= 2) return ;
	      vCap = v->cap;
	      continue;
	    }

	  wHere = 1;
	  sum4 = arcCap(vw);

	  assert(vCap == v->cap);
      
	  g->context = PR34;
	  ForAllIncidentArcs(w,wu) 
	    {
	      g->edgeScanCnt[PR34]++;

	      u = wu->head;
	      a = u->auxArc;

	      if ((a != NULL) && (a->head == v)) 
		{
		  sum4 += MIN(arcCap(wu), arcCap(a)); /* accounting for test4*/
		  
		  /* test 3 */
		  if ((vCap <= 2*(arcCap(vw) + arcCap(a))) &&
		      (w->cap <= 2*(arcCap(vw) + arcCap(wu)))) 
		    {
		      g->PR3Cnt++;	
		      compactContract(g, v, w, MERGE_INTO_V);		      
		      wHere = 0; /* w gone */		      
		      vCap = v->cap; 
		      if (g->currentN <= 2)  return;
		      break;
		    }  
		}
	    }

	  if (wHere && (sum4 >= g->minCap)) 
	    {
	      g->PR4Cnt++;
	      compactContract(g, v, w, MERGE_INTO_V);
	      vCap = v->cap; 
	      if (g->currentN <= 2)  return;
	    }
	}
    } while(g->currentN < oldN*bailout && g->currentN >= 3);
}


/* iteratePR and print stats */
void PRpreprocess(graph *g, double bailout, double bailout12, 
		  double bailout34)
{
  float tt = timer();

  PRpass(g, bailout, bailout34, bailout34, 1);

  /* print stats */
  printf("c pnodes:  %14d    ptime: %14.2f\n", g->currentN, timer() - tt); 
  printf("c initial  PR: %-6d PR 1: %-6d PR 2: %-6d PR 3: %-6d PR 4: %-6d\n", 
	 g->PR1Cnt+g->PR2Cnt+g->PR3Cnt+g->PR4Cnt, g->PR1Cnt, g->PR2Cnt, 
	 g->PR3Cnt, g->PR4Cnt);
  printf("c PR initial edge scans  1+2: %-8ld 3+4: %-8ld\n", 
	 g->edgeScanCnt[PRPRE12], g->edgeScanCnt[PRPRE34]);
  g->PR1Cnt = g->PR2Cnt = g->PR3Cnt = g->PR4Cnt = 0;
}

#endif /* HO_INTERNAL */
