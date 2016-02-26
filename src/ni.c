/* Nagomoci et. al. version of the Ibaraki-Nagomoci min-cut algorithm */
/* see Math. Prog. A VOL. 67, 1994, 325--342 */
/* Written by Andrew V. Goldberg - avg@research.nj.nec.edu */
/* NEC Research Institute, Inc. */

/* We attempted to replicate the PR strategy used by NOI in hybrid, but 
 * we got confused by differences between the paper and their code and 
 * eventually dropped it. Our code is at the end of this file. It should be 
 * ignored unless you want to get the NOI paper and code and figure out what 
 * it should be. */

/* 3/4  deleted all things related to outputing a cut  CS */
/* $Id: ni.c,v 1.10 1997/05/16 03:41:53 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define NI

#ifndef PR_FREQ
#define PR_FREQ 2
#endif

#include "graph.h"
#include "pr.h"
#include "timer.h"
#include "heap.h"
#include "ni_sparse_cert.h"
#include "contract.h"

/************************************* global variables **************/
heap   h;

/************************************* definitions *************/

int main (int argc, char *argv[])
{
  graph *g;
  node *v;
  float t;
  int numPhases=0;
  
  /* read input */
  g = dimacsParse(stdin);
  printf("c nodes:   %14d    arcs: %15d\n", g->currentN, g->m);

  t = timer ();

  /* initialization */
  makeHeap(h, g->n);
  compact(g);  /* clean up multiple edges and/or self-loops in input */

#ifndef NO_PR
  PRpreprocess(g, 0.75, 0.5, 0.5);  /* preprocess with PR heuristics */
#endif  

  /* main loop */
  while (g->currentN > 2)  
    {
      ForAllNodes(g, v) 
	{
	  v->key = 0;  /* remove junk left behind by PR */
#ifdef PR_34
	  v->age = 0;  /* clear contract counts */
#endif      
	}

      g->context = MAIN;
      v = sparseCertContract(g, g->minCap);
      if ( g->currentN <= 1 || v==NULL) break;

      compact(g);

      numPhases++;
#ifndef NO_PR
#if PR_FREQ > 1
      if (!(numPhases % PR_FREQ))
#endif
	{
	  sourcePR(g, v, 0.9);
	  PRpass(g, 0.75, 0.5, 0.5, 0);
	}
#endif
    } 
  
  freeHeap(h);

  g->dtime = g->dtime - t;
  t = timer() - t;

  /* print out stats */

#ifndef NO_PR
  printf("c internal PR: %-6d PR 1: %-6d PR 2: %-6d PR 3: %-6d PR 4: %-6d\n", 
	 g->PR1Cnt+g->PR2Cnt+g->PR3Cnt+g->PR4Cnt, g->PR1Cnt, g->PR2Cnt, 
	 g->PR3Cnt, g->PR4Cnt);
  printf("c PR internal edge scans 1+2: %-8ld 3+4: %-8ld\n", 
	 g->edgeScanCnt[PR12], g->edgeScanCnt[PR34]);
#endif
  printf("c phases:  %14d    scans: %14d\n", numPhases, g->numScans);
  printf("c edge scans: %11ld\n", g->edgeScanCnt[MAIN]);
  printf("c MinCuts discovered:%4d\n", g->cutCount);
  printf("c ttime: %16.2f    capacity: ", t);
  fprintf_wt(stdout, "%11", g->minCap);
  printf("\nc dtime: %16.2f\n", g->dtime);

  printf("\n");

#ifdef SAVECUT
  printCut(g);
#endif

  freeGraph(g);
  
  return 0;
}

