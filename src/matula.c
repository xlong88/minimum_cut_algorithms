/* Matula's 2+epsilon approximation algorithm */

/* $Id: matula.c,v 1.5 1997/04/14 03:22:27 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>

#define MATULA

#include "graph.h"
#include "contract.h"
#include "pr.h"
#include "timer.h"
#include "heap.h"
#include "ni_sparse_cert.h"

heap   h;

int main (int argc, char *argv[])
{
  graph *g, *backup;
  node *v;
  float t;
  double epsilon=1.;
  int i, numPhases = 0;

  while ((i=getopt(argc, argv, "e:")) != -1)
    switch (i)
      {
      case '?': 
	fprintf(stderr, "usage: %s [-e epsilon]\n", argv[0]);
	exit(1);
      case 'e': 
	if (*optarg != '=') epsilon = atof(optarg); 
	else epsilon = atof(optarg+1); 
	break;
      }

  /* read input */
  g = dimacsParse(stdin);
  printf("c nodes:   %14d    arcs: %15d\n", g->currentN, g->m);

  t = timer ();

  /* initialization */
  makeHeap ( h, g->currentN);
  compact(g);  /* clean up multiple edges and/or self-loops in input */

#ifndef NO_PR
  PRpreprocess(g, 0.75, 0.5, 0.5);  /* preprocess with PR heuristics */
#endif  

  backup = makeBackupGraph(g);
  
  if (g->currentN > 2)
    {
      numPhases++;
      matulaApprox(g, epsilon);
      v = sparseCertContract(g, g->minCap);
      if (g->currentN > 1)
	{
	  compact(g);
	  sourcePR(g, v, 0.9);
	  PRpass(g, 0.75, 0.5, 0.5, 0);
	}
    }
 
  restoreBackupGraph(g, backup);
  sparseCertContract(g, g->minCap);

  freeBackupGraph(backup);
  freeHeap(h);

  g->dtime = g->dtime - t;
  t = timer() - t;

  printf("c sparse_nodes: %d  sparse_arcs: %d\n", g->currentN, g->currentM/2);

  /* print out stats */
#ifndef NO_PR
  printf("c internal PR: %-6d PR 1: %-6d PR 2: %-6d PR 3: %-6d PR 4: %-6d\n", 
	 g->PR1Cnt+g->PR2Cnt+g->PR3Cnt+g->PR4Cnt, g->PR1Cnt, g->PR2Cnt, 
	 g->PR3Cnt, g->PR4Cnt);
  printf("c PR internal edge scans 1+2: %-8ld 3+4: %-8ld\n", 
	 g->edgeScanCnt[PR12], g->edgeScanCnt[PR34]);
#endif

  printf("c phases:  %14d    scans: %14d\n", numPhases, g->numScans);
  printf("c edge scans: %11ld\n", totalEdgeScans(g));
  printf("c MinCuts discovered:%4d\n", g->cutCount);
  printf("c ttime: %16.2f    capacity: ", t);
  fprintf_wt(stdout, "%11", g->minCap);
  printf("\nc dtime: %16.2f\n", g->dtime);

  printf("\n");

  freeGraph(g);
  
  return 0;
}
