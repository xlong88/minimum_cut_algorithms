/* Read in a graph in DIMACS format, PRpreprocesses() until stuck and
 * print out the remaining graph. 
 *
 * $Id: pr_only.c,v 1.6 1997/04/13 05:54:57 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define PR_ONLY

#include "pr.h"
#include "contract.h"
#include "timer.h"
#include "heap.h"

heap h;

int main (void)
{
  float t;
  graph *g;
  
  /* read input */
  g = dimacsParse(stdin);

  /* initialization */
  makeHeap ( h, g->n );
  printf("c nodes:   %12d    arcs:     %15d\n", g->n, g->m);

  t = timer ();

  compact(g);
  PRpreprocess(g, 1, 0.5, 0.5);
  
  g->dtime = g->dtime - t;
  t = timer() - t;

  printf("c ttime: %16.2f    capacity: ", t);
  fprintf_wt(stdout, "%11", g->minCap);
  printf("\nc dtime: %16.2f\nc\n", g->dtime);

  printGraph(g); 

  freeGraph(g);
  freeHeap(h);

  return 0;
}

