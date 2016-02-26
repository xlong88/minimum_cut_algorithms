/* Gabow's minimum cut algorithm
 *
 * Usage: g
 *
 * Description: Takes a directed graph in DIMACS format and runs
 * Gabow's minimum cut algorithm on it. Prints the cut value and time
 * used to stdout. 
 *
 * Reference: Harold N. Gabow, "A Matroid Approach to Finding Edge
 * Connectivity and Packing Arborescences", Journal of Computer and
 * System Sciences 50, 259-273 (1995)
 * 
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 * $Id: g.c,v 1.3 1997/01/03 02:32:57 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include "memassert.h"
#define K
#define GABOW
#include "gabow_pack.h"
#include "graphutil.h"
#include "random.h"
#include "timer.h"
#include "pr.h"

#ifndef NO_PR
heap h;
#endif

int main(int argc, char *argv[])
{
  int n, m, maxk;
  weight_t cut;
  float start, stop;
  packing pack;
  node *v;
  arc *a;
  
  parse(&n, &m);

  printf("c nodes: %14d    arcs: %16d\n", n, m);

  start = timer();
  
  minCap = MAXWEIGHT;
  ForAllNodes(v) deleteExtras(v);  

#ifndef NO_PR
  makeHeap(h, n);
  n = PRpreprocess(n); 
#endif  

  maxk = minCap;

  ForAllLeaderNodes(v)  /* put input in right place in data structures */
    {
      v->s_first = v->first;
      ForAllArcs(v, a)
	{
	  a->s_prev = a->prev;
	  a->s_next = a->next;
	  a->s_cap = a->cap;
	}
    }

  pack = gabow_pack(findLeader(nodes), n, nodes, maxk);
  cut = pack->k;
  free_packing(pack);

  stop = timer();

  printf("c ttime: %14.2f    capacity: " stop-start);
  fprintf_wt(stdout, "%12", cut);
  printf("\n\n");

  freegraph();

  return 0;
}
