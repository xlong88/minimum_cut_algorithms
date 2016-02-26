/* Karger's mincut algorithm
 *
 * Usage: k [-s S]
 *
 * Description: Takes a directed graph in DIMACS format and runs
 * Karger's mincut algorithm on it. The random number generator is
 * seeded with S if specified; otherwise it uses 123456. The value of
 * the cut and the run time are printed to stdout.
 * 
 * Reference: David R. Karger, "Minimum Cuts in Near Linear Time", in 
 * proceedings of the ACM Symposium on the Theory of Computation 1996.
 *
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu) 
 * $Id: k.c,v 1.9 1997/05/16 03:41:52 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#ifdef BROKEN_HDRS   
extern char *optarg;  /* stupid SUNs don't have getopt stuff in headers */
#else
#include <getopt.h>
#endif
#include <time.h>
#include "karger_mincut.h"
#ifdef GABOW
#include "gabow_pack.h"
#endif
#ifdef PST
#include "pst.h"
#endif
#include "pr.h"
#include "random.h"
#include "timer.h"

int main(int argc, char *argv[])
{
  int seed=123456, i;
  double successprob = SUCCESSPROB;
  double gabowbailout = GABOWBAILOUT;
  double maxdyndensity = MAXDYNDENSITY;
  double disrespect = 1.;
  weight_t cut;
  float start, stop;
  graph *g;

  while ((i=getopt(argc, argv, "s:p:r:d:b:")) != -1)
    switch (i)
      {
      case '?': 
	fprintf(stderr, "usage: %s [-s seedvalue] [-p success probability] [-r disrespect fraction] [-b gabow bailout] [-d max dyn density]\n",
		argv[0]);
	exit(1);
      case 's': 
	if (*optarg != '=') seed = atoi(optarg); 
	else seed = atoi(optarg+1); 
	break;
      case 'p':
	if (*optarg != '=') successprob = atof(optarg); 
	else successprob = atof(optarg+1); 
	break;
      case 'r':
	if (*optarg != '=') disrespect = atof(optarg); 
	else disrespect = atof(optarg+1); 
	break;
      case 'd':
	if (*optarg != '=') maxdyndensity = atof(optarg); 
	else maxdyndensity = atof(optarg+1); 
	break;
      case 'b':
	if (*optarg != '=') gabowbailout = atof(optarg); 
	else gabowbailout = atof(optarg+1); 
	break;
      }

  SetRandom(seed);

  /* this leaves the arcs in edge lists in the global nodes, and
     sets the global currentN */
  g = dimacsParse(stdin);

  printf("c nodes: %16d    arcs: %15d\n", g->n, g->m);

  start = timer();

  /* graph input comes from globals set by parse */
  cut = karger_mincut(g, successprob, disrespect, 
		      gabowbailout, maxdyndensity); 

  stop = timer();

#ifndef NO_PR
  printf("c internal PR: %-6d PR 1: %-6d PR 2: %-6d PR 3: %-6d PR 4: %-6d\n", 
	 g->PR1Cnt+g->PR2Cnt+g->PR3Cnt+g->PR4Cnt, g->PR1Cnt, g->PR2Cnt, 
	 g->PR3Cnt, g->PR4Cnt);
  printf("c PR internal edge scans 1+2: %-8ld 3+4: %-8ld\n", 
	 g->edgeScanCnt[PR12], g->edgeScanCnt[PR34]);
#endif

  printf("c edge scans: %11ld\n", totalEdgeScans(g));
  printf("c Packing edge scans:%4ld\n", g->edgeScanCnt[PACKING]);
  printf("c MinCuts discovered:%4d\n", g->cutCount);
  printf("c ttime: %16.2f    capacity: ", stop-start);
  fprintf_wt(stdout, "%11", g->minCap);
  printf("\nc dtime: %16.2f\n", g->dtime-start);
  printf("\n");  
  
#ifdef SAVECUT
  printCut(g);
#endif

  freeGraph(g);

  return 0;
}
