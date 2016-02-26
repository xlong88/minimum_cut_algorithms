/* min cut problem generator
 * makes random regular graphs, either unions of random cycles or unions of 
 * random matchings
 * 
 * $Id: regulargen.c,v 1.3 1997/01/03 02:33:01 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include "random.h"

#define DASH '-'

void error (int error_no);
  
void error (int error_no)
{
  switch ( error_no ) {

  case 1: 
    fprintf ( stderr, "\nUsage: regulargen n k [-sS] [-m] [-wW] \n");
    fprintf ( stderr, "where n is the number of nodes\n");
    fprintf ( stderr, "      k is the number of cycles/matchings\n");
    fprintf ( stderr, "      S is a seed\n");
    fprintf ( stderr, "      W is the maximum edge weight\n");
    fprintf ( stderr, "-m means matchings instead of cycles\n");
    break;
  case 2: 
    fprintf (stderr, "error: number of nodes must exceed 1\n");
    break;
  case 3:
    fprintf (stderr, "error: number of nodes must be even\n");
    break;
  }
  exit(17);
}

int main ( int argc, char *argv[])
{

  char   args[30];
  long   n;
  int fMatchings = 0;
  long   i;
  long   seed=214365;
  int d;
  int    *rgi;
  int wt, wtMax=1;
  int *rgwt;

  if (argc < 2) error(1);

  /* first parameter - number of nodes */
  if (( n = atoi ( argv[1] ) )  <  2  ) error (2);

  rgi = malloc(n*sizeof(int));
  
  /* second parameter - cycles */
  d = atof ( argv[2] );

  argv +=3; argc-=3;
  
  while(argc--){
    strcpy ( args, *(argv++));
    if ( args[0] != DASH ) error (1);
    switch(args[1])
      {
      case 'w':
	wtMax= atoi( &args[2] );
	break;
      case 'm':
	fMatchings=1;
	break;
      case 's':
	seed  =  atoi ( &args[2] );
	break;
      default:
	error(1);
      }
  }

  if (fMatchings && (n&1))
    error(3);

  SetRandom(seed);

  rgwt = malloc(d*sizeof(int));
  wt = 0;
  for (i=0; i< d; i++){
    rgwt[i] = 1+Random(1,wtMax);
    wt += rgwt[i];
  }
/*  for (i=0; i< d; i++)
    rgwt[i] /= wt; */

  printf ("c REGULARGEN min-cut problem\n");
  printf ("p cut %8ld %8ld\n", n, fMatchings? n*d/2: n*d );

  while(d--){
    for(i=0; i<n; i++)
      rgi[i] = i+1;
    RandomPermute(n,rgi);
    for(i=0; i<n-1; i+= 1+fMatchings)
      {
	printf("a %d %d %d\n",rgi[i],rgi[i+1],wtMax==1? 1:rgwt[d]);
      }
    if (!fMatchings)
      {
	printf("a %d %d %d\n", rgi[n-1], rgi[0], wtMax==1? 1:rgwt[d]);
      }
  }
  return 0;
}
