/*
  min cut problem generator

  Usage: dblcyclegen n [-U capacity] 

  written by Andrew V. Goldberg, Matthew Levine
  NECI
  avg@research.nj.nec.com, mslevine@theory.lcs.mit.edu
  $Id: dblcyclegen.c,v 1.2 1997/05/16 03:41:51 mslevine Exp $
*/

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

  case 1: {
    fprintf ( stderr, "\nUsage: dblcyclegen n\n");
    fprintf ( stderr, "                  [-U rim capacity]\n");
    break;
  }

  case 2: {
    fprintf ( stderr, "\nError: unknown option\n");
    break;
  }

  case 3: {
    fprintf ( stderr, "\nError: n must be >5");
    break;
  }

  }

  exit(error_no);
}

int main (int argc, char* argv[])
{

  char   args[30];
  long   n, m;
  long   i, x;
  long   np;
  long   U=1000;

  if (( argc < 2 ) || (argc > 4)) error (1);

  np = 0;
  strcpy ( args, argv[1] );

  /* first parameter - number of nodes */
  if ((n = atoi ( argv[1] ))  <  6) error (3);

  for ( np = 2; np < argc; np++ )  {
    strcpy ( args, argv[np]);
    if ( args[0] != DASH ) error (1);

    switch ( args[1] ) {
      
    case 'U' : {
      U  =  (long) atof ( &args[2] );
      break;
    }

    default: {
      error (2);
      break;
    }
    }
  }

  m = 2*n;

  printf ("c \"dblcycle\" min-cut problem\n");
  printf ("p cut %8ld %8ld\n", n, m );

  /* generate cycle arcs */
  for ( i = 1; i <= n; i++ ) 
    {
      if ((i % (n/2)) == 50)
	printf ("a %ld %ld %ld\n", i, ( i % n) + 1, U-3);
      else
	printf ("a %ld %ld %ld\n", i, ( i % n) + 1, U);

      if ((((i-1) % (n/2)) == 50) || ((i+3) % (n/2)) == 50)
	printf ("a %ld %ld %d\n", i, ( (i+2) % n) + 1, 4);
      else
	printf ("a %ld %ld %d\n", i, ( (i+2) % n) + 1, 1);
    }

  return(0);
}



