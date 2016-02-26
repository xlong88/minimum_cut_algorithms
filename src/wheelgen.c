/*
  min cut problem generator
  written by Andrew V. Goldberg
  NECI
  avg@research.nj.nec.com
  $Id: wheelgen.c,v 1.3 1997/01/03 02:33:01 mslevine Exp $
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
    fprintf ( stderr, "\nUsage: cyclegen n\n");
    fprintf ( stderr, "                  [-r range]\n");
    fprintf ( stderr, "                  [-s seed]\n");
    
    break;
  }

  case 2: {
    fprintf ( stderr, "\nError: unknown option\n");
    break;
  }

  case 3: {
    fprintf ( stderr, "\nError: r must be between 0 and n-3\n");
    break;
  }

  }

  exit(error_no);
}

int main (int argc, char* argv[])
{

  char   args[30];
  long   n, m;
  long   i;
  long   seed;
  long   R;
  long   np;
  long   U, u;

  if (( argc < 2 ) || (argc > 4)) error (1);

  np = 0;
  strcpy ( args, argv[1] );

  /* first parameter - number of nodes */
  if (( n = atoi ( argv[1] ) )  <  3  ) error (2);

  /* optional parameters */
  /* set default values */
  seed = 214365;
  R = 0;

  for ( np = 2; np < argc; np++ )  {
    strcpy ( args, argv[np]);
    if ( args[0] != DASH ) error (1);

    switch ( args[1] ) {
      
    case 'r' : {
      R  =  (long) atof ( &args[2] );
      break;
    }

    case 's' : {
      seed  =  (long) atof ( &args[2] );
      break;
    }

    default: {
      error (2);
      break;
    }
    }
  }
  /* sanity check */

  if ((R < 0) || ( R > n - 2)) error(3);

  SetRandom(seed);
  m = 2 * (n-1);
  u = 2;
  U = (n - 2) - R;

  printf ("c \"wheelgen\" min-cut problem\n");
  printf ("p cut %8ld %8ld\n", n, m );

  /* generate cycle arcs */
  for ( i = 1; i < n; i++ ) {
    printf ("a %ld %ld %ld\n", i, ( i % (n-1)) + 1, U + Random ( 0, 2*R ));
  }

  /* generate spoke arcs */
  for ( i = 1; i < n; i++ ) {
    printf ("a %ld %ld %ld\n", n, i, u );
  }

  return(0);
}



