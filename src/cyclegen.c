/*
  min cut problem generator
  written by Andrew V. Goldberg
  Computer Science Department, Stanford University
  goldberg@cs.stanford.edu
  $Id: cyclegen.c,v 1.3 1997/01/03 02:32:57 mslevine Exp $
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
    fprintf ( stderr, "\nUsage: cyclegen numer_of_nodes\n");
    fprintf ( stderr, "         [-s seed]\n");
    fprintf ( stderr, "         [-m number_of_arcs]\n");
    fprintf ( stderr, "         [-U cycle_capacity_upper_bound]\n");
    fprintf ( stderr, "         [-L cycle_capacity_lower_bound]\n");
    fprintf ( stderr, "         [-u nom_cycle_capacity_upper_bound]\n");
    fprintf ( stderr, "         [-l non_cycle_capacity_lower_bound]\n");
    break;
  }

  case 2: {
    fprintf ( stderr, "\nError: number_of_nodes must be int > 2\n");
    break;
  }

  case 3: {
    fprintf ( stderr, "\nError: number_of_arcs < number_or_nodes\n");
    break;
  }

  case 4: {
    fprintf ( stderr, "\nError: lower_bound > upper_bound\n");
    break;
  }

  case 5: {
    fprintf ( stderr, "\nError: negative bound\n");
    break;
  }

  case 6: {
    fprintf ( stderr, "\nError: unknown option\n");
    break;
  }

  }

  exit(error_no);
}

int main (int argc, char* argv[])
{

  char   args[30];
  long   n, m;
  long   i, j, k;
  long   seed;
  long   np;
  long   l, u, L, U;

  if (( argc < 2 ) || (argc > 8)) error (1);

  np = 0;
  strcpy ( args, argv[1] );

  /* first parameter - number of nodes */
  if (( n = atoi ( argv[1] ) )  <  3  ) error (2);

  /* optional parameters */
  /* set default values */
  seed = 214365;
  m = n;
  U = 10000;
  L = 1;
  u = 10;
  l = 1;

  for ( np = 2; np < argc; np++ )  {
    strcpy ( args, argv[np]);
    if ( args[0] != DASH ) error (1);

    switch ( args[1] ) {
      
    case 'm' : {
      m  =  (long) atof ( &args[2] );
      break;
    }

    case 's' : {
      seed  =  (long) atof ( &args[2] );
      break;
    }

    case 'l' : {
      l  =  (long) atof ( &args[2] );
      break;
    }

    case 'u' : {
      u  =  (long) atof ( &args[2] );
      break;
    }

    case 'L' : {
      L  =  (long) atof ( &args[2] );
      break;
    }

    case 'U' : {
      U  =  (long) atof ( &args[2] );
      break;
    }
    default: {
      error (6);
      break;
    }
    }
  }
  /* sanity check */

  if ( m < n ) error (3);
  if ( l > u ) error (4);
  if ( L > U ) error (4);
  if (( l < 0 ) || ( u < 0 ) || ( L < 0 ) || ( U < 0 )) error (5);

  SetRandom(seed);

  printf ("c \"cyclegen\" min-cut problem\n");
  printf ("p cut %8ld %8ld\n", n, m );

  /* generate cycle arcs */
  for ( i = 1; i <= n; i++ ) {
    printf ("a %ld %ld %ld\n", i, ( i % n ) + 1, Random ( L, U ));
  }

  /* generate cross arcs */
  for ( k = n+1; k <= m; k++ ) {
    do {
      i = Random ( 1, n );
      j = Random ( 1, n );
    } while ( i == j );
    printf ("a %ld %ld %ld\n", i, j, Random ( l, u ));
  }

  return(0);
}



