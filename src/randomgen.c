/* min cut problem generator 
 *
 * makes almost regular graphs. The intention was to make graphs where
 * all the vertex degrees are the same, but counting a self loop as 2,
 * so there would be a slight break in regularity. It turns out there
 * is also a bug that moves a few edges around and breaks regularity a
 * little more. oh well. it doesn't make much difference and we
 * already tested with this.
 * 
 * $Id: randomgen.c,v 1.4 1997/01/03 02:33:01 mslevine Exp $ */

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
    fprintf ( stderr, "\nUsage: randomgen n d [-sS]\n");
    fprintf ( stderr, "where n is the number of nodes\n");
    fprintf ( stderr, "      d is the degree\n");
    fprintf ( stderr, "      S is a seed\n");
    break;
  }
  
  case 2: {
    fprintf ( stderr, "\nError: number_of_nodes must be even int");
    break;
  }

  case 3: {
    fprintf ( stderr, "\nError: degree out of range\n");
    break;
  }

  }

  exit(error_no);
}


void main(int argc, char* argv[])
{
  char args[30];
  int n,d,seed;
  int *pv, *rgv, v;
  int v1,v2,i1,i2,i,cv;

  if (( argc < 3 ) || ( argc > 4 )) error (1);
  
  /* first parameter - number of nodes */
  n = atoi (argv[1]);
  if ( (n<2) || (n&1) ) error (2);
 
  /* second parameter - degree */
  d = atoi( argv[2]);
  if ( d > n-1 ) error(3);

  /* optional parameters */
  /* set default values */
  seed = 214365;

  if ( argc == 4 ) {
    strcpy ( args, argv[3]);
    if (( args[0] != DASH ) || ( args[1] != 's')) error (1);
    seed  =  atoi ( &args[2] );
  }

  SetRandom(seed);

  pv = rgv = calloc(n*d,sizeof(int));

  for (v=1; v<=n; v++)
    for (i=0; i<d; i++)
      *(pv++) = v;

  printf ("c RANDOMGEN min-cut problem\n");
  printf ("p cut %8d %8d\n", n, n*d/2 );
    
  for (cv = n*d; cv; cv-= 2)
    {
      /*generate i1 not equal i2*/
      i1 = Random(0,cv-1);
      i2 = Random(0,cv-2);
      if (i2 >= i1) i2++;
      
      printf ("a %d %d %d\n", rgv[i1], rgv[i2], 1);

      /* move (cv-1,cv-2)  to (i1,i2)
         being careful about overlap*/

      v1 = rgv[cv-1];
      v2 = rgv[cv-2];
      rgv[i1] = v1;
      rgv[i2] = v2;
    }

}
