/* min cut problem generator 
   based on the generator described by Padberg and Rinaldi
   $Id: prgen.c,v 1.3 1997/01/03 02:33:00 mslevine Exp $
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include "random.h"

void error(int error_no);

#define RANGE 100
#define DASH '-'
void error (int error_no)
{
  switch ( error_no ) {

  case 1: {
    fprintf ( stderr, "\nUsage: prgen n d c [-sS]\n");
    fprintf ( stderr, "where n is the number of nodes\n");
    fprintf ( stderr, "      d is the density (%%)(double)\n");
    fprintf ( stderr, "      c is the scenario (1 or 2)\n");
    break;
  }

  case 2: {
    fprintf ( stderr, "\nError: number_of_nodes must be int > 1\n");
    break;
  }

  case 3: {
    fprintf ( stderr, "\nError: density out of range\n");
    break;
  }

  case 4: {
    fprintf ( stderr, "\nError: scenario c out of range\n");
    break;
  }

  }

  exit(error_no);
}

int main ( int argc, char* argv[])
{

  char   args[30];
  long   n, m;
  long   i, j;
  long   seed;
  double d;
  long c;
  double r;
  long cap;
  int hack;


  if (( argc < 4 ) || ( argc > 5 )) error (1);

  strcpy ( args, argv[1] );

  /* first parameter - number of nodes */
  if (( n = atoi ( argv[1] ) )  <  2  ) error (2);
 
  /* second parameter - density */
  d = atof ( argv[2] );
  if ( (d > 100) || (d < 0))
    error (3);

  /* third parameter - scenario number */
  c = atoi ( argv[3] );
  if (( c < 1 ) || ( c > 2 ))
    error (4);
 
  /* optional parameters */
  /* set default values */
  seed = 214365;

  if ( argc == 5 ) {
    strcpy ( args, argv[4]);
    if (( args[0] != DASH ) || ( args[1] != 's')) error (1);
    seed  =  atoi ( &args[2] );
  }
  m=0;
  for(hack = 1;hack <=2; hack++) {
    SetRandom(seed);

    if (hack == 2) 
      {
	printf ("c NOIGEN min-cut problem\n");
	printf ("p cut %8ld %8ld\n", n, m );
      }

    for (i=1;i <= n-1; i++) {
      for(j= i+1;j <= n; j++) {
	r = Random (1, 100);
	if (j == i+1)
	  cap = r;
	else
	  {
	    if (r < (100 - d) )
	      cap = 0;
	    else
	      cap = 100 * (1.0 - (double) r/ 100.0) / (d/100.0);
	  }
	if ( ( c==2) && ( (i*2 > n) || (j*2 <= n)))
	  cap *= n;

	if (cap !=0) {
	  if (hack == 1)
	    m++;
	  else
	    printf("a %ld %ld %ld\n",i,j,cap);
	}
      }
    }
  }
  return 0;
}


