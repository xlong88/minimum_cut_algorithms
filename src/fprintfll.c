/* fprintfll: fprintf a long long int
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 * $Id: fprintfll.c,v 1.5 1997/04/13 05:54:54 mslevine Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __CHECKER__
#include "memassert.h"
#else
#include <assert.h>
#endif

#ifdef __GNUC__   /* this file can't be compiled without it! */

#include "fprintfll.h"

int fprintfll(FILE *fp, char *format, long long int val)
{
  char s[22], myformat[20];
  int i=21, j=0;  

  s[21] = 0;
  if (val < 0)
    {
      j = 1;
      val = -val;
    }
  else if (val == 0)
    s[--i] = '0';
  while(val)
    {
      s[--i] = '0'+val%10;
      val /= 10;
    }
  if (j) s[--i] = '-';
  
  j = strlen(format);
  assert(j <= 18);
  strcpy(myformat, format);
  myformat[j] = 's';
  myformat[j+1] = 0;
  return fprintf(fp, myformat, s+i);
}

#endif
