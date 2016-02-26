/* A portable function to print long long ints. 
 *
 * gcc knows about long long ints, so the code involving them is
 * mostly portable, but unfortunately printf is part of libc, not cc,
 * so we get left at the mercy of the system when we go to
 * print. Hence a simple home-grown function to print long long
 * ints...
 * 
 * This is only capable of one long long int at a time, and the
 * "format" string can only be used to specify placement.
 * For example,
 * fprintfll(stdout, "%20", 123456789)
 * gives
 *            123456789
 * fprintfll(stdout, "%-20", 123456789)
 * gives
 * 123456789
 *
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 * $Id: fprintfll.h,v 1.3 1997/01/14 21:10:12 mslevine Exp $
 */

#ifndef FPRINTFLL_H
#define FPRINTFLL_H

#include <stdio.h>

int fprintfll(FILE *fp, char *format, long long int val);

#endif
