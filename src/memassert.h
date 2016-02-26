/* homegrown assert specifically for checking allocations.

 * The point is that some assertions are gratitous sanity checks and
 * should be removed when not debugging, but some are
 * important. Defining NDEBUG removes assert() statements, so they are
 * good for the former; memassert() is intended for checking
 * allocations, which is an instance of the latter.
 *
 * Sample use: memassert(arry = (int *)malloc(10*sizeof(int)))
 *
 * $Id: memassert.h,v 1.3 1997/04/13 05:54:56 mslevine Exp $ */

#ifndef MEMASSERT_H
#define MEMASSERT_H

#include <stdio.h>
#include <stdlib.h>

#define memassert(x) { if ((x)==NULL) {fprintf(stderr, "Memory exhausted at " __FILE__ ":%d " #x "\n", __LINE__); exit(1); } }

/* hack to deal with fact that checker doesn't like assert */
#ifdef __CHECKER__
#define assert(x) { if (!(x)) {fprintf(stderr, "Assertion failed at " __FILE__ ":%d " #x "\n", __LINE__); exit(1); } }
#endif

#endif
