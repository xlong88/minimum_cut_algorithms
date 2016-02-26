/*** This is a portable random number generator whose origins are
 *** unknown.  As far as can be told, this is public domain software. */
/* $Id: random.h,v 1.4 1997/02/03 04:37:41 mslevine Exp $ */

/*** SetRandom - initialize constants and seed */
void SetRandom(long seed);

/*** Random - generate a random integer in the interval [a,b] (b >= a >= 0) */
long Random(long a, long b);

/***  DblUnitRandom - generate a random double in the interval [0,1] */
double DblUnitRandom(void);

/*** RandomPermute - randomly permute the iMac elements of rgi */
void RandomPermute(int iMac, int rgi[/*iMac*/]);

/*** RandomPermute - randomly permute the iMac elements of rgi, 
  each of size objsize */
void RandomPermuteObj(int iMac, char rgi[/*iMac*/], int objsize);
