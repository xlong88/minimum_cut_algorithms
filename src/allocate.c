/* 
 * Description: Generic "node" allocator (used by gabow_pack)
 *
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 *
 * $Id: allocate.c,v 1.3 1997/01/14 21:10:12 mslevine Exp $
 */

#include "allocate.h"

#ifdef __CHECKER__  
/* lobotomize homegrown allocator if using checker, or else we lose the 
   benefits of using checker! */
allocator make_allocator(int nodecount, int nodesize)
{
  allocator a;

  memassert(a=(allocator)malloc(sizeof(struct allocator_struct)));
  a->nodesize = nodesize;
  return a;
}

char *get_node(allocator a)
{
  return malloc(a->nodesize);
}

void recycle_node(allocator a, char *node)
{
  free(node);
}

void recycle_all(allocator a)
{
  return;
}

void free_allocator(allocator a)
{
  free(a);
}
#else
allocator make_allocator(int nodecount, int nodesize)
{
  allocator a;
  int i;

  memassert(a=(allocator)malloc(sizeof(struct allocator_struct)));
  a->nodesize = nodesize;

  a->blocksize = (nodecount+1)*nodesize;     /* determine min blocksize */
  for(i=1; a->blocksize/=2; i++);            /* find first power of 2 > min */
  a->blocksize = 1<<i;       

  memassert(a->space = (char *)malloc(a->blocksize));  /* get first block */
  a->current = a->space+nodesize;         /* current leaves room for nxt ptr */
  *(char **)a->space = NULL;              /* init linked list */
  a->recycled = NULL;
  return a;
}

char *get_node(allocator a)
{
  char *t;
  
  if (a->recycled != NULL)
    {
      t = a->recycled;
      a->recycled = *(char **)(a->recycled);
      return t;
    }

  if (a->current - a->space + a->nodesize >= a->blocksize) /* out of room? */
    {
      memassert(t=(char *)malloc(a->blocksize));  /* get space */
      *(char **)t = a->space;                  /* link new block into list */
      a->space = t;
      a->current = t+a->nodesize;
    }
  t = a->current;                              /* give next block away */
  a->current += a->nodesize;
  return t;
}

void recycle_node(allocator a, char *node)
{
  *(char **)node = a->recycled;
  a->recycled = node;
}

void recycle_all(allocator a)
{
  char *t;
  /* free extra blocks */
  while (*(char**)a->space != NULL)     /* while there are extra blocks... */
    {
      t = a->space;                     /* walk along */
      a->space = *(char **)a->space;
      free(t);                          /* free front */
    }
  /* reset */
  a->recycled = NULL;
  a->current = a->space+a->nodesize;
}


void free_allocator(allocator a)
{
  char *t;
  
  while (a->space != NULL)              /* while there are blocks... */
    {
      t = a->space;                     /* walk along */
      a->space = *(char **)a->space;
      free(t);                          /* free front */
    }
  free(a);
}

#endif
