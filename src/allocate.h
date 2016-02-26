/*
 * Description: Generic "node" allocator (used by gabow_pack)
 *
 * This module provides an allocator object, a data structure that
 * dispenses blocks of memory of a given size. It works by allocating
 * a large block of memory, dispensing it incrementally until it runs
 * out, and getting new large block as needed. It is useful for
 * allocating memory for data structures like linked lists and trees
 * that need many "nodes" of the same size. 
 * It is better than malloc for 2 reasons:  
 *
 *    1. It is generally faster, because knowing that all blocks are
 *    the same size cuts overhead, and it doesn't even try to give memory
 *    back to the system when blocks are recycled.
 * 
 *    2. It is possible to free an entire linked list without traversing
 *    it by freeing the blocks all the nodes were derived from. This is a
 *    fast and foolproof way to avoid memory leaks.
 * 
 * It has the following disadvantage with respect to malloc:
 *
 *    1. You have to keep track of the allocator objects in your code.
 *
 *    2. Memory return to the system is all or nothing. 
 *
 * An allocator is created by calling make_allocator() with a
 * "suggested" number of nodes to allocate memory for in each block
 * and the size of each node in bytes (as returned by sizeof()). (It
 * will use as blocksize the first power of 2 greater than what is needed
 * to satisfy the hint)
 * Nodes are dispensed by get_node().
 * All the memory dispensed by a given allocator can be freed by 
 * calling free_allocator().
 *
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 *
 * $Id: allocate.h,v 1.3 1997/01/14 21:10:12 mslevine Exp $ */

#ifndef ALLOCATE_H
#define ALLOCATE_H

#include <stdlib.h>
#include <stdio.h>
#include "memassert.h"

/*****************************************************************
 * DATA TYPES   
 *****************************************************************/

typedef struct allocator_struct
{
  int nodesize;                   /* size of node to dispense */
  int blocksize;                  /* size of memory blocks */
  char *current;                  /* pointer to next node to give out */
  char *space;                    /* pointer to block of memory 
				     (implicity a linked list of blocks 
				     of memory, the link pointer being 
				     in the space of the first node) */
  char *recycled;
} *allocator;

/*****************************************************************
 * FUNCTION PROTOTYPES   
 *****************************************************************/

allocator make_allocator(int nodecount, int nodesize);  /* constructor */
char *get_node(allocator a);                            /* dispenser */
void recycle_node(allocator a, char *node);
void recycle_all(allocator a);
void free_allocator(allocator a);                       /* destructor */

/*****************************************************************
 * FUNCTION DEFINITIONS
 * (These are here for inline purposes only)
 *****************************************************************/

#if defined(__GNUC__) && !defined(__CHECKER__)
extern inline allocator make_allocator(int nodecount, int nodesize)
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

extern inline char *get_node(allocator a)
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

extern inline void recycle_node(allocator a, char *node)
{
  *(char **)node = a->recycled;
  a->recycled = node;
}

extern inline void recycle_all(allocator a)
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

extern inline void free_allocator(allocator a)
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

#endif
