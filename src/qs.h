/* queue and stack macros for HO
 * $Id: qs.h,v 1.5 1997/04/13 05:54:57 mslevine Exp $*/

#ifndef QS_H
#define QS_H

typedef
   struct queueSt
{
  node             **storage;
  node             **head;
  node             **tail;
  node             **end;
}  queue;

typedef
   struct stackSt
{
  node             **storage;
  node             **tOS;
}  stack;

#define qCreate( size, q ) \
{\
  memassert(q = ( queue * ) calloc ( 1, sizeof ( queue )));\
  memassert(q -> storage = (node**) calloc ( size, sizeof (node*) + 1));\
  q -> end = q -> storage + ( size - 1 );\
}

#define qFree(q) \
{\
   free(q->storage); \
   free(q); \
}

#define qReset( q ) \
{\
  q -> head = q -> tail = q -> storage;\
}

#define qEmpty( q ) ( q -> head == q -> tail )

#define qEnqueue( i, q ) \
{\
  *( q -> tail ) = i;\
  if ( q -> tail == q -> end ) q -> tail = q -> storage;\
  else q -> tail++;\
}

#define qDequeue( i, q ) \
{\
  i = *( q -> head );\
  if ( q -> head == q -> end ) q -> head = q -> storage;\
  else q -> head++;\
}

#define sCreate( size, S ) \
{\
  memassert(S = ( stack * ) calloc ( 1, sizeof ( stack )));\
  memassert(S -> storage = (node**) calloc ( size, sizeof (node*) + 1 ));\
}

#define sFree(s) \
{\
  free(s->storage); \
  free(s); \
}

#define sReset( S ) \
{\
  S -> tOS = S -> storage;\
}

#define sPush( i, S ) \
{\
  *( S -> tOS ) = i;\
  S -> tOS++;\
}

#define sPop( i, S ) \
{\
  S -> tOS--;\
  i = *( S -> tOS );\
}

#define sEmpty( S ) ( S -> tOS == S -> storage )

/* not worrying about empty, etc. because of parallel stack for nodes */

#define nPush(x) \
{\
   nStack[0] += 1; \
   nStack[(int) nStack[0]] = x; \
}


#define nPop(x) \
{\
   x = nStack[ (int) nStack[0] ]; \
   nStack[0] -= 1; \
}


#endif
