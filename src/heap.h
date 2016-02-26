/* heap.h
 * 
 * k-ary heap of nodes. uses the key and heap_pos fields of nodes. 
 *
 * $Id: heap.h,v 1.4 1997/04/11 16:17:08 mslevine Exp $
 */

#ifndef HEAP_H
#define HEAP_H

typedef struct heap_st
{
  long              size;          /* the number of the last heap element */
  node            **node;         /* heap of the pointers to nodes       */ 

  /* temporaries (here because everything macros) */
  long              current_pos,
                    new_pos,
                    pos,
                    last_pos;
  node             *node_k;  
  heap_weight       key_k, key_min;
} heap;

#define HEAP_DEGREE  3
#define NILL        -1

/* internal definition */
#define PUT_TO_POS_IN_HEAP( h, node_i, pos )\
{\
   h.node[pos]        = node_i;\
   node_i -> heap_pos = pos;\
}

#define makeHeap( h, n )\
{\
   h.size = 0;\
   h.node = (node **) calloc (( n+1 ), sizeof(node*));\
}

#define freeHeap(h) \
{ \
    free(h.node); \
}

#define clearHeap( h ) (h.size=0);

#define nonEmptyH( h )  ( h.size > 0 )

#define NODE_IN_HEAP( node_i ) ( node_i -> heap_pos != NILL )

#define increaseKey( h, node_i, key_i ) \
{\
   for ( h.current_pos =  node_i -> heap_pos;\
	h.current_pos > 0;\
	h.current_pos = h.new_pos\
	)\
   { \
      h.new_pos = ( h.current_pos - 1 ) / HEAP_DEGREE;\
      h.node_k = h.node[h.new_pos];\
      if ( key_i  <=  h.node_k -> key ) break;\
      PUT_TO_POS_IN_HEAP ( h, h.node_k, h.current_pos )\
   } \
   PUT_TO_POS_IN_HEAP ( h, node_i, h.current_pos )\
}

#define hInsert( h, node_i )\
{\
   PUT_TO_POS_IN_HEAP ( h, node_i, h.size )\
   h.size ++;\
   increaseKey ( h, node_i, node_i -> key );\
}

#define extractMax( h, node_0 ) \
{\
   node_0             = h.node[0];\
   node_0 -> heap_pos = NILL;\
   \
   h.size -- ;\
   \
   if ( h.size > 0 )\
   {\
      h.node_k =  h.node [ h.size ];\
      h.key_k =  h.node_k -> key;\
      \
      h.current_pos = 0;\
      \
      while ( 1 )\
      {\
         h.new_pos = h.current_pos * HEAP_DEGREE  +  1;\
         if ( h.new_pos >= h.size ) break;\
	 \
         h.key_min  = h.node[h.new_pos] -> key;\
	 \
         h.last_pos  = h.new_pos + HEAP_DEGREE;\
	 if ( h.last_pos > h.size ) h.last_pos = h.size;\
	 \
         for ( h.pos = h.new_pos + 1; h.pos < h.last_pos; h.pos ++ )\
	 {\
	    if ( h.node[h.pos] -> key > h.key_min )\
	    {\
	       h.new_pos = h.pos;\
	       h.key_min  = h.node[h.pos] -> key;\
	   }\
	}\
	 \
         if ( h.key_k >= h.key_min ) break;\
	 \
         PUT_TO_POS_IN_HEAP ( h, h.node[h.new_pos], h.current_pos )\
	 \
         h.current_pos = h.new_pos;\
     }\
      \
      PUT_TO_POS_IN_HEAP ( h, h.node_k, h.current_pos )\
  }\
}

#endif /* ! HEAP_H */
