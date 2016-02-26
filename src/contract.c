/* contract.h: routines for graph contraction 
* $Id: contract.c,v 1.5 1997/05/16 03:41:51 mslevine Exp $ */

#include <stdlib.h>
#ifdef __CHECKER__
#include "memassert.h"
#else
#include <assert.h>
#endif
#include "contract.h"

#ifdef HO
#include "ho.h"
#endif

/*************************************************** prototypes *******/

/* delete self-loops and merge parallel edges */
static void deleteExtras(graph *g, node *v);

/*************************************************** definitions ******/

/* DeleteArc is very conservative.
   I am not sure it is be necessary to set head to NULL 
   or to update v->current.
   But this helps debugging, and since an arc is deleted only once,
   this is not too inefficient.
*/
void DeleteArc(graph *g, node *v, arc *a) 
{
  if ( a == v -> first ) v -> first = a -> next;
  if ( a == v -> last ) v -> last = a -> prev;
  if ( a -> next != NULL ) a -> next -> prev = a -> prev;
  if ( a -> prev != NULL ) a -> prev -> next = a -> next;
#ifdef HO
  if (v->current == a)
    if (a->next != NULL)
      v->current = a->next;
    else
      v->current = a ->prev;
#endif
#ifdef ENABLE_GRAPH_BACKUPS
  if (a->index != --g->currentM) 
    {
      g->currentArcs[a->index] = g->currentArcs[g->currentM]; 
      g->currentArcs[a->index]->index = a->index;
    }
  g->currentArcs[g->currentM] = NULL; 
#endif
  a->head = NULL;
}	


/* set-union contraction */
void setUnionContract(graph *g, node *v, node *w)
{
#ifndef NO_PR
  v->age += 1;
#endif

  w->leader = v;
  if (INDEX(w) != --g->currentN)
    {
      g->nodes[INDEX(w)] = g->nodes[g->currentN];
      INDEX(g->nodes[INDEX(w)]) = INDEX(w);
    }
  g->nodes[g->currentN] = w;
}

/* find the node this vertex is part of (and do path compression) */
node *findLeader (node *v)
{
  if ( v != v->leader )
    v->leader = findLeader ( v->leader );
  return ( v->leader );
}

/* compactify set-union representation
   do path compression on all vertices, merge parallel arcs, and delete 
   self-loops. */
void compact(graph *g)
{
  node *v, **w;

  for(w=g->nodes+g->currentN; *w != NULL; w++) 
    {
      v=findLeader(*w);
      v->last->next = (*w)->first;
      (*w)->first->prev = v->last;
      v->last = (*w)->last;	
    }
  g->nodes[g->currentN] = NULL;
  
  /* now there are no leader chains until contract is called */

  ForAllNodes(g, v)
    deleteExtras(g, v);
}

static void deleteExtras(graph *g, node *v)
     /* delete parallel edges adjacent to v */
     /* delete self-loops */
{
  node *w;
  arc *a, *b;

  ForAllIncidentArcs(v, a) 
    {
      a->head = a->head->leader;
      a->head->auxArc = NULL;
    }

  v->cap = 0;
  ForAllIncidentArcs (v, a) 
    {
      g->edgeScanCnt[g->context]++;
      w = a->head;
      if ( w == v ) {
	/* delete self-loop */
	DeleteArc (g, v, a );
      }
      else {
	v->cap += arcCap ( a );
	if ( w->auxArc == NULL ) {
	  w->auxArc = a;
	}
	else {
	  b = w->auxArc;
	  b->CAPFIELD += a->CAPFIELD;
	  Reverse ( b )->CAPFIELD +=  Reverse ( a )->CAPFIELD;
	  DeleteArc (g, a->head, Reverse ( a ) );
	  DeleteArc (g, v, a );
	}
      }
    }
  
  if ( v->cap < g->minCap ) 
    saveTrivialCut(g, v, v->cap);
}

/* reset auxArc field if necessary */
void setAuxArcs(graph *g, node *v)
{
  static node *last = NULL;
  arc *a;
 
  if (v==NULL) last=NULL;
  else 
    if (last != v) 
      {
	last = v;
	ForAllIncidentArcs(v,a) 
	  {
	    a->head->auxArc = Reverse(a);
	    g->edgeScanCnt[g->context]++;
	  }
      }
#ifndef NDEBUG
    else
      ForAllIncidentArcs(v,a) {
	assert(a->head->auxArc == Reverse(a));
      }
#endif
}


/* compact contraction 

   WARNING: Although it is not at all obvious, the order of operations
   in this function is extremely important. contract is called
   inside ForAllArcs(v,a) loops. These loops rely on a->next always
   pointing to the next edge in the list, even though a call to
   contract will delete the current edge a from the list!
   DeleteArc() does not damage the next pointer of the deleted arc,
   but it is important that this next pointer point to a valid
   edge. In particular, it is NOT ok if we delete a and then delete
   a->next, because then a->next will then be invalid. Furthermore, if
   a is the last edge of v's list, we will want a->next to be the
   first non-deleted edge of w's list. This means that deleting the
   edge originally from v->w MUST BE THE VERY LAST THING. 

   decrement currentN
   sets leaders 
   merges edgelists*/

void compactContract(graph *g, node *v, node *w, int type)
{
  arc *a, *to_delete;
  arc *ra, *rd, *self_loop=NULL;
  node *x;

  assert(v->leader == v);
  assert(w->leader == w);
  assert(g->currentN > 2);

#ifdef HO
  if (g->source != NULL)  /* flow has started */
    {
      fixFlowForContraction(g, v, w);
      if (g->currentN <= 2) return;
    }
#endif

  w->leader = v;
  if (INDEX(w) != --g->currentN)
    {
      g->nodes[INDEX(w)] = g->nodes[g->currentN];
      INDEX(g->nodes[INDEX(w)]) = INDEX(w);
    }
  g->nodes[g->currentN] = NULL;

#ifdef PR_34
  v->age += 1;
#endif

  setAuxArcs(g, v);  

  /*  scan w's arcs looking for parallel edges and self loops */
  ForAllIncidentArcs(w, a)
    {
      g->edgeScanCnt[g->context]++;
      ra = Reverse(a);
      x = a-> head;
      if (x == v)
	{ /* self loop */
	  /* self_loop must be deleted last, so we save it for later */
	  assert(self_loop == NULL);
	  self_loop = a;   /* pointing w to v */ 
	}
      else 
	{
	  if (type == MERGE_INTO_W)   
	    if ((x->auxArc != NULL) && (x->auxArc->head == v)) 
	      { 
		/* catch for a vertex that never had auxArc */ 
		/* parallel edge, merge earlier edge into this one */
		to_delete = x->auxArc;
		rd = Reverse(to_delete);
		
		a->CAPFIELD += rd->CAPFIELD;
		ra -> CAPFIELD += to_delete->CAPFIELD;
		ra -> head = v;
		
		x->auxArc = ra;  /* reset auxArc to new one */
		
		DeleteArc(g, x,to_delete);
		DeleteArc(g, v,rd);
	      }
	    else 
	      {  /* not a parallel arc */
		x->auxArc = ra;
		ra -> head = v;
	      }
	  else  /* (type == MERGE_INTO_V) */
	    if ((x->auxArc != NULL) && (x->auxArc->head == v))
	      { /* catch for a vertex that never had auxArc */ 
		/* parallel edge, merge this edge into earlier one */
		x->auxArc->CAPFIELD += ra->CAPFIELD;
		Reverse(x->auxArc)->CAPFIELD += a->CAPFIELD;
		
		DeleteArc(g, w,a);
		DeleteArc(g, x,ra);
		/* we use the fact that a->next is not changed by DeleteArc*/
	      }
	    else 
	      {  /* not a parallel arc */
		x->auxArc = ra;
		ra -> head = v;
	      }
	}
    }

  /* must merge lists before deleting self_loop */
  if ( v->first == NULL ) { /* catch case when v's neighbors are all w's 
				neighbors too */
    v->first = w->first;
    v->last = w->last;
  }
  else {
    v -> last -> next = w -> first;
    w -> first -> prev = v -> last;
    v -> last = w -> last;
  }

  /* finally we can delete the self_loop */
  if (self_loop != NULL)
    {
      ra = Reverse(self_loop);
      /* in the next line we subtract 2* the capacity of the self loop
	 from v's capacity. we write it starngely so that when using the
	 residual capaities in HO, teh right thing still happens. */ 
      v -> cap += w->cap - self_loop->CAPFIELD - ra->CAPFIELD;
      if (v->cap < g->minCap) 
	saveTrivialCut(g, v, v->cap);  
            
      DeleteArc(g,v,self_loop);    /* v and w already merged ! */
      DeleteArc(g,v,ra);           /* this is the deletion that MUST BE LAST */
    }
  else {
    v -> cap  += w -> cap;   /* note no chance for new mincap */
  }
}


#ifdef K
/* special form of compactContract for use during K's 2-respect
   computations.  does NOT delete self-loops. does NOT update vertex
   capacities or leaders.  does NOT support being called inside a
   ForAllIncidentArcs() loop, as does comapctContract. */
void compactContractLite(graph *g, node *v, node *w) 
{ 
  arc *a, *to_delete=NULL; 
  arc *ra, *rd;
  node *x;

  v->auxArc = NULL;
  ForAllIncidentArcs(v,a) 
    {   /* we preseve self-loops. they get merged here */
      if (a->head == w)
	{ /* new self loop */
	  if (v->auxArc == NULL)
	    {
	      v->auxArc = a;
	      a->head = v;
	      a->CAPFIELD += a->CAPFIELD;
	    }
	  else
	    {
	      v->auxArc->CAPFIELD += 2*a->CAPFIELD;
	      DeleteArc(g, v, a);
	    }
	  DeleteArc(g, w, Reverse(a));
	}
      else if (a->head == v)
	{ /* old self loop */
	  if (v->auxArc == NULL)
	    v->auxArc = a;
	  else
	    {
	      v->auxArc->CAPFIELD += a->CAPFIELD;
	      DeleteArc(g, v, a);
	    }
	}
      else  /* non-self loop, just set auxArc */
	a->head->auxArc = Reverse(a);
      g->edgeScanCnt[g->context]++;
    }

  /*  scan w's arcs looking for parallel edges and self loops */
  
  ForAllIncidentArcs(w, a)
    {
      g->edgeScanCnt[g->context]++;
      ra = Reverse(a);
      x = a->head;
      if (x==w)
	{
	  to_delete = v->auxArc;
	  a->CAPFIELD += to_delete->CAPFIELD;
	  a->head = v;
	  v->auxArc = a;
	  DeleteArc(g, v, to_delete);
	}
      else if ((x->auxArc != NULL) && (x->auxArc->head == v)) 
	{ 
	  /* catch for a vertex that never had auxArc */ 
	  /* parallel edge, merge earlier edge into this one */
	  to_delete = x->auxArc;
	  rd = Reverse(to_delete);
	  
	  a->CAPFIELD += rd->CAPFIELD;
	  ra -> CAPFIELD += to_delete->CAPFIELD;
	  ra -> head = v;
	  
	  x->auxArc = ra;  /* reset auxArc to new one */
	  
	  DeleteArc(g,x,to_delete);
	  DeleteArc(g,v,rd);
	}
      else 
	{  /* not a parallel arc */
	  x->auxArc = ra;
	  ra -> head = v;
	}
    }
  
  /* concatenate edge lists */
  if ( v->first == NULL ) 
    { /* catch case when v's neighbors are all w's 
	 neighbors too */
      v->first = w->first;
      v->last = w->last;
    }
  else if (w->first != NULL)
    {
      v -> last -> next = w -> first;
      w -> first -> prev = v -> last;
      v -> last = w -> last;
    }
}
#endif

/* print out current contracted  graph */
void printGraph(graph *g)
{ 
  int ca=0;
  node *v;
  arc *a;

  ForAllNodes(g, v) 
    ForAllIncidentArcs(v, a) ++ca;
      
  ca /= 2;  /*double counted twin arcs*/
  
  printf("c reduced problem instance\n");
  printf("p cut %d %d\n", g->currentN, ca);
  ForAllNodes(g, v)
    ForAllIncidentArcs(v,a) 
      {
	if (INDEX(a->head) < INDEX(a->rev->head))
	  {
	    printf("a %d %d ", INDEX(a->head)+1, INDEX(a->rev->head)+1);
	    fprintf_wt(stdout, "%", a->CAPFIELD);
	    printf("\n");
	  }
      }
}


#ifdef ENABLE_GRAPH_BACKUPS
graph *makeBackupGraph(graph *g)
{
  arc *a, *ac;
  graph *ans;
  int i;

  memassert(ans = (graph *)calloc(1, sizeof(graph)));
  
  ans->currentN = g->currentN;
  ans->currentM = g->currentM;

  memassert(ans->nodes = (node**)malloc((g->currentN+1)*sizeof(node *)));
  memcpy(ans->nodes, g->nodes, (g->currentN+1)*sizeof(node *));

  memassert(ans->arcs = (arc*)calloc(g->currentM, sizeof(arc)));
  memassert(ans->currentArcs = (arc**)malloc((g->currentM+1)*sizeof(arc *)));

  ForAllArcs(g, a)
    {
      g->edgeScanCnt[g->context]++;
      i = a->index;
      ans->currentArcs[i] = ac = ans->arcs+i;
      ac->CAPFIELD = a->CAPFIELD;
      ac->head = a->head;
      ac->next = a->next != NULL ? g->arcs+a->next->index : NULL;
      ac->prev = a->prev != NULL ? g->arcs+a->prev->index : NULL;
      ac->rev = g->arcs+a->rev->index;
      ac->index = i;
    }
  ans->currentArcs[ans->currentM] = NULL;
  return ans;
}

void restoreBackupGraph(graph *g, graph *backup)
{
  arc *a, *ac;
  node *v;
  int i = 0;

  g->currentN = backup->currentN;
  g->currentM = backup->currentM;

  memcpy(g->nodes, backup->nodes, (backup->currentN+1)*sizeof(node *));
  ForAllNodes(g, v)
    {
      v->leader = v;
      v->cap = 0;
      v->auxArc = NULL;
      v->index =  i++;
      v->key = 0;
      /* v->{first,last,cap} restored below, v->current lost!*/
    }
  
  ForAllArcs(backup, ac)
    {
      g->edgeScanCnt[g->context]++;
      i = ac->index;
      g->arcs[i] = *ac;
      g->currentArcs[i] = a = g->arcs+i;
      a->head->cap += a->CAPFIELD;
      if (a->rev < a)
	{
	  v = a->rev->head;
	  if (a->prev == NULL)
	    v->first = a;
	  if (a->next == NULL)
	    v->last = a;
	  v = a->head;
	  if (a->rev->prev == NULL)
	    v->first = a->rev;
	  if (a->rev->next == NULL)
	    v->last = a->rev;
	}
    }
  g->currentArcs[g->currentM] = NULL;
}

void freeBackupGraph(graph *backup)
{
  free(backup->nodes);
  free(backup->arcs);
  free(backup->currentArcs);
  free(backup);
}
#endif /* ENABLE_GRAPH_BACKUPS */
