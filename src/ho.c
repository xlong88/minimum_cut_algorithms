/* Hao-Orlin min cut algorithm 
 *
 * Reference: J. Hao and J. B. Orlin, "A Faster Algorithm for Finding
 * the Minimum Cut in a Directed Graph". Journal of Algorithms, vol
 * 17, pp 424-446, 1994.
 *
 * $Id: ho.c,v 1.10 1997/05/16 03:41:52 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __CHECKER__
#include "memassert.h"
#else
#include <assert.h>
#endif

#define HO

/************************** constants ********************************/

#define GLOB_UPDATE_FREQ 2.0
#define INT_PR_FREQ_12  18.0  /* this should be played with some more */
#define INT_PR_FREQ_34   9.0  /* this should be played with some more */

/************************** includes *******************************/

#include "graph.h"
#include "pr.h"
#include "timer.h"
#include "ho.h"
#define HO_INTERNAL
#include "pr.h"
#include "heap.h"

/**************************** global variables **********************/

#ifdef PR_34
heap h;
#endif

/************************************** prototypes *******************/

static void insertActive(graph *g, node *j, bucket *l);
static void insertNode(node *i, bucket *l);
static void deleteNode(node *i, bucket *l);
static node *FirstRealActive(bucket *l);
static node *NextRealActive(node *i);

static void pushFlow(graph *g, arc *a, node *i, node *j, int jD,
		     weight_t delta);
static int relabel(graph *g, node *i);
static int discharge(graph *g, node *i);
static void saturateOutgoing(graph *g, node *w);

static void gap(graph *g, bucket *gapB);
static void globalUpdate(graph *g);
static void saveFlowCut(graph *g);
static int STCut(graph *g);

static void oneNodeLayer (graph *g, node *i);
#ifdef EXCESS_DETECTION
static void excessGap (graph *g, node *i);
static int excessCheck(graph *g, node *j);
#endif

static int restartCutComp(graph *g);
static void initCutComp(graph *g);

static void mainInit(graph *g);
static void mainCleanup(graph *g);

void computeCuts(graph *g);

/************************************** function definitions *********/

/* ================================================================ */

#ifdef CONTRACT
void fixFlowForContraction(graph *g, node *v, node *w)
{
  bucket *bb, *l;

  assert(v->status != FROZEN);
  assert(v->leader == v);
  assert(w->leader == w);
  assert(w->status != SOURCE);
  
  w->toContract = 1;
  
  if (v->status == SOURCE)
    saturateOutgoing(g, w);

  assert(w->leader == w);  

  if (g->currentN > 2) 
    {
      if (w->status == REGULAR)
	g->regN--;

      if ((w->status != FROZEN) && (w != g->sink)) 
	{
	  bb = g->buckets + w->d;
	  deleteNode(w,bb);
	}
      
      if (( v != g->sink) && (v != g->source) && 
	  (v->excess == 0) && (w->excess > 0)) 
	{
	  l = g->buckets + v->d;
	  insertActive(g,v,l);
	}
      
      v->excess += w->excess;
    }

  w->status = CONTRACTED;
}
#endif

/* ================================================================ */
/* insert node in bucket's active list */
static void insertActive(graph *g, node *j, bucket *l) 
{
  j->nextA = l->firstActive;
  l->firstActive  = j;
  if (j->d > g->aMax)
    g->aMax = j->d;
}

/* ================================================================ */
/* insert node in bucket */
static void insertNode(node *i, bucket *l)
{
  node *tmpNode; 

  tmpNode = l->first;
  i->bNext = tmpNode;
  if ( tmpNode != NULL ) tmpNode->bPrev = i;
  l->first = i;
}

/* ================================================================ */
/* delete node from bucket */
static void deleteNode(node *i, bucket *l)
{
  node *tmpNode;

  if ( (l)->first == (i) )
    { 
       (l)->first = (i)->bNext;
	if ((i)->bNext != NULL) (i)->bNext->bPrev = NULL; 
    } 
  else 
    {
      tmpNode = (i)->bNext;
      (i)->bPrev->bNext = tmpNode;
      if ( tmpNode != NULL ) tmpNode->bPrev = (i)->bPrev;
    }
}

/* ================================================================ */
/* finds the first non-frozen active node */
static node *FirstRealActive(bucket *l)
{
  node *i;

  for (i = l->firstActive; i != NULL; i = i->nextA) 
    if (i->status == REGULAR)
      {
	assert(i->toContract == 0);
	l->firstActive = i;
	return i;
      }
  l->firstActive = NULL;
  return NULL;
}

/* ================================================================ */
/* finds the next non-frozen active node */
static node *NextRealActive(node *i)
{
  node *j;

  for (j = i->nextA; j != NULL; j = j-> nextA) 
    if (j->status == REGULAR)
      {
	assert(j->toContract == 0);
	return j;
      }
  return NULL;
}


/* ================================================================ */
/* push flow and take care of bucket updates */

static void pushFlow(graph *g, arc *a, node *i, node *j, int jD,
		     weight_t delta)
{
  a->resCap -= delta;
  Reverse(a)->resCap += delta;

  if ( j != g->sink && j->excess == 0 )  
    insertActive (g, j, g->buckets+jD );
  
  j->excess += delta;
  i->excess -= delta;
}

/* ================================================================ */
/*--- relabelling node i */
/* assumes i->excess > 0 */
static int relabel(graph *g, node *i)
{
  node  *j;
  int  minD;       /* minimum d of a node reachable from i */
  arc   *minA=NULL; /* an arc which leads to the node with minimal d */
  arc   *a;
  bucket *l;        /* bucket containing node i */

  g->relabelCnt++; 
  assert(i->status == REGULAR);

  l = g->buckets + i->d;
  deleteNode ( i, l );

  minD = g->sourceD;

  /* find the minimum */
  ForAllIncidentArcs(i,a) 
    {
      g->edgeScanCnt[MAIN]++;
      if ( a->resCap > 0 ) 
	{
	  j = a->head;
	  if ((j->status == REGULAR ) && ( j->d < minD )) 
	    {
	      minD = j->d;
	      minA = a;
	    }
	}
    }
  
  minD++;
  
  if ( minD < g->sourceD ) 
    {
      i->d = minD;
      i->current = minA;
      
      l = g->buckets + minD;
      
      insertNode(i, l);
      if ( minD > g->dMax ) g->dMax = minD;
      
      insertActive(g,i , l);
    }
  else
    oneNodeLayer(g, i);

  return(minD);

} /* end of relabel */

/* ================================================================ */
/* discharge:  push flow out of i until i becomes inactive or i is relabeled */
static int discharge (graph *g, node  *i)
{
  node  *j;                /* successor of i */
  int  jD;                /* d of the next bucket */
  arc   *a;                /* current arc (i,j) */
  weight_t  delta;           /* flow to push through the arc */
  int flag;
#ifdef EXCESS_DETECTION
  int zeroExcess=0;
#endif 

  jD = (i->d) - 1;

  /* scanning arcs out of  i  */
  assert(Reverse(i->current)->head == i);
  a = i->current;
  while (a != NULL) {

    assert(Reverse(a)->head == i);
    flag = 0;

    g->edgeScanCnt[g->context]++;
    if ( a->resCap > 0 ) {
      j = a->head;

      if (( j->status == REGULAR ) && ( j->d == jD )) {
 	g->pushCnt ++; 
	delta = MIN ( i->excess, a->resCap );

	pushFlow (g, a, i, j, jD, delta );

	if (j->excess >= g->minCap)
	  if (j == g->sink ) {
	    i->current = a;
	    return (2);
	  }
#ifdef EXCESS_DETECTION
	/* we want to excessCheck(j) and, if excessCheck succeeds,
	   return appropriate value.
	   It is possible to continue if i was not reinserted into the
	   active list, but it is possible that aMax increased and
	   i no longer has max label among active nodes
	   */
	  else {
	    zeroExcess = (i->excess == 0);
	    if (excessCheck(g, j)) {
#ifndef NO_PR
	      if (!isLeader(i))
#else
	      if (i->status == SOURCE)
#endif
		return(4); /* i got contacted */
	      if (!zeroExcess) {
		assert(i->excess > 0);
		insertActive(g,i, (g->buckets+i->d));
		return(3); /* i just got reinserted */
	      }
	      else if (i->excess > 0)
		return(3); /* i got reinserted during excessCheck(j) */
	    } 
	  } 

#endif

	if (i->excess == 0) break;

      } /* j belongs to the next bucket */
    } /* a  is not saturated */

    a = a->next;
    
  } /* end of scanning arcs from  i */


  if (a == NULL) {
    return (1);
  }
  else {
#ifdef CONTRACT
    if (a->head != NULL)
      i->current = a;
    else { /* a was deleted */
      if (a->next != NULL)
	i->current = a->next;
      else
	i->current = a->prev; 
    }
#else
    i->current = a;   /* a never deleted if no contractions */
#endif
    return (0);
  }
}



/* ================================================================ */
/* saturateOutgoing is called only before contracting w
   into the source */
static void saturateOutgoing(graph *g, node * w)
{
  arc *a;
  weight_t delta;
  int jD;
  node *j;

  ForAllIncidentArcs ( w, a ) {
    g->edgeScanCnt[g->context]++;
    j = a->head;
    /* some nodes may get contracted by excess detection */
#ifdef EXCESS_DETECTION
    if (j == NULL) continue;
    if (j->toContract) continue;
#endif
    
    delta = a->resCap;
    if ( delta > 0 ) {
      switch ( j->status ) {
	
      case REGULAR: {
	jD = j->d;
	pushFlow (g, a, w, j, jD, delta );

#ifdef EXCESS_DETECTION
	excessCheck(g, j);       
#endif

	break;
      }
      case FROZEN: {
	a->resCap = 0;
	Reverse ( a )->resCap += delta;
	      
	w->excess -= delta;
	j->excess += delta;

#ifdef EXCESS_DETECTION
	excessCheck(g, j);       
#endif

	break;
      }
      case SOURCE: 
	{ 
	  break; 
	}
      } /* switch */
    }
  }
}

/* ================================================================ */
/* gap relabeling */
static void gap (graph *g, bucket *gapB)
{
  bucket *b;
  node  *i;
  int emptyLayer = 1;
  int  r;           /* index of the bucket before the gap */

  g->gapCnt++;

  r = ( gapB - g->buckets ) - 1;

  /* freeze nodes beyond the gap */

  if (( gapB + 1 )->first == NULL ) {
    i = gapB->first;
    gapB->firstActive = gapB->first = NULL;
    oneNodeLayer(g, i);
    emptyLayer = 0;
  }
  else  { 
    sPush ( NULL, g->freezer );
    for ( b = gapB; b <= g->buckets + g->dMax; b++ ) {
      for ( i = b->first; i != NULL; i = i->bNext ) {
	i->status = FROZEN;
	sPush ( i, g->freezer );
	emptyLayer = 0;
	g->regN--;
	g->gNodeCnt++;
      }
      b->firstActive = b->first = NULL;
    }
  }

  assert(!emptyLayer);

  g->dMax = r;
  g->aMax = r;

}
/* ================================================================ */


/* ================================================================ */
/* global update via backward breadth first search from the sink */
/* update time proportional to the number of arcs adjacent to regular nodes */

static void globalUpdate(graph *g)
{
  node  *i, *j;
  arc   *a;
  bucket *b, *bS;
  int  dist;
  int  fCnt = 0; /* count of freezed nodes */

  g->updateCnt++;
  g->relsSinceUpdate = 0;

  /* initialization */

  bS = g->buckets + g->dMax;
  for ( b = g->buckets + g->sinkD; b <= bS; b++ ) {
    for ( i = b->first; i != NULL; i = i->bNext ) {
      i->d += g->infD; /* this way we can reconstruct i->d */
      sPush ( i, g->regNodes );
    }
    b->first = b->firstActive = NULL;
  }

  g->dMax = g->aMax = g->sinkD;

  /* backward breadth first search from the sink */

  qEnqueue ( g->sink, g->BFSqueue );  
  g->sink->d = g->sinkD;
  do {
    qDequeue ( i, g->BFSqueue );

    dist = i->d;
    b = g->buckets + dist;
    insertNode( i, b );
    if ( dist > g->dMax ) g->dMax = dist;

    if (( i->excess  > 0 ) && ( i != g->sink )) {
      insertActive (g,i, b);
    }

    /* scanning arcs incident to node i */
    dist++;
    ForAllIncidentArcs(i, a ) { 
      g->edgeScanCnt[GLUPDATE]++;
      if ( Reverse ( a )->resCap > 0 ) {
	j = a->head;
	
	if (( j->status == REGULAR ) && ( j->d >= g->infD )) {
	  j->d = dist;
	  j->current = j->first;
    	  qEnqueue ( j, g->BFSqueue );
    	}
      }
    }
  } while ( !qEmpty ( g->BFSqueue ) );

  /* freeze gapped nodes */
  sPush ( NULL, g->freezer );
  do {
    sPop (i, g->regNodes);
    if (i->d >= g->infD) {
      i->d -= g->infD; /* restore label */
      i->status = FROZEN;
      g->regN--;
      sPush ( i, g->freezer );
      fCnt++;
    }
  } while ( !sEmpty ( g->regNodes ));

  if (fCnt == 1) 
    {
      sPop(i, g->freezer);
      fCnt--;
      i->status = REGULAR;
      g->regN++;
      oneNodeLayer(g, i);
    }

  if (fCnt == 0) sPop(i, g->freezer);
}

/* ================================================================ */
/* save minimum cut value found from a flow computation. */
static void saveFlowCut(graph *g) 
{
#ifdef SAVECUT
  node *v;
#endif

  g->minCap = g->sink->excess;
  g->cutCount++;
#ifdef VERBOSE
  printf("c New mc valu = ");
  fprintf_wt(stdout, "%", g->minCap);
  printf("\n");
#endif
  g->dtime = timer();

#ifdef SAVECUT
  ForAllNodes(g, v) 
    v->in_cut = (v->status == REGULAR);
  g->sink->in_cut = 1;
#ifdef CONTRACT
  ForAllVertices(g, v)
    v->in_cut = findLeader(v)->in_cut;
#endif
#endif
}

/* ================================================================ */
/* find a minimum s-t cut */
/* based on the push-relabel method */
static int STCut(graph *g)
{
  node   *i;
  bucket  *l;             /* current bucket */

  /* main loop */
  g->STCutCnt++;
  g->totalSize += g->regN + 1;
  g->relsSinceUpdate = 0;
  while (g->aMax >= g->sinkD) {
    l = g->buckets + g->aMax;

    i = FirstRealActive(l);   /* skip over contracted vertices */
    if ( i == NULL ) {
      g->aMax--;
    }
    else {
      l->firstActive = NextRealActive(i); /* delete i */

      assert(i->status == REGULAR);

      switch (discharge(g, i))
	{
	case 4: 
	  {
	    /* i was contracted into source */
	    break;
	  }
	  
	case 3: {
	  /* i got reinserted into the active list */
	  break;
	}
	case 2: {
	  /* current preflow value at least minCapValue */
	  /* early termination */
	  if (i->excess > 0) {
	    insertActive(g,i,l);
	  }
	  return (1); 
	}
	case 1:
	  {
	    /* i must be relabeled */

	    if ( l->first->bNext == NULL ) {
	      /* relabeling would create a gap */
	      gap (g, l);
	    }
	    else {
	      relabel(g, i);
	      g->relsSinceUpdate ++;

	      /* is it time for global update? */
#ifdef GLOBAL_UPDATES
	      if ( g->relsSinceUpdate > GLOB_UPDATE_FREQ * (double)g->regN ) {
		globalUpdate(g);
	      }
#endif

	    }
	    break;
	  }
	case 0:
	  {
	    break;
	  }
	}
    }
  }
  return (0);
}

/* ================================================================ */
/* 
   a common one-node layer case
   handle in a special way for efficiency

   relabeled node with no residual arcs to regular nodes
   must be handled in a special way in any case
   */
static void oneNodeLayer(graph *g, node *i)
{
#ifndef CONTRACT
  weight_t delta;
  node *j;
#endif
  arc *a;
  weight_t cutCap;

  g->oneCnt++;
  /* 
     when i would have become sink, all currently regular nodes
     would be sources and would push to i
     when i becomes source, it pushes flow to the currently 
     frozen nodes
     */

#ifdef CONTRACT
  g->regN--;
  i->status = FROZEN;  /* so that we don't try to delete it from buckets */
  cutCap = i->excess;
  ForAllIncidentArcs(i, a ) 
    {
      g->edgeScanCnt[g->context]++;
      if (a->head->status == REGULAR)
	cutCap += Reverse (a)->resCap;
    }
  compactContract(g, g->source, i, MERGE_INTO_W);
#else
  assert(i->status == REGULAR);
  i->status = SOURCE;
  i->d = g->sourceD;
  g->regN--;
  g->currentN--;
  g->gNodeCnt++;

  cutCap = i->excess;
  ForAllIncidentArcs(i, a ) {
    g->edgeScanCnt[g->context]++;
    j = a->head;

    /* sum capacity from regular nodes */
    if (j->status == REGULAR)
      {
	/*	assert((j->toContract != 0) || (a->resCap == 0)); */
	cutCap += Reverse ( a )->resCap;
      }

    /* saturate outgoing to frozen nodes */
    if ( j->status == FROZEN ) {
      delta = a->resCap;
      a->resCap = 0;
      Reverse ( a )->resCap += delta;
      j->excess += delta;

#ifdef EXCESS_DETECTION
      excessCheck(g, j);
#endif
    }
  }
#endif
  if ( cutCap < g->minCap ) 
    saveTrivialCut(g, i, cutCap);
}
/* ================================================================ */

#ifdef EXCESS_DETECTION

/* ================================================================ */
/*
   When node with large excess is contrracted into the source
   we check if this created a gap
*/
static void excessGap (graph *g, node *i)
{
  bucket *b;
  node  *j;
  int  r;           /* index of the bucket before the gap */
  int fCnt = 0;     /* number of frozen nodes */

  g->gapCnt++;
  b = g->buckets + (i->d + 1);
  r = i->d;

  sPush (NULL, g->freezer);
  for ( ; b <= g->buckets + g->dMax; b++ ) {
    for (j = b->first; j != NULL; j = j->bNext) {
      assert(j != g->sink);
      j->status = FROZEN;
      sPush (j, g->freezer);
      fCnt++;
      g->regN--;
      g->gNodeCnt++;
    }
    b->firstActive = b->first = NULL;
  }

  if (fCnt == 1) {
    sPop(j, g->freezer);
    if (j->toContract) {
      sPush(j, g->freezer);  /* WHY?? */
    }
    else {
      fCnt--;
      j->status = REGULAR;
      g->regN++;
      oneNodeLayer(g, j);
    }
  }

  if (fCnt == 0) sPop(i, g->freezer);

  g->dMax = r;
  g->aMax = r;

}

/* ================================================================ */

/* perform the excess heuristic test
 
   returns 0 if nothing happened
   1 if j was contracted
*/
static int excessCheck(graph *g, node *j)
{
  node *v;
  int retv = 0;
  static int doing_excess_contraction = 0;

  if (j->excess < g->minCap || j->toContract || 
      g->currentN <= 2 || j == g->sink )
    return 0;
  
  j->toContract = 1;
  sPush(j, g->exStack);
  
  if (doing_excess_contraction == 0)  /* don't allow this while 
					 loop to happen recursively */
    {
      doing_excess_contraction = 1;
      while (!sEmpty(g->exStack)) 
	{
	  sPop(v, g->exStack);
#ifndef NO_PR
	  if (isLeader(v))
#else
	  if (v->status != SOURCE)
#endif
	    {
	      if (g->currentN > 2) 
		{
		  /* node v is now safe to contract into source */
		  
		  g->excessCnt++;
#ifndef NO_PR
		  if ((v->status != FROZEN) && (v != g->sink)) {
		    if ((v->bNext == NULL) && (v->bPrev == NULL)) {
		      excessGap(g, v);
		    }
		  }
		  compactContract(g, g->source, v, MERGE_INTO_W);
#else
		  saturateOutgoing(g, v);
		  if ((v->status != FROZEN) && (v != g->sink)) {
		    if ((v->bNext == NULL) && (v->bPrev == NULL))
		      excessGap(g, v);
		    deleteNode(v, g->buckets + v->d);
		  }
		  if (v->status == REGULAR)
		    g->regN--;
		  g->currentN--;
		  v->status = SOURCE;
		  v->d = g->sourceD;
#endif
		  retv = 1;
		}
	      else break;
	    }
	} 
      doing_excess_contraction = 0;
    }
  
  return(retv);
  
}
#endif



/* ================================================================ */
/* do another max-flow */
static int restartCutComp(graph *g)
{
  node *i, *j;
  int newD, iD;
  node *newSink=NULL;
  bucket *b;

  /* aMax = 0; can have left-over ones from early termination! */

  b = g->buckets + g->sink->d;
  assert(b->first != NULL);
  deleteNode ( g->sink, b );

  assert(g->sink->toContract == 0);

#ifdef NO_PR
  /* saturate all arcs out of the sink */
  /* which will become a source */
  g->sink->toContract = 1;
  saturateOutgoing(g, g->sink);
  g->source = g->sink;
  g->source->status = SOURCE;
  g->source->d = g->sourceD;
  g->regN--;
#else
  if (g->currentN <= 2) return 0;
  compactContract(g,  g->source, g->sink, MERGE_INTO_W); 
  if (g->currentN <= 2) return(0);     /* catch case when excess contracts 
					  bring currentN down to 2 */
#endif

  if ( g->regN > 0 ) {
    newD = g->sinkD;
    for ( b = g->buckets + newD; b->first == NULL; b++ )
      newD++;
    newSink = b->first;
  }
  else {
    /* if no regular nodes, unfreeze the last layer */
    newD = g->infD;
    g->dMax = 0; /* need to compute new value */

    do {
      do {
	if (sEmpty(g->freezer))
	  return ( 0 );
	sPop ( i, g->freezer );
      } while (i == NULL);

      /* now i is not NULL */
      do {
	if (i->status == FROZEN) {

	  if (g->currentN <= 2) return(0);

/* avg: this breaks. Does not seem to help much in any case 
   but it would be nice to have it fixed */
	  /*	 if (!excessCheck(i)) */ {
	    i->status = REGULAR;
	    g->regN++;

	    iD = i->d;

	    b = g->buckets + iD;
	    insertNode ( i, b );
	    if ( iD > g->dMax ) g->dMax = iD;
	    if (( i->excess > 0 )/* && ( iD < sourceD )*/) {
	      insertActive (g,i, b);
	    }

	    if ( iD < newD ) {
	      newD = iD;
	      newSink = i;
	    }

	  }
	}
	sPop ( i, g->freezer );
      } while ( i != NULL );
    } while (g->regN == 0);
  }

  if (g->currentN <= 2) return(0);
  assert(g->regN > 0);

  g->sink = newSink;
  g->sinkD = newD;

  if (g->sink->excess > 0) {
    /* delete sink from active list */
    b = g->buckets + g->sink->d;
    i = FirstRealActive(b);
    if (i == g->sink)
      b->firstActive = NextRealActive(g->sink);
    else {
      while ( 1 ) {
	j = NextRealActive(i);
	if ( j == g->sink ) break;
	else i = j;
      }
      /* now i points to sink */
      i->nextA = NextRealActive(g->sink);
    }
  }

  g->source->d = g->sourceD;

#ifndef NO_PR
  if (g->edgeScanCnt[PR12] * INT_PR_FREQ_12  < g->edgeScanCnt[MAIN]) 
    {
#ifndef DONT_DO_INTERNAL
      sourcePR12(g, g->source);
#endif
      if (g->currentN <= 2) return (0);
    }
  
  if (g->edgeScanCnt[PR34] * INT_PR_FREQ_34  < g->edgeScanCnt[MAIN]) 
    {
#ifndef DONT_DO_INTERNAL
      sourcePR34(g, g->source);
#endif

      if (g->currentN <= 2) return (0);
    }
#endif
  g->context = MAIN;

  return ( 1 );
}


/* ================================================================ */
/* initialize maxflow computation */
static void initCutComp(graph *g)
{
  node *v;
  bucket *b;
#ifdef NO_PR
  arc *a;
#endif
  weight_t vCap;
  weight_t maxV;

  for ( b = g->buckets + g->sourceD; b >= g->buckets; b-- )
    b->first = b->firstActive = NULL;
  
  /* make the largest capacity node initial sink */
  maxV = 0;
  ForAllNodes (g, v) 
    {
#ifdef NO_PR
      vCap = 0;
      ForAllIncidentArcs(v, a )
	vCap += a->resCap;
#else
      vCap = v->cap;
#endif
      if ( vCap > maxV ) 
	{
	  maxV = vCap;
	  g->sink = v;
	}
      v->excess = 0;
      v->status = REGULAR;
      v->current = v->first;
    }

  /* assume graph is compacted */
  g->source = g->sink->first->head;
  g->source->status = SOURCE;

  /* initialize labels */
  g->source->d = g->sourceD;
  g->sink->d = 0;
  g->sinkD = 0;

  insertNode(g->sink, g->buckets);

  b = g->buckets + 1;
  ForAllNodes (g, v)
    if ((v != g->source) && (v != g->sink)) 
      {
	v->d = 1;
	insertNode(v, b);
      }

  g->dMax = 1;
  g->aMax = 0; /* no active nodes at this point */

  saturateOutgoing(g, g->source);
  
#ifdef INIT_UPDATES
  globalUpdate(g);
#endif
}

/* ================================================================ */
/* initialization procedure */

static void mainInit (graph *g)
{
  g->sourceD = 2*g->currentN - 1;
  g->infD = g->sourceD + 1;

#if defined(GLOBAL_UPDATES) || defined(INIT_UPDATES)
  qCreate ( g->currentN, g->BFSqueue );
  qReset ( g->BFSqueue );
#endif

  sCreate ( g->infD, g->freezer );
  sReset ( g->freezer );

  sCreate ( g->currentN, g->regNodes );
  sReset ( g->regNodes );

#ifdef EXCESS_DETECTION
  sCreate ( g->currentN, g->exStack );
  sReset ( g->exStack );
#endif

  memassert(g->buckets = (bucket*)calloc(g->sourceD+2, sizeof(bucket)));
  (g->buckets + g->infD)->first = NULL;
}

/* ================================================================ */

static void mainCleanup(graph *g)
{
#ifdef PR_34
  freeHeap(h);
#endif
  free(g->buckets); 
  qFree(g->BFSqueue); 
  sFree(g->freezer);
  sFree(g->regNodes);

  g->buckets = NULL;
  g->BFSqueue = NULL;
  g->freezer = NULL;
  g->regNodes = NULL;

#ifdef EXCESS_DETECTION
  sFree(g->exStack);
  g->exStack = NULL;
#endif
}


/* ================================================================ */
/* main loop of Hao-Orlin */
void computeCuts(graph *g)
{
  node *v;

  /* Initial graph compaction */

  ForAllNodes(g, v) 
    v->status = REGULAR;

  compact(g);  
  
#ifndef NO_PR
#ifdef PR_34
  makeHeap(h, g->currentN);
#endif
  PRpreprocess(g, 0.75, 0.5, 0.5);
#endif

  mainInit(g);
  g->regN = g->currentN - 1;

  /* initial PR tests done */

  initCutComp(g);
  while (g->currentN > 1) 
    {
      STCut(g);
      
      if (g->sink->excess < g->minCap)   
	saveFlowCut(g);

      if (!restartCutComp(g))
	break;
      assert(g->sink->status == REGULAR);
    }

  mainCleanup(g);
}

/* ================================================================ */

int main (int argc, char *argv[])

{
  float t;
  graph *g;

  g = dimacsParse(stdin);

  printf("c nodes: %14d    arcs: %15d\n", g->n, g->m);

  t = timer();

  computeCuts(g);

  g->dtime = g->dtime - t;
  t = timer() - t;

  printf("c relabels:  %10d    pushes:   %15d\n", g->relabelCnt, g->pushCnt); 

#ifndef NO_PR
  printf("c internal PR: %-6d PR 1: %-6d PR 2: %-6d PR 3: %-6d PR 4: %-6d\n", 
	 g->PR1Cnt+g->PR2Cnt+g->PR3Cnt+g->PR4Cnt, g->PR1Cnt, g->PR2Cnt, 
	 g->PR3Cnt, g->PR4Cnt);
  printf("c PR internal edge scans 1+2: %-8ld 3+4: %-8ld\n", 
	 g->edgeScanCnt[PR12], g->edgeScanCnt[PR34]);
#endif

  printf("c Ave. size: %12.2f\n", (double) g->totalSize / g->STCutCnt);

#ifdef GLOBAL_UPDATES
  printf("c update scans:  %8ld\n", g->edgeScanCnt[GLUPDATE]);
#endif

#ifdef EXCESS_DETECTION
  printf("c Excess contracts:%6d\n", g->excessCnt);
#endif

  printf("c edge scans: %11ld\n", totalEdgeScans(g));
  printf("c s-t Cuts:   %11d\n", g->STCutCnt);
  printf("c one node layers:%7d\n", g->oneCnt);
  printf("c MinCuts discovered:%4d\n",g->cutCount);
#ifndef NO_PR

#endif
  printf("c ttime: %16.2f    capacity: ", t);
  fprintf_wt(stdout, "%11", g->minCap);
  printf("\nc dtime: %16.2f\n", g->dtime);

  printf("\n");

#ifdef SAVECUT
  printCut(g);
#endif

  freeGraph(g);

  return 0;
}

/* ================================================================ */

