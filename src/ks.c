/*----------------------------------------------------------------------
$Id: ks.c,v 1.10 1997/05/16 14:44:07 mslevine Exp $

COMPILATION FLAGS:

NDEBUG compiles out asserts, etc.

FLOATWEIGHTS  (default) Assumes input graph weights are FLOATsSeyon Copyright
INTWEIGHTS Assumes input graph weights are INTs

PR    (default) activates PR (Padberg Rinaldi) heuristics
NO_PR deactivates PR  heuristics.

CUTVAL  (default) Only want value of cut---all PR tests OK
ONECUT  Want cut as well as value.  Limits PR tests
ALLCUTS Compile to identify all (rather than 1) minimum cuts.  Limits
        PR tests.  Presumably involves user implementation of
	FSaveMinCut.  NOT IMPLEMENTED.
APPROXCUTS Compile to identify approximately minimum cuts.  Presumably
        involves user implementation of FSaveMinCut.  Eliminates
        PR tests.  Requires changes in contraction probabilities.
	NOT IMPLEMENTED. 

DEPTHTIMES Accumulates info about how much time spent at each recursion depth
VERBOSE   Define flag VERBOSE for additional output about cuts.
DEBUG   set to various levels, gives increased debugging info (0=none)

----------------------------------------------------------------------*/

#define KS

#ifndef DEBUG
#define DEBUG 0
#endif

/*#define ABSURD .75   /* don't ask */

/* must be defined for ks.h: (not any more -ml)
#if  (!defined(INTWEIGHTS) && !defined(LONGLONGINTWEIGHTS))
#define FLOATWEIGHTS 1
#endif */

#define OFF 0
#define STRICT 1
#define TIGHT 2

#ifdef SAVECUT
#define ONECUT 1
#endif

#if defined(APPROXCUTS)
#define ALLCUTS 0
#define ONECUT 0
#define CUTVAL 0
#elif defined(ALLCUTS)
#define APPROXCUTS 0
#define ONECUT 0
#define CUTVAL 0
#elif defined(ONECUT)
#define APPROXCUTS 0
#define ALLCUTS 0
#define CUTVAL 0
#else /*default CUTVAL*/
#define CUTVAL 1
#define ONECUT 0
#define ALLCUTS 0
#define APPROXCUTS 0
#endif

#if (APPROXCUTS || NO_PR)
#define PR1 OFF
#define PR2 OFF
#elif ALLCUTS
#define PR1 STRICT 
#define PR2 STRICT
#elif ONECUT 
#define PR1 TIGHT
#define PR2 TIGHT
#elif CUTVAL
#define PR1 TIGHT
#define PR2 TIGHT
#endif

#include <stdio.h>
#ifdef BROKEN_HDRS  
extern int printf();
extern int fprintf();
extern int fscanf();
extern int fgetc();
#endif
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <malloc.h>
#include <search.h>
#include <math.h>
#include <values.h>

#ifdef NO_PR
#undef NO_PR
#endif

#include "random.h"
#include "timer.h"
#include "ks.h"

#include "graph.h"
#include "contract.h"
#include "pr.h"
#include "heap.h"

/* System dependecies are such fun. I don't think it matters that this 
dumb substitution isn't as fast a possible. -ml*/
#ifdef linux
#define log2(x) (log(x)*M_LOG2E)
#endif

/* prototypes */

char *safe_malloc(int c);
void InitRoots(metavertex **hmv, int cmv);
void FreeMetagraph(metagraph *pmg);
void FreeMetagraph(metagraph *pmg);
metavertex *pmvSafeFindRoot(metavertex *pmv);
metavertex *pmvFindRoot(metavertex *pmv);
int fUnion(metavertex *pmv1, metavertex *pmv2);
void MergeVertices(metavertex *rgpmv[], int cmv, 
		   metavertex *rgpmvNew[], int cmvNew);
void CountingSort(metaedge rgmeIn[], int imeMax, int icMax, int iWhich,
		  metaedge rgmeOut[]);
void CountingSortVertices(metavertex *rgpmv[], int ipmvMax, int icMax);
metagraph * pmgCompact(metagraph *pmg, int cmvNew);
void ComputeDegrees(metagraph *pmg);
bool FWalkDegrees(metagraph *pmg, 
		  bool (*FNewSmallCut)(metagraph *pmg,int imvOpt,
				       weight_t wt));
bool NewTrivialCut(metagraph *pmg, int imvOpt, weight_t wt);
metagraph *pmgIteratePR(metagraph *pmg, metagraph (*PRAlg(metagraph *pmg, 
							  int cDepth)), 
			double dbl, 
			bool (*FNewSmallCut)(metagraph *pmg, int imvOpt, 
					     weight_t wt),
			int cDepth);
metagraph *pmgPR12Pass(metagraph *pmg, int cDepth);
metagraph *pmgPR34Pass(metagraph *pmg, int cDepth);
metagraph *pmgContract(metagraph *pmg, bool (*FContract)(metaedge *pme), 
		       bool (*FNewSmallCut)(metagraph *pmg, int imvOpt,
					    weight_t wt));
bool fRecursiveContract(metagraph *pmg, bool (*FContract)(metaedge *pme), 
			bool (*FNewSmallCut)(metagraph *pmg, int imvOpt, 
					     weight_t wt), int cDepth);
bool FContract(metaedge *pme);
void WalkMinCuts(metagraph *pmg, double dFailProb,
		 bool (*FNewSmallCut)(metagraph *pmg, int imvOpt, 
				      weight_t wt));
/* metagraph *pmgReadMetagraph(char * filename); */
metagraph *pmgConvertToMetaGraph(graph *g);

/*----------------------------------------------------------------------
           ABOUT THE PADBERG-RINALDI TESTS

We have 2 options for the PR1 tests.
  1) STRICT PR1 tests contract edges bigger than the min-cut bound.
  This can never damage a minimum cut, so we needn't worry about
  storing the cut that yielded the bound: if it was the minimum cut,
  we will encounter it later.  
  2) TIGHT PR1 tests contract edges greater than or _equal_ to the
  upper bound.  This is more dangerous: it can cause us to contract
  the minimum cut (if it is made of one edge).  Therefore, if we plan
  to identify the minimum cut (and not just its value), if we
  encounter a low-degree vertex and use it to update the lower bound,
  we have to make sure to store the cut that it defines since we may
  never see that cut again.

The PR2 tests are also tricky.  
  1) STRICT PR2 test contract edges exceeding half their endpoint's
  degree.  This can eliminate the minimum cut if it is a single
  vertex.  
  2) TIGHT PR2 tests contract edges exactly half one endpoint's
  degree.  If we use tight PR2 tests, we must ensure that every vertex
  is involved in at most one tight PR2 contraction, since otherwise
  our algorithm might never find the minimum cut.

The output objectives have impact on what PR tests are valid:
  1) Minimum cut value.  All PR tests are OK without storing anything,
     (if we contract at most one tight PR2 edge per vertex).
  2) One minimum cut.  All tests are fine iff we store trivial cuts we
     encounter.  If we don't store trivial cuts, we can only use
     strict PR1 tests.
  3) All minimum cuts.  Only "strict" PR tests (ie strictly more than
     mincut for PR1, strictly more than 1/2 degree for PR2) are valid,
     since "tight" PR tests (with equality) might kill various
     mincuts.  We still have to store trivial cuts.
  4) Approximately minimum cuts.  PR tests don't respect these, so
     they shouldn't be done.  For very near minimum cuts, there are
     weaker PR test variants that are OK, but they are not implemented.

----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
ASSORTED NOTES

memory management could be much improved; I did the most 
straightforward malloc/free stuff but since we have good bounds on 
space needed it would probably be faster to manage our own memory, 
but i don't know if it's worth the effort

Variable naming convention may seem quite odd.  it's a scheme called
hungarian; full text available on request, but to summarize what's
used here: 
base typenames:
  integer i, double dbl, count c, function fn, boolean f
  metavertex mv, metaedge me, metagraph mg
  weight wt
prefixes:
  p pointer, h handle (=pp), rg array
sufixes:
  Max: upper limit on value.  not a legal value
  Out: output 
  Tmp: temporary storage
Functions get capitalized.
----------------------------------------------------------------------*/

#define PR12RATE 0.95
#define PR34RATE 0.95
#define DEFAULTSEED (long)3
#define DEFAULTEPSILON 0.05
#define DASH '-'
#define min(a,b)  (((a)<(b))?(a):(b))

/*  GLOBAL VARIABLES */

/*for input*/
metavertex *rgmv;		/* vertices in input graph */
int cmvIn;			/*   how many */
heap h;                         /* for PR preprocessing */

/*for output*/
weight_t  wtMax;			/* best known upper bound on minimum
				   cut.  Might have been discovered 
				   as a small degree that hasn't been
				   visited as a cut yet, so keep
				   separate from wtOpt*/
weight_t  wtOpt;			/* best cut value found so far*/
#ifdef SAVECUT
metavertex **rgpmvOpt;		/* best cut vertices */
int cmvOpt;			/*   how many */
#endif


/*for statistics*/
double tmStart;			/*starting time*/
double tmFound;			/* time mincut was found */
int cTrees;			/*number of trees enumerated*/
int cLeaves;			/*number of leaves enumerated in this tree*/
int cLeavesTot;			/* and all trees */
int cPR12;
int cPR34;
int cScans;                    /* edge scans */


#if (PR2==TIGHT)
int iTime;			/* For distinguishing phases of PR2
				   tests.  Increments each time we
				   start a new series of tests. */
#endif

#ifdef DEPTHTIMES
double depthTimes[200];
int depthCount[200];
int depthNodes[200];
int pr12Count[200];
int pr34Count[200];
#endif

int cDepthMax;                 /*track maximum recursion depth*/

#define Malloc(n,type) ((type *)safe_malloc((n)*sizeof(type)))

#define FREE(x) (free ((char *) (x)))

char *safe_malloc(int c)
    {
    char *pch;
    if ((pch = malloc((unsigned)c)))
	{
	return pch;
	}
    printf("error: out of memory\n");
    fprintf(stderr,"error: out of memory\n");
    exit(17);
    return NULL; /* unreachable, but makes lcc happy */
    }

void FreeMetagraph(metagraph *pmg)
    {
    FREE(pmg->rgpmv);
    FREE(pmg->rgme);
    FREE(pmg);
    }


/*------------------------------------------------------------------------
fRoot
returns true if argument points to root in disjoint sets structure
-----------------------------------------------------------------------*/
#define fRoot(pmv) ((pmv)->pmv == (pmv))

/*----------------------------------------------------------------------
InitRoots

Initializes root pointers for union find, so all sets are singletons

----------------------------------------------------------------------*/

void InitRoots(metavertex **hmv, int cmv)
{
  metavertex *pmv;
  
  while (cmv--) {
    pmv = *(hmv++);
    pmv->pmv = pmv;
    pmv->c = 1;
  }
}


/*----------------------------------------------------------------------
pmgSafeFindRoot

Find operation for union-find (without path compression).  Returns
pointer to root.
----------------------------------------------------------------------*/
metavertex *pmvSafeFindRoot(metavertex *pmv)
    {
    if (fRoot(pmv))
	return pmv;
    else
	return (pmvSafeFindRoot(pmv->pmv));
    }

/*----------------------------------------------------------------------
pmgFindRoot

Find operation for union-find (with path compression).  Returns
pointer to root. 
----------------------------------------------------------------------*/

metavertex *pmvFindRoot(metavertex *pmv)
    {
    if (fRoot(pmv))
	return pmv;
    else
	return (pmv->pmv = pmvFindRoot(pmv->pmv));
    }


/*-------------------------------------------------------------------
fUnion    
Input: two metavertices (sets).

union operation for union-find.  maintains ranks.  

Returns: fFalse iff inputs are already part of same metavertex, 
         otherwise fTrue.
--------------------------------------------------------------------*/

int fUnion(metavertex *pmv1, metavertex *pmv2)
    {
    pmv1 = pmvFindRoot(pmv1);
    pmv2 = pmvFindRoot(pmv2);
    
    if (pmv1==pmv2)
	return fFalse;
    if (pmv1->c > pmv2->c)
	pmv2->pmv = pmv1;
    else 
	{
	pmv1->pmv = pmv2;
	if (pmv1->c == pmv2->c)
	    ++(pmv2->c);
	}
    return fTrue;
    }

/*----------------------------------------------------------------------
MergeVertices 
Input: array of cmv metavertex pointers in rgpmv.
       count cmvNew of number of output vertices.

Results: the cmvNew roots of the current list of vertices rgpmv become the
vertices of the new list of vertices rgpmvNew.  Each vertex in the old
list rgpmv is assigned a number equal to its root's so that CountingSort
can identify them by index.

If PR34, uses global rgc to reorder vertices
----------------------------------------------------------------------*/
    int *rgc = NULL;  /*counts occurences of various verts*/

void MergeVertices(metavertex *rgpmv[], int cmv,
                   metavertex *rgpmvNew[], int cmvNew)
    {
#ifdef PR34
      int i;
#endif
    metavertex *pmv;
    metavertex **hmv;
    metavertex **hmvMax = rgpmv+cmv;

#ifdef PR34
    for (i=0; i<cmvNew; i++)
      rgc[i]=i;
    RandomPermute(rgc,cmvNew);
#endif

    for (hmv = rgpmv; hmv<hmvMax; ++hmv) {
      pmv = (*hmv);
      if (fRoot(pmv))
	    {
#ifdef PR34
	      i = rgc[--cmvNew];
	      pmv->i = i;
	      rgpmvNew[i] = pmv;
#else
	    pmv->i = --cmvNew;
	    rgpmvNew[cmvNew] = pmv;
#endif
	    }
      }

    assert(cmvNew==0);  /*check really were right number of roots*/

    for (hmv = rgpmv; hmv<hmvMax; hmv++)   /*name other vertices by root*/
	(*hmv)->i = pmvFindRoot(*hmv)->i;  
    }

/*----------------------------------------------------------------------
CountingSort

BEFORE CALLING: global rgc needs to hold an array for accumulating counts.

Input: array rgmeIn of imeMax edges pointing to vertices with name in
range 0...icMax-1.  iWhich=0 or 1 specifying an endpoint.

Output:
In array rgmeOut, edges of rgmeIn are reordered accoring to name of
iWhich endpoint.  
Called twice, gives a radix sort of edges by endpoint names.
------------------------------------------------------------------------*/

void CountingSort(metaedge rgmeIn[], int imeMax, int icMax, int iWhich,
                  metaedge rgmeOut[])
    {
    int c;
    int *pc;
    int *pcMax=rgc+icMax;
    metaedge *pme;
    metaedge *pmeMax = rgmeIn+imeMax;

    assert(rgc != NULL);
    
    /*count occurences of various verts*/
    for (pc=rgc; pc<pcMax; pc++)
	*pc = 0;
    for (pme=rgmeIn; pme<pmeMax; pme++)
	++rgc[pme->rgpmv[iWhich]->i];

    /*now prefix sums tell where in rgmeOut the edges for each */
    /*particular vertex name should begin*/
    for (pc=rgc, c=0; pc < pcMax; pc++) /*set rgc[ic] to total number of */
	{			        /*edges with indices <= ic*/
	c += *pc;           
	*pc = c-(*pc);     /*old c*/
	}

    for (pme = rgmeIn; pme < pmeMax; pme++)
	rgmeOut[rgc[pme->rgpmv[iWhich]->i]++] = *pme;
    }

	      
#ifdef 0

void CountingSortVertices(metavertex *rgpmv[], int ipmvMax, int icMax)
    {
    int c;
    int *pc;
    int *pcMax=rgc+icMax;
    metavertex **hmv;
    metavertex **hmvMax = rgpmv+ipmvMax;

    assert(rgc != NULL);
    
    /*count occurences of various scores*/
    for (pc=rgc; pc<pcMax; pc++)
	*pc = 0;
    for (hmv=rgpmv; hmv<hmvMax; hmv++)
	++rgc[(*hmv)->iScore];

    /*now prefix sums tell how many verts have a given score or less*/
    for (pc=rgc, c=0; pc < pcMax; pc++) /*set rgc[ic] to total number of */
	{			        /*verts with soores <= ic*/
	c += *pc;           
	*pc = c-(*pc);     /*old c*/
	}

    for (hmv = rgpmv; hmv < hmvMax; hmv++)
	(*hmv)->i = rgc[(*hmv)->iScore]++;
    }
	      
#endif	  

/*----------------------------------------------------------------------
pmgCompact

BEFORE CALLING: setup rgmeTmp to hold a list of edges

Input: graph pmg with some vertices merged, cmvNew counting number of
merged groups

Output: A new graph, whose metavertices are the roots of the various
sets of metavertices in the old graph, and with a new edge list
pointing only to the new metavertices.

Removes contracted edges, and merges edge with the same endpoints.
Merger is carried out by using two calls to CountingSort to radix-sort
the edges; then duplicates are adjacent and can easily be merged.
Nondestructive of original graph.
----------------------------------------------------------------------*/
metaedge *rgmeTmp = NULL;

metagraph * pmgCompact(metagraph *pmg,
		       int cmvNew)  /*number of verts in compact graph*/
    {
    metaedge *pme,*pmeOut;
    int cme = pmg->cme;
    metaedge *rgme = pmg->rgme;
    metaedge *rgmeOut;
    metavertex **rgpmvOut = Malloc(cmvNew,metavertex *);
    metagraph *pmgOut = Malloc(1,metagraph);
    metavertex *pmv0, *pmv1;

    assert(rgmeTmp != NULL);

    MergeVertices(pmg->rgpmv,pmg->cmv,rgpmvOut,cmvNew);
    /*now each metavertex points directly at its root*/
    /*and roots are named 0,1,...,cmvNew-1*/

    /*collect non-contracted edges, and dereference their endpoints*/
    /*put low ordered vertex first in copied edges, so we don't */
    /*keep two distinct edges for each pair of endpoints*/
    for (pme=rgme+cme-1, pmeOut = rgmeTmp; pme>=rgme; pme--)
	{
	if ((pmv0=pme->rgpmv[0]->pmv) != (pmv1=pme->rgpmv[1]->pmv))
	    {
	    if (pmv0 < pmv1)
		{
		pmeOut->rgpmv[0] = pmv0;
		pmeOut->rgpmv[1] = pmv1;
		}
	    else
		{
		pmeOut->rgpmv[0] = pmv1;
		pmeOut->rgpmv[1] = pmv0;
		}
	    pmeOut->wt = pme->wt;
	    ++pmeOut;
	    }
	}
    cme = pmeOut-rgmeTmp;
    rgmeOut = Malloc(cme,metaedge);

    CountingSort(rgmeTmp,cme,cmvNew,1,rgmeOut);
    CountingSort(rgmeOut,cme,cmvNew,0,rgmeTmp);

    /*merge same-endpoint edges*/
    pmeOut = rgmeOut;
    rgmeOut[0] = rgmeTmp[cme-1];
    for (pme = rgmeTmp+cme-2; pme >= rgmeTmp; pme--)
	if ((pmeOut->rgpmv[0] == pme->rgpmv[0]) &&
	    (pmeOut->rgpmv[1] == pme->rgpmv[1])) 
	    {
	    pmeOut->wt +=pme->wt;
	    }
	else
	    *++pmeOut = *pme;

    pmgOut->rgme = rgmeOut;
    pmgOut->cme = pmeOut-rgmeOut+1;
    pmgOut->rgpmv = rgpmvOut;
    pmgOut->cmv = cmvNew;

    return pmgOut;
    }

/* -------------------------------------------------------
ComputeDegrees

Store degree of each vertex in its "wt" field.
----------------------------------------- */
void ComputeDegrees(metagraph *pmg)
    {
    metavertex **rgpmv = pmg->rgpmv;
    metaedge *rgme = pmg->rgme;
    int cme = pmg->cme;
    int cmv = pmg->cmv;
    int i;

    for (i=0; i < cmv; i++)
      rgpmv[i]->wt = 0;

    for (i=0; i<cme; i++)
      {
	rgme[i].rgpmv[0]->wt += rgme[i].wt;
	rgme[i].rgpmv[1]->wt += rgme[i].wt;
      }
    }


bool FWalkDegrees( metagraph *pmg, 
		   bool (*FNewSmallCut)(metagraph *pmg, int imvOpt, 
					weight_t wt))
{
  metavertex **rgpmv = pmg->rgpmv;
    int cmv = pmg->cmv;
    int i;
    weight_t wt=wtMax;
    int imv=0;
    bool f=0; /*default retval*/

    ComputeDegrees(pmg);

    for (i=0; i < cmv; i++)
	{
 	if (wt > rgpmv[i]->wt)
	  {
	    wt = rgpmv[i]->wt;
            imv = i;
	  }
	}
    
    if (wt < wtMax)
      {
      f=(*FNewSmallCut)(pmg,imv,wt);
      wtMax=wt;
      }

    return f;
}

     

/*----------------------------------------------------------------------
NewTrivialCut

This takes a graph with a single vertex identified by imv whose weight
(degree) may be smaller than any previously encountered cut.  This is
checked; if it is a new minimum, the cut is saved as the new best
cut.  This only needs to be done if we want the cut itself and not
just its value.

This is a client function.  it can be replaced by anything the user
wants without affecting the performance of the mincut algorithm.  That
is, the rest of the code will walk through all the minimum cuts and is
guaranteed to call this function on the minimum cut (and maybe
others).

----------------------------------------------------------------------*/

bool NewTrivialCut(metagraph *pmg, int imvOpt, weight_t wt)
    {    
#ifdef SAVECUT
    int i,c;
    bool fSaveNonTrivialSide;
    metavertex *pmvOpt = pmvSafeFindRoot(pmg->rgpmv+imvOpt);
#endif

    if (wt >= wtOpt)  /* not a better cut */
	return fFalse;  /*keep going*/
    
    wtOpt = wt;
    

    tmFound = timer();
#if (DEBUG >= 2)
    printf("leaf %5d %7d %8d %12lf %12f\n",
	   cTrees,cLeaves,cLeavesTot, tmFound-tmStart, (double)wtOpt);
#endif
        
#ifdef SAVECUT

    /*arrange to save 'smaller' side of cut*/
    /*set fSaveNonTrivialSide to 0 or 1 according to prevalence of
      trivial- or non-trivial- side verts.*/
    c=0;
    for (i=0; i<cmvIn; i++)
	{
	if (pmvSafeFindRoot(rgmv+i) == pmvOpt)
	    ++c;
	}
    fSaveNonTrivialSide = ((2*c)>cmvIn);
    
    /*collect vertices on smaller side*/
    cmvOpt=0;

    for (i=0; i<cmvIn; i++)
	if ((pmvSafeFindRoot(rgmv+i) == pmvOpt) ^ fSaveNonTrivialSide)
	    rgpmvOpt[cmvOpt++] = &(rgmv[i]);
#endif

    return fFalse;  /*signal to keep going*/
  }

/*----------------------------------------------------------------------
pmgIteratePR

Input should be a graph with more than 2 vertices.

Takes a PR test and iterates it until it stops reducing vertices by
more than the specified threshold dbl.

A bit messy because we want to discard all intermediate graphs but not
the input graph.

----------------------------------------------------------------------*/

metagraph *pmgIteratePR(metagraph *pmg, metagraph (*PRAlg(metagraph *pmg, 
							  int cDepth)), 
			double dbl, 
			bool (*FNewSmallCut)(metagraph *pmg, int imvOpt, 
					     weight_t wt),
			int cDepth)
{
  bool fAgain=1;
  metagraph *pmgOld, *pmgOrig=pmg;
  
  if (pmg->cmv <= 2)
    return pmgCompact(pmg,pmg->cmv);/* a copy, for consistent malloc/free */
  
  do {
      
    if (FWalkDegrees(pmg,FNewSmallCut)) {
      if (pmg != pmgOrig)
	FreeMetagraph(pmg);
      return NULL; /*user terminated*/
    }
    
    InitRoots(pmg->rgpmv,pmg->cmv);

    pmgOld = pmg;
    pmg = (*PRAlg)(pmg, cDepth);
    fAgain = (pmg->cmv < dbl * pmgOld->cmv);
    if (pmgOld != pmgOrig)
      FreeMetagraph(pmgOld);  /*dump intermediate versions*/
  }
  while (fAgain && (pmg->cmv >2));
  
  return pmg;
}


/*----------------------------------------------------------------------
pmgPR12Pass

Takes a graph (sorted edge list, no contractions yet) and performs PR12 tests.

Assumes pmv->wt degrees and pmv->pmv roots have been set up.
----------------------------------------------------------------------*/

#if (PR1 || PR2)

metagraph *pmgPR12Pass(metagraph *pmg, int cDepth)
{
  metaedge *pme;
  int cme = pmg->cme;
  int cmv = pmg->cmv;
  bool fDoit;
  weight_t wt;

#if (PR2==TIGHT)
  ++iTime;
#endif

  for (pme = pmg->rgme; cme--; ++pme) {
    wt = pme->wt;
    cScans++;

#if (DEBUG >= 3)    
    printf("consider %d %d: ", pme->rgpmv[0]->name,
	   pme->rgpmv[1]->name);
#endif    
    
    fDoit = 
      (
#if (PR1==TIGHT)
       (wt >= wtMax) ||
#elif (PR1==STRICT)
       (wt > wtMax) ||
#endif

#if ((PR2==STRICT) || PR2==TIGHT)
       ((pme->rgpmv[0]->wt < 2*wt) || (pme->rgpmv[1]->wt < 2*wt)) ||
#if (PR2==TIGHT)
       ((pme->rgpmv[0]->wt == 2*wt) &&
	(pme->rgpmv[0]->iTime != iTime) && /*no prior PR this time around*/
	(pme->rgpmv[0]->iTime = iTime)) /*NOT A BUG!
					  Mark doing PR this time around*/
       ||
       ((pme->rgpmv[1]->wt == 2*wt) &&
	(pme->rgpmv[1]->iTime != iTime) && /*no prior PR this time around*/
	(pme->rgpmv[1]->iTime = iTime)) /*mark doing PR this time around*/
       );
#endif
#endif

#if (DEBUG >= 3)    
    printf("%d \n",fDoit);
#endif    

    if (fDoit) {
      if (fUnion(pme->rgpmv[0],pme->rgpmv[1]))
	{
	  ++cPR12;
#ifdef DEPTHTIMES
	  ++pr12Count[cDepth];
#endif	  
	  --cmv;
	  if (cmv==2) {
	    ++cLeavesTot;
	    ++cLeaves;
	    break;
	  }
	}
    }
  }

  return pmgCompact(pmg,cmv);
}

#endif

#ifdef PR34
/*---------------------------------------------------------------
pmgPR34Pass

Takes a graph (sorted edge list, no contractions yet) and performs
PR34 tests.

Assumes pmv->wt degrees and pmv->pmv roots set up.

Doesn't (yet) distinguish between STRICT and TIGHT options for PR34
tests.  
----------------------------------------------------------------------*/

metagraph * pmgPR34Pass(metagraph *pmg, int cDepth)
{
  metavertex **rgpmv = pmg->rgpmv;
  metaedge *rgme = pmg->rgme;
  int cmv = pmg->cmv;
  int cmvOut = cmv;
  int cme = pmg->cme;
  metavertex *pmv, *pmv1, *pmv2, *pmvOld;
  metaedge *pme, *pme1, *pmeMax, *pme1Max;
  metavertex **hmv;
  metavertex **hmvMax = rgpmv+cmv;
  bool fDoit;
  weight_t wt;

  iTime++;
  
  for (hmv = rgpmv; hmv<hmvMax; hmv++) {
    pmv = *hmv;
    pmv->pmv34 = NULL;
    /*set default edge list to be empty */
    pmv->pmeFirst = NULL;
    ++pmv->pmeFirst;
    pmv->pmeMax = NULL;
  }

  /* Set edge list boundaries for all vertices */
  pmvOld = NULL;
  for (pme=rgme, pmeMax = rgme+pmg->cme; pme<pmeMax; pme++)    {
    pmv = pme->rgpmv[0];
    cScans++;
    if (pmv != pmvOld) {      /* new vertex*/
      if (pmvOld != NULL) 
	pmvOld->pmeMax = pme;
      pmv->pmeFirst = pme;    /* head of adjacency list for v*/
      pmvOld = pmv;
    }
  }
  pmvOld->pmeMax = rgme+cme; /*last vertex's edge list never ran out
			       in loop*/
    
  /* scan vertices one at a time.  a vertex can be scanned as the
     distance-one neighbor of at most one vertex */
  for (hmv = rgpmv; hmv < hmvMax; hmv++) {
    pmv = *hmv;
    /* initialize my neighbors*/
    for (pme = pmv->pmeFirst, pmeMax = pmv->pmeMax; pme < pmeMax; pme++) {
      pmv1 = pme->rgpmv[1];
      pmv1->wt34 = pmv1->wt34Tot = pme->wt;
      pmv1->pmv34 = pmv; 
    }
    /* apply edge test to my neighbors */
    for (pme = pmv->pmeFirst, pmeMax = pmv->pmeMax; pme < pmeMax; pme++) { 
      pmv1 = pme->rgpmv[1];
      wt = pme->wt;
      cScans++;

      /*PR12 tests are cheap.  squeeze them in*/
      fDoit = 
	( 
#if (PR1==TIGHT)
	 (wt >= wtMax) ||
#elif (PR1==STRICT)
	 (wt > wtMax) ||
#endif

#if ((PR2==STRICT) || PR2==TIGHT)
	 ((pmv->wt < 2*wt) || (pmv1->wt < 2*wt)) ||
#if (PR2==TIGHT)
	 ((pmv->wt == 2*wt) &&
	  (pmv->iTime != iTime) && /*no prior PR this time around*/
	  (pmv->iTime = iTime)) /*NOT A BUG!
				  Mark doing PR this time around*/
	 ||
	 ((pmv1->wt == 2*wt) &&
	  (pmv1->iTime != iTime) && /*no prior PR this time around*/
	  (pmv1->iTime = iTime)) /*mark doing PR this time around*/
	 );
#endif
#endif
      if (fDoit) {
	if (fUnion(pmv, pmv1)) {
#ifdef DEPTHTIMES
	  ++pr12Count[iDepth];
#endif	  
	 ++cPR12;
	 if  (--cmvOut ==2) {
	  ++cLeaves;
	  ++cLeavesTot;
	  goto exit;
	 }
	}
      }
      
      /*PR34 tests*/
      else if (!pmv1->fScanned) {
	pmv1->fScanned = fTrue;
	for (pme1 = pmv1->pmeFirst, pme1Max = pmv1->pmeMax; 
	     pme1 < pme1Max; pme1++) {
	  pmv2 = pme1->rgpmv[1];
	  wt = min(pme->wt,pme1->wt);
	  if (pmv2->pmv34 != pmv) { /* first visit from v */
	    pmv2->pmv34 = pmv;
	    pmv2->wt34 = 0;
	    pmv2->wt34Tot = wt;
	    /*since PR12 tests were already applied, it is impossible
	      for a PR34 test to succeed at this point, so don't bother*/
	  }
	  else {
	    pmv2->wt34Tot += wt;
	    wt = 2*(wt+pmv2->wt34); /*twice triangle weight for PR3*/
	    fDoit = (
		     (pmv2->wt34Tot >= wtMax) || /*PR4*/
		     (pmv->wt < wt) || 
		     (pmv2->wt < wt) ||
		     ((pmv->wt == wt) &&
		      (pmv->iTime != iTime) && /*no prior PR this time*/
		      (pmv->iTime = iTime))    /*Mark doing PR this time*/
		     ||
		     ((pmv2->wt == wt) &&
		      (pmv2->iTime != iTime) && /*no prev PR this time*/
		      (pmv2->iTime = iTime))    /*mark doing PR this time*/
		     );
	    if (fDoit) {
	      if (fUnion(pmv, pmv1) && ++cPR34 && (--cmvOut ==2)) {
		++cLeaves;
		++cLeavesTot;
		goto exit;
	      }
	    }
	    if (pmv2->wt34 != 0) {  /* v2 neighbors v.  accumulate to v1*/
	      wt =min(pme1->wt,pmv2->wt34);
	      pmv1->wt34Tot += wt;
	      wt = 2*(pme->wt + wt);  /*triangle weight for PR3*/
	      fDoit = (
		       (pmv1->wt34Tot >= wtMax) || /*PR4*/
		       (pmv->wt < wt) || /*PR3*/
		       (pmv1->wt < wt) ||
		       ((pmv->wt == wt) &&
			(pmv->iTime != iTime) && /*no prior PR this time*/
			(pmv->iTime = iTime))    /*Mark doing PR this time*/
		       ||
		       ((pmv1->wt == 2*wt) &&
			(pmv1->iTime != iTime) && /*no prev PR this time*/
			(pmv1->iTime = iTime))    /*mark doing PR this time*/
		       );
	      if (fDoit) {
		if (fUnion(pmv, pmv1)) {
		  ++cPR34;
#ifdef DEPTHTIMES
		  ++pr34Count[iDepth];
#endif	  
		  if (--cmvOut ==2) {
		    ++cLeaves;
		    ++cLeavesTot;
		    goto exit;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
      
exit:
  return pmgCompact(pmg,cmvOut);
}

#endif 

/*----------------------------------------------------------------------
pmgContract

BEFORE: Malloc rgwt to size cmeIn

Input: metagraph pmg,
       boolean function FContract, called to see whether each edge should
          be contracted.

Basic contraction algorithm of Karger:SODA93.

Output: New graph.  NULL if user requested termination.
-----------------------------------------------------------------------*/
weight_t *rgwt = NULL;

metagraph *pmgContract(metagraph *pmg, bool (*FContract)(metaedge *pme), 
		       bool (*FNewSmallCut)(metagraph *pmg, int imvOpt, 
					    weight_t wt))
    {
    metavertex **rgpmv = pmg->rgpmv;
    metaedge *rgme = pmg->rgme;
    int i;
    int cme = pmg->cme;
    int cmv = pmg->cmv;
    
#if (PR1 || PR2)
    if (FWalkDegrees(pmg,FNewSmallCut))
      return NULL;
#endif    

    InitRoots(rgpmv,cmv);
    for (i=0; i<cme; i++)
	{
	  cScans++;
	  if((*FContract)(rgme+i))
	    {
	    if (fUnion(rgme[i].rgpmv[0],rgme[i].rgpmv[1]))
		{
		--cmv;
		if (cmv==2)
		    break;
		}
	    }
	}

    return pmgCompact(pmg,cmv);

    }

/*----------------------------------------------------------------------
fRecursiveContract

Recursive Contraction Algorithm of Karger/Stein:STOC93

input: graph to contract, contraction function (decides edges
to contract to in each step).

Output: True if should stop generating cuts.
----------------------------------------------------------------------*/

bool fRecursiveContract(metagraph *pmg, bool (*FContract)(metaedge *pme), 
			bool (*FNewSmallCut)(metagraph *pmg, int imvOpt, 
					     weight_t wt), 
			int cDepth)
{
  metagraph *pmgNew;
  bool f;
#ifdef DEPTHTIMES
  double stime=timer();
#endif

  assert(pmg->cmv >= 2);
  assert(pmg->cme > 0);         /*if fails, graph is disconnected*/

  if (cDepth > cDepthMax)
    cDepthMax=cDepth;
  if (pmg->cmv == 2) 
    {	
      cLeaves++;
      cLeavesTot++;
      if (pmg->rgme[0].wt < wtMax)
	wtMax=pmg->rgme[0].wt;
      return (*FNewSmallCut)(pmg,0,pmg->rgme[0].wt);
    }			
  else {

    if ((pmgNew = pmgContract(pmg,FContract,FNewSmallCut)) == NULL)
      return fTrue;

#if (PR1 || PR2)
    pmg = pmgNew;
    pmgNew = pmgIteratePR(pmg,pmgPR12Pass,PR12RATE,FNewSmallCut,cDepth);
    FreeMetagraph(pmg);
    if (pmgNew ==NULL) {
      return fTrue;
    }
#endif    
#ifdef PR34    
    pmg = pmgNew;
    pmgNew = pmgIteratePR(pmg,pmgPR34Pass,PR34RATE,FNewSmallCut,cDepth);
    FreeMetagraph(pmg);
    if (pmgNew ==NULL) {
      return fTrue;
    }
#endif

#ifdef DEPTHTIMES
    /*    ++depthCount[(int)(log((float)pmgNew->cmv)/log(2.0))];
    depthTimes[(int)(log((float)pmgNew->cmv)/log(2.0))] += timer()-stime; */
    ++depthCount[cDepth];
    depthTimes[cDepth] += timer()-stime;
    depthNodes[cDepth] += pmg->cmv;
#endif
    f = fRecursiveContract(pmgNew,FContract,FNewSmallCut,cDepth+1)
      || fRecursiveContract(pmgNew,FContract,FNewSmallCut,cDepth+1);
    FreeMetagraph(pmgNew);
    return f;
  }
}

/*----------------------------------------------------------------------
FContract

Function deciding whether to contract an edge.  Returns true if a PR
test applies or if a coin flip says to. 
----------------------------------------------------------------------*/

bool FContract(metaedge *pme)
    {
    double wt=pme->wt;
    bool fDoit;
    
#if (DEBUG >= 3)    
    printf("consider %d %d: ", pme->rgpmv[0]->name,
	   pme->rgpmv[1]->name);
#endif    

    fDoit = (
#if (PR1==TIGHT)
      (wt >= wtMax) ||
#elif (PR1==STRICT)
      (wt > wtMax) ||
#endif

#if ((PR2==STRICT) || PR2==TIGHT)
      ((pme->rgpmv[0]->wt < 2*wt) || (pme->rgpmv[1]->wt < 2*wt)) ||
#if (PR2==TIGHT)
      ((pme->rgpmv[0]->wt == 2*wt) &&
       (pme->rgpmv[0]->iTime != iTime) && /*no prior PR this time around*/
       (pme->rgpmv[0]->iTime = iTime)) /*NOT A BUG!
					 Mark doing PR this time around*/
      ||
      ((pme->rgpmv[1]->wt == 2*wt) &&
       (pme->rgpmv[1]->iTime != iTime) && /*no prior PR this time around*/
       (pme->rgpmv[1]->iTime = iTime)) /*mark doing PR this time around*/
      ||
#endif
#endif
#ifdef ABSURD
      (DblUnitRandom() > pow(2.0,(double) -wt*ABSURD/ wtMax)));
#else
      (DblUnitRandom() > pow(2.0,(double) -wt/ wtMax)));
#endif

#if (DEBUG >= 3)    
    printf("%d \n",fDoit);
#endif    

    return fDoit;
    
    }
  
/*----------------------------------------------------------------------
WalkMinCuts

Enumerates cuts so that the probability of failing to encounter any
particular minimum cut is at most dFailProb.  Calls FNewCut on each cut it
encounters.  Aborts enumeration if FNewCut returns fTrue on some cut.
----------------------------------------------------------------------*/

void WalkMinCuts(metagraph *pmg, double dFailProb, 
		 bool (*FNewSmallCut)(metagraph *pmg, int imvOpt, 
				      weight_t wt))
    {
    int c;
    double dblTreeFindProb, dblTreeFindProbSav, dbllgn;
    int iTreeDepth, iTreeDepthSav=0;
    int cIter;

    /*setup temporary storage*/
    rgc = Malloc(pmg->cmv,int);	/*for CountingSort*/
    rgmeTmp = Malloc(pmg->cme,metaedge); /*for pmgCompact*/
    rgwt = Malloc(pmg->cme+1,weight_t); /*for pmgContract*/

    /*ComputeDegrees(pmg); computed already by deleteExtras() (after input)*/ 
    wtMax= 1+pmg->rgpmv[0]->wt;
    if (FWalkDegrees(pmg,FNewSmallCut))
      goto exit;
    
    cPR12=cPR34=0;

    if (pmg->cmv == 2) {
      cIter = 1;
      iTreeDepth = 0;
    }
    else {
      dbllgn = log2(pmg->cmv);
      dblTreeFindProbSav = 0;
      for (iTreeDepth = ceil(2*dbllgn); 
	   iTreeDepth < 5*dbllgn; iTreeDepth++)	{
	/*	iTreeDepth = 2*dbllgn+2*log2(dbllgn);*/
	dblTreeFindProb = 2/(5/2.0+iTreeDepth+log(2+iTreeDepth))*
	  (1-pmg->cmv*pow(2.0,-iTreeDepth/2.0));
	if (dblTreeFindProb > dblTreeFindProbSav) {
	  iTreeDepthSav = iTreeDepth;
	  dblTreeFindProbSav = dblTreeFindProb;
	}
      }
      dblTreeFindProb = dblTreeFindProbSav;
      iTreeDepth = iTreeDepthSav;
#ifdef ABSURD      
      dblTreeFindProb = pow(2.0,ABSURD+1)-pow(2.0, 2*ABSURD);
#endif
      cIter = ceil(log(dFailProb) / log(1.0 - dblTreeFindProb));
      printf("c tree success prob. %f\n",dblTreeFindProb);
    }

#if (DEBUG >= 1)
    fprintf(stdout,"c %d depth %d trees for success probability %f\n",
	    cIter,(int)iTreeDepth,1-dFailProb);
    cDepthMax=0;
#endif
 
#if (DEBUG >= 2)
    /* header to track encountered mincuts */
    printf("     %5s %7s %8s %12s %12s\n",
	   "tree","leaf","total","time (sec)","value");
#endif

    cLeavesTot = 0;
    for (c = 0; c<cIter; c++)
	{
	cTrees++;  
	cLeaves=0;
	if (fRecursiveContract(pmg,FContract,FNewSmallCut,0))
	    {
	    break;
	    }
	}
#if (DEBUG >= 1)
    if (cDepthMax > 0) {
      dblTreeFindProb = 2.0 / ( 5/2.0 + cDepthMax + log(2+cDepthMax));
      cIter = ceil(log(dFailProb) / log(1.0 - dblTreeFindProb));
    }
    else {
      cIter = 1;
    }
    fprintf(stdout,"c Actually sufficient: %d trees, depth %d\n",
	    cIter,(int)cDepthMax);
#endif
    exit:
#if (DEBUG >= 1)
    printf("c internal PR12:   %7d (%6d)   PR34:     %6d (%5d)\n", 
	   cPR12, (int)(cPR12/(float)cTrees), 
	   cPR34, (int)(cPR34/(float)cTrees)); 
#endif
    FREE(pmg);
    FREE(rgc);
    FREE(rgmeTmp);
    FREE(rgwt);
    }

#ifdef 0  /* we now use Andrew's parse() */
/*------------------------------------------------------------------------
  pmgReadMetagraph 

input: file name
output: pointer to a metagraph representing the input graph

this function takes a file name, representing a file in the DIMACS format
and reads the data into a metagraph data structure.  

The input is not robust.  There can be one (or no) comment line,
followed by p cut num_verts num_edges and the rest must be the "a"
lines.  Andrew's code catches 23 different types of errors in the
input file.  I will write a more robust input routine soon, it didn't
seem that important for now.

NOTE:  We translate vertex names by subtracting 

-------------------------------------------------------------------*/
metagraph *pmgReadMetagraph(char * filename)
{
    FILE *f;
    metagraph *pmg = Malloc(1,metagraph);
    metavertex **rgpmv;
    metaedge *rgme;

    int i;
    int end1, end2;
    char intype;
    char probtype[255];
    int cmeIn;
    weight_t wt;

    if (strcmp (filename, "-")) {
	f = fopen (filename, "r");
	if (!f) {
	    perror (filename);
	    fprintf (stderr, "Unable to open %s for input\n",filename);
	    exit (1);
	}
    } else {
	f = stdin;
    }

    /* QUICK AND DIRTY INPUT FORMAT.  Allows one optional comment line,
       followed by a line  beginning with p, followed by lines beginning
       with a.
       SHOULD FIX */

    /* Allow Line 1 to have an optional comment */
    if (fscanf(f,"%c",&intype) != 1)
	{
	
	fprintf(stderr,"%s: First line is blank, exiting\n",filename);
	exit(1);
	}
    
    if (intype == 'c')
	{
	while (fgetc(f) != '\n');  /* skip rest of line */
	
	if (fscanf (f, "%c %s %d %d",&intype,probtype,&cmvIn, &cmeIn) != 4) {
	fprintf (stderr, "%s: Second line has wrong number of args. \n",filename);
	exit (1);
	}
	
       if (intype != 'p') 
	{
	fprintf (stderr, "%s: Second line is not of right form \n",filename);
	exit (1);
	}
	}
    
    
    else if (intype == 'p') 
	{
	if (fscanf (f, "%s %d %d",probtype,&cmvIn, &cmeIn) != 3) {
	fprintf (stderr, "%s: First line is not of right form \n",filename);
	exit (1);
	}
	
	}
    else 
	{
	fprintf (stderr, "%s: First line is not of right form \n",filename);
	exit (1);
	}
    
	

    rgpmv = Malloc(cmvIn, metavertex *);
    rgmv = Malloc(cmvIn, metavertex );
    rgme = Malloc(cmeIn, metaedge);
    pmg->rgme = rgme;
    pmg->rgpmv = rgpmv;
    pmg->cme = cmeIn;
    pmg->cmv = cmvIn;

    /* initialize node list */
    for (i=0; i < cmvIn; i++)
	{
	rgmv[i].pmv = &(rgmv[i]);
	rgpmv[i] = &(rgmv[i]);
	rgmv[i].name = i+1;
	}
    
    /* read edges and set up metaedge data structure */

    for (i=0;i < cmeIn; i++)
	{
	if (fscanf(f,"%s %d %d %" WT_RD_FORMAT,&intype,&end1,&end2,&wt) != 4)
	    {
	    fprintf (stderr, "%s: incorrect format for edge file in line %d\n",
		     filename,i+1); 
	    exit (1);
	    }
	if (intype != 'a')
	    {
	    fprintf (stderr, "%s: incorrect format for edge file in line %d\n",
		     filename,cmeIn+1); 
	    exit (1);
	    }
	rgme[i].rgpmv[0] = &(rgmv[end1-1]);
	rgme[i].rgpmv[1] = &(rgmv[end2-1]);
	rgme[i].wt = wt;
	}
    return pmg;
    }
    
#endif /* 0 */

/*------------------------------------------------------------------------
  pmgConvertToMetagraph 

input: info in global vars (nodes, currentN)
output: pointer to a metagraph representing the input graph

NOTE:  We translate vertex names by subtracting 

-------------------------------------------------------------------*/
metagraph *pmgConvertToMetaGraph(graph *g)
{
  float t = timer();
  metagraph *pmg = Malloc(1,metagraph);
  metavertex **rgpmv;
  metaedge *rgme;
  
  int i;
  int cmeIn=0;
  node *v;
  arc *a;
  
  cmvIn = g->currentN;
  
  rgpmv = Malloc(cmvIn, metavertex *);
  rgmv = Malloc(cmvIn, metavertex );
  
  /* count edges, initialize node list */
  i = 0;
  ForAllNodes(g, v)
    {
      rgmv[i].pmv = &(rgmv[i]);
      rgpmv[i] = &(rgmv[i]);
      rgmv[i].name = i+1;
      rgmv[i].wt = v->cap;
      ForAllIncidentArcs(v,a) cmeIn++;
      i++;
    }
  cmeIn = cmeIn>>1;
  
  rgme = Malloc(cmeIn, metaedge);
  pmg->rgme = rgme;
  pmg->rgpmv = rgpmv;
  pmg->cme = cmeIn;
  pmg->cmv = cmvIn;
  
  /* set up metaedge data structure */
  i = 0;
  ForAllNodes(g, v)
    ForAllIncidentArcs(v,a)
    if (INDEX(v) < INDEX(a->head))
      {
	rgme[i].rgpmv[0] = &(rgmv[INDEX(v)]);
	rgme[i].rgpmv[1] = &(rgmv[INDEX(a->head)]);
	rgme[i].wt = a->cap;
	i++;
      }
  
  printf("c data structure conversion in %.2f time\n", timer()-t);
  return pmg;
}

int main(int argc, char ** argv)
{
    double dFailProb;
    char *filename;
    char args[50];
    metagraph *pmg;
    int i;
    int firstoptional;
    graph *g;  /* graph type from graph.h; for input and PRpreprocess only */

#ifdef DEPTHTIMES
    int iDepthMax=0;

    for (i=0; i<200; i++)
      depthTimes[i] = pr12Count[i] = pr34Count[i] = 0;
#endif
    cScans = 0;
    /* set up defaults in case not specified in the input */
    SetRandom(DEFAULTSEED);
    dFailProb = DEFAULTEPSILON;
    
    if (argc > 4) 
	{	
	  fprintf(stderr, "Usage: %s graph -e epsilon -s seed \n",argv[0]);
	  exit(1);
	}

    if (argc == 1)
      filename = "-";
    else
      {
    
      strcpy(args,argv[1]); /* check first argument */

    if ( args[0] != DASH) /* we have a filename */
    {
      filename = argv[1];
      firstoptional = 2;    
    }
    else
      {
      filename = "-";
      firstoptional = 1;
    };

    for (i=firstoptional; i<argc;i++)  /* now must have dash */
      {
	strcpy(args,argv[i]);
	if (args[0] != DASH) 
	  {printf("optional arguments must start with -\n");
	   exit(4);
	 }

	switch (args[1]) {
	case 's' : {
	  SetRandom(atol(&args[2]));
	  break;
	}

	case 'e' : {
	  dFailProb = atof(&args[2]);
	  break;
	}
    
	default:{
	  printf("bad input option \n");
	  exit(5);
	  break;
	}
	}
      }
}    
/* end of input */
    tmStart = timer();  /* only to time input reading */
    
    /* read input */
    /* WARNING: we use the same data-structures and code as the other 
       algorithms to read in the input and do the PR preprocessing. This 
       means that two different variable/function naming conventions are mixed 
       in the rest of this function. Sorry. */
    g = dimacsParse(stdin);

    printf("c nodes:   %14d    arcs: %15d\n", g->n, g->m);
#ifdef LONGOUTPUT
    printf("graph read in %.2f seconds\n\n",timer()- tmStart);
#endif

    tmStart = timer ();

    /* initialization */
    cLeaves = cTrees = 0;
    
    /* clean up multiple edges */
    compact(g);
    
    makeHeap ( h, g->currentN );
    PRpreprocess(g, 0.95, 0.5, 0.5);
    tmFound = g->dtime;
    
    wtOpt = g->minCap;

    if (g->currentN > 2)
      {
	/* move into different data-structures */
	pmg = pmgConvertToMetaGraph(g);

#ifdef DEPTHTIMES
	iDepthMax = (int)(log((float)pmg->cmv)/log(2.0));
#endif
	
#ifdef SAVECUT
	rgpmvOpt = Malloc(cmvIn,metavertex*);  /*for FSaveMinCut*/
#endif
	
#if (PR2==TIGHT)
	iTime=1;
	for (i=0; i<pmg->cmv; i++)
	  {
	    pmg->rgpmv[i]->iTime = 0;
	  }
#endif
	
	WalkMinCuts(pmg,dFailProb,NewTrivialCut);
      }

/* output in format like noi */

  printf("c ttime: %14.2f    capacity: ", timer()-tmStart); 
  fprintf_wt(stdout, "%15", wtOpt);
  printf("\n");

  printf("c PR12: %15d    PR34:   %15d\n",cPR12, cPR34);
#if (DEBUG >=1)
  cScans += g->edgeScanCnt[PRPRE12] + g->edgeScanCnt[PRPRE34];
  printf("c escans: %15d\n",cScans);
#endif
  printf("c dtime:%15.2f    leaves:   %12d    \n",tmFound-tmStart, cLeavesTot);

#if (ONECUT && defined(VERBOSE))
  printf("\nVertices on one side of the cut:\n");
  printf(" %d",rgpmvOpt[0]->name);
  for(i=1;i<cmvOpt;i++)
    printf(", %d",rgpmvOpt[i]->name);
  printf("\n");
#endif
    
#ifdef DEPTHTIMES
  for (i=0; i<=cDepthMax; i++) {
    if (depthCount[i]) {
      printf("c depth %2d",i);
      printf("%9d avg size: %6.1f tm: %4.2f/%4.2f PR12: %7d/%5.2f\n",
	     depthCount[i], depthNodes[i]/(float)depthCount[i], depthTimes[i], 
	     depthTimes[i]/depthCount[i],pr12Count[i], 
	     pr12Count[i]/(float)depthCount[i]);
    }
  }
#endif

  printf("\n");
    
  return(0);
}

    
