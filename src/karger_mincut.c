/* Karger's sampling+tree packing minimum cut algorithm. 
 * See karger_mincut.h
 * This file is in good shape to use use a fractional packing algorithm of
 * Plotkin, Shmoys, and Tardos in place of Gabow's algorithm, but I'm not yet 
 * confident that it is right. I apologize for the (temporarily) useless goo
 * this introduces.
 *
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 * $Id: karger_mincut.c,v 1.20 1997/05/16 15:30:32 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#ifndef __CHECKER__
#include <assert.h>
#endif
#include "karger_mincut.h"
#ifdef GABOW
#include "gabow_pack.h"
#endif
#ifdef PST
#include "pst.h"
#endif
#include "pr.h"
#include "contract.h"
#include "random.h"
#include "timer.h"
#include "heap.h"
#include "dyn_tree/dyn_tree.h"
#include "ni_sparse_cert.h"
#include "memassert.h"

/***********************************************************************
 * Prototypes
 ***********************************************************************/

/* see below for descriptions of the functions */
static node *find(node *v);
static void saveKargerCut(graph *g, node *v1, node *v2, char descendant, 
			  weight_t value);
void sampleGraph(graph *g, double sampleprob);
static int binomial(weight_t n, double p);
static int treecmp(tree **t1, tree **t2);
#ifdef CHECK_GABOW_CUT
static weight_t dfs0(graph *g, node *v);
#endif
#if defined(GABOW) && (defined(INTWEIGHTS) || defined(CHECK_GABOW_CUT))
static void saveGabowCut(graph *g, packing pack, weight_t value);
#endif
static void dfs1_resp(graph *g, node *v);
static void dfs1_resp_desc(graph *g, node *v);
static void dfs2_respA(graph *g, node *v);
static void dfs2_respB(graph *g, node *v);
static void make_dyn_tree(graph *g, node *v);
static node *check_boughs(graph *g, node *v);
#ifdef SAVECUT
static void dfs_mark(graph *g, node *v, int label);
#endif

/***********************************************************************
 * Global Variables
 ***********************************************************************/

heap h;                        /* heap for PR 34 tests and/or PST */
static graph *gg;      /* treecmp needs g, but can only have 2 args 
			  'cause passed to qsort */

/***********************************************************************
 * Function definitions
 ***********************************************************************/

weight_t karger_mincut(graph *g, double successprob, double disrespect,  
		       double gabowbailout, double maxdyndensity)
{
  double logfail, epsilon, probconst, sampleprob;    /* sampling probability */
  float tm, tim;
  weight_t value=0, savedMinCap;
  int i, k=0, lim, l=0;
  int vname;
  int treetotal = 0;
  packing pack=NULL;
  node *v;
  tree **trees;
  graph *backup, *backup2;
  double matula_epsilon = 0.1;

  assert(successprob < 1);
  assert(disrespect*3. > 2.);
  logfail = log(1/(1-successprob));

  tm = timer();

  /* clean up input */
  compact(g);
  if (g->minCap == 0 || g->currentN <= 2) /* if disconnected or small, quit */
    {
      printf("c trees total: %10d    final: %14d\n", treetotal, k);      
      return g->minCap; 
    }

  makeHeap(h, g->currentN);  /* initialize heap */

  /* Preprocess with PR heauristics */
#ifndef NO_PR
  PRpreprocess(g, 0.95, 0.5, 0.5);
  if (g->currentN <= 2 )
    {
      freeHeap(h);
      printf("c trees total: %10d    final: %14d\n", treetotal, k);      
      return g->minCap; 
    }
#endif  

  /* get a 2+matula_epsilon approx to mincut and sparsify */
  backup = makeBackupGraph(g);
  matulaApprox(g, matula_epsilon);
  restoreBackupGraph(g, backup);
  freeBackupGraph(backup);
  v = sparseCertContract(g, g->minCap);
  if (g->currentN > 1)
    {
      compact(g);
#ifndef NO_PR
      sourcePR(g, v, 0.95);
      PRpass(g, 0.95, 0.5, 0.5, 0);
#endif
    }
  if (g->currentN <= 2 )
    {
      freeHeap(h);
      printf("c solved by sparsification\n");
      printf("c trees total: %10d    final: %14d\n", treetotal, k);      
      return g->minCap; 
    }

  printf("c post-matula nodes: %d   total preprocesstime: %f \n", 
	 g->currentN, timer()-tm);
  
  /* compute sampling probability */
  tm = timer();
  epsilon = (3*disrespect-2)/(3.*disrespect+sqrt(logfail/(2.*log(g->currentN)+logfail)));
  probconst = (4*log(g->currentN)+2*logfail)/(epsilon*epsilon); 
  /*  probconst = 4*log(g->currentN)+2*logfail; */
  sampleprob = probconst*(2+matula_epsilon)/(double)g->minCap; 
#ifdef VERBOSE
  printf("c epsilon: %f\n", epsilon);
  printf("c expected #trees if estimate is off: %.0f (= %.2f log n)\n", 
	 probconst,  probconst/log(g->currentN));
  printf("c expected #trees if estimate right: %.0f (= %.2f log n)\n",
	 probconst*(2+matula_epsilon),  
	 probconst*(2+matula_epsilon)/log(g->currentN)); 
  printf("c cut estimate: "); fprintf_wt(stdout, "%", g->minCap); printf("\n");
#endif
  
#ifdef CHEAT
  sampleprob = 6*log(g->currentN)*(2+matula_epsilon)/(double)g->minCap;
#endif

  /* if the mincut is a small integer (ie, sampleprob is large and the
     weights are integers), then we should just bailout into Gabow's
     mincut algorithm */
#if defined(INTWEIGHTS) && defined(GABOW)
  /* if integer weights we can bail out when when sampling prob >= const */
  if (sampleprob >= gabowbailout)
    {
#ifdef VERBOSE
      printf("c bailed to straight gabow\n");
#endif /* VERBOSE */
      g->root = NODE(g,0);
      pack = gabow_pack(g, (int)g->minCap, pack);
      if (pack->k < g->minCap)
	saveGabowCut(g, pack, (weight_t)pack->k);
      treetotal += pack->k;
      k = 0;
      goto cleanup;
    }
#endif /* INTWEIGHTS && GABOW */

  /****** sample the graph *****/
  
  backup = makeBackupGraph(g);
  savedMinCap = g->minCap;
  sampleGraph(g, sampleprob);
#ifdef VERBOSE
  printf("c sampled graph has %d edges, and min trivial cut %d\n",
	 g->currentM/2, (int)g->minCap);
#endif

  /****** approx mincut of sample (if we have reason to care) ******/
  if (disrespect < 1)
    {
      backup2 = makeBackupGraph(g);
      matulaApprox(g, matula_epsilon);
      restoreBackupGraph(g, backup2);
      freeBackupGraph(backup2);
#ifdef VERBOSE
      printf("c MatulaApprox says %d <= mincut of sample <= %d\n", 
	     (int)g->minCap, (int)(g->minCap/(2.+matula_epsilon)));
#endif  
    }

  /******  find a tree packing ******/
  g->root = NODE(g,0);
  g->context = PACKING;
#ifdef GABOW
  tim = timer();
  lim = g->minCap*disrespect; 
  printf("c safe to quit after %d trees\n", lim);
  pack = gabow_pack(g, lim, pack);
  k = pack->k;
#ifdef VERBOSE
  printf("c gabow found %d trees in %f\n", k, timer()-tim);
#endif
  treetotal += k;
#endif /* GABOW */

#ifdef PST
  pack = pst_pack(root, g->currentN, sample_m, (int)(probconst), .16, pack);
  cont = (pack->lambda > 4);
#endif /* PST */

  restoreBackupGraph(g, backup);
  freeBackupGraph(backup);
  g->minCap = savedMinCap;

#ifdef CHECK_GABOW_CUT
  if (pack->onesideroot != NULL)
    {
      /* check vertex partition returned by gabow_pack in real graph */
      /* this gets us a (1+epsilon) approx */
      ForAllNodes(g, v) v->active = 0;
      g->currentTree = pack->last_tree;
      if ((value=dfs0(g, pack->onesideroot)) < g->minCap)
	saveGabowCut(g, pack, value);
      if (g->minCap == 0) { k = 0; goto cleanup;}
    }
#endif

  printf("c pack time: %f\n", timer()-tm);
  
  /*** find min cut 1 or 2-respecting each tree ***/
  g->context = MAIN;

  /* pick the trees to check for 2-respect */
  tm = timer();
  k = pack->k;

  if (lim > k && disrespect < 1) 
    lim = -logfail/log(disrespect);

#ifdef CHEAT
  lim = 10;
#endif

  tim = timer();
  /* prune duplicate trees */
  memassert(trees = (tree **)malloc(pack->k*sizeof(tree *)));
  i = 0;
  ForAllTrees(pack, g->currentTree)
     trees[i++] = g->currentTree;
  assert(i == k);
  gg = g;
  qsort(trees, k, sizeof(tree **), 
	(int (*)(const void *, const void *))treecmp);

  l=1;
  for(i=1; i<k; i++)
    if (treecmp(trees+i-1, trees+i))
      trees[l++] = trees[i]; 

#ifdef VERBOSE
  printf("c pruned %d redundant trees in %f time\n", k-l, timer()-tim);
#endif
  pack->k = k = l;
      
  if (lim > k) lim = k;

  printf("c have to check %d trees for 2-resp cuts\n", lim);

  if (lim < k)
    RandomPermuteObj(k, (char *)trees, sizeof(tree *));  

  if (g->currentM/2 > g->currentN*(g->currentN-1)/2*maxdyndensity)
    {
      /* use simple O(n^2) algorithm to check for 2-respecting cuts */
      memassert(g->subtreecuts=(weight_t*)
		malloc(g->currentN*g->currentN*sizeof(weight_t)));
      memassert(g->is_descendant=(char *)
		malloc(g->currentN*g->currentN*sizeof(char)));
      
      for(i=0, g->currentTree=trees[0]; i<k; i++, g->currentTree=trees[i])
	{
	  /* initialize data structures used in computing cut */
	  ForAllNodes(g, v)
	    {
	      v->rho = 0;
	      v->rhodown = 0;
	      v->capdown = v->cap;
	      v->cutdown = 0;
	      v->parent = NULL;
	    }
	  
	  if (i<lim)   /* clear old values */
	    {
	      memset(g->subtreecuts, 0, 
		     g->currentN*g->currentN*sizeof(weight_t));
	      /* n instead of n^2 ok for descendant */
	      memset(g->is_descendant, 0, g->currentN*sizeof(char)); 
	    }
	  
	  if (i >= lim)
	    /* just compute 1-respecting cuts */
	    dfs1_resp(g, g->root);
	  else
	    {
	      /* compute cuts 1-respecting the tree and tree descendant 
		 relation */
	      dfs1_resp_desc(g, g->root);
	      /* compute cuts between all pairs of subtrees */
	      dfs2_respA(g, g->root);
	      dfs2_respB(g, g->root);
	    }
	  
	  /* put together information from dfs[123] to find the minumum */
	  if (i >= lim)
	    {
	      ForAllNodes(g, v)
		if (v->cutdown < g->minCap && v!= g->root) 
		  saveKargerCut(g, v, NULL, 0, v->cutdown);  
	    }
	  else
	    ForAllNodes(g, v)
	    {
	      if (v == g->root) continue;
	      
	      if (v->cutdown < g->minCap) 
		saveKargerCut(g, v, NULL, 0, v->cutdown); 
	      
	      for(l=(vname=INDEX(v))+1; l<g->currentN; l++)
		{
		  if (g->is_descendant[vname*g->currentN+l])
		    {
		      if ((value = NODE(g,l)->cutdown - v->cutdown +
			   2*g->subtreecuts[vname*g->currentN+l]) < g->minCap)
			saveKargerCut(g, NODE(g,l), v, 1, value);
		    }
		  else if (g->is_descendant[l*g->currentN+vname])
		    {
		      if ((value = v->cutdown - NODE(g,l)->cutdown +
			   2*g->subtreecuts[l*g->currentN+vname]) < g->minCap)
			saveKargerCut(g, v, NODE(g,l), 1, value);
		    }
		  else
		    {
		      if ((value = v->cutdown + NODE(g,l)->cutdown -
			   2*g->subtreecuts[vname*g->currentN+l]) <g->minCap)
			saveKargerCut(g, v, NODE(g,l), 0, value);
		    }
		}
	    }
	}
      free(g->subtreecuts);
      free(g->is_descendant);
    }
  else
    {
      /* use fancy O(m log^2 n) algorithm (using dynamic trees) to check for 
	 2-respecting cuts */
    
      dyn_init(g->currentN);    /* alloc memory for dynamic trees */
      backup = makeBackupGraph(g);
      
      for(i=0, g->currentTree=trees[0]; i<k; i++, g->currentTree=trees[i])
	{
	  /* initialize data structures used in computing cut */
	  ForAllNodes(g, v)
	    {
	      v->rho = 0;
	      v->rhodown = 0;
	      v->capdown = v->cap;
	      v->cutdown = 0;
	      v->parent = NULL;
#ifdef SAVECUT
	      savefirst(g->currentTree, v) = NULL;
#endif
	    }
	  
	  /* compute cuts 1-respecting the tree and tree descendant relation */
	  dfs1_resp(g, g->root);
	  
	  /* see if any 1-respecting cuts are a new minimum */
	  ForAllNodes(g, v)
	    {
	      if (v->cutdown < g->minCap && v!= g->root) 
		saveKargerCut(g, v, NULL, 0, v->cutdown); 	  
	    }
	  
	  if (i < lim)
	    {
	      /* build dynamic tree */
	      dyn_reset();
	      make_dyn_tree(g, g->root);
	      
	      /* check all boughs and contract them away */
	      while (HasChildren(g->currentTree, g->root))
		/* be careful about root = top of bough */
		if (HasOneChild(g->currentTree, g->root)) 
		  {
		    if (check_boughs(g, first(g->currentTree, g->root)->head)
			!= NULL)
		      {
#ifdef SAVECUT
			first(g->currentTree, g->root)->next =
			  savefirst(g->currentTree, g->root);
			savefirst(g->currentTree, g->root) =
			  first(g->currentTree, g->root);
			first(g->currentTree, g->root) = NULL;
#endif
			break;
		      }
		  }
		else
		  check_boughs(g, g->root);

	      /* undo damage done by contractions */
	      restoreBackupGraph(g, backup);
	    }
	}
      
      freeBackupGraph(backup);
      dyn_free();
    }

#ifdef SAVECUT
  /* set in_cut fields for a non-trivial cut */
  g->currentTree = g->savet;
  if (g->save1 != NULL)
    {
      ForAllNodes(g, v)
	{
	  v->in_cut = 0;
	  if (first(g->currentTree, v) == NULL)
	    first(g->currentTree, v) = savefirst(g->currentTree, v);
	}
      if (g->save2 != NULL)
	if (g->saveDescendant)
	  {
	    dfs_mark(g, g->save1, 1);
	    dfs_mark(g, g->save2, 0);
	  }
	else
	  {
	    dfs_mark(g, g->save1, 1);
	    dfs_mark(g, g->save2, 1);
	  }
      else
	dfs_mark(g, g->save1, 1);
#ifdef CONTRACT
      ForAllVertices(g, v)
	v->in_cut = findLeader(v)->in_cut;
#endif /* CONTRACT */
    }
#endif /* SAVECUT */

  printf("c respect time: %f\n", timer()-tm);

#if defined(GABOW) && (defined(CHECK_GABOW_CUT) || defined(INTWEIGHTS))
 cleanup:
#endif
  printf("c trees total: %10d    final: %14d\n", treetotal, k);
  free_packing(pack);
  freeHeap(h);
  return g->minCap;
}

/******************************************************************/

/* save a new non-trivial cut
   Note that this only takes constant time! */
static void saveKargerCut(graph *g, node *v1, node *v2, char descendant, 
			  weight_t value)
{
  g->minCap = value;
  g->dtime = timer();
  g->cutCount++;

#ifdef VERBOSE
  if (v2 == NULL)
    printf("c new min in 1-respect (tree %d): ", g->currentTree->name);
  else
    printf("c new min in 2-respect (tree %d): ", g->currentTree->name);
  fprintf_wt(stdout, "%", g->minCap);
  printf("\n");
#endif

#ifdef SAVECUT
  g->save1 = v1;
  g->save2 = v2;
  g->savet = g->currentTree;
  g->saveDescendant = descendant;
#endif
}

#if defined(GABOW) && (defined(INTWEIGHTS) || defined(CHECK_GABOW_CUT))
/* save a new non-trivial cut determined by a tree-packing */
static void saveGabowCut(graph *g, packing pack, weight_t value)
{
#ifdef SAVECUT
  node *v;
#endif

  g->minCap = value;
  g->dtime = timer();
  g->cutCount++;

#ifdef VERBOSE
  printf("c new min in gabow: ");
  fprintf_wt(stdout, "%", g->minCap);
  printf("\n");
#endif

#ifdef SAVECUT
  ForAllNodes(g, v)
    v->in_cut = 0;
  dfs_mark(g, pack->onesideroot, 1);
#ifdef CONTRACT
  ForAllVertices(g, v)
    v->in_cut = findLeader(v)->in_cut;
#endif
#endif
}
#endif

/* find, of union-find used for finding least common ancestors */
static node *find(node *v)
{
  node *p=v, *q;
  while (p->parent != NULL) p = p->parent;
  if (p != v)    /* compress path */
    while((q=v)->parent != p) 
      {
	v = v->parent;
	q->parent = p;
      }
  return p;
}

/* sample the graph. DESTRUCTIVE */
void sampleGraph(graph *g, double sampleprob)
{
  arc *a, *ra, **i;
  node *v, *w;
  int newCap;

  ForAllNodes(g, v) 
    v->cap = 0;

  /*  ForAllArcs(g, a) */
  for(i=g->currentArcs; (a=*i) != NULL; i++)
    {
      ra = a->rev;
      if (a->index < ra->index)  /* only sample each edge once */
	{
	  g->edgeScanCnt[g->context]++; 
	  v = ra->head;
	  w = a->head;
	  newCap = binomial(a->cap, sampleprob);
	  if (newCap == 0)
	    {
	      DeleteArc(g, v, a);
	      DeleteArc(g, w, ra);
	      i--;
	    }
	  else
	    {
	      a->cap = ra->cap = newCap;
	      v->cap += newCap;
	      w->cap += newCap;
	    }
	}
      
    }
  
  g->minCap = MAXWEIGHT;
  ForAllNodes(g, v)
    if (v->cap < g->minCap) g->minCap = v->cap;
}

/* generate a random number according to the binomial distribution 
 * NOTE: this function approximates the binomial distribution by the 
 * poisson distribution, as implemented by Knuth. 
 *
 * UGLY NOTE: The efficient way to do this is to multiply q by a
 * random number each time. Sadly, it is easy for exp(-p*n) to
 * underflow, causing this code to loop forever. So unless you're sure
 * that p*n is small, leave it paying a log each iteration. 
 */
static int binomial(weight_t n, double p)
{
  /*  double q=1., emu = exp(-p*n);  */
  double q=0., emu = -p*n; 
  int ans=-1;
  
  do 
    {
      ans++;
      /*      q *= DblUnitRandom(); */
      q += log(DblUnitRandom()); 
    } while (q >= emu); 
  return ans;
}

static int treecmp(tree **t1, tree **t2)
{
  int i;

  for(i=1; i<gg->currentN && parent_edge((*t1), NODE(gg,i))->head == 
	parent_edge((*t2), NODE(gg,i))->head;
      i++);
  if (i == gg->currentN) return 0;
  return INDEX(parent_edge((*t1), NODE(gg,i))->head) - 
    INDEX(parent_edge((*t2), NODE(gg,i))->head);
}

#ifdef CHECK_GABOW_CUT
/* determine value of cut in full graph induced by vertex partition of
   min cut in sample */
static weight_t dfs0(graph *g, node *v)
{
  node *w;
  tedge *e;
  weight_t ans=0, wt;
  arc *a;

  ForAllIncidentArcs(v, a)
    {                 /* note this part is NOT recursive */
      g->edgeScanCnt[g->context]++;
      wt = a->cap;
      if (a->head->active) ans -= wt;
      else ans += wt;
    }

  v->active = 1;
  ForAllTreeEdges(g->currentTree, v, e)
    {
      g->edgeScanCnt[g->context]++;
      if (!((w=HEAD(e))->active))
	ans += dfs0(g, w);
    }

  return ans;
}
#endif /* CHECK_GABOW_CUT */

/* treefix to compute 1-respecting cuts */
static void dfs1_resp(graph *g, node *v)
{
  node *w, *lca;
  tedge *e;
  weight_t wt;
  arc *a;

  ForAllTreeEdges(g->currentTree, v, e)
    {
      g->edgeScanCnt[g->context]++;
      if (IsParentEdge(g->currentTree, v, e)) 
	{
	  DeleteTreeEdge(g->currentTree, v, e);
	  continue; /* don't go up */
	}      
      w=HEAD(e);
	  
      dfs1_resp(g, w);
      
      /* most things are computed up from the bottom */ 
      w->parent = v;     /* union, of union-find used below */
      v->rhodown += w->rhodown;
      v->capdown += w->capdown;
    }

  ForAllIncidentArcs(v, a)
    {                 /* note this part is NOT recursive */
      g->edgeScanCnt[g->context]++;
      wt = a->cap;
      w = a->head;

      /* this is offline LCA done by union-find (path compression)
	 (R.E. Tarjan, JACM 26:4, p. 698)
	 basically, whenever we're done with a vertex we union it with its 
	 parent (done above), and the the second time we see an edge we do 
	 a find on the other end to get the LCA */
      if (w->parent != NULL)
	{
	  lca = find(w);
	  lca->rho += wt;
	  a->lca = a->rev->lca = lca;
	}
    }

  v->rhodown += v->rho;
  v->cutdown = v->capdown - 2*v->rhodown;
}

/* treefix to compute cuts 1-respecting tree and is_descendant */
static void dfs1_resp_desc(graph *g, node *v)
{
  node *w;
  int vname, wname;
  int i;
  tedge *e;
  weight_t wt;
  arc *a;

  vname = INDEX(v);
  g->is_descendant[vname*g->currentN+vname] = 1;
  ForAllTreeEdges(g->currentTree, v, e)
    {
      g->edgeScanCnt[g->context]++;
      if (IsParentEdge(g->currentTree, v, e)) 
	{
	  DeleteTreeEdge(g->currentTree, v, e);
	  continue; /* don't go up */
	}

      w = HEAD(e); wname = INDEX(w);
      
      /* descendant is computed down from the top */
      for (i=0; i<g->currentN; i++)
	g->is_descendant[wname*g->currentN + i] = 
	  g->is_descendant[vname*g->currentN + i];
      g->is_descendant[wname*g->currentN+vname] = 1;
      
      dfs1_resp_desc(g, w);
      
      /* most things are computed up from the bottom */ 
      w->parent = v;     /* union, of union-find used below */
      v->rhodown += w->rhodown;
      v->capdown += w->capdown;
    }

  ForAllIncidentArcs(v, a)
    {                 /* note this part is NOT recursive */
      g->edgeScanCnt[g->context]++;
      wt = a->cap;
      w = a->head;
      
      /* init subtreecuts with edge weights */
      g->subtreecuts[vname*g->currentN + INDEX(w)] += wt; 

      /* this is offline LCA done by union-find (path compression)
	 (R.E. Tarjan, JACM 26:4, p. 698)
	 basically, whenever we're done with a vertex we union it with its 
	 parent (done above), and the the second time we see an edge we do 
	 a find on the other end to get the LCA */
      if (w->parent != NULL)
	find(w)->rho += wt;
    }

  v->rhodown += v->rho;
  v->cutdown = v->capdown - 2*v->rhodown;
}

/* treefix to compute C(v, w^down) */
static void dfs2_respA(graph *g, node *v)
{
  node *w;
  int i, vname, wname;
  tedge *e;

  vname = INDEX(v);
  ForAllTreeEdges(g->currentTree, v, e)
    {
      g->edgeScanCnt[g->context]++;
      w=HEAD(e);
      wname = INDEX(w);
      dfs2_respA(g, w);
      for (i=0; i<g->currentN; i++)
	if (!g->is_descendant[i*g->currentN+vname])
	  g->subtreecuts[vname*g->currentN + i] += 
	    g->subtreecuts[wname*g->currentN + i];
    }
  /* wipe out any initial values that got filled in by dfs1_resp_desc but 
     don't make sense given the descendant relation */
  for (i=0; i<g->currentN; i++)   
    if (g->is_descendant[i*g->currentN+vname])
      g->subtreecuts[vname*g->currentN + i] = 0;
}

/* treefix to compute C(v^down, w^down) given C(v, w^down) (from dfs2_respA) */
static void dfs2_respB(graph *g, node *v)
{
  node *w;
  int i, vname, wname;
  tedge *e;

  vname = INDEX(v);
  ForAllTreeEdges(g->currentTree, v, e)
    {
      g->edgeScanCnt[g->context]++;
      w=HEAD(e);
      wname = INDEX(w);
      dfs2_respB(g, w);
      for (i=0; i<g->currentN; i++)
	if (!g->is_descendant[vname*g->currentN+i])
	  g->subtreecuts[i*g->currentN + vname] += 
	    g->subtreecuts[i*g->currentN + wname];
    }
}


/* build up a dynamic tree version of g->currentTree */
static void make_dyn_tree(graph *g, node *v)
{
  node *w;
  tedge *e;
  
  dyn_make_tree(v, 0);
  ForAllTreeEdges(g->currentTree, v, e)
    {
      g->edgeScanCnt[g->context]++;
      w=HEAD(e);
      make_dyn_tree(g, w);
      dyn_link(w, v, w->cutdown);
    }
}

/* check all of g->currentTree's boughs for better 2-respecting cuts */
static node *check_boughs(graph *g, node *v)
{
  weight_t val, tmpval;
  tedge *te, *e; 
  node *tmppartner, *partner=NULL;
  arc *a;

  e = first(g->currentTree, v);
  first(g->currentTree, v) = NULL;
  
  if (e != NULL) /* v has children */
    {
      if (e->next != NULL) /* v has many children => not on a bough */
	{
	  while(e != NULL)
	    {
	      g->edgeScanCnt[g->context]++;
	      if (check_boughs(g, e->head) != NULL) 
		{  /* e->head is the top of a bough */

		  /* clean up e->head */
		  ForAllIncidentArcs(e->head, a)
		    {
		      g->edgeScanCnt[g->context]++;
		      if (a->lca != a->head && a->lca != e->head)
			{
			  dyn_add_value(a->head, 2*a->cap);
			  dyn_add_value(a->lca, -4*a->cap);
			}
		      else
			dyn_add_value(a->head, -2*a->cap);
		    }
		  /* contract it away */
		  compactContractLite(g, v, e->head);

#ifdef SAVECUT
		  te = e;
#endif
		  e = e->next;
#ifdef SAVECUT
		  te->next = savefirst(g->currentTree, v);
		  savefirst(g->currentTree, v) = te;
#endif

		}
	      else
		{   /* child not top of bough, put back on edge list */
		  te = e;
		  e = e->next;
		  te->next = first(g->currentTree, v);
		  first(g->currentTree, v) = te;
		}
	    }
	}
      else  /* v has one child => might be on bough */
	{
	  if ((partner=check_boughs(g, e->head)) != NULL)
	    { /* on a bough */
	      /* accumulate weights */
	      ForAllIncidentArcs(v, a)
		{
		  g->edgeScanCnt[g->context]++;
		  if (a->lca != a->head && a->lca != v)
		    {
		      dyn_add_value(a->head, -2*a->cap);
		      dyn_add_value(a->lca, 4*a->cap);
		    }
		  else
		    dyn_add_value(a->head, 2*a->cap);
		}
	      /* look for best incomparable partner */	      
	      val = partner->rho;
	      ForAllIncidentArcs(v, a)
		{
		  g->edgeScanCnt[g->context]++;
		  if (a->lca != a->head && a->lca != v)
		    {
		      tmpval=dyn_find_value(tmppartner=dyn_find_min(a->head));
		      if (tmpval < val)
			{
			  val = tmpval;
			  partner = tmppartner;
			}
		    } 
		}
	      
	      /* save value to return w/ partner */
	      partner->rho = val;
	      
	      /* found a new mincut ? */
	      if ((tmpval = val+v->cutdown) < g->minCap)
		saveKargerCut(g, v, partner, 0, tmpval);
	      
	      tmppartner=dyn_find_min(parent_edge(g->currentTree, v)->head);
	      if (tmppartner != v && tmppartner != g->root &&
		  ((val=dyn_find_value(tmppartner)
		    -v->cutdown-4*v->rhodown) < g->minCap))
		saveKargerCut(g, tmppartner, v, 1, val);
	      
	      compactContractLite(g, v, e->head);

#ifdef SAVECUT
	      e->next = savefirst(g->currentTree, v);
	      savefirst(g->currentTree, v) = e;
#endif
	    }
	  else  /* not on a bough, put child back on edge list */
	    first(g->currentTree, v) = e;
	}
    }
  else  /* v is a leaf */ 
    {
      /* accumulate weights */
      ForAllIncidentArcs(v, a)
	{
	  g->edgeScanCnt[g->context]++;
	  if (a->lca != a->head && a->lca != v)
	    {
	      dyn_add_value(a->head, -2*a->cap);
	      dyn_add_value(a->lca, 4*a->cap);
	    }
	  else
	    dyn_add_value(a->head, 2*a->cap);
	}
      /* look for best comparable partner */
      partner=dyn_find_min(parent_edge(g->currentTree, v)->head);

      /* found a new mincut ? */
      if (partner != v && partner != g->root && 
	  ((val=dyn_find_value(partner)-v->cutdown-4*v->rhodown) 
	   < g->minCap))
	saveKargerCut(g, partner, v, 1, val);

      /* look for best incomparable partner */
      val = MAXWEIGHT;
      partner = NULL;
      ForAllIncidentArcs(v, a)
	{
	  g->edgeScanCnt[g->context]++;
	  if (a->lca != a->head && a->lca != v)
	    {
	      tmpval = dyn_find_value(tmppartner=dyn_find_min(a->head));
	      if (tmpval < val)
		{
		  val = tmpval;
		  partner = tmppartner;
		}
	    }
	}

      if (partner != NULL)
	{
	  /* save value to return w/ partner */
	  partner->rho = val;
	  
	  /* found a new mincut ? */
	  if ((tmpval = val+v->cutdown) < g->minCap)
	    saveKargerCut(g, v, partner, 0, tmpval);
	}
      else
	{  /* we have to return non-NULL for a leaf, even if there isn't
	      an incomparable partner, so we (ab)use root */
	  partner = g->root;
	  partner->rho = g->minCap;
	}
    }
  
  return partner;
}

#ifdef SAVECUT
/* dfs to mark nodes in cut */
static void dfs_mark(graph *g, node *v, int label)
{
  tedge *e;

  v->in_cut = label;
  ForAllTreeEdges(g->currentTree, v, e)
    {
      g->edgeScanCnt[g->context]++;
      if (IsParentEdge(g->currentTree, v, e)) continue; /* don't go up */
      dfs_mark(g, HEAD(e), label);
    }
}
#endif /* SAVECUT */


