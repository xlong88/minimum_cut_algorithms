/* Gabow's round robin tree packing algorithm. See gabow_pack.h
 *
 * Author: Matthew Levine (mslevine@theory.lcs.mit.edu)
 * $Id: gabow_pack.c,v 1.9 1997/04/15 15:18:17 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include "gabow_pack.h"
#include "graph.h"
#include "allocate.h"
#include "timer.h"

/***********************************************************************
 * Prototypes
 ***********************************************************************/

static void dfs_init(node *v, tree *t, tedge *parent_edge, int depth, 
		     node *label);
#ifdef __GNUC__
static inline int label(tree *t, tedge *e, tedge *label);
#else
static int label(tree *t, tedge *e, tedge *label);
#endif
static int label_unused_edges(node *v, tedge *label);
static int search(node *vertex);

static void dfs_init_pack(graph *g, node *root);
static void find_dfs_tree(graph *g, node *v);

/***********************************************************************
 * Variables global to this file
 ***********************************************************************/

static int dfs_visited;          /* for dfs_init_packing */

/* lists (in array form) */
static tedge *path_end_edges;    /* last edges of augmenting paths */
static tedge *q;                 /* queue of edges for labelling */
static tedge *q_last;             /* last item in q */

static long *edgeScanCnt;

/* the answer (it is more convenient not to pass them around) */
static packing ans;

/***********************************************************************
 * Macros
 ***********************************************************************/

#define depth(t, v) (t->info[INDEX(v)].depth)
#define is_labeled(t, v) (t->info[INDEX(v)].is_labeled)

#define AddTreeEdge(t, v, e) \
{ \
  e->next = t->info[INDEX(v)].first; \
  if (t->info[INDEX(v)].first != NULL) t->info[INDEX(v)].first->prev = e; \
  t->info[INDEX(v)].first = e; \
  e->prev = NULL; \
}

/***********************************************************************
 * Function definitions
 ***********************************************************************/

packing gabow_pack(graph *g, int c, packing ans2)
{ 
  node *ftree_roots, *new_ftree_roots;    /* current roots of active ftrees */
  node *v, *w;
  tedge *e, *e2, *re;
  tree *t;
  arc *ta;
  /*  float tm;*/
  
  /* initialize */
  ans = ans2;
  if (ans == NULL)
    {
      memassert(ans = (packing)malloc(sizeof(struct packing)));
      ans->spare_trees = NULL;
      ans->tedge_allocator = make_allocator(g->currentN, sizeof(tedge));
    }
  else 
    {
      if (ans->last_tree != NULL)
	{
	  ans->last_tree->next_tree = ans->spare_trees;
	  ans->spare_trees = ans->first_tree;
	}
    }
  ans->k = 0;
  ans->first_tree = NULL;
  ans->last_tree = NULL;
  ans->onesideroot = NULL;
  
  edgeScanCnt = g->edgeScanCnt+g->context;
  
  q = (tedge *)get_node(ans->tedge_allocator);  /* dummy node (so q_last has 
				  something to point to when q is empty) */
  q->head = NULL;
  q->auxnext = NULL;
  q_last = q;

  ForAllNodes(g, v)
    v->used_first = NULL; 

  dfs_init_pack(g, g->root);
  if (ans->k == 0) return ans; 

  /*  tm = timer(); */

  while(ans->k <c) 
    {
      ans->k++;
      /* printf("trying tree %d. %f since last tree\n", ans->k, timer()-tm);
      tm = timer(); */

      /* get new tree */
      if (ans->spare_trees != NULL)
	{
	  t = ans->spare_trees;
	  ans->spare_trees = t->next_tree;
	}
      else
	{
	  memassert(t= (tree *)malloc(sizeof(tree)));
	  memassert(t->info = (struct treeinfo *)
		    malloc(g->currentN*sizeof(struct treeinfo)));
	}
      if (ans->k == 1)  /* first is special case */
	{
	  ans->last_tree = t;
	  ans->first_tree = ans->last_tree;
	}
      else
	{
	  ans->last_tree->next_tree = t;
	  ans->last_tree = ans->last_tree->next_tree;
	}
      t->name = ans->k-1;	  
      t->next_tree = ans->first_tree;
      t->was_changed = 1;
      
      /* init kth forest. (each vertex a tree) */
      new_ftree_roots = NULL;
      ForAllNodes(g,v)
	{
	  first(t,v) = NULL;
	  if (v==g->root) continue;
	  v->next_ftree = new_ftree_roots;
	  if (new_ftree_roots != NULL) new_ftree_roots->prev_ftree = v;
	  new_ftree_roots = v;
	}
      new_ftree_roots->prev_ftree = NULL;

      /* while there are active ftrees... */
      while(new_ftree_roots != NULL)
	{
	  ftree_roots = new_ftree_roots;
	  new_ftree_roots = NULL;

	  ForAllTrees(ans, t) /* (re)set depths, parents, vertex labels */
	    if (t->was_changed)	       
	      {
		t->was_changed = 0;
		dfs_init(g->root, t, NULL, 0, NULL);
	      }

	  /* this is expensive; should try to be lazier... */
	  ForAllNodes(g,v)
	    {
	      ForAllTrees(ans, t)     /* clear marks */
		is_labeled(t, v) = 0;
	      is_labeled(t, v) = 0;
	      v->has_no_labeled_edges = 1;
	    }
	  
	  /* init the ftrees too... */
	  t = ans->last_tree;
	  dfs_init(g->root, t, NULL, 0, g->root);
	  for(v=ftree_roots; v!=NULL; v=v->next_ftree)  
	    {
	      dfs_init(v, t, NULL, 0, v);
	      v->active = 1;
	    }
	  
	  path_end_edges = NULL;     /* no augmenting paths yet */ 

	  /* for each active ftree... */
	  for(v=ftree_roots; v!=NULL; v=v->next_ftree)  
	    /* search for an augmenting path */
	    if (!search(v)) 
	      {
		ans->onesideroot = v; /* no path means ftree defines cut */
		ans->k--;  /* this k failed, so succeeded on one less */
		return ans;
	      }
	    else
	      {               /* connect ftrees  */
		v->active = 0;
		
		/* if other ftree is still active then not yet added to
		   new_ftree_roots, so add it and mark inactive */ 
		if ((w=path_end_edges->head->ftree_label)->active)
		  {
		    /* note w can't be first on list (it is after v) */
		    w->prev_ftree->next_ftree = w->next_ftree;
		    if (w->next_ftree != NULL)
		      w->next_ftree->prev_ftree = w->prev_ftree;
		    w->next_ftree = new_ftree_roots;
		    if (new_ftree_roots!=NULL) new_ftree_roots->prev_ftree = w;
		    new_ftree_roots = w;
		    w->active = 0;
		  }  
	      }

	  /* do augmentations */
	  for(e2=path_end_edges; e2!=NULL;  e2 = e2->auxnext)
	    {
	      t = ans->last_tree;  /* t holds the tree to insert into */
	      for(e=e2; e != NULL; e=e->label)
		{
		  re = e->rev;
		  v = re->head;
		  w = e->head;

		  /* delete edge from wherever it is */
		  if (e->tree != NULL)
		    {
		      e->tree->was_changed = 1;
		      DeleteTreeEdge(e->tree, v, e);
		      DeleteTreeEdge(e->tree, w, re);
		    }
		  
		  /* insert edge where it belongs */
		  if (t != NULL)
		    {
		      t->was_changed = 1;
		      
		      AddTreeEdge(t, v, e);
		      AddTreeEdge(t, w, re);
		    }

		  e->tree = t;   /* update e->tree and t */
		  t = re->tree;
		  re->tree = e->tree;
		}
	      
	      e2->label = e2->rev->label = NULL;
	      e2->startlabel = e2->rev->startlabel = 0;
	    }

	  /* all created and labeled edges that aren't in path_end_edges
	     are in the q, so we use it to clear labels and reclaim memory
	     for unused edges */
	  for(e=q->auxnext; e!=NULL; e=e2) 
	    {
	      e2 = e->auxnext;     /* careful! e might get free()'ed */
	      if (e->tree == NULL)      /* clean up unused tedges */
		{
		  if ((++(e->e->cap)) == 1)
		    {
		      /* just regained capacity -> move to unused list */
		      ta = e->e;
		      v = ta->rev->head;
		      /* delete from used list */
		      if (ta->prev != NULL) ta->prev->next = ta->next;
		      else v->used_first = ta->next;
		      if (ta->next != NULL) ta->next->prev = ta->prev;
		      /* add to unused list */
		      ta->next = v->first;
		      if (v->first != NULL) v->first->prev = ta;
		      ta->prev = NULL;
		      v->first = ta;
		    }
		  recycle_node(ans->tedge_allocator, (char *)(e->rev));
		  recycle_node(ans->tedge_allocator, (char *)e);
		}
	      else                    /* clear edge labels of used edges */
		{
		  e->label = e->rev->label = NULL;
		  e->startlabel = e->rev->startlabel = 0;
		}
	    }
	  
	  q->auxnext = NULL;   /* empty q */
	  q_last = q;
	}
    }


  ans->last_tree->next_tree = NULL;
  ans->last_tree = NULL;
  ForAllTrees(ans, t) /* (re)set depths, parents, vertex labels */
    if (t->was_changed)	       
      {
	t->was_changed = 0;
	dfs_init(g->root, t, NULL, 0, NULL);
      }
  return ans;  
}

/* clean up */
void free_packing(packing p)
{
  tree *t1, *t2;
  if (p->last_tree != NULL) p->last_tree->next_tree = NULL;
  for(t2=NULL, t1 = p->first_tree; t1 != NULL; t1 = t1->next_tree)
    {
      if (t2 != NULL)
	{
	  free(t2->info);
	  free(t2);
	}
      t2 = t1;
    }
  for(t1 = p->spare_trees; t1 != NULL; t1 = t1->next_tree)
    {
      if (t2 != NULL)
	{
	  free(t2->info);
	  free(t2);
	}
      t2 = t1;
    }
  if (t2 != NULL)
    {
      free(t2->info);
      free(t2);
    }
  free_allocator((allocator)p->tedge_allocator);
  free(p);
}

/* this does a depth first search on spanning tree tree, starting at vertex,
   and sets parent_edge, depth, and label as appropriate */
static void dfs_init(node *v, tree *t, tedge *parent_edge, int depth, 
		     node *label)
{
  tedge *e;

  v->ftree_label = label;
  parent_edge(t,v) = parent_edge;
  depth(t,v) = depth++;
  ForAllTreeEdges(t, v, e)
    {
      (*edgeScanCnt)++;
      if (parent_edge == e) continue;  /* don't redo parent */

      dfs_init(e->head, t, e->rev, depth, label);
    }
}

#ifndef NDEBUG
void printtree(node *v, tree *t, node *parent)
{
  tedge *e;

  ForAllTreeEdges(t, v, e)
    {
      if (parent == e->head)
	continue;  /* don't redo parent */

      printf("%d %d\n", INDEX(e->e->rev->head), INDEX(e->e->head));
      if (e->label != NULL || e->startlabel)
	printf("\t\t\tlabelled edge!! %d %d\n", (int)e->label, e->startlabel);
      printtree(e->head, t, v);
    }
}
#endif

/* label an edge */
#ifdef __GNUC__
static inline int label(tree *t, tedge *e, tedge *label)
#else
static int label(tree *t, tedge *e, tedge *label)
#endif
{
  (*edgeScanCnt)++;
  e->label = e->rev->label = label;
  q_last->auxnext = e;
  q_last = e;
  e->auxnext = NULL;
  is_labeled(t, e->head) = 1;
  is_labeled(t, e->rev->head) = 1;
  return (e->e->rev->head->has_no_labeled_edges && 
	  label_unused_edges(e->e->rev->head, e));
}

/* label unused edges adjacent to a vertex */
static int label_unused_edges(node *v, tedge *label)
{
  tedge *e, *re;
  arc *a, *a2;

  v->has_no_labeled_edges = 0;

  a2=v->first;
  while( (a=a2) != NULL )
    {
      (*edgeScanCnt)++;
      a2 = a->next;   /* careful! a might get moved to used list */
      e = (tedge *)get_node(ans->tedge_allocator);
      re = (tedge *)get_node(ans->tedge_allocator);
      e->rev = re;
      re->rev = e;
      e->head = a->head;
      re->head = v;

      if (! (--(a->cap)))
	{
	  /* a is out of capacity -> move to used list */
	  /* delete from unused list */
	  if (a->prev != NULL) a->prev->next = a->next;
	  else v->first = a->next;
	  if (a->next != NULL) a->next->prev = a->prev;
	  /* add to used list */
	  a->next = v->used_first;
	  if (v->used_first != NULL) v->used_first->prev = a;
	  a->prev = NULL;
	  v->used_first = a;
	}

      e->e = re->e = a;      
      e->tree = re->tree = NULL;
      e->label = re->label = label;
      e->startlabel = re->startlabel = (label == NULL);
      re->auxnext = NULL;

      if (e->head->ftree_label != v->ftree_label)
	{
	  e->auxnext = path_end_edges;
	  path_end_edges = e;
	  return 1;  /* found a joining edge (whew!) */
	}

      q_last->auxnext = e;   /* add to q */
      q_last = e;
      e->auxnext = NULL;
    }
  return 0;
}

/* search for an augmenting path */
static int search(node *vertex)
{
  node *u, *r, *v;
  tedge *e, *f;
  tree *t;       /* t walks thru the trees (i in the paper) */
  tedge *fund_cycle_up;     /* unlabeled edges on fundamental cycle 
			       up from labeled subtree */
  tedge *fund_cycle_down;   /* unlabeled edges on fundamental cycle 
			       down to other end of edge */
  tedge *rev_fund, *te;

  ForAllTrees(ans, t)
    {
      t->root_of_labeled_subtree = vertex;
      is_labeled(t, vertex) = 1;
    }
  t->root_of_labeled_subtree = vertex;
  is_labeled(t, vertex) = 1;

  t = ans->last_tree;

  e = q_last;

  /* label all unused edges from vertex */
  if (label_unused_edges(vertex, NULL)) return 1;

  /* while Q is not empty... */
  for(e=e->auxnext; e!=NULL; e=e->auxnext)
    {
      if (e->tree == t)  t = t->next_tree;

      r = t->root_of_labeled_subtree;
      if (is_labeled(t,e->head))
	if (is_labeled(t,e->rev->head)) continue;
	else u = e->rev->head;
      else u = e->head;
      
      fund_cycle_up = NULL;
      fund_cycle_down = NULL;
      while(1)  /* traverse fundamental cycle:
		   walk up from unlabeled end of e and 
		   root of labeled subtree, deeper first 
		   */
	{
	  (*edgeScanCnt)++;
	  /* v = deeper of u,r */
	  if (depth(t, r) > depth(t, u)) v = r;
	  else v = u;    
	  
	  f = parent_edge(t, v);

	  if (parent_edge(t, u) == parent_edge(t, r))
	    { /* found a common edge (at lca(u, r)) 
		 label up from r to here and back down to u */
	      rev_fund = NULL;
	      while(fund_cycle_up != NULL)  /* reverse fund_cycle_up */
		{
		  te = fund_cycle_up->auxnext;
		  fund_cycle_up->auxnext = rev_fund;
		  rev_fund = fund_cycle_up;
		  fund_cycle_up = te;
		}
	      fund_cycle_up = rev_fund;
	      while (fund_cycle_up != NULL)  /* label up */
		{
		  te = fund_cycle_up->auxnext;
		  if (label(t, fund_cycle_up, e)) return 1;  
		  fund_cycle_up = te;
		}
	      while (fund_cycle_down != NULL) /* label down */
		{
		  te = fund_cycle_down->auxnext;
		  if (label(t, fund_cycle_down, e)) return 1;  
		  fund_cycle_down = te;
		}
	      t->root_of_labeled_subtree = v;
	      break;     /* stop labeling (breaks while(1)) */
	    }
	  else if (f->head->ftree_label != f->rev->head->ftree_label)
	    { /* found a joining edge */
	      f->label = f->rev->label = e;
	      
	      f->auxnext = path_end_edges;
	      path_end_edges = f;
	      return 1;
	    }
	  else if (f->label != NULL || f->startlabel)
	    { /* found a labeled edge, label down to u */
	      while (fund_cycle_down != NULL)
		{
		  te = fund_cycle_down->auxnext;
		  if (label(t, fund_cycle_down, e)) return 1;  
		  fund_cycle_down = te;
		}
	      break;     /* stop labeling (breaks while(1)) */
	    }
	  else  /* nothing doing, add this edge to the list of 
		   fundamental cycle edges, which will be labeled later */
	    {
	      if (v == r) 
		{
		  f->auxnext = fund_cycle_up;
		  fund_cycle_up = f;
		  r = f->head;
		}
	      else 
		{ 
		  f->auxnext = fund_cycle_down;
		  fund_cycle_down = f;
		  u = f->head;
		}
	    }
	}
    }
  
  return 0;  /* fall thru means no augmenting path */
}

static void dfs_init_pack(graph *g, node *root)
{
  node *v;
  arc *a, *list;
  tedge *e;
  float tm = timer();
  tree *lastT = NULL;

  ans->first_tree = NULL;  
  while (1)
    {
      /* make new tree */
      if (ans->spare_trees != NULL)
	{
	  g->currentTree = ans->spare_trees;
	  ans->spare_trees = g->currentTree->next_tree;
	}
      else
	{
	  memassert(g->currentTree= (tree *)malloc(sizeof(tree)));
	  memassert(g->currentTree->info = (struct treeinfo *)
		    malloc(g->currentN*sizeof(struct treeinfo)));
	}
      if (ans->k == 0)  /* first is special case */
	{
	  ans->last_tree = g->currentTree;;
	  ans->first_tree = ans->last_tree;
	}
      else
	{
	  lastT = ans->last_tree;
	  lastT->next_tree = g->currentTree;
	  ans->last_tree = g->currentTree;
	}
      g->currentTree->next_tree = ans->first_tree;
      g->currentTree->name = ans->k++;
      g->currentTree->was_changed = 0;
      depth(g->currentTree, root) = 0;
      parent_edge(g->currentTree, root) = NULL;
      ForAllNodes(g,v) 
	{
	  first(g->currentTree, v) = NULL;
	  parent_edge(g->currentTree, v) = NULL;
	  v->active = 1;
	}
      
      /* find a tree */
      dfs_visited = 0;
      find_dfs_tree(g, root);

      /* are we done? */
      if (dfs_visited != g->currentN)  
	{
	  /* clear last tree */
	  ans->last_tree = lastT;
	  if (lastT != NULL)
	    lastT->next_tree = ans->first_tree;
	  else
	    ans->first_tree = NULL;
	  g->currentTree->next_tree = ans->spare_trees;
	  ans->spare_trees = g->currentTree;
	  ans->k--;

	  ForAllNodes(g,v)
	    {
	      /* free edges of last tree*/
	      e = parent_edge(g->currentTree, v);
	      if (e != NULL)
		{
		  e->e->cap++;
		  recycle_node(ans->tedge_allocator, (char *)(e->rev));
		  recycle_node(ans->tedge_allocator, (char *)e);
		}
	    }
	  
	  ForAllNodes(g,v)
	    {
	      /* move 0 capacity edges to used list*/
	      list = v->first;
	      v->first = a = NULL;
	      while (list != NULL)
		{
		  a = list;
		  list = list->next;
		  if (a->cap > 0)
		    {
		      a->next = v->first;
		      if (v->first != NULL)
			v->first->prev = a;
		      v->first = a;
		    }
		  else
		    {
		      a->next = v->used_first;
		      if (v->used_first != NULL)
			v->used_first->prev = a;
		      v->used_first = a;
		    }
		  a->prev = NULL;
		}
	    }
  
	  printf("c dfs init packing found %d trees in %f time\n", ans->k, 
		 timer()-tm);
	  return;
	}
    }
}

static void find_dfs_tree(graph *g, node *v)
{
  arc *a, *ra;
  tedge *e, *re;

  v->active = 0;
  dfs_visited++;
  ForAllIncidentArcs(v, a)
    {
      if ((ra=a->rev)->cap > 0 && a->head->active)
	{
	  e = (tedge *)get_node(ans->tedge_allocator);
	  re = (tedge *)get_node(ans->tedge_allocator);
	  e->rev = re;
	  re->rev = e;
	  e->e = re->e = ra;
	  e->head = a->head;
	  re->head = v;
	  e->tree = re->tree = g->currentTree;
	  e->label = re->label = NULL;
	  e->startlabel = re->startlabel = 0;
	  e->auxnext = re->auxnext = NULL;
	  AddTreeEdge(g->currentTree, v, e);
	  AddTreeEdge(g->currentTree, a->head, re);
	  ra->cap--;
	  parent_edge(g->currentTree, a->head) = re;
	  depth(g->currentTree, a->head) = depth(g->currentTree, v)+1;
	  
	  find_dfs_tree(g, a->head);
	}
      if (dfs_visited == g->currentN) return;
    }
}
