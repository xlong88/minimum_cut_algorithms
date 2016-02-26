/* General purpose graph functions
 * $Id: graph.c,v 1.4 1997/05/16 15:30:31 mslevine Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph.h"
#include "timer.h"
#include "contract.h"

/*********************************************** function definitions */

graph *makeGraph(int n, int m)
{
  graph *g;

  memassert(g = (graph *)calloc(1, sizeof(graph)));
  
  g->n = g->currentN = n;
  g->m = m; g->currentM = 2*m;

  memassert(g->vertexalloc    = (node*)calloc(g->n+1, sizeof(node)));
  g->vertices = g->vertexalloc;
  g->sentinelVertex = g->vertices+n;
#ifdef CONTRACT
  memassert(g->nodes = (node**)calloc(g->n+1, sizeof(node *)));
#ifdef ENABLE_GRAPH_BACKUPS
  memassert(g->currentArcs = (arc**)calloc(g->currentM+1, sizeof(arc *)));
#endif
#endif
  memassert(g->arcs = (arc*)calloc( 2*g->m+1, sizeof(arc) ));
  g->sentinelArc = g->arcs + 2*g->m;

  g->minCap = MAXWEIGHT;
  g->cutCount = 0;
  g->context = MAIN;
  return g;
}

long totalEdgeScans(graph *g)
{
  int i;
  long ans=0;

  for(i=0; i<CONTEXTS; i++)
    ans += g->edgeScanCnt[i];
  return ans;
}


#ifdef SAVECUT
void printCut(graph *g)
{
  node *w;
  int count=0, in=1;

  printf("c one side of the mincut\n");
  ForAllVertices(g, w)
    if ( w -> in_cut == 0 ) count++;
  if (count < g->n/2) in=0;
  ForAllVertices(g, w)
    if ( w -> in_cut == in )   
      printf("n %d\n", VERTEX_NAME(g, w));
}
#endif

void saveTrivialCut(graph *g, node *v, weight_t val)
{
#ifdef SAVECUT
  node *w;
#endif

  g->minCap = val;
  g->dtime = timer () ;
  g->cutCount++;

#ifdef VERBOSE
  printf("c new trivial cut: ");
  fprintf_wt(stdout, "%", g->minCap);
  printf("\n");
#endif

#ifdef SAVECUT
#ifdef K
  g->save1 = NULL;
#endif
#ifdef CONTRACT
  v->in_cut = 1;
  ForAllVertices(g, w) 
    if (findLeader(w) == v)
      w->in_cut = 1;
    else
      w->in_cut = 0;
#else
  ForAllVertices(g,w)
    w->in_cut = 0;
  v->in_cut = 1;
#endif /* CONTRACT */
#endif /* SAVECUT */
}

void freeGraph(graph *g)
{
  free(g->vertexalloc);
  free(g->arcs);
#ifdef CONTRACT
  free(g->nodes);
#endif
  free(g);
}

/*
   parse (...) :
   1. Reads maximal flow problem in extended DIMACS format.
   2. Prepares internal data representation.
   
   types: 'arc' and 'node' must be predefined
   
   type arc  must contain fields 'head', 'rev', 'next', CAPFIELD
   
   typedef 
   struct arcSt
   {
   weight_t         resCap;     .. residual capacity
   struct nodeSt   *head;     .. head node 
   struct arcSt    *rev;   .. opposite arc 
   struct arcSt    *next;     .. next arc with the same tail 
   ..................
   }
   arc;
   
   type   node   must contain the field 'first': 
   
   typedef
   struct nodeSt
   {
   arcSt        *first;    ..  first outgoing arc 
   ....................
   }
   node;
   */

/* ----------------------------------------------------------------- */


graph  *dimacsParse(FILE *fp)
{

#define MAXLINE       100	/* max line length in the input file */
#define ARC_FIELDS      3	/* no of fields in arc line  */
#define NODE_FIELDS     2	/* no of fields in node line  */
#define P_FIELDS        3       /* no of fields in problem line */
#define PROBLEM_TYPE "cut"      /* name of problem type*/

  int n,m;                      /* number of vertices, arcs */
  graph *g=NULL;                /* the return value */
  int node_min=INT_MAX,		/* minimal no of node  */
      node_max=0,		/* maximal no of nodes */
     *arc_first=NULL,		/* internal array for holding
				   - node degree
				   - position of the first outgoing arc */
      *arc_tail=NULL,		/* internal array: tails of the arcs */
  /* temporary variables carrying no of nodes */
      head, tail, i,

  /* temporary variables carrying no of arcs */
       last, arc_num, arc_new_num;

  node   *head_p, *ndp; /* pointers to the node structure */

  arc    *arc_current=NULL, /* pointer to the arc structure */
         *arc_new,
         *arc_tmp;

  weight_t cap;			/* capacity of the current arc */
#ifdef LONGLONGINTWEIGHTS
  long capr;                    /* we read long longs as longs */
#else
  weight_t capr;                /* we read most types as themselves */
#endif 

  int    no_lines=0,		/* no of current input line */
  no_plines=0,			/* no of problem-lines */
  /*no_nslines=0,            no of node-source-lines */
  /*no_nklines=0,            no of node-source-lines */
  no_alines=0,			/* no of arc-lines */
  pos_current=0;		/* 2*no_alines */

  char    in_line[MAXLINE],	/* for reading input line */
    pr_type[3];			/* for reading type of the problem */
  /*nd;                      source (s) or sink (t) */

  int     err_no;		/* no of detected error */

  /* -------------- error numbers & error messages ---------------- */
#define EN1   0
#define EN2   1
#define EN3   2
#define EN4   3
#define EN6   4
#define EN10  5
#define EN7   6
#define EN8   7
#define EN9   8
#define EN11  9
#define EN12 10
#define EN13 11
#define EN14 12
#define EN16 13
#define EN15 14
#define EN17 15
#define EN18 16
#define EN21 17
#define EN19 18
#define EN20 19
#define EN22 20

  static char *err_message[] = 
    { 
      /* 0*/    "more than one problem line.",
      /* 1*/    "wrong number of parameters in the problem line.",
      /* 2*/    "it is not a Max Flow problem line.",
      /* 3*/    "bad value of a parameter in the problem line.",
      /* 4*/    "can't obtain enough memory to solve this problem.",
      /* 5*/    "more than one line with the problem name.",
      /* 6*/    "can't read problem name.",
      /* 7*/    "problem description must be before arc description.",
      /* 8*/    "this parser doesn't support multiply sources and sinks.",
      /* 9*/    "wrong number of parameters in the node line.",
      /*10*/    "wrong value of parameters in the node line.",
      /*11*/    " ",
      /*12*/    "source and sink descriptions must be before arc descriptions.",
      /*13*/    "too many arcs in the input.",
      /*14*/    "wrong number of parameters in the arc line.",
      /*15*/    "wrong value of parameters in the arc line.",
      /*16*/    "unknown line type in the input.",
      /*17*/    "reading error.",
      /*18*/    "not enough arcs in the input.",
      /*19*/    "source or sink doesn't have incident arcs.",
      /*20*/    "can't read anything from the input file."
	      };
  /* --------------------------------------------------------------- */

  /* The main loop:
     -  reads the line of the input,
     -  analyses its type,
     -  checks correctness of parameters,
     -  puts data to the arrays,
     -  does service functions
     */

  while ( fgets ( in_line, MAXLINE, fp ) != NULL )
    {
      no_lines ++;

      switch (in_line[0])
	{
	case 'c':		/* skip lines with comments */
	case '\n':		/* skip empty lines   */
	case '\0':		/* skip empty lines at the end of file */
	  break;

	case 'p':		/* problem description      */
	  if ( no_plines > 0 )
	    /* more than one problem line */
	    { err_no = EN1 ; goto error; }

	  no_plines = 1;
   
	  if (
	      /* reading problem line: 
		 type of problem, no of vertices, no of arcs */
	      sscanf ( in_line, "%*c %3s %d %d", pr_type, &n, &m )
	      != P_FIELDS
	      )
	    /*wrong number of parameters in the problem line*/
	    { err_no = EN2; goto error; }

	  if ( strcmp ( pr_type, PROBLEM_TYPE ) )
	    /*wrong problem type*/
	    { err_no = EN3; goto error; }

	  if ( n <= 0  || m <= 0 )
	    /*wrong value of no of arcs or vertices*/
	    { err_no = EN4; goto error; }

	  g = makeGraph(n+1, m);
	  /* allocating memory for 'vertices', 'arcs'  and internal arrays */
	  memassert(arc_tail = (int*) calloc ( 2*g->m,   sizeof(int) )); 
	  memassert(arc_first= (int*) calloc ( g->n+1, sizeof(int) ));
	  /* arc_first [ 0 .. n+1 ] = 0 - initialized by calloc */

	  if ( g->vertexalloc == NULL || g->arcs == NULL || 
	      arc_first == NULL || arc_tail == NULL )
	    /* memory is not allocated */
	    { err_no = EN6; goto error; }
		     
	  /* setting pointer to the first arc */
	  arc_current = g->arcs;

	  node_max = 0;
	  node_min = g->n;

	  break;

	case 'a':		/* arc description */
	  if ( no_plines < 1)
	    { err_no = EN7; goto error; }
		
	  if ( no_alines >= g->m )
	    /*too many arcs on input*/
	    { err_no = EN16; goto error; }
		
	  if (
	      /* reading an arc description */
	      sscanf ( in_line, "%*c %d %d %" WT_RD_FORMAT, 
		      &tail, &head, &capr )
	      != ARC_FIELDS
	      ) 
	    /* arc description is not correct */
	    { err_no = EN15; goto error; }

	  if ( tail < 0  ||  tail > n  ||
	      head < 0  ||  head > n  
	      )
	    /* wrong value of nodes */
	    { err_no = EN17; goto error; }

	  /* no of arcs incident to node i is stored in arc_first[i+1] */
	  arc_first[tail + 1] ++; 
	  arc_first[head + 1] ++;

	  /* storing information about the arc */
	  arc_tail[pos_current]        = tail;
	  arc_tail[pos_current+1]      = head;
	  arc_current       -> head    = g->vertexalloc + head;
	  arc_current       -> CAPFIELD    = capr;
	  arc_current       -> rev  = arc_current + 1;
	  ( arc_current + 1 ) -> head    = g->vertexalloc + tail;
	  ( arc_current + 1 ) -> CAPFIELD    = capr;
	  ( arc_current + 1 ) -> rev  = arc_current;

	  /* update minimum and maximum node */
	  if ( head < node_min ) node_min = head;
	  if ( tail < node_min ) node_min = tail;
	  if ( head > node_max ) node_max = head;
	  if ( tail > node_max ) node_max = tail;

	  no_alines   ++;
	  arc_current += 2;
	  pos_current += 2;
	  break;

	default:
	  /* unknown type of line */
	  err_no = EN18; goto error;
	  break;

	} /* end of switch */
    } /* end of input loop */

  /* ----- all is read  or  error while reading ----- */ 

  if ( feof (stdin) == 0 )	/* reading error */
    { err_no=EN21; goto error; } 

  if ( no_lines == 0 )		/* empty input */
    { err_no = EN22; goto error; } 

  if ( no_alines < m )		/* not enough arcs */
    { err_no = EN19; goto error; } 

  /********** ordering arcs - linear time algorithm ***********/

  /* first arc from the first node */
  ( g->vertexalloc + node_min ) -> first = g->arcs;

  /* before below loop arc_first[i+1] is the number of arcs outgoing from i;
     after this loop arc_first[i] is the position of the first 
     outgoing from node i arcs after they would be ordered;
     this value is transformed to pointer and written to node.first[i]
     */
 
  for ( i = node_min + 1; i <= node_max + 1; i ++ ) 
    {
      arc_first[i]          += arc_first[i-1];
      ( g->vertexalloc + i ) -> first = g->arcs + arc_first[i];
    }


  for ( i = node_min; i < node_max; i ++ ) /* scanning all the vertices  
					      exept the last*/
    {

      last = ( ( g->vertexalloc + i + 1 ) -> first ) - g->arcs;
      /* arcs outgoing from i must be cited    
	 from position arc_first[i] to the position
	 equal to initial value of arc_first[i+1]-1  */

      for ( arc_num = arc_first[i]; arc_num < last; arc_num ++ )
	{ tail = arc_tail[arc_num];

	  while ( tail != i )
	    /* the arc no  arc_num  is not in place because arc cited here
	       must go out from i;
	       we'll put it to its place and continue this process
	       until an arc in this position would go out from i */

	    { arc_new_num  = arc_first[tail];
	      arc_current  = g->arcs + arc_num;
	      arc_new      = g->arcs + arc_new_num;
	    
	      /* arc_current must be cited in the position arc_new    
		 swapping these arcs:                                 */

	      head_p               = arc_new -> head;
	      arc_new -> head      = arc_current -> head;
	      arc_current -> head  = head_p;

	      cap                 = arc_new -> CAPFIELD;
	      arc_new -> CAPFIELD     = arc_current -> CAPFIELD;
	      arc_current -> CAPFIELD = cap;

	      if ( arc_new != arc_current -> rev )
		{
		  arc_tmp                = arc_new -> rev;
		  arc_new  -> rev     = arc_current -> rev;
		  arc_current -> rev  = arc_tmp;

		  ( arc_current -> rev ) -> rev = arc_current;
		  ( arc_new     -> rev ) -> rev = arc_new;
		}

	      arc_tail[arc_num] = arc_tail[arc_new_num];
	      arc_tail[arc_new_num] = tail;

	      /* we increase arc_first[tail]  */
	      arc_first[tail] ++ ;

	      tail = arc_tail[arc_num];
	    }
	}
      /* all arcs outgoing from  i  are in place */
    }       

  /* -----------------------  arcs are ordered  ------------------------- */

  /*----------- constructing lists ---------------*/


  for ( ndp=g->vertexalloc+node_min; ndp <= g->vertexalloc + node_max; ndp++)
    ndp -> first = (arc*) NULL;

  for (arc_current =g->arcs+(2*g->m-1); arc_current >= g->arcs; arc_current--)
    {
      arc_num = arc_current - g->arcs;
      tail = arc_tail [arc_num];
      ndp = g->vertexalloc + tail;
      arc_current -> next = ndp -> first;
      ndp -> first = arc_current;
    }


  /* ----------- assigning output values ------------*/
  g->currentN = g->n = node_max - node_min + 1;
  g->vertices = g->vertexalloc + node_min;
  g->sentinelVertex = g->vertices + g->n;
  g->minCap = MAXWEIGHT;

  /* free internal memory */
  free ( arc_first ); free ( arc_tail );

#ifdef CONTRACT
  /* set v->last and a->prev pointers */
  ForAllVertices(g, ndp) 
    if (ndp->first != NULL)
      {
	ndp->first->prev = NULL;
	for (arc_tmp=ndp->first; arc_tmp->next!=NULL; arc_tmp=arc_tmp->next)
	  arc_tmp->next->prev = arc_tmp;
	ndp->last = arc_tmp;
      }
    else
      ndp->last = NULL;

  for(i=0; i<g->n; i++)
    {
      g->nodes[i] = g->vertices+i;
      g->vertices[i].leader = g->vertices+i;
      g->vertices[i].index = i;
    }
  g->nodes[g->n] = NULL;
#ifdef ENABLE_GRAPH_BACKUPS
  for(i=0; i<g->currentM; i++)
    {
      g->currentArcs[i] = g->arcs+i;
      g->arcs[i].index = i;
    }
  g->currentArcs[g->currentM] = NULL;
#endif /* ENABLE_GRAPH_BACKUPS */
#endif /* CONTRACT */

  /* Thanks God! all is done */
  return g;

  /* ---------------------------------- */
 error:				/* error found reading input */

  printf ( "\nline %d of input - %s\n", 
	  no_lines, err_message[err_no] );

  exit(1);
}
/* --------------------   end of parser  -------------------*/


#ifndef CONTRACT

/* we don't have the fields we want, so we'll (ab)use some others */
#ifdef HO
#define AUXARC current
#define VCAP vCap
#endif
#ifdef K
#define AUXARC s_first
#define VCAP vCap
#endif

void compact(graph *g)
{ 
  node *v, *w;
  arc *a, *b, *prevA;
#ifdef HO
  weight_t vCap;
#endif

  ForAllNodes (g, v ) 
    {
      ForAllIncidentArcs ( v, a ) 
	if (a->head != NULL) a->head->AUXARC = NULL;
      
      VCAP = 0;
      prevA = NULL;
      ForAllIncidentArcs ( v, a ) 
	{
	  g->edgeScanCnt[g->context]++;
	  w = a->head;
	  if ( w == v || w== NULL )  /* delete self-loops, 
					arcs marked for deletion */
	    if ( prevA != NULL )  prevA->next = a->next;
	    else v->first = a->next;
	  else 
	    {
	      VCAP += a->CAPFIELD;
	      if (w->AUXARC == NULL) 
		{
		  w->AUXARC = a;
		  prevA = a;
		}
	      else 
		{
		  b = w->AUXARC;
		  b->CAPFIELD += a->CAPFIELD;
		  Reverse(b)->CAPFIELD +=  Reverse (a)->CAPFIELD;
		  Reverse(a)->head = NULL; /* mark for deletion */
		  if ( prevA != NULL ) prevA->next = a->next;
		  else v->first = a->next;
		}
	    }
	}
  
      if ( VCAP < g->minCap ) 
	saveTrivialCut(g, v, VCAP);
    }
}
#endif /* ! CONTRACT */

