# -*-makefile-*- for programs and generators
# $Id: makefile,v 1.18 1997/05/16 15:30:32 mslevine Exp $
######################################

CC = gcc
# uncomment -DBROKEN_HDRS is missing getopt.h
#CFLAGS = -g -Wall -DSAVECUT #-DBROKEN_HDRS
CFLAGS = -O4 -Wall -DNDEBUG -DSAVECUT #-DBROKEN_HDRS
LDFLAGS = -O4 
LOADLIBES = -lm

# place to copy binaries to
BINDIR=../bin

PROGS=ho ho_nopr ho_noxs ho_noprxs ni ni_nopr ks ks_nopr k_g

GENS=cyclegen noigen randomgen regulargen wheelgen bikewheelgen prgen \
  pr_only irregulargen dblcyclegen

all: progs gens

clean:
	rm -f $(PROGS) $(GENS) *.o dyn_tree/*.o

install: $(PROGS) $(GENS)
	cp $(PROGS) $(GENS) $(BINDIR)

# utilities
random.o: random.c random.h
	$(CC) $(CFLAGS) -c -o random.o random.c

timer.o: timer.c timer.h
	$(CC) $(CFLAGS) -c -o timer.o timer.c

graph_K_GABOW.o: graph.c graph.h memassert.h fprintfll.h timer.h 
	$(CC) $(CFLAGS) -DK -DGABOW -c graph.c -o graph_K_GABOW.o

graph_K_PST.o: graph.c graph.h memassert.h fprintfll.h timer.h 
	$(CC) $(CFLAGS) -DK -DPST -c graph.c -o graph_K_PST.o

graph_HO_nopr.o: graph.c graph.h memassert.h fprintfll.h timer.h 
	$(CC) $(CFLAGS) -DHO -DNO_PR -c graph.c -o graph_HO_nopr.o

graph_NI_nopr.o: graph.c graph.h memassert.h fprintfll.h timer.h 
	$(CC) $(CFLAGS) -DNI -DNO_PR -c graph.c -o graph_NI_nopr.o

graph_MATULA_nopr.o: graph.c graph.h memassert.h fprintfll.h timer.h 
	$(CC) $(CFLAGS) -DPR_ONLY -DNO_PR -c graph.c -o graph_MATULA_nopr.o

graph_HO.o: graph.c graph.h memassert.h fprintfll.h timer.h 
	$(CC) $(CFLAGS) -DHO -c graph.c -o graph_HO.o

graph_NI.o: graph.c graph.h memassert.h fprintfll.h timer.h 
	$(CC) $(CFLAGS) -DNI -c graph.c -o graph_NI.o

graph_MATULA.o: graph.c graph.h memassert.h fprintfll.h timer.h 
	$(CC) $(CFLAGS) -DPR_ONLY -c graph.c -o graph_MATULA.o

contract_K_GABOW.o: contract.c contract.h graph.h memassert.h fprintfll.h 
	$(CC) $(CFLAGS) -DGABOW -DK -c contract.c -o contract_K_GABOW.o

contract_K_PST.o: contract.c contract.h graph.h memassert.h fprintfll.h 
	$(CC) $(CFLAGS) -DPST -DK -c contract.c -o contract_K_PST.o

contract_HO_nopr.o: contract.c contract.h graph.h memassert.h fprintfll.h 
	$(CC) $(CFLAGS) -DHO -DNO_PR -c contract.c -o contract_HO_nopr.o

contract_NI_nopr.o: contract.c contract.h graph.h memassert.h fprintfll.h 
	$(CC) $(CFLAGS) -DNI -DNO_PR -c contract.c -o contract_NI_nopr.o

contract_MATULA_nopr.o: contract.c contract.h graph.h memassert.h fprintfll.h 
	$(CC) $(CFLAGS) -DPR_ONLY -DNO_PR -c contract.c -o contract_MATULA_nopr.o

contract_HO.o: contract.c contract.h graph.h memassert.h fprintfll.h 
	$(CC) $(CFLAGS) -DHO -c contract.c -o contract_HO.o

contract_NI.o: contract.c contract.h graph.h memassert.h fprintfll.h 
	$(CC) $(CFLAGS) -DNI -c contract.c -o contract_NI.o

contract_MATULA.o: contract.c contract.h graph.h memassert.h fprintfll.h 
	$(CC) $(CFLAGS) -DPR_ONLY -c contract.c -o contract_MATULA.o

ni_sparse_cert_K_GABOW.o:  ni_sparse_cert.c ni_sparse_cert.h graph.h \
 memassert.h fprintfll.h heap.h timer.h pr.h 
	$(CC) $(CFLAGS) -DK -DGABOW -c ni_sparse_cert.c -o ni_sparse_cert_K_GABOW.o

ni_sparse_cert_K_PST.o:  ni_sparse_cert.c ni_sparse_cert.h graph.h \
 memassert.h fprintfll.h heap.h timer.h pr.h 
	$(CC) $(CFLAGS) -DK -DPST -c ni_sparse_cert.c -o ni_sparse_cert_K_PST.o

ni_sparse_cert_HO_nopr.o:  ni_sparse_cert.c ni_sparse_cert.h graph.h \
 memassert.h fprintfll.h heap.h timer.h pr.h contract.h
	$(CC) $(CFLAGS) -DHO -DNO_PR -c ni_sparse_cert.c -o ni_sparse_cert_HO_nopr.o

ni_sparse_cert_NI_nopr.o:  ni_sparse_cert.c ni_sparse_cert.h graph.h \
 memassert.h fprintfll.h heap.h timer.h pr.h contract.h
	$(CC) $(CFLAGS) -DNI -DNO_PR -c ni_sparse_cert.c -o ni_sparse_cert_NI_nopr.o

ni_sparse_cert_MATULA_nopr.o:  ni_sparse_cert.c ni_sparse_cert.h graph.h \
 memassert.h fprintfll.h heap.h timer.h pr.h contract.h
	$(CC) $(CFLAGS) -DPR_ONLY -DNO_PR -c ni_sparse_cert.c -o ni_sparse_cert_MATULA_nopr.o

ni_sparse_cert_HO.o:  ni_sparse_cert.c ni_sparse_cert.h graph.h \
 memassert.h fprintfll.h heap.h timer.h pr.h contract.h
	$(CC) $(CFLAGS) -DHO -c ni_sparse_cert.c -o ni_sparse_cert_HO.o

ni_sparse_cert_NI.o:  ni_sparse_cert.c ni_sparse_cert.h graph.h \
 memassert.h fprintfll.h heap.h timer.h pr.h contract.h
	$(CC) $(CFLAGS) -DNI -c ni_sparse_cert.c -o ni_sparse_cert_NI.o

ni_sparse_cert_MATULA.o:  ni_sparse_cert.c ni_sparse_cert.h graph.h \
 memassert.h fprintfll.h heap.h timer.h pr.h contract.h
	$(CC) $(CFLAGS) -DPR_ONLY -c ni_sparse_cert.c -o ni_sparse_cert_MATULA.o

pr_K_GABOW.o: pr.c pr.h graph.h memassert.h fprintfll.h heap.h timer.h contract.h
	$(CC) $(CFLAGS) -DK -DGABOW -c pr.c -o pr_K_GABOW.o

pr_K_PST.o: pr.c pr.h graph.h memassert.h fprintfll.h heap.h timer.h contract.h
	$(CC) $(CFLAGS) -DK -DPST -c pr.c -o pr_K_PST.o

pr_HO_INTERNAL.o: pr.c pr.h graph.h memassert.h fprintfll.h heap.h timer.h \
 contract.h
	$(CC) $(CFLAGS) -DHO_INTERNAL -c pr.c -o pr_HO_INTERNAL.o

pr_PR_ONLY.o: pr.c pr.h graph.h memassert.h fprintfll.h heap.h timer.h \
 contract.h
	$(CC) $(CFLAGS) -DPR_ONLY -c pr.c -o pr_PR_ONLY.o

pr_HO.o: pr.c pr.h graph.h memassert.h fprintfll.h heap.h timer.h contract.h
	$(CC) $(CFLAGS) -DHO -c pr.c -o pr_HO.o

pr_NI.o: pr.c pr.h graph.h memassert.h fprintfll.h heap.h timer.h contract.h
	$(CC) $(CFLAGS) -DNI -c pr.c -o pr_NI.o

pr_MATULA.o: pr.c pr.h graph.h memassert.h fprintfll.h heap.h timer.h contract.h
	$(CC) $(CFLAGS) -DPR_ONLY -c pr.c -o pr_MATULA.o

# problem generators
gens: $(GENS)

sgbgen: sgbgen.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(SGBINC) -o sgbgen  sgbgen.c $(LOADLIBES) $(SGBLIB) -lgb

noigen: noigen.c random.o
	$(CC) $(CFLAGS) -o noigen  noigen.c random.o

cyclegen: cyclegen.c random.o
	$(CC) $(CFLAGS) -o cyclegen  cyclegen.c random.o

randomgen: randomgen.c random.o
	$(CC) $(CFLAGS) -o randomgen  randomgen.c random.o

regulargen: regulargen.c random.o
	$(CC) $(CFLAGS) -o regulargen  regulargen.c random.o

irregulargen: irregulargen.c random.o
	$(CC) $(CFLAGS) -o irregulargen  irregulargen.c random.o

wheelgen: wheelgen.c random.o
	$(CC) $(CFLAGS) -o wheelgen  wheelgen.c random.o

bikewheelgen: bikewheelgen.c random.o
	$(CC) $(CFLAGS) -o bikewheelgen  bikewheelgen.c random.o

prgen: prgen.c random.o
	$(CC) $(CFLAGS) -o prgen  prgen.c random.o

dblcyclegen: dblcyclegen.c random.o
	$(CC) $(CFLAGS) -o dblcyclegen  dblcyclegen.c random.o

# programs
progs: $(PROGS)

###### HO

ho: ho.o random.o timer.o graph_HO.o pr_HO.o pr_HO_INTERNAL.o contract_HO.o \
fprintfll.o 
	$(CC) $(LDFLAGS)  ho.o random.o timer.o graph_HO.o pr_HO.o pr_HO_INTERNAL.o contract_HO.o fprintfll.o $(LOADLIBES) -o ho

ho_noint: ho.o random.o timer.o graph_HO.o pr_HO.o contract_HO.o fprintfll.o 
	$(CC) $(LDFLAGS)  ho.o random.o timer.o graph_HO.o pr_HO.o contract_HO.o fprintfll.o $(LOADLIBES) -o ho_noint

ho_nopr: ho_nopr.o random.o timer.o graph_HO_nopr.o fprintfll.o
	$(CC) $(LDFLAGS)  ho_nopr.o random.o timer.o graph_HO_nopr.o fprintfll.o $(LOADLIBES) -o ho_nopr

ho_noxs: ho_noxs.o random.o timer.o graph_HO.o pr_HO.o pr_HO_INTERNAL.o \
  fprintfll.o contract_HO.o
	$(CC) $(LDFLAGS)  ho_noxs.o random.o timer.o graph_HO.o pr_HO.o pr_HO_INTERNAL.o fprintfll.o contract_HO.o $(LOADLIBES) -o ho_noxs

ho_noprxs: ho_noprxs.o random.o timer.o graph_HO_nopr.o fprintfll.o
	$(CC) $(LDFLAGS)  ho_noprxs.o random.o timer.o graph_HO_nopr.o fprintfll.o $(LOADLIBES) -o ho_noprxs

ho.o: ho.c graph.h memassert.h fprintfll.h pr.h contract.h timer.h \
 ho.h qs.h  
	$(CC)  $(CFLAGS) -c ho.c -o ho.o

ho_noint.o: ho.c graph.h memassert.h fprintfll.h pr.h contract.h timer.h \
 ho.h qs.h  
	$(CC)  $(CFLAGS) -DDONT_DO_INTERNAL -c ho.c -o ho_noint.o

ho_nopr.o: ho.c graph.h memassert.h fprintfll.h timer.h ho.h qs.h  
	$(CC)  $(CFLAGS) -DNO_PR -c ho.c -o ho_nopr.o

ho_noxs.o: ho.c graph.h memassert.h fprintfll.h pr.h contract.h timer.h \
 ho.h qs.h  
	$(CC)  $(CFLAGS) -DNO_EXCESS_DETECTION -c ho.c -o ho_noxs.o

ho_noprxs.o: ho.c graph.h memassert.h fprintfll.h timer.h ho.h qs.h   
	$(CC)  $(CFLAGS) -DNO_PR -DNO_EXCESS_DETECTION -c ho.c -o ho_noprxs.o

###### NI

ni: ni.o random.o timer.o graph_NI.o pr_NI.o ni_sparse_cert_NI.o \
 fprintfll.o contract_NI.o
	$(CC) $(LDFLAGS)  ni.o random.o timer.o graph_NI.o pr_NI.o ni_sparse_cert_NI.o fprintfll.o contract_NI.o $(LOADLIBES) -o ni

ni_nopr: ni_nopr.o random.o timer.o graph_NI_nopr.o fprintfll.o \
 ni_sparse_cert_NI_nopr.o contract_NI_nopr.o
	$(CC) $(LDFLAGS)  ni_nopr.o random.o timer.o graph_NI_nopr.o fprintfll.o ni_sparse_cert_NI_nopr.o contract_NI_nopr.o $(LOADLIBES) -o ni_nopr

ni.o: ni.c graph.h memassert.h fprintfll.h pr.h heap.h timer.h \
 ni_sparse_cert.h  
	$(CC)  $(CFLAGS) -c ni.c -o ni.o

ni_nopr.o: ni.c graph.h memassert.h fprintfll.h heap.h timer.h \
 ni_sparse_cert.h  
	$(CC) $(CFLAGS) -DNO_PR -o ni_nopr.o -c ni.c

hybrid_main: hybrid.o timer.o
	$(CC) $(LDFLAGS)  hybrid.o timer.o -lF77 -lI77 $(LOADLIBES) -o hybrid_main

hybrid.o: hybrid.c
	$(CC) $(CFLAGS) -c hybrid.c -o hybrid.o

hybrid.c: hybrid.F
	f2c -A hybrid.F

##

matula: matula.o random.o timer.o contract_MATULA.o graph_MATULA.o \
 pr_MATULA.o fprintfll.o ni_sparse_cert_MATULA.o
	$(CC) $(LDFLAGS)  matula.o random.o timer.o contract_MATULA.o graph_MATULA.o pr_MATULA.o fprintfll.o ni_sparse_cert_MATULA.o $(LOADLIBES) -o matula

matula.o: matula.c graph.h memassert.h fprintfll.h contract.h pr.h \
 heap.h timer.h ni_sparse_cert.h  
	$(CC)  $(CFLAGS) -c matula.c -o matula.o

###### KS

ks: ks.o random.o timer.o graph_NI.o pr_NI.o fprintfll.o contract_NI.o
	$(CC) $(LDFLAGS)  ks.o random.o timer.o graph_NI.o pr_NI.o fprintfll.o contract_NI.o $(LOADLIBES) -o ks

ks_nopr: ks_nopr.o random.o timer.o graph_NI.o fprintfll.o pr_NI.o \
 contract_NI.o
	$(CC) $(LDFLAGS)  ks_nopr.o random.o timer.o graph_NI.o fprintfll.o pr_NI.o contract_NI.o $(LOADLIBES) -o ks_nopr

ks_pr34: ks_pr34.o random.o timer.o graph_NI.o pr_NI.o fprintfll.o \
 contract_NI.o
	$(CC) $(LDFLAGS)  ks_pr34.o random.o timer.o graph_NI.o pr_NI.o fprintfll.o contract_NI.o $(LOADLIBES) -o ks_pr34

ks.o: ks.c random.h ks.h graph.h memassert.h fprintfll.h pr.h heap.h \
 timer.h 
	$(CC)  $(CFLAGS) -DDEBUG=1 -DDEPTHTIMES -c ks.c -o ks.o

ks_nopr.o: ks.c random.h ks.h graph.h memassert.h fprintfll.h timer.h 
	$(CC)  $(CFLAGS) -DDEBUG=1 -DDEPTHTIMES -DNO_PR -c ks.c -o ks_nopr.o

ks_pr34.o: ks.c random.h ks.h graph.h memassert.h fprintfll.h pr.h heap.h \
 timer.h 
	$(CC)  $(CFLAGS) -DDEBUG=1 -DDEPTHTIMES -DPR34 -c ks.c -o ks_pr34.o

###### K 

k_g: k__GABOW.o karger_mincut__GABOW.o gabow_pack.o allocate.o random.o \
 timer.o graph_K_GABOW.o pr_K_GABOW.o dyn_tree/dyn_splay__GABOW.o \
 dyn_tree/dyn_tree__GABOW.o ni_sparse_cert_K_GABOW.o \
 contract_K_GABOW.o fprintfll.o
	$(CC) $(LDFLAGS)  k__GABOW.o karger_mincut__GABOW.o gabow_pack.o allocate.o random.o timer.o graph_K_GABOW.o pr_K_GABOW.o dyn_tree/dyn_splay__GABOW.o dyn_tree/dyn_tree__GABOW.o ni_sparse_cert_K_GABOW.o contract_K_GABOW.o fprintfll.o $(LOADLIBES) -o k_g

k_pst: k__PST.o karger_mincut__PST.o pst.o allocate.o random.o timer.o \
 graph_K_PST.o pr_K_PST.o  fprintfll.o
	$(CC) $(LDFLAGS)  k__PST.o karger_mincut__PST.o pst.o allocate.o random.o timer.o graph_K_PST.o pr_K_PST.o fprintfll.o $(LOADLIBES) -o k_pst

k__GABOW.o: k.c karger_mincut.h graph.h memassert.h fprintfll.h pr.h \
 contract.h random.h timer.h
	$(CC)  $(CFLAGS) -DGABOW -DVERBOSE -c k.c -o k__GABOW.o

k__PST.o: k.c karger_mincut.h graph.h memassert.h fprintfll.h pr.h \
 contract.h random.h timer.h
	$(CC)  $(CFLAGS) -DPST -DVERBOSE -c k.c -o k__PST.o

allocate.o: allocate.c allocate.h
	$(CC)  $(CFLAGS) -c allocate.c -o allocate.o

gabow_pack.o:  gabow_pack.c gabow_pack.h graph.h memassert.h \
 fprintfll.h allocate.h timer.h 
	$(CC)  $(CFLAGS) -c gabow_pack.c -o gabow_pack.o

pst.o: pst.c pst.h  graph.h memassert.h heap.h random.h 
	$(CC)  $(CFLAGS) -DVERBOSE -c pst.c -o pst.o

dyn_tree/dyn_splay__GABOW.o: dyn_tree/dyn_splay.c dyn_tree/_dyn_tree.h \
 dyn_tree/dyn_tree.h  fprintfll.h dyn_tree/macros.h graph.h
	$(CC)  $(CFLAGS) -DGABOW -c dyn_tree/dyn_splay.c -o dyn_tree/dyn_splay__GABOW.o

dyn_tree/dyn_splay__PST.o: dyn_tree/dyn_splay.c dyn_tree/_dyn_tree.h \
 dyn_tree/dyn_tree.h  fprintfll.h dyn_tree/macros.h graph.h
	$(CC)  $(CFLAGS) -DPST -c dyn_tree/dyn_splay.c -o dyn_tree/dyn_splay__PST.o

dyn_tree/dyn_tree__GABOW.o: dyn_tree/dyn_tree.c dyn_tree/_dyn_tree.h \
 dyn_tree/dyn_tree.h  fprintfll.h dyn_tree/macros.h allocate.h \
 memassert.h  graph.h
	$(CC)  $(CFLAGS) -DGABOW -c dyn_tree/dyn_tree.c -o dyn_tree/dyn_tree__GABOW.o

dyn_tree/dyn_tree__PST.o: dyn_tree/dyn_tree.c dyn_tree/_dyn_tree.h \
 dyn_tree/dyn_tree.h  fprintfll.h dyn_tree/macros.h allocate.h \
 memassert.h  graph.h
	$(CC)  $(CFLAGS) -DPST -c dyn_tree/dyn_tree.c -o dyn_tree/dyn_tree__PST.o

karger_mincut__GABOW.o: karger_mincut.c karger_mincut.h graph.h memassert.h \
 fprintfll.h pr.h contract.h random.h timer.h heap.h \
 dyn_tree/dyn_tree.h dyn_tree/../graph.h ni_sparse_cert.h 
	$(CC)  $(CFLAGS) -O -DVERBOSE -DGABOW -DCHEAT -c karger_mincut.c -o karger_mincut__GABOW.o

karger_mincut__PST.o: karger_mincut.c karger_mincut.h graph.h memassert.h \
 fprintfll.h pr.h contract.h random.h timer.h heap.h \
 dyn_tree/dyn_tree.h dyn_tree/../graph.h ni_sparse_cert.h 
	$(CC)  $(CFLAGS) -O -DVERBOSE -DPST -DCHEAT -c karger_mincut.c -o karger_mincut__PST.o

###### PR_ONLY

pr_only: pr_only.o random.o timer.o graph_NI.o pr_PR_ONLY.o contract_NI.o \
 fprintfll.o
	$(CC) $(LDFLAGS)  pr_only.o random.o timer.o graph_NI.o pr_PR_ONLY.o contract_NI.o fprintfll.o $(LOADLIBES) -o pr_only

pr_only.o:  pr_only.c pr.h graph.h memassert.h fprintfll.h heap.h \
 timer.h 	
	$(CC) $(CFLAGS) -c pr_only.c -o pr_only.o

