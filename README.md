# Implementations of Minimum Cut Algorithms

## Description

This package contains implementations of four recent algorithms for
the minimum cut problem and associated problem generators, as
described in "[**Experimental Study of Minimum Cut Algorithms**](http://dl.acm.org/citation.cfm?id=314161.314315)", in
proceedings of the 8th Annual Symposium on Discrete Algorithms (SODA), 1997.


The cut programs all take graphs in DIMACS format on stdin and output
the value of the minimum cut, the time taken to find it, and various
stats on stdout. It shouldn't be hard to modify the code to
output the actual cut. (In some cases the code for it is there, but it
is untested.) See the README in the source directory for more info.

## Authors

* [Chandra Chekuri](http://theory.stanford.edu/people/chekuri)		
* [Andrew Goldberg](http://www.neci.nj.nec.com/homepages/avg.html)	
* [David Karger](http://theory.lcs.mit.edu/~karger)	
* [Matthew Levine](http://theory.lcs.mit.edu/~mslevine)
* [Cliff Stein](http://www.cs.dartmouth.edu/~cliff)	
			
## Files/Subdirectories

~~~
README		this file
COPYRIGHT	copying policy
bin		dir for compiled binaries
inputs		some TSP inputs courtesy of Applegate and Cook
results		dir for results
scripts		dir with scripts for testing codes
src		dir with C source code
src/dyn_tree	dir with Dynamic Tree code from Tamas Badics
~~~

## Install

The makefile in the src subdirectory builds the 4 programs with
various options, as well as the generators. The command 'make install'
will make all the binaries and copy them to the bin subdirectory.


The codes take advantage of features present in gcc, but they should
compile with any ANSI C compiler. If your system doesn't have getopt.h
(eg, it is a Sun), you'll need to uncomment -DBROKEN_HDRS on the
CFLAGS line in the makefile in order to get some of the programs to
compile.


See the [README](./src/README) in the src directory for more info.

## Availability

This package is available at [http://www.columbia.edu/~cs2035/code.html](http://www.columbia.edu/~cs2035/code.html).

It is for **NON-COMMERCIAL USE ONLY**.

## Notice

This package **COPIED** from [http://www.columbia.edu/~cs2035/code.html](http://www.columbia.edu/~cs2035/code.html).

