# DSSSP

This repository stores the source code whe ran through our compiled executions for the project present on the paper "LP-based Heuristics for the Distinguishing Substring Selection Problem"[1]. 

General information about the dataset we have artificially generated for the Distinguishing String Selection Problem (DSSP) and the Distinguishing Substring Selection Problem (DSSSP) may be found on the repository jeanpttorres/dssp [2].

An alternative implementation of our algorithms is also presented in a more user oriented implementation using the AIMMS language on the repository jeanpttorres/aimms-dsssp [3].

The source code we present in here is written mainly in the C language and we use the SCIPT Optimization Suite [4]. We use the SCIP version 3.2.1. This version can be downloaded from here: https://www.scipopt.org/download.php?fname=scipoptsuite-3.2.1.tgz

Notice that in order to run SCIP a license is required, but SCIP offers an Academic License.

## Installation


## Usage

### Compiling

During compilation, the make command executed builds a object depending on the arguments passed trhough command line. They are:

```
TRACE         =       NOTRACE
# NOTRACE : shows no debug information
# DEBUG   : display debug prints trhough execution

EXECUTE       =       GLOBAL
# ONLYROOT  : execute a Linear Relaxation of the problem, solving only the root node from the branch-and-bound tree
# GLOBAL    : executes the entire branch-and-bound tree

PROBLEM       =       SUB
# SUB       : indicates if the target string is a substring 
# STRING    : indicates if the target string has length equal to the strings found in the sets Sc and Sf

HEUR       =       NONE
# RA      : runs the Roundg Algorithm as heuristic
# VNS     : runs the Variable Neighbourhood Search heuristic
# BCPA    : runs the Basic Core Problem Algorithm heuristic
# NONE    : runs no heuristic

TYPE       = 	   DSP
# DSP    : Considers the sets Sc and Sf
# CSP    : Considers only the set Sc
# FSP    : Considers only the set Sf

VNS_TYPE = NOTYPE
# VNS implementation to run
# BLPL 
# SBPL
# VBPL
# NOTYPE

TIMEH = SS
# Heuristic's execution time
# SSS = 5
# SS=30
# S=60
# M=120
# L=180
# LL=300

NSIZE = NSIZEM
# Number of neighbourhoods
# S=2
# M=3
# L=4
# LL=5

OMEGA = OMEGAM
#S=2, M=3, L=4, LL=5

IMPROVEMENT=BEST
#BEST, FIRST
```

Example

```
      make EXECUTE=ONLYROOT PROBLEM=SUB HEUR=RA TIMEH=M
```

## Execution

In order to execute the code, the object generated trhough compilation should be executed followed by a -f flag and argument, which indicates the instance file.

Example:

```
      ./dsssp.out -f instance.in
```

# References

[1] Torres, J.P.T., Hoshino, E.A. LP-based heuristics for the distinguishing string and substring selection problems. Ann Oper Res 316, 1205–1234 (2022). https://doi.org/10.1007/s10479-021-04138-5

[2] https://github.com/jeanpttorres/dssp

[3] https://github.com/jeanpttorres/aimms-dsssp

[4] Gamrath, G., Fischer, T., Gally, T., Gleixner, A.M., Hendel, G., Koch, T., Maher, S.J., Miltenberger, M., Müller, B., Pfetsch, M.E., Puchert, C., Rehfeldt, D., Schenker, S., Schwarz, R., Serrano, F., Shinano, Y., Vigerske, S., Weninger, D., Winkler, M., Witt, J.T., & Witzig, J. (2016). The scip optimization suite 3.2. Tech. Rep. 15-60, ZIB, Berlin
