#
#  This file is part of MUMPS 4.8.0, built on Fri Jul 25 14:46:02 2008
#
#Begin orderings

LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

LMETISDIR = $(topdir)/../metis-5.1.0/installed-libs/lib
IMETIS    = # Metis doesn't need include files (Fortran interface avail.)
LMETIS    = -L$(LMETISDIR) -lmetis


# The following variables will be used in the compilation process.
#ORDERINGSF = -Dscotch -Dmetis -Dpord -Dparmetis
ORDERINGSF = -Dmetis -Dpord
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD)
#IORDERINGS = $(IMETIS) $(IPORD) $(ISCOTCH)
#IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD)


#End orderings
########################################################################
################################################################################

PLAT    =
LIBEXT  = .a
OUTC    = -o
OUTF    = -o
RM = /bin/rm -f
CC = gcc
FC = gfortran
FL = gfortran
AR = ar vr 
#RANLIB = ranlib
RANLIB  = echo
INCSEQ = -I$(topdir)/libseq
LIBSEQ  =  -L$(topdir)/libseq -lmpiseq

#LIBPAR = $(SCALAP)

OS := $(shell uname)

ifeq ($(OS),Darwin)
	LIBBLAS = /System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libBLAS.dylib
else
	 BLDR = /opt/OpenBLAS/lib
	 LIBBLAS = $(BLDR)/libopenblas_haswellp-r0.2.14.a
endif

LIBOTHERS = -lpthread
#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
#OPTF    = -O3 -xHost -vec-report0 -assume buffered_io -warn noalignments -Dintel_ -DALLOW_NON_INIT 
#OPTF    = -O3 -vec-report0 -assume buffered_io -warn noalignments -Dintel_ -DALLOW_NON_INIT 

OPTF    = -O3 -fPIC -fdefault-integer-8

OPTL    = -O3   $(SCALAP) 
OPTC    = -O3 -fPIC -DINTSIZE64
#End Optimized options
INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded

