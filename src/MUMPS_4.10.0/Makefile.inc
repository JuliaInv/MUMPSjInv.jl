#
#  This file is part of MUMPS 4.8.0, built on Fri Jul 25 14:46:02 2008
#
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. If you would like to
# use other orderings, you need to obtain the corresponding package and modify
# the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis is now available as an internal ordering for MUMPS.
#



LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

LMETISDIR = $(topdir)/../metis-4.0.3/
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

LIBOTHERS = -lpthread
#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
#OPTF    = -O3 -xHost -vec-report0 -assume buffered_io -warn noalignments -Dintel_ -DALLOW_NON_INIT 
#OPTF    = -O3 -vec-report0 -assume buffered_io -warn noalignments -Dintel_ -DALLOW_NON_INIT 

OPTF    = -O2 -fPIC

OPTL    = -O2   $(SCALAP) 
OPTC    = -O1 -fPIC
#End Optimized options
INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded

