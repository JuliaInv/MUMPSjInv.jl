
topdir = MUMPS_4.10.0/
libdir = $(topdir)/lib

include $(topdir)/Makefile.inc
OS := $(shell uname)

ifeq ($(OS),Darwin)
	LIBBLAS = /System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libBLAS.dylib
else
	 LIBBLAS = /usr/lib/atlas-base/atlas/libblas.so.3 
endif

SRC = mumps_cmplx_p.f90  mumps_p.f90 call_mumps_cmplx_p.f90 call_mumps_p.f90  
LIB = -L$(topdir)/lib -ldmumps -lzmumps -lmumps_common -lpord  -L$(LMETISDIR) -lmetis  $(LIBSEQ)  $(LIBBLAS)


all: $(SRC)
	  (cd metis-4.0.3; make)
	 (cd MUMPS_4.10.0; make; make z)
	 $(FL) -shared -fPIC -fcray-pointer $(SRC) -o ../lib/MUMPS -I$(topdir)/include $(LIB)
clean:
	(cd metis-4.0.3; make realclean)
	(cd MUMPS_4.10.0; make clean)
	rm ../lib/MUMPS