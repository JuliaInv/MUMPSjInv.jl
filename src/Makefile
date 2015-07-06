srcdir = /home/patrick/Julia/v0.4/MUMPS.jl/src
topdir = $(srcdir)/MUMPS_5.0.0


include $(topdir)/Makefile.inc
OS := $(shell uname)


SRC = mumps_cmplx_p.f90  mumps_p.f90 call_mumps_cmplx_p.f90 call_mumps_p.f90  
LIB = -L$(topdir)/lib -ldmumps -lzmumps -lmumps_common -lpord  -L$(LMETISDIR) -lmetis  $(LIBSEQ)  $(LIBBLAS)

all: $(SRC)
	(cd metis-5.1.0; make config prefix=$(srcdir)/metis-5.1.0/installed-libs ; make; make install)
	(cd MUMPS_5.0.0; make d z)
	$(FL) -shared -fPIC -fcray-pointer -O2 $(SRC) -o ../lib/MUMPS -I$(topdir)/include $(LIB)
clean:
	(cd MUMPS_5.0.0; make clean)
	rm ../lib/MUMPS
	(cd metis-5.1.0; make clean; make uninstall; make distclean)