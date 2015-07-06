/*
 *
 *  This file is part of MUMPS 5.0.0, released
 *  on Fri Feb 20 08:19:56 UTC 2015
 *
 *
 *  Copyright 1991-2015 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license,
 *  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html, 
 *  except for the external and optional ordering PORD, 
 *  in separate directory PORD, which is public domain (see PORD/README).
 *
 *  You can acknowledge (using references [1] and [2]) the contribution of
 *  this package in any scientific publication dependent upon the use of
 *  the package. Please use reasonable endeavours to notify the authors
 *  of the package of this publication.
 *
 *   [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
 *   A fully asynchronous multifrontal solver using distributed dynamic
 *   scheduling, SIAM Journal of Matrix Analysis and Applications,
 *   Vol 23, No 1, pp 15-41 (2001).
 *
 *   [2] P. R. Amestoy, A. Guermouche, J.-Y. L'Excellent and
 *   S. Pralet, Hybrid scheduling for the parallel solution of linear
 *   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
 *
 *  As a counterpart to the access to the source code and rights to copy,
 *  modify and redistribute granted by the license, users are provided only
 *  with a limited warranty  and the software's author,  the holder of the
 *  economic rights,  and the successive licensors  have only  limited
 *  liability. 
 *
 *  In this respect, the user's attention is drawn to the risks associated
 *  with loading,  using,  modifying and/or developing or reproducing the
 *  software by the user in light of its specific status of free software,
 *  that may mean  that it is complicated to manipulate,  and  that  also
 *  therefore means  that it is reserved for developers  and  experienced
 *  professionals having in-depth computer knowledge. Users are therefore
 *  encouraged to load and test the software's suitability as regards their
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and,  more generally, to use and operate it in the 
 *  same conditions as regards security. 
 *
 *  The fact that you are presently reading this means that you have had
 *  knowledge of the CeCILL-C license and that you accept its terms.
 *
 */
#ifndef MUMPS_ORDERINGS_H
#define MUMPS_ORDERINGS_H
#include "mumps_common.h"
#include "mumps_c_types.h"
#if defined(pord)
#include <space.h>
MUMPS_INT mumps_pord( MUMPS_INT, MUMPS_INT, MUMPS_INT *, MUMPS_INT *, MUMPS_INT * );
#define MUMPS_PORDF \
    F_SYMBOL(pordf,PORDF)
void MUMPS_CALL
MUMPS_PORDF( MUMPS_INT *nvtx, MUMPS_INT *nedges,
             MUMPS_INT *xadj, MUMPS_INT *adjncy,
             MUMPS_INT *nv, MUMPS_INT *ncmpa );
MUMPS_INT mumps_pord_wnd( MUMPS_INT, MUMPS_INT, MUMPS_INT *, MUMPS_INT *, MUMPS_INT *, MUMPS_INT * );
#define MUMPS_PORDF_WND          \
    F_SYMBOL(pordf_wnd,PORDF_WND)
void MUMPS_CALL
MUMPS_PORDF_WND( MUMPS_INT *nvtx, MUMPS_INT *nedges,
                 MUMPS_INT *xadj, MUMPS_INT *adjncy,
                 MUMPS_INT *nv, MUMPS_INT *ncmpa, MUMPS_INT *totw );
#endif /*PORD*/
#if defined(scotch) || defined(ptscotch)
MUMPS_INT esmumps( const MUMPS_INT n, const MUMPS_INT iwlen, MUMPS_INT * const pe, const MUMPS_INT pfree,
             MUMPS_INT * const len, MUMPS_INT * const iw, MUMPS_INT * const nv, MUMPS_INT * const elen,
             MUMPS_INT * const last);
#define MUMPS_SCOTCH        \
    F_SYMBOL(scotch,SCOTCH)
void MUMPS_CALL
MUMPS_SCOTCH( const MUMPS_INT * const  n,
              const MUMPS_INT * const  iwlen,
              MUMPS_INT * const        petab,
              const MUMPS_INT * const  pfree,
              MUMPS_INT * const        lentab,
              MUMPS_INT * const        iwtab,
              MUMPS_INT * const        nvtab,
              MUMPS_INT * const        elentab,
              MUMPS_INT * const        lasttab,
              MUMPS_INT * const        ncmpa );
#endif /*scotch or ptscotch*/
#if defined(ptscotch)
#include "mpi.h"
#include <stdio.h>
#include "ptscotch.h"
#define MUMPS_DGRAPHINIT \
  F_SYMBOL(dgraphinit,DGRAPHINIT)
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr);
#endif /*ptscotch*/
#if defined(parmetis) || defined(parmetis3)
#include "mpi.h"
#include "parmetis.h"
#define MUMPS_PARMETIS \
  F_SYMBOL(parmetis,PARMETIS)
void MUMPS_CALL
MUMPS_PARMETIS(MUMPS_INT *first,      MUMPS_INT *vertloctab, 
               MUMPS_INT *edgeloctab, MUMPS_INT *numflag, 
               MUMPS_INT *options,    MUMPS_INT *order, 
               MUMPS_INT *sizes,      MUMPS_INT *comm,
               MUMPS_INT *ierr);
#endif /*PARMETIS*/
#endif /* MUMPS_ORDERINGS_H */
