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
#ifdef INTSIZE64
#include <inttypes.h>
#define LIBSEQ_INT int64_t
#else
#define LIBSEQ_INT int
#endif

#ifndef MUMPS_MPI_H
#define MUMPS_MPI_H

/* We define all symbols as extern "C" for users who call MUMPS with its
   libseq from a C++ driver. */
#ifdef __cplusplus
extern "C" {
#endif

/* This is the minimum to have the C interface of MUMPS work.
 * Most of the time, users who need this file have no call to MPI functions in
 * their own code. Hence it is not worth declaring all MPI functions here.
 * However if some users come to request some more stub functions of the MPI
 * standards, we may add them. But it is not worth doing it until then. */

typedef LIBSEQ_INT MPI_Comm; /* Simple type for MPI communicator */
static MPI_Comm MPI_COMM_WORLD=(MPI_Comm)0;

LIBSEQ_INT MPI_Init(LIBSEQ_INT *pargc, char ***pargv);
LIBSEQ_INT MPI_Comm_rank(LIBSEQ_INT  comm, LIBSEQ_INT  *rank);
LIBSEQ_INT MPI_Finalize(void);

#ifdef __cplusplus
}
#endif

#endif /* MUMPS_MPI_H */
