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
#include "mumps_common.h"
/* Special case of mapping and pivnul_list -- allocated from MUMPS */
static MUMPS_INT * MUMPS_MAPPING;
static MUMPS_INT * MUMPS_PIVNUL_LIST;
/* as uns_perm and sym_perm */
static MUMPS_INT * MUMPS_SYM_PERM;
static MUMPS_INT * MUMPS_UNS_PERM;
MUMPS_INT*
mumps_get_mapping()
{
    return MUMPS_MAPPING;
}
void MUMPS_CALL
MUMPS_ASSIGN_MAPPING(MUMPS_INT * f77mapping)
{
    MUMPS_MAPPING = f77mapping;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_MAPPING()
{
    MUMPS_MAPPING = 0;
}
MUMPS_INT*
mumps_get_pivnul_list()
{
    return MUMPS_PIVNUL_LIST;
}
void MUMPS_CALL
MUMPS_ASSIGN_PIVNUL_LIST(MUMPS_INT * f77pivnul_list)
{
    MUMPS_PIVNUL_LIST = f77pivnul_list;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_PIVNUL_LIST()
{
    MUMPS_PIVNUL_LIST = 0;
}
MUMPS_INT*
mumps_get_sym_perm()
{
    return MUMPS_SYM_PERM;
}
void MUMPS_CALL
MUMPS_ASSIGN_SYM_PERM(MUMPS_INT * f77sym_perm)
{
    MUMPS_SYM_PERM = f77sym_perm;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_SYM_PERM()
{
    MUMPS_SYM_PERM = 0;
}
MUMPS_INT*
mumps_get_uns_perm()
{
    return MUMPS_UNS_PERM;
}
void MUMPS_CALL
MUMPS_ASSIGN_UNS_PERM(MUMPS_INT * f77uns_perm)
{
    MUMPS_UNS_PERM = f77uns_perm;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_UNS_PERM()
{
    MUMPS_UNS_PERM = 0;
}
