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
#ifndef MUMPS_COMMON_H
#define MUMPS_COMMON_H
#include "mumps_compat.h"
#include "mumps_c_types.h"
/**
 * F_SYMBOL is a macro that converts a couple (lower case symbol, upper
 * case symbol) into the symbol defined by the compiler convention.
 * Example: For MUMPS_XXX, first define
 *   #define MUMPS_XXX F_SYMBOL(xxx,XXX) and then use
 *   MUMPS_XXX in the code to get rid of any symbol convention annoyance.
 *
 * NB: We need to provide both upper and lower case versions because to our
 *     knowledge, there is no way to perform the conversion with CPP
 *     directives only.
 */
#if defined(UPPER) || defined(MUMPS_WIN32)
# define F_SYMBOL(lower_case,upper_case) MUMPS_##upper_case
#elif defined(Add_)
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case##_
#elif defined(Add__)
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case##__
#else
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case
#endif
MUMPS_INT*
mumps_get_mapping();
#define MUMPS_ASSIGN_MAPPING \
    F_SYMBOL(assign_mapping,ASSIGN_MAPPING)
void MUMPS_CALL
MUMPS_ASSIGN_MAPPING(MUMPS_INT *f77mapping);
#define MUMPS_NULLIFY_C_MAPPING F_SYMBOL(nullify_c_mapping,NULLIFY_C_MAPPING)
void MUMPS_CALL
MUMPS_NULLIFY_C_MAPPING();
MUMPS_INT*
mumps_get_pivnul_list();
#define MUMPS_ASSIGN_PIVNUL_LIST \
    F_SYMBOL(assign_pivnul_list,ASSIGN_PIVNUL_LIST)
void MUMPS_CALL
MUMPS_ASSIGN_PIVNUL_LIST(MUMPS_INT *f77pivnul_list);
#define MUMPS_NULLIFY_C_PIVNUL_LIST \
    F_SYMBOL(nullify_c_pivnul_list,NULLIFY_C_PIVNUL_LIST)
void MUMPS_CALL
MUMPS_NULLIFY_C_PIVNUL_LIST();
MUMPS_INT*
mumps_get_uns_perm();
#define MUMPS_ASSIGN_UNS_PERM \
    F_SYMBOL(assign_uns_perm,ASSIGN_UNS_PERM)
void MUMPS_CALL
MUMPS_ASSIGN_UNS_PERM(MUMPS_INT *f77sym_perm);
#define MUMPS_NULLIFY_C_UNS_PERM \
    F_SYMBOL(nullify_c_uns_perm,NULLIFY_C_UNS_PERM)
void MUMPS_CALL
MUMPS_NULLIFY_C_UNS_PERM();
MUMPS_INT*
mumps_get_sym_perm();
#define MUMPS_ASSIGN_SYM_PERM \
    F_SYMBOL(assign_sym_perm,ASSIGN_SYM_PERM)
void MUMPS_CALL
MUMPS_ASSIGN_SYM_PERM(MUMPS_INT * f77sym_perm);
#define MUMPS_NULLIFY_C_SYM_PERM \
    F_SYMBOL(nullify_c_sym_perm,NULLIFY_C_SYM_PERM)
void MUMPS_CALL
MUMPS_NULLIFY_C_SYM_PERM();
#endif /* MUMPS_COMMON_H */
