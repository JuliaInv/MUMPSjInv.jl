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
#include "mumps_io_err.h"
#include "mumps_io_basic.h"
#include "mumps_c_types.h"
#if defined( MUMPS_WIN32 )
# include <string.h>
#endif
/* Exported global variables */
char* mumps_err;
MUMPS_INT* dim_mumps_err;
MUMPS_INT mumps_err_max_len;
MUMPS_INT err_flag;
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
pthread_mutex_t err_mutex;
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
/* Functions */
/* Keeps a C pointer to store error description string that will be
   displayed by the Fortran layers.
   * dim contains the size of the Fortran character array to store the
   description.
*/
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_ERR_STR(MUMPS_INT *dim, char* err_str, mumps_ftnlen l1){
  mumps_err = err_str;
  dim_mumps_err = (MUMPS_INT *) dim;
  mumps_err_max_len = (MUMPS_INT) *dim;
  err_flag = 0;
  return;
}
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
MUMPS_INLINE MUMPS_INT
mumps_io_protect_err()
{
  if(mumps_io_flag_async==IO_ASYNC_TH){
    pthread_mutex_lock(&err_mutex);
  }
  return 0;
}
MUMPS_INLINE MUMPS_INT
mumps_io_unprotect_err()
{
  if(mumps_io_flag_async==IO_ASYNC_TH){
    pthread_mutex_unlock(&err_mutex);
  }
  return 0;
}
MUMPS_INT
mumps_io_init_err_lock()
{
  pthread_mutex_init(&err_mutex,NULL);
  return 0;
}
MUMPS_INT
mumps_io_destroy_err_lock()
{
  pthread_mutex_destroy(&err_mutex);
  return 0;
}
MUMPS_INT
mumps_check_error_th()
{
  /* If err_flag != 0, then error_str is set */
  return err_flag;
}
#endif /* MUMPS_WIN32 && WITHOUT_PTHREAD */
MUMPS_INT
mumps_io_error(MUMPS_INT mumps_errno, const char* desc)
{
    MUMPS_INT len;
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_protect_err();
#endif
  if(err_flag == 0){
    strncpy(mumps_err, desc, mumps_err_max_len);
    /* mumps_err is a FORTRAN string, we do not care about adding a final 0 */
    len = (MUMPS_INT) strlen(desc);
    *dim_mumps_err = (len <= mumps_err_max_len ) ? len : mumps_err_max_len;
    err_flag = mumps_errno;
  }
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_unprotect_err();
#endif
  return mumps_errno;
}
MUMPS_INT
mumps_io_sys_error(MUMPS_INT mumps_errno, const char* desc)
{
  MUMPS_INT len = 2; /* length of ": " */
  const char* _desc;
  char* _err;
#if defined( MUMPS_WIN32 )
  MUMPS_INT _err_len;
#endif
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_protect_err();
#endif
  if(err_flag==0){
    if(desc == NULL) {
      _desc = "";
    } else {
        len += (MUMPS_INT) strlen(desc);
      _desc = desc;
    }
#if ! defined( MUMPS_WIN32 )
    _err = strerror(errno);
    len += (MUMPS_INT) strlen(_err);
    snprintf(mumps_err, mumps_err_max_len, "%s: %s", _desc, _err);
    /* mumps_err is a FORTRAN string, we do not care about adding a final 0 */
#else
    /* This a VERY UGLY workaround for snprintf: this function has been
     * integrated quite lately into the ANSI stdio: some windows compilers are
     * not up-to-date yet. */
    if( len >= mumps_err_max_len - 1 ) { /* then do not print sys error msg at all */
      len -= 2;
      len = (len >= mumps_err_max_len ) ? mumps_err_max_len - 1 : len;
      _err = strdup( _desc );
      _err[len] = '\0';
      sprintf(mumps_err, "%s", _err);
    } else {
      _err = strdup(strerror(errno));
      _err_len = (MUMPS_INT) strlen(_err);
      /* We will use sprintf, so make space for the final '\0' ! */
      if((len + _err_len) >= mumps_err_max_len) {
        /* truncate _err, not to overtake mumps_err_max_len at the end. */
        _err[mumps_err_max_len - len - 1] = '\0';
        len = mumps_err_max_len - 1;
      } else {
        len += _err_len;
      }
      sprintf(mumps_err, "%s: %s", _desc, _err);
    }
    free(_err);
#endif
    *dim_mumps_err = (len <= mumps_err_max_len ) ? len : mumps_err_max_len;
    err_flag = mumps_errno;
  }
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_unprotect_err();
#endif
  return mumps_errno;
}
