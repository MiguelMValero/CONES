
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/
/*----------------------------------------------------------------------------
 * Standard C/C++ library headers
 *----------------------------------------------------------------------------*/

#include <mpi.h>

#include <cassert>
#include <cstring>
#include <cstdlib>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

// #include "bftc_printf.h"
#include "pdm_printf.h"

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fortran/new/cwp_printfort_cf.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Output listing File (C printing)
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *
 * Set bftc_printf proxy for Fortran interface
 *
 *----------------------------------------------------------------------------*/

#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
static int
_cwp_print_with_fortran
(
 const char     *const format,
       va_list         arg_ptr
)
{
  int  msgsize;

 /* Tampon pour impressions depuis du code C : on imprime dans un chaîne
    de caractères, qui sera imprimée vers un fichier par du code Fortran.
    Une fois les impressions Fortran totalement remplacées par des impressions
    C, on pourra supprimer cette étape, mais elle est nécessaire pour l'instant
    afin de pouvoir utiliser les mêmes fichiers de sortie */

#undef BUF_PRINT_F_SIZE
#define BUF_PRINT_F_SIZE 16384

  static char buf_print_f[BUF_PRINT_F_SIZE];

 /* Impression dans le tampon */

#if defined  __STDC_VERSION__
  msgsize = vsprintf (buf_print_f, format, arg_ptr);
#else
  msgsize = vsnprintf (buf_print_f, BUF_PRINT_F_SIZE, format, arg_ptr);
#endif

  if (msgsize == -1 || msgsize > BUF_PRINT_F_SIZE - 1) {
    fprintf(stderr,
            "Fatal error: printf() called on a message of size %d\n"
            "whereas the print buffer is of size %d.",
            msgsize, BUF_PRINT_F_SIZE);

    /* Try to force segmentation fault (to call signal handlers);
       as stack has most likely been corrupted, this is the most
       "similar" error that allows for portable handling. */
    {
      int *_force_err = NULL;
      *_force_err = 0;
    }
    exit(1);
  }

  /* Impression effective par le code Fortran */

  printfortran (buf_print_f, &msgsize);
  return msgsize;
}
#endif

#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
static int
_cwp_flush_with_fortran(void)
{
  flushfortran();
  return 42; // TO DO: how to know fortran flush failed?
}
#endif

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Set up the file used for the output listing
 *
 *----------------------------------------------------------------------------*/

#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
void cwp_set_output_listing_cf ()
{
  PDM_printf_proxy_set(_cwp_print_with_fortran);
  PDM_printf_flush_proxy_set(_cwp_flush_with_fortran);

}
#endif


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
