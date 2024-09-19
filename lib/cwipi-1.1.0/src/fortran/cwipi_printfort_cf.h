#ifndef __CWIPI_PRINTFORT_CF_H__
#define __CWIPI_PRINTFORT_CF_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2021-2023  ONERA

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

#include <mpi.h>
#include <stdio.h>

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

/*----------------------------------------------------------------------------
 * Macro used to handle automatic "Fortran string length" arguments
 * (not used by CWIPI, but set by many compilers).
 * Some compilers, like the Fujitsu VPP 5000 compiler, may not
 * support the variable length lists in mixed C/Fortran calls.
 *----------------------------------------------------------------------------*/

#if defined (__uxpv__)  /* Fujitsu VPP 5000 case */
#define ARGF_SUPP_CHAINE
#else
#define ARGF_SUPP_CHAINE , ...
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Set up the file used for the output listing
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_output_listing_cf, CWIPI_SET_OUTPUT_LISTING_CF)
     ();


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWIPI_H__ */
