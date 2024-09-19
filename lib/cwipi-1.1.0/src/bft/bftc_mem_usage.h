#ifndef __BFTC_MEM_USAGE_H__
#define __BFTC_MEM_USAGE_H__

/*============================================================================
 * Base memory usage information (System and Library dependent)
 *============================================================================*/

/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2004-2009  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

/*
 * Obtain definitions such as that of size_t through stddef.h (C99 standard)
 * if available (preferred method), or through stdlib.h (which defines
 * malloc() and family and so must define size_t some way) otherwise.
 */

#if defined(__STDC_VERSION__)
#  if (__STDC_VERSION__ == 199901L)
#    include <stddef.h>
#  else
#    include <stdlib.h>
#  endif
#else
#  include <stdlib.h>
#endif

/* BFT headers */

#include "bftc_config.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Initialize memory usage count depending on system.
 *
 * This functions checks if it has already been called, so
 * it is safe to call more than once (though it is not
 * thread-safe). Only the first call is effective.
 */

void
bftc_mem_usage_init(void);

/*
 * End memory usage count depending on system.
 */

void
bftc_mem_usage_end(void);

/*
 * Indicates if bftc_mem_usage_...() functions are initialized.
 *
 * returns:
 *   1 if bftc_mem_usage_init has been called, 0 otherwise.
 */

int
bftc_mem_usage_initialized(void);

/*
 * Return current process memory use (in kB) depending on OS.
 */

size_t
bftc_mem_usage_pr_size(void);

/*
 * Return maximum process memory use (in kB) depending on OS.
 *
 * The returned value is the maximum returned by bftc_mem_usage_pr_size()
 * during the program's lifetime. With memory allocations which return
 * memory to the system (such as the GNU glibc on Linux systems),
 * this value will be correct only if allocation is tracked. This should
 * be the case if malloc hooks are used with the glibc allocation
 * functions (BFT library's default configuration/installation option),
 * but may give results lower than the true maximum in other cases.
 */

size_t
bftc_mem_usage_max_pr_size(void);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BFTC_MEM_USAGE_H__ */
