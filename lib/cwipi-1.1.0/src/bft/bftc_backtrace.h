#ifndef __BFTC_BACKTRACE_H__
#define __BFTC_BACKTRACE_H__

/*============================================================================
 * Obtaining a stack backtrace
 *============================================================================*/

/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2006  EDF

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

/*-----------------------------------------------------------------------------*/

/* Standard C library headers */

#include <stdarg.h>

/* BFT library headers */

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

/* BFT backtrace descriptor */

typedef struct _bftc_backtrace_t bftc_backtrace_t;

/* Pointers for backtrace print functions */

typedef void (bftc_backtrace_print_t) (int  start_depth);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Build a backtrace description structure.
 *
 * returns:
 *   pointer to bftc_backtrace_t backtrace descriptor (NULL in case of
 *   error, or if backtracing is unavailable on this architecture).
 */

bftc_backtrace_t *
bftc_backtrace_create(void);

/*
 * Free a backtrace description structure.
 *
 * parameters:
 *   bt: <-> pointer to backtrace description structure.
 *
 * returns:
 *   NULL pointer.
 */

bftc_backtrace_t *
bftc_backtrace_destroy(bftc_backtrace_t  *bt);

/*!
 * Demangle a backtrace description structure (for C++).
 *
 * parameters:
 *   bt: <-> pointer to backtrace description structure.
 */

void
bftc_backtrace_demangle(bftc_backtrace_t  *bt);

/*
 * Return the total depth of a backtrace.
 *
 * parameters:
 *   bt: <-- pointer to backtrace description structure.
 *
 * returns:
 *   total backtrace depth.
 */

int
bftc_backtrace_size(const bftc_backtrace_t  *bt);

/*
 * Return file name associated with a backtrace at a given depth.
 *
 * parameters:
 *   bt:    <-- pointer to backtrace description structure.
 *   depth: <-- index in backtrace structure (< bftc_backtrace_size(bt)).
 *
 * returns:
 *   file name at the given depth, or NULL.
 */

const char *
bftc_backtrace_file(bftc_backtrace_t  *bt,
                   int               depth);

/*
 * Return function name associated with a backtrace at a given depth.
 *
 * parameters:
 *   bt:    <-- pointer to backtrace description structure.
 *   depth: <-- index in backtrace structure (< bftc_backtrace_size(bt)).
 *
 * returns:
 *   function name at the given depth, or NULL.
 */

const char *
bftc_backtrace_function(bftc_backtrace_t  *bt,
                       int               depth);

/*
 * Return address associated with a backtrace at a given depth.
 *
 * parameters:
 *   bt:    <-- pointer to backtrace description structure.
 *   depth: <-- index in backtrace structure (< bftc_backtrace_size(bt)).
 *
 * returns:
 *   address at the given depth, or NULL.
 */

const char *
bftc_backtrace_address(bftc_backtrace_t  *bt,
                      int               depth);

/*
 * Print a backtrace.
 *
 * parameters:
 *   start_depth: <-- depth of backtrace at which to start printing
 *                    (0 for all, including backtrace print function)
 */

void
bftc_backtrace_print(int  start_depth);

/*
 * Returns backtrace print function.
 *
 * returns:
 *   pointer to the backtrace print function.
 */

bftc_backtrace_print_t *
bftc_backtrace_print_get(void);

/*
 * Sets a backtrace print function.
 *
 * parameters:
 *   fct: <-- pointer to a bftc_backtrace_print_t type function.
 */

void
bftc_backtrace_print_set(bftc_backtrace_print_t  *const fct);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BFTC_BACKTRACE_H__ */
