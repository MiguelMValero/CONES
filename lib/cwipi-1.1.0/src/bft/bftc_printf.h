#ifndef __BFTC_PRINTF_H__
#define __BFTC_PRINTF_H__

/*============================================================================
 * Base user-definable printf() wrapper or replacement
 *============================================================================*/

/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2004  EDF

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

/* Function pointers for printf() and fflush(stdout) type functions */

typedef int (bftc_printf_proxy_t) (const char     *const format,
                                  va_list               arg_ptr);

typedef int (bftc_printf_flush_proxy_t) (void);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Replacement for printf() with modifiable behavior.
 *
 * This function calls vprintf() by default, or a function with similar
 * arguments indicated by bftc_printf_proxy_set().
 *
 * parameters:
 *   format: <-- format string, as printf() and family.
 *   ... :   <-- variable arguments based on format string.
 *
 * returns:
 *   number of characters printed, not counting the trailing '\0' used
 *   to end output strings
 */

int
bftc_printf(const char  *const format,
           ...);

/*
 * Flush for output of bftc_printf() with modifiable behavior.
 *
 * This function calls fflush(stdout) if bftc_printf()'s default behavior is
 * used. If bftc_printf's behavior is modified with bftc_printf_proxy_set(),
 * bftc_printf_flush()'s behavior may have to be also adjusted with
 * bftc_printf_flush_proxy_set().
 *
 * returns:
 *   using the default behavior, the return value is that of
 *   fflush(stdout): O upon successful completion, EOF otherwise
 *   (with errno set to indicate the error).
 */

int
bftc_printf_flush(void);

/*
 * Returns function associated with the bftc_printf() function.
 *
 * returns:
 *   pointer to the vprintf() or replacement function.
 */

bftc_printf_proxy_t *
bftc_printf_proxy_get(void);

/*
 * Associates a vprintf() type function with the bftc_printf() function.
 *
 * parameters:
 *   fct: <-- pointer to a vprintf() type function.
 */

void
bftc_printf_proxy_set(bftc_printf_proxy_t  *const fct);

/*
 * Returns function associated with bftc_printf_flush().
 *
 * returns:
 *   pointer to the bftc_printf_flush() proxy.
 */

bftc_printf_flush_proxy_t *
bftc_printf_flush_proxy_get(void);

/*
 * Associates a proxy function with bftc_printf_flush().
 *
 * warning:
 *   bftc_printf() is called by the default bftc_error() error handler
 *   (so as to ensure that the error text appears at the end of the
 *   program output), so a bftc_print_flush replacement must not itself
 *   call (directly or indirectly) bftc_error() if the default error
 *   handler is used.
 *
 * parameter:
 *   fct <-- pointer to a function similar to {return fflush(stdout)}.
 */

void
bftc_printf_flush_proxy_set(bftc_printf_flush_proxy_t  *const fct);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BFTC_PRINTF_H__ */
