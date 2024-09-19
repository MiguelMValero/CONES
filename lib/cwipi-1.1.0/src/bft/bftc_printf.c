/*============================================================================
 * Base user-definable printf() wrapper or replacement.
 *============================================================================*/

/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2004-2006  EDF

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

#include "bftc_config_defs.h"

/*
 * Standard C library headers
 */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Optional library and BFT headers
 */

#include "bftc_printf.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/* Associated typedef documentation (for bftc_printf.h) */

/*!
 * \typedef bftc_printf_proxy_t
 *
 * \brief Function pointer for printf() type functions.
 *
 * \param [in] format       format string, as printf() and family.
 * \param [in, out] arg_ptr pointer to variable argument list based on format
 *                          string.
 */

/*!
 * \typedef bftc_printf_flush_proxy_t
 *
 * \brief Function pointer for fflush(stdout) type functions.
 */

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default bftc_printf_flush() proxy.
 *
 * returns:
 *   return code of fflush(stdout).
 */

static int
_bftc_printf_flush_proxy_default(void);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static bftc_printf_proxy_t        *_bftc_printf_proxy = vprintf;
static bftc_printf_flush_proxy_t  *_bftc_printf_flush_proxy
                                    = _bftc_printf_flush_proxy_default;

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Default bftc_printf_flush() proxy.
 *
 * returns:
 *   return code of fflush(stdout).
 */

static int
_bftc_printf_flush_proxy_default(void)
{
  return fflush(stdout);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Replacement for printf() with modifiable behavior.
 *
 * This function calls vprintf() by default, or a function with similar
 * arguments indicated by bftc_printf_proxy_set().
 *
 * \param [in] format format string, as printf() and family.
 * \param [in] ...    variable arguments based on format string.
 *
 * \return number of characters printed, not counting the trailing '\\0'
 *         used to end output strings
 */

int
bftc_printf(const char *const format,
           ...)
{
  int  retval;
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  retval = _bftc_printf_proxy(format, arg_ptr);

  va_end(arg_ptr);

  return retval;
}

/*!
 * \brief Flush for output of bftc_printf() with modifiable behavior.
 *
 * This function calls fflush(stdout) if bftc_printf()'s default behavior is
 * used. If bftc_printf's behavior is modified with bftc_printf_proxy_set(),
 * bftc_printf_flush()'s behavior may have to be also adjusted with
 * bftc_printf_flush_proxy_set().
 *
 * \return using the default behavior, the return value is that of
 *         fflush(stdout): O upon successful completion, EOF otherwise
 *         (with errno set to indicate the error).
 */

int
bftc_printf_flush(void)
{
  return _bftc_printf_flush_proxy();
}

/*!
 * \brief Returns function associated with the bftc_printf() function.
 *
 * \return pointer to the vprintf() or replacement function.
 */

bftc_printf_proxy_t *
bftc_printf_proxy_get(void)
{
  return _bftc_printf_proxy;
}

/*!
 * \brief Associates a vprintf() type function with the bftc_printf() function.
 *
 * \param [in] fct pointer to a vprintf() type function.
 */

void
bftc_printf_proxy_set(bftc_printf_proxy_t *const fct)
{
  _bftc_printf_proxy = fct;
}

/*!
 * \brief Returns function associated with bftc_printf_flush().
 *
 * \return pointer to the bftc_printf_flush() proxy.
 */

bftc_printf_flush_proxy_t *
bftc_printf_flush_proxy_get(void)
{
  return _bftc_printf_flush_proxy;
}

/*!
 * \brief Associates a proxy function with bftc_printf_flush().
 *
 * \warning
 *   bftc_printf() is called by the default bftc_error() error handler
 *   (so as to ensure that the error text appears at the end of the
 *   program output), so a bftc_print_flush replacement must not itself
 *   call (directly or indirectly) bftc_error() if the default error
 *   handler is used.
 *
 * \param [in] fct pointer to a function similar to {return fflush(stdout)}.
 */

void
bftc_printf_flush_proxy_set(bftc_printf_flush_proxy_t *const fct)
{
  _bftc_printf_flush_proxy = fct;
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
