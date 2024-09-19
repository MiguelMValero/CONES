#ifndef __BFTC_CONFIG_DEFS_H__
#define __BFTC_CONFIG_DEFS_H__

/*============================================================================
 * Base macro and typedef definitions for system portability
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

#include "bftc_config.h"
#include "config_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
#ifndef BFTC_CPPCALLER
extern "C" {
#endif
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Basic macros
 *============================================================================*/


/*============================================================================
 * Suppress warning
 *============================================================================*/

#if defined(__INTEL_COMPILER)
#define BFTC_PRAGMA_TO_STR(x) _Pragma(#x)
#define BFTC_INTEL_SUPPRESS_WARNING_PUSH _Pragma("warning(push)")
#define BFTC_INTEL_SUPPRESS_WARNING(w) BFTC_PRAGMA_TO_STR(warning(disable:w))
#define BFTC_INTEL_SUPPRESS_WARNING_POP _Pragma("warning(pop)")
#define BFTC_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    BFTC_INTEL_SUPPRESS_WARNING_PUSH BFTC_INTEL_SUPPRESS_WARNING(w)
#else // BFTC_INTEL
#define BFTC_INTEL_SUPPRESS_WARNING_PUSH
#define BFTC_INTEL_SUPPRESS_WARNING(w)
#define BFTC_INTEL_SUPPRESS_WARNING_POP
#define BFTC_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // BFTC_INTEL

#if defined(__clang__)
#define BFTC_PRAGMA_TO_STR(x) _Pragma(#x)
#define BFTC_CLANG_SUPPRESS_WARNING_PUSH _Pragma("clang diagnostic push")
#define BFTC_CLANG_SUPPRESS_WARNING(w) BFTC_PRAGMA_TO_STR(clang diagnostic ignored w)
#define BFTC_CLANG_SUPPRESS_WARNING_POP _Pragma("clang diagnostic pop")
#define BFTC_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    BFTC_CLANG_SUPPRESS_WARNING_PUSH BFTC_CLANG_SUPPRESS_WARNING(w)
#else // BFTC_CLANG
#define BFTC_CLANG_SUPPRESS_WARNING_PUSH
#define BFTC_CLANG_SUPPRESS_WARNING(w)
#define BFTC_CLANG_SUPPRESS_WARNING_POP
#define BFTC_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // BFTC_CLANG

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#define BFTC_PRAGMA_TO_STR(x) _Pragma(#x)
#define BFTC_GCC_SUPPRESS_WARNING_PUSH _Pragma("GCC diagnostic push")
#define BFTC_GCC_SUPPRESS_WARNING(w) BFTC_PRAGMA_TO_STR(GCC diagnostic ignored w)
#define BFTC_GCC_SUPPRESS_WARNING_POP _Pragma("GCC diagnostic pop")
#define BFTC_GCC_SUPPRESS_WARNING_WITH_PUSH(w)                                                  \
    BFTC_GCC_SUPPRESS_WARNING_PUSH BFTC_GCC_SUPPRESS_WARNING(w)
#else // BFTC_GCC
#define BFTC_GCC_SUPPRESS_WARNING_PUSH
#define BFTC_GCC_SUPPRESS_WARNING(w)
#define BFTC_GCC_SUPPRESS_WARNING_POP
#define BFTC_GCC_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // DOCTEST_GCC

#define BFTC_UNUSED(x) (void)(x)

BFTC_GCC_SUPPRESS_WARNING("-Wcast-qual")
BFTC_GCC_SUPPRESS_WARNING("-Wunknown-pragmas")


/*============================================================================
 * Internationalization (future)
 *============================================================================*/

#if defined(ENABLE_NLS)

#  include <libintl.h>
#  define _(String) dgettext(PACKAGE,String)
#  ifdef gettext_noop
#    define N_(String) gettext_noop(String)
#  else
#    define N_(String) String
#  endif /* gettext_noop */

#else

#  define _(String) (String)
#  define N_(String) String
#  define textdomain(String) (String)
#  define gettext(String) (String)
#  define dgettext(Domain,String) (String)
#  define dcgettext(Domain,String,Type) (String)
#  define bindtextdomain(Domain,Directory) (Domain)

#endif

/*============================================================================
 * C99 Qualifiers
 *============================================================================*/

/* inline provided by bftc_config.h */

/* restrict type qualifier (standard in C99) */

#if defined(__GNUC__)
#define restrict __restrict
#else
#define restrict
#endif

/*============================================================================
 * Definitions that may not always be provided directly by the system
 *============================================================================*/

/*
 * Obtain definitions such as that of size_t through stddef.h (C99 standard)
 * if available (preferred method), or through stdlib.h (which defines
 * malloc() and family and so must define size_t some way) otherwise.
 */

#if HAVE_STDDEF_H
# include <stddef.h>
#else
# include <stdlib.h>
#endif

/*
 * Usually stdint.h is included by inttypes.h, but only inttypes.h exists
 * on certain systems, such as Tru64Unix
 */

#if HAVE_STDINT_H
# include <stdint.h>
#elif HAVE_INTTYPES_H
# include <inttypes.h>
#endif

/* C99 _Bool type */

#ifdef __cplusplus
typedef bool _Bool;
#else
# if HAVE_STDBOOL_H
#  include <stdbool.h>
# else
#  if !HAVE__BOOL
typedef unsigned char _Bool;
#  endif
# endif
# define bool _Bool
# define false 0
# define true 1
#endif

/* int32_t type */

#if !defined(HAVE_INT32_T)
# if (BFTC_SIZEOF_INT == 4)
typedef int int32_t;
# elif (BFTC_SIZEOF_SHORT == 4)
typedef short int32_t;
# else
#  error
# endif
#endif

/* Directory name separator ('/' for Unix/Linux, '\' for Windows, ':' for Mac */

#define BFTC_DIR_SEPARATOR '/'

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
#ifndef BFTC_CPPCALLER
} 
#endif
#endif /* __cplusplus */

#endif /* __BFTC_CONFIG_DEFS_H__ */
