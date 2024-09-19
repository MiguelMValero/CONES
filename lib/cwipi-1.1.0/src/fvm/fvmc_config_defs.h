#ifndef __FVMC_CONFIG_DEFS_H__
#define __FVMC_CONFIG_DEFS_H__

/*============================================================================
 * Base macro and typedef definitions for system portability
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2008  EDF

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

#include "fvmc_config.h"
#include "config_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Basic macros
 *============================================================================*/


#if defined(__INTEL_COMPILER)
#define FVMC_PRAGMA_TO_STR(x) _Pragma(#x)
#define FVMC_INTEL_SUPPRESS_WARNING_PUSH _Pragma("warning(push)")
#define FVMC_INTEL_SUPPRESS_WARNING(w) FVMC_PRAGMA_TO_STR(warning(disable:w))
#define FVMC_INTEL_SUPPRESS_WARNING_POP _Pragma("warning(pop)")
#define FVMC_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    FVMC_INTEL_SUPPRESS_WARNING_PUSH FVMC_INTEL_SUPPRESS_WARNING(w)
#else // FVMC_INTEL
#define FVMC_INTEL_SUPPRESS_WARNING_PUSH
#define FVMC_INTEL_SUPPRESS_WARNING(w)
#define FVMC_INTEL_SUPPRESS_WARNING_POP
#define FVMC_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // FVMC_INTEL

#if defined(__clang__)
#define FVMC_PRAGMA_TO_STR(x) _Pragma(#x)
#define FVMC_CLANG_SUPPRESS_WARNING_PUSH _Pragma("clang diagnostic push")
#define FVMC_CLANG_SUPPRESS_WARNING(w) FVMC_PRAGMA_TO_STR(clang diagnostic ignored w)
#define FVMC_CLANG_SUPPRESS_WARNING_POP _Pragma("clang diagnostic pop")
#define FVMC_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    FVMC_CLANG_SUPPRESS_WARNING_PUSH FVMC_CLANG_SUPPRESS_WARNING(w)
#else // FVMC_CLANG
#define FVMC_CLANG_SUPPRESS_WARNING_PUSH
#define FVMC_CLANG_SUPPRESS_WARNING(w)
#define FVMC_CLANG_SUPPRESS_WARNING_POP
#define FVMC_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // FVMC_CLANG

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#define FVMC_PRAGMA_TO_STR(x) _Pragma(#x)
#define FVMC_GCC_SUPPRESS_WARNING_PUSH _Pragma("GCC diagnostic push")
#define FVMC_GCC_SUPPRESS_WARNING(w) FVMC_PRAGMA_TO_STR(GCC diagnostic ignored w)
#define FVMC_GCC_SUPPRESS_WARNING_POP _Pragma("GCC diagnostic pop")
#define FVMC_GCC_SUPPRESS_WARNING_WITH_PUSH(w)                                                  \
    FVMC_GCC_SUPPRESS_WARNING_PUSH FVMC_GCC_SUPPRESS_WARNING(w)
#else // FVMC_GCC
#define FVMC_GCC_SUPPRESS_WARNING_PUSH
#define FVMC_GCC_SUPPRESS_WARNING(w)
#define FVMC_GCC_SUPPRESS_WARNING_POP
#define FVMC_GCC_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // DOCTEST_GCC

#define FVMC_UNUSED(x) (void)(x)


FVMC_GCC_SUPPRESS_WARNING("-Wcast-qual")
FVMC_GCC_SUPPRESS_WARNING("-Wunknown-pragmas")
FVMC_GCC_SUPPRESS_WARNING("-Wsign-compare")
FVMC_GCC_SUPPRESS_WARNING("-Wfloat-equal")


/*============================================================================
 * Internationalization (future)
 *============================================================================*/

#ifdef FVMC_HAVE_GETTEXT
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
#endif /* FVMC_HAVE_GETTEXT */

/*============================================================================
 * C99 Qualifiers
 *============================================================================*/

/* inline provided by fvmc_config.h if necessary */

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
# if (FVMC_SIZEOF_INT == 4)
typedef int int32_t;
# elif (FVMC_SIZEOF_SHORT == 4)
typedef short int32_t;
# else
#  error
# endif
#endif

/* int64_t type */

#if !defined(HAVE_INT64_T)
# if (FVMC_SIZEOF_INT == 8)
typedef int int64_t;
# elif (FVMC_SIZEOF_LONG == 8)
typedef long int64_t;
# elif (HAVE_LONG_LONG == 8)  /* FVMC_SIZEOF_LONG_LONG not generally available */
typedef long long int64_t;
# else
#  error
# endif
#endif

/* uint32_t type */

#if !defined(HAVE_UINT32_T)
# if (FVMC_SIZEOF_INT == 4)
typedef unsigned uint32_t;
# elif (FVMC_SIZEOF_SHORT == 4)
typedef unsigned short uint32_t;
# else
#  error
# endif
#endif

/* uint64_t type */

#if !defined(HAVE_UINT64_T)
# if (FVMC_SIZEOF_INT == 8)
typedef unsigned uint64_t;
# elif (FVMC_SIZEOF_LONG == 8)
typedef unsigned long uint64_t;
# elif (HAVE_LONG_LONG) /* FVMC_SIZEOF_LONG_LONG not generally available */
typedef unsigned long long uint64_t;
# else
#  error
# endif
#endif

/* Directory name separator ('/' for Unix/Linux, '\' for Windows, ':' for Mac */

#define FVMC_DIR_SEPARATOR '/'

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_CONFIG_DEFS_H__ */
