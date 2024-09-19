#ifndef __CWIPI_PRIV_H__
#define __CWIPI_PRIV_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-20  ONERA

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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Suppress warning
 *============================================================================*/

#if defined(__INTEL_COMPILER)
#define CWIPI_PRAGMA_TO_STR(x) _Pragma(#x)
#define CWIPI_INTEL_SUPPRESS_WARNING_PUSH _Pragma("warning(push)")
#define CWIPI_INTEL_SUPPRESS_WARNING(w) CWIPI_PRAGMA_TO_STR(warning(disable:w))
#define CWIPI_INTEL_SUPPRESS_WARNING_POP _Pragma("warning(pop)")
#define CWIPI_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    CWIPI_INTEL_SUPPRESS_WARNING_PUSH CWIPI_INTEL_SUPPRESS_WARNING(w)
#else // CWIPI_INTEL
#define CWIPI_INTEL_SUPPRESS_WARNING_PUSH
#define CWIPI_INTEL_SUPPRESS_WARNING(w)
#define CWIPI_INTEL_SUPPRESS_WARNING_POP
#define CWIPI_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // CWIPI_INTEL

#if defined(__clang__)
#define CWIPI_PRAGMA_TO_STR(x) _Pragma(#x)
#define CWIPI_CLANG_SUPPRESS_WARNING_PUSH _Pragma("clang diagnostic push")
#define CWIPI_CLANG_SUPPRESS_WARNING(w) CWIPI_PRAGMA_TO_STR(clang diagnostic ignored w)
#define CWIPI_CLANG_SUPPRESS_WARNING_POP _Pragma("clang diagnostic pop")
#define CWIPI_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    CWIPI_CLANG_SUPPRESS_WARNING_PUSH CWIPI_CLANG_SUPPRESS_WARNING(w)
#else // CWIPI_CLANG
#define CWIPI_CLANG_SUPPRESS_WARNING_PUSH
#define CWIPI_CLANG_SUPPRESS_WARNING(w)
#define CWIPI_CLANG_SUPPRESS_WARNING_POP
#define CWIPI_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // CWIPI_CLANG

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#define CWIPI_PRAGMA_TO_STR(x) _Pragma(#x)
#define CWIPI_GCC_SUPPRESS_WARNING_PUSH _Pragma("GCC diagnostic push")
#define CWIPI_GCC_SUPPRESS_WARNING(w) CWIPI_PRAGMA_TO_STR(GCC diagnostic ignored w)
#define CWIPI_GCC_SUPPRESS_WARNING_POP _Pragma("GCC diagnostic pop")
#define CWIPI_GCC_SUPPRESS_WARNING_WITH_PUSH(w)                                                  \
    CWIPI_GCC_SUPPRESS_WARNING_PUSH CWIPI_GCC_SUPPRESS_WARNING(w)
#else // CWIPI_GCC
#define CWIPI_GCC_SUPPRESS_WARNING_PUSH
#define CWIPI_GCC_SUPPRESS_WARNING(w)
#define CWIPI_GCC_SUPPRESS_WARNING_POP
#define CWIPI_GCC_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // DOCTEST_GCC

#define CWIPI_UNUSED(x) (void)(x)


CWIPI_GCC_SUPPRESS_WARNING("-Wcast-qual")
CWIPI_GCC_SUPPRESS_WARNING("-Wunknown-pragmas")


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWIPI_H__ */
