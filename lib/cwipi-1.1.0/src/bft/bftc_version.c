/*============================================================================
 * BFT version and build information
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

#include <string.h>

/*
 * Optional library and BFT headers
 */

#if defined(HAVE_ZLIB)
#include <zlib.h>
#endif

#include "bftc_version.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-------------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------
 * Local static strings
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Return BFT version version string
 *
 * \return Pointer to static string BFT version version string.
 */

const char *
bftc_version(void)
{
  static const char buf[] = PACKAGE_VERSION;

  return buf;
}

/*!
 * \brief Indicate Zlib version available at run time.
 *
 * It may be useful to compare the Zlib version used at compile
 * and run time in case we use dynamic libraries.
 *
 * \return pointer to string indicating Zlib version in use, or NULL
 *         if Zlib support is not available.
 */

const char *
bftc_version_zlib(void)
{
#if defined(HAVE_ZLIB)
  return zlibVersion();
#else
  return NULL;
#endif
}

/*!
 * \brief Indicate Zlib version available at compilation time.
 *
 * It may be useful to compare the Zlib version used at compile
 * and link time in case we use dynamic libraries.
 *
 * \return pointer to string indicating Zlib version at compilation, or NULL
 *         if Zlib support is not available.
 */

const char *
bftc_version_build_zlib(void)
{
#if defined(HAVE_ZLIB)
#if defined(ZLIB_VERSION)
  return ZLIB_VERSION;
#else
  return _("unknown");
#endif
#else
  return NULL;
#endif
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
