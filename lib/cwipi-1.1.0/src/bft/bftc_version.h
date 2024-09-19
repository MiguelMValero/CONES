#ifndef __BFTC_VERSION_H__
#define __BFTC_VERSION_H__

/*============================================================================
 * BFT library information
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

/*----------------------------------------------------------------------------*/

#include "bftc_config.h"

/*----------------------------------------------------------------------------*/

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
 * Return BFT library version
 *
 * returns:
 *   Pointer to static string BFT library version string.
 */

const char *
bftc_version(void);

/*
 * Indicate Zlib version available at run time.
 *
 * It may be useful to compare the Zlib version used at compile
 * and run time in case we use dynamic libraries.
 *
 * returns:
 *   pointer to string indicating Zlib version in use, or NULL
 *   if Zlib support is not available.
 */

const char *
bftc_version_zlib(void);

/*
 * Indicate Zlib version available at compilation time.
 *
 * It may be useful to compare the Zlib version used at compile
 * and link time in case we use dynamic libraries.
 *
 * returns:
 *   pointer to string indicating Zlib version at compilation, or NULL
 *   if Zlib support is not available.
 */

const char *
bftc_version_build_zlib(void);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BFTC_VERSION_H__ */
