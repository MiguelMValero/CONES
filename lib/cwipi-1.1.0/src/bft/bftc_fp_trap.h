#ifndef __BFTC_FP_TRAP_H__
#define __BFTC_FP_TRAP_H__

/*============================================================================
 * Trap floating-point exceptions on some systems
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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Enable trapping of floating point exceptions if possible.
 *
 * This function is deprecated, and does nothing.
 * Use C99's feenableexcept() instead.
 */

void
bftc_fp_trap_set(void);

/*
 * Revert to initial enable trapping behavior of floating point
 * exceptions if possible.
 *
 * This function is deprecated, and does nothing.
 */

void
bftc_fp_trap_unset(void);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BFTC_FP_TRAP_H__ */
