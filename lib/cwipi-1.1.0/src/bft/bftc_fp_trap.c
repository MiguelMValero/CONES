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

#include "bftc_config_defs.h"

/*
 * Standard C library headers
 */

/*
 * Optional library and BFT headers
 */

#include "bftc_fp_trap.h"

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

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Enable trapping of floating point exceptions if possible.
 *
 * This function is deprecated, and does nothing.
 * Use C99's feenableexcept() instead.
 */

void
bftc_fp_trap_set(void)
{
}

/*!
 * \brief Revert to initial enable trapping behavior of floating point
 *  exceptions if possible.
 *
 * This function is deprecated, and does nothing.
 */

void
bftc_fp_trap_unset(void)
{
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

