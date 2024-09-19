#ifndef __PDM_MEMORY_STATS_PRIV_H__
#define __PDM_MEMORY_STATS_PRIV_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_extract_part_t
 * \brief  Define a partition mesh. Arrays are shared
 *
 */

struct _pdm_memory_stats_t
{
  PDM_MPI_Comm           comm;
  int                    n_memory_snapshot;

  char**                 snapshot_name;
  long*                  curr_real_mem;
  long*                  peak_real_mem;
  long*                  curr_virt_mem;
  long*                  peak_virt_mem;

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MEMORY_STATS_PRIV_H__ */
