/*
 * \file
 */

#ifndef __PDM_GNUM_FROM_HASH_VALUES_H__
#define __PDM_GNUM_FROM_HASH_VALUES_H__

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

typedef struct _pdm_gnum_from_hv_t PDM_gnum_from_hv_t;

typedef int (*gnum_from_hv_compare)(const void* a, const void* b, void* );
typedef int (*gnum_from_hv_equal  )(const void* a, const void* b, void* );

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a global numbering structure
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   equilibrate  Use algorithm to equilibrate the block treatment (hash value is not a priori equi-reparti)
 * \param [in]   tolerance    Geometric tolerance (used if merge double points is activated)
 * \param [in]   comm         PDM_MPI communicator
 * \param [in]   owner        Owner
 *
 * \return     Pointer to \ref PDM_gnum_from_hv_t object
 */

PDM_gnum_from_hv_t *
PDM_gnum_from_hash_values_create
(
 const int            n_part,
 const PDM_bool_t     equilibrate,
 const size_t         s_data,
 gnum_from_hv_compare fcompare,
 gnum_from_hv_equal   fequal,
 const PDM_MPI_Comm   comm,
 const PDM_ownership_t owner
);


/**
 *
 * \brief Set hash values for one partition
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   part_hkey    For each elements the hash value associated
 * \param [in]   part_strid   Stride between each data in part_hdata
 * \param [in]   part_hdata   Partition data which compute the hash value, we need it to setup in a block way
 *
 */

void
PDM_gnum_set_hash_values
(
 PDM_gnum_from_hv_t  *gnum_from_hv,
 const int            i_part,
 const int            n_elts,
 const size_t        *part_hkeys,
 const int           *part_hstri,
 const unsigned char *part_hdata
);


/**
 *
 * \brief Compute
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 *
 */

void
PDM_gnum_from_hv_compute
(
 PDM_gnum_from_hv_t *gnum_from_hv
);


/**
 *
 * \brief Get the global ids for the current partition
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 * \param [in]   i_part       Current partition
 *
 * \return  Array of global ids for the current partition
 *
 */

PDM_g_num_t *
PDM_gnum_from_hv_get
(
 PDM_gnum_from_hv_t *gnum_from_hv,
 const int           i_part
);


/**
 *
 * \brief Dump elapsed an CPU time
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 *
 */

void
PDM_gnum_from_hv_dump_times
(
 PDM_gnum_from_hv_t *gnum_from_hv
);


/**
 *
 * \brief Free
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 * \param [in]   partial      1 to free partially, 0 else
 *
 */

void
PDM_gnum_from_hv_free
(
 PDM_gnum_from_hv_t *gnum_from_hv
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GNUM_FROM_HASH_VALUES_H__ */
