/*
 * \file
 */

#ifndef __PDM_gen_GNUM_H__
#define __PDM_gen_GNUM_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017-2023       ONERA

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

typedef struct _pdm_gen_gnum_t PDM_gen_gnum_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a structure to build a global numbering
 *
 * \param [in]   dim          Spatial dimension
 * \param [in]   n_part       Number of local partitions
 * \param [in]   merge        Merge double points or not
 * \param [in]   tolerance    Geometric tolerance (used if \p merge is set to \ref PDM_TRUE)
 * \param [in]   comm         PDM_MPI communicator
 * \param [in]   owner        Ownership of results
 *
 * \return     Pointer to \ref PDM_gen_gnum_t object
 */

PDM_gen_gnum_t *
PDM_gnum_create
(
 const int             dim,
 const int             n_part,
 const PDM_bool_t      merge,
 const double          tolerance,
 const PDM_MPI_Comm    comm,
 const PDM_ownership_t owner
);



/**
 *
 * \brief Set from coordinates
 *
 * \note The ordering is based on the <a href="https://en.wikipedia.org/wiki/Z-order_curve">Morton space-filling curve</a>.
 *       Elements are expected to be unique (i.e. not duplicated on 2 or more ranks).
 *       If two elements share the same Morton code, their global
 *       id will be determined by lexicographical ordering of coordinates.
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum_t object
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   coords       Coordinates (size = 3 * \p n_elts)
 * \param [in]   char_length  Characteristic length (or NULL)
 *                            (used only if \p merge is set to \ref PDM_TRUE)
 *
 */

void
PDM_gnum_set_from_coords
(
       PDM_gen_gnum_t *gen_gnum,
 const int             i_part,
 const int             n_elts,
 const double         *coords,
 const double         *char_length
);



/**
 *
 * \brief Set parent global numbering
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum_t object
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   parent_gnum  Parent global ids (size = \p n_elts)
 *
 */

void
PDM_gnum_set_from_parents
(
       PDM_gen_gnum_t  *gen_gnum,
 const int              i_part,
 const int              n_elts,
 const PDM_g_num_t     *parent_gnum
);

/**
 *
 * \brief Set size of tuple for nuplet
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum_t object
 * \param [in]   nuplet       Size of tuple
 *
 */
void
PDM_gnum_set_parents_nuplet
(
       PDM_gen_gnum_t  *gen_gnum,
 const int              nuplet
);

/**
 *
 * \brief Build global numbering
 *
 * \param [in]   gen_gnum         Pointer to \ref PDM_gen_gnum_t object
 *
 */

void
PDM_gnum_compute
(
 PDM_gen_gnum_t  *gen_gnum
);


/**
 *
 * \brief Get global ids for a given partition
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum_t object
 * \param [in]   i_part       Current partition
 *
 * \return     Array of global ids
 *
 */

PDM_g_num_t *
PDM_gnum_get
(
       PDM_gen_gnum_t *gen_gnum,
 const int             i_part
);


/**
 *
 * \brief Free a \ref PDM_gen_gnum_t object
 *
 * \param [in]   gen_gnum         Pointer to \ref PDM_gen_gnum_t object
 *
 */

void
PDM_gnum_free
(
PDM_gen_gnum_t *gen_gnum
);


/**
 *
 * \brief Get number of elements in a partition
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum_t object
 * \param [in]   i_part       Current partition
 *
 * \return     Number of elements in current partition
 *
 */

int
PDM_gnum_n_elt_get
(
       PDM_gen_gnum_t *gen_gnum,
 const int             i_part
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_gen_GNUM_H__ */
