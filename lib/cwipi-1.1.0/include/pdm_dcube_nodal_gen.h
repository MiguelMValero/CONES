/*
 * \file
 */

#ifndef __PDM_DCUBE_NODAL_GEN_H__
#define __PDM_DCUBE_NODAL_GEN_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_ho_ordering.h"
#include "pdm_domain_interface.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_dcube_nodal_t PDM_dcube_nodal_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Initialize a \ref PDM_dcube_nodal_t structure
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   n_vtx_x        Number of vertices on segments in x-direction
 * \param [in]   n_vtx_y        Number of vertices on segments in y-direction
 * \param [in]   n_vtx_z        Number of vertices on segments in z-direction
 * \param [in]   length         Segment length
 * \param [in]   zero_x         x-coordinate of the origin
 * \param [in]   zero_y         y-coordinate of the origin
 * \param [in]   zero_z         z-coordinate of the origin
 * \param [in]   t_elt          Element type
 * \param [in]   order          Element order
 * \param [in]   owner          Ownership
 *
 * \return   Pointer to new \ref PDM_dcube_nodal_t object
 *
 */

PDM_dcube_nodal_t *
PDM_dcube_nodal_gen_create
(
 PDM_MPI_Comm          comm,
 const PDM_g_num_t     n_vtx_x,
 const PDM_g_num_t     n_vtx_y,
 const PDM_g_num_t     n_vtx_z,
 const double          length,
 const double          zero_x,
 const double          zero_y,
 const double          zero_z,
 PDM_Mesh_nodal_elt_t  t_elt,
 const int             order,
 PDM_ownership_t       owner
);


/**
 *
 * \brief Free a \ref PDM_dcube_nodal_t structure
 *
 * \param [in]  dcube      Pointer to \ref PDM_dcube_nodal_t object
 *
 */

void
PDM_dcube_nodal_gen_free
(
 PDM_dcube_nodal_t *dcube
);


/**
 *
 * \brief Set the HO-ordering for a \ref PDM_dcube_nodal_t structure
 *
 * \param [in]  dcube      Pointer to \ref PDM_dcube_nodal_t object
 * \param [in]  ordering   Name of the HO-ordering
 *
 */

void PDM_dcube_nodal_gen_ordering_set
(
       PDM_dcube_nodal_t *dcube,
 const char              *ordering
);


/**
 *
 * \brief Build a \ref PDM_dcube_nodal_t structure
 *
 * \param [in]  dcube      Pointer to \ref PDM_dcube_nodal_t object
 *
 * \return   Pointer to the associated \ref PDM_dmesh_nodal_t object
 *
 */

PDM_dmesh_nodal_t *
PDM_dcube_nodal_gen_build
(
 PDM_dcube_nodal_t *dcube
);


/**
 *
 * \brief Get the \ref PDM_dmesh_nodal_t associated to a \ref PDM_dcube_nodal_t
 *
 * \param [in]  dcube      Pointer to \ref PDM_dcube_nodal_t object
 *
 * \return   Pointer to the associated \ref PDM_dmesh_nodal_t object
 *
 */

PDM_dmesh_nodal_t *
PDM_dcube_nodal_gen_dmesh_nodal_get
(
 PDM_dcube_nodal_t *dcube
);


/**
 *
 * \brief Generate a multi-domain, 3D array of \ref PDM_dcube_nodal_t
 *        with (possibly periodic) inter-domain joints
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   n_dom_i        Number of domains in 1st dimension
 * \param [in]   n_dom_j        Number of domains in 2nd dimension
 * \param [in]   n_dom_k        Number of domains in 3rd dimension
 * \param [in]   periodic_i     Enable periodic joint in 1st dimension (0 or 1)
 * \param [in]   periodic_j     Enable periodic joint in 2nd dimension (0 or 1)
 * \param [in]   periodic_k     Enable periodic joint in 3rd dimension (0 or 1)
 * \param [in]   n_vtx_x_in     Number of vertices on segments in x-direction (for each domain)
 * \param [in]   n_vtx_y_in     Number of vertices on segments in y-direction (for each domain)
 * \param [in]   n_vtx_z_in     Number of vertices on segments in z-direction (for each domain)
 * \param [in]   length         Segment length (for each domain)
 * \param [in]   zero_x         x-coordinate of the origin
 * \param [in]   zero_y         y-coordinate of the origin
 * \param [in]   zero_z         z-coordinate of the origin
 * \param [in]   t_elt          Element type
 * \param [in]   order          Element order
 * \param [in]   owner          Ownership
 * \param [out]  dcube          Set of \ref PDM_dcube_nodal_t objects
 * \param [out]  dom_intrf      Inter-domain joints
 * \param [in]   owner          Ownership
 *
 */

void
PDM_dcube_nodal_cart_topo
(
       PDM_MPI_Comm              comm,
       int                       n_dom_i,
       int                       n_dom_j,
       int                       n_dom_k,
       int                       periodic_i,
       int                       periodic_j,
       int                       periodic_k,
 const PDM_g_num_t               n_vtx_x_in,
 const PDM_g_num_t               n_vtx_y_in,
 const PDM_g_num_t               n_vtx_z_in,
 const double                    length,
 const double                    zero_x,
 const double                    zero_y,
 const double                    zero_z,
       PDM_Mesh_nodal_elt_t      t_elt,
 const int                       order,
       PDM_dcube_nodal_t      ***dcube,
       PDM_domain_interface_t  **dom_intrf,
       PDM_ownership_t           owner
);


void
PDM_generate_lines
(
  PDM_MPI_Comm  comm,
  double        zero_x,
  double        zero_y,
  double        zero_z,
  double        length,
  PDM_g_num_t   n_g_pts,
  PDM_g_num_t **distrib_edge_out,
  PDM_g_num_t **distrib_vtx_out,
  PDM_g_num_t **dedge_vtx_out,
  double      **dvtx_coord_out
);


void
PDM_generate_cart_topo_lines
(
 PDM_MPI_Comm              comm,
 int                       n_dom_i,
 int                       periodic_i,
 double                    zero_x,
 double                    zero_y,
 double                    zero_z,
 double                    length,
 PDM_g_num_t               n_g_pts,
 PDM_g_num_t            ***distrib_edge_out,
 PDM_g_num_t            ***distrib_vtx_out,
 PDM_g_num_t            ***dedge_vtx_out,
 double                 ***dvtx_coord_out,
 PDM_domain_interface_t  **dom_intrf
);


/**
 * \brief Set randomization factor
 *
 * \param [in]  dcube          Pointer to \ref PDM_dcube_nodal_t object
 * \param [in]  random_factor  Randomization factor (between 0 and 1)
 *
 */

void
PDM_dcube_nodal_gen_random_factor_set
(
 PDM_dcube_nodal_t *dcube,
 double             random_factor
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DCUBE_NODAL_GEN_H__ */
