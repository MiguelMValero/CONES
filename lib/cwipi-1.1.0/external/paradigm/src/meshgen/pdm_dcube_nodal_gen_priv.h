#ifndef __PDM_DCUBE_NODAL_GEN_PRIV_H__
#define __PDM_DCUBE_NODAL_GEN_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_ho_ordering.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_dcube_t
 * \brief  Distributed cube
 *
 * _dcube_t define a distributed mesh of a cube
 *
 */

struct _pdm_dcube_nodal_t {
  PDM_MPI_Comm          comm;                   /*!< MPI communicator                          */
  PDM_ownership_t       owner;                  /*!< Which have the responsabilities of results*/

  PDM_g_num_t           nx;                     /*!< Number of elements in segments along x    */
  PDM_g_num_t           ny;                     /*!< Number of elements in segments along y    */
  PDM_g_num_t           nz;                     /*!< Number of elements in segments along z    */
  double                length;                 /*!< Segment length                            */
  double                zero_x;                 /*!< Coordinates of the origin                 */
  double                zero_y;                 /*!< Coordinates of the origin                 */
  double                zero_z;                 /*!< Coordinates of the origin                 */

  PDM_Mesh_nodal_elt_t  t_elt;                  /*!< Type of elements to generate              */
  int                   order;                  /*!< Order of elements                         */

  double                random_factor;          /*!< Randomization factor                      */

  PDM_g_num_t          *distrib_hexa;
  PDM_g_num_t          *distrib_quad;
  PDM_g_num_t          *distrib_bar;

  int                   dn_hexa;
  int                   dn_quad;
  int                   dn_bar;


  PDM_dmesh_nodal_t    *dmesh_nodal;            /*!< Results                                   */

  char                 *ordering;
} ;


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DCUBE_NODAL_GEN_PRIV_H__ */
