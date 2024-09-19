#ifndef __PDM_DCUBE_GEN_PRIV_H__
#define __PDM_DCUBE_GEN_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"

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

struct _pdm_dcube_t {
  PDM_MPI_Comm     comm;            /*!< MPI communicator                          */
  PDM_ownership_t  owner;           /*!< Which have the responsabilities of results*/
  PDM_g_num_t      n_vtx_seg;       /*!< Number of vertices in segments            */
  double           length;          /*!< Segment length                            */
  double           zero_x;          /*!< Coordinates of the origin                 */
  double           zero_y;          /*!< Coordinates of the origin                 */
  double           zero_z;          /*!< Coordinates of the origin                 */
  int              n_face_group;    /*!< Number of faces groups                    */
  int              dn_cell;         /*!< Number of cells stored in this process    */
  int              dn_face;         /*!< Number of faces stored in this process    */
  int              dn_vtx;          /*!< Number of vertices stored in this process */
  PDM_g_num_t     *dface_cell;      /*!< Faces from cells connectivity             */
  int             *dface_vtx_idx;   /*!< Faces from vertices connectivity index    */
  PDM_g_num_t     *dface_vtx;       /*!< Faces from vertices connectivity          */
  double          *dvtx_coord;      /*!< Vertices coordinates                      */
  int             *dface_group_idx; /*!< Faces groups index                        */
  PDM_g_num_t     *dface_group;     /*!< Faces groups                              */
} ;


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DCUBE_GEN_PRIV_H__ */
