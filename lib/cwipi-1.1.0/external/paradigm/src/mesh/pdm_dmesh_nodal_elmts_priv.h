#ifndef __PDM_DMESH_NODAL_ELMTS_PRIV_H__
#define __PDM_DMESH_NODAL_ELMTS_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */


/**
 * \struct PDM_Mesh_nodal_section_std_t
 * \brief  Standard geometric section
 *
 */

typedef struct PDM_DMesh_nodal_section_std_t {

  PDM_Mesh_nodal_elt_t    t_elt;       /*!< Element type */
  PDM_l_num_t             n_elt;       /*!< Number elements */
  PDM_g_num_t            *_connec;     /*!< Connectivity (Memory mapping)
                                        *   (size = Number of vertices per element * \ref n_elt)  */
  PDM_g_num_t            *distrib;     /*!< Distribution on the processes (size = \ref n_rank + 1) */
  int                     order;       /*!< Element order */
  const char             *ho_ordering; /*!< HO node ordering */
  PDM_ownership_t         owner;

} PDM_DMesh_nodal_section_std_t;


/**
 * \struct PDM_Mesh_nodal_section_poly2d_t
 * \brief  Polygon geometric section
 *
 */

typedef struct PDM_DMesh_nodal_section_poly2d_t {

  PDM_l_num_t     n_elt;         /*!< Number of elements of each partition */
  PDM_l_num_t     *_connec_idx;  /*!< Index of elements connectivity
                                  *  (Memory mapping) (size = \ref n_elt + 1) */
  PDM_g_num_t     *_connec;      /*!< Elements connectivity
                                  * (Memory mapping) (size = \ref connec_idx[\ref n_elt]) */
  PDM_g_num_t     *distrib;      /*!< Distribution on the processes (size = \ref n_rank + 1) */
  PDM_ownership_t  owner;

} PDM_DMesh_nodal_section_poly2d_t;


/**
 * \struct PDM_Mesh_nodal_section_poly3d_t
 * \brief  Polyhedron geometric section
 *
 */

typedef struct PDM_DMesh_nodal_section_poly3d_t{

  PDM_l_num_t     n_elt;          /*!< Number of elements */
  PDM_l_num_t     n_face;         /*!< Number of faces */
  PDM_l_num_t     *_face_vtx_idx; /*!< Index of faces connectivity
                                   * (Memory mapping) (Size = \ref n_face + 1) */

  PDM_g_num_t     *_face_vtx;      /*!< Faces connectivity
                                    * (Memory mapping) (Size = \ref _face_vtx_idx[\ref n_face])*/
  PDM_l_num_t     *_cell_face_idx; /*!< Index of cell->face connectivity
                                    * (Memory mapping) (Size = \ref n_cell + 1) */

  PDM_g_num_t     *_cell_face;     /*!< cell->face connectivity
                                    * (Memory mapping) (Size = \ref _cell_face_idx[\ref n_cell]) */
  PDM_g_num_t     *distrib;        /*!< Distribution on the processes (size = \ref n_rank + 1) */
  PDM_ownership_t  owner;

} PDM_DMesh_nodal_section_poly3d_t;


typedef struct _pdm_dmesh_nodal_elts_t {

  PDM_MPI_Comm                       comm;                   /*!< MPI Communicator             */
  int                                n_rank;                 /*!< Number of processes          */
  int                                i_rank;                 /*!< Number of processes          */
  int                                mesh_dimension;         /*! Principal dimension of meshes */
  PDM_g_num_t                        n_g_elmts;              /*!< Global number of elements    */

  int                                n_section;              /*!< Total number of sections */
  int                                n_section_std;          /*!< Total number of standard sections   */
  int                                n_section_poly2d;       /*!< Total number of polyhedron sections */
  int                                n_section_poly3d;       /*!< Total number of polyhedron sections */
  int                               *sections_id;

  PDM_DMesh_nodal_section_std_t    **sections_std;           /*!< Standard sections            */
  PDM_DMesh_nodal_section_poly2d_t **sections_poly2d;        /*!< Polygon sections             */
  PDM_DMesh_nodal_section_poly3d_t **sections_poly3d;        /*!< Polyhedron sections          */
  PDM_g_num_t                       *section_distribution;   /*!< Element distribution         */

  int              n_group_elmt;
  int             *dgroup_elmt_idx;
  PDM_g_num_t     *dgroup_elmt;
  PDM_ownership_t  dgroup_elmt_owner;

  PDM_g_num_t     *dparent_gnum;
  int             *dparent_sign;
  int             *dparent_idx;
  PDM_g_num_t     *delmt_child_distrib;


} _pdm_dmesh_nodal_elts_t;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_NODAL_ELMTS_PRIV_H__ */
