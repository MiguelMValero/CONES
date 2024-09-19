#ifndef __PDM_MULTIPART_PRIV_H__
#define __PDM_MULTIPART_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multipart.h"
#include "pdm_dmesh_priv.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_part_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/


/**
 * \struct _part_mesh_t
 * \brief  This private structure stores partionned meshes obtained on a given
 *         domain.
 */
typedef struct {
  /* Data shared between all the partitions */
  int  tn_part;          // Total number of partitions created for this mesh
  int  n_bounds;         // Global number of boundary groups
  int  n_joins;          // Global number of interface groups
  int *joins_ids;        // Global id of each interface (size=n_joins)

  /* Partitions -- see pdm_part_priv.h for struct definition */
  // _part_t  **parts;
  PDM_part_mesh_t *pmesh;
  int        renum_cell_method;     // Choice of renumbering method for cells
  const int *renum_cell_properties; // Parameters used by some renumbering methods
  int        renum_face_method;     // Choice of renumbering method for faces
  int        renum_edge_method;     // Choice of renumbering method for faces
  int        renum_vtx_method;      // Choice of renumbering method for faces

  /*
   * Additional info
   */
  int        **vtx_ghost_information;
  PDM_bool_t   is_owner_vtx_ghost_information;

  int        **hyperplane_color;
  int        **thread_color;
  PDM_bool_t   is_owner_hyperplane_color;
  PDM_bool_t   is_owner_thread_color;

} _part_mesh_t;


/**
 * \struct _pdm_multipart_t
 * \brief  This structure describe a multipart.
 *         It includes distributed meshes, partioned meshes and
 *         partitioning parameters.
 */
struct _pdm_multipart_t {
  /* Multipart description */
  int           n_domain;            // Number of initial domains

  PDM_dmesh_t **dmeshes;             // Ids of dmesh structure storing
                                     // distributed meshes (size = n_domain)
  PDM_bool_t   *is_owner_dmeshes;
  PDM_dmesh_nodal_t          **dmeshes_nodal;
  PDM_dmesh_nodal_to_dmesh_t **dmn_to_dm;

  int           n_total_joins;       // Total number of joins between domains (each counts twice)
  const int    *join_to_opposite;    // For each global joinId, give the globalId of
                                     //   the opposite join (size = n_total_joins)

  PDM_MPI_Comm     comm;             // MPI communicator
  PDM_ownership_t  owner;            // Which have the responsabilities of results

  /* Partitioning parameters */
  PDM_bool_t       merge_blocks;     // Merge before partitionning or not
  PDM_split_dual_t split_method;     // Partitioning method (Metis or Scotch)
  PDM_part_size_t  part_size_method; // Procude homogeneous or heterogeneous partitions
  int             *n_part;           // Number of wanted partitions per proc
                                     // in each domain (size = n_domain)
  const double    *part_fraction;    // Weight (in %) of each partition, in each domain
                                     //   (size = sum n_part[i]), if heterogeneous
  /* Partitioned meshes */
  _part_mesh_t *pmeshes;             // Partitioned meshes structures (size=n_domain)

};



/*============================================================================
 * Private functions
 *============================================================================*/
void
_run_ppart_domain_nodal
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  _part_mesh_t      *pmesh,
  PDM_split_dual_t   split_method,
  int                dn_part,
  PDM_MPI_Comm       comm
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTIPART_PRIV_H__ */
