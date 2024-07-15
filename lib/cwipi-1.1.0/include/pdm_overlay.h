/*
 * \file
 */

#ifndef __PDM_OVERLAY_H__
#define __PDM_OVERLAY_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"

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
 * Type
 *============================================================================*/

typedef struct _PDM_ol_t PDM_ol_t;

/**
 * \enum PDM_ol_mesh_t
 * \brief 3 Meshes
 *
 */

typedef enum {

  PDM_OL_MESH_A       = 0,  /*!< First mesh to overlay */
  PDM_OL_MESH_B       = 1,  /*!< Second mesh to overlay */

} PDM_ol_mesh_t;


/**
 * \enum PDM_ol_parameter_t
 * \brief Parameters for ovelay meshes building
 *
 */

typedef enum {

  PDM_OL_CAR_LENGTH_TOL = 0,  /*!< Absolute tolerance for caracteristic length */
  PDM_OL_EXTENTS_TOL    = 1,  /*!< Absolute tolerance for extents */
  PDM_OL_SAME_PLANE_TOL = 2,  /*!< Absolute tolerance for check if 2 surfaces are
                                   the same plane surface*/

} PDM_ol_parameter_t;

/**
 * \enum PDM_ol_mv_t
 * \brief Type of moving mesh
 *
 */

typedef enum {

  PDM_OL_MV_TRANSFORMATION  = 0,  /*!< Moving with combination of geometric transformations */
  PDM_OL_MV_UNKNOWN         = 1,  /*!< Unknown moving type */

} PDM_ol_mv_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Build and initialize an overlaying object
 *
 * This function builds an initializes an overlaying surface meshes object
 *
 * \param [in]  n_partMeshA   Number of local partitions of the meshA input
 * \param [in]  n_partMeshB   Number of local partitions of the meshB input
 * \param [in]  projectCoeff Projection coefficient to define the overlay surface projection
 *                           If value == 0, the surface projection is MeshA
 *                           If value == 1, the surface projection is MeshB
 *                           If 0 < value < 1 , the projection surface is an intermediate surface
 * \param [in]  comm         MPI communicator.
 *
 * \return      id           Overlay object identifier.
 *
 */

PDM_ol_t *
PDM_ol_create
(
 const int          n_partMeshA,
 const int          n_partMeshB,
 const double       projectCoeff,
 const PDM_MPI_Comm comm
);


/**
 * \brief Set an overlay parameter
 *
 * This function sets en overlay parameter
 *
 * \param [in]  id          PDM_ol identifier
 * \param [in]  parameter   Parameter to define
 * \param [in]  value       Parameter value
 *
 */

void
PDM_ol_parameter_set
(
       PDM_ol_t           *ol,
 const PDM_ol_parameter_t  parameter,
 const double              value
);


/**
 * \brief Define input meshes properties
 *
 * This function defines the input meshes properties
 *
 * \param [in]  id          PDM_ol identifier
 * \param [in]  mesh        Input mesh to define
 *                          (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  i_part       Partition to define
 * \param [in]  n_face       Number of faces
 * \param [in]  face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]  face_vtx_idx  face -> vertex connectivity
 * \param [in]  face_ln_to_gn  Local face numbering to global face numbering
 * \param [in]  n_vtx        Number of vertices
 * \param [in]  coords      Coordinates
 * \param [in]  vtx_ln_to_gn   Local vertex numbering to global vertex numbering
 *
 */

void
PDM_ol_input_mesh_set
(
       PDM_ol_t      *ol,
 const PDM_ol_mesh_t  mesh,
 const int            i_part,
 const int            n_face,
 const int           *face_vtx_idx,
 const int           *face_vtx,
 const PDM_g_num_t   *face_ln_to_gn,
 const int            n_vtx,
 const double        *coords,
 const PDM_g_num_t   *vtx_ln_to_gn
);



/**
 * \brief Define the type of a mesh moving
 *
 * This function defines the type of a mesh moving.
 * Only a mesh can move
 *
 * \param [in]  id       PDM_ol identifier
 * \param [in]  mesh     Moving mesh
 * \param [in]  mv       Type of moving
 *
 */

void
PDM_ol_moving_type_set
(
       PDM_ol_t      *ol,
 const PDM_ol_mesh_t  mesh,
 const PDM_ol_mv_t    mv
);


/**
 * \brief Define a translation
 *
 * This function defines a translation for the moving mesh
 *
 * \param [in]  id       PDM_overlay identifier
 * \param [in]  vect     Translation vector
 * \param [in]  center   Translation center
 *
 */

void
PDM_ol_translation_set
(
       PDM_ol_t     *ol,
 const double       *vect,
 const double       *center
);



/**
 * \brief Define a rotation
 *
 * This function defines a rotation for the moving mesh
 *
 * \param [in]  id        PDM_ol identifier
 * \param [in]  direction Rotation direction
 * \param [in]  center    Rotation center
 * \param [in]  angle     Rotation center (degrees)
 *
 */

void
PDM_ol_rotation_set
(
       PDM_ol_t *ol,
 const double   *direction,
 const double   *center,
 const double    angle
);


/**
 * \brief Overlaying the input surface meshes
 *
 * This function overlays the input surface meshes
 *
 * \param [in]  id       PDM_ol identifier
 *
 */

void
PDM_ol_compute
(
 PDM_ol_t *ol
);


/**
 * \brief Return the entitie sizes of the overlay mesh
 *
 * This function returns the entities sizes of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id        PDM_ol identifier
 * \param [in]  mesh      Input mesh
 *                        (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [out] nGOlFace  Global number of faces of the overlay mesh
 * \param [out] nGOlVtx   Global number of vertices of the overlay mesh
 *
 */


void
PDM_ol_mesh_dim_get
(
 const PDM_ol_t      *ol,
 const PDM_ol_mesh_t  mesh,
       PDM_g_num_t   *nGOlFace,
       PDM_g_num_t   *nGOlVtx
);


/**
 * \brief Return the entitie sizes of the overlay mesh
 *
 * This function returns the entities sizes of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id            PDM_ol identifier
 * \param [in]  mesh          Input mesh
 *                            (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  i_part         Partition to define
 * \param [out] nOlFace        Number of faces of the overlay mesh
 * \param [out] nOlLinkedFace  Number of linked faces
 * \param [out] nOlVtx         Number of vertices of the overlay mesh
 * \param [out] sOlface_vtx    Size of olface_vtx for each partition
 * \param [out] sInitToOlFace  Size of initToOlFace for each partition
 *
 */

void
PDM_ol_part_mesh_dim_get
(
 const PDM_ol_t      *ol,
 const PDM_ol_mesh_t  mesh,
 const int            i_part,
       int           *nOlFace,
       int           *nOlLinkedFace,
       int           *nOlVtx,
       int           *sOlFaceIniVtx,
       int           *sOlface_vtx,
       int           *sInitToOlFace
);


/**
 * \brief Return the entitie of the overlay mesh
 *
 * This function returns the entities of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id              PDM_ol identifier
 * \param [in]  mesh            Input mesh
 *                              (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  i_part           Mesh partition identifier
 * \param [out] olface_vtx_idx    Array adress of \ref olface_vtx index
 *                              (size : \ref nOlFace + 1)
 * \param [out] olface_vtx       Array adress of face vertex connectivity
 *                              (size : \ref sOlface_vtx[\ref i_part])
 * \param [out] olLinkedface_procIdx olLinkedFace Index (size = n_proc + 1)
 * \param [out] olLinkedFace    Array adress of linked face in other mesh
 *                              For each face, 4 link properties :
 *                                    - local face number
 *                                    - linked process,
 *                                    - linked part number,
 *                                    - linked local face number
 *                              (size : \ref 4 * nOlLinkedFace)
 * \param [out] olface_ln_to_gn    Array adress of local to global face numbering
 *                              (size : \ref nOlFace)
 * \param [out] olCoords        Array adress of vertex coodinates
 *                              (size : 3 * \ref nOlVtx)
 * \param [out] olvtx_ln_to_gn     Array adress of local to global vertex numbering array
 *                              (size : \ref nOlVtx)
 * \param [out] initToOlFaceIdx Array adress of \ref initToOlFace index
 *                              (size : \ref nOlVtx + 1)
 * \param [out] initToOlFace    Array adress of initial to overlay faces
 * \param [out] initToOlVtx     Array adress of initial to overlay vertices
 *
 */

void
PDM_ol_mesh_entities_get
(
 const PDM_ol_t        *ol,
 const PDM_ol_mesh_t    mesh,
 const int              i_part,
 int                  **olFaceIniVtxIdx,
 int                  **olFaceIniVtx,
 int                  **olface_vtx_idx,
 int                  **olface_vtx,
 int                  **olLinkedface_procIdx,
 int                  **olLinkedFace,
 PDM_g_num_t          **olface_ln_to_gn,
 double               **olCoords,
 PDM_g_num_t          **olvtx_ln_to_gn,
 int                  **initToOlFaceIdx,
 int                  **initToOlFace
);


/**
 * \brief Delete an overlay object
 *
 * This function deletes an overlay object
 *
 * \param [in]  id                PDM_ol identifier.
 *
 */

void
PDM_ol_del
(
 PDM_ol_t *ol
);


/**
 * \brief Dump elapsed an CPU time
 *
 *
 * \param [in]  id                PDM_ol identifier.
 *
 */

void
PDM_ol_dump_times
(
 PDM_ol_t *ol
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_OVERLAY_H__ */
