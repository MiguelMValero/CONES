/*
 * \file
 */

#ifndef __PDM_SURF_MESH_H__
#define __PDM_SURF_MESH_H__

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

/**
 * \struct PDM_surf_mesh_t
 * \brief  Surface mesh
 *
 *  PDM_surf_mesh_t defines a surface mesh
 *
 */

typedef struct _pdm_surf_mesh_t PDM_surf_mesh_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_surf_mesh_t structure
 *
 * This function returns an initialized \ref PDM_surf_mesh_t structure
 *
 * \param [in]  n_part        Number of partition
 * \param [in]  comm         MSG communicator of mesh
 *
 * \return      A new initialized \ref PDM_surf_mesh_t structure
 *
 */

PDM_surf_mesh_t *
PDM_surf_mesh_create
(
const int         n_part,
PDM_MPI_Comm          comm
 );


/**
 * \brief Delete a \ref _mesh_t structure
 *
 * This function returns an initialized \ref _mesh_t structure
 *
 * \param [in]  mesh         Mesh to delete
 *
 * \return     Null pointer
 *
 */

PDM_surf_mesh_t *
PDM_surf_mesh_free
(
PDM_surf_mesh_t *mesh
 );


/**
 * \brief Compute edge global numbering
 *
 * This function computes edge global numbering
 *
 * \param [in]  mesh    mesh to compute global numbering
 *
 */

void
PDM_surf_mesh_build_edges_gn_and_edge_part_bound
(
 PDM_surf_mesh_t *mesh
 );


/**
 * \brief Compute inter partition communication graph between vertices
 *
 * This function computes edge global numbering
 *
 * \param [in]  mesh    mesh to compute global numbering
 *
 */

void
PDM_surf_mesh_build_vtx_part_bound
(
PDM_surf_mesh_t *mesh
 );


/**
 * \brief Build ghost faces and edges
 *
 * This function computes ghost edges and ghost faces
 *
 * \param [in]  mesh    mesh to compute global numbering
 *
 */

void
PDM_surf_mesh_build_ghost_element
(
 PDM_surf_mesh_t *mesh
 );


/**
 * \brief Build communication graph between internal partitions of any
 * initial mesh
 *
 * This function builds the communication graph between internal partitions
 * of each initial mesh
 *
 * \param [in]  mesh    mesh to compute global numbering
 *
 */

void
PDM_surf_mesh_build_exchange_graph
(
 PDM_surf_mesh_t *mesh
 );


/**
 * \brief Build ghost faces and ghost edges
 *
 * This function builds ghost faces and ghost edges
 * of each initial mesh
 *
 * \param [in]  mesh        Mesh
 *
 */

void
PDM_surf_mesh_compute_carLgthVtx
(
PDM_surf_mesh_t *mesh
);


/**
 * \brief Return face normal
 *
 * This function returns face normal after computation
 *
 * \param [in]  mesh       mesh
 *
 */

const double *
PDM_surf_mesh_face_normal_get
(
 PDM_surf_mesh_t  *mesh,
 int              i_part
);


/**
 * \brief Chek if mesh is a plane surfece
 *
 * This function cheks if the mesh is a plane surface
 * and returns plane equation ant vertices barycenter
 *
 * \param [in]  mesh           Mesh object
 * \param [in]  tolerance      Tolerance to accept surface as plane
 * \param [out] planeEquation  Plane equation
 * \param [out] barycenter     Vertices barycenter
 *
 */

int
PDM_surf_mesh_is_plane_surface
(
 PDM_surf_mesh_t  *mesh,
 double            tolerance,
 double            planeEquation[4],
 double            barycenter[3]
);


/**
 * \brief Compute face extents
 *
 * This function computes face extents
 *
 * \param [in]  mesh       Mesh object
 *
 */

void
PDM_surf_mesh_compute_faceExtentsMesh
(
 PDM_surf_mesh_t *mesh,
 double           tolerance
);


/**
 * \brief Compute face extents
 *
 * This function computes face extents
 *
 * \param [in]  mesh       Mesh object
 *
 */

void
PDM_surf_mesh_build_edges
(
 PDM_surf_mesh_t *mesh
);


/**
 * \brief Return global minimum of caracteristic length vertex
 *
 * This function returns global minimum of caracteristic length vertex
 *
 * \param [in]  mesh       Mesh object
 *
 */

double
PDM_surf_mesh_gMinCarLgthVtx_get
(
 PDM_surf_mesh_t *mesh
 );


/**
 * \brief Return Vertex caracteristic length
 *
 * This function returns global minimum of caracteristic length of vertex
 *
 * \param [in]  mesh       Mesh object
 *
 * \return  Vertex caracteristic length
 */

double *
PDM_surf_mesh_part_carLgthVtx_get
(
 PDM_surf_mesh_t *mesh,
 int              i_part
);


/**
 * \brief Input a partition
 *
 * This function inputs a partition
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part       Partition to define
 * \param [in]  n_face       Number of faces
 * \param [in]  face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]  face_vtx     face -> vertex connectivity
 * \param [in]  face_ln_to_gn  Local face numbering to global face numbering
 * \param [in]  n_vtx        Number of vertices
 * \param [in]  coords      Coordinates
 * \param [in]  vtx_ln_to_gn   Local vertex numbering to global vertex numbering
 *
 */

void
PDM_surf_mesh_part_input
(
 PDM_surf_mesh_t      *mesh,
 const int            i_part,
 const int            n_face,
 const int           *face_vtx_idx,
 const int           *face_vtx,
 const PDM_g_num_t    *face_ln_to_gn,
 const int            n_vtx,
 const double        *coords,
 const PDM_g_num_t    *vtx_ln_to_gn
);


/**
 * \brief Return number of partitions
 *
 * \param [in]  mesh       Mesh object
 *
 * \return    Number of partitions
 */

int
PDM_surf_mesh_n_part_get
(
 PDM_surf_mesh_t      *mesh
);


/**
 * \brief Return number of faces
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Number of faces
 */

int
PDM_surf_mesh_part_n_face_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return number of vertices
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Number of faces
 */

int
PDM_surf_mesh_part_n_vtx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);

/**
 * \brief Return extents for any face
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Extents
 */

const double *
PDM_surf_mesh_part_extents_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return face global number
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return     Face global number
 */

const PDM_g_num_t *
PDM_surf_mesh_part_face_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return Vertex global number
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Vertex global number
 */

const PDM_g_num_t *
PDM_surf_mesh_part_vtx_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return Edge global number
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return  Edge global number
 */

const PDM_g_num_t *
PDM_surf_mesh_part_edge_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return Face to edge connectivity
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Face to edge connectivity
 */

const int *
PDM_surf_mesh_part_face_edge_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return Face to vertex connectivity
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Face to vertex connectivity
 */

const int *
PDM_surf_mesh_part_face_vtx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return Face to vertex connectivity index
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Face to vertex connectivity index
 */

const int *
PDM_surf_mesh_part_face_vtx_idx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return Face to edge connectivity index
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Face to edge connectivity index
 */

const int *
PDM_surf_mesh_part_face_edge_idx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return vertex coordinates
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Vertex coordinates
 */

const double *
PDM_surf_mesh_part_vtx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
);


/**
 * \brief Return global number of edges
 *
 * \param [in]  mesh       Mesh object
 *
 * \return   Global number of edges
 */

PDM_g_num_t
PDM_surf_mesh_n_g_edge_get
(
 PDM_surf_mesh_t      *mesh
);


/**
 * \brief Return global number of vertices
 *
 * \param [in]  mesh       Mesh object
 *
 * \return   Global number of vertices
 */

PDM_g_num_t
PDM_surf_mesh_n_g_vtx_get
(
 PDM_surf_mesh_t      *mesh
);


/**
 * \brief Return global number of faces
 *
 * \param [in]  mesh       Mesh object
 *
 * \return   Global number of faces
 */

PDM_g_num_t
PDM_surf_mesh_n_g_face_get
(
 PDM_surf_mesh_t      *mesh
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_MESH_H__ */
