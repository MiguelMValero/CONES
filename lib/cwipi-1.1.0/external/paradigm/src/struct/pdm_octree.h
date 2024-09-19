/*
 * \file
 * \author equemera
 *
 * \date November 8, 2017, 11:27 AM
 */

#ifndef PDM_OCTREE_H
#define	PDM_OCTREE_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_octree_seq.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_octree_t PDM_octree_t;

/**
 * \enum PDM_octree_child_t
 * \brief Names of 8 children of a node
 *
 */

typedef enum {
  PDM_NADIR,
  PDM_ZENITH,
  PDM_WEST,
  PDM_EAST,
  PDM_NORTH,
  PDM_SOUTH,
} PDM_octree_direction_t;

/**
 * \enum PDM_octree_child_t
 * \brief Names of 8 children of a node
 *
 */

typedef enum {
  PDM_NORTH_WEST_NADIR,
  PDM_NORTH_WEST_ZENITH,
  PDM_NORTH_EAST_NADIR,
  PDM_NORTH_EAST_ZENITH,
  PDM_SOUTH_WEST_NADIR,
  PDM_SOUTH_WEST_ZENITH,
  PDM_SOUTH_EAST_NADIR,
  PDM_SOUTH_EAST_ZENITH,
} PDM_octree_child_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create an octree structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 * \param [in]   comm               MPI communicator
 *
 * \return     Pointer to \ref PDM_octree object
 */

PDM_octree_t *
PDM_octree_create
(
 const int          n_point_cloud,
 const int          depth_max,
 const int          points_in_leaf_max,
 const double       tolerance,
 const PDM_MPI_Comm comm
);


/**
 *
 * \brief Create an octree structure from a sequential octree
 *
 * \param [in]   octree_seq         Pointer to sequential octree
 * \param [in]   comm               MPI communicator
 *
 * \return     Pointer to \ref PDM_octree object
 */

PDM_octree_t *
PDM_octree_from_octree_seq_create
(
 PDM_octree_seq_t   *octree_seq,
 const PDM_MPI_Comm  comm
);


/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   octree     Pointer to \ref PDM_octree object
 *
 */

void
PDM_octree_free
(
 PDM_octree_t *octree
);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [in]   i_point_cloud      Number of point cloud
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 * \param [in]   g_num              Point global number or NULL
 *
 */

void
PDM_octree_point_cloud_set
(
 PDM_octree_t      *octree,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords,
 const PDM_g_num_t *g_num
);


/**
 *
 * \brief Build octree
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 *
 */

void
PDM_octree_build
(
 PDM_octree_t      *octree
);


/**
 *
 * \brief Used processes extents
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [out]  used_ranks         Used ranks
 * \param [out]  extents            Used ranks extents
 *
 * \return Number of used ranks
 */

int
PDM_octree_processes_extents_get
(
 PDM_octree_t     *octree,
 int              *used_ranks[],
 double           *extents[]
);


/**
 *
 * \brief Get root node id
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 *
 * \return     Root node identifier (-1 if octree is not built)
 *
 */

int
PDM_octree_root_node_id_get
(
 PDM_octree_t      *octree
);


/**
 *
 * \brief Get ancestor node id
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [in]   node_id            Node identifier
 *
 * \return     Ancestor node identifier
 *
 */

int
PDM_octree_ancestor_node_id_get
(
 PDM_octree_t      *octree,
 const int          node_id
);


/**
 *
 * \brief Get node extents
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [in]   node_id            Node identifier
 *
 * \return     Extents
 *
 */

const double *
PDM_octree_node_extents_get
(
 PDM_octree_t      *octree,
 const int          node_id
);


/**
 *
 * \brief Get children of a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [in]   node_id            Node identifier
 * \param [in]   child              Children
 *
 * \return     Children node id
 *
 */

int
PDM_octree_children_get
(
 PDM_octree_t             *octree,
 const int                 node_id,
 const PDM_octree_child_t  child
);


/**
 *
 * \brief Get Neighbor of node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [in]   node_id            Node identifier
 * \param [in]   direction          Neighbor direction
 *
 * \return     Neighbor node id (-1 if no neighbor)
 *
 */

int
PDM_octree_neighbor_get
(
 PDM_octree_t                 *octree,
 const int                     node_id,
 const PDM_octree_direction_t  direction
);

/**
 *
 * \brief Get the number of point inside a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [in]   node_id            Node identifier
 *
 * \return   Number of points
 *
 */

int
PDM_octree_n_points_get
(
 PDM_octree_t            *octree,
 const int                node_id
);


/**
 *
 * \brief Get indexes of points inside a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [in]   node_id            Node identifier
 * \param [out]  point_clouds_id    Point clouds number
 *                                  (size = Number of points inside the node)
 * \param [out]  point_indexes      Point indexes
 *                                  (size = Number of points inside the node)
 *
 */

void
PDM_octree_points_get
(
 PDM_octree_t            *octree,
 const int                node_id,
 int                    **point_clouds_id,
 int                    **point_indexes
);


/**
 *
 * \brief Is it a leaf
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 * \param [in]   node_id            Node identifier
 *
 * \return   1 or 0
 *
 */

int
PDM_octree_leaf_is
(
 PDM_octree_t            *octree,
 const int                node_id
);


/**
 *
 * \brief Get extents
 *
 * \param [in]   octree             Pointer to \ref PDM_octree object
 *
 * \return     Extents
 *
 */

double *
PDM_octree_extents_get
(
 PDM_octree_t            *octree
);


/**
 *
 * Look for closest points stored inside an octree
 *
 *
 * \param [in]   octree                 Pointer to \ref PDM_octree object
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [in]   pts_g_num              Point global numbers
 * \param [out]  closest_octree_pt_id   Closest point in octree global number
 * \param [out]  closest_octree_pt_dist Closest point in octree distance
 *
 */

void
PDM_octree_closest_point
(
PDM_octree_t     *octree,
const int         n_pts,
double           *pts,
PDM_g_num_t      *pts_g_num,
PDM_g_num_t      *closest_octree_pt_g_num,
double           *closest_octree_pt_dist2
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_OCTREE_H */

