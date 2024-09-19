#ifndef PDM_PARA_OCTREE_CUH
#define	PDM_PARA_OCTREE_CUH

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_para_octree.h"

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

/*  */

/*============================================================================
 * Global variable
 *============================================================================*/

 //extern MANAGED PDM_Handles_t *_octrees;

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
 * \return     Identifier
 */

int
PDM_para_octree_create_GPU
(
 const int n_point_cloud,
 const int depth_max,
 const int points_in_leaf_max,
 const int build_leaf_neighbours,
 const PDM_MPI_Comm comm
);

/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_free_GPU
(
 const int          id
);

/**
 *
 * Look for closest points stored inside an octree
 *
 * \param [in]   id                     Identifier
 * \param [in]   n_closest_points       Number of closest points to find
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [in]   pts_g_num              Point global numbers
 * \param [out]  closest_octree_pt_id   Closest points in octree global number
 * \param [out]  closest_octree_pt_dist Closest points in octree distance
 *
 */

void
PDM_para_octree_closest_point_GPU
(
 PDM_octree_t *octree,
 const int    n_closest_points,
 const int    n_pts,
 double      *pts,
 PDM_g_num_t *pts_g_num,
 PDM_g_num_t *closest_octree_pt_g_num,
 double      *closest_octree_pt_dist2
 );

//test
__global__
void
print_from_gpu
(
PDM_octree_t *octree
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PARA_OCTREE_CUH */
