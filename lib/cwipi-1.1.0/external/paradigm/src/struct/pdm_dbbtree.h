/*
 * \file
 */

#ifndef __PDM_DBBTREE_H__
#define __PDM_DBBTREE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_box.h"
#include "pdm_part_to_block.h"

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
 * Type
 *============================================================================*/

/**
 * \struct PDM_dbbtree_t
 * \brief  Distributed boundary box tree
 *
 *  PDM_dbbtree_t defines a distributed boundary box tree
 *
 */

typedef struct _PDM_dbbtree_t PDM_dbbtree_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_dbbtree_t structure
 *
 * This function returns an initialized \ref PDM_dbbtree_t structure
 *
 * \param [in]  comm             Associated communicator
 * \param [in]  dim              boxes dimension
 * \param [in]  global_extents   Globals of elements to storage into the tree
 *                               (automatic computation if NULL)
 *
 * \return      A new initialized \ref PDM_dbbtree_t structure
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_create
(
 PDM_MPI_Comm  comm,
 int           dim,
 double       *global_extents
 );



/**
 * \brief Free a \ref PDM_dbbtree_t structure
 *
 * \param [in]  dbbt   Pointer to a distributed bounding box tree
 *
 * \return      NULL
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_free
(
 PDM_dbbtree_t     *dbbt
 );


/**
 * \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  n_part    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according to
 * the tree location
 *
 */

PDM_box_set_t *
PDM_dbbtree_boxes_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum
);


PDM_box_set_t *
PDM_dbbtree_boxes_set_with_init_location
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const int         **init_location,
 const double      **extents,
 const PDM_g_num_t **gNum
);

PDM_box_set_t *
PDM_dbbtree_boxes_set_for_intersect_line
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum,
 const int           n_line,
 double             *line_coord
);

/**
 * \brief Assign boxes to intersect to the tree.
 *
 * This function  assigns boxes to intersect to the tree.
 *
 * \param [in]  dbbt       Pointer to a distributed bounding box tree
 * \param [in]  n_part     Number of partitions
 * \param [in]  nElts      Number of elements of each partition
 * \param [in]  extents    Extents of each element of each partition
 * \param [in]  gNum       Global number of each element of each partition
 * \param [out] box_index  Pointer to the index array on associated tree bounding boxeq
 * \param [out] box_g_num  Pointer to the list of intersecting bounding boxes
 *
 * \return associated \ref PDM_box_set_t structure distributed according
 * to the tree intersection
 *
 */
PDM_box_set_t *
PDM_dbbtree_intersect_boxes_with_init_location_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const int         **init_location,
 const double      **extents,
 const PDM_g_num_t **gNum,
 int                *box_index[],
 int                *box_l_num[]
);


/**
 * \brief Assign boxes to intersect to the tree.
 *
 * This function  assigns boxes to intersect to the tree.
 *
 * \param [in]  dbbt           Pointer to a distributed bounding box tree
 * \param [in]  n_part         Number of partitions
 * \param [in]  nElts          Number of elements of each partition
 * \param [in]  init_location  Init location of each element of each partition (triplet rank/ipart/ielt)
 * \param [in]  extents        Extents of each element of each partition
 * \param [in]  gNum           Global number of each element of each partition
 * \param [out] box_index      Pointer to the index array on associated tree bounding boxeq
 * \param [out] box_g_num      Pointer to the list of intersecting bounding boxes
 *
 * \return associated \ref PDM_box_set_t structure distributed according
 * to the tree intersection
 *
 */

PDM_box_set_t *
PDM_dbbtree_intersect_boxes_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum,
 int                *box_index[],
 int                *box_l_num[]
 );


/**
 *
 * \brief Get the boxes closer than the upper bound distance
 *
 *   \param [in]  bt                 Pointer to box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts                Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  pts_g_num          Point global ids
 *   \param [in]  upper_bound_dist2  Upper bound of the square of the distance (size = \ref n_pts)
 *   \param [out] box_index          Index of boxes (size = \ref n_pts + 1)
 *   \param [out] box_g_num          Global ids of boxes (size = \ref i_boxes[\ref n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_get
(
 PDM_dbbtree_t   *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *box_index[],
 PDM_g_num_t     *box_g_num[]
);

/**
 *
 * \brief Get the boxes closer than the upper bound distance
 *
 *   \param [in]  bt                 Pointer to box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts                Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  pts_g_num          Point global ids
 *   \param [in]  upper_bound_dist2  Upper bound of the square of the distance (size = \ref n_pts)
 *   \param [out] n_extract_boxes    Number of extracted box
 *   \param [out] box_l_num          Index of boxes (size = \ref n_extract_boxes )
 *   \param [out] box_pts_idx        Index of boxes (size = \ref n_extract_boxes + 1)
 *   \param [out] box_g_num          Global ids of boxes (size = \ref box_pts_idx[\ref n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_pts_shared_get
(
 PDM_dbbtree_t   *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *out_n_extract_boxes,
 PDM_g_num_t     *out_box_gnum[],
 int             *out_box_init_location[],
 int             *out_dbox_pts_idx[],
 PDM_g_num_t     *out_dbox_pts_g_num[],
 double          *out_dbox_pts_coord[]
);

/**
 *
 * \brief Get the boxes closer than the upper bound distance (Asynchronous)
 *
 *   \param [in]  bt                 Pointer to box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts                Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  pts_g_num          Point global ids
 *   \param [in]  upper_bound_dist2  Upper bound of the square of the distance (size = \ref n_pts)
 *   \param [out] box_index          Index of boxes (size = \ref n_pts + 1)
 *   \param [out] box_g_num          Global ids of boxes (size = \ref i_boxes[\ref n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_get_async
(
 PDM_dbbtree_t   *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *box_index[],
 PDM_g_num_t     *box_g_num[]
);


/**
 *
 * \brief Get an indexed list of all points inside the boxes of a distributed box tree
 *
 *   \param [in]  dbbt               Pointer to distributed box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts_g_num          Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord          Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  n_boxes            Number of boxes
 *   \param [in]  box_g_num          Global ids of boxes (size = \ref n_boxes)
 *   \param [out] pts_in_box_idx     Index of points in boxes (size = \ref n_boxes + 1, allocated inside function)
 *   \param [out] pts_in_box_g_num   Global ids of points in boxes (size = \ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [out] pts_in_box_coord   Coordinates of points in boxes (size = 3 *\ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [in]  ellipsoids         Consider boxes as axis-aligned ellipsoids (1 or 0)
 *
 */

void
PDM_dbbtree_points_inside_boxes
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 const int           n_boxes,
 const PDM_g_num_t   box_g_num[],
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord,
 const int           ellipsoids
 );

void
PDM_dbbtree_points_inside_boxes_block_frame
(
 PDM_dbbtree_t        *dbbt,
 const int             n_pts,
 PDM_g_num_t           pts_g_num[],
 double                pts_coord[],
 PDM_part_to_block_t **ptb_out,
 int                 **dbox_pts_n,
 PDM_g_num_t         **dbox_pts_g_num,
 double              **dbox_pts_coord,
 const int             ellipsoids
 );


/**
 *
 * \brief Get an indexed list of all points inside the boxes of a distributed box tree
 *
 *   \param [in]  dbbt               Pointer to distributed box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts_g_num          Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord          Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  n_boxes            Number of boxes
 *   \param [in]  box_g_num          Global ids of boxes (size = \ref n_boxes)
 *   \param [out] pts_in_box_idx     Index of points in boxes (size = \ref n_boxes + 1, allocated inside function)
 *   \param [out] pts_in_box_g_num   Global ids of points in boxes (size = \ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [out] pts_in_box_coord   Coordinates of points in boxes (size = 3 *\ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [in]  ellipsoids         Consider boxes as axis-aligned ellipsoids (1 or 0)
 *
 */
void
PDM_dbbtree_points_inside_boxes_shared
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 const int           n_boxes,
 const PDM_g_num_t   box_g_num[],
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord,
 const int           ellipsoids
);

/**
 *
 * \brief Get an indexed list of all boxes containing points
 *
 *   \param [in]  dbbt        Pointer to distributed box tree structure
 *   \param [in]  n_pts       Number of points
 *   \param [in]  pts_g_num   Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord   Point coordinates (size = 3 * \ref n_pts)
 *   \param [out] box_idx     Index of boxes (size = \ref n_pts + 1, allocated inside function)
 *   \param [out] box_g_num   Global ids of boxes (size = \ref box_idx[\ref n_pts], allocated inside function)
 *   \param [in]  ellipsoids  Consider boxes as axis-aligned ellipsoids (1 or 0)
 *
 */

void
PDM_dbbtree_boxes_containing_points
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 int               **box_idx,
 PDM_g_num_t       **box_g_num,
 const int           ellipsoids
 );


/**
 *
 * \brief Get an indexed list of all boxes intersecting lines
 *
 *   \param [in]  dbbt        Pointer to distributed box tree structure
 *   \param [in]  n_line      Number of points
 *   \param [in]  line_g_num  Line global ids (size = \ref n_line)
 *   \param [in]  line_coord  Line coordinates (size = 6 * \ref n_line)
 *   \param [out] box_idx     Index of boxes (size = \ref n_line + 1, allocated inside function)
 *   \param [out] box_g_num   Global ids of boxes (size = \ref box_idx[\ref n_line], allocated inside function)
 *
 */

void
PDM_dbbtree_lines_intersect_boxes
(
 PDM_dbbtree_t  *dbbt,
 const int       n_line,
 PDM_g_num_t    *line_g_num,
 double         *line_coord,
 int           **box_index,
 PDM_g_num_t   **box_g_num
);


/**
 *
 * \brief Export boxes and points to ASCII VTK format for debugging purposes
 *
 * A set of files (normalized/unnormalized points and boxes) is written
 * for each MPI rank.
 *
 *   \param [in]  filename_pattern  File name pattern (prefix)
 *   \param [in]  dbbt              Pointer to distributed box tree structure
 *   \param [in]  n_pts             Number of points
 *   \param [in]  pts_g_num         Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord         Point coordinates (size = 3 * \ref n_pts)
 *
 */

void
PDM_dbbtree_points_debug
(
 char*               filename_pattern,
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 double             *pts_coord,
 PDM_g_num_t        *pts_g_num
);


/**
 *
 * \brief Get an indexed list of all lines intersecting boxes
 * /!\ Results conform to the dbbtree's box partitionning
 *
 *   \param [in]  dbbt             Pointer to distributed box tree structure
 *   \param [in]  n_line           Number of points
 *   \param [in]  line_g_num       Line global ids (size = \ref n_line)
 *   \param [in]  line_coord       Line coordinates (size = 6 * \ref n_line)
 *   \param [out] box_line_idx     Index of lines (size = \ref n_box + 1, allocated inside function)
 *   \param [out] box_line_g_num   Global ids of lines (size = \ref box_idx[\ref n_box], allocated inside function)
 *
 */

void
PDM_dbbtree_lines_intersect_boxes2
(
 PDM_dbbtree_t  *dbbt,
 const int       n_line,
 PDM_g_num_t    *line_g_num,
 double         *line_coord,
 int            *n_part,
 int           **redistrib_n_box,
 PDM_g_num_t  ***redistrib_box_ln_to_gn,
 int          ***redistrib_box_init_location,
 int          ***box_line_idx,
 PDM_g_num_t  ***box_line_g_num
 );

void
PDM_dbbtree_lines_intersect_boxes2_shared
(
 PDM_dbbtree_t   *dbbt,
 const int        n_line,
 PDM_g_num_t     *line_g_num,
 double          *line_coord,
 int             *n_part,
 int            **redistrib_n_box,
 PDM_g_num_t   ***redistrib_box_ln_to_gn,
 int           ***redistrib_box_init_location,
 int           ***box_line_idx,
 PDM_g_num_t   ***box_line_g_num
);

/**
 *
 * \brief Get an indexed list of all boxes intersecting volumes
 *
 * \param [in]   dbbt                  Pointer to distributed box tree structure
 * \param [in]   n_volumes             Number of volumes
 * \param [in]   volume_g_num          Global number of volumes
 * \param [in]   volume_plane_idx      Index of the number of planes per volume
 * \param [in]   plane_normal          Oriented normal vector for a given plane (oriented toward the interior of the volume)
 * \param [in]   plane_pt_coord        Point on plane coordinates (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [out]  volume_box_idx        Index of boxes (size = \ref n_line + 1, allocated inside function)
 * \param [out]  volume_box_g_num      Global ids of boxes (size = \ref box_line_idx[\ref n_line], allocated inside function)
 *
 */

void
PDM_dbbtree_volumes_intersect_boxes
(
 PDM_dbbtree_t  *dbbt,
 const int       n_volumes,
 PDM_g_num_t    *volume_g_num,
 int            *volume_plane_idx,
 double         *plane_normal,
 double         *plane_pt_coord,
 int           **volume_box_idx,
 PDM_g_num_t   **out_volume_box_g_num
);


void
PDM_dbbtree_box_tree_write_vtk
(
 const char    *filename,
 PDM_dbbtree_t *dbbt,
 const int      i_copied_rank,
 const int      normalized
 );


void
PDM_dbbtree_box_tree_write_vtk2
(
 const char    *filename,
 PDM_dbbtree_t *dbbt,
 const int      i_copied_rank
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DBBTREE_H__ */
