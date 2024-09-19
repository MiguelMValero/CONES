/*
 * \file
 */

#ifndef __PDM_BOX_H__
#define __PDM_BOX_H__

/*============================================================================
 * Handle boxes aligned with Cartesian axes.
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#include <stdbool.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_morton.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

typedef struct _PDM_boxes_t PDM_boxes_t;

/* Collection of boxes */

typedef struct _PDM_box_set_t PDM_box_set_t;

/* Distribution on octree or quadtree */

typedef struct _PDM_box_distrib_t PDM_box_distrib_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Normalize the coordinates of a point according to a box set
 *
 * \param [in]   boxes              Pointer to box set structure
 * \param [in]   pt_origin          Coordinates (size = 3)
 * \param [out]  pt_normalized      Normalized coordinates (size = 3)
 *
 */

void
PDM_box_set_normalize
(
 PDM_box_set_t *boxes,
 const double  *pt_origin,
 double        *pt_nomalized
 );


/**
 *
 * \brief De-normalize the coordinates of a point according to a box set
 *
 * \param [in]   boxes              Pointer to box set structure
 * \param [in]   pt_normalized      Normalized coordinates (size = 3)
 * \param [out]  pt_origin          Coordinates (size = 3)
 *
 */

void
PDM_box_set_normalize_inv
(
 PDM_box_set_t *boxes,
 const double  *pt_nomalized,
 double        *pt_origin
 );


/**
 *
 * \brief Normalize a set of coordinates according to a box set
 *
 * This implementation prevents division by zero.
 *
 * \param [in]   boxes              Pointer to box set structure
 * \param [in]   n_pts              Number of coordinates
 * \param [in]   pts_origin         Coordinates (size = 3 * \ref n_pts)
 * \param [out]  pts_normalized     Normalized coordinates (size = 3 * \ref n_pts)
 *
 */

void
PDM_box_set_normalize_robust
(
 PDM_box_set_t  *boxes,
 const int       n_pts,
 double         *pts_origin,
 double         *pts_normalized
 );

/**
 *
 * \brief Normalize a set of normal vectors according to a box set

 *
 * \param [in]   boxes              Pointer to box set structure
 * \param [in]   n_pts              Number of coordinates
 * \param [in]   pts_origin         Coordinates (size = 3 * \ref n_pts)
 * \param [out]  pts_normalized     Normalized coordinates (size = 3 * \ref n_pts)
 *
 */

void
PDM_box_set_normalize_normal_vector
(
 PDM_box_set_t  *boxes,
 const int       n_pts,
 double         *pts_origin,
 double         *pts_normalized
 );


/**
 *
 * \brief Remove duplicated boxes in a box set
 *
 * \param [in,out]  boxes     Pointer to box set structure
 *
 */

void
PDM_box_set_remove_duplicate
(
 PDM_box_set_t  *boxes
 );


/**
 * \brief Create a \ref PDM_boxes_t structure and initialize it
 *
 * \param [in] dim          Spatial dimension
 * \param [in] n_boxes      Number of elements to create
 * \param [in] box_gnum     Global ids of boxes
 * \param [in] box_extents  Coordinate extents (size = 2 * \ref n_boxes * \ref dim, as
 *                          xmin1, ymin1, .. xmax1, ymax1, ..., xmin2, ...)
 * \param [in] origin       Initial location (size = 3 * \ref n_boxes, as
 *                           iproc, i_part, local num, ...)
 *
 * \return   A new allocated pointer to a \ref PDM_boxes_t structure.
 */

PDM_boxes_t *
PDM_boxes_create(const int          dim,
                 int                n_boxes,
                 const PDM_g_num_t *box_gnum,
                 const double      *box_extents,
                 const int          n_part_orig,
                 const int         *n_boxes_orig,
                 const int         *origin);


/**
 * \brief Create a \ref PDM_box_set_t structure and initialize it
 *
 * \param[in] dim               Spatial dimension
 * \param[in] normalize         1 if boxes are to be normalized, 0 otherwize
 * \param[in] allow_projection  If 1, project to lower dimension if all boxes
 *                              are cut by the median plane of the set
 * \param[in] n_boxes           Number of elements to create
 * \param[in] box_gnum          Global ids of boxes
 * \param[in] extents           Coordinate extents (size = 2 * \ref n_boxes * \ref dim, as
 *                              xmin1, ymin1, .. xmax1, ymax1, ..., xmin2, ...)
 * \param[in] origin            Initial location (size = 3 * \ref n_boxes, as
 *                              iproc, i_part, local num, ...)
 * \param[in] comm              Associated MPI communicator
 *
 * \return   A new allocated pointer to a \ref PDM_box_set_t structure
 */

PDM_box_set_t *
PDM_box_set_create(int                dim,
                   int                normalize,
                   int                allow_projection,
                   int                n_boxes,
                   const PDM_g_num_t *box_gnum,
                   const double      *box_extents,
                   const int          n_part_orig,
                   const int         *n_boxes_orig,
                   const int         *origin,
                   PDM_MPI_Comm       comm);


/**
 * \brief Delete a \ref PDM_boxes_t structure
 *
 * \param [in,out]  boxes  Pointer to the \ref PDM_boxes_t structure to delete
 */

void
PDM_boxes_destroy(PDM_boxes_t  *boxes);


/**
 * \brief Destroy a \ref PDM_box_set_t structure
 *
 * \param [in,out] boxes Pointer to pointer to the \ref PDM_box_set_t structure to delete
 */

void
PDM_box_set_destroy(PDM_box_set_t  **boxes);


/**
 * \brief Return the dimension associated with a set of boxes.
 *
 * \param [in] boxes  Pointer to set of boxes
 *
 * \return   Associated spatial dimension
 */

int
PDM_box_set_get_dim(const PDM_box_set_t  *boxes);


/**
 * \brief Return the local number of boxes in a set
 *
 * \param [in] boxes  Pointer to set of boxes
 *
 * \return   Local number of boxes
 */

int
PDM_box_set_get_size(const PDM_box_set_t  *boxes);


/**
 * \brief Return the global number of boxes in a set.
 *
 * \param [in] boxes  Pointer to set of boxes
 *
 * \return  Global number of boxes
 */

PDM_g_num_t
PDM_box_set_get_global_size(const PDM_box_set_t  *boxes);


/**
 * \brief Return extents associated with a set of boxes.
 *
 * The extents array is organized in the following fashion:
 * {x_min_0, y_min_0, ..., x_max_0, y_max_0, ...
 *  x_min_n, y_min_n, ..., x_max_n, y_max_n, ...}
 *
 * Its size is thus: \ref n_boxes * \ref dim * 2.
 *
 * \param [in]   boxes  Pointer to set of boxes
 *
 * \return   Pointer to extents array
 */

const double *
PDM_box_set_get_extents(PDM_box_set_t  *boxes);


/**
 * \brief Return global ids associated with a set of boxes.
 *
 * \param [in]  boxes  Pointer to set of boxes
 *
 * \returns  Pointer to global box ids array
 */

const PDM_g_num_t *
PDM_box_set_get_g_num(PDM_box_set_t  *boxes);


/**
 * \brief Return global ids associated with a set of boxes (copied from another rank).
 *
 * \param [in]   boxes   Pointer to set of boxes
 * \param [in]   i_rank  Copied rank
 *
 * \return   Pointer to global box ids array
 */

PDM_g_num_t *
PDM_box_set_get_rank_boxes_g_num(PDM_box_set_t  *boxes,
                                 const int       i_rank);


/**
 * \brief Return initial location associated with a set of boxes.
 *
 * \param [in]   boxes  Pointer to set of boxes
 *
 * \return   Pointer to initial location array
 */

const int *
PDM_box_set_origin_get(PDM_box_set_t  *boxes);


/**
 * \brief Build a Morton_index to get a well-balanced distribution of the boxes.
 *
 * \param [in]     boxes       Pointer to associated \ref PDM_box_set_t structure
 * \param [in,out] distrib     Pointer to a \ref PDM_box_distrib_t structure
 * \param [in]     n_leaves    Number of leaves with weight > 0
 * \param [in]     leaf_codes  Morton code for each leaf
 * \param [in]     weight      Number of boxes related to each leaf
 */

void
PDM_box_set_build_morton_index(const PDM_box_set_t *boxes,
                               PDM_box_distrib_t   *distrib,
                               int                  n_leaves,
                               PDM_morton_code_t   *leaf_codes,
                               double              *weight);


/**
 * \brief Redistribute boxes over the ranks according to the Morton index to
 * assume a better balanced distribution of the boxes.
 *
 *  \param[in]     box_distrib Data structure on box distribution
 *  \param[in,out] box_set     Pointer to the structure to redistribute
 */

void
PDM_box_set_redistribute(const PDM_box_distrib_t  *box_distrib,
                         PDM_box_set_t            *boxes);


/**
 * \brief Dump a \ref PDM_box_set_t structure.
 *
 * \param [in] box_set   Pointer to the PDM_box_t structure
 * \param [in] verbosity Verbosity level (0 or 1)
 */

void
PDM_box_set_dump(const PDM_box_set_t  *boxes,
                 int                   verbosity);


/**
 * \brief Receive data from origin for any box
 *
 * \param[in]     box_set                Pointer to the PDM_box_t structure
 * \param[in]     t_stride               Type of stride
 * \param[in]     stride_cst             Constant stride
 * \param[in]     data_size              Size of data
 * \param[in]     origin_distrib_stride  Origin stride distribution
 * \param[in]     origin_distrib_data    Origin data distribution
 * \param[in,out] current_distrib_stride Current stride distribution (Allocate and compute if input is NULL,
 *                                       otherwise nothing)
 * \param[out]    current_distrib_data   Current data distribution
 *
 * \return  Size of current_distrib_data
 */

void
PDM_box_set_recv_data_from_origin_distrib
(
 PDM_box_set_t  *boxes,
 PDM_stride_t    t_stride,
 int             stride_cst,
 size_t          data_size,
 int           **origin_distrib_stride,
 void          **origin_distrib_data,
 int           **current_distrib_stride,
 void          **current_distrib_data
 );


/**
 * \brief Send data to origin for any box
 *
 * \param [in]  box_set                 pointer to the \ref PDM_box_t structure
 * \param [in]  t_stride                Type of stride
 * \param [in]  stride_cst              Constant stride
 * \param [in]  data_size               Size of data
 * \param [in]  current_distrib_stride  Current stride distribution
 * \param [in]  current_distrib_data    Current data distribution
 * \param [out] origin_distrib_stride   Origin stride distribution
 * \param [out] origin_distrib_data     Origin data distribution
 *
 * \return   Size of origin_distrib_data
 */

void
PDM_box_set_send_data_to_origin_distrib
(
 PDM_box_set_t *boxes,
 PDM_stride_t   t_stride,
 int            stride_cst,
 size_t         data_size,
 int           *current_distrib_stride,
 void          *current_distrib_data,
 int          **origin_distrib_stride,
 void         **origin_distrib_data
);


/**
 * \brief Send copies of boxes from selected ranks to all other ranks for better load balancing
 *
 * \param [in] boxes            Pointer to the PDM_box_t structure
 * \param [in] n_copied_ranks   Number of copied ranks
 * \param [in] copied_ranks     List of copied ranks
 */

void
PDM_box_copy_boxes_to_ranks
(
 PDM_box_set_t  *boxes,
 const int       n_copied_ranks,
 int            *copied_ranks
);

/**
 * \brief Setup a shared structure among nodes
 *
 * \param [in] boxes            Pointer to the PDM_box_t structure
 */
void
PDM_box_copy_boxes_to_shm
(
 PDM_box_set_t  *boxes
);

void
PDM_box_set_free_copies
(
 PDM_box_set_t  **boxes
);


/**
 * \brief Create a \ref PDM_box_distrib_t structure.
 *
 * \param [in] n_boxes    Number of boxes
 * \param [in] n_g_boxes  Global number of boxes
 * \param [in] max_level  Max level reached locally in the related tree
 * \param [in] comm       MPI communicator. on which the distribution takes place
 *
 * \return  A pointer to a new allocated \ref PDM_box_distrib_t structure.
 */

PDM_box_distrib_t *
PDM_box_distrib_create(int          n_boxes,
                       PDM_g_num_t  n_g_boxes,
                       int          max_level,
                       PDM_MPI_Comm comm);


PDM_box_distrib_t *
PDM_box_distrib_shared_create(int          n_boxes,
                              PDM_g_num_t  n_g_boxes,
                              int          gmax_level,
                              PDM_MPI_Comm comm);

/**
 * \brief Destroy a \ref PDM_box_distrib_t structure.
 *
 * \param [in,out]  distrib  Pointer to pointer to the structure to destroy
 */

void
PDM_box_distrib_destroy(PDM_box_distrib_t  **distrib);


/**
 * \brief Delete redundancies in box distribution
 *
 * \param [in,out]  distrib  Pointer to the \ref PDM_box_distrib_t structure
 */

void
PDM_box_distrib_clean(PDM_box_distrib_t  *distrib);


/**
 * \brief Display a histogramm on leaves associated to the boxes and
 * several other pieces of information (min, max, ...)
 *
 * \param [in] distrib  Pointer to the \ref PDM_box_distrib_t structure
 * \param [in] comm     Associated MPI communicator
 */

void
PDM_box_distrib_dump_statistics(const PDM_box_distrib_t *distrib,
                                PDM_MPI_Comm             comm);


/**
 * \brief Dump a \ref PDM_box_distrib_t structure
 *
 * \param [in] distrib  Pointer to the \ref PDM_box_distrib_t structure
 */

void
PDM_box_distrib_dump(const PDM_box_distrib_t *distrib);


void
PDM_box_set_normalization_get
(
 PDM_box_set_t  *boxes,
 double        **s,
 double        **d
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_BOX_H__ */
