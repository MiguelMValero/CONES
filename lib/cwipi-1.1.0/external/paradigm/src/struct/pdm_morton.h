/*
 * \file
 */

#ifndef __PDM_MORTON_H__
#define __PDM_MORTON_H__

/*============================================================================
 * Morton encoding for 2D or 3D coordinates.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA
  Copyright (C) 2008-2010  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

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

typedef enum {

  PDM_MORTON_EQUAL_ID,
  PDM_MORTON_SAME_ANCHOR,
  PDM_MORTON_DIFFERENT_ID

} PDM_morton_compare_t;

typedef unsigned int   PDM_morton_int_t;

typedef struct {

  PDM_morton_int_t   L;     /* Level in the tree structure */
  PDM_morton_int_t   X[3];  /* X, Y, Z coordinates in Cartesian grid */

} PDM_morton_code_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

static const PDM_morton_int_t PDM_morton_max_level = 31u;

/**
 * \brief Determine the global extents associated with a set of coordinates
 *
 * \param [in]   dim        Spatial dimension
 * \param [in]   n_coords   Local number of coordinates
 * \param [in]   coords     Coordinates, interleaved (size = \ref dim * \ref n_coords)
 * \param [out]  g_extents  Global extents (size = 2 * \ref dim)
 * \param [in]   comm       Associated MPI communicator
 *
 */

void
PDM_morton_get_coord_extents(int            dim,
                             size_t         n_coords,
                             const double   coords[],
                             double         g_extents[],
                             PDM_MPI_Comm   comm);


/**
 * \brief Determine the global extents associated with a set of local extents
 *
 * \param [in]   dim        Spatial dimension
 * \param [in]   n_extents  Local number of extents
 * \param [in]   extents    Local extents (size = 2 * \ref dim * \ref n_extents)
 * \param [out]  g_extents  Global extents (size = 2 * \ref dim)
 * \param [in]   comm       Associated MPI communicator
 *
 */

void
PDM_morton_get_global_extents(int               dim,
                              size_t            n_extents,
                              const double      extents[],
                              double            g_extents[],
                              PDM_MPI_Comm      comm);


/**
 * \brief Build a Morton code according to the level in an octree grid and its coordinates in the grid.
 *
 * \param [in]   dim        Spatial dimension (1, 2 or 3)
 * \param [in]   level      Level in the grid
 * \param [in]   coords     Coordinates in the grid (normalized) (size = \ref dim)
 *
 * \return                  a Morton code
 *
 */

PDM_morton_code_t
PDM_morton_encode(int                dim,
                  PDM_morton_int_t   level,
                  const double       coords[]);


/**
 * \brief Encode an array of coordinates.
 *
 * The caller is responsible for freeing the returned array once it is
 * no longer useful.
 *
 * \param [in]   dim        Spatial dimension (1, 2 or 3)
 * \param [in]   level      Level in the grid
 * \param [in]   extents    Coordinate extents for normalization (size = 2 * \ref dim)
 * \param [in]   n_coords   Local number of coordinates
 * \param [in]   coords     Coordinates, interleaved (size = \ref dim * \ref n_coords)
 * \param [out]  m_code     Array of corresponding Morton codes (size = \ref n_coords)
 * \param [out]  d          Normalization vector (dilatation component)
 * \param [out]  s          Normalization vector (translation component)
 *
 */

void
PDM_morton_encode_coords(int                dim,
                         PDM_morton_int_t   level,
                         const double       extents[],
                         size_t             n_coords,
                         const double       coords[],
                         PDM_morton_code_t  m_code[],
                         double             d[3],
                         double             s[3]);


/**
 * \brief Compute the Morton codes of the children of a node
 *
 * Given a Morton code in the grid, compute the Morton codes of its
 * children when refining the grid by one level.
 *
 * \param [in]   dim        Spatial dimension (1, 2 or 3)
 * \param [in]   parent     Morton code associated with the parent
 * \param [out]  children   Array of children Morton codes (size = 2^\ref dim)
 *
 */

void
PDM_morton_get_children(int                dim,
                        PDM_morton_code_t  parent,
                        PDM_morton_code_t  children[]);


/**
 * \brief Compare two Morton codes
 *
 * Compare two Morton encoding and check if these two codes are equal,
 * different or shared the same anchor.
 *
 * \param [in]   dim        Spatial dimension (2 or 3)
 * \param [in]   code_a     First Morton code to compare
 * \param [in]   code_b     Second Morton code to compare
 *
 * \return                  A type on the kind of relation between the two Morton codes.
 *
 */

PDM_morton_compare_t
PDM_morton_compare(int                dim,
                   PDM_morton_code_t  code_a,
                   PDM_morton_code_t  code_b);


/**
 * \brief Get local order in a list of Morton codes
 *
 * \param [in]   n_codes       Number of Morton codes to order
 * \param [in]   morton_codes  Array of Morton codes to order
 * \param [out]  order         Pointer to pre-allocated ordering table
 *
 */

void
PDM_morton_local_order(int                      n_codes,
                       const PDM_morton_code_t  morton_codes[],
                       int                      order[]);


/**
 * \brief Sort a local list of Morton codes
 *
 * \param [in]     n_codes       Number of Morton codes to sort
 * \param [in,out] morton_codes  Array of Morton codes to sort
 *
 */

void
PDM_morton_local_sort(int                n_codes,
                      PDM_morton_code_t  morton_codes[]);



#ifndef __cplusplus
/**
 * \brief Test if Morton code "a" is greater than Morton code "b"
 *
 * \param [in]   a     First Morton code to compare
 * \param [in]   b     Second Morton code to compare
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_a_gt_b(PDM_morton_code_t  a,
                  PDM_morton_code_t  b);

/**
 * \brief Test if Morton code "a" is greater than Morton code "b" (compare \em anchors)
 *
 * \param [in]   a     First Morton code to compare
 * \param [in]   b     Second Morton code to compare
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_a_gtmin_b(PDM_morton_code_t  a,
                     PDM_morton_code_t  b);
#endif


/**
 * \brief Copy Morton code "a" into Morton code "b"
 *
 * \param [in]   a     Morton code to copy
 * \param [out]  b     Morton code receiving the copy
 *
 */

void
PDM_morton_copy (PDM_morton_code_t  a,
                 PDM_morton_code_t *b);


/**
 * \brief Get the nearest common ancestor of two Morton codes
 *
 * \param [in]   code_a     First Morton code
 * \param [in]   code_b     Second Morton code
 * \param [out]  c          Nearest common ancestor of the two codes
 *
 */

void
PDM_morton_nearest_common_ancestor (PDM_morton_code_t  code_a,
                                    PDM_morton_code_t  code_b,
                                    PDM_morton_code_t *c);



#ifndef __cplusplus
/**
 * \brief Test if Morton code 'a' is an ancestor of Morton code 'b'
 *
 * \param [in]   a     First Morton code
 * \param [in]   b     Second Morton code
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_ancestor_is (PDM_morton_code_t a,
                        PDM_morton_code_t b);


/**
 * \brief Test if Morton code "a" is greater than or equal to Morton code "b"
 *
 * \param [in]   a     First Morton code to compare
 * \param [in]   b     Second Morton code to compare
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_a_ge_b(PDM_morton_code_t a,
                  PDM_morton_code_t b);


/**
 * \brief Test if Morton code "a" is equal to Morton code "b" at the level = max (level a , level b)
 *
 * \param [in]   a     First Morton code to compare
 * \param [in]   b     Second Morton code to compare
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_a_eq_b(PDM_morton_code_t  a,
                  PDM_morton_code_t  b);
#endif


/**
 * \brief Assign a level to a Morton code
 *
 * \param [in,out]  a     Morton code
 * \param [in]      l     Level to assign
 *
 */

void
PDM_morton_assign_level (PDM_morton_code_t  *a,
                         int                l);


/**
 * \brief Get the id associated to a Morton code in an array using a binary search.
 *
 * \param [in]  size   Size of the array
 * \param [in]  code   Morton code we are searching for
 * \param [in]  codes  Array of Morton codes
 *
 * \return             Id associated to the given code in the codes array.
 *
 */

int
PDM_morton_binary_search(int                 size,
                         PDM_morton_code_t   code,
                         PDM_morton_code_t  *codes);


/**
 * \brief Get the quantile associated to a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the quantiles.
 *
 * \param [in]  n_quantiles     Number of quantiles
 * \param [in]  code            Morton code we are searching for
 * \param [in]  quantile_start  First Morton code in each quantile (size = \ref n_quantiles)
 *
 * \return                      Quantile associated to the given code in the codes array.
 *
 */

size_t
PDM_morton_quantile_search(size_t              n_quantiles,
                           PDM_morton_code_t   code,
                           PDM_morton_code_t  *quantile_start);


/**
 * \brief Get the quantiles intersected by a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the quantiles.
 *
 * \param [in]  n_quantiles     Number of quantiles
 * \param [in]  code            Morton code we are searching for
 * \param [in]  quantile_start  First Morton code in each quantile (size = \ref n_quantiles)
 * \param [out] n_intersect     Number of intersected quantiles
 * \param [out] intersect       Intersected quantiles (size = \ref n_quantiles)
 *
 */

void
PDM_morton_quantile_intersect(size_t              n_quantiles,
                              PDM_morton_code_t   code,
                              PDM_morton_code_t  *quantile_start,
                              size_t             *start,
                              size_t             *end);


/**
 * \brief Get the Morton codes intersected by a given Morton code using a binary search.
 *
 * The array of Morton codes \emph must be sorted in ascending order.
 *
 * \param [in]  n_codes       Number of codes
 * \param [in]  code          Morton code we are searching for
 * \param [in]  codes         First Morton code in each code (size = \ref n_codes)
 * \param [out] start         Id of the first intersected Morton code
 * \param [out] end           Id of the last intersected Morton code + 1
 *
 */

void
PDM_morton_list_intersect(size_t              n_codes,
                          PDM_morton_code_t   code,
                          PDM_morton_code_t  *codes,
                          size_t             *start,
                          size_t             *end);


/**
 * \brief Build a global Morton encoding rank index from sorted Morton codes
 *
 * The rank_index[i] contains the first Morton code assigned to rank i.
 *
 * \param [in]  dim           Spatial dimension (1, 2 or 3)
 * \param [in]  gmax_level    Level in octree used to build the Morton encoding
 * \param [in]  n_codes       Number of Morton codes to be indexed
 * \param [in]  ordered_code  Array of Morton codes to be indexed (size = \ref n_codes)
 * \param [in]  weight        Weight associated to each Morton code (size = \ref n_codes)
 * \param [out] rank_index    Pointer to the global Morton encoding rank index (size = n_rank + 1)
 * \param [in]  comm          MPI communicator on which we build the global index
 *
 */

void
PDM_morton_ordered_build_rank_index
(
 int                      dim,
 int                      gmax_level,
 PDM_l_num_t              n_codes,
 const PDM_morton_code_t  ordered_code[],
 const int                weight[],
 PDM_morton_code_t        rank_index[],
 PDM_MPI_Comm             comm
 );


/**
 * \brief Build a global Morton encoding rank index
 *
 * The rank_index[i] contains the first Morton code assigned to rank i.
 *
 * \param [in]  dim           Spatial dimension (1, 2 or 3)
 * \param [in]  gmax_level    Level in octree used to build the Morton encoding
 * \param [in]  n_codes       Number of Morton codes to be indexed
 * \param [in]  ordered_code  Array of Morton codes to be indexed (size = \ref n_codes)
 * \param [in]  weight        Weight associated to each Morton code (size = \ref n_codes)
 * \param [in]  order         Ordering array (size = \ref n_codes)
 * \param [out] rank_index    Pointer to the global Morton encoding rank index (size = n_rank + 1)
 * \param [in]  comm          MPI communicator on which we build the global index
 *
 */

double
PDM_morton_build_rank_index(int                      dim,
                            int                      gmax_level,
                            PDM_l_num_t              n_codes,
                            const PDM_morton_code_t  code[],
                            const double             weight[],
                            const int                order[],
                            PDM_morton_code_t        rank_index[],
                            PDM_MPI_Comm             comm);


/**
 * \brief Dump a Morton to standard output or to a file
 *
 * \param [in]  dim           Spatial dimension (2 or 3)
 * \param [in]  code          Morton code to dump
 *
 */

void
PDM_morton_dump(int                 dim,
                PDM_morton_code_t   code);


/**
 * \brief Intersect a box with an array of sorted Morton codes
 *
 * A recursive top-down approach is used, starting from the deepest common
 * ancestor of all the Morton codes in the array.
 *
 * \param [in]     dim           Spatial dimension
 * \param [in]     node          Morton code of current node
 * \param [in]     box_min       Morton code of the box's lower corner
 * \param [in]     box_max       Morton code of the box's upper corner
 * \param [in]     nodes         Array of Morton codes
 * \param [in]     n_points      Number of points stored inside each node
 * \param [in]     start         Id of the first descendant of the current node
 * \param [in]     end           Id of the last descendant of the current node + 1
 * \param [in,out] n_intersect   Number of intersected nodes
 * \param [in,out] intersect     Intersected nodes
 *
 */

void
PDM_morton_intersect_box
(
 const int                dim,
 const PDM_morton_code_t  node,
 const PDM_morton_code_t  box_min,
 const PDM_morton_code_t  box_max,
 const PDM_morton_code_t  nodes[],
 int                     *n_points,
 const size_t             start,
 const size_t             end,
 size_t                  *n_intersect,
 int                     *intersect
 );


/**
 * \brief Get the closest node (described by Morton codes) of a given point.
 *
 * A recursive top-down approach is used, starting from the deepest common
 * ancestor of all the Morton codes in the array.
 * (NOT USED ANYMORE)
 *
 * \param [in]     dim            Spatial dimension
 * \param [in]     node           Morton code of current node
 * \param [in]     nodes          Array of Morton codes
 * \param [in]     point          Coordinates of the point
 * \param [in]     d              Normalization vector (dilatation component)
 * \param [in]     start          Id of the first descendant of the current node
 * \param [in]     end            Id of the last descendant of the current node + 1
 * \param [in,out] closest_node   Id of the closest node
 * \param [in,out] closest_dist2  Squared distance between the point and the closest node
 *
 */

void
PDM_morton_closest_node
(
 const int                dim,
 const PDM_morton_code_t  node,
 const PDM_morton_code_t  nodes[],
 const double             point[],
 const double             d[],
 const size_t             start,
 const size_t             end,
 int                     *closest_node,
 double                  *closest_dist2
 );

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MORTON_H__ */
