/*
 * \file
 */

#ifndef __PDM_PART_H__
#define __PDM_PART_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"
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
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_PART_SPLIT_PARMETIS = 1,
  PDM_PART_SPLIT_PTSCOTCH = 2,
  PDM_PART_SPLIT_HILBERT  = 3
} PDM_part_split_t;

typedef struct _PDM_part_t PDM_part_t;
typedef struct _part_t     part_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a initial partitioning
 *  Build a initial partitioning from :
 *      - Cell block distribution with implicit global numbering
 *         (the first cell is the first cell of the first process and
 *          the latest cell is the latest cell of the latest process)
 *      - Face block distribution with implicit global numbering
 *      - Vertex block distribution with implicit global numbering
 *  To repart an existing partition use \ref PDM_part_repart function
 *
 * \param [in]   comm                   MPI Comminicator
 * \param [in]   split_method           Split method
 * \param [in]   renum_cell_method      Cell renumbering method
 * \param [in]   renum_face_method      Face renumbering method
 * \param [in]   n_property_cell        Number of cell properties
 * \param [in]   renum_properties_cell  For cache blocking [ n_cell_per_cache_wanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
 * \param [in]   n_property_face        Number of face properties
 * \param [in]   renum_properties_face  NOT USED
 * \param [in]   n_part                 Number of partition to build on this process
 * \param [in]   dn_cell                Number of distributed cells
 * \param [in]   dn_face                Number of distributed faces
 * \param [in]   dn_vtx                 Number of distributed vertices
 * \param [in]   n_face_group           Number of face groups
 * \param [in]   dcell_face_idx         Distributed cell face connectivity index or NULL
 *                                      (size : dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face             Distributed cell face connectivity or NULL
 *                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 * \param [in]   dcell_tag              Cell tag (size : n_cell) or NULL
 * \param [in]   dcell_weight           Cell weight (size : n_cell) or NULL
 * \param [in]   have_dcell_part        Presence of an array of cell part id
 * \param [in]   dcell_part             Distributed cell partitioning
 *                                      (size = dn_cell) or NULL (No partitioning if != NULL)
 * \param [in]   dface_cell             Distributed face cell connectivity or NULL
 *                                      (size : 2 * dn_face, numbering : 1 to n)
 * \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
 *                                      (size : dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx              Distributed face to vertex connectivity
 *                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 * \param [in]   dface_tag              Distributed face tag (size : dn_face)
 *                                      or NULL
 * \param [in]   dvtx_coord             Distributed vertex coordinates
 *                                      (size : 3*dn_vtx)
 * \param [in]   dvtx_tag               Distributed vertex tag (size : dn_vtx) or NULL
 * \param [in]   dface_group_idx        Index of distributed faces list of each group
 *                                      (size = n_face_group + 1) or NULL
 * \param [in]   dface_group            Distributed faces list of each group
 *                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
 *                                      or NULL
 *
 * \return    Pointer to \ref PDM_part object
 *
 */

PDM_part_t *
PDM_part_create
(
 const PDM_MPI_Comm           comm,
 const PDM_part_split_t       split_method,
 const char                  *renum_cell_method,
 const char                  *renum_face_method,
 const int                    n_property_cell,
 const int                   *renum_properties_cell,
 const int                    n_property_face,
 const int                   *renum_properties_face,
 const int                    n_part,
 const int                    dn_cell,
 const int                    dn_face,
 const int                    dn_vtx,
 const int                    n_face_group,
 const int                   *dcell_face_idx,
 const PDM_g_num_t           *dcell_face,
 const int                   *dcell_tag,
 const int                   *dcell_weight,
 const int                    have_dcell_part,
       int                   *dcell_part,
 const PDM_g_num_t           *dface_cell,
 const int                   *dface_vtx_idx,
 const PDM_g_num_t           *dface_vtx,
 const int                   *dface_tag,
 const double                *dvtx_coord,
 const int                   *dvtx_tag,
 const int                   *dface_group_idx,
 const PDM_g_num_t           *dface_group
 );

/**
 *
 * \brief Return a mesh partition dimensions
 *
 * \param [in]   ppart               Pointer to \ref PDM_part object
 * \param [in]   i_part              Current partition
 * \param [out]  n_cell              Number of cells
 * \param [out]  n_face              Number of faces
 * \param [out]  n_face_part_bound   Number of partitioning boundary faces
 * \param [out]  n_vtx               Number of vertices
 * \param [out]  n_proc              Number of processus
 * \param [out]  n_total_part        Number of partitions
 * \param [out]  scell_face          Size of cell-face connectivity
 * \param [out]  sface_vtx           Size of face-vertex connectivity
 * \param [out]  sFacePartBound      Size of face_part_bound array
 * \param [out]  sface_group         Size of face_group array
 * \param [out]  n_face_group        Number of face groups
 *
 */

void
PDM_part_part_dim_get
(
PDM_part_t   *ppart,
const  int    i_part,
       int   *n_cell,
       int   *n_face,
       int   *n_face_part_bound,
       int   *n_vtx,
       int   *n_proc,
       int   *n_total_part,
       int   *scell_face,
       int   *sface_vtx,
       int   *sface_group,
       int   *n_face_group
);


/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   ppart                     Pointer to \ref PDM_part object
 * \param [in]   i_part                    Current partition
 * \param [out]  cell_tag                  Cell tag (size = n_cell)
 * \param [out]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 * \param [out]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 *                                                                   numbering : 1 to n)
 * \param [out]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 * \param [out]  face_tag                  Face tag (size = n_face)
 * \param [out]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 * \param [out]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 * \param [out]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 * \param [out]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 * \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
 * \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
 * \param [out]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
 *                                          sorted by processus, sorted by partition in each processus, and
 *                                          sorted by absolute face number in each partition
 *                                         For each face :
 *                                           - Face local number (numbering : 1 to n)
 *                                           - Connected process (numbering : 0 to n-1)
 *                                           - Connected Partition
 *                                             on the connected process (numbering :1 to n)
 *                                           - Connected face local number
 *                                             in the connected partition (numbering :1 to n)
 * \param [out]  vtx_tag                   Vertex tag (size = nVertex)
 * \param [out]  vtx                       Vertex coordinates (size = 3 * nVertex)
 * \param [out]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 * \param [out]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 * \param [out]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 * \param [out]  face_group_ln_to_gn       Faces global numbering for each group
 *                                         (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 *
 */

void PDM_part_part_val_get
(
PDM_part_t          *ppart,
const int            i_part,
      int          **cell_tag,
      int          **cell_face_idx,
      int          **cell_face,
      PDM_g_num_t  **cell_ln_to_gn,
      int          **face_tag,
      int          **face_cell,
      int          **face_vtx_idx,
      int          **face_vtx,
      PDM_g_num_t  **face_ln_to_gn,
      int          **face_part_bound_proc_idx,
      int          **face_part_bound_part_idx,
      int          **face_part_bound,
      int          **vtx_tag,
      double       **vtx,
      PDM_g_num_t  **vtx_ln_to_gn,
      int          **face_group_idx,
      int          **face_group,
      PDM_g_num_t  **face_group_ln_to_gn
);


/**
 *
 * \brief Return the coloring of a mesh partition
 *
 * \param [in]   ppart               Pointer to \ref PDM_part object
 * \param [in]   i_part              Current partition
 * \param [out]  cell_color          Cell color (size = n_cell)
 * \param [out]  face_color          Face color (size = n_face)
 * \param [out]  thread_color        Thread color (size = n_cell)
 * \param [out]  hyperplane_color    Hyperplane color (size = n_cell)
 *
 */

void PDM_part_part_color_get
(
PDM_part_t          *ppart,
const int            i_part,
      int          **cell_color,
      int          **face_color,
      int          **thread_color,
      int          **hyperplane_color
);


/**
 *
 * \brief Free ppart
 *
 * \param [in]   ppart               Pointer to \ref PDM_part object
 *
 */

void
PDM_part_free
(
PDM_part_t  *ppart
);


/**
 *
 * \brief Return times
 *
 * \param [in]   ppart       Pointer to \ref PDM_part object
 * \param [out]  elapsed     Elapsed times (size = 4)
 * \param [out]  cpu         Cpu times (size = 4)
 * \param [out]  cpu_user    User cpu times (size = 4)
 * \param [out]  cpu_sys     System cpu times (size = 4)
 *
 */

void
PDM_part_time_get
(
 PDM_part_t  *ppart,
 double     **elapsed,
 double     **cpu,
 double     **cpu_user,
 double     **cpu_sys
 );


/**
 *
 * \brief Return statistics
 *
 * \param [in]   ppart                          Pointer to \ref PDM_part object
 * \param [out]  cells_average                  average of cells number
 * \param [out]  cells_median                   median of cells number
 * \param [out]  cells_std_deviation            standard deviation of cells number
 * \param [out]  cells_min                      minimum of cells nummber
 * \param [out]  cells_max                      maximum of cells nummber
 * \param [out]  bound_part_faces_average       average of partitioning boundary faces
 * \param [out]  bound_part_faces_median        median of partitioning boundary faces
 * \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 * \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 * \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 *
 */

void
PDM_part_stat_get
(
PDM_part_t  *ppart,
int         *cells_average,
int         *cells_median,
double      *cells_std_deviation,
int         *cells_min,
int         *cells_max,
int         *bound_part_faces_average,
int         *bound_part_faces_median,
double      *bound_part_faces_std_deviation,
int         *bound_part_faces_min,
int         *bound_part_faces_max,
int         *bound_part_faces_sum
);


/**
 *
 * \brief Free partially a \ref PDM_part_object
 *
 * \param [in]   ppart       Pointer to \ref PDM_part object
 *
 */

void
PDM_part_partial_free
(
PDM_part_t  *ppart
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_part_H__ */
