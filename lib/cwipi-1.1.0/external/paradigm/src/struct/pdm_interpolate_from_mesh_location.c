/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_interpolate_from_mesh_location_priv.h"
#include "pdm_interpolate_from_mesh_location.h"
#include "pdm_timer.h"
#include "pdm_logging.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that compute interpolation from mesh_location information
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

PDM_interpolate_from_mesh_location_t*
PDM_interpolate_from_mesh_location_create
(
 const int                    n_cloud_target,
       PDM_interpolate_kind_t interp_kind,
 const PDM_MPI_Comm           comm
)
{
  PDM_interpolate_from_mesh_location_t *interp_from_ml = (PDM_interpolate_from_mesh_location_t *) malloc (sizeof(PDM_interpolate_from_mesh_location_t));

  interp_from_ml->comm = comm;

  PDM_MPI_Comm_size (comm, &interp_from_ml->n_rank);
  PDM_MPI_Comm_rank (comm, &interp_from_ml->i_rank);

  interp_from_ml->n_cloud_target     = n_cloud_target;

  interp_from_ml->point_clouds       = malloc( sizeof(_point_cloud_t      ) * interp_from_ml->n_cloud_target );
  interp_from_ml->points_in_elements = malloc( sizeof(_points_in_element_t) * interp_from_ml->n_cloud_target );

  PDM_UNUSED(interp_kind);

  return interp_from_ml;
}

void
PDM_interpolate_from_mesh_location_free
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml
)
{

  //printf("PDM_interpolate_from_mesh_location_free \n");

  for (int icloud = 0; icloud < interp_from_ml->n_cloud_target; icloud++) {

    _points_in_element_t *_points_in_elements = interp_from_ml->points_in_elements + icloud;

    free (_points_in_elements->n_elts);
    free (_points_in_elements->pts_inside_idx);
    free (_points_in_elements->gnum);
    free (_points_in_elements->uvw);
    free (_points_in_elements->coords);
    free (_points_in_elements->projected_coords);
    free (_points_in_elements->weights_idx);
    free (_points_in_elements->weights);
    free (_points_in_elements->dist2);

    _point_cloud_t *pcloud = interp_from_ml->point_clouds + icloud;

    if (pcloud->n_points != NULL) {
      free (pcloud->n_points);
    }

    if (pcloud->coords != NULL) {
      free (pcloud->coords);
    }

    if (pcloud->gnum != NULL) {
      free (pcloud->gnum);
    }

  }

  free(interp_from_ml->n_cell       );
  free(interp_from_ml->cell_face_idx);
  free(interp_from_ml->cell_face    );
  free(interp_from_ml->cell_ln_to_gn);
  free(interp_from_ml->n_face       );
  free(interp_from_ml->face_vtx_idx );
  free(interp_from_ml->face_vtx     );
  free(interp_from_ml->face_ln_to_gn);
  free(interp_from_ml->n_vtx        );
  free(interp_from_ml->coords       );
  free(interp_from_ml->vtx_ln_to_gn );


  free(interp_from_ml->point_clouds);
  free(interp_from_ml->points_in_elements);
  free(interp_from_ml);
}


void
PDM_interpolate_from_mesh_location_compute
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml
)
{
  PDM_UNUSED(interp_from_ml);

  /*
   * On prepare les échanges pour faire le part_to_part
   */

}

void
PDM_interpolate_from_mesh_location_exch
(
 PDM_interpolate_from_mesh_location_t   *interp_from_ml,
 int                                     i_point_cloud,
 size_t                                  s_data,
 double                                **part_data_in,
 double                               ***cloud_data_out
)
{
  PDM_UNUSED(interp_from_ml);
  PDM_UNUSED(s_data);
  PDM_UNUSED(part_data_in);
  PDM_UNUSED(cloud_data_out);

  assert (interp_from_ml->points_in_elements != NULL);

  _points_in_element_t *_points_in_elements = interp_from_ml->points_in_elements + i_point_cloud;

  /*
   * For now only first order with
   */
  double** cloud_data_in_current_src = (double **) malloc( interp_from_ml->n_part_src * sizeof(double *));
  int    **cloud_data_in_current_src_n = (int    **) malloc( interp_from_ml->n_part_src * sizeof(int    *));
  assert(_points_in_elements->n_part == interp_from_ml->n_part_src );

  _point_cloud_t *pcloud = interp_from_ml->point_clouds + i_point_cloud;

  int* n_point_tot = (int *) malloc( interp_from_ml->n_part_src * sizeof(int));

  for(int i_part = 0; i_part < interp_from_ml->n_part_src; ++i_part){

    int n_cell = interp_from_ml->n_cell[i_part];
    int n_elmt = n_cell; // By construction of _points_in_element_t
    int* _elt_pts_inside_idx = _points_in_elements->pts_inside_idx[i_part];

    n_point_tot[i_part] = _elt_pts_inside_idx[n_elmt];

    //printf("n_cell[%i] = %i \n", i_part, n_cell);
    //printf("_elt_pts_inside_idx[%i] = %i \n", i_part, _elt_pts_inside_idx[1]);
    //printf("n_point_tot[%i] = %i \n", i_part, n_point_tot[i_part]);
    // printf("pcloud->n_points[%i] = %i \n", i_part, pcloud->n_points[i_part]);

    cloud_data_in_current_src  [i_part] = (double *) malloc( _elt_pts_inside_idx[n_elmt] * sizeof(double));
    cloud_data_in_current_src_n[i_part] = (int    *) malloc( _elt_pts_inside_idx[n_elmt] * sizeof(int   ));

    for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
      for (int i_point = _elt_pts_inside_idx[i_cell]; i_point < _elt_pts_inside_idx[i_cell+1]; i_point++) {
        cloud_data_in_current_src  [i_part][i_point] = part_data_in[i_part][i_cell]; // Simple extrapolation
        cloud_data_in_current_src_n[i_part][i_point] = 1;
        // printf(" cloud_data_in_current_src[%i][%i] = %12.5e (from cell = %i) | gnum = %i \n", i_part, i_point, part_data_in[i_part][i_cell], i_cell, (int)_points_in_elements->gnum[i_part][i_point] );
      }
    }
  }

  PDM_g_num_t n_g_cloud = 0;

  for(int i_part = 0; i_part < pcloud->n_part; ++i_part){
    for(int i = 0; i < pcloud->n_points[i_part]; ++i) {
      n_g_cloud = PDM_MAX(pcloud->gnum[i_part][i], n_g_cloud);
    }
  }
  PDM_g_num_t _n_g_cloud = 0;
  PDM_MPI_Allreduce (&n_g_cloud, &_n_g_cloud, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, interp_from_ml->comm);
  n_g_cloud = _n_g_cloud;

  /*
   *  Create the part_to_block to have block of cloud point
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP, /* Dans notre cas on a une seule localisation possible */
                                                       1.,
                                                       _points_in_elements->gnum,      /* Numero absolu des points (donc les target ) */
                                                       NULL,
                                                       n_point_tot,
                                                       interp_from_ml->n_part_src,
                                                       interp_from_ml->comm);


  int stride_one = 1;
  int    *block_strid;
  double *block_data;
  int s_block_data = PDM_part_to_block_exch(ptb,
                                            sizeof(double),
                                            PDM_STRIDE_VAR_INTERLACED,
                                            stride_one,
                                            cloud_data_in_current_src_n,
                                  (void **) cloud_data_in_current_src,
                                           &block_strid,
                                  (void **)&block_data);

  if(0 == 1) {
    printf(" s_block_data = %i \n", s_block_data);
    for(int i_point = 0; i_point < s_block_data; ++i_point) {
      printf("block_data[%i] = %12.5e \n", i_point, block_data[i_point]);
    }
  }


  PDM_g_num_t* point_block_distrib_idx = PDM_part_to_block_adapt_partial_block_to_block(ptb, &block_strid, n_g_cloud);

  if(1 == 1) {
    int tmp = point_block_distrib_idx[interp_from_ml->i_rank+1] - point_block_distrib_idx[interp_from_ml->i_rank];
    printf("tmp = %i \n", tmp);
    for(int i_point = 0; i_point < tmp; ++i_point) {
      printf("block_strid[%i] = %i \n", i_point, block_strid[i_point]);
    }
  }

  PDM_block_to_part_t *btp = PDM_block_to_part_create(point_block_distrib_idx,
                               (const PDM_g_num_t **) pcloud->gnum,
                                                      pcloud->n_points,
                                                      pcloud->n_part,
                                                      interp_from_ml->comm);

  int** part_strid = NULL;
  PDM_block_to_part_exch(btp,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_strid,
                          block_data,
               (int ***)  &part_strid,
               (void ***) cloud_data_out);

  PDM_part_to_block_free(ptb);
  free(block_data);
  free(block_strid);

  for(int i_part = 0; i_part < pcloud->n_part; ++i_part){
    PDM_log_trace_array_double((*cloud_data_out)[i_part], pcloud->n_points[i_part], "cloud_data_out :: ");
    PDM_log_trace_array_long(pcloud->gnum[i_part], pcloud->n_points[i_part], "cloud_gnum :: ");
  }

  PDM_block_to_part_free(btp);

  for(int i_part = 0; i_part < pcloud->n_part; ++i_part){
    free(part_strid[i_part]);
  }
  free(part_strid);

  for(int i_part = 0; i_part < interp_from_ml->n_part_src; ++i_part){
    free(cloud_data_in_current_src[i_part]);
    free(cloud_data_in_current_src_n[i_part]);
  }
  free(cloud_data_in_current_src);
  free(cloud_data_in_current_src_n);
  free(n_point_tot);

}



void
PDM_interpolate_from_mesh_location_exch_inplace
(
 PDM_interpolate_from_mesh_location_t   *interp_from_ml,
 int                                     i_point_cloud,
 size_t                                  s_data,
 double                                **part_data_in,
 double                                **cloud_data_out
)
{
  PDM_UNUSED(interp_from_ml);
  PDM_UNUSED(s_data);
  PDM_UNUSED(part_data_in);
  PDM_UNUSED(cloud_data_out);

  assert (interp_from_ml->points_in_elements != NULL);

  _points_in_element_t *_points_in_elements = interp_from_ml->points_in_elements + i_point_cloud;

  /*
   * For now only first order with
   */
  double **cloud_data_in_current_src   = (double **) malloc( interp_from_ml->n_part_src * sizeof(double *));
  int    **cloud_data_in_current_src_n = (int    **) malloc( interp_from_ml->n_part_src * sizeof(int    *));
  assert(_points_in_elements->n_part == interp_from_ml->n_part_src );

  _point_cloud_t *pcloud = interp_from_ml->point_clouds + i_point_cloud;

  int* n_point_tot = (int *) malloc( interp_from_ml->n_part_src * sizeof(int));

  for(int i_part = 0; i_part < interp_from_ml->n_part_src; ++i_part){

    int n_elmt = _points_in_elements->n_elts[i_part]; // By construction of _points_in_element_t
    int* _elt_pts_inside_idx = _points_in_elements->pts_inside_idx[i_part];

    n_point_tot[i_part] = _elt_pts_inside_idx[n_elmt];

    // printf("n_point_tot[%i] = %i \n", i_part, n_point_tot[i_part]);

    cloud_data_in_current_src  [i_part] = (double *) malloc( _elt_pts_inside_idx[n_elmt] * sizeof(double));
    cloud_data_in_current_src_n[i_part] = (int    *) malloc( _elt_pts_inside_idx[n_elmt] * sizeof(int   ));

    for(int i_cell = 0; i_cell < n_elmt; ++i_cell) {
      for (int i_point = _elt_pts_inside_idx[i_cell]; i_point < _elt_pts_inside_idx[i_cell+1]; i_point++) {
        cloud_data_in_current_src  [i_part][i_point] = part_data_in[i_part][i_cell]; // Simple extrapolation
        // printf(" cloud_data_in_current_src[%i][%i] = %12.5e (from cell = %i) | gnum = %i \n", i_part, i_point, part_data_in[i_part][i_cell], i_cell, (int)_points_in_elements->gnum[i_part][i_point] );
        cloud_data_in_current_src_n[i_part][i_point] = 1;
      }
    }
  }

  PDM_g_num_t n_g_cloud = 0;

  for(int i_part = 0; i_part < pcloud->n_part; ++i_part){
    for(int i = 0; i < pcloud->n_points[i_part]; ++i) {
      n_g_cloud = PDM_MAX(pcloud->gnum[i_part][i], n_g_cloud);
    }
  }
  PDM_g_num_t _n_g_cloud = 0;
  PDM_MPI_Allreduce (&n_g_cloud, &_n_g_cloud, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, interp_from_ml->comm);
  n_g_cloud = _n_g_cloud;

  /*
   *  Create the part_to_block to have block of cloud point
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP, /* Dans notre cas on a une seule localisation possible */
                                                       1.,
                                                       _points_in_elements->gnum,      /* Numero absolu des points (donc les target ) */
                                                       NULL,
                                                       n_point_tot,
                                                       interp_from_ml->n_part_src,
                                                       interp_from_ml->comm);


  int stride_one = 1;
  int    *block_strid;
  double *block_data;
  int s_block_data = PDM_part_to_block_exch(ptb,
                                            sizeof(double),
                                            PDM_STRIDE_VAR_INTERLACED,
                                            stride_one,
                                            cloud_data_in_current_src_n,
                                  (void **) cloud_data_in_current_src,
                                           &block_strid,
                                  (void **)&block_data);

  if(0 == 1) {
    printf(" s_block_data = %i \n", s_block_data);
    for(int i_point = 0; i_point < s_block_data; ++i_point) {
      printf("block_data[%i] = %12.5e \n", i_point, block_data[i_point]);
      printf("block_strid[%i] = %i \n", i_point, block_strid[i_point]);
    }
  }

  PDM_g_num_t* point_block_distrib_idx = PDM_part_to_block_adapt_partial_block_to_block(ptb, &block_strid, n_g_cloud);

  if(0 == 1) {
    int tmp = point_block_distrib_idx[interp_from_ml->i_rank+1] - point_block_distrib_idx[interp_from_ml->i_rank];
    printf("tmp = %i \n", tmp);
    for(int i_point = 0; i_point < tmp; ++i_point) {
      printf("block_strid[%i] = %i \n", i_point, block_strid[i_point]);
    }
  }
  PDM_block_to_part_t *btp = PDM_block_to_part_create(point_block_distrib_idx,
                               (const PDM_g_num_t **) pcloud->gnum,
                                                      pcloud->n_points,
                                                      pcloud->n_part,
                                                      interp_from_ml->comm);

  // for(int i_part = 0; i_part < pcloud->n_part; ++i_part){
  //   PDM_log_trace_array_double(cloud_data_out[i_part], pcloud->n_points[i_part], "cloud_data_out :: ");
  //   PDM_log_trace_array_long(pcloud->gnum[i_part], pcloud->n_points[i_part], "cloud_gnum :: ");
  // }

  int** part_strid;
  double** tmp_cloud_data_out;
  PDM_block_to_part_exch(btp,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_strid,
                          block_data,
                          &part_strid,
               (void ***) &tmp_cloud_data_out);

  /*
   * Recopie dans le vrai tableau
   */
  for(int i_part = 0; i_part < pcloud->n_part; ++i_part){
    int idx = 0;
    for(int i = 0; i < pcloud->n_points[i_part]; ++i) {
      int n_interp = part_strid[i_part][i];
      if(n_interp == 1) { // Normalement ne peux valoir que 1 ou 0
        cloud_data_out[i_part][i] = tmp_cloud_data_out[i_part][idx++];
      }
      assert(n_interp <= 1);
    }
  }

  PDM_part_to_block_free(ptb);
  free(block_data);
  free(block_strid);

  PDM_block_to_part_free(btp);

  for(int i_part = 0; i_part < pcloud->n_part; ++i_part){
    free(part_strid[i_part]);
    free(tmp_cloud_data_out[i_part]);
  }
  free(part_strid);
  free(tmp_cloud_data_out);

  for(int i_part = 0; i_part < interp_from_ml->n_part_src; ++i_part){
    free(cloud_data_in_current_src[i_part]);
    free(cloud_data_in_current_src_n[i_part]);
  }
  free(cloud_data_in_current_src);
  free(cloud_data_in_current_src_n);
  free(n_point_tot);
  free(point_block_distrib_idx);

}

void
PDM_interpolate_from_mesh_location_send
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml,
 size_t                                 s_data,
 void                                  **part_data_in,
 void                                 ***cloud_data_out
)
{
  PDM_UNUSED(interp_from_ml);
  PDM_UNUSED(s_data);
  PDM_UNUSED(part_data_in);
  PDM_UNUSED(cloud_data_out);
}



void
PDM_interpolate_from_mesh_location_recv
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml,
 size_t                                 s_data,
 void                                  **part_data_in,
 void                                 ***cloud_data_out
)
{
  PDM_UNUSED(interp_from_ml);
  PDM_UNUSED(s_data);
  PDM_UNUSED(part_data_in);
  PDM_UNUSED(cloud_data_out);
}



/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_interpolate_from_mesh_location_mesh_global_data_set
(
       PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                                   n_part
)
{
  interp_from_ml->n_part_src          = n_part;

  for (int icloud = 0; icloud < interp_from_ml->n_cloud_target; icloud++) {

    _points_in_element_t *_points_in_elements = interp_from_ml->points_in_elements + icloud;

    _points_in_elements->n_part           = n_part;
    _points_in_elements->n_elts           = (int          *) malloc( interp_from_ml->n_part_src * sizeof(int          ));
    _points_in_elements->pts_inside_idx   = (int         **) malloc( interp_from_ml->n_part_src * sizeof(int         *));
    _points_in_elements->gnum             = (PDM_g_num_t **) malloc( interp_from_ml->n_part_src * sizeof(PDM_g_num_t *));
    _points_in_elements->uvw              = (double      **) malloc( interp_from_ml->n_part_src * sizeof(double      *));
    _points_in_elements->coords           = (double      **) malloc( interp_from_ml->n_part_src * sizeof(double      *));
    _points_in_elements->projected_coords = (double      **) malloc( interp_from_ml->n_part_src * sizeof(double      *));
    _points_in_elements->weights_idx      = (int         **) malloc( interp_from_ml->n_part_src * sizeof(int         *));
    _points_in_elements->weights          = (double      **) malloc( interp_from_ml->n_part_src * sizeof(double      *));
    _points_in_elements->dist2            = (double      **) malloc( interp_from_ml->n_part_src * sizeof(double      *));

  }

  interp_from_ml->n_cell        = (int          *) malloc(n_part * sizeof(int     ));
  interp_from_ml->cell_face_idx = (int         **) malloc(n_part * sizeof(int    *));
  interp_from_ml->cell_face     = (int         **) malloc(n_part * sizeof(int    *));
  interp_from_ml->cell_ln_to_gn = (PDM_g_num_t **) malloc(n_part * sizeof(int    *));
  interp_from_ml->n_face        = (int          *) malloc(n_part * sizeof(int     ));
  interp_from_ml->face_vtx_idx  = (int         **) malloc(n_part * sizeof(int    *));
  interp_from_ml->face_vtx      = (int         **) malloc(n_part * sizeof(int    *));
  interp_from_ml->face_ln_to_gn = (PDM_g_num_t **) malloc(n_part * sizeof(int    *));
  interp_from_ml->n_vtx         = (int          *) malloc(n_part * sizeof(int     ));
  interp_from_ml->coords        = (double      **) malloc(n_part * sizeof(double *));
  interp_from_ml->vtx_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(int    *));

}

void
PDM_interpolate_from_mesh_location_n_part_cloud_set
(
       PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                                   i_point_cloud,
 const int                                   n_part
)
{
  // printf("PDM_interpolate_from_mesh_location_n_part_cloud_set : %i --> %i \n", i_point_cloud, n_part);

  _point_cloud_t *pcloud = interp_from_ml->point_clouds + i_point_cloud;

  pcloud->n_part   = n_part;
  pcloud->n_points = malloc( n_part * sizeof(int          ));
  pcloud->coords   = malloc( n_part * sizeof(double *     ));
  pcloud->gnum     = malloc( n_part * sizeof(PDM_g_num_t *));

  for (int i_part = 0; i_part < n_part; i_part++) {
    pcloud->n_points[i_part] = -1;
    pcloud->coords  [i_part] = NULL;
    pcloud->gnum    [i_part] = NULL;
  }
}

/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_interpolate_from_mesh_location_cloud_set
(
       PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                                   i_point_cloud,
 const int                                   i_part,
 const int                                   n_points,
       double                               *coords,
       PDM_g_num_t                          *gnum
)
{
  interp_from_ml->point_clouds[i_point_cloud].n_points[i_part] = n_points;
  interp_from_ml->point_clouds[i_point_cloud].coords  [i_part] = coords;
  interp_from_ml->point_clouds[i_point_cloud].gnum    [i_part] = gnum;

}

/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_face_idx Index in the cell -> face connectivity
 * \param [in]   cell_face     cell -> face connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */
void
PDM_interpolate_from_mesh_location_part_set
(
       PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                                   i_part,
 const int                                   n_cell,
 const int                                  *cell_face_idx,
 const int                                  *cell_face,
 const PDM_g_num_t                          *cell_ln_to_gn,
 const int                                   n_face,
 const int                                  *face_vtx_idx,
 const int                                  *face_vtx,
 const PDM_g_num_t                          *face_ln_to_gn,
 const int                                   n_vtx,
 const double                               *coords,
 const PDM_g_num_t                          *vtx_ln_to_gn
)
{
  interp_from_ml->n_cell       [i_part] = n_cell;
  interp_from_ml->cell_face_idx[i_part] = (int         *) cell_face_idx;
  interp_from_ml->cell_face    [i_part] = (int         *) cell_face;
  interp_from_ml->cell_ln_to_gn[i_part] = (PDM_g_num_t *) cell_ln_to_gn;
  interp_from_ml->n_face       [i_part] = n_face;
  interp_from_ml->face_vtx_idx [i_part] = (int         *) face_vtx_idx;
  interp_from_ml->face_vtx     [i_part] = (int         *) face_vtx;
  interp_from_ml->face_ln_to_gn[i_part] = (PDM_g_num_t *) face_ln_to_gn;
  interp_from_ml->n_vtx        [i_part] = n_vtx;
  interp_from_ml->coords       [i_part] = (double *) coords;
  interp_from_ml->vtx_ln_to_gn [i_part] = (PDM_g_num_t *) vtx_ln_to_gn;
}


/**
 *
 * \brief Set point list located in elements
 *
 * \param [in]   id                      Identifier
 * \param [in]   i_part                  Index of partition of the mesh
 * \param [in]   i_point_cloud           Index of cloud
 * \param [in]   elt_pts_inside_idx      Points index (size = n_elt + 1)
 * \param [in]   points_gnum             Points global number
 * \param [in]   points_coords           Points coordinates
 * \param [in]   points_uvw              Points parametric coordinates in elements
 * \param [in]   points_weights_idx      Interpolation weights index (size = elt_pts_inside_idx[n_elt] + 1)
 * \param [in]   points_weights          Interpolation weights
 * \param [in]   points_dist2            Distance element-points (dist < 0 if the point is inside)
 * \param [in]   points_projected_coords Point projection on element if the point is outside
 *
 */
void
PDM_interpolate_from_mesh_location_points_in_elt_set
(
 PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                             i_part,
 const int                             i_point_cloud,
 const int                             n_elts,
 int                                  *elt_pts_inside_idx,
 PDM_g_num_t                          *points_gnum,
 double                               *points_coords,
 double                               *points_uvw,
 int                                  *points_weights_idx,
 double                               *points_weights,
 double                               *points_dist2,
 double                               *points_projected_coords
)
{
  assert (i_point_cloud < interp_from_ml->n_cloud_target);

  _points_in_element_t *_points_in_elements = interp_from_ml->points_in_elements + i_point_cloud;

  // printf("PDM_interpolate_from_mesh_location_points_in_elt_set : i_part = %i | _points_in_elements->n_part = %i \n", i_part, _points_in_elements->n_part);

  assert (interp_from_ml->points_in_elements != NULL);
  assert (i_part < _points_in_elements->n_part);

  _points_in_elements->n_elts          [i_part] = n_elts;
  _points_in_elements->pts_inside_idx  [i_part] = elt_pts_inside_idx;
  _points_in_elements->gnum            [i_part] = points_gnum;
  _points_in_elements->coords          [i_part] = points_coords;
  _points_in_elements->uvw             [i_part] = points_uvw;
  _points_in_elements->weights_idx     [i_part] = points_weights_idx;
  _points_in_elements->weights         [i_part] = points_weights;
  _points_in_elements->dist2           [i_part] = points_dist2;
  _points_in_elements->projected_coords[i_part] = points_projected_coords;

}

#ifdef  __cplusplus
}
#endif
