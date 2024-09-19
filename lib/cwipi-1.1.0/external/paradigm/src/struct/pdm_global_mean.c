
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_global_mean_priv.h"
#include "pdm_global_mean.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure that compute a global mean
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to \ref PDM_global_mean object
 */

PDM_global_point_mean_t*
PDM_global_mean_create
(
 const int          n_part,
 const PDM_MPI_Comm comm
)
{

  PDM_global_point_mean_t *gmean = (PDM_global_point_mean_t *) malloc(sizeof(PDM_global_point_mean_t));

  gmean->n_part  = n_part;
  gmean->comm    = comm;
  gmean->g_nums  = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);
  gmean->n_elts  = (int          *) malloc (sizeof(int          ) * n_part);
  gmean->strides = (int         **) malloc (sizeof(int         *) * n_part);
  gmean->ptb     = NULL;
  gmean->btp     = NULL;

  gmean->s_weight = NULL;

  for (int i = 0; i < n_part; i++) {
    gmean->g_nums [i] = NULL;
    gmean->strides[i] = NULL;
  }

  gmean->local_field       = (double **) malloc (sizeof(double * ) * n_part);
  gmean->local_weight      = (double **) malloc (sizeof(double * ) * n_part);
  gmean->global_mean_field = (double **) malloc (sizeof(double * ) * n_part);

  for (int i = 0; i < n_part; i++) {
    gmean->local_field      [i] = NULL;
    gmean->local_weight     [i] = NULL;
    gmean->global_mean_field[i] = NULL;
  }

  return gmean;
}


/**
 *
 * \brief Set absolute number
 *
 * \param [in]   gmean        Pointer to \ref PDM_global_mean object
 * \param [in]   i_part       Current partition
 * \param [in]   n_point      Number of points in the partition
 * \param [in]   numabs       Absolute number of points
 *
 */

void
PDM_global_mean_set
(
       PDM_global_point_mean_t *gmean,
 const int                      i_part,
 const int                      n_point,
 const PDM_g_num_t             *numabs
)
{
  gmean->g_nums [i_part] = (PDM_g_num_t *) numabs;
  gmean->n_elts [i_part] = n_point;
  gmean->strides[i_part] = malloc (sizeof(int) * n_point);
}

/**
 *
 * \brief Free a global point mean structure
 *
 * \param [in]   gmean           Pointer to \ref PDM_global_mean object
 *
 */

void
PDM_global_mean_free
(
 PDM_global_point_mean_t *gmean
)
{

  if (gmean->ptb != NULL) {
    gmean->ptb = PDM_part_to_block_free (gmean->ptb);
  }

  if (gmean->btp != NULL) {
    gmean->btp = PDM_block_to_part_free (gmean->btp);
  }

  free (gmean->g_nums);
  free (gmean->n_elts);
  free (gmean->local_field);
  free (gmean->local_weight);
  free (gmean->global_mean_field);

  for (int i = 0; i < gmean->n_part; i++) {
    if (gmean->strides[i] != NULL) {
      free (gmean->strides[i]);
    }
  }

  if (gmean->strides != NULL) {
    free (gmean->strides);
  }

  if (gmean->s_weight != NULL) {
    free (gmean->s_weight);
  }

  free(gmean);

}

/**
 *
 * \brief Set local field and it associated weight
 *
 * \param [in]   gmean                 Pointer to \ref PDM_global_mean object
 * \param [in]   i_part                Current partition
 * \param [in]   stride                Stride of the field
 * \param [in]   local_field           Local value of field
 * \param [in]   local_weight          Local weight used to compute the mean value
 * \param [in]   global_mean_field_ptr Pointer where global mean field
 *                                     will be stored after computing
 */

void
PDM_global_mean_field_set
(
 PDM_global_point_mean_t *gmean,
 const int                i_part,
 const int                stride,
 const double            *local_field,
 const double            *local_weight,
 double                  *global_mean_field_ptr
)
{
  gmean->local_field      [i_part] = (double *) local_field;
  gmean->local_weight     [i_part] = (double *) local_weight;
  gmean->global_mean_field[i_part] = (double *) global_mean_field_ptr;
  gmean->stride                    = stride;

}

/**
 *
 * \brief Compute the global average field
 *
 * \param [in]   gmean        Pointer to \ref PDM_global_mean object
 *
 */

void
PDM_global_mean_field_compute
(
 PDM_global_point_mean_t *gmean
)
{

  if (gmean->ptb == NULL) {
    gmean->ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                          1.,
                                          gmean->g_nums,
                                          NULL,
                                          gmean->n_elts,
                                          gmean->n_part,
                                          gmean->comm);

    int n_elt_block = PDM_part_to_block_n_elt_block_get(gmean->ptb);

    gmean->s_weight = malloc (sizeof(double) * n_elt_block);

  }

  if (gmean->btp == NULL) {
    PDM_g_num_t *distrib = PDM_part_to_block_distrib_index_get (gmean->ptb);
    gmean->btp = PDM_block_to_part_create (distrib,
                                          (const PDM_g_num_t **) gmean->g_nums,
                                          gmean->n_elts,
                                          gmean->n_part,
                                          gmean->comm);
  }

  int    *block_field_stride  = NULL;
  double *block_field         = NULL;
  int    *block_weight_stride = NULL;
  double *block_weight        = NULL;

  for (int i = 0; i < gmean->n_part; i++) {
    for (int j = 0; j < gmean->n_elts[i]; j++) {
      gmean->strides[i][j] = gmean->stride;
    }
  }

  PDM_part_to_block_exch (gmean->ptb,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                         gmean->strides,
                         (void **) gmean->local_field,
                                   &block_field_stride,
                         (void **) &block_field);

  int **_stride_w = NULL;
  if (gmean->local_weight[0] != NULL) {

    _stride_w = malloc (sizeof(int *) * gmean->n_part);
    for (int i = 0; i < gmean->n_part; i++) {
      _stride_w[i] = PDM_array_const_int(gmean->n_elts[i], 1);
    }

    PDM_part_to_block_exch (gmean->ptb,
                            sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            _stride_w,
                           (void **) gmean->local_weight,
                                     &block_weight_stride,
                           (void **) &block_weight);
    free (block_weight_stride);
  }

  //TODO: Remplisage du tableau moyenne

  int n_elt_block = PDM_part_to_block_n_elt_block_get(gmean->ptb);

  for (int i = 0; i < n_elt_block; i++) {
    gmean->s_weight[i] = 0.;
  }

  // Attention la definition de block_field_stride n'est pas la bonne
  // il faut creer un tableau stride_idx

  int *stride_idx = malloc(sizeof(int) * (n_elt_block + 1));
  stride_idx[0] = 0;

  for (int i = 0; i < n_elt_block; i++) {
    stride_idx[i+1] = stride_idx[i] + block_field_stride[i]/gmean->stride;
  }

  free (block_field_stride);

  for (int i = 0; i < n_elt_block; i++) {
    for (int j = stride_idx[i]; j < stride_idx[i+1]; j++) {
      double weight = 1.;
      if (block_weight != NULL) {
        weight = block_weight[j];
      }
      gmean->s_weight[i] += weight;
      for (int k = 0; k < gmean->stride; k++) {
        double val = block_field[j*gmean->stride + k];
        block_field[j*gmean->stride + k]  = 0;
        block_field[i*gmean->stride + k] += weight * val;
      }
    }
  }

  for (int i = 0; i <  n_elt_block; i++) {
    if (PDM_ABS(gmean->s_weight[i]) < 1e-15) {
      PDM_error (__FILE__, __LINE__, 0, "Sum of weights < 1e-15\n");
    }
    for (int k = 0; k < gmean->stride; k++) {
       block_field[gmean->stride * i + k] =
       block_field[gmean->stride * i + k] / gmean->s_weight[i];
    }
  }

  PDM_block_to_part_exch_in_place (gmean->btp,
                          sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          &gmean->stride,
                          block_field,
                          NULL,
                          (void **) gmean->global_mean_field);

  free (block_field);
  if (block_weight != NULL) {
    free (block_weight);
  }

  for (int i = 0; i < gmean->n_part; i++) {
    gmean->local_field      [i] = NULL;
    gmean->local_weight     [i] = NULL;
    gmean->global_mean_field[i] = NULL;

  }

  if (_stride_w != NULL) {
    for (int i = 0; i < gmean->n_part; i++) {
      free (_stride_w[i]);
    }
    free (_stride_w);
    _stride_w = NULL;
  }

  free (stride_idx);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
