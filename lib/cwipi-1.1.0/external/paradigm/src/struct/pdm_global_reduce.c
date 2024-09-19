
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

#include "pdm_global_reduce_priv.h"
#include "pdm_global_reduce.h"

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

static inline double _min_a_b (const double a, const double b) {
  return PDM_MIN (a, b);
}

static inline double _max_a_b (const double a, const double b) {
  return PDM_MAX (a, b);
}

static inline double _sum_a_b (const double a, const double b) {
  return a + b;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that computes a global reduction
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to \ref PDM_global_reduce object
 */

PDM_global_reduce_t *
PDM_global_reduce_create
(
 const int          n_part,
 const PDM_MPI_Comm comm
)
{
  PDM_global_reduce_t *gre = (PDM_global_reduce_t *) malloc (sizeof(PDM_global_reduce_t));

  gre->n_part  = n_part;
  gre->comm    = comm;
  gre->g_nums  = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);
  gre->n_elts  = (int          *) malloc (sizeof(int          ) * n_part);
  gre->strides = (int         **) malloc (sizeof(int         *) * n_part);
  gre->ptb     = NULL;
  gre->btp     = NULL;

  gre->operation = PDM_REDUCE_OP_MIN;

  for (int i = 0; i < n_part; i++) {
    gre->g_nums [i] = NULL;
    gre->strides[i] = NULL;
  }

  gre->local_field          = (double **) malloc (sizeof(double * ) * n_part);
  gre->global_reduced_field = (double **) malloc (sizeof(double * ) * n_part);

  for (int i = 0; i < n_part; i++) {
    gre->local_field         [i] = NULL;
    gre->global_reduced_field[i] = NULL;
  }

  return gre;
}



/**
 *
 * \brief Free a global point mean structure
 *
 * \param [in]   gre          Pointer to \ref PDM_global_reduce object
 *
 */

void
PDM_global_reduce_free
(
 PDM_global_reduce_t *gre
)
{
  if (gre->ptb != NULL) {
    gre->ptb = PDM_part_to_block_free (gre->ptb);
  }

  if (gre->btp != NULL) {
    gre->btp = PDM_block_to_part_free (gre->btp);
  }

  free (gre->g_nums);
  free (gre->n_elts);
  free (gre->local_field);
  free (gre->global_reduced_field);

  for (int i = 0; i < gre->n_part; i++) {
    if (gre->strides[i] != NULL) {
      free (gre->strides[i]);
    }
  }

  if (gre->strides != NULL) {
    free (gre->strides);
  }

  free(gre);
}


/**
 *
 * \brief Set global ids
 *
 * \param [in]   gre           Pointer to \ref PDM_global_reduce object
 * \param [in]   i_part        Current partition
 * \param [in]   n_pts         Number of points in the partition
 * \param [in]   pts_ln_to_gn  Global ids of points in the partition
 *
 */

void
PDM_global_reduce_g_num_set
(
 PDM_global_reduce_t *gre,
 const int            i_part,
 const int            n_pts,
 const PDM_g_num_t   *pts_ln_to_gn
 )
 {
  gre->g_nums [i_part] = (PDM_g_num_t *) pts_ln_to_gn;
  gre->n_elts [i_part] = n_pts;
  gre->strides[i_part] = malloc (sizeof(int) * n_pts);
 }


/**
 *
 * \brief Set reduction operation
 *
 * \param [in]   gre                       Pointer to \ref PDM_global_reduce object
 * \param [in]   operation                 Type of reduction operation
 */

void
PDM_global_reduce_operation_set
(
 PDM_global_reduce_t   *gre,
 const PDM_reduce_op_t  operation
)
{
  gre->operation = operation;
}


 /**
 *
 * \brief Set local field
 *
 * \param [in]   gre                       Pointer to \ref PDM_global_reduce object
 * \param [in]   i_part                    Current partition
 * \param [in]   stride                    Stride of the field
 * \param [in]   local_field               Local value of field
 * \param [in]   global_reduced_field_ptr  Pointer where global reduced field
 *                                         will be stored after computing
 */

void
PDM_global_reduce_field_set
(
 PDM_global_reduce_t *gre,
 const int            i_part,
 const int            stride,
 const double        *local_field,
 double              *global_reduced_field_ptr
)
{
  gre->local_field         [i_part] = (double *) local_field;
  gre->global_reduced_field[i_part] = (double *) global_reduced_field_ptr;
  gre->stride                       = stride;
}



/**
 *
 * \brief Compute the global reduced field
 *
 * \param [in]   gre          Pointer to \ref PDM_global_reduce object
 *
 */

void
PDM_global_reduce_field_compute
(
 PDM_global_reduce_t *gre
)
{
  if (gre->ptb == NULL) {
    gre->ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                         PDM_PART_TO_BLOCK_POST_MERGE,
                                         1.,
                                         gre->g_nums,
                                         NULL,
                                         gre->n_elts,
                                         gre->n_part,
                                         gre->comm);
  }

  if (gre->btp == NULL) {
    PDM_g_num_t *distrib = PDM_part_to_block_distrib_index_get (gre->ptb);

    gre->btp = PDM_block_to_part_create (distrib,
                                         (const PDM_g_num_t **) gre->g_nums,
                                         gre->n_elts,
                                         gre->n_part,
                                         gre->comm);
  }

  int    *block_field_stride  = NULL;
  double *block_field         = NULL;

  for (int i = 0; i < gre->n_part; i++) {
    for (int j = 0; j < gre->n_elts[i]; j++) {
      gre->strides[i][j] = gre->stride;
    }
  }

  PDM_part_to_block_exch (gre->ptb,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          gre->strides,
                          (void **) gre->local_field,
                          &block_field_stride,
                          (void **) &block_field);

  int n_elt_block = PDM_part_to_block_n_elt_block_get(gre->ptb);

  int *stride_idx = malloc(sizeof(int) * (n_elt_block + 1));
  stride_idx[0] = 0;

  for (int i = 0; i < n_elt_block; i++) {
    stride_idx[i+1] = stride_idx[i] + block_field_stride[i]/gre->stride;
  }
  //free (block_field_stride);


  double (*_reduce) (const double, const double) = NULL;

  switch (gre->operation) {
    case PDM_REDUCE_OP_MIN:
    _reduce = &_min_a_b;
    break;

    case PDM_REDUCE_OP_MAX:
    _reduce = &_max_a_b;
    break;

    case PDM_REDUCE_OP_SUM:
    _reduce = &_sum_a_b;
    break;

    default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid reduction operation\n");
  }


  if (1) {//gre->operation == PDM_REDUCE_OP_SUM) {
    int idx = 0;
    for (int i = 0; i < n_elt_block; i++) {
      for (int j = 0; j < gre->stride; j++) {
        block_field[gre->stride*i + j] = block_field[gre->stride*idx + j];
      }
      idx++;

      for (int k = 1; k < block_field_stride[i]/gre->stride; k++) {
        for (int j = 0; j < gre->stride; j++) {
          //block_field[gre->stride*i + j] += block_field[gre->stride*idx + j];
          block_field[gre->stride*i + j] = _reduce(block_field[gre->stride*i + j],
                                                   block_field[gre->stride*idx + j]);
        }
        idx++;
      }
    }

  }
  else {
    for (int i = 0; i < n_elt_block; i++) {
      for (int j = stride_idx[i]+1; j < stride_idx[i+1]; j++) {
        for (int k = 0; k < gre->stride; k++) {
          double val = block_field[j*gre->stride + k];
          block_field[i*gre->stride + k] = _reduce(block_field[i*gre->stride + k],
                                                   val);
        }
      }
    }
  }
  free (stride_idx);
  free (block_field_stride);


  PDM_block_to_part_exch_in_place (gre->btp,
                          sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          &gre->stride,
                          block_field,
                          NULL,
                          (void **) gre->global_reduced_field);
  free (block_field);

  for (int i = 0; i < gre->n_part; i++) {
    gre->local_field         [i] = NULL;
    gre->global_reduced_field[i] = NULL;
  }
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
