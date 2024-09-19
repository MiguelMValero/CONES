/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

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

/*============================================================================
 * Look for the location of a global numbering element. The location has three
 * propertie : process, partition, number of element in this partition.
 * A global numbering can be located in multiple partitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Standard headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum_location.h"
#include "pdm_gnum_location_priv.h"
#include "pdm_priv.h"
#include "pdm_logging.h"

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
 * \brief Build a global numbering location structure
 *
 * \param [in]   n_part_in      Number of local partitions for elements
 * \param [in]   n_part_out     Number of local partitions for requested locations
 * \param [in]   comm           PDM_MPI communicator
 * \param [in]   owner          Owner
 *
 * \return     Pointer to \ref PDM_gnum_locaion object
 */

PDM_gnum_location_t*
PDM_gnum_location_create
(
 const int          n_part_in,
 const int          n_part_out,
 const PDM_MPI_Comm comm,
 const PDM_ownership_t owner
)
{

  PDM_gnum_location_t *gnum_loc = (PDM_gnum_location_t *) malloc(sizeof(PDM_gnum_location_t));

  gnum_loc->n_part_in  = n_part_in;
  gnum_loc->n_part_out = n_part_out;

  gnum_loc->n_elts_in = (      int          *) malloc (sizeof(      int          ) * n_part_in);
  gnum_loc->g_nums_in = (const PDM_g_num_t **) malloc (sizeof(const PDM_g_num_t *) * n_part_in);
  for (int i = 0; i < n_part_in; i++) {
    gnum_loc->g_nums_in[i] = NULL;
  }
  gnum_loc->n_elts_out = (      int          *) malloc (sizeof(      int          ) * n_part_out);
  gnum_loc->g_nums_out = (const PDM_g_num_t **) malloc (sizeof(const PDM_g_num_t *) * n_part_out);
  for (int i = 0; i < n_part_out; i++) {
    gnum_loc->g_nums_out[i] = NULL;
  }

  gnum_loc->location_idx = NULL;
  gnum_loc->location     = NULL;
  gnum_loc->comm         = comm;
  gnum_loc->tag_results_get = 0;
  gnum_loc->owner           = owner;

  return gnum_loc;
}



/**
 *
 * \brief Set global numbering
 *
 * \param [in]   gnum_loc    Pointer to \ref PDM_gnum_locaion object
 * \param [in]   i_part_in   Current partition
 * \param [in]   n_elts_in   Number of elements
 * \param [in]   gnum_in     Global numbering
 *
 */

void
PDM_gnum_location_elements_set
(
       PDM_gnum_location_t *gnum_loc,
 const int                  i_part_in,
 const int                  n_elts_in,
 const PDM_g_num_t         *gnum_in
)
{
  gnum_loc->n_elts_in[i_part_in] = n_elts_in;
  gnum_loc->g_nums_in[i_part_in] = gnum_in;
}



/**
 *
 * \brief Set requested elements
 *
 * \param [in]   gnum_loc     Pointer to \ref PDM_gnum_locaion object
 * \param [in]   i_part_out   Current partition
 * \param [in]   n_elts_out   Number of elements
 * \param [in]   gnum_out     Global numbering
 *
 */

void
PDM_gnum_location_requested_elements_set
(
       PDM_gnum_location_t *gnum_loc,
 const int                  i_part_out,
 const int                  n_elts_out,
 const PDM_g_num_t         *gnum_out
)
{
  gnum_loc->n_elts_out[i_part_out] = n_elts_out;
  gnum_loc->g_nums_out[i_part_out] = gnum_out;
}

/**
 *
 * \brief Compute the location (processus, partittion, local number in the partition)
 *
 * \param [in]   gnum_loc     Pointer to \ref PDM_gnum_locaion object
 *
 */

void
PDM_gnum_location_compute
(
 PDM_gnum_location_t *gnum_loc
)
{
  // PDM_MPI_Barrier(gnum_loc->comm);
  // double t1 = PDM_MPI_Wtime();

  int rank;
  PDM_MPI_Comm_rank (gnum_loc->comm, &rank);

  int n_rank;
  PDM_MPI_Comm_size (gnum_loc->comm, &n_rank);

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1,
                                                       (PDM_g_num_t **) gnum_loc->g_nums_in,
                                                       NULL,
                                                       gnum_loc->n_elts_in,
                                                       gnum_loc->n_part_in,
                                                       gnum_loc->comm);

  // PDM_g_num_t *block_distrib_index = PDM_part_to_block_distrib_index_get (ptb);

  const int s_data = sizeof(int);
  const PDM_stride_t t_stride = PDM_STRIDE_VAR_INTERLACED;
  const int cst_stride = 3;

  int  **part_stride = (int **) malloc (sizeof(int *) * gnum_loc->n_part_in);
  int  **part_data = (int **) malloc (sizeof(int *) * gnum_loc->n_part_in);

  for (int i = 0; i < gnum_loc->n_part_in; i++) {
    part_stride[i] = malloc (sizeof(int) * gnum_loc->n_elts_in[i]);
    part_data[i] = malloc (sizeof(int) * 3 * gnum_loc->n_elts_in[i]);
    for (int j = 0; j < gnum_loc->n_elts_in[i]; j++) {
      part_stride[i][j]   = 3;
      part_data[i][3*j]   = rank;
      part_data[i][3*j+1] = i;
      part_data[i][3*j+2] = j;
    }
  }

  int  *block_stride = NULL;
  int  *block_data = NULL;

  PDM_part_to_block_exch (ptb,
                          s_data,
                          t_stride,
                          cst_stride,
                          part_stride,
                          (void **) part_data,
                          &block_stride,
                          (void **) &block_data);

  // int request_id = -1;
  // PDM_part_to_block_iexch (ptb,
  //                          PDM_MPI_COMM_KIND_COLLECTIVE,
  //                          s_data,
  //                          t_stride,
  //                          cst_stride,
  //                          part_stride,
  //                (void **) part_data,
  //                          &block_stride,
  //                (void **) &block_data,
  //                          &request_id);

  for (int i = 0; i < gnum_loc->n_part_in; i++) {
    free (part_stride[i]);
    free (part_data[i]);
  }
  free (part_data);
  free (part_stride);


  PDM_g_num_t* block_g_num = PDM_part_to_block_block_gnum_get (ptb);
  int          dn_block    = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(block_g_num,
                                                                        dn_block,
                                                                        gnum_loc->g_nums_out,
                                                                        gnum_loc->n_elts_out,
                                                                        gnum_loc->n_part_out,
                                                                        gnum_loc->comm);

  // PDM_part_to_block_iexch_wait(ptb, request_id);

  PDM_block_to_part_exch (btp,
                          s_data,
                          PDM_STRIDE_VAR_INTERLACED,
                          block_stride,
                          block_data,
                          &part_stride,
               (void ***) &gnum_loc->location);

  gnum_loc->location_idx = (int **) malloc (sizeof(int *) * gnum_loc->n_part_out);
  for (int i = 0; i < gnum_loc->n_part_out; i++) {
    gnum_loc->location_idx[i] = malloc (sizeof(int) * (gnum_loc->n_elts_out[i] + 1));
    gnum_loc->location_idx[i][0] = 0;
    for (int j = 0; j < gnum_loc->n_elts_out[i]; j++) {
      gnum_loc->location_idx[i][j+1] = gnum_loc->location_idx[i][j] + part_stride[i][j];
    }
  }
  free (block_stride);
  free (block_data);
  // free (block_distrib_index_correct);

  for (int i = 0; i < gnum_loc->n_part_out; i++) {
    free (part_stride[i]);
  }

  free (part_stride);

  PDM_part_to_block_free (ptb);
  PDM_block_to_part_free (btp);

  // PDM_MPI_Barrier(gnum_loc->comm);
  // double dt = PDM_MPI_Wtime()-t1;
  // if(rank == 0) {
  //   printf("Time gnum_location compute : %12.5e \n", dt);
  // }

}


/**
 *
 * \brief Get location
 *
 * \param [in]    gnum_loc       Pointer to \ref PDM_gnum_locaion object
 * \param [in]    i_part_out     Current partition
 * \param [out]   location_idx   Index in the location arrays (size = \ref n_elts + 1)
 * \param [out]   location       Locations of each element
 *                                (Three informations : process, partition, element)
 *
 */

void
PDM_gnum_location_get
(
       PDM_gnum_location_t  *gnum_loc,
 const int                   i_part_out,
       int                 **location_idx,
       int                 **location
)
{
  *location_idx = gnum_loc->location_idx[i_part_out];
  *location     = gnum_loc->location    [i_part_out];
  gnum_loc->tag_results_get = 1;
}


/**
 *
 * \brief Free
 *
 * \param [in]   gnum_loc      Pointer to \ref PDM_gnum_locaion object
 *
 */

void
PDM_gnum_location_free
(
  PDM_gnum_location_t *gnum_loc
)
{

  free (gnum_loc->n_elts_in);
  free (gnum_loc->g_nums_in);

  free (gnum_loc->n_elts_out);
  free (gnum_loc->g_nums_out);

  if(( gnum_loc->owner == PDM_OWNERSHIP_KEEP ) ||
     ( gnum_loc->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !gnum_loc->tag_results_get)) {
    for (int i = 0; i < gnum_loc->n_part_out; i++) {
      free (gnum_loc->location_idx[i]);
    }
    // free (gnum_loc->location_idx);

    for (int i = 0; i < gnum_loc->n_part_out; i++) {
      free (gnum_loc->location[i]);
    }
    // free (gnum_loc->location);
  }

  free (gnum_loc->location_idx);
  free (gnum_loc->location);

  free (gnum_loc);

}


/**
 *
 * \brief Get the number of requested elements in a given partition
 *
 * \param [in]  gnum_loc      Pointer to \ref PDM_gnum_locaion object
 * \param [in]  i_part_out    Id of current partition
 *
 * \return  Number of requested elements in the current partition
 *
 */

int
PDM_gnum_location_n_requested_elt_get
(
       PDM_gnum_location_t *gnum_loc,
 const int                  i_part_out
 )
{
  assert (gnum_loc != NULL);
  assert (i_part_out < gnum_loc->n_part_out);

  return gnum_loc->n_elts_out[i_part_out];
}



/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
