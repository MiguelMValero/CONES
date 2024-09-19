/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_block_to_block.h"
#include "pdm_block_to_block_priv.h"
#include "pdm_binary_search.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a block to partitions redistribution
 *
 * \param [in]   blockDistribIdx Block distribution (size : \ref size of \ref comm + 1)
 *                               C numbering (blockDistribIdx[0] = 0)
 * \param [in]   gnum_elt        Element global number (size : \ref n_part)
 * \param [in]   n_elt           Local number of elements (size : \ref n_part)
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_block instance
 *
 */

PDM_block_to_block_t *
PDM_block_to_block_create
(
 PDM_g_num_t   *block_distrib_ini_idx,
 PDM_g_num_t   *block_distrib_end_idx,
 PDM_MPI_Comm   comm
)
{

  _pdm_block_to_block_t *btb =
    (_pdm_block_to_block_t *) malloc (sizeof(_pdm_block_to_block_t));

  btb->comm = comm;
  PDM_MPI_Comm_size (comm, &btb->n_rank);
  PDM_MPI_Comm_rank (comm, &btb->i_rank);

  /*
   * Define requested data for each process
   */

  btb->block_distrib_ini_idx = malloc (sizeof(PDM_g_num_t) * (btb->n_rank + 1));
  for (int i = 0; i < btb->n_rank + 1; i++) {
    btb->block_distrib_ini_idx[i] = block_distrib_ini_idx[i];
  }

  btb->block_distrib_end_idx = malloc (sizeof(PDM_g_num_t) * (btb->n_rank + 1));
  for (int i = 0; i < btb->n_rank + 1; i++) {
    btb->block_distrib_end_idx[i] = block_distrib_end_idx[i];
  }

  /*
   * Verbose
   */

  if(0 == 1){

    PDM_printf("block_distrib_ini_idx : ");
    for(int i = 0; i < btb->n_rank+1; i++){
      PDM_printf("%i ", btb->block_distrib_ini_idx[i]);
    }
    PDM_printf("\n");

    PDM_printf("block_distrib_end_idx : ");
    for(int i = 0; i < btb->n_rank+1; i++){
      PDM_printf("%i ", btb->block_distrib_end_idx[i]);
    }
    PDM_printf("\n");
  }

  return (PDM_block_to_block_t *) btb;

}


/**
 *
 * \brief Initialize an exchange (Variable stride is not yet available)
 *
 * \param [in]   btb          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Constant stride
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride or NULL
 * \param [out]  part_data    Partition data
 *
 */

int
PDM_block_to_block_exch
(
 PDM_block_to_block_t *btb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                 *block_stride_ini,
 void                *block_data_ini,
 int                 *block_stride_end,
 void                **block_data_end
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;

  int *i_send_buffer = (int *) malloc (sizeof(int) * _btb->n_rank);
  int *i_recv_buffer = (int *) malloc (sizeof(int) * _btb->n_rank);
  int *n_send_buffer = (int *) malloc (sizeof(int) * _btb->n_rank);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * _btb->n_rank);

  for (int i = 0; i < _btb->n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  unsigned char *sendBuffer = (unsigned char *) block_data_ini;

  int s_recv_buffer = 0;

  /*
   * Exchange Stride and build buffer properties
   */

  for (PDM_g_num_t i = _btb->block_distrib_ini_idx[_btb->i_rank]; i < _btb->block_distrib_ini_idx[_btb->i_rank+1]; i++) {

    int send_rank = PDM_binary_search_gap_long (i,
                                               _btb->block_distrib_end_idx,
                                               _btb->n_rank + 1);
    n_send_buffer[send_rank] += 1;
  }

  for (PDM_g_num_t i = _btb->block_distrib_end_idx[_btb->i_rank]; i < _btb->block_distrib_end_idx[_btb->i_rank+1]; i++) {

    int recv_rank = PDM_binary_search_gap_long (i,
                                               _btb->block_distrib_ini_idx,
                                               _btb->n_rank + 1);
    n_recv_buffer[recv_rank] += 1;
  }

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    for(int i = 1; i < _btb->n_rank; i++){
      i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
    }

    PDM_MPI_Alltoallv(block_stride_ini,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_INT,
                      block_stride_end,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_INT,
                      _btb->comm);

    for (int i = 0; i < _btb->n_rank; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
      i_send_buffer[i] = 0;
      i_recv_buffer[i] = 0;
    }

    int k = 0;
    for (PDM_g_num_t i = _btb->block_distrib_ini_idx[_btb->i_rank]; i < _btb->block_distrib_ini_idx[_btb->i_rank+1]; i++) {

      int send_rank = PDM_binary_search_gap_long (i,
                                                 _btb->block_distrib_end_idx,
                                                 _btb->n_rank + 1);
      n_send_buffer[send_rank] += block_stride_ini[k++] * s_data;
    }

    k = 0;
    for (PDM_g_num_t i = _btb->block_distrib_end_idx[_btb->i_rank]; i < _btb->block_distrib_end_idx[_btb->i_rank+1]; i++) {

      int recv_rank = PDM_binary_search_gap_long (i,
                                                 _btb->block_distrib_ini_idx,
                                                 _btb->n_rank + 1);
      n_recv_buffer[recv_rank] += block_stride_end[k++] * s_data;
    }

    for(int i = 1; i < _btb->n_rank; i++){
      i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
    }

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    for(int i = 0; i < _btb->n_rank; i++){
      n_send_buffer[i] *= cst_stride * s_data;
      n_recv_buffer[i] *= cst_stride * s_data;
    }
    for(int i = 1; i < _btb->n_rank; i++){
      // i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1] * s_data * cst_stride;
      // i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1] * s_data * cst_stride;
      i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
    }

  }

  s_recv_buffer = i_recv_buffer[_btb->n_rank-1] + n_recv_buffer[_btb->n_rank-1];

  unsigned char *recv_buffer =
    (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

  /*
   * Data exchange
   */

  PDM_MPI_Alltoallv(sendBuffer,
                    n_send_buffer,
                    i_send_buffer,
                    PDM_MPI_BYTE,
                    recv_buffer,
                    n_recv_buffer,
                    i_recv_buffer,
                    PDM_MPI_BYTE,
                    _btb->comm);

  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  *block_data_end = recv_buffer;

  return s_recv_buffer/( (int)s_data );

}



int
PDM_block_to_block_exch_with_mpi_type
(
 PDM_block_to_block_t  *btb,
 PDM_stride_t           t_stride,
 PDM_MPI_Datatype       mpi_type,
 int                   *block_stride_ini,
 void                  *block_data_ini,
 int                   *block_stride_end,
 void                 **block_data_end
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;

  int *i_send_buffer = (int *) malloc (sizeof(int) * _btb->n_rank);
  int *i_recv_buffer = (int *) malloc (sizeof(int) * _btb->n_rank);
  int *n_send_buffer = (int *) malloc (sizeof(int) * _btb->n_rank);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * _btb->n_rank);

  for (int i = 0; i < _btb->n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  unsigned char *sendBuffer = (unsigned char *) block_data_ini;

  int s_recv_buffer = 0;

  /*
   * Exchange Stride and build buffer properties
   */

  for (PDM_g_num_t i = _btb->block_distrib_ini_idx[_btb->i_rank]; i < _btb->block_distrib_ini_idx[_btb->i_rank+1]; i++) {

    int send_rank = PDM_binary_search_gap_long (i,
                                               _btb->block_distrib_end_idx,
                                               _btb->n_rank + 1);
    n_send_buffer[send_rank] += 1;
  }

  for (PDM_g_num_t i = _btb->block_distrib_end_idx[_btb->i_rank]; i < _btb->block_distrib_end_idx[_btb->i_rank+1]; i++) {

    int recv_rank = PDM_binary_search_gap_long (i,
                                               _btb->block_distrib_ini_idx,
                                               _btb->n_rank + 1);
    n_recv_buffer[recv_rank] += 1;
  }

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    for(int i = 1; i < _btb->n_rank; i++){
      i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
    }

    PDM_MPI_Alltoallv(block_stride_ini,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_INT,
                      block_stride_end,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_INT,
                      _btb->comm);

    for (int i = 0; i < _btb->n_rank; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
      i_send_buffer[i] = 0;
      i_recv_buffer[i] = 0;
    }

    int k = 0;
    for (PDM_g_num_t i = _btb->block_distrib_ini_idx[_btb->i_rank]; i < _btb->block_distrib_ini_idx[_btb->i_rank+1]; i++) {

      int send_rank = PDM_binary_search_gap_long (i,
                                                 _btb->block_distrib_end_idx,
                                                 _btb->n_rank + 1);
      n_send_buffer[send_rank] += block_stride_ini[k++];
    }

    k = 0;
    for (PDM_g_num_t i = _btb->block_distrib_end_idx[_btb->i_rank]; i < _btb->block_distrib_end_idx[_btb->i_rank+1]; i++) {

      int recv_rank = PDM_binary_search_gap_long (i,
                                                 _btb->block_distrib_ini_idx,
                                                 _btb->n_rank + 1);
      n_recv_buffer[recv_rank] += block_stride_end[k++];
    }

    for(int i = 1; i < _btb->n_rank; i++){
      i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
    }

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {
    for(int i = 1; i < _btb->n_rank; i++){
      i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
    }

  }

  s_recv_buffer = i_recv_buffer[_btb->n_rank-1] + n_recv_buffer[_btb->n_rank-1];

  int s_type = 0;
  PDM_MPI_Type_size(mpi_type, &s_type);

  size_t tot_size = ((size_t) s_type ) * ((size_t) s_recv_buffer);
  unsigned char *recv_buffer = (unsigned char *) malloc(tot_size);

  /*
   * Data exchange
   */

  PDM_MPI_Alltoallv(sendBuffer,
                    n_send_buffer,
                    i_send_buffer,
                    mpi_type,
                    recv_buffer,
                    n_recv_buffer,
                    i_recv_buffer,
                    mpi_type,
                    _btb->comm);

  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  *block_data_end = recv_buffer;

  return s_recv_buffer/( (int)s_type );

}



/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] btb  Block to part structure
 *
 * \return       NULL
 */

PDM_block_to_block_t *
PDM_block_to_block_free
(
 PDM_block_to_block_t *btb
)
{
  _pdm_block_to_block_t *_btb = (_pdm_block_to_block_t *) btb;

  free (_btb->block_distrib_ini_idx);
  free (_btb->block_distrib_end_idx);

  free (_btb);

  return NULL;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
