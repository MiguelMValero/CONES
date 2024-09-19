#ifndef __PDM_writer_PART_TO_BLOCK_PRIV_H__
#define __PDM_writer_PART_TO_BLOCK_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

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
 * \struct _pdm_part_to_block_t
 *
 * \brief  Data transfer from partitions to blocks
 *
 */

struct _pdm_part_to_block_t {

  /*
   * Block distribution properties
   */

  PDM_part_to_block_distrib_t  t_distrib;          /*!< Distribution type */
  PDM_part_to_block_post_t     t_post;             /*!< post processing type */
  int                          n_active_ranks;     /*!< Number of active ranks */
  int                         *active_ranks;       /*!< List of active ranks */
  PDM_MPI_Comm                 comm;               /*!< MSG communicator */
  int                          s_comm;             /*!< Communicator size */
  int                          i_rank;             /*!< Current rank in comm */
  int                          is_my_rank_active;  /*!< Is active current rank */
  double                       part_active_node;   /*!< Part of active nodes */

  /*
   * General exchange data
   */

  int                         n_part;               /*!< Number of parts */
  int                        *n_elt;                /*!< Number of elements for any part */
  int                         n_elt_proc;           /*!< Number of elements on the current processus */
  double                    **weight;               /*!< Weight of elements */
  PDM_g_num_t               **gnum_elt;             /*!< Global numbering of elements for any part */
  int                        *dest_proc;            /*!< Destination process for any element (size = n_elt_proc) */

  PDM_g_num_t                *data_distrib_index;  /*!< Data distribution on ranks
                                                         (size = s_comm + 1) */
  int                         s_block_min;          /*!< Minimum block size */
  int                         s_block_max;          /*!< Maximum block size */

  int                        *i_send_data;          /*!< Data to send to other processes index
                                                     (size = s_comm) */
  int                        *i_recv_data;          /*!< Received Data from other processes index
                                                     (size = s_comm) */
  int                        *n_send_data;          /*!< Number of data to send to other processes
                                                     (size = s_comm) */
  int                        *n_recv_data;          /*!< Number of received Data from other processes
                                                     (size = s_comm) */

  int                         tn_send_data;         /*!< Total number of sended data */
  int                         tn_recv_data;         /*!< Total number of received data */
  PDM_g_num_t                *sorted_recv_gnum;     /*!< Sorted Global number of
                                                        reveived data (size = tn_recvData) */
  int                        *order;                /*!< Order of sorted_recvGnum
                                                      (size = tn_recvData) */
  int                         n_elt_block ;         /*!< Number of element in current block */
  PDM_g_num_t                *block_gnum;           /*!< Sorted Global number of
                                                         reveived data (size = block_n_elt) */
  int                        *block_gnum_count;     /*!< Number of occurence of each gnum in partitions
                                                         (size = block_n_elt) */

  double                    **weight_g;             /*!< Global weights of elements for any part */

  int                         enable_reverse;
  int                        *idx_partial;          /*! Index of gnum in partial block */

  /* Asynchrone */
  int                      max_exch_request;
  int                      next_request;
  size_t                  *s_data;
  PDM_stride_t*            t_stride;
  int                     *cst_stride;
  int                     *wait_status; /* 0 - Send / 1 - Is recv / 2 - Is treated */
  PDM_MPI_Request         *request_mpi;
  unsigned char          **send_buffer;
  unsigned char          **recv_buffer;
  int                    **recv_stride;
  int                    **n_send_buffer; // TODO -> size_t if long
  int                    **i_send_buffer;
  int                    **n_recv_buffer;
  int                    **i_recv_buffer;

  int                   ***block_stride;
  void                  ***block_data;

  int                  ****part_stride;
  void                 ****part_data;

  PDM_mpi_comm_kind_t     *comm_kind;
  PDM_MPI_Win             *win_send;
  PDM_MPI_Win             *win_recv;

} ;

/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_writer_PART_TO_BLOCK_PRIV_H__ */
