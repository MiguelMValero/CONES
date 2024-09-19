/*
 * \file
 */

#ifndef __PDM_PART_TO_BLOCK_H__
#define __PDM_PART_TO_BLOCK_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_geom.h"

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
 * \enum PDM_part_to_block_distrib_t
 * \brief Type of block distribution
 *
 */

typedef enum {

  PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC          = 0,  /*!< Distribute block on all processes */
  PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE = 1,  /*!< Distribute block on a single process per node */
  PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE      = 2   /*!< Distribute block on part of nodes */

} PDM_part_to_block_distrib_t;


/**
 * \enum PDM_part_to_block_post_t
 * \brief Type of post-processing performed on blocks
 *
 */

typedef enum {

  PDM_PART_TO_BLOCK_POST_NOTHING       = 0,  /*!< No post processing                 */
  PDM_PART_TO_BLOCK_POST_CLEANUP       = 1,  /*!< Cleanup multi-elements             */
  PDM_PART_TO_BLOCK_POST_MERGE         = 2,  /*!< Merge multi-elements               */
  PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM = 3   /*!< Merge but keep the original buffer */

} PDM_part_to_block_post_t;


/**
 * \struct PDM_part_to_block_t
 * \brief  Partition-to-Block redistribution
 *
 */

typedef struct _pdm_part_to_block_t PDM_part_to_block_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief PDM_extents_conformize
 *
 * Correction extents to manage singular cases and breaks symmetry
 * eps = 1.e-3 is a standard value
 *
 * \param [in]    dim      Spatial dimension
 * \param [inout] extents  Spatial extents (size = 2 * \p dim : [*x_min*, *y_min*, ..., *x_max*, *y_max*, ...])
 * \param [in]    eps      Geometric tolerance
 *
 */
void
PDM_extents_conformize(int    dim,
                       double extents[],
                       double eps);


/**
 *
 * \brief Reset global part-to-block statistics
 *
 */

void
PDM_part_to_block_global_statistic_reset
(
void
);


/**
 *
 * \brief Get global part-to-block statistics
 *
 * \param [in]   comm                 MPI communicator
 * \param [out]  min_exch_rank_send   Global min part of ranks used to send
 * \param [out]  min_exch_rank_recv   Global min part of ranks used to receive
 * \param [out]  max_exch_rank_send   Global max part of ranks used to send
 * \param [out]  max_exch_rank_recv   Global max part of ranks used to receive
 * \param [out]  min_exch_data_send   Global min sent data for a rank
 * \param [out]  min_exch_data_recv   Global min received data for a rank
 * \param [out]  max_exch_data_send   Global max sent data for a rank
 * \param [out]  max_exch_data_recv   Global max received data for a rank
 *
 */

void
PDM_part_to_block_global_statistic_get
(
 PDM_MPI_Comm comm,
 int *min_exch_rank_send,
 int *min_exch_rank_recv,
 int *max_exch_rank_send,
 int *max_exch_rank_recv,
 unsigned long long *min_exch_data_send,
 unsigned long long *min_exch_data_recv,
 unsigned long long *max_exch_data_send,
 unsigned long long *max_exch_data_recv
);


/**
 *
 * \brief Get global part-to-block timer
 *
 * \param [in]   comm              MPI communicator
 * \param [out]  min_elaps         Min elapsed time
 * \param [out]  max_elaps         Max elapsed time
 * \param [out]  min_cpu           Min cpu time
 * \param [out]  max_cpu           Max cpu time
 * \param [out]  min_elaps_create  Global min elapsed for create function
 * \param [out]  max_elaps_create  Global max elapsed for create function
 * \param [out]  min_cpu_create    Global min cpu for create function
 * \param [out]  max_cpu_create    Global max cpu for create function
 * \param [out]  min_elaps_create2 Global min elapsed for create2 function
 * \param [out]  max_elaps_create2 Global max elapsed for create2 function
 * \param [out]  min_cpu_create2   Global min cpu for create2 function
 * \param [out]  max_cpu_create2   Global max cpu for create2 function
 * \param [out]  min_elaps_exch    Global min elapsed for exch function
 * \param [out]  max_elaps_exch    Global max elapsed for exch function
 * \param [out]  min_cpu_exch      Global min cpu for exch function
 * \param [out]  max_cpu_exch      Global max cpu for exch function
 *
 */

void
PDM_part_to_block_global_timer_get
(
 PDM_MPI_Comm comm,
 double       *min_elaps_create,
 double       *max_elaps_create,
 double       *min_cpu_create,
 double       *max_cpu_create,
 double       *min_elaps_create2,
 double       *max_elaps_create2,
 double       *min_cpu_create2,
 double       *max_cpu_create2,
 double       *min_elaps_exch,
 double       *max_elaps_exch,
 double       *min_cpu_exch,
 double       *max_cpu_exch
);

/**
 *
 * \brief Create a part-to-block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   partActiveNode  Part of active nodes (\ref PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global numbers
 * \param [in]   weight          Element weights of elements (or ``NULL``)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partitions
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_block_t object
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        partActiveNode,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
);


/**
 *
 * \brief Create a part-to-block redistribution from a given distribution index
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   partActiveNode  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global numbers
 * \param [in]   distrib         Distribution index (\p distrib[0] = 0 and size = *n_rank* + 1)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partitions
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_block_t object
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create_from_distrib
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        partActiveNode,
 PDM_g_num_t                 **gnum_elt,
 const PDM_g_num_t            *distrib,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
);


/**
 *
 * \brief Create a part-to-block redistribution from geometric renumbering
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   partActiveNode  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   geom_renum      Geometric renumbering
 * \param [in]   gnum_elt        Element global numbers
 * \param [in]   weight          Element weights of elements (or ``NULL``)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partitions
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_block_t object
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_geom_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_part_geom_t               geom_renum,
 double                      **coords,
 PDM_g_num_t                 **gnum,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
);

/**
 *
 * \brief Return number of active ranks
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return Number of active ranks
 *
 */

int
PDM_part_to_block_n_active_ranks_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return active ranks
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return  List of active ranks
 *
 */

int *
PDM_part_to_block_active_ranks_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return if current rank is active
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return  1 if current rank is active, 0 otherwise
 *
 */

int
PDM_part_to_block_is_active_rank
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return number of elements on current process (in the block-distributed frame)
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return Number of elements in the current process
 *
 */

int
PDM_part_to_block_n_elt_block_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return global numbers of elements on current process (in the block-distributed frame)
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return  Global numbers
 *
 */

PDM_g_num_t *
PDM_part_to_block_block_gnum_get
(
 PDM_part_to_block_t *ptb
);

/**
 *
 * \brief Return numbers of occurrence of each gnum element on current process
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return  Global numbers counter
 *
 */

int *
PDM_part_to_block_block_gnum_count_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Exchange data from *partitioned* frame to *block-distributed* frame
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride (only for \ref PDM_STRIDE_CST_INTERLACED or \ref PDM_STRIDE_CST_INTERLEAVED)
 * \param [in]   part_stride  Variable stride (size = \p n_part) (only for \ref PDM_STRIDE_VAR_INTERLACED)
 * \param [in]   part_data    Partitioned data (size = \p n_part)
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 * \return       Size of \p block_data array
 *
 */

int
PDM_part_to_block_exch
(
 PDM_part_to_block_t       *ptb,
 size_t                     s_data,
 PDM_stride_t               t_stride,
 int                        cst_stride,
 int                      **part_stride,
 void                     **part_data,
 int                      **block_stride,
 void                     **block_data
);


/**
 *
 * \brief Exchange data from *block-distributed* frame to *partitioned* frame
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride (only for \ref PDM_STRIDE_CST_INTERLACED or \ref PDM_STRIDE_CST_INTERLEAVED)
 * \param [in]   block_stride Block stride
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Variable stride (size = \p n_part) (only for \ref PDM_STRIDE_VAR_INTERLACED)
 * \param [out]  part_data    Partitioned data (size = \p n_part)
 *
 */

void
PDM_part_to_block_reverse_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                 *block_stride,
 void                *block_data,
 int               ***part_stride,
 void              ***part_data
);


/**
 *
 * \brief Initiate a data exchange
 * (TODO: Replace by \ref PDM_part_to_block_iexch ?)
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 * \param [in]   part_stride  Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   part_data    partitioned data
 *
 * \return       Exchange request
 *
 */

int
PDM_part_to_block_async_exch
(
 PDM_part_to_block_t       *ptb,
 size_t                     s_data,
 PDM_stride_t               t_stride,
 int                        cst_stride,
 int                      **part_stride,
 void                     **part_data
);

/**
 *
 * \brief Initiate a non-blocking data exchange from *partitioned* frame
 *        to *block-distributed* frame
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   k_comm       Kind of MPI communication
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride (only for \ref PDM_STRIDE_CST_INTERLACED or \ref PDM_STRIDE_CST_INTERLEAVED)
 * \param [in]   part_stride  Variable stride (size = \p n_part) (only for \ref PDM_STRIDE_VAR_INTERLACED)
 * \param [in]   part_data    Partitioned data (size = \p n_part)
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 * \param [out]  request      Request
 *
 */

void
PDM_part_to_block_iexch
(
       PDM_part_to_block_t  *ptb,
 const PDM_mpi_comm_kind_t   k_comm,
       size_t                s_data,
       PDM_stride_t          t_stride,
       int                   cst_stride,
       int                 **part_stride,
       void                **part_data,
       int                 **block_stride,
       void                **block_data,
       int                  *request
);

/**
 *
 * \brief Finalize and post-process a non-blocking data exchange
 *        from *partitioned* frame to *block-distributed* frame
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   request_id   Request to wait / post-process
 *
 * \return       Size of block data
 *
 */

int
PDM_part_to_block_iexch_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
);


/**
 *
 * \brief Initiate a non-blocking data exchange from *block-distributed* frame
 *        to *partitioned* frame
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   k_comm       Kind of MPI communication
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride (only for \ref PDM_STRIDE_CST_INTERLACED or \ref PDM_STRIDE_CST_INTERLEAVED)
 * \param [in]   block_stride Block stride
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Variable stride (size = \p n_part) (only for \ref PDM_STRIDE_VAR_INTERLACED)
 * \param [out]  part_data    Partitioned data (size = \p n_part)
 * \param [out]  request      Request
 *
 */
void
PDM_part_to_block_reverse_iexch
(
       PDM_part_to_block_t  *ptb,
 const PDM_mpi_comm_kind_t   k_comm,
       size_t                s_data,
       PDM_stride_t          t_stride,
       int                   cst_stride,
       int                  *block_stride,
       void                 *block_data,
       int                ***part_stride,
       void               ***part_data,
       int                  *request
);

/**
 *
 * \brief Finalize a non-blocking data exchange from *block-distributed* frame to *partitioned* frame
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   request_id   Request to wait
 *
 */

void
PDM_part_to_block_reverse_iexch_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
);

/**
 *
 * \brief Wait for an exchange
 * (TODO: Replace by \ref PDM_part_to_block_iexch_wait ?)
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 */
void
PDM_part_to_block_async_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
);

/**
 *
 * \brief Get the raw exchange buffer and stride and free memory
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 * \return Size of raw received data array
 */
int
PDM_part_to_block_asyn_get_raw
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
);

/**
 *
 * \brief Get the raw exchange buffer and stride and free memory
 *
 * \param [in]   ptb          Part-to-Block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 * \return       Size of block data
 */
int
PDM_part_to_block_asyn_post_treatment
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
);

/**
 *
 * \brief Free a Part-to-Block structure
 *
 * \param [inout] ptb         Part-to-Block structure
 *
 * \return       Null pointer
 */

PDM_part_to_block_t *
PDM_part_to_block_free
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Return block distribution
 *
 * \param [in] ptb         Part-to-Block structure
 *
 * \return  Distribution (size = *n_rank* + 1)
 */

PDM_g_num_t *
PDM_part_to_block_distrib_index_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Return destination process
 *
 * \param [in] ptb         Part-to-Block structure
 *
 * \return  Destination (size = sum of partition elements)
 */

PDM_l_num_t *
PDM_part_to_block_destination_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Adapt a partial block (padd stride and distribution arrays)
 *
 * \param [in]    ptb         Part-to-Block structure
 * \param [inout] block_n     Block stride
 * \param [in]    n_g_block   Total number of elements (across all processes)
 *
 * \return New block distribution (size = *n_rank* + 1)
 *
 */

PDM_g_num_t *
PDM_part_to_block_adapt_partial_block_to_block
(
 PDM_part_to_block_t  *ptb,
 int                 **block_n,
 PDM_g_num_t           n_g_block
);



/**
 *
 * \brief Return global weights of elements on current process
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return  Global weights (size = *n_part*)
 *
 */

double **
PDM_part_to_block_global_weight_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Get number of MPI ranks
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return  Number of MPI ranks
 *
 */

int
PDM_part_to_block_n_ranks_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Return total number of elements on current process (summed over all partitions)
 *
 * \param [in]   ptb          Part-to-Block structure
 *
 * \return Total number of elements on current process (in *partition* frame)
 *
 */

int
PDM_part_to_block_n_elt_proc_get
(
 PDM_part_to_block_t *ptb
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_writer_PART_TO_BLOCK_H__ */
