/*
 * \file
 */

#ifndef __PDM_PART_TO_PART_H__
#define	__PDM_PART_TO_PART_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_part_to_part_t
 * \brief  Partition-to-Partition redistribution
 *
 */

typedef struct _pdm_part_to_part_t PDM_part_to_part_t;


/**
 * \enum PDM_part_to_part_data_def_t
 * \brief Kind of data definition
 *
 */

typedef enum {

  PDM_PART_TO_PART_DATA_DEF_ORDER_PART1           = 0, /*!< Data defined according
                                                         to the part1 arrays order*/
  PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2  = 1, /*!< Data defined according
                                                         to the part1_to_part2 arrays order */
  PDM_PART_TO_PART_DATA_DEF_ORDER_PART2           = 2, /*!< Data defined according
                                                         to the part2 arrays order*/
  PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM = 3, /*!< Data defined according
                                                         to the gnum1_come_from arrays order */

} PDM_part_to_part_data_def_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a Partition-to-Partition redistribution from global ids
 *
 * \param [in]   gnum_elt1          Element global ids in Part1 (size : \p n_part1)
 * \param [in]   n_elt1             Local number of elements in Part1 (size : \p n_part1)
 * \param [in]   n_part1            Number of partitions in Part1
 * \param [in]   gnum_elt2          Element global ids in Part2 (size : \p n_part2)
 * \param [in]   n_elt2             Local number of elements in Part2 (size : \p n_part2)
 * \param [in]   n_part2            Number of partitions in Part2
 * \param [in]   part1_to_part2_idx Index for Part1→Part2 mapping <br>
 *                                  (for each part, size : \p n_elt1 + 1)
 * \param [in]   part1_to_part2     Part1→Part2 mapping (global ids) <br>
 *                                  (for each part, size : \p part1_to_part2_idx[\p n_elt1])
 * \param [in]   comm               MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_part instance
 *
 */

PDM_part_to_part_t *
PDM_part_to_part_create
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t   **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **part1_to_part2_idx,
 const PDM_g_num_t   **part1_to_part2,
 const PDM_MPI_Comm    comm
);


/**
 *
 * \brief Create a Partition-to-Partition redistribution from location triplets
 *
 * \param [in]   gnum_elt1                   Element global numbers in Part1 (size : \p n_part1)
 * \param [in]   n_elt1                      Local number of elements in Part1 (size : \p n_part1)
 * \param [in]   n_part1                     Number of partitions in Part1
 * \param [in]   n_elt2                      Local number of elements in Part2 (size : \p n_part2)
 * \param [in]   n_part2                     Number of partitions in Part2
 * \param [in]   part1_to_part2_idx          Index for Part1→Part2 mapping <br>
 *                                           (for each part, size : \p n_elt1 + 1)
 * \param [in]   part1_to_part2_triplet_idx  Index for multiple locations in Part2 <br>
 *                                           (for each part, size : \p part1_to_part2_idx[\p n_elt1] + 1)
 * \param [in]   part1_to_part2_triplet      Part1→Part2 mapping (location triplets: (*irank2*, *ipart2*, *ielt2*)) <br>
 *                                           (for each part, size : \p part1_to_part2_triplet_idx[\p part1_to_part2_idx[\p n_elt1]] + 1)
 * \param [in]   comm                        MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_part instance
 *
 */

PDM_part_to_part_t *
PDM_part_to_part_create_from_num2_triplet
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const int            *n_elt2,
 const int             n_part2,
 const int           **part1_to_part2_idx,
 const int           **part1_to_part2_triplet_idx,
 const int           **part1_to_part2_triplet,
 const PDM_MPI_Comm    comm
);


/**
 *
 * \brief Initiate a collective Part1→Part2 exchange (based on `MPI_ialltoall`)
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Part1 data to send (ordered as in *part1_to_part2*)
 * \param [out]  ref_part2_data      Part2 data to receive (only referenced Part2 elements)
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_ialltoall
(
PDM_part_to_part_t *ptp,
 const size_t       s_data,
 const int          cst_stride,
 void             **part1_to_part2_data,
 void             **ref_part2_data,
 int               *request
);



/**
 *
 * \brief Initiate a collective Part2→Part1 exchange (based on `MPI_ialltoall`)
 *
 * \param [in]   ptp                      Part-to-Part structure
 * \param [in]   s_data                   Data size
 * \param [in]   cst_stride               Constant stride
 * \param [in]   ref_part2_to_part1_data  Part2 data to send (ordered as *ref_lnum2*)
 * \param [in]   part1_to_part2_data      Part1 data to receive (ordered as *part1_to_part2*)
 * \param [out]  request                  Request
 *
 */

void
PDM_part_to_part_reserve_ialltoall
(
PDM_part_to_part_t *ptp,
 const size_t       s_data,
 const int          cst_stride,
 void             **ref_part2_to_part1_data,
 void             **part1_part2_data,
 int               *request
);


/**
 *
 * \brief Finalize a collective Part1→Part2 exchange
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_ialltoall_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
);


/**
 *
 * \brief Finalize a collective Part2→Part1 exchange
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_ialltoall_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
);


/**
 *
 * \brief Initiate a collective Part1→Part2 exchange (based on `MPI_ineighbor_alltoall`)
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Part1 data to send (ordered as \ref part1_to_part2)
 * \param [out]  ref_part2_data      Part2 data to receive (only referenced Part2 elements)
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_ineighbor_alltoall
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part1_to_part2_data,
 void              **ref_part2_data,
 int                *request
);


/**
 *
 * \brief Finalize a collective Part1→Part2 exchange
 *
 * \param [in]  ptp                 Part-to-Part structure
 * \param [in]  request             Request
 *
 */


void
PDM_part_to_part_ineighbor_alltoall_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
);

/**
 *
 * \brief Get selected numbers of part2 index
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  n_elt1              Number of gnum1 element
 * \param [out]  part1_to_part2_idx  Index of data to send to gnum2 from gnum1
 *                                  (for each part size : \ref n_elt1+1)
 */

void
PDM_part_to_part_part1_to_part2_idx_get
(
 PDM_part_to_part_t *ptp,
 int               **n_elt1,
 int              ***part1_to_part2_idx
);


/**
 *
 * \brief Get Part1→Part2 mapping (global ids)
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [out]  n_elt1              Number of Part1 elements (size = *n_part1*)
 * \param [out]  part1_to_part2_idx  Index for Part1→Part2 mapping <br>
 *                                   (for each part, size : \p n_elt1 + 1)
 * \param [out]  part1_to_part2      Part1→Part2 mapping (global ids) <br>
 *                                   (for each part, size : \p part1_to_part2_idx[\p n_elt1] + 1)
 */

void
PDM_part_to_part_part1_to_part2_get
(
 PDM_part_to_part_t *ptp,
 int               **n_elt1,
 int              ***part1_to_part2_idx,
 PDM_g_num_t      ***part1_to_part2
);



/**
 *
 * \brief Get referenced Part2 elements
 *
 * \param [in]   ptp           Part-to-Part structure
 * \param [out]  n_ref_lnum2   Number of referenced Part2 elements
 * \param [out]  ref_lnum2     Referenced Part2 elements (1-based local ids)
 *
 */

void
PDM_part_to_part_ref_lnum2_get
(
 PDM_part_to_part_t *ptp,
 int               **n_ref_lnum2,
 int              ***ref_lnum2
);


/**
 *
 * \brief Get unreferenced Part2 elements
 *
 * \param [in]   ptp             Part-to-Part structure
 * \param [out]  n_unref_lnum2   Number of referenced Part2 elements
 * \param [out]  unref_lnum2     Unreferenced Part2 elements (1-based local ids)
 *
 */

void
PDM_part_to_part_unref_lnum2_get
(
 PDM_part_to_part_t  *ptp,
 int                **n_unref_lnum2,
 int               ***unref_lnum2
);


/**
 *
 * \brief Get Part2→Part1 mapping for referenced Part2 elements
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [out]  gnum1_come_from_idx Index for Part2→Part1 mapping (size = *n_part2*)
 * \param [out]  gnum1_come_from     Part2→Part1 mapping (global ids) <br>
 *                                   (for each part, size = \p gnum1_come_from_idx[*n_ref_lnum2*])
 *
 */

void
PDM_part_to_part_gnum1_come_from_get
(
 PDM_part_to_part_t *ptp,
 int              ***gnum1_come_from_idx,
 PDM_g_num_t      ***gnum1_come_from
);


/**
 *
 * \brief Initiate a non-blocking send from Part1 to Part2 (based on `MPI_issend`)
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Part1 data to send (order as in *part1_to_part2*)
 * \param [in]   tag                 Exchange tag
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_issend
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 const void        **part1_to_part2_data,
 const int           tag,
 int                *request
);


/**
 *
 * \brief Initiate a *raw*, non-blocking send from Part1 to Part2 (based on `MPI_issend`)
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Part1 data to send (order as in *part1_to_part2*)
 * \param [in]   tag                 Exchange tag
 * \param [out]  request             Request
 *
 */
void
PDM_part_to_part_issend_raw
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 const void         *raw_buffer,
 const int           tag,
 int                *request
);


/**
 *
 * \brief Finalize a non-blocking send from Part1 to Part2
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_issend_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
);


/**
 *
 * \brief Initiate a non-blocking send from Part2 to Part1 (based on `MPI_issend`)
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part2_to_part1_data Part2 data to send (order given by *gnum1_come_from* and *ref_lnum2* arrays)
 * \param [in]   tag                 Exchange tag
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_reverse_issend
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 const void        **part2_to_part1_data,
 const int           tag,
 int                *request
);


/**
 *
 * \brief Finalize a non-blocking send from Part2 to Part1
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_issend_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
);


/**
 *
 * \brief Test a non-blocking send from Part2 to Part1
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

int
PDM_part_to_part_reverse_issend_test
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Post (after test completion is OK) a non-blocking send from Part2 to Part1
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_issend_post
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Initiate a non-blocking recv (Part1→Part2) (based on `MPI_irecv`)
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part2_data          Part2 data (order given by *gnum1_come_from* and *ref_lnum2* arrays)
 * \param [in]   tag                 Exchange tag
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_irecv
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part2_data,
 int                 tag,
 int                *request
);

/**
 *
 * \brief Initiate a *raw*, non-blocking recv (Part1→Part2) (based on `MPI_irecv`)
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  raw_buffer    Buffer given by user
 * \param [in]  tag           Tag of the exchange
 * \param [out] request       Request
 *
 */

void
PDM_part_to_part_irecv_raw
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void               *raw_buffer,
 const int           tag,
 int                *request
);

/**
 *
 * \brief Finalize a non-blocking recv (Part1→Part2)
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_irecv_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Finalize a *raw*, non-blocking recv (Part1→Part2)
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */


void
PDM_part_to_part_irecv_wait_raw
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Test a non-blocking send (Part1→Part2)
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */


int
PDM_part_to_part_issend_test
(
 PDM_part_to_part_t *ptp,
 const int           request
);


/**
 *
 * \brief Test a non-blocking recv (Part1→Part2)
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

int
PDM_part_to_part_irecv_test
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Post (after test completion is OK) a non-blocking send from Part1 to Part2
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */
void
PDM_part_to_part_issend_post
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Post (after test completion is OK) a non-blocking recv (Part1→Part2)
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */
void
PDM_part_to_part_irecv_post
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Initiate a non-blocking recv (Part2→Part1) (based on `MPI_irecv`)
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_data          Part1 data to receive (ordered as in *part1_to_part2*)
 * \param [in]   tag                 Exchange tag
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_reverse_irecv
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part1_data,
 const int           tag,
 int                *request
);


/**
 *
 * \brief Finalize a non-blocking recv (Part2→Part1)
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_irecv_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Finalize a non-blocking recv (Part2→Part1)
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */

int
PDM_part_to_part_reverse_irecv_test
(
 PDM_part_to_part_t *ptp,
 const int           request
);

/**
 *
 * \brief Post (after test completion is OK) a non-blocking recv (Part2→Part1)
 *
 * \param [in]  ptp           Part-to-Part structure
 * \param [in]  request       Request
 *
 */
void
PDM_part_to_part_reverse_irecv_post
(
 PDM_part_to_part_t *ptp,
 const int           request
);


/**
 *
 * \brief Initiate a non-blocking exchange (Part1→Part2)
 *
 * \param [in]   ptp              Part-to-Part structure
 * \param [in]   k_comm           Kind of MPI communication
 * \param [in]   t_stride         Kind of stride
 * \param [in]   t_part1_data_def Kind of Part1 data definition
 * \param [in]   cst_stride       Constant stride
 * \param [in]   s_data           Data size
 * \param [in]   part1_stride     Stride of Part1 data (according to \p t_part1_data_def)
 * \param [in]   part1_data       Part1 data           (according to \p t_part1_data_def)
 * \param [out]  part2_stride     Stride of Part2 data (order given by *gnum1_come_from* and *ref_lnum2* arrays)
 * \param [out]  part2_data       Part2 data           (order given by *gnum1_come_from* and *ref_lnum2* arrays)
 * \param [out]  request          Request
 *
 */

void
PDM_part_to_part_iexch
(
 PDM_part_to_part_t                *ptp,
 const PDM_mpi_comm_kind_t          k_comm,
 const PDM_stride_t                 t_stride,
 const PDM_part_to_part_data_def_t  t_part1_data_def,
 const int                          cst_stride,
 const size_t                       s_data,
 const int                        **part1_stride,
 const void                       **part1_data,
 int                             ***part2_stride,
 void                            ***part2_data,
 int                               *request
);


/**
 *
 * \brief Finalize a non-blocking exchange (Part1→Part2)
 *
 * \param [in]  ptp      Part-to-Part structure
 * \param [in]  request  Request
 *
 */

void
PDM_part_to_part_iexch_wait
(
 PDM_part_to_part_t                *ptp,
 const int                          request
);


/**
 *
 * \brief Initiate a non-blocking exchange (Part2→Part1)
 *
 * \param [in]   ptp              Part-to-Part structure
 * \param [in]   k_comm           Kind of MPI communication
 * \param [in]   t_stride         Kind of stride
 * \param [in]   t_part2_data_def Kind of Part2 data definition
 * \param [in]   cst_stride       Constant stride
 * \param [in]   s_data           Data size
 * \param [in]   part2_stride     Stride of Part1 data (according to \p t_part2_data_def)
 * \param [in]   part2_data       Part1 data           (according to \p t_part2_data_def)
 * \param [out]  part1_stride     Stride of Part2 data (order given by *part1_to_part2*)
 * \param [out]  part1_data       Part2 data           (order given by *part1_to_part2*)
 * \param [out]  request          Request
 *
 */

void
PDM_part_to_part_reverse_iexch
(
 PDM_part_to_part_t                *ptp,
 const PDM_mpi_comm_kind_t          k_comm,
 const PDM_stride_t                 t_stride,
 const PDM_part_to_part_data_def_t  t_part2_data_def,
 const int                          cst_stride,
 const size_t                       s_data,
 const int                        **part2_stride,
 const void                       **part2_data,
 int                             ***part1_stride,
 void                            ***part1_data,
 int                               *request
);


/**
 *
 * \brief Finalize a non-blocking exchange (Part2→Part1)
 *
 * \param [in]  ptp      Part-to-Part structure
 * \param [in]  request  Request
 *
 */

void
PDM_part_to_part_reverse_iexch_wait
(
 PDM_part_to_part_t                *ptp,
 const int                          request
);


/**
 *
 * \brief Free a Part-to-Part structure
 *
 * \param [inout] ptp  Part-to-Part structure
 *
 * \return       NULL
 */

PDM_part_to_part_t *
PDM_part_to_part_free
(
 PDM_part_to_part_t *ptp
);


/**
 *
 * \brief Get number of partitions
 *
 * \param [in]  ptp       Pointer to \ref PDM_part_to_part_t object
 * \param [out] n_part1   Number of partitions in Part1
 * \param [out] n_part2   Number of partitions in Part2
 *
 */

void
PDM_part_to_part_n_part_get
(
 PDM_part_to_part_t *ptp,
 int                *n_part1,
 int                *n_part2
 );


/**
 *
 * \brief Get number of partitions and number of elements
 *
 * \param [in]  ptp       Pointer to \ref PDM_part_to_part_t object
 * \param [out] n_part1   Number of partitions in Part1
 * \param [out] n_part2   Number of partitions in Part2
 * \param [out] n_elt1    Number of Part1 elements (size = \p n_part1)
 * \param [out] n_elt2    Number of Part2 elements (size = \p n_part2)
 *
 */
void
PDM_part_to_part_n_part_and_n_elt_get
(
 PDM_part_to_part_t *ptp,
 int                *n_part1,
 int                *n_part2,
 int               **n_elt1,
 int               **n_elt2
 );


/**
 *
 * \brief Get referenced Part2 elements in current partition
 *
 * \param [in]   ptp           Part-to-Part structure
 * \param [in]   i_part        Id of current partition
 * \param [out]  n_ref_lnum2   Number of referenced Part2 elements
 * \param [out]  ref_lnum2     Referenced Part2 elements (zero-based local ids)
 *
 */

void
PDM_part_to_part_ref_lnum2_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                 *n_ref_lnum2,
       int                **ref_lnum2
);


/**
 *
 * \brief Get unreferenced Part2 elements in current partition
 *
 * \param [in]   ptp           Part-to-Part structure
 * \param [in]   i_part        Id of partition
 * \param [out]  n_unref_lnum2 Number of unreferenced Part2 elements
 * \param [out]  unref_lnum2   Unreferenced Part2 elements (zero-based local ids)
 *
 */

void
PDM_part_to_part_unref_lnum2_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                 *n_unref_lnum2,
       int                **unref_lnum2
);


/**
 *
 * \brief Get Part2→Part1 mapping for referenced Part2 elements in current partition
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   i_part              Id of current partition
 * \param [out]  gnum1_come_from_idx Index for Part2→Part1 mapping
 * \param [out]  gnum1_come_from     Part2→Part1 mapping
 *
 */

void
PDM_part_to_part_gnum1_come_from_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                **gnum1_come_from_idx,
       PDM_g_num_t        **gnum1_come_from
);

/**
 *
 * \brief Get selected numbers of part2 (only index)
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   i_part              Id of partition
 * \param [out]  n_elt1              Number of gnum1 element
 * \param [out]  part1_to_part2_idx  Index of data to send to gnum2 from gnum1
 *                                  (for each part size : \ref n_elt1+1)
 *
 */

void
PDM_part_to_part_part1_to_part2_idx_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                 *n_elt1,
       int                **part1_to_part2_idx
);

/**
 *
 * \brief Get Part1→Part2 mapping for current partition
 *
 * \param [in]   ptp                 Part-to-Part structure
 * \param [in]   i_part              Id of current partition
 * \param [out]  n_elt1              Number of Part1 elements
 * \param [out]  part1_to_part2_idx  Index for Part1→Part2 mapping (size : \p n_elt1 + 1)
 * \param [out]  part1_to_part2      Part1→Part2 mapping (global ids)
 *
 */

void
PDM_part_to_part_part1_to_part2_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                 *n_elt1,
       int                **part1_to_part2_idx,
       PDM_g_num_t        **part1_to_part2
);


/**
 *
 * \brief Get indirection from part1_to_part2 to send buffer (useful to setup buffer outside ptp) (??)
 *
 * \param [in]   ptp                       Part-to-Part structure
 * \param [out]  gnum1_to_send_buffer_idx  Index of data to send to gnum2 from gnum1 (for each part size : *n_elt1* + 1)
 * \param [out]  gnum1_to_send_buffer      For each gnum1 the position in send buffer
 *
 */
void
PDM_part_to_part_gnum1_to_send_buffer_get
(
 PDM_part_to_part_t    *ptp,
 int                 ***gnum1_to_send_buffer_idx,
 int                 ***gnum1_to_send_buffer
);

/**
 *
 * \brief Get indirection from ref_lnum2 to recv buffer  (useful to setup buffer outside ptp) (??)
 *
 * \param [in]   ptp                       Part-to-Part structure
 * \param [out]  recv_buffer_to_ref_lnum2  For each gnum2 the position in recv buffer (size = *gnum1_come_from_idx*[*n_ref_lnum2*])
 *
 */
void
PDM_part_to_part_recv_buffer_to_ref_lnum2_get
(
 PDM_part_to_part_t    *ptp,
 int                 ***recv_buffer_to_ref_lnum2
);


/**
 *
 * \brief Get buffer size and stride for send
 *
 * \param [in]   ptp                       Part-to-Part structure
 * \param [out]  default_n_send_buffer     Number of entities to send (size = *n_rank*)
 * \param [out]  default_i_send_buffer     Index (size = *n_rank* + 1)
 *
 */
void
PDM_part_to_part_default_send_buffer_get
(
 PDM_part_to_part_t    *ptp,
 int                  **default_n_send_buffer,
 int                  **default_i_send_buffer
);

/**
 *
 * \brief Get buffer size and stride for recv
 *
 * \param [in]   ptp                       Part-to-Part structure
 * \param [out]  default_n_recv_buffer     Number of entities to recv (size = *n_rank*)
 * \param [out]  default_i_recv_buffer     Index (size = *n_rank* + 1)
 *
 */
void
PDM_part_to_part_default_recv_buffer_get
(
 PDM_part_to_part_t    *ptp,
 int                  **default_n_recv_buffer,
 int                  **default_i_recv_buffer
);

/**
 *
 * \brief Get number of MPI ranks
 *
 * \param [in]   ptp          Part to part structure
 *
 * \return  Number of MPI ranks
 *
 */

int
PDM_part_to_part_n_ranks_get
(
 PDM_part_to_part_t    *ptp
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_tO_pART_H */
