/*
 * \file
 */

#ifndef PDM_DISTRIB_H
#define PDM_DISTRIB_H

/*----------------------------------------------------------------------------
 *  Local headers
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


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Compute distribution from dn_elt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     offset               Can be -1 or 0 to shift the first elements of distribution (-1 is the standard and the the distrib begin at 0 )
 * \param [in]     comm                 MPI Communicator
 */
void
PDM_distrib_compute
(
 const int           dn_elt,
       PDM_g_num_t  *elt_distrib,
       int           offset,
 const PDM_MPI_Comm  comm
);

/**
 * \brief Compute distribution from dNelmt
 *
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 * \return elt_distrib, Distribution of elements on processes (size = n_rank+1)
 */
PDM_g_num_t*
PDM_compute_entity_distribution
(
 const PDM_MPI_Comm     comm,
 const int              dn_entity
);


/**
 * \brief Compute an uniform size for all rank with step and reminder from the total number of global entity
 *        (All ranks should have the same n_g_entity)
 *
 * \param [in]     comm        MPI Communicator
 * \param [in]     n_g_entity  Global number of entity
 * \return dn_elmt, Number of element on current process
 */
int
PDM_compute_uniform_dn_entity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t      n_g_entity
);


/**
 * \brief Compute an uniform distribution array for all rank with step and reminder from the total number of global entity
 *        (All ranks should have the same n_g_entity)
 *
 * \param [in]     comm        MPI Communicator
 * \param [in]     n_g_entity  Global number of entity
 * \return elt_distrib, Distribution of elements on processes (size = n_rank+1)
 */
PDM_g_num_t*
PDM_compute_uniform_entity_distribution
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t      n_g_entity
);


/**
 * \brief Compute an uniform distribution array for all rank with step and reminder from the total number of global entity
 *        (All ranks should have the same n_g_entity). This function automaticly compute the total number of entity and setp a uniform distribution
 *
 * \param [in]     comm        MPI Communicator
 * \param [in]     n_part      Number of partition in current process
 * \param [in]     n_elmts     Number of elements for each partition
 * \param [in]     ln_to_gn    Local to global numbering for each partition (size = n_part, and each component have size pn_elmt[i_part])
 * \return elt_distrib, Distribution of elements on processes (size = n_rank+1)
 */
PDM_g_num_t*
PDM_compute_uniform_entity_distribution_from_partition
(
 const PDM_MPI_Comm     comm,
 const int              n_part,
 const int             *n_elmts,
 const PDM_g_num_t    **ln_to_gn
);


/**
 * \brief Compute an equilibrate distribution array for all rank.
 *        Algorithm can take weight to equilibrate among all process.
 * \param [in]     sampling_factor   Size of the sampling of distribution. Typical value are in range [1, 4].
 * \param [in]     n_active_ranks    Number of ranks actives to computes samplings
 * \param [in]     n_part            Number of partition in current process
 * \param [in]     n_elmts           Number of elements for each partition
 * \param [in]     ln_to_gn          Local to global numbering for each partition (size = n_part, and each component have size pn_elmt[i_part])
 * \param [in]     weight            Weight associte to each elements
 * \param [in]     n_iter_max        Maximum iteration of refinement
 * \param [in]     tolerance         Tolerance for load imbalance
 * \param [in]     comm              MPI Communicator
 * \param [in]     n_g_entity        Global number of entity
 * \param [out]    rank_index        Distributation among n_active_ranks (size = n_active_ranks)
 */
void
PDM_distrib_weight
(
  const int            sampling_factor,
  const int            n_active_ranks,
  const int            n_part,
  const int           *n_elmts,
  const PDM_g_num_t  **ln_to_gn,
  const double       **weight,
  const int            n_iter_max,
  const double         tolerance,
  const PDM_MPI_Comm   comm,
        PDM_g_num_t  **rank_index
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif // PDM_DISTRIB_H
