/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#include "pdm_distrib.h"
#include "pdm_binary_search.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 * \brief   Evaluate a distribution array.
 *
 * \param [in]  n_ranges     Number of ranges in the distribution
 * \param [in]  distribution Number of elements associated to each range of the distribution
 * \param [in]  optim        Optimal count in each range
 *
 * \return  a fit associated to the distribution. If fit = 0, distribution is perfect.
 *
 */
static double
_evaluate_distribution(int          n_ranges,
                       double      *distribution,
                       double       optim)
{
  int  i;
  double  d_low = 0, d_up = 0, fit = 0;

  /*
     d_low is the max gap between the distribution count and the optimum when
     distribution is lower than optimum.
     d_up is the max gap between the distribution count and the optimum when
     distribution is greater than optimum.
  */

  for (i = 0; i < n_ranges; i++) {

    if (distribution[i] > optim)
      d_up = PDM_MAX(d_up, distribution[i] - optim);
    else
      d_low = PDM_MAX(d_low, optim - distribution[i]);

  }

  fit = (d_up + d_low) / optim;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0)
    PDM_printf( "<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n",
               optim, fit);
#endif

  return  fit;
}


/**
 * \brief Define a global distribution associated to a sampling array i.e. count
 * the number of elements in each range.
 *
 *   \param [in]    dim           2D or 3D
 *   \param [in]    n_ranks       number of ranks (= number of ranges)
 *   \param [in]    gsum_weight   global sum of all weightings
 *   \param [in]    n_codes       local number of Hilbert codes
 *   \param [in]    hilbert_codes local list of Hilbert codes to distribute
 *   \param [in]    weight        weighting related to each code
 *   \param [in]    order         ordering array
 *   \param [in]    sampling      sampling array
 *   \param [inout] c_freq        pointer to the cumulative frequency array
 *   \param [inout] g_distrib     pointer to a distribution array
 *   \param [in]    comm          mpi communicator
 */

static void
_define_rank_distrib(const int             sampling_factor,
                     const int             n_part,
                     const int            *n_elt,
                     const PDM_g_num_t   **gnum_elt,
                     const double        **weight,
                     const int             n_ranks,
                     const double          gsum_weight,
                     const PDM_g_num_t     sampling[],
                           double          cfreq[],
                           double          g_distrib[],
                           PDM_MPI_Comm    comm)
{
  int  id, rank_id;

  const int  n_samples       = sampling_factor * n_ranks;

  /* Initialization */
  double   *l_distrib = (double  *) malloc (n_samples * sizeof(double));

  for (id = 0; id < n_samples; id++) {
    l_distrib[id] = 0;
    g_distrib[id] = 0;
  }

  if(weight != NULL) {
    for (int i = 0; i < n_part; i++) {
      for (int j = 0; j < n_elt[i]; j++) {
        PDM_g_num_t _gnum_elt = PDM_ABS(gnum_elt[i][j]) - 1;
        int i_sample = PDM_binary_search_gap_long(_gnum_elt, sampling, n_samples + 1);
        l_distrib[i_sample] += weight[i][j];
      }
    }
  } else {
    for (int i = 0; i < n_part; i++) {
      for (int j = 0; j < n_elt[i]; j++) {
        PDM_g_num_t _gnum_elt = PDM_ABS(gnum_elt[i][j]) - 1;
        int i_sample = PDM_binary_search_gap_long(_gnum_elt, sampling, n_samples + 1);
        l_distrib[i_sample] += 1.;
      }
    }
  }

  /* Define the global distribution */
  PDM_MPI_Allreduce(l_distrib, g_distrib, n_samples,
                    PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  free(l_distrib);

  /* Define the cumulative frequency related to g_distribution */
  cfreq[0] = 0.;
  for (id = 0; id < n_samples; id++) {
    double _g_distrib  = (double)g_distrib[id];
    double _gsum_weight = (double)gsum_weight;
    cfreq[id+1] = cfreq[id] + _g_distrib/_gsum_weight;
  }
  cfreq[n_samples] = 1.0;

#if 0 && defined(DEBUG) && !defined(DEBUG) /* For debugging purpose only */

  if (cs_glob_rank_id <= 0) {

    FILE  *dbg_file = NULL;
    int  len;
    static int  loop_id1 = 0;

    len = strlen("DistribOutput_l.dat")+1+2;
    char  *rfilename = (char *) malloc (len * sizeof(char));
    sprintf(rfilename, "DistribOutput_l%02d.dat", loop_id1);

    loop_id1++;

    dbg_file = fopen(rfilename, "w");

    fprintf(dbg_file,
            "# Sample_id  |  OptCfreq  |  Cfreq  |  Sampling  |"
            "Global Distrib\n");
    for (i = 0; i < n_samples; i++)
      fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
              i, (double)i/(double)n_samples, cfreq[i],
              (double)(sampling[i]), distrib[i]);
    fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
            i, 1.0, 1.0, 1.0, 0);

    fclose(dbg_file);
    free(rfilename);

  }

#endif /* debugging output */

  /* Convert global distribution from n_samples to n_ranks */

  for (rank_id = 0; rank_id < n_ranks; rank_id++) {

    double   sum = 0.;
    int   shift = rank_id * sampling_factor;

    for (id = 0; id < sampling_factor; id++) {
      sum += g_distrib[shift + id];
    }
    g_distrib[rank_id] = sum;

  } /* End of loop on ranks */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Sanity check in debug */
  {
    PDM_g_num_t   sum = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      sum += g_distrib[rank_id];

    if (sum != gsum_weight)
      PDM_error(__FILE__, __LINE__, 0,
                "Error while computing global distribution.\n"
                "sum = %u and gsum_weight = %u\n",
                sum, gsum_weight);
    exit(1);
  }
#endif /* sanity check */

}


/**
 * \brief Update a distribution associated to sampling to assume a well-balanced
 * distribution of the leaves of the tree.
 *
 *   \param [in]    dim      1D, 2D or 3D
 *   \param [in]    n_ranks  number of ranks (= number of ranges)
 *   \param [inout] c_freq   cumulative frequency array
 *   \param [inout] sampling pointer to pointer to a sampling array
 *   \param [in]    comm     mpi communicator
 */

static void
_update_sampling(int            sampling_factor,
                 int            n_ranks,
                 double         c_freq[],
                 PDM_g_num_t  *sampling[])
{
  int  i, j, next_id;
  double  target_freq, f_high, f_low, delta;
  PDM_g_num_t  s_low, s_high;

  PDM_g_num_t  *new_sampling = NULL, *_sampling = *sampling;

  const int  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute new_sampling */

  new_sampling = ( PDM_g_num_t  *) malloc (sizeof(PDM_g_num_t) * (n_samples + 1));

  new_sampling[0] = _sampling[0];

  next_id = 1;

  for (i = 0; i < n_samples; i++) {

    target_freq = (i+1)*unit;

    /* Find the next id such as c_freq[next_id] >= target_freq */

    for (j = next_id; j < n_samples + 1; j++) {
      if (c_freq[j] >= target_freq) {
        next_id = j;
        break;
      }
    }

    /* Find new s such as new_s is equal to target_freq by
       a linear interpolation */

    f_low = c_freq[next_id-1];
    f_high = c_freq[next_id];

    s_low = _sampling[next_id-1];
    s_high = _sampling[next_id];

    if (f_high - f_low > 0) {
      delta = (target_freq - f_low) * (s_high - s_low) / (f_high - f_low);
      new_sampling[i+1] = (PDM_g_num_t) (s_low + delta);
    }
    else /* f_high = f_low */
      new_sampling[i+1] = (PDM_g_num_t) (s_low + 0.5 * (s_low + s_high));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    PDM_printf( " <_update_distrib> (rank: %d) delta: %g, target: %g,"
               " next_id: %d, f_low: %g, f_high: %g, s_low: %g, s_high: %g\n"
               "\t => new_sampling: %g\n",
               cs_glob_rank_id, delta, target_freq, next_id,
               f_low, f_high, s_low, s_high, new_sampling[i+1]);
#endif

  } /* End of loop on samples */

  new_sampling[n_samples] = _sampling[n_samples];

  free(_sampling);

  /* Return pointers */

  *sampling = new_sampling;
}


/*=============================================================================
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
 const int           dnelt,
       PDM_g_num_t  *elt_distrib,
       int           offset,
 const PDM_MPI_Comm  comm
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Compute distribution for element */

  // PDM_g_num_t* elt_distrib = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t  _dnelt      = (PDM_g_num_t) dnelt;

  PDM_MPI_Allgather((void *) &_dnelt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&elt_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  elt_distrib[0] = 1+offset;
  for (int i = 1; i < n_rank+1; i++) {
    elt_distrib[i] +=  elt_distrib[i-1];
  }

  /* Verbose */
  if (1 == 0) {
    PDM_printf("elt_distrib : "PDM_FMT_G_NUM,  elt_distrib[0]);
    for (int i = 1; i < n_rank+1; i++) {
      PDM_printf(" "PDM_FMT_G_NUM, elt_distrib[i]);
    }
    PDM_printf("\n");
  }
}

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
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t* dentity_proc = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));

  /*
   * Exchange
   */

  PDM_g_num_t _dn_entity = (PDM_g_num_t) dn_entity;
  PDM_MPI_Allgather((void *) &_dn_entity,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&dentity_proc[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  dentity_proc[0] = 0;
  for (int i = 1; i < n_rank+1; i++) {
    dentity_proc[i] = dentity_proc[i] + dentity_proc[i-1];
  }

  return dentity_proc;
}


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
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t step      = n_g_entity/n_rank;
  PDM_g_num_t remainder = n_g_entity%n_rank;

  int dn_elmt = step;
  if (i_rank < remainder) { /* Distribute the remainder */
    dn_elmt += 1;
  }

  return dn_elmt;
}

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
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t* dentity_proc = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));

  PDM_g_num_t step      = n_g_entity/n_rank;
  PDM_g_num_t remainder = n_g_entity%n_rank;

  dentity_proc[0] = 0;
  for (int i = 1; i < n_rank + 1; i++) {
    dentity_proc[i] = step;
    const int i1 = i - 1;
    if (i1 < remainder) { /* Distribute the remainder */
      dentity_proc[i]  += 1;
    }
  }

  for (int i = 1; i < n_rank + 1; i++) {
    dentity_proc[i] += dentity_proc[i-1];
  }

  return dentity_proc;
}


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
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * Compute the max
   */
  PDM_g_num_t _id_max = 0;
  PDM_g_num_t n_g_entity = 0;

  for(int i_part = 0; i_part < n_part; ++i_part) {
    for(int i_elmt = 0; i_elmt < n_elmts[i_part]; ++i_elmt) {
      _id_max = PDM_MAX (_id_max, ln_to_gn[i_part][i_elmt]);
    }
  }

  // double t1 = PDM_MPI_Wtime();
  PDM_MPI_Allreduce (&_id_max,
                     &n_g_entity,
                     1,
                     PDM__PDM_MPI_G_NUM,
                     PDM_MPI_MAX,
                     comm);
  // double t2 = PDM_MPI_Wtime();
  // double dt = t2-t1;
  // printf("[%i] dt = %12.5e \n ", i_rank, dt);

  return PDM_compute_uniform_entity_distribution(comm, n_g_entity);
}



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

)
{

  PDM_g_num_t _id_max     = 0;
  PDM_g_num_t _id_max_max = 0;

  for (int i = 0; i < n_part; i++) {
    for (int j = 0; j < n_elmts[i]; j++) {
      PDM_g_num_t gnum = PDM_ABS(ln_to_gn[i][j]);
      _id_max = PDM_MAX (_id_max, gnum);
    }
  }

  PDM_MPI_Allreduce (&_id_max,
                     &_id_max_max,
                     1,
                     PDM__PDM_MPI_G_NUM,
                     PDM_MPI_MAX,
                     comm);


  const int n_samples = sampling_factor * n_active_ranks;

  PDM_g_num_t *sampling = malloc(sizeof(PDM_g_num_t) * (n_samples + 1));

  double  lsum_weight = 0.;
  if(weight != NULL) {
    for (int i = 0; i < n_part; i++) {
      for (int j = 0; j < n_elmts[i]; j++) {
        lsum_weight += weight[i][j];
      }
    }
  } else {
    for (int i = 0; i < n_part; i++) {
      lsum_weight += n_elmts[i];
    }
  }

  double  gsum_weight = 0.;
  PDM_MPI_Allreduce(&lsum_weight, &gsum_weight, 1,
                    PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  double optim = gsum_weight / n_active_ranks;

  /* Define a naive sampling (uniform distribution) */
  PDM_g_num_t _n_sample_data = _id_max_max / n_samples;
  PDM_g_num_t _samplerest    = _id_max_max % n_samples;

  sampling[0] = 0;
  int k = 0;
  for (int i = 0; i < n_samples; i++) {
    sampling[i+1] = sampling[i];
    sampling[i+1] += _n_sample_data;
    if (k < _samplerest) {
      sampling[i+1] += 1;
      k += 1;
    }
  }

  /* Define the distribution associated to the current sampling array */

  double *distrib = (double *) malloc (sizeof(double) * n_samples);
  double  *cfreq = (double *) malloc (sizeof(double) * (n_samples + 1));

  _define_rank_distrib(sampling_factor,
                       n_part,
                       n_elmts,
                       ln_to_gn,
                       weight,
                       n_active_ranks,
                       gsum_weight,
                       sampling,
                       cfreq,
                       distrib,
                       comm);

  /* Initialize best choice */

  double fit = _evaluate_distribution(n_active_ranks, distrib, optim);
  double best_fit = fit;

  PDM_g_num_t  *best_sampling = (PDM_g_num_t  *) malloc (sizeof(PDM_g_num_t) * (n_samples + 1));

  for (int i = 0; i < (n_samples + 1); i++) {
    best_sampling[i] = sampling[i];
  }

  /* Loop to get a better sampling array */
  for (int n_iters = 0;
       (   n_iters < n_iter_max
           && fit > tolerance);
       n_iters++)  {

    _update_sampling(sampling_factor, n_active_ranks, cfreq, &sampling);

    /* Compute the new distribution associated to the new sampling */
    _define_rank_distrib(sampling_factor,
                         n_part,
                         n_elmts,
                         ln_to_gn,
                         weight,
                         n_active_ranks,
                         gsum_weight,
                         sampling,
                         cfreq,
                         distrib,
                         comm);


    fit = _evaluate_distribution(n_active_ranks, distrib, optim);

    /* Save the best sampling array and its fit */

    if (fit < best_fit) {

      best_fit = fit;
      for (int i = 0; i < (n_samples + 1); i++) {
        best_sampling[i] = sampling[i];
      }
    }
  } /* End of while */

  free (distrib);
  free (cfreq);
  free (sampling);

  sampling = best_sampling;

  PDM_g_num_t *_rank_index = malloc (sizeof(PDM_g_num_t) * (n_active_ranks + 1));
  for (int i = 0; i < n_active_ranks + 1; i++) {
    int id = i * sampling_factor;
    _rank_index[i] = sampling[id];
  }
  free(sampling);

  *rank_index = _rank_index;
}



#ifdef __cplusplus
}
#endif /* __cplusplus */
