/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_unique.h"
#include "pdm_array.h"
#include "pdm_quick_sort.h"
#include "pdm_radix_sort.h"
#include "pdm_binary_search.h"

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


/**
 *
 * \brief Unique
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
int
PDM_inplace_unique
(
 int a[],
 int l,
 int r
)
{
  int array_size = r - l + 1;
  if (array_size == 0) {
    return 0;
  }
  PDM_sort_int(&a[l], NULL, array_size);

  int new_size  = 1;
  int idx_write = l;
  PDM_g_num_t last_value = a[l];
  a[idx_write++] = last_value;
  for (int idx = l+1; idx <= r; idx++) {
    if(last_value != a[idx]){
      last_value = a[idx];
      a[idx_write++] = a[idx];
      new_size++;
    }
  }

  return new_size;
}

/**
 *
 * \brief Unique
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
int
PDM_inplace_unique_long
(
 PDM_g_num_t a[],
 int         order[],
 int l,
 int r
)
{
  // PDM_quick_sort_long(a, l, r); /* Less optimal than PDM_sort_long */
  int array_size = r - l + 1;
  if(array_size == 0) {
    return 0;
  }
  // printf("PDM_inplace_unique_long::array_size::%d\n", array_size);
  PDM_sort_long(&a[l], order, array_size);

  int new_size  = 1;
  int idx_write = l;
  PDM_g_num_t last_value = a[l];
  // if(order != NULL) {
  //   order[idx_write] = order[l];
  // }
  a[idx_write++] = last_value;
  for (int idx = l+1; idx <= r; idx++) {
    if(last_value != a[idx]){
      // if(order != NULL) {
      //   order[idx_write] = order[idx];
      // }
      last_value = a[idx];
      a[idx_write++] = a[idx];
      new_size++;
    }
  }

  return new_size;
}


/**
 *
 * \brief Same as unique but apply unique to order to know for each element the place in original array
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
int
PDM_inplace_unique_long_and_order
(
 PDM_g_num_t a[],
 int         order[],
 int l,
 int r
)
{
  // PDM_quick_sort_long(a, l, r); /* Less optimal than PDM_sort_long */
  int array_size = r - l + 1;
  if(array_size == 0) {
    return 0;
  }
  // printf("PDM_inplace_unique_long::array_size::%d\n", array_size);
  PDM_sort_long(&a[l], order, array_size);

  int new_size  = 1;
  int idx_write = l;
  PDM_g_num_t last_value = a[l];
  // if(order != NULL) {
  //   order[idx_write] = order[l];
  // }
  a[idx_write++] = last_value;
  for (int idx = l+1; idx <= r; idx++) {
    if(last_value != a[idx]){
      // if(order != NULL) {
      //   order[idx_write] = order[idx];
      // }
      last_value = a[idx];
      a    [idx_write] = a    [idx];
      order[idx_write] = order[idx];
      idx_write++;
      new_size++;
    }
  }

  return new_size;
}



/**
 *
 * \brief Unique
 *
 * \param [inout]   a             Array to sort
 * \param [inout]   unique_order  Unique index in old numbering
 * \param [in]      l             First element
 * \param [in]      r             Last  element
 *
 */
int
PDM_inplace_unique_long2
(
 PDM_g_num_t a[],
 int unique_order[],
 int l,
 int r
)
{
  int array_size = r - l + 1;
  if(array_size == 0) {
    return array_size;
  }
  // printf("PDM_inplace_unique_long::array_size::%d\n", array_size);
  int* order = (int *) malloc( (array_size) * sizeof(int));

  for(int i = 0; i < array_size; ++i){
    order[i] = i;
  }
  // PDM_radix_sort_long(&a[l], order, array_size);
  PDM_sort_long(&a[l], order, array_size);

  int new_size  = 1;
  int idx_write = l;
  PDM_g_num_t last_value = a[l];
  int idx_save = l;
  unique_order[order[0]] = idx_save;
  a[idx_write++] = last_value;
  for (int idx = l+1; idx <= r; idx++) {
    if(last_value != a[idx]){
      last_value = a[idx];
      // printf(" order[%d] = %d\n", idx-l, order[idx-l]);
      idx_save = idx_write;
      unique_order[order[idx-l]] = idx_save;
      a[idx_write++] = a[idx];
      new_size++;
    }
    unique_order[order[idx-l]] = idx_save;
  }

  free(order);

  return new_size;
}


int
PDM_unique_long_with_distrib
(
  PDM_MPI_Comm   comm,
  PDM_g_num_t   *dentity1_entity2_gnum,
  PDM_g_num_t   *distrib_entity2,
  int            array_size,
  int          **unique_order,
  PDM_g_num_t  **unique_dentity1_entity2_gnum
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  int *_unique_order = malloc( array_size * sizeof(int));
  int *order         = malloc( array_size * sizeof(int));
  int *rank_id_n     = PDM_array_zeros_int(n_rank+1);

  for(int i = 0; i < array_size; ++i) {
    int t_rank = PDM_binary_search_gap_long(PDM_ABS(dentity1_entity2_gnum[i]) - 1, distrib_entity2, n_rank + 1);
    rank_id_n[t_rank]++;
    _unique_order[i] = t_rank;
  }

  int *rank_id_idx = PDM_array_new_idx_from_sizes_int(rank_id_n, n_rank);
  int max_dblock = 0;
  for(int i = 0; i < n_rank; ++i) {
    max_dblock = PDM_MAX(max_dblock, distrib_entity2[i+1] - distrib_entity2[i]);
    rank_id_n[i] = 0;
  }

  // Swap by rank
  PDM_g_num_t *tmp_dentity1_entity2_gnum     = malloc( array_size * sizeof(PDM_g_num_t));
  PDM_g_num_t *_unique_dentity1_entity2_gnum = malloc( array_size * sizeof(PDM_g_num_t));
  for(int i = 0; i < array_size; ++i) {
    int t_rank    = _unique_order[i];
    int idx_write = rank_id_idx[t_rank]+rank_id_n[t_rank]++;
    tmp_dentity1_entity2_gnum[idx_write] = dentity1_entity2_gnum[i];
    order[idx_write] = i;
  }

  int* unique_tag = malloc(max_dblock * sizeof(int));

  int unique_size = 0;
  for(int i = 0; i < n_rank; ++i) {

    int dn_block =  distrib_entity2[i+1] - distrib_entity2[i];
    for(int j = 0; j < dn_block; ++j) {
      unique_tag[j] = -1;
    }

    // PDM_log_trace_array_long(&tmp_dentity_edge[ rank_id_idx[i]], rank_id_idx[i+1]-rank_id_idx[i], "tmp_dentity_edge :: ");
    for(int j = rank_id_idx[i]; j < rank_id_idx[i+1]; ++j) {

      int lnum = PDM_ABS(tmp_dentity1_entity2_gnum[j]) - distrib_entity2[i] - 1;
      // printf("lnum = %i \n", lnum);
      if(unique_tag[lnum] == -1) {
        unique_tag[lnum] = unique_size;
        _unique_dentity1_entity2_gnum[unique_size++] = tmp_dentity1_entity2_gnum[j];
      }
      _unique_order[order[j]] = unique_tag[lnum];
    }
  }

  _unique_dentity1_entity2_gnum = realloc(_unique_dentity1_entity2_gnum, unique_size * sizeof(PDM_g_num_t));

  free(rank_id_idx);
  free(rank_id_n);
  free(order);
  free(unique_tag);

  free(tmp_dentity1_entity2_gnum);

  *unique_order                 = _unique_order;
  *unique_dentity1_entity2_gnum = _unique_dentity1_entity2_gnum;

  return unique_size;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
