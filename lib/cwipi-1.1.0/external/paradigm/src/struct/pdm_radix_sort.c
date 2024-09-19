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
#include "pdm_array.h"
#include "pdm_radix_sort.h"

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
// const int _radix_shift = 7;
// const int _radix_base  = 128;
static const int _radix_shift = 10;
static const int _radix_base  = 1024;
// const int _radix_shift = 11;
// const int _radix_base  = 2048;
// const int _radix_shift = 6;
// const int _radix_base  = 64;
// const int _radix_shift = 8;
// const int _radix_base  = 256;
// const int base2_power[11] = {0, 2, 4, 8, 16, 32, 64, _radix_base, 256, 512, 1024};
// const int _radix_base = base2_power[radix_shift];

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static inline
PDM_g_num_t
_get_max_value_in_array_long
(
 PDM_g_num_t* v,
 int          size
)
{
  PDM_g_num_t max = v[0];
  for(int i = 1; i < size; ++i){
    if( v[i] > max){
      max = v[i];
    }
  }
  return max;
}


static inline
int*
_counting_sort
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int beg,
 int end,
 int place
)
{
  /* First step - Count */
  int n_buckets = _radix_base;
  int* count = PDM_array_zeros_int(n_buckets + 1);

  /* Count */
  for(int i = beg; i < end; ++i) {
    count[(v[i]/place)%_radix_base]++;
  }

  /* Rebuild array of displacement */
  for(int i = 1; i < n_buckets+1; ++i) {
    count[i] += count[i-1];
  }

  for(int i = end-1; i >= beg; i--){
    tmp[beg+count[(v[i]/place)%_radix_base]-1] = v[i];
    count[(v[i]/place)%_radix_base]--;
  }

  for(int i = beg; i < end; ++i) {
    v[i] = tmp[i];
  }

  return count;
}


// static inline
// int*
// _counting_sort_with_order
// (
//  PDM_g_num_t* v,
//  PDM_g_num_t* tmp,
//  int*         order,
//  int*         order_tmp,
//  int beg,
//  int end,
//  int place
// )
// {
//   /* First step - Count */
//   int n_buckets = _radix_base;
//   int* count = malloc( (n_buckets + 1) * sizeof(int));

//   /* Set to zero */
//   for(int i = 0; i < n_buckets+1; ++i){
//     count[i] = 0;
//   }

//   /* Count */
//   for(int i = beg; i < end; ++i) {
//     count[(v[i]/place)%_radix_base]++;
//   }

//   /* Rebuild array of displacement */
//   for(int i = 1; i < n_buckets+1; ++i) {
//     count[i] += count[i-1];
//   }

//   for(int i = end-1; i >= beg; i--){
//     tmp[beg+count[(v[i]/place)%_radix_base]-1] = v[i];
//     order_tmp[beg+count[(v[i]/place)%_radix_base]-1] = order[i];
//     count[(v[i]/place)%_radix_base]--;
//   }

//   for(int i = beg; i < end; ++i) {
//     v[i] = tmp[i];
//     order[i] = order_tmp[i];
//   }

//   return count;
// }


static inline
int*
_counting_sort_with_order
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int*         order,
 int*         order_tmp,
 int beg,
 int end,
 int place,
 int place_power
)
{
  PDM_UNUSED(place);
  /* First step - Count */
  int n_buckets = _radix_base;
  int* count = PDM_array_zeros_int(n_buckets + 1);

  /* Count */
  for(int i = beg; i < end; ++i) {
    count[(v[i] >> place_power) % _radix_base]++;
  }

  /* Rebuild array of displacement */
  for(int i = 1; i < n_buckets+1; ++i) {
    count[i] += count[i-1];
  }

  for(int i = end-1; i >= beg; i--){
    int idx = (v[i] >> place_power) % _radix_base;
    tmp[beg+count[idx]-1] = v[i];
    order_tmp[beg+count[idx]-1] = order[i];
    count[idx]--;
  }

  for(int i = beg; i < end; ++i) {
    v[i] = tmp[i];
    order[i] = order_tmp[i];
  }

  return count;
}

static inline
void
_counting_sort_with_order_recurse
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int*         order,
 int*         order_tmp,
 int          beg,
 int          end,
 int          place,
 int          place_power,
 int*         count_curr,
 int*         count_next,
 int          iteration
)
{
  PDM_UNUSED(place);
  /* Set to zero */
  if(iteration == 0){
    for(int i = 0; i < _radix_base+1; ++i){
      count_curr[i] = 0;
      count_next[i] = 0;
    }

    /* Count */
    for(int i = beg; i < end; ++i) {
      count_curr[(v[i] >> place_power) % _radix_base]++;
    }

  } else {
    PDM_array_reset_int(count_next, _radix_base+1, 0);
  }

  /* Rebuild array of displacement */
  for(int i = 1; i < _radix_base+1; ++i) {
    count_curr[i] += count_curr[i-1];
  }

  // printf("count_curr[%d/%d]:: ", iteration, place_power);
  // for(int i = 1; i < _radix_base+1; ++i) {
  //   printf("%d ", count_curr[i]);
  // }
  // printf("\n");
  int next_shift = place_power+_radix_shift;

  for(int i = end-1; i >= beg; i--){
    int idx   = (v[i] >> place_power) % _radix_base;
    int idx_t = beg+count_curr[idx]-1;

    tmp[idx_t]       = v[i];
    order_tmp[idx_t] = order[i];

    count_curr[idx]--;

    /* Prepare for next step */
    count_next[(v[i] >> next_shift ) % _radix_base]++;
  }

  for(int i = beg; i < end; ++i) {
    v[i]     = tmp[i];
    order[i] = order_tmp[i];
  }
}



static inline
void
_cc_radix_sort
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int beg,
 int end,
 int place
)
{
  int place_init = place*_radix_base;
  place = 1;
  while(place < place_init){
    int* range = _counting_sort(v, tmp, beg, end, place);
    free(range);
    place *= _radix_base;
  }

}

// static inline
// void
// _cc_radix_sort_with_order
// (
//  PDM_g_num_t* v,
//  PDM_g_num_t* tmp,
//  int*         order,
//  int*         order_tmp,
//  int beg,
//  int end,
//  int place,
//  int place_power
// )
// {
//   int place_init = place*_radix_base;
//   place = 1;
//   place_power = 0;
//   printf("_cc_radix_sort_with_order::%d --> %d :: %d \n", end-beg, place_init, place_power);
//   while(place < place_init){
//     int* range = _counting_sort_with_order(v, tmp, order, order_tmp, beg, end, place, place_power);
//     free(range);
//     place *= _radix_base;
//     place_power += _radix_shift;
//   }
//   printf("_cc_radix_sort_with_order end ::%d --> %d :: %d \n", end-beg, place_init, place_power);

// }

static inline
void
_cc_radix_sort_with_order
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int*         order,
 int*         order_tmp,
 int          beg,
 int          end,
 int          place,
 int          place_power
)
{
  int place_init = place*_radix_base;
  place = 1;
  place_power = 0;

  int iteration = 0;

  int* count_curr = malloc( (_radix_base + 1) * sizeof(int));
  int* count_next = malloc( (_radix_base + 1) * sizeof(int));

  int* count_tmp;
  int* _count_curr = count_curr;
  int* _count_next = count_next;

  // printf("_cc_radix_sort_with_order::%d --> %d :: %d \n", end-beg, place_init, place_power);
  while(place < place_init){

    // printf("_cc_radix_sort_with_order in while :: %p | %p \n", _count_curr, _count_next);

    // printf("_cc_radix_sort_with_order::%d --> %d :: %d \n", end-beg, place_init, place_power);

    _counting_sort_with_order_recurse(v, tmp, order, order_tmp, beg, end,
                                      place, place_power,
                                      _count_curr, _count_next, iteration);
    place       *= _radix_base;
    place_power += _radix_shift;
    iteration   += 1;

    /* Swap ptr */
    count_tmp   = _count_curr;
    _count_curr = _count_next;
    _count_next = count_tmp;

  }
  // printf("_cc_radix_sort_with_order end ::%d --> %d :: %d \n", end-beg, place_init, place_power);

  free(count_curr);
  free(count_next);

}


static inline
void
_std_radix_sort
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int beg,
 int end,
 int place
)
{
  // int cache_size = 20000;
  // int cache_size = -1;//(int) pow(_radix_base, 2);
  int cache_size = -1;//(int) pow(_radix_base, 2);

  if(beg == end){
    return;
  }

  int* range = _counting_sort(v, tmp, beg, end, place);

  place /= _radix_base;

  if(place == 0) {
    free(range);
    return;
  }

  /* Il faut aglomerer les ranges succesives pour gagner si elle sont trop petite */
  for(int i = 0; i < _radix_base; ++i){
    // fmt::print("Manage {0} = {1} / {2}\n", i, range[i], range[i+1]);
    if( (range[i+1] - range[i]) == 0){
    } else if( (range[i+1] - range[i]) > cache_size){
      // fmt::print("Case 1 :: {0} - {1} \n ", range[i], range[i+1]);
      _std_radix_sort(v, tmp, beg+range[i], beg+range[i+1] , place);
    } else if( (range[i+1] - range[i]) <= cache_size) {
      // fmt::print("Case 3 :: {0} - {1} \n ", range[i], range[i+1]);
      _cc_radix_sort(v, tmp, beg+range[i], beg+range[i+1], place);
    }
  }

  free(range);

}


static inline
void
_std_radix_sort_with_order
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int*         order,
 int*         order_tmp,
 int          beg,
 int          end,
 int          place,
 int          place_power
)
{
  // int cache_size = 900000;
  // int cache_size = (int) pow(_radix_base, 10);
  // int cache_size = (int) pow(_radix_base, 3);
  // int cache_size = (int) pow(256, 2);
  // int cache_size = -1;//(int) pow(_radix_base, 2);
  // int cache_size = (int) pow(256, 3);
  int cache_size = -1; //(int) pow(64, 2);

  if(beg == end){
    return;
  }

  // printf("_std_radix_sort_with_order::%d --> %d | place_power::%d\n", end-beg, place, place_power);
  int* range = _counting_sort_with_order(v, tmp, order, order_tmp, beg, end, place, place_power);

  place /= _radix_base;
  place_power -= _radix_shift;
  if(place_power < 0){
    place_power = 0;
  }

  if(place == 0) {
    free(range);
    return;
  }

  /* Il faut aglomerer les ranges succesives pour gagner si elle sont trop petite */
  for(int i = 0; i < _radix_base; ++i){
    if( (range[i+1] - range[i]) == 0){
    } else if( (range[i+1] - range[i]) > cache_size){
      _std_radix_sort_with_order(v, tmp, order, order_tmp, beg+range[i], beg+range[i+1] , place, place_power);
    } else if( (range[i+1] - range[i]) <= cache_size) {
      _cc_radix_sort_with_order(v, tmp, order, order_tmp, beg+range[i], beg+range[i+1], place, place_power);
    }
  }

  free(range);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout] array        Array to sort
 * \param [inout] order        new indice to old indice (or NULL)
 * \param [in]    lArray       Array length
 *
 */
void
PDM_radix_sort_long
(
 PDM_g_num_t *array,
 int         *order,
 int          lArray
)
{
  if(lArray == 0) {
    return;
  }
  PDM_g_num_t max = _get_max_value_in_array_long(array, lArray);

  // int n_step = -1;
  int n_step = 0;
  PDM_g_num_t lar = max;
  while(lar > 0){
    n_step++;
    lar /= _radix_base;
  }

  // Il faut trouver le max dans la base pour reprensenter le min et le max
  //  --> si beaucoup d'écart on tente la moyenne
  PDM_g_num_t* tmp = (PDM_g_num_t *) malloc( (lArray+1) * sizeof(PDM_g_num_t));

  int place = (int) pow(_radix_base, (n_step-1));
  int place_power = _radix_shift*(n_step-1);

  // printf("PDM_radix_sort_long::place ::%d \n", place);
  // printf("PDM_radix_sort_long::max   ::"PDM_FMT_G_NUM" \n", max);
  // printf("PDM_radix_sort_long::n_step::%d \n", n_step);
  // abort();

  if(order == NULL){
    _std_radix_sort(array, tmp, 0, lArray, place);
  } else {
    int* order_tmp = (int *) malloc( (lArray+1) * sizeof(int));
    _std_radix_sort_with_order(array, tmp, order, order_tmp, 0, lArray, place, place_power);
    // _std_radix_sort(array, tmp, 0, lArray, place);
    free(order_tmp);
  }


  free(tmp);

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
