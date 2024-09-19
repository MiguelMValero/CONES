/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_order.h"

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
 * This function is part of Code_Saturne, a general-purpose CFD tool.
 *  Copyright (C) 1998-2014 EDF S.A.
 *
 * \brief Descend binary tree for the lexicographical ordering of a strided array
 *
 * \param [in]     number pointer to numbers of entities that should be ordered.
 * \param [in]            (if NULL, a default 1 to n numbering is considered)
 * \param [in]     stride stride of array (number of values to compare)
 * \param [in]     level  level of the binary tree to descend
 * \param [in]     nb_ent number of entities in the binary tree to descend
 * \param [in,out] order  ordering array
 */

static inline void
_order_lnum_descend_tree_s
(
const int    number[],
size_t       stride,
size_t       level,
const size_t nb_ent,
int          order[]
)
{
  size_t i_save, i1, i2, j, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      for (j = 0; j < stride; j++) {
        if (number[i1*stride + j] != number[i2*stride + j])
          break;
      }

      if (j < stride) {
        if (number[i1*stride + j] > number[i2*stride + j])
          lv_cur++;
      }

    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    for (j = 0; j < stride; j++) {
      if (number[i1*stride + j] != number[i2*stride + j])
        break;
    }

    if (j == stride) break;
    if (number[i1*stride + j] >= number[i2*stride + j]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = (int) i_save;
}



static inline void
_order_gnum_descend_tree_s
(
const PDM_g_num_t number[],
size_t            stride,
size_t            level,
const size_t      nb_ent,
int               order[]
)
{
  size_t i_save, i1, i2, j, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      for (j = 0; j < stride; j++) {
        if (number[i1*stride + j] != number[i2*stride + j])
          break;
      }

      if (j < stride) {
        if (number[i1*stride + j] > number[i2*stride + j])
          lv_cur++;
      }

    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    for (j = 0; j < stride; j++) {
      if (number[i1*stride + j] != number[i2*stride + j])
        break;
    }

    if (j == stride) break;
    if (number[i1*stride + j] >= number[i2*stride + j]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = (int) i_save;
}

inline
static
int
_lexicographic_compare_long
(
  const PDM_g_num_t *x,
  const PDM_g_num_t *y,
  const int          stride
)
{
  int res = x[0] == y[0];
  if(res == 1 && stride > 1) {
    return _lexicographic_compare_long(&x[1], &y[1], stride-1);
  }
  return x[0] < y[0];
}

inline
static
int
_lexicographic_equal_long
(
  const PDM_g_num_t *x,
  const PDM_g_num_t *y,
  const int          stride
)
{
  int res = x[0] == y[0];
  if(res == 1 && stride > 1) {
    return _lexicographic_equal_long(&x[1], &y[1], stride-1);
  }
  return x[0] == y[0];
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Order an array
 *
 * \param [in]      size_array       Number of elements
 * \param [in]      new_to_old_order   New order (size = \ref nElt
 * \param [in, out] Array           Array to renumber
 *
 */

void
PDM_order_array
(
const int     size_array,
const size_t  elt_size,
const int    *new_to_old_order,
void         *array
)
{
  unsigned char *old_array = (unsigned char *) malloc (size_array * elt_size);
  unsigned char *_array    = (unsigned char *) array;

  for (int i = 0; i < size_array; ++i) {
    for (int j = 0; j < (int) elt_size; ++j) {
      old_array[elt_size * i + j] = _array[elt_size * i + j];
    }
  }

  for (int i = 0; i < size_array; ++i) {
    for (int j = 0; j < (int) elt_size; ++j) {
      _array[elt_size * i + j] = old_array[elt_size * new_to_old_order[i] +j];
    }
  }

  free(old_array);
}

/**
 * This function is part of Code_Saturne, a general-purpose CFD tool.
 *  Copyright (C) 1998-2014 EDF S.A.
 *
 * \brief Order a strided array of global numbers lexicographically.
 *
 * \param [in]     number array of entity numbers (if NULL, a default 1 to n numbering is considered)
 * \param [in]     stride stride of array (number of values to compare)
 * \param [in,out] order  pre-allocated ordering table
 * \param [in]     nb_ent number of entities considered
 */

void
PDM_order_lnum_s
(
const int    number[],
size_t       stride,
int          order[],
const size_t nb_ent
)
{
  size_t i;
  int o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = (int) i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2) ;
  do {
    i--;
    _order_lnum_descend_tree_s(number, stride, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_lnum_descend_tree_s(number, stride, 0, i, order);
  }
}




/**
 * This function is part of Code_Saturne, a general-purpose CFD tool.
 *  Copyright (C) 1998-2014 EDF S.A.
 *
 * \brief Order a strided array of global numbers lexicographically.
 *
 * \param [in]     number array of entity numbers (if NULL, a default 1 to n numbering is considered)
 * \param [in]     stride stride of array (number of values to compare)
 * \param [in,out] order  pre-allocated ordering table
 * \param [in]     nb_ent number of entities considered
 */

void
PDM_order_gnum_s
(
const PDM_g_num_t number[],
size_t            stride,
int               order[],
const size_t      nb_ent
)
{
  size_t i;
  int o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = (int) i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2) ;
  do {
    i--;
    _order_gnum_descend_tree_s(number, stride, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_gnum_descend_tree_s(number, stride, 0, i, order);
  }
}



int
PDM_order_binary_search_long
(
 const PDM_g_num_t elt   [],
 const PDM_g_num_t array [],
 const size_t      stride,
 const size_t      nb_ent
)
{
  int left  = 0;
  int right = nb_ent - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if(_lexicographic_compare_long(elt, &array[stride*ind], stride)) {
      right = ind;
    } else {
      left = ind;
    }

    ind = (left + right) / 2;
  }

  if (_lexicographic_equal_long(elt, &array[stride*ind], stride)) {
    return ind;
  }

  if (_lexicographic_equal_long(elt, &array[stride*right], stride)) {
    return right;
  } else {
    return -1;
  }

}


int
PDM_order_inplace_unique_long
(
const int              n_entity,
const size_t           stride,
      PDM_g_num_t     *array,
      int             *order
)
{
  if(n_entity == 0) {
    return 0;
  }

  PDM_order_gnum_s(array, stride, order, n_entity);

  /* Apply sort */
  PDM_order_array (n_entity, stride * sizeof(PDM_g_num_t), order, array);

  PDM_g_num_t *last_value = malloc(stride * sizeof(PDM_g_num_t));

  int new_size  = 1;
  int idx_write = 1;
  for(int j = 0; j < (int) stride; ++j) {
    last_value[j] = array[j];
  }

  // PDM_log_trace_array_long(array, n_entity * stride, "array :: ");

  for(int i = 1; i < n_entity; i++){
    int is_same = _lexicographic_equal_long(last_value, &array[stride*i], stride);
    if(is_same == 0){ // N'est pas le meme
      for(int j = 0; j < (int) stride; ++j) {
        last_value[j] = array[stride*i+j];
        array[stride*idx_write+j] = last_value[j];
      }
      new_size++;
      idx_write++;
    }
  }

  free(last_value);

  return new_size;
}


int
PDM_order_inplace_unique_and_order_long
(
const int              n_entity,
const size_t           stride,
      PDM_g_num_t     *array,
      int             *order
)
{
  if(n_entity == 0) {
    return 0;
  }

  PDM_order_gnum_s(array, stride, order, n_entity);

  /* Apply sort */
  PDM_order_array (n_entity, stride * sizeof(PDM_g_num_t), order, array);

  PDM_g_num_t *last_value = malloc(stride * sizeof(PDM_g_num_t));

  int new_size  = 1;
  int idx_write = 1;
  for(int j = 0; j < (int) stride; ++j) {
    last_value[j] = array[j];
  }

  // PDM_log_trace_array_long(array, n_entity * stride, "array :: ");

  for(int i = 1; i < n_entity; i++){
    int is_same = _lexicographic_equal_long(last_value, &array[stride*i], stride);
    if(is_same == 0){ // N'est pas le meme
      for(int j = 0; j < (int) stride; ++j) {
        last_value[j] = array[stride*i+j];
        array[stride*idx_write+j] = last_value[j];
      }
      order[idx_write] = order[i];
      new_size++;
      idx_write++;
    }
  }

  free(last_value);

  return new_size;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */


