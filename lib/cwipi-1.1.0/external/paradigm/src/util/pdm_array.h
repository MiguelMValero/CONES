/*
 * \file
 */

#ifndef __PDM_ARRAY_H__
#define __PDM_ARRAY_H__

/*----------------------------------------------------------------------------*/
#include "pdm.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*============================================================================
 * Interfaces des fonctions publiques
 *============================================================================*/

/* Utils functions creating new arrays */

// Create an array and fill it with 0
int* PDM_array_zeros_int(const int size);

// Create an array and fill it with given value
int*         PDM_array_const_int (const int size, const int         value);
PDM_g_num_t* PDM_array_const_gnum(const int size, const PDM_g_num_t value);
double*      PDM_array_const_double(const int size, const double value);

// Create an index array from a size array
int*         PDM_array_new_idx_from_sizes_int (const int *size_array, const int size);
PDM_g_num_t* PDM_array_new_idx_from_sizes_gnum(const int *size_array, const int size);
int*         PDM_array_new_idx_from_const_stride_int(const int stride, const int size);

/* Utils functions compararing arrays */

// Return 1 if the two arrays are equal, 0 otherwise
int PDM_array_are_equal_int(const int *array1, const int *array2, const int size);
int PDM_array_are_equal_gnum(const PDM_g_num_t *array1, const PDM_g_num_t *array2, const int size);

/* Utils functions modifying arrays*/

// Fill an array with the given value
void PDM_array_reset_int (int         *array, const int size, const int         value);
void PDM_array_reset_gnum(PDM_g_num_t *array, const int size, const PDM_g_num_t value);

// Compute an index array from a size array
void PDM_array_idx_from_sizes_int (const int *size_array, const int size, int         *idx_array);
void PDM_array_idx_from_sizes_gnum(const int *size_array, const int size, PDM_g_num_t *idx_array);

// Accumulate the values of an array
void PDM_array_accumulate_int (int         *array, const int size);
void PDM_array_accumulate_gnum(PDM_g_num_t *array, const int size);

/* Arrays algorithms */

// Count the number of occurences in a colored array
void PDM_array_count_per_col_int
(
 const int  n_col,
 const int  n_elem,
 const int *elem_col,
 int       *n_per_col
);
// Compute an index order to be used to read a colored array in increasing color order
void PDM_array_repart_per_col_int
(
 const int   n_col,
 const int   n_elem,
 const int *elem_col,
 int       *ordered_idx,
 int       *ordered
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_SEQ_H__ */
