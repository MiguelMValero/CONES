/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_quick_sort.h"
#include "pdm_sort.h"
#include "pdm_compare_operator.h"

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
 * \brief Equal operator for connectivities
 *
 */
int
PDM_operator_equal_string
(
const void* a,
const void* b,
      void* ctxt)
{
  int i = *(const int *) a;
  int j = *(const int *) b;

  PDM_user_defined_sort* us = (PDM_user_defined_sort*) ctxt;

  int ni = us->idx[i+1] - us->idx[i];
  int nj = us->idx[j+1] - us->idx[j];

  // log_trace("PDM_operator_equal_string:: %d %d - %d %d - %d %d \n", i, j, ni, nj, us->idx[i], us->idx[j]);

  char* arr_i = (char *) &us->arr[us->idx[i]*sizeof(char)];
  char* arr_j = (char *) &us->arr[us->idx[j]*sizeof(char)];

  if(ni != nj){
    return 0;
  } else {
    int i_comp = strncmp(arr_i, arr_j, ni);
    if(i_comp == 0) {
      return 1;
    } else {
      return 0;
    }
  }
}

/**
 *
 * \brief Equal operator for connectivities
 *
 */
int
PDM_operator_equal_connectivity
(
const void* a,
const void* b,
      void* ctxt)
{
  int i = *(const int *) a;
  int j = *(const int *) b;


  PDM_user_defined_sort* us = (PDM_user_defined_sort*) ctxt;
  // printf("PDM_operator_equal_connectivity:: %d %d - %d %d - %d %d \n", i, j, ni, nj, us->idx[i], us->idx[j]);

  if( us->key[i] == us->key[j]){
    int ni = us->idx[i+1] - us->idx[i];
    int nj = us->idx[j+1] - us->idx[j];

    if(ni == nj){

      int* arr_i = (int*) &us->arr[us->idx[i]*sizeof(int)];
      int* arr_j = (int*) &us->arr[us->idx[j]*sizeof(int)];

      /* Dans notre cas on veut sort les entiers avant de les comparers */
      int* sort_arr_i = (int*) malloc( ni * sizeof(int));
      int* sort_arr_j = (int*) malloc( ni * sizeof(int));

      for(int k = 0; k < ni; ++k){
        sort_arr_i[k] = arr_i[k];
        sort_arr_j[k] = arr_j[k];
      }
      PDM_quick_sort_int(sort_arr_i, 0, ni-1);
      PDM_quick_sort_int(sort_arr_j, 0, ni-1);

      for(int k = 0; k < ni; ++k){
        // printf(" \t sort_arr_i[%d] = %d | sort_arr_j[%d] = %d \n", k, sort_arr_i[k], k, sort_arr_j[k]);
        if(sort_arr_i[k] != sort_arr_j[k]) {
          free(sort_arr_i);
          free(sort_arr_j);
          return 0;
        }
      }

      free(sort_arr_i);
      free(sort_arr_j);
    } else {
      return 0;
    }
  } else {
    return 0;
  }

  return 1;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
