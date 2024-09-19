/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_printf.h"

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
 * \brief Search gap in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_gap_long
(
 const PDM_g_num_t   elt,
 const PDM_g_num_t  *array,
 const int          lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if ((elt >= array[ind]) && (elt < array[right])) {
    return ind;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search gap in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_gap_size_t
(
 const size_t  elt,
 const size_t *array,
 const int     lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if ((elt >= array[ind]) && (elt < array[right])) {
    return ind;
  }
  else {
    return -1;
  }
}

/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_long
(
 const PDM_g_num_t   elt,
 const PDM_g_num_t  *array,
 const int          lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if (elt == array[ind]) {
    return ind;
  }
  if (elt == array[right]) {
    return right;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search gap in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_gap_int
(
 const int   elt,
 const int  *array,
 const int   lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if ((elt >= array[ind]) && (elt < array[right])) {
    return ind;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_int
(
 const int          elt,
 const int         *array,
 const int          lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if (elt == array[ind]) {
    return ind;
  }
  if (elt == array[right]) {
    return right;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_uint32t
(
 const uint32_t     elt,
 const uint32_t    *array,
 const int          lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if (elt == array[ind]) {
    return ind;
  }
  if (elt == array[right]) {
    return right;
  }
  else {
    return -1;
  }
}


/**
 *
 * \brief Search gap in a sorted array
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   lArray       Array length
 *
 * \return       Index where element is stored
 */

int
PDM_binary_search_gap_double
(
 const double   elt,
 const double  *array,
 const int      lArray
)
{
  int left  = 0;
  int right = lArray - 1;
  int ind   = (left + right) / 2;

  while ((right - left) > 1) {

    if (elt < array[ind]) {
      right = ind;
    }
    else if (elt >= array[ind]) {
      left = ind;
    }

    ind = (left + right) / 2;

  }

  if ((elt >= array[ind]) && (elt < array[right])) {
    return ind;
  }
  else {
    return -1;
  }
}

/**
 *
 * \brief Search the rank where element of distributed array is storage
 *
 * \param [in]   elt          Element to find
 * \param [in]   array        Array where to search
 * \param [in]   id1          First index into array
 * \param [in]   id2          Last index into array
 *
 * \return       Rank where the element is stored
 */

int
PDM_search_rank
(
 PDM_g_num_t   elt,
 PDM_g_num_t  *array,
 int            id1,
 int            id2
)
{
  if (elt >= array[id2]) {
    PDM_printf("PDM_search_rank error : Element not in initial distributed array "
           PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
           elt, array[id1], array[id2]);
    abort();
  }

  if (elt < array[id1]) {
    PDM_printf("PDM_search_rank error : Element not in initial distributed array "
           PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
           elt, array[id1], array[id2]);
    abort();
  }

  if (id2 == id1 + 1) {
    return id1;
  }

  else {

    while(array[id1] == array[id1+1]) id1++;

    int midId = (id2 + id1) / 2;

    if (elt == array[id1])
      return id1;
    else if (elt < array[midId])
      return PDM_search_rank(elt, array, id1, midId);
    else if (elt >= array[midId])
      return PDM_search_rank(elt, array, midId, id2);
  }
  return -1;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */


