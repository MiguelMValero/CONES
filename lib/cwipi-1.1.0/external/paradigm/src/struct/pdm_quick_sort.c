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
#include "pdm_quick_sort.h"

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
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
void
PDM_quick_sort_long
(
 PDM_g_num_t a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    PDM_g_num_t t;
    PDM_g_num_t pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    PDM_quick_sort_long(a, l  , j-1);
    PDM_quick_sort_long(a, j+1,   r);
  }
}


/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
void
PDM_quick_sort_int
(
 int a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    int t;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    PDM_quick_sort_int(a, l  , j-1);
    PDM_quick_sort_int(a, j+1,   r);
  }
}



/**
 *
 * \brief Quick sort 2
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 * \param [inout]   c     Array sorted as a
 *
 */
void
PDM_quick_sort_int2
(
 int          a[],
 int          l,
 int          r,
 int          c[]
)
{
  if (l < r) {
    int j = r+1;
    int  t, v;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

      v    = c[i];
      c[i] = c[j];
      c[j] = v;
    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    v    = c[l];
    c[l] = c[j];
    c[j] = v;

    PDM_quick_sort_int2(a, l  , j-1, c);
    PDM_quick_sort_int2(a, j+1,   r, c);
  }
}

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 * \param [inout]   c     Array sorted as a
 *
 */

void
PDM_quick_sort_long2
(
 PDM_g_num_t a[],
 int          l,
 int          r,
 int          c[]
)
{
  if (l < r) {
    int j = r+1;

    PDM_g_num_t t;
    int v;
    PDM_g_num_t pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

      v    = c[i];
      c[i] = c[j];
      c[j] = v;
    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    v    = c[l];
    c[l] = c[j];
    c[j] = v;

    PDM_quick_sort_long2(a, l  , j-1, c);
    PDM_quick_sort_long2(a, j+1,   r, c);
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
