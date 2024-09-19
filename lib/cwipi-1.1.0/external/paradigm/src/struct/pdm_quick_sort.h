/*
 * \file
 */

#ifndef __PDM_QUICK_SORT_H__
#define __PDM_QUICK_SORT_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
);


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
);

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
);

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
PDM_quick_sort_int2
(
 int          a[],
 int          l,
 int          r,
 int          c[]
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_QUICK_SORT_H__ */
