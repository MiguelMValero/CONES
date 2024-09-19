/*
 * \file
 */

#ifndef __PDM_UNIQUE_H__
#define __PDM_UNIQUE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_sort.h"
#include "pdm_quick_sort.h"

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
 * \brief Unique in place
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
);

/**
 *
 * \brief Unique in place
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
);



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
);

/**
 *
 * \brief Unique in place
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
int
PDM_inplace_unique_long2
(
 PDM_g_num_t a[],
 int unique_order[],
 int l,
 int r
);


int
PDM_unique_long_with_distrib
(
  PDM_MPI_Comm   comm,
  PDM_g_num_t   *dentity1_entity2_gnum,
  PDM_g_num_t   *distrib_entity2,
  int            array_size,
  int          **unique_order,
  PDM_g_num_t  **unique_dentity1_entity2_gnum
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_SORT_H__ */
