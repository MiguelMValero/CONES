/*
 * \file
 */

#ifndef __PDM_SORT_H__
#define __PDM_SORT_H__

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
typedef int (*pdm_operator_compare)(const void* a, const void* b, void* );

typedef struct  {
  int           *idx;
  unsigned char *arr; /* Not void* to avoid arithmetic pointer */
  size_t        *key;
} PDM_user_defined_sort;

// typedef struct  {
//   int*  idx;
//   unsigned char* arr;
//   size_t*       hkey;
// } PDM_user_defined_sort2;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout] array        Array to sort
 * \param [inout] order        new indice to old indice  (or NULL)
 * \param [in]    lArray       Array length
 *
 */

void
PDM_sort_long
(
 PDM_g_num_t  *array,
 int         *order,
 int          lArray
);

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
PDM_sort_int
(
 int         *array,
 int         *order,
 int          lArray
);

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
PDM_sort_double
(
 double     *array,
 int        *order,
 int         lArray
);


/**
 *
 * \brief Compare operator for connectivities
 *
 */
int
PDM_operator_compare_string
(
const void* a,
const void* b,
      void* ctxt
);

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
      void* ctxt
);

/**
 *
 * \brief Compare operator for connectivities
 *
 */
int
PDM_operator_compare_connectivity
(
const void* a,
const void* b,
      void* ctxt
);

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
      void* ctxt
);


/**
 *
 * \brief Indirect sort of structure according to a user comparison function
 *
 * \param [inout] array        Array to sort
 * \param [in]    lArray       Array length
 * \param [in]    lArray       User compare function (return int )
 * \param [in]    lArray       Context anonymous pointer to a struct to perfomr user compare
 *
 */
void
PDM_sort_int_special
(
 int                  *array,
 int                   lArray,
 pdm_operator_compare  comp,
 void*                 context
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_SORT_H__ */
