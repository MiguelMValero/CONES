/*
 * \file
 */

#ifndef __PDM_BINARY_SEARCH_H__
#define __PDM_BINARY_SEARCH_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdint.h>

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
 const PDM_g_num_t  elt,
 const PDM_g_num_t *array,
 int         lArray
);

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
);



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
);


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
);


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
);


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
 );


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
);


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
PDM_search_rank
(
 PDM_g_num_t   elt,
 PDM_g_num_t  *array,
 int            id1,
 int            id2
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_BINARY_SEARCH_H__ */
