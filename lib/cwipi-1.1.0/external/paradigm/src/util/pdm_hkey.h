/*
 * \file
 */

#ifndef __PDM_HKEY_H__
#define __PDM_HKEY_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*----------------------------------------------------------------------------*/

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
 * Public function prototypes
 *============================================================================*/


/**
 *
 * \brief Build a hash from n int components
 *
 * \param [in]  n_key          Number of keys
 * \param [in]  components_idx Index of components in components (size = n_component)
 * \param [in]  components     List of components
 * \param [in]  max_key        Max authorized key 
 * \param [out] hkeys          Hash keys 
 *
 */

void
PDM_hkey_int 
(
  const int  n_key, 
  const int *components_idx,  
  const int *components,
  const int  max_key,
        int *keys
);


/**
 *
 * \brief Build a hash from n gnum components
 *
 * \param [in]  n_key          Number of keys
 * \param [in]  components_idx Index of components in components (size = n_component)
 * \param [in]  components     List of components
 * \param [in]  max_key        Max authorized key
 * \param [out] hkeys          Hash keys 
 *
 */

void
PDM_hkey_gnum 
(
  const int          n_key, 
  const int         *components_idx,  
  const PDM_g_num_t *components,
  const PDM_g_num_t  max_key,
        PDM_g_num_t *keys
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HKEY_H__ */
