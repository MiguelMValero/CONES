/*
 * \file
 */

#ifndef __PDM_HO_ORDERING_H_
#define __PDM_HO_ORDERING_H_


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_hash_tab.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/


/*=============================================================================
 * Global variables
 *============================================================================*/


/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Initialize the structure that stores HO orderings
 */

void
PDM_ho_ordering_init
(
void
);


/**
 * \brief Free the structure that stores HO orderings
 *
 * This function is automatically called upon exit of
 * the program which called \ref PDM_ho_ordering_init.
 *
 */

void
PDM_ho_ordering_free
(
void
);


/**
 * \brief Add a user-defined HO ordering from the locations
 * in the reference uvw-grid of a given element type of a given order
 *
 * \param[in] name         Name of the HO ordering
 * \param[in] t_elt        Element type
 * \param[in] order        Element order
 * \param[in] n_nodes      Number of nodes in the high-order element
 * \param[in] user_to_ijk  IJK-coordinates of the nodes in the high-order element
 *
 * \return                 Id of the HO ordering
 */

int
PDM_ho_ordering_user_to_ijk_add
(
 const char                 *name,
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                   n_nodes,
 const int                  *user_to_ijk
);


/**
 * \brief Get the node locations in the reference uvw-grid
 * for a user-defined HO ordering of a given element type of a given order
 *
 * \param[in] name         Name of the HO ordering
 * \param[in] t_elt        Element type
 * \param[in] order        Element order
 *
 * \return                 IJK-coordinates of the nodes in the high-order element
 */

int *
PDM_ho_ordering_user_to_ijk_get
(
 const char                 *name,
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order
);


/**
 * \brief Get the map from ijk HO ordering to a user-defined HO ordering
 * for a given element type of a given order
 *
 * \param[in] name         Name of the HO ordering
 * \param[in] t_elt        Element type
 * \param[in] order        Element order
 *
 * \return                 User-defined ordering of the nodes in the high-order element
 */

int *
PDM_ho_ordering_ijk_to_user_get
(
 const char                 *name,
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order
);


/**
 * \brief Compute the map from ijk HO ordering to a user-defined HO ordering
 * for a given element type of a given order
 *
 * \param[in] name         Name of the HO ordering
 * \param[in] t_elt        Element type
 * \param[in] order        Element order
 *
 * \return                 User-defined ordering of the nodes in the high-order element
 */

int *
PDM_ho_ordering_compute_ijk_to_user
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                  *user_to_ijk
);

/**
 * \brief Return the ID of an HO ordering
 *
 * \param[in] name         Name of the HO ordering
 *
 * \return   ID of the HO ordering (-1 if not found)
 */

int
PDM_ho_ordering_id_get
(
 const char *name
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_ORDERING_H__ */
