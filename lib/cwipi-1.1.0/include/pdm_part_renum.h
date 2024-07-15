/*
 * \file
 */

#ifndef __PDM_PART_RENUM_H__
#define	__PDM_PART_RENUM_H__

/*============================================================================
 * Mesh entities renumbering
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/**
 * \struct PDM_part_renum_fct_t
 *
 * \brief  Function pointer used to define a renumbering method
 *
 */

typedef void (*PDM_part_renum_fct_t) (part_t  **ppart, int n_part, void* specific_data);

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Add a new method for cell renumbering
 *
 * \param [in]      name           Mesh entity to renumber
 * \param [in]      renum_fct      Renumbering function
 *
 */

int
PDM_part_renum_method_cell_add
(
 const char                 *name,     /*!< Name          */
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_cell function for the format */
);

/**
 *
 * \brief Add a new method for face renumbering
 *
 * \param [in]      name           Mesh entity to renumber
 * \param [in]      renum_fct      Renumbering function
 *
 */

int
PDM_part_renum_method_face_add
(
 const char                 *name,     /*!< Name          */
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_face function for the format */
);

/**
 *
 * \brief Add a new method for edge renumbering
 *
 * \param [in]      name           Mesh entity to renumber
 * \param [in]      renum_fct      Renumbering function
 *
 */

int
PDM_part_renum_method_edge_add
(
 const char                 *name,     /*!< Name          */
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_face function for the format */
);


/**
 *
 * \brief Add a new method for face renumbering
 *
 * \param [in]      name           Mesh entity to renumber
 * \param [in]      renum_fct      Renumbering function
 *
 */

int
PDM_part_renum_method_vtx_add
(
 const char                 *name,     /*!< Name          */
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_face function for the format */
);

/**
 *
 * \brief Get index of a renumbering cell method
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

int
PDM_part_renum_method_cell_idx_get
(
const char *name
);


/**
 *
 * \brief Get index of a renumbering face method
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

int
PDM_part_renum_method_face_idx_get
(
const char *name
);

/**
 *
 * \brief Get index of a renumbering edge method
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

int
PDM_part_renum_method_edge_idx_get
(
const char *name
);

/**
 *
 * \brief Get index of a renumbering vtx method
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

int
PDM_part_renum_method_vtx_idx_get
(
const char *name
);



/**
 *
 * \brief Get name of the cell renumbering method
 *
 * \param [in]  idx     Index of the method
 *
 * \return Name of the method
 *
 */

void
PDM_part_renum_method_cell_name_get_cf
(
 const int  idx,
 char      *name,
 int       *l_name
 );

const char *
PDM_part_renum_method_cell_name_get
(
const int idx
);


/**
 *
 * \brief Get name of the face renumbering method
 *
 * \param [in]  idx     Index of the method
 *
 * \return Name of the method
 *
 */

void
PDM_part_renum_method_face_name_get_cf
(
 const int  idx,
 char      *name,
 int       *l_name
 );

const char *
PDM_part_renum_method_face_name_get
(
const int idx
);


/**
 *
 * \brief Get the number of renumbering face methods
 *
 * \return Name of the method
 *
 */

int
PDM_part_n_renum_method_cell_get
(
void
);


/**
 *
 * \brief Get the number of renumbering face methods
 *
 * \return Number of methods
 *
 */

int
PDM_part_n_renum_method_face_get
(
void
);


/**
 *
 * \brief Get the number of renumbering face methods
 *
 */

int
PDM_part_n_renum_method_cell_get
(
void
);


/**
 *
 * \brief Get the number of renumbering face methods
 *
 * \return Name of the method
 *
 */

int
PDM_part_n_renum_method_face_get
(
void
);


/**
 *
 * \brief Purge renumbering methods
 *
 */

void
PDM_part_renum_method_purge
(
void
);

/**
 *
 * \brief Load local renumbering methods
 *
 */

void
PDM_part_renum_method_load_local
(
void
);


/**
 *
 * \brief Perform cell renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void
PDM_part_renum_cell
(
 part_t **part,
 int       n_part,
 int       renum_cell_method,
 void     *specific_data
);


/**
 *
 * \brief Perform face renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void
PDM_part_renum_face
(
 part_t **part,
 int       n_part,
 int       renum_face_method,
 void     *specific_data
);


/**
 *
 * \brief Perform vtx renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void
PDM_part_renum_edge
(
 part_t **part,
 int       n_part,
 int       renum_vtx_method,
 void     *specific_data
);

/**
 *
 * \brief Perform vtx renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void
PDM_part_renum_vtx
(
 part_t **part,
 int       n_part,
 int       renum_vtx_method,
 void     *specific_data
);


/**
 *
 * \brief Perform cells renumbering from a new order
 *        Actualise all cells array according to the new numbering
 *        Connectivities/cell_tag/cell_color/cell_ln_to_gn
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */

void
PDM_part_reorder_cell
(
 part_t *part,
 int     *new_to_old_order
);


/**
 *
 * \brief Perform faces renumbering from a new order
 *        Actualise all cells array according to the new numbering
 *        Connectivities/face_tag/face_color/face_ln_to_gn
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */

void
PDM_part_reorder_face
(
 part_t *part,
 int     *new_to_old_order
);

/**
 *
 * \brief Perform vtx renumbering from a new order
 *        Actualise all faces array according to the new numbering
 *        Connectivities/face_tag/face_color/face_ln_to_gn
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */

void
PDM_part_reorder_edge
(
 part_t *part,
 int     *new_to_old_order
);

/**
 *
 * \brief Perform vtx renumbering from a new order
 *        Actualise all faces array according to the new numbering
 *        Connectivities/face_tag/face_color/face_ln_to_gn
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */

void
PDM_part_reorder_vtx
(
 part_t *part,
 int     *new_to_old_order
);


/**
 *
 * \brief Perform faces renumbering from a new order
 *        Actualise all cells array according to the new numbering
 *        Connectivities/face_tag/face_color/face_ln_to_gn
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */

void
PDM_part_renum_connectivities
(
  const int nElt,
  const int *new_to_old_order,
  int       *connectivity_idx,
  int       *connectivities
);

/**
 * \brief Order face_cell or edge_vtx array
 *
 * \param [in]      n_face              Number of elements
 * \param [in]      new_to_old_order    New order (size = \ref nElt
 * \param [in, out] face_cell           Connectivity ( size = 2 * n_face)
 *
 */
void
PDM_order_face_cell
(
int          n_face,
int         *new_to_old_order,
int         *face_cell
);

/**
 * \brief Order an array
 *
 * \param [in]      sizeArray       Number of elements
 * \param [in]      new_to_old_order        New order (size = \ref nElt
 * \param [in, out] Array           Array to renumber
 *
 */

void
PDM_part_renum_array
(
const int  sizeArray,
const int *old_to_new_order,
int       *array
);

/**
 * \brief Order an array
 *
 * \param [in]      sizeArray       Number of elements
 * \param [in]      new_to_old_order        New order (size = \ref n_elmt
 * \param [in, out] Array           Array to renumber
 *
 */

void
PDM_part_renum_graph
(
const int   n_entity1,
      int  *entity1_entity1_idx,
      int  *entity1_entity1,
const int  *new_to_old_order,
int         start
);

/**
 * \brief Order an array for face_cell
 *
 * \param [in]      sizeArray       Number of elements
 * \param [in]      new_to_old_order        New order (size = \ref nElt
 * \param [in, out] Array           Array to renumber
 *
 */

void
PDM_part_renum_array_face_cell
(
const int  sizeArray,
const int *old_to_new_order,
int       *array
);

/**
 *
 * \brief Perform faces renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */
void
PDM_part_reorder_face_bound
(
 part_t *part,
 int    *new_to_old_order
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_RENUM_H */

