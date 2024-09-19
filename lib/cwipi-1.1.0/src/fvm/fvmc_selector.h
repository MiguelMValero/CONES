#ifndef __FVMC_SELECTOR_H__
#define __FVMC_SELECTOR_H__

/*============================================================================
 * Mechanism for entity selection based on groups or attributes
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2007  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_group.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _fvmc_selector_t  fvmc_selector_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a selector object.
 *
 * parameters:
 *   dim                  <-- spatial dimension (coordinates and normals)
 *   n_elements           <-- number of selectable elements
 *   group_class_set      <-- pointer to group class set definition
 *   group_class_id       <-- group class id associated with each element
 *                            (size: n_elements)
 *   group_class_id_base; <-- Starting group class id base (usually 0 or 1)
 *   coords               <-- coordinates (interlaced) associated with each
 *                            element, whether vertex, face or cell center, ...
 *                            (size: n_elements * dim)
 *   normals              <-- normals (interlaced) associated with each element
 *                            if applicable (such as for face normals), or NULL
 *
 * returns:
 *   pointer to new selector
 *----------------------------------------------------------------------------*/

fvmc_selector_t *
fvmc_selector_create(int                           dim,
                    fvmc_lnum_t                    n_elements,
                    const fvmc_group_class_set_t  *group_class_set,
                    const int                     group_class_id[],
                    int                           group_class_id_base,
                    const double                  coords[],
                    const double                  normals[]);

/*----------------------------------------------------------------------------
 * Destruction of a selector structure.
 *
 * parameters:
 *   this_selector <-> selector to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_selector_t *
fvmc_selector_destroy(fvmc_selector_t  *this_selector);

/*----------------------------------------------------------------------------
 * Define the list of the elements verifying the criteria described
 * by a character string
 *
 * The selected_element[] array must be pre-allocated, and be of sufficient
 * size to contain all elements associated with the selector.
 *
 * parameters:
 *   this_selector       <-> pointer to selector
 *   str                 <-- string defining selection criteria
 *   n_selected_elements <-- number of elements selected
 *   selected_elements   <-> selected elements list (1 to n numbering)
 *
 * returns:
 *   criteria id associated by selector with str
 *----------------------------------------------------------------------------*/

int
fvmc_selector_get_list(fvmc_selector_t  *this_selector,
                      const char      *str,
                      fvmc_lnum_t      *n_selected_elements,
                      fvmc_lnum_t      *selected_elements);

/*----------------------------------------------------------------------------
 * Return the number of operands associated with a selection criteria
 * which are missing in the selector's associated group class set.
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *   criteria_id   <-- id of criteria returned by fvmc_selector_get_list()
 *
 * returns:
 *   number of missing operands
 *----------------------------------------------------------------------------*/

int
fvmc_selector_n_missing(const fvmc_selector_t  *this_selector,
                       int                    criteria_id);

/*----------------------------------------------------------------------------
 * Return a pointer to the name of an of operand associated with a selection
 * criteria which is missing in the selector's associated group class set.
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *   criteria_id   <-- id of criteria returned by fvmc_selector_get_list()
 *   missing_id    <-- id of missing operand for this criteria
 *
 * returns:
 *   pointer to name of missing operand
 *----------------------------------------------------------------------------*/

const char *
fvmc_selector_get_missing(const fvmc_selector_t  *this_selector,
                         int                    criteria_id,
                         int                    missing_id);

/*----------------------------------------------------------------------------
 * Dump the contents of a selector structure in human readable form
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *----------------------------------------------------------------------------*/

void
fvmc_selector_dump(const fvmc_selector_t  *this_selector);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_SELECTOR_H__ */
