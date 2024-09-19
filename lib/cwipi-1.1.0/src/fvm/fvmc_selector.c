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

  This library is distibuted in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_error.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_config_defs.h"
#include "fvmc_group.h"
#include "fvmc_selector_postfix.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_selector.h"

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
 * Local type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local structure types
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Description (caching) of all previously interpreted criteria
 *----------------------------------------------------------------------------*/

typedef struct {

  int     n_operations;                /* Number of previously interpreted
                                          operations */
  int     n_max_operations;            /* Maximum size before reallocation */

  fvmc_selector_postfix_t **postfix;    /* Array of postfix operations */

  size_t *n_calls;                     /* Number of calls per operation */

  int    *n_group_classes;             /* Array of group class numbers
                                          for each operation */
  int   **group_class_set;             /* Array of group class lists
                                          for each operation */
} _operation_list_t;

/*----------------------------------------------------------------------------
 * Opaque type for management of all necessary info for elements selection
 *----------------------------------------------------------------------------*/

struct _fvmc_selector_t {

  int                dim;                      /* Spatial dimension */
  fvmc_lnum_t         n_elements;               /* Number of elements */

  const int         *group_class_id;           /* Element group class ids */
  int               *_group_class_id;          /* private group_class_id,
                                                  or NULL */
  int                group_class_id_base;      /* Starting base (usually
                                                  0 or 1) of group class ids */

  int                n_group_classes;          /* Number of group classes */

  int                n_groups;                 /* Total number of groups */
  int                n_attributes;             /* Total number of attributes */

  char             **group_name;               /* Ordered group names */
  int               *attribute;                /* Ordered attributes */

  int               *n_class_groups;           /* Number of groups per class */
  int              **group_ids;                /* Id of groups per class in
                                                  group_name */

  int               *n_class_attributes;       /* Number of attrs. per class */
  int              **attribute_ids;            /* Id of attributes per class in
                                                  attribute */

  const double      *coords;                   /* Element coordinates
                                                  (i.e. centers), interlaced */
  double            *_coords;                  /* private coords, or NULL */

  const double      *normals;                  /* Element normals, interlaced */
  double            *_normals;                 /* private normals, or NULL */

  _operation_list_t *_operations;              /* Array which caches all
                                                  previously interpreted
                                                  strings (operations) */

  fvmc_lnum_t        *_n_group_class_elements;  /* Number of elements per
                                                  group class */
  fvmc_lnum_t        **_group_class_elements;   /* Group class elements array */
};

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compare groups (qsort wrapper to strcmp).
 *
 * parameters:
 *   x <-> pointer to first group name
 *   y <-> pointer to second group name
 *
 * returns:
 *   result of strcmp() on group names
 *----------------------------------------------------------------------------*/

static int _compare_groups(const void *x, const void *y)
{
  return strcmp(*(char * const *)x, *(char * const *)y);
}

/*----------------------------------------------------------------------------
 * Add group information to a selector.
 *
 * A global list of group names from all group classes is built and sorted,
 * and the id's of groups of each class are determined.
 *
 * parameters:
 *   this_selector   <-> pointer to selector
 *   group_class_set <-- group classes description
 *----------------------------------------------------------------------------*/

static void
_assign_groups(fvmc_selector_t               *this_selector,
               const fvmc_group_class_set_t  *group_class_set)
{
  int i, j;
  char **_set_group_names = NULL;
  const char **set_group_names = NULL;
  int n_groups_tot = 0;

  /* Build basic arrays which exist in any case */

  BFTC_MALLOC(this_selector->n_class_groups,
             this_selector->n_group_classes,
             int);
  BFTC_MALLOC(this_selector->group_ids, this_selector->n_group_classes, int *);

  for (i = 0; i < this_selector->n_group_classes; i++) {

    const fvmc_group_class_t *gc = fvmc_group_class_set_get(group_class_set, i);
    const int n_groups = fvmc_group_class_get_n_groups(gc);

    n_groups_tot += n_groups;

    this_selector->n_class_groups[i] = n_groups;
    if (n_groups > 0)
      BFTC_MALLOC(this_selector->group_ids[i], n_groups, int);
    else
      this_selector->group_ids[i] = NULL;
  }

  if (n_groups_tot == 0)
    return;

  /* Fill arrays with unsorted group class info */

  BFTC_MALLOC(_set_group_names, n_groups_tot, char *);
  set_group_names = (const char **)_set_group_names;

  n_groups_tot = 0;

  for (i = 0; i < this_selector->n_group_classes; i++) {

    const fvmc_group_class_t *gc = fvmc_group_class_set_get(group_class_set, i);
    const int n_groups = fvmc_group_class_get_n_groups(gc);
    const char **group_names = fvmc_group_class_get_group_names(gc);

    for (j = 0; j < n_groups; j++)
      set_group_names[n_groups_tot + j] = group_names[j];

    n_groups_tot += n_groups;
  }

  qsort(set_group_names, n_groups_tot, sizeof(char *), &_compare_groups);

  BFTC_MALLOC(this_selector->group_name, n_groups_tot, char *);

  BFTC_MALLOC(this_selector->group_name[0],
             strlen(set_group_names[0]) + 1,
             char);
  strcpy(this_selector->group_name[0], set_group_names[0]);
  for (i = 1, j = 1; i < n_groups_tot; i++) {
    const char *name = set_group_names[i];
    if (strcmp(name, set_group_names[i-1]) != 0) {
      BFTC_MALLOC(this_selector->group_name[j], strlen(name) + 1, char);
      strcpy(this_selector->group_name[j], name);
      j++;
    }
  }

  set_group_names = NULL;
  BFTC_FREE(_set_group_names);

  this_selector->n_groups = j;
  BFTC_REALLOC(this_selector->group_name, this_selector->n_groups, char *);

  /* Now build group_id arrays */

  for (i = 0; i < this_selector->n_group_classes; i++) {

    int mid_id;
    const char *name;
    const fvmc_group_class_t *gc = fvmc_group_class_set_get(group_class_set, i);
    const int n_groups = fvmc_group_class_get_n_groups(gc);
    const char **group_names = fvmc_group_class_get_group_names(gc);

    for (j = 0; j < n_groups; j++) {

      /* use binary search */

      int start_id = 0;
      int end_id = this_selector->n_groups - 1;
      mid_id = (end_id -start_id) / 2;
      name = group_names[j];

      while (start_id <= end_id) {
        int att_cmp = strcmp(this_selector->group_name[mid_id], name);
        if (att_cmp < 0)
          start_id = mid_id + 1;
        else if (att_cmp > 0)
          end_id = mid_id - 1;
        else
          break;
        mid_id = start_id + ((end_id -start_id) / 2);
      }

      assert(strcmp(this_selector->group_name[mid_id], name) == 0);
      this_selector->group_ids[i][j] = mid_id;
    }
  }
}

/*----------------------------------------------------------------------------
 * Compare attributes (for qsort).
 *
 * parameters:
 *   x <-> pointer to first attribute
 *   y <-> pointer to second attribute
 *
 * returns:
 *   result of strcmp() on group names
 *----------------------------------------------------------------------------*/

static int _compare_attributes(const void *x, const void *y)
{
  return (*(const int *)x - *(const int *)y);
}

/*----------------------------------------------------------------------------
 * Add attribute information to a selector.
 *
 * A global list of attributes from all group classes is built and sorted,
 * and the id's of attributes of each class are determined.
 *
 * parameters:
 *   this_selector   <-> pointer to selector
 *   group_class_set <-- group classes description
 *----------------------------------------------------------------------------*/

static void
_assign_attributes(fvmc_selector_t               *this_selector,
                   const fvmc_group_class_set_t  *group_class_set)
{
  int i, j;
  int *set_attributes = NULL;
  int n_attributes_tot = 0;

  /* Build basic arrays which exist in any case */

  BFTC_MALLOC(this_selector->n_class_attributes,
             this_selector->n_group_classes,
             int);
  BFTC_MALLOC(this_selector->attribute_ids,
             this_selector->n_group_classes,
             int *);

  for (i = 0; i < this_selector->n_group_classes; i++) {

    const fvmc_group_class_t *gc = fvmc_group_class_set_get(group_class_set, i);
    const int n_attributes = fvmc_group_class_get_n_attributes(gc);

    n_attributes_tot += n_attributes;

    this_selector->n_class_attributes[i] = n_attributes;
    if (n_attributes > 0)
      BFTC_MALLOC(this_selector->attribute_ids[i], n_attributes, int);
    else
      this_selector->attribute_ids[i] = NULL;
  }

  if (n_attributes_tot == 0)
    return;

  /* Fill arrays with unsorted group class info */

  for (i = 0; i < this_selector->n_group_classes; i++) {
    const fvmc_group_class_t *gc = fvmc_group_class_set_get(group_class_set, i);
    n_attributes_tot += fvmc_group_class_get_n_attributes(gc);
  }

  if (n_attributes_tot == 0)
    return;

  BFTC_MALLOC(set_attributes, n_attributes_tot, int);

  n_attributes_tot = 0;

  for (i = 0; i < this_selector->n_group_classes; i++) {

    const fvmc_group_class_t *gc = fvmc_group_class_set_get(group_class_set, i);
    const int n_attributes = fvmc_group_class_get_n_attributes(gc);
    const int *attributes = fvmc_group_class_get_attributes(gc);

    for (j = 0; j < n_attributes; j++)
      set_attributes[n_attributes_tot + j] = attributes[j];

    n_attributes_tot += n_attributes;
  }

  qsort(set_attributes, n_attributes_tot, sizeof(int), &_compare_attributes);

  BFTC_MALLOC(this_selector->attribute, n_attributes_tot, int);

  this_selector->attribute[0] = set_attributes[0];
  for (i = 1, j = 1; i < n_attributes_tot; i++) {
    if (set_attributes[i] != set_attributes[i-1])
      this_selector->attribute[j++] = set_attributes[i];
  }

  BFTC_FREE(set_attributes);

  this_selector->n_attributes = j;
  BFTC_REALLOC(this_selector->attribute, this_selector->n_attributes, int);

  /* Now build attribute_id arrays */

  for (i = 0; i < this_selector->n_group_classes; i++) {

    int mid_id;
    const fvmc_group_class_t *gc = fvmc_group_class_set_get(group_class_set, i);
    const int n_attributes = fvmc_group_class_get_n_attributes(gc);
    const int *attributes = fvmc_group_class_get_attributes(gc);

    for (j = 0; j < n_attributes; j++) {

      /* use binary search */

      int start_id = 0;
      int end_id = this_selector->n_attributes - 1;
      int val = attributes[j];
      mid_id = (end_id -start_id) / 2;

      /* use binary search */

      while (start_id <= end_id) {
        int att_cmp = this_selector->attribute[mid_id];
        if (att_cmp < val)
          start_id = mid_id + 1;
        else if (att_cmp > val)
          end_id = mid_id - 1;
        else
          break;
        mid_id = start_id + ((end_id -start_id) / 2);
      }

      assert(this_selector->attribute[mid_id] == val);
      this_selector->attribute_ids[i][j] = mid_id;

    }
  }
}

/*----------------------------------------------------------------------------
 * Create an operations list.
 *
 * The default number of operations is 30.
 *
 * returns:
 *   pointer to the operations list
 *----------------------------------------------------------------------------*/

static _operation_list_t *
_operation_list_allocate(void)
{
  _operation_list_t *ops;
  int i;
  const int n_operations = 16;

  BFTC_MALLOC(ops, 1, _operation_list_t);

  /*  Definitions */

  ops->n_operations = 0;
  ops->n_max_operations = n_operations;

  BFTC_MALLOC(ops->postfix,
             ops->n_max_operations,
             fvmc_selector_postfix_t *);

  BFTC_MALLOC(ops->n_calls, ops->n_max_operations, size_t);

  BFTC_MALLOC(ops->n_group_classes, ops->n_max_operations, int);
  BFTC_MALLOC(ops->group_class_set, ops->n_max_operations, int *);

  for (i = 0; i < ops->n_max_operations; i++) {
    ops->postfix[i] = NULL;
    ops->group_class_set[i] = NULL;
    ops->n_calls[i] = 0;
    ops->n_group_classes[i] = 0;
  }

  return ops;
}

/*----------------------------------------------------------------------------
 * Increase operations list length.
 *
 * Length is multiplied by 2.
 *
 * parameters:
 *   ops <-> operations list to be updated
 *----------------------------------------------------------------------------*/

static void
_operation_list_reallocate(_operation_list_t  *ops)
{
  int old_size;

  int i = 0;

  old_size = ops->n_max_operations;
  ops->n_max_operations *= 2;

  /* Reallocation */

  BFTC_REALLOC(ops->postfix,
              ops->n_max_operations,
              fvmc_selector_postfix_t *);

  BFTC_REALLOC(ops->n_calls, ops->n_max_operations, size_t);

  BFTC_REALLOC(ops->n_group_classes, ops->n_max_operations, int);
  BFTC_REALLOC(ops->group_class_set, ops->n_max_operations, int *);

  for (i = old_size; i < ops->n_max_operations; i++) {
    ops->postfix[i] = NULL;
    ops->group_class_set[i] = NULL;
    ops->n_calls[i] = 0;
    ops->n_group_classes[i] = 0;
  }
}

/*----------------------------------------------------------------------------
 * Delete operations list.
 *
 * parameters:
 *   ops <-> operations list to be deleted
 *
 * returns:
 *   NULL pointer;
 *----------------------------------------------------------------------------*/

static _operation_list_t *
_operation_list_free(_operation_list_t  *ops)
{
  int i = 0;

  if (ops != NULL) {
    BFTC_FREE(ops->n_calls);
    BFTC_FREE(ops->n_group_classes);
    for (i = 0; i < ops->n_max_operations; i++) {
      if (ops->group_class_set[i] != NULL)
        BFTC_FREE(ops->group_class_set[i]);
      if (ops->postfix[i] != NULL)
        fvmc_selector_postfix_destroy(ops->postfix + i);
    }
    BFTC_FREE(ops->postfix);
    BFTC_FREE(ops->group_class_set);
    BFTC_FREE(ops);
  }

  return NULL;
}

/*----------------------------------------------------------------------------
 * Interpret the postfix string for the last operation of operations list to
 * build the group class list which corresponds to this operation.
 *
 * parameters:
 *   this_selector <-- selector structure
 *   operations    <-> operations list to be updated
 *----------------------------------------------------------------------------*/

static void
_create_operation_group_class_set(const fvmc_selector_t  *this_selector,
                                  _operation_list_t     *operations)
{
  int gc_id;
  int *group_class_set;

  int n_group_classes = 0;

  const fvmc_selector_postfix_t  *pf
    = operations->postfix[operations->n_operations -1];

  BFTC_MALLOC(operations->group_class_set[operations->n_operations - 1],
             this_selector->n_group_classes,
             int);

  group_class_set
    = operations->group_class_set[operations->n_operations - 1];

  for (gc_id = 0; gc_id < this_selector->n_group_classes; gc_id++) {

    /* update group class list for current operation */

    if (fvmc_selector_postfix_eval(pf,
                                  this_selector->n_class_groups[gc_id],
                                  this_selector->n_class_attributes[gc_id],
                                  this_selector->group_ids[gc_id],
                                  this_selector->attribute_ids[gc_id],
                                  NULL,
                                  NULL))
      group_class_set[n_group_classes++] = gc_id;
  }

  operations->n_group_classes[operations->n_operations-1] = n_group_classes;

  BFTC_REALLOC(operations->group_class_set[operations->n_operations-1],
              n_group_classes,
              int);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("  - group_classes list: ");
  {
    int i ;
    /* fvmc_group_class_set_dump(class_defs); */
    for (i = 0; i < n_group_classes; i++)
      bftc_printf("%i ",
                 operations->group_class_set[operations->n_operations-1][i]);
    bftc_printf("\n");
  }
#endif /* defined(DEBUG) && !defined(NDEBUG) */
}

/*----------------------------------------------------------------------------
 * Add a new operation in the operation_list :
 *
 * - infix string is parsed into postfix string
 * - interpret the postfix string
 *
 * parameters:
 *   selector     <-> selector to be updated
 *   infix_string <-- string parsed
 *----------------------------------------------------------------------------*/

static void
_add_new_operation(fvmc_selector_t  *selector,
                   const char      *infix_string)
{
  fvmc_selector_postfix_t *pf = NULL;

  /* reallocation  */

  if (   selector->_operations->n_operations
      >= selector->_operations->n_max_operations)
    _operation_list_reallocate(selector->_operations);

  /* Parse infix_string */

  pf = fvmc_selector_postfix_create(infix_string,
                                   selector->n_groups,
                                   selector->n_attributes,
                                   (const char **)selector->group_name,
                                   selector->attribute);

  /* update n_operations */

  selector->_operations->postfix[selector->_operations->n_operations] = pf;
  selector->_operations->n_operations++;


  /* Create group class list if there are no numerical tests in postfix */

  if (   fvmc_selector_postfix_coords_dep(pf) == false
      && fvmc_selector_postfix_normals_dep(pf) == false)
    _create_operation_group_class_set(selector,
                                      selector->_operations);
}

/*----------------------------------------------------------------------------
 * Get the test number in the operation_list which corresponds to
 * the infix string "test_str". If this string doesn't correspond to a an
 * operation already parsed, a new operation is added.
 *
 * parameters:
 *   teststr  <-- string parsed
 *   selector <-> current selector
 *
 * returns:
 *   the number of the operation which corresponds to "test_str" in the
 *   "selector" operations list
 *----------------------------------------------------------------------------*/

static int
_get_criteria_id(fvmc_selector_t  *selector,
                 const char      *teststr)
{
  int op = 0;

  assert(teststr != NULL);

  /* Search for teststr in the operations list */

  if (selector->_operations == NULL)
    selector->_operations = _operation_list_allocate();

  for (op = 0; op < selector->_operations->n_operations; op++) {
    const fvmc_selector_postfix_t *pf = selector->_operations->postfix[op];
    if (!strcmp(fvmc_selector_postfix_get_infix(pf), teststr))
      break;
  }

  /* if teststr is not in the list : add teststrcpy in the list */
  if (op == selector->_operations->n_operations)
    _add_new_operation(selector, teststr);

  return op;
}

/*============================================================================
 * Public function definitions
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
                    const double                  normals[])
{
  int i;
  fvmc_lnum_t j;
  fvmc_selector_t *selector;

  int n_group_classes = fvmc_group_class_set_size(group_class_set);

  BFTC_MALLOC(selector, 1, fvmc_selector_t);

  selector->dim = dim;
  selector->n_elements = n_elements;

  selector->group_class_id = group_class_id;
  selector->_group_class_id = NULL;
  selector->group_class_id_base = group_class_id_base;

  selector->n_group_classes = fvmc_group_class_set_size(group_class_set);

  selector->n_groups = 0;
  selector->n_attributes = 0;
  selector->group_name = NULL;
  selector->attribute = NULL;

  selector->n_class_groups = NULL;
  selector->group_ids = NULL;
  selector->n_class_attributes = NULL;
  selector->attribute_ids = NULL;

  selector->coords = coords;
  selector->_coords = NULL;
  selector->normals = normals;
  selector->_normals = NULL;

  selector->_operations = NULL;

  selector->_n_group_class_elements = NULL;
  selector->_group_class_elements = NULL;

  _assign_groups(selector, group_class_set);
  _assign_attributes(selector, group_class_set);

  if (group_class_id != NULL && n_group_classes > 0) {

    BFTC_MALLOC(selector->_n_group_class_elements, n_group_classes, fvmc_lnum_t);
    BFTC_MALLOC(selector->_group_class_elements, n_group_classes, fvmc_lnum_t *);

    /* Counting loop and allocation */

    for (i = 0; i < n_group_classes; i++)
      selector->_n_group_class_elements[i] = 0;

    for (j = 0; j < n_elements; j++)
      selector->_n_group_class_elements[  group_class_id[j]
                                        - group_class_id_base] += 1;

    for (i = 0; i < n_group_classes; i++)
      BFTC_MALLOC(selector->_group_class_elements[i],
                 selector->_n_group_class_elements[i],
                 int);

    /* Definition loop */

    for (i = 0; i < n_group_classes; i++)
      selector->_n_group_class_elements[i] = 0;

    for (j = 0; j < n_elements; j++) {
      selector->_group_class_elements
                   [group_class_id[j]-1]
                   [selector->_n_group_class_elements
                                [  group_class_id[j]
                                 - group_class_id_base]++]
        = j+1;

    }
  }

  return selector;
}

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
fvmc_selector_destroy(fvmc_selector_t  *this_selector)
{
  int i;

  /* Delete private arrays */

  _operation_list_free(this_selector->_operations);

  if (this_selector->_coords != NULL)
    BFTC_FREE(this_selector->_coords);
  if (this_selector->_normals != NULL)
    BFTC_FREE(this_selector->_normals);

  for (i = 0; i < this_selector->n_groups; i++)
    BFTC_FREE(this_selector->group_name[i]);
  BFTC_FREE(this_selector->group_name);

  BFTC_FREE(this_selector->attribute);

  BFTC_FREE(this_selector->n_class_groups);
  BFTC_FREE(this_selector->n_class_attributes);

  for (i = 0; i < this_selector->n_group_classes; i++) {
    if (this_selector->group_ids[i] != NULL)
      BFTC_FREE(this_selector->group_ids[i]);
    if (this_selector->attribute_ids[i] != NULL)
      BFTC_FREE(this_selector->attribute_ids[i]);
  }

  BFTC_FREE(this_selector->group_ids);
  BFTC_FREE(this_selector->attribute_ids);

  if (this_selector->_group_class_elements != NULL) {
    for (i = 0; i < this_selector->n_group_classes; i++)
      BFTC_FREE(this_selector->_group_class_elements[i]);

    BFTC_FREE(this_selector->_n_group_class_elements);
    BFTC_FREE(this_selector->_group_class_elements);
  }

  /* Delete selector */

  BFTC_FREE(this_selector);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Define the list of the elements verifying the criteria described
 * by a character string.
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
                      fvmc_lnum_t      *selected_elements)

{
  int  c_id, gc_id;
  const fvmc_selector_postfix_t *pf = NULL;
  fvmc_selector_t  *ts = this_selector;
  fvmc_lnum_t  i;

  assert(this_selector != NULL);

  *n_selected_elements = 0;

  /* Add or find the test number in the the cached operations list */

  c_id = _get_criteria_id(ts, str);

  ts->_operations->n_calls[c_id] += 1;
  pf = ts->_operations->postfix[c_id];

  /* Case without geometrical test: get group class list without the
     interpretration of the postfix writing of the criteria string */

  if (   fvmc_selector_postfix_coords_dep(pf) == false
      && fvmc_selector_postfix_normals_dep(pf) == false) {

    if (ts->_operations->group_class_set[c_id] != NULL) {

      int n_criteria_group_classes
        = ts->_operations->n_group_classes[c_id];
      const int *_criteria_group_class_set
        = ts->_operations->group_class_set[c_id];

      if (ts->_n_group_class_elements != NULL) {

        for (gc_id = 0; gc_id < n_criteria_group_classes; gc_id++) {
          for (i = 0;
               i < ts->_n_group_class_elements
                         [_criteria_group_class_set[gc_id]];
               i++) {
            selected_elements[(*n_selected_elements)++]
              = ts->_group_class_elements
                      [_criteria_group_class_set[gc_id]][i];
          }
        }

      }
    }

  }

  /* Case with geometrical test:
     evaluation of the postfix expression for each element */

  else if (ts->n_elements > 0) {

    const int dim = ts->dim;

    assert(ts->group_class_id != NULL);

    if (fvmc_selector_postfix_coords_dep(pf) == true && ts->coords == NULL)
      bftc_error(__FILE__, __LINE__, 0,
                _("Selection criteria:\n"
                  "\"%s\"\n"
                  "depends on coordinates, but the current selector\n"
                  "has no associated coordinates."),
                str);
    else if (   fvmc_selector_postfix_normals_dep(pf) == true
             && ts->normals == NULL)
      bftc_error(__FILE__, __LINE__, 0,
                _("Selection criteria:\n"
                  "\"%s\"\n"
                  "depends on normals, but the current selector\n"
                  "has no associated normals."),
                str);
    if (dim != 3)
        bftc_error(__FILE__, __LINE__, 0,
                  _("Selection criteria:\n"
                  "\"%s\"\n"
                    "is associated with %d spatial dimensions, but\n"
                    "geometric conditions are only currently implemented\n"
                    "for 3 spatial dimension."),
                  str, dim);

    /* Loop on all elements */

    for (i = 0; i < ts->n_elements; i++) {

      gc_id = ts->group_class_id[i] - ts->group_class_id_base;

      if (fvmc_selector_postfix_eval(pf,
                                    ts->n_class_groups[gc_id],
                                    ts->n_class_attributes[gc_id],
                                    ts->group_ids[gc_id],
                                    ts->attribute_ids[gc_id],
                                    ts->coords + (i*dim),
                                    ts->normals + (i*dim)))
        selected_elements[(*n_selected_elements)++] = i+1;

    }
  }

  return c_id;
}

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
                       int                    criteria_id)
{
  int retval = 0;

  if (this_selector != NULL && criteria_id >= 0) {
    if (   this_selector->_operations != NULL
        && this_selector->_operations->n_operations > criteria_id) {
      const fvmc_selector_postfix_t *pf
        = this_selector->_operations->postfix[criteria_id];

      retval = fvmc_selector_postfix_n_missing(pf);
    }
  }

  return retval;
}

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
                         int                    missing_id)
{
  const char *retval = NULL;

  if (this_selector != NULL && criteria_id >= 0) {
    if (   this_selector->_operations != NULL
        && this_selector->_operations->n_operations > criteria_id) {
      const fvmc_selector_postfix_t *pf
        = this_selector->_operations->postfix[criteria_id];

      retval = fvmc_selector_postfix_get_missing(pf, missing_id);
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Dump the contents of a selector structure in human readable form
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *----------------------------------------------------------------------------*/

void
fvmc_selector_dump(const fvmc_selector_t  *this_selector)
{
  int i, j;
  const fvmc_selector_t  *ts = this_selector;

  if (ts == NULL) {
    bftc_printf("\n"
               "Null selector dump:\n");
    return;
  }

  bftc_printf("\n"
             "Selector dump:\n"
             "  Dimension:                          %d\n"
             "  Number of selectable elements:      %d\n"
             "  Shared group class id's:            %p\n"
             "  Private group class id's:           %p\n"
             "  Group class id base:                %d\n"
             "  Number of associated group classes: %d\n"
             "  Number of associated groups:        %d\n"
             "  Number of associated attributes:    %d\n",
             ts->dim, (int)ts->n_elements,
             ts->group_class_id, ts->_group_class_id,
             ts->group_class_id_base,
             ts->n_group_classes, ts->n_groups, ts->n_attributes);

  if (ts->n_groups > 0) {
    bftc_printf("  Group names:\n");
    for (i = 0; i < ts->n_groups; i++)
      bftc_printf("    \"%s\"\n", ts->group_name[i]);
  }
  if (ts->n_attributes > 0) {
    bftc_printf("  Attributes:\n");
    for (i = 0; i < ts->n_attributes; i++)
      bftc_printf("    %d\n", ts->attribute[i]);
  }

  if (ts->n_group_classes > 0) {
    bftc_printf("  Group classes:\n");
    for (i = 0; i < ts->n_group_classes; i++) {
      bftc_printf("    Group class %d\n", (int)i);
      if (ts->n_groups > 0) {
        bftc_printf("      Number of groups: %d\n",
                   ts->n_class_groups[i]);
        for (j = 0; j < ts->n_class_groups[i]; j++)
          bftc_printf("        %d\n", ts->group_ids[i][j]);
      }
      if (ts->n_attributes > 0) {
        bftc_printf("      Number of attributes: %d\n",
                   ts->n_class_attributes[i]);
        for (j = 0; j < ts->n_class_attributes[i]; j++)
          bftc_printf("        %d\n", ts->attribute_ids[i][j]);
      }
    }
  }

  bftc_printf("  Shared coordinates:                 %p\n"
             "  Private coordinates:                %p\n"
             "  Shared normals;                     %p\n"
             "  Private normals:                    %p\n"
             "  Operations list:                    %p\n",
             ts->coords, ts->_coords, ts->normals, ts->_normals,
             ts->_operations);

  if (ts->n_group_classes > 0) {
    bftc_printf("  Number of elements per group class:\n");
    for (i = 0; i < ts->n_group_classes; i++) {
      bftc_printf("    %d (%p)\n",
                 (int)ts->_n_group_class_elements[i],
                 ts->_group_class_elements[i]);
    }
  }

  if (ts->_operations != NULL) {

    bftc_printf("\n");

    for (i = 0; i < ts->_operations->n_operations; i++) {
      bftc_printf("  Operation %d (cached, n_calls = %lu)\n",
                 i, (unsigned long) ts->_operations->n_calls[i]);
      fvmc_selector_postfix_dump(ts->_operations->postfix[i],
                                ts->n_groups, ts->n_attributes,
                                (const char **)ts->group_name,
                                ts->attribute);
    }

  }

  bftc_printf("\n");
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
