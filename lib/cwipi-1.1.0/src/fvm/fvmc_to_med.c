/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to MED files
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2007  EDF

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

/*----------------------------------------------------------------------------*/

#include "fvmc_config.h"
#include "config_priv.h"

/* We undefine the following macros to avoid conflicts with those defined
   by MED's config.h */

#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_VERSION
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef VERSION

/*----------------------------------------------------------------------------*/

#if defined(HAVE_MED)

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * MED library headers
 *----------------------------------------------------------------------------*/

#include <med.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_error.h>
#include <bftc_mem.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_convert_array.h"
#include "fvmc_gather.h"
#include "fvmc_io_num.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_parall.h"
#include "fvmc_writer_helper.h"
#include "fvmc_writer_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_to_med.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * MED field structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char    name[MED_TAILLE_NOM + 1];  /* MED field name */

  int     id;                        /* MED field id */
  int     n_components;              /* Number of components */
  med_type_champ  datatype;          /* Field datatype */

} fvmc_to_med_field_t;

/*----------------------------------------------------------------------------
 * MED mesh structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char    name[MED_TAILLE_NOM + 1]; /* Med_Mesh name */
  int     num;                      /* MED mesh number */

  med_int  entity_dim;  /* 3 for a volume mesh, 2 for a surface mesh */
  med_int  space_dim;   /* Number of coordinates to define a vertex */

} fvmc_to_med_mesh_t;

/*----------------------------------------------------------------------------
 * MED writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char       *name;          /* Writer name */
  char       *filename;      /* MED file name */
  med_idt     fid;           /* MED file id */

  int                  n_med_meshes;  /* Number of MED meshes */
  fvmc_to_med_mesh_t  **med_meshes;    /* Array of pointers to MED mesh
                                         structure */

  fvmc_writer_time_dep_t time_dependency; /* Mesh time dependency */

  int                  n_fields;      /* Number of fields */
  fvmc_to_med_field_t **fields;        /* Array of pointers to MED field
                                         structure */

  int         n_time_steps;    /* Number of meshes time steps */
  int        *time_steps;      /* Array of meshes time steps */
  double     *time_values;     /* Array of meshes time values */

  _Bool       is_open;            /* True if MED file is open, else false */

  _Bool       discard_polygons;   /* Option to discard polygonal elements */
  _Bool       discard_polyhedra;  /* Option to discard polyhedral elements */
  _Bool       divide_polygons;    /* Option to tesselate polygonal elements */
  _Bool       divide_polyhedra;   /* Option to tesselate polyhedral elements */

  int         rank;            /* Rank of current process in communicator */
  int         n_ranks;         /* Number of processes in communicator */

#if defined(FVMC_HAVE_MPI)
  MPI_Comm          comm;     /* Associated MPI communicator */
#endif

} fvmc_to_med_writer_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

#define FVMC_MED_MAX_N_NODES  8

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convert FVM float type into MED float datatype according to its storage.
 *
 * parameters:
 *   fvmc_data     <-- FVM data array to convert.
 *   fvmc_datatype <-- FVM datatype.
 *   med_data     --> MED data array converted.
 *   n_elems      <-- Number of elements to convert.
 *
 *----------------------------------------------------------------------------*/

static void
_convert_float_fvmc_to_med(const void      *fvmc_data,
                          fvmc_datatype_t   fvmc_datatype,
                          med_float       *med_data,
                          int              n_elems)
{
  int  i_elem;

  const float   *data_f = (const float *)fvmc_data;
  const double  *data_d = (const double *)fvmc_data;

  assert(fvmc_data != NULL);

  /* Type conversion adaptation */

  if (fvmc_datatype == FVMC_DOUBLE) {

    for (i_elem = 0; i_elem < n_elems; i_elem++)
      med_data[i_elem] = (med_float)data_d[i_elem];

  }
  else if (fvmc_datatype == FVMC_FLOAT)  {

    for (i_elem = 0; i_elem < n_elems; i_elem++)
      med_data[i_elem] = (med_float)data_f[i_elem];

  }
  else
    bftc_error(__FILE__, __LINE__, 0,
              "_convert_float_fvmc_to_med() incorrect datatype\n");

  return;
}

/*----------------------------------------------------------------------------
 * Convert FVM datatype fvmc_gnum_t into MED datatype med_int.
 *
 * parameters:
 *   fvmc_data    <-- FVM data array to convert.
 *   med_data    <-> MED data array converted.
 *   n_vals      <-- Number of values to convert.
 *
 * returns:
 *----------------------------------------------------------------------------*/

static void
_convert_fvmc_gnum_to_med_int(fvmc_gnum_t        *fvmc_data,
                             med_int           *med_data,
                             const fvmc_gnum_t   n_vals)
{
  fvmc_gnum_t i_val;

  for (i_val = 0; i_val < n_vals; i_val++) {

    if (sizeof(med_int) > sizeof(fvmc_gnum_t))
      med_data[n_vals - 1 - i_val] = (med_int)fvmc_data[n_vals - 1 -i_val];
    else
      med_data[i_val] = (med_int)fvmc_data[i_val];

  }

  return;
}

/*----------------------------------------------------------------------------
 * Return the MED mesh number associated with a given MED mesh name,
 * or 0 if no association found.
 *
 * parameters:
 *   writer         <-- MED writer structure
 *   med_mesh_name  <-- MED mesh name
 *
 * returns:
 *    MED mesh number, or 0 if MED mesh name is not associated with this
 *    MED writer structure in FVM
 *----------------------------------------------------------------------------*/

static int
_get_med_mesh_num(fvmc_to_med_writer_t  *writer,
                  const char           *med_mesh_name)
{
  int i;
  int retval = 0;

  fvmc_to_med_mesh_t  **med_meshes = writer->med_meshes;

  assert(writer != NULL);

  for (i = 0; i < writer->n_med_meshes; i++) {
    if (strcmp(med_mesh_name, med_meshes[i]->name) == 0)
      break;
  }

  if (i == writer->n_med_meshes)
    retval = 0;
  else
    retval = med_meshes[i]->num;

  return retval;
}

/*----------------------------------------------------------------------------
 * Create a MED mesh structure in MED writer structure.
 *
 * parameters:
 *   writer          <-- MED writer structure.
 *   med_mesh_name   <-- name in MED format of the mesh .
 *   mesh            <-- FVM mesh  structure.
 *
 * returns:
 *   MED mesh number, or 0 if any problem has been encountered.
 *----------------------------------------------------------------------------*/

static int
_add_med_mesh(fvmc_to_med_writer_t  *writer,
              char                 *med_mesh_name,
              const fvmc_nodal_t    *mesh)
{
  int id;
  int rank = writer->rank;
  char med_info[MED_TAILLE_DESC + 1] = "Generated by FVM.";
  char family_name[MED_TAILLE_NOM + 1] = "";

  med_err  retval = 0;

  assert(writer != NULL);

  /* Add a new MED mesh structure */

  writer->n_med_meshes += 1;
  id = writer->n_med_meshes - 1;

  BFTC_REALLOC(writer->med_meshes, writer->n_med_meshes, fvmc_to_med_mesh_t *);
  BFTC_MALLOC(writer->med_meshes[id], 1, fvmc_to_med_mesh_t);

  strncpy(writer->med_meshes[id]->name, med_mesh_name, MED_TAILLE_NOM);
  writer->med_meshes[id]->name[MED_TAILLE_NOM] = '\0';

  /* BUG in med: entity_dim not well treated.
     writer->med_meshes[id]->entity_dim = (med_int)entity_dim;
     So, this variable is set to space_dim */
  writer->med_meshes[id]->entity_dim = (med_int)mesh->dim;
  writer->med_meshes[id]->space_dim  = (med_int)mesh->dim;

  if (rank == 0) {

    retval = MEDmaaCr(writer->fid,
                      med_mesh_name,
                      writer->med_meshes[id]->entity_dim,
                      MED_NON_STRUCTURE,
                      med_info);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDmaaCr() failed to create a new med_mesh.\n"
                  "Associated med_mesh name: \"%s\"\n"
                  "Associated writer name: \"%s\"\n"),
                med_mesh_name, writer->name);

    retval = MEDdimEspaceCr(writer->fid,
                            med_mesh_name,
                            writer->med_meshes[id]->space_dim);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDdimEspaceCr() failed to set spacial dimension "
                  "for a new med_mesh.\n"
                  "Associated med_mesh name: \"%s\"\n"
                  "Associated writer name: \"%s\"\n"),
                med_mesh_name, writer->name);

    /* Add family 0 */

    {
      int i;
      strncpy(family_name, _("FAMILY_0"), MED_TAILLE_NOM);
      family_name[MED_TAILLE_NOM] = '\0';
      for (i = strlen(family_name) + 1; i < MED_TAILLE_NOM; i++)
        family_name[i] = '\0';
    }

    retval = MEDfamCr(writer->fid,
                      med_mesh_name,
                      family_name,
                      0,
                      NULL,
                      NULL,
                      NULL,
                      0,
                      NULL,
                      0);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDfamCr() failed to create family 0"
                  "for a new med_mesh.\n"
                  "Associated med_mesh name: \"%s\"\n"
                  "Associated writer name: \"%s\"\n"),
                med_mesh_name, writer->name);

  }

  writer->med_meshes[id]->num = writer->n_med_meshes;

  return writer->n_med_meshes;
}

/*----------------------------------------------------------------------------
 * Count number of extra vertices when tesselations are present
 *
 * parameters:
 *   this_writer         <-- pointer to associated writer
 *   mesh                <-- pointer to nodal mesh structure
 *   n_extra_vertices_g  --> global number of extra vertices (optional)
 *   n_extra_vertices    --> local number of extra vertices (optional)
 *----------------------------------------------------------------------------*/

static void
_count_extra_vertices(const fvmc_to_med_writer_t  *this_writer,
                      const fvmc_nodal_t          *mesh,
                      fvmc_gnum_t                 *n_extra_vertices_g,
                      fvmc_lnum_t                 *n_extra_vertices)
{
  int  i;

  const int  export_dim = fvmc_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  if (n_extra_vertices_g != NULL)
    *n_extra_vertices_g = 0;
  if (n_extra_vertices != NULL)
    *n_extra_vertices   = 0;

  for (i = 0 ; i < mesh->n_sections ; i++) {

    const fvmc_nodal_section_t  *section = mesh->sections[i];

    /* Output if entity dimension equal to highest in mesh
       (i.e. no output of faces if cells present, or edges
       if cells or faces) */

    if (   section->entity_dim == export_dim
        && section->type == FVMC_CELL_POLY
        && section->tesselation != NULL
        && this_writer->divide_polyhedra == true) {

      if (n_extra_vertices_g != NULL)
        *n_extra_vertices_g
          += fvmc_tesselation_n_g_vertices_add(section->tesselation);

      if (n_extra_vertices != NULL)
        *n_extra_vertices
          += fvmc_tesselation_n_vertices_add(section->tesselation);

    }

  }

}

/*----------------------------------------------------------------------------
 * Return extra vertex coordinates when tesselations are present
 *
 * parameters:
 *   this_writer <-- pointer to associated writer
 *   mesh        <-- pointer to nodal mesh structure that should be written
 *
 * returns:
 *   array containing all extra vertex coordinates
 *----------------------------------------------------------------------------*/

static fvmc_coord_t *
_extra_vertex_coords(const fvmc_to_med_writer_t  *this_writer,
                     const fvmc_nodal_t          *mesh)
{
  int  i;
  fvmc_lnum_t  n_extra_vertices_section;

  fvmc_lnum_t  n_extra_vertices = 0;
  size_t  coord_shift = 0;
  fvmc_coord_t  *coords = NULL;

  _count_extra_vertices(this_writer,
                        mesh,
                        NULL,
                        &n_extra_vertices);

  if (n_extra_vertices > 0) {

    const int  export_dim = fvmc_nodal_get_max_entity_dim(mesh);

    BFTC_MALLOC(coords, n_extra_vertices * 3, fvmc_coord_t);

    for (i = 0 ; i < mesh->n_sections ; i++) {

      const fvmc_nodal_section_t  *section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (   section->entity_dim == export_dim
          && section->type == FVMC_CELL_POLY
          && section->tesselation != NULL
          && this_writer->divide_polyhedra == true) {

        n_extra_vertices_section
          = fvmc_tesselation_n_vertices_add(section->tesselation);

        if (n_extra_vertices_section > 0) {

          fvmc_tesselation_vertex_coords(section->tesselation,
                                        coords + coord_shift);

          coord_shift += n_extra_vertices_section * 3;

        }

      }
    }
  }

  return coords;
}

/*----------------------------------------------------------------------------
 * Define MED geometrical element type according to FVM element type
 *
 * parameters:
 *   fvmc_elt_type  <-- pointer to fvm element type.
 *
 * return:
 *   med geometrical element type.
 *----------------------------------------------------------------------------*/

static med_geometrie_element
_get_med_elt_type(const fvmc_element_t fvmc_elt_type)
{
  med_geometrie_element  med_elt_type;

    switch (fvmc_elt_type) {

    case FVMC_EDGE:
      med_elt_type = MED_SEG2;
      break;

    case FVMC_FACE_TRIA:
      med_elt_type = MED_TRIA3;
      break;

    case FVMC_FACE_QUAD:
      med_elt_type = MED_QUAD4;
      break;

    case FVMC_FACE_POLY:
      med_elt_type = MED_POLYGONE;
      break;

    case FVMC_CELL_TETRA:
      med_elt_type = MED_TETRA4;
      break;

    case FVMC_CELL_PYRAM:
      med_elt_type = MED_PYRA5;
      break;

    case FVMC_CELL_PRISM:
      med_elt_type = MED_PENTA6;
      break;

    case FVMC_CELL_HEXA:
      med_elt_type = MED_HEXA8;
      break;

    case FVMC_CELL_POLY:
      med_elt_type = MED_POLYEDRE;
      break;

    default:
      med_elt_type = MED_NONE;
      bftc_error(__FILE__, __LINE__, 0,
                "_get_med_elt_type(): "
                "No association with MED element type has been found\n"
                "FVM element type: \"%i\"\n",
                fvmc_elt_type);

    } /* End of switch on element type */

    return med_elt_type;
}

/*----------------------------------------------------------------------------
 * Get vertex order to describe MED element type.
 *
 * parameters:
 *   med_elt_type  <-- MED element type.
 *   vertex_order  --> Pointer to vertex order array.
 *
 *----------------------------------------------------------------------------*/

static void
_get_vertex_order(const med_geometrie_element med_elt_type,
                  int  *vertex_order)
{
  switch(med_elt_type) {

  case MED_SEG2:
    vertex_order[0] = 1;
    vertex_order[1] = 2;
    break;

  case MED_TRIA3:
    vertex_order[0] = 1;
    vertex_order[1] = 2;
    vertex_order[2] = 3;
    break;

  case MED_QUAD4:
    vertex_order[0] = 1;
    vertex_order[1] = 2;
    vertex_order[2] = 3;
    vertex_order[3] = 4;
    break;

  case MED_TETRA4:
    vertex_order[0] = 1;
    vertex_order[1] = 3;
    vertex_order[2] = 2;
    vertex_order[3] = 4;
    break;

  case MED_PYRA5:
    vertex_order[0] = 1;
    vertex_order[1] = 4;
    vertex_order[2] = 3;
    vertex_order[3] = 2;
    vertex_order[4] = 5;
    break;

  case MED_PENTA6:
    vertex_order[0] = 1;
    vertex_order[1] = 3;
    vertex_order[2] = 2;
    vertex_order[3] = 4;
    vertex_order[4] = 6;
    vertex_order[5] = 5;
    break;

  case MED_HEXA8:
    vertex_order[0] = 1;
    vertex_order[1] = 4;
    vertex_order[2] = 3;
    vertex_order[3] = 2;
    vertex_order[4] = 5;
    vertex_order[5] = 8;
    vertex_order[6] = 7;
    vertex_order[7] = 6;
    break;

  case MED_POLYGONE:
    vertex_order[0] = 0;
    break;

  case MED_POLYEDRE:
    vertex_order[0] = 0;
    break;

  default:
    bftc_error(__FILE__, __LINE__, 0,
              "_get_vertex_order(): No associated MED element type known\n"
              "MED element type: \"%i\"\n",
              med_elt_type);
  }

  return;
}

/*----------------------------------------------------------------------------
 * Compute the connectivity size of the current section.
 *
 * parameters:
 *   writer           <-- pointer to MED writer structure.
 *   export_sections  <-> pointer to a list of MED section structures
 *
 * returns:
 * connect_section_size <-> connectivity size of the current section.
 *----------------------------------------------------------------------------*/

static size_t
_get_connect_section_size(const fvmc_to_med_writer_t   *writer,
                          const fvmc_writer_section_t  *export_sections)
{
  fvmc_gnum_t n_g_sub_elements = 0;
  fvmc_gnum_t n_g_elements = 0;

  size_t connect_section_size = 0;

  const fvmc_nodal_section_t *section = export_sections->section;

  if (export_sections->type == section->type) {

    /* Ordinary section */
    /*------------------*/

    if (section->stride > 0) {

      n_g_elements = fvmc_nodal_section_n_g_elements(section);
      connect_section_size = n_g_elements * section->stride;

    }
    else { /* section->stride == 0 */

      fvmc_lnum_t  i, j, f_id;

      size_t l_connect_size = 0;

      const fvmc_lnum_t  *f_idx = section->face_index;
      const fvmc_lnum_t  *f_num = section->face_num;
      const fvmc_lnum_t  *v_idx = section->vertex_index;

      assert(   section->type == FVMC_FACE_POLY
             || section->type == FVMC_CELL_POLY);

      /* Compute local connectivity size */

      if (section->type == FVMC_CELL_POLY) {

        for (i = 0; i < section->n_elements; i++) {
          for (j = f_idx[i]; j < f_idx[i+1]; j++) {
            f_id = FVMC_ABS(f_num[j]) - 1;
            l_connect_size += v_idx[f_id + 1] - v_idx[f_id];
          }
        }

      }
      else /* if (section->type == FVMC_FACE_POLY) */
        l_connect_size = section->connectivity_size;

      if (writer->n_ranks > 1) {
#if defined(FVMC_HAVE_MPI)
        MPI_Allreduce(&l_connect_size,
                      &connect_section_size,
                      1,
                      MPI_UNSIGNED_LONG,
                      MPI_SUM,
                      writer->comm);
#endif
      }
      else  /* Serial mode */
        connect_section_size = l_connect_size;

    }

  }
  else {

    /* Tesselated section */
    /*--------------------*/

    fvmc_tesselation_get_global_size(export_sections->section->tesselation,
                                    export_sections->type,
                                    &n_g_sub_elements,
                                    NULL);

    connect_section_size =
      n_g_sub_elements * fvmc_nodal_n_vertices_element[export_sections->type];

  }

  return connect_section_size;
}

/*----------------------------------------------------------------------------
 * Compute connectivity buffer size for MED writer file.
 *
 * parameters:
 *   writer           <-- pointer to MED writer structure.
 *   export_sections  <-> pointer to a list of MED section structures
 *
 * returns:
 *   global_connect_buffer_size: maximum buffer size useful to export
 *   connectivity.
 *----------------------------------------------------------------------------*/

static size_t
_get_connect_buffer_size(const fvmc_to_med_writer_t   *writer,
                         const fvmc_writer_section_t  *export_sections)

{
  med_geometrie_element  med_type;

  size_t  connect_section_size = 0;
  size_t  global_connect_buffer_size = 0;

  const fvmc_writer_section_t *current_section = NULL;

  med_geometrie_element  previous_med_type
    = _get_med_elt_type(export_sections->type);

  current_section = export_sections;

  while (current_section != NULL) {

    med_type = _get_med_elt_type(current_section->type);

    if (med_type != previous_med_type) {

      previous_med_type = med_type;

      connect_section_size = _get_connect_section_size(writer,
                                                       current_section);

    }
    else
      connect_section_size += _get_connect_section_size(writer,
                                                        current_section);

    global_connect_buffer_size = FVMC_MAX(connect_section_size,
                                         global_connect_buffer_size);

    current_section = current_section->next;

  } /* End of loop on sections */

  /* When mesh has only one section */

  if (  global_connect_buffer_size == 0
     && connect_section_size > 0)
    global_connect_buffer_size = connect_section_size;

  return global_connect_buffer_size;
}

/*----------------------------------------------------------------------------
 * Adapt datatype to be MED compliant and associate its own MED datatype
 *
 * parameters:
 *   input_fvmc_datatype   <-- input field datatype in FVM.
 *   output_fvmc_datatype  --> output field datatype in FVM.
 *   med_datatype         --> associated MED field datatype.
 *   data_sizeof          --> size of datatype to be exported.
 *----------------------------------------------------------------------------*/

static void
_get_datatypes(const fvmc_datatype_t    input_fvmc_datatype,
               fvmc_datatype_t         *output_fvmc_datatype,
               med_type_champ         *med_datatype,
               int                    *data_sizeof)
{
  int med_int_sizeof, med_float_sizeof;

  assert(sizeof(med_float) == 8);

  /* Verify if double define over 8 bytes */

#if (FVMC_SIZEOF_DOUBLE != 8)
#error
#endif

  med_int_sizeof = sizeof(med_int);
  med_float_sizeof = sizeof(med_float);

  /* Define output_fvmc_datatype and med_datatype by input_fvmc_datatype */

  switch(input_fvmc_datatype) {
  case FVMC_DOUBLE:
    *output_fvmc_datatype = FVMC_DOUBLE;
    *med_datatype = MED_FLOAT64;
    *data_sizeof = med_float_sizeof;
    break;

  case FVMC_FLOAT:
    *output_fvmc_datatype = FVMC_DOUBLE;
    *med_datatype = MED_FLOAT64;
    *data_sizeof = med_float_sizeof;
    break;

  case FVMC_INT32:
    if (med_int_sizeof == 4) {
      *output_fvmc_datatype = FVMC_INT32;
      *med_datatype = MED_INT32;
      *data_sizeof = med_int_sizeof;
    }
    else if (med_int_sizeof == 8) {
      *output_fvmc_datatype = FVMC_INT64;
      *med_datatype = MED_INT64;
      *data_sizeof = med_int_sizeof;
    }
    else
      assert(0);
    break;

  case FVMC_INT64:
    if (med_int_sizeof == 4) {
      *output_fvmc_datatype = FVMC_INT32;
      *med_datatype = MED_INT32;
      *data_sizeof = med_int_sizeof;
    }
    else if (med_int_sizeof == 8) {
      *output_fvmc_datatype = FVMC_INT64;
      *med_datatype = MED_INT64;
      *data_sizeof = med_int_sizeof;
    }
    else
      assert(0);
    break;

  case FVMC_UINT32:
    if (med_int_sizeof == 4) {
      *output_fvmc_datatype = FVMC_INT32;
      *med_datatype = MED_INT32;
      *data_sizeof = med_int_sizeof;
    }
    else if (med_int_sizeof == 8) {
      *output_fvmc_datatype = FVMC_INT64;
      *med_datatype = MED_INT64;
      *data_sizeof = med_int_sizeof;
    }
    else
      assert(0);
    break;

  case FVMC_UINT64:
    if (med_int_sizeof == 4) {
      *output_fvmc_datatype = FVMC_INT32;
      *med_datatype = MED_INT32;
      *data_sizeof = med_int_sizeof;
    }
    else if (med_int_sizeof == 8) {
      *output_fvmc_datatype = FVMC_INT64;
      *med_datatype = MED_INT64;
      *data_sizeof = med_int_sizeof;
    }
    else
      assert(0);
    break;

  default:
    assert(0);
  }

  return;
}

/*----------------------------------------------------------------------------
 * Get med fieldname, update field struture and call MEDchampCr() if necessary
 *
 * parameters:
 *   writer          <-- MED writer structure.
 *   fieldname       <-- input fieldname.
 *   datatype_med    <-- associated MED field datatype.
 *   dimension       <-- dimension of field to export.
 *   med_fieldname   --> MED name of the field.
 *
 *----------------------------------------------------------------------------*/

static void
_get_med_fieldname(fvmc_to_med_writer_t    *writer,
                   const char             *fieldname,
                   med_type_champ          datatype_med,
                   int                     dimension,
                   char                   *med_fieldname)
{
  int   i, i_char, i_dim, n_chars;
  int   n_fields, i_field, name_size;
  med_int        n_components;

  char *component_name = NULL;
  char *units_name = NULL;

  const int     rank = writer->rank;

  med_err retval = 0;

  /* Fieldname adaptation */

  strncpy(med_fieldname, fieldname, MED_TAILLE_NOM);
  n_chars = strlen(med_fieldname);

  for (i_char = n_chars + 1; i_char < MED_TAILLE_NOM; i_char++)
    med_fieldname[i_char] = ' ';

  med_fieldname[MED_TAILLE_NOM] = '\0';

  /* Loop on fields to know if field has already been created */

  n_fields = writer->n_fields;

  for (i_field = 0; i_field < n_fields; i_field++) {

    fvmc_to_med_field_t *field = writer->fields[i_field];

    if (strcmp(med_fieldname, field->name) == 0)
      break;

  }

  if (i_field == n_fields) { /* Create a new field for this writer */

    writer->n_fields++;
    BFTC_REALLOC(writer->fields, writer->n_fields, fvmc_to_med_field_t *);

    BFTC_MALLOC(writer->fields[n_fields], 1, fvmc_to_med_field_t);

    for (i = 0; i < MED_TAILLE_NOM; i++)
      writer->fields[n_fields]->name[i] = med_fieldname[i];
    writer->fields[n_fields]->name[MED_TAILLE_NOM] = '\0';
    writer->fields[n_fields]->id = n_fields;
    writer->fields[n_fields]->n_components = dimension;
    writer->fields[n_fields]->datatype = datatype_med;

    if (rank == 0) {

      /* Component and unit names */

      name_size = MED_TAILLE_PNOM * dimension;
      BFTC_MALLOC(component_name, name_size + 1, char);
      BFTC_MALLOC(units_name, name_size + 1, char);

      component_name[0] = '\0';
      units_name[0] = '\0';

      for (i = 1; i < name_size; i++) {
        component_name[i] = ' ';
        units_name[i] = ' ';
      }

      if (dimension == 1)
        sprintf(&component_name[0], "Scalar");

      else if (dimension == 3) {
        const char  *const xyz[3] = {"X", "Y", "Z"};
        for (i_dim = 0; i_dim < dimension; i_dim++)
          sprintf(&component_name[i_dim * MED_TAILLE_PNOM],
                  "Component %s", xyz[i_dim]);
      }

      else if (dimension == 6) {
        const char  *const ij_sym[6] = {"11", "22", "33", "12", "13", "23"};
        for (i_dim = 0; i_dim < dimension; i_dim++)
          sprintf(&component_name[i_dim * MED_TAILLE_PNOM],
                  "Component %s", ij_sym[i_dim]);
      }

      else if (dimension == 9) {
        const char  *const ij_asym[9] = {"11", "12", "13",
                                         "21", "22", "23",
                                         "31", "32", "33" };
        for (i_dim = 0; i_dim < dimension; i_dim++)
          sprintf(&component_name[i_dim * MED_TAILLE_PNOM],
                  "Component %s", ij_asym[i_dim]);
      }

      component_name[name_size] = '\0';
      units_name[name_size] = '\0';

      /* Creation of a MED field */

      assert(writer->is_open == true);
      n_components = (med_int)dimension;

      retval = MEDchampCr(writer->fid,
                          med_fieldname,
                          datatype_med,
                          component_name,
                          units_name,
                          n_components);

      if (retval < 0)
        bftc_error(__FILE__, __LINE__, 0,
                  "MEDchampCr() failed to create a field.\n"
                  "Associated writer name: \"%s\"\n"
                  "Associated fieldname: \"%s\"\n",
                  writer->name, med_fieldname);

      BFTC_FREE(units_name);
      BFTC_FREE(component_name);

    } /* End if rank = 0 */

  } /* End of field creation */

  /* If field exists, check that dimensions and type are compatible */

  else { /*  if (i_field < n_fields) */

    fvmc_to_med_field_t *field = writer->fields[i_field];

    if (dimension != field->n_components)
      bftc_error(__FILE__, __LINE__, 0,
                _("MED field \"%s\" already defined\n"
                  "for writer \"%s\" with %d components,\n"
                  "but re-defined with %d components."),
                med_fieldname, writer->name,
                (int)field->n_components, (int)dimension);

    if (datatype_med != field->datatype )
      bftc_error(__FILE__, __LINE__, 0,
                _("MED field \"%s\" already defined\n"
                  "for writer \"%s\" with datatype %d,\n"
                  "but re-defined with datatype %d."),
                med_fieldname, writer->name,
                (int)field->datatype, datatype_med);

  }

  return;
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write vertex coordinates to a MED file in parallel mode
 *
 * parameters:
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *   med_mesh      <-- pointer to med_mesh structure associated with the mesh.
 *   writer        <-- pointer to associated writer.
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_g(const fvmc_nodal_t    *mesh,
                        fvmc_to_med_mesh_t    *med_mesh,
                        fvmc_to_med_writer_t  *writer)
{
  int  i, i_dim, name_size;
  fvmc_lnum_t   i_lnod;
  MPI_Datatype  mpi_datatype;
  fvmc_datatype_t  fvmc_datatype;

  fvmc_lnum_t  n_extra_vertices = 0;
  fvmc_gnum_t  n_g_extra_vertices = 0;
  fvmc_gnum_t  global_num_start = 0, global_num_end = 0;

  char  *unit_name = NULL;
  char  *coord_name = NULL;
  fvmc_coord_t  *extra_vertex_coords = NULL;
  med_float  *med_coords = NULL;
  med_float  *global_coords_buffer = NULL;
  fvmc_gather_slice_t  *vertices_slice = NULL;
  fvmc_gather_slice_t  *extra_vertices_slice = NULL;

  const int          dim = mesh->dim;
  const int          rank = writer->rank;
  const MPI_Comm     comm = writer->comm;
  const fvmc_lnum_t   n_vertices = mesh->n_vertices;
  const fvmc_gnum_t   n_g_vertices =
    fvmc_io_num_get_global_count(mesh->global_vertex_num);

  const char          *const XYZ[3] = {"X", "Y", "Z"};
  const fvmc_lnum_t    *parent_vertex_num = mesh->parent_vertex_num;
  const fvmc_coord_t   *vertex_coords = mesh->vertex_coords;
  const fvmc_io_num_t  *global_vertex_num = mesh->global_vertex_num;

  med_err   retval = 0;

  assert(writer->is_open == true);
  assert(med_mesh != NULL);

  /* MED does not accept an empty mesh */

  if (n_g_vertices == 0)
    bftc_error(__FILE__, __LINE__, 0,
              "MED does not allow to export an empty mesh,\n"
              "Mesh: \"%s\" has no vertex.\n"
              "Associated file: \"%s\".",
              mesh->name, writer->filename);

  /* Define MPI and FVM datatype */

#if (FVMC_COORD == FVMC_DOUBLE)
  fvmc_datatype = FVMC_DOUBLE;
  mpi_datatype = MPI_DOUBLE;
#else
  if (sizeof(med_float) == sizeof(fvmc_coord)) {
    mpi_datatype = MPI_FLOAT;
    fvmc_datatype = FVMC_FLOAT;
  }
  else {
    bftc_error(__FILE__, __LINE__, 0 ,
              _("Incompatible datatype sizes between MPI, FVM and MED."));
  }
#endif

  if (rank == 0) {

    /* Coordinates name and unit names */
    /*---------------------------------*/

    name_size = MED_TAILLE_PNOM * dim;
    BFTC_MALLOC(coord_name, name_size + 1, char);
    BFTC_MALLOC(unit_name, name_size + 1, char);

    coord_name[0] = '\0';
    unit_name[0] = '\0';

    for (i = 1; i < name_size; i++) {
      coord_name[i] = ' ';
      unit_name[i] = ' ';
    }

    for (i_dim = 0; i_dim < dim; i_dim++)
      sprintf(&coord_name[i_dim * MED_TAILLE_PNOM],"Coord%s",XYZ[i_dim]);

    coord_name[name_size] = '\0';
    unit_name[name_size] = '\0';

  } /* End if rank == 0 */

  /* Compute extra vertices buffer size if necessary */

  _count_extra_vertices(writer,
                        mesh,
                        &n_g_extra_vertices,
                        &n_extra_vertices);

  /* Get local extra vertex coords */

  extra_vertex_coords = _extra_vertex_coords(writer,
                                             mesh);

  /* Export vertex coordinates to a MED file */
  /*-----------------------------------------*/

  vertices_slice = fvmc_gather_slice_create(global_vertex_num,
                                           n_g_vertices,
                                           comm);

  BFTC_MALLOC(med_coords,
             FVMC_MAX(n_vertices, n_extra_vertices) * dim,
             med_float);

  if (parent_vertex_num != NULL) {

    fvmc_lnum_t  idx = 0;

    for (i_lnod = 0; i_lnod < n_vertices; i_lnod++) {
      for (i_dim = 0; i_dim < dim; i_dim++)
        med_coords[idx++]
          = (med_float)vertex_coords[(parent_vertex_num[i_lnod]-1)*dim + i_dim];
    }

  }
  else
    _convert_float_fvmc_to_med(vertex_coords,
                              fvmc_datatype,
                              med_coords,
                              n_vertices * dim);

  BFTC_MALLOC(global_coords_buffer,
             (n_g_vertices + n_g_extra_vertices) * dim,
             med_float);

  /* Gather slices from other ranks to rank 0 */

  while (fvmc_gather_slice_advance(vertices_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    fvmc_gather_array(med_coords,
                     global_coords_buffer,
                     mpi_datatype,
                     dim,
                     global_vertex_num,
                     comm,
                     vertices_slice);

  } /* End of slice advance */

  fvmc_gather_slice_destroy(vertices_slice);

  if (n_g_extra_vertices > 0) {

    /* Handle extra vertices */
    /*-----------------------*/

    int section_id;
    fvmc_lnum_t extra_vertices_count = 0;
    fvmc_gnum_t extra_vertices_count_g = 0;

    for (section_id = 0 ; section_id < mesh->n_sections ; section_id++) {

      const fvmc_nodal_section_t  *section = mesh->sections[section_id];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (   section->entity_dim == mesh->dim
          && section->type == FVMC_CELL_POLY
          && section->tesselation != NULL
          && writer->divide_polyhedra == true) {

        const fvmc_io_num_t *extra_vertex_num
          = fvmc_tesselation_global_vertex_num(section->tesselation);
        const fvmc_gnum_t n_extra_vertices_section
          = fvmc_tesselation_n_vertices_add(section->tesselation);
        const fvmc_gnum_t n_g_extra_vertices_section
          = fvmc_tesselation_n_g_vertices_add(section->tesselation);

        fvmc_lnum_t  idx = 0;

        for (i_lnod = 0 ; i_lnod < n_extra_vertices ; i_lnod++) {
          for (i_dim = 0; i_dim < dim; i_dim++)
            med_coords[idx++] = (med_float)
              extra_vertex_coords[(i_lnod + extra_vertices_count)*dim + i_dim];
        }

        extra_vertices_slice = fvmc_gather_slice_create(extra_vertex_num,
                                                       n_g_extra_vertices_section,
                                                       comm);

        /* loop on slices in parallel mode */

        while (fvmc_gather_slice_advance(extra_vertices_slice,
                                        &global_num_start,
                                        &global_num_end) == 0) {

          fvmc_gather_array(med_coords,
                           global_coords_buffer
                           + (n_g_vertices + extra_vertices_count_g) * dim ,
                           mpi_datatype,
                           dim,
                           extra_vertex_num,
                           comm,
                           extra_vertices_slice);

        }

        fvmc_gather_slice_destroy(extra_vertices_slice);

        extra_vertices_count_g += n_g_extra_vertices_section;
        extra_vertices_count   += n_extra_vertices_section;

      }

    } /* end of loop on sections for extra vertices */

  } /* n_g_extra_vertices > 0 */

  if (rank == 0) {

    /* Write all the coordinates */
    /*---------------------------*/

    if (global_coords_buffer != NULL)
      retval = MEDcoordEcr(writer->fid,
                           med_mesh->name,
                           /* according to MED documentation */
                           med_mesh->entity_dim,
                           global_coords_buffer,
                           MED_FULL_INTERLACE,
                           (med_int)(n_g_vertices + n_g_extra_vertices),
                           MED_CART,
                           coord_name,
                           unit_name);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDcoordEcr() failed to write coords:\"%s\"\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                coord_name, writer->name, med_mesh->name);

    BFTC_FREE(coord_name);
    BFTC_FREE(unit_name);

  } /* End if rank == 0 */

  /* Free buffers */

  BFTC_FREE(med_coords);
  BFTC_FREE(global_coords_buffer);

  if (extra_vertex_coords != NULL)
    BFTC_FREE(extra_vertex_coords);

}

#endif

/*----------------------------------------------------------------------------
 * Write vertex coordinates to a MED file in serial mode
 *
 * parameters:
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *   med_mesh      <-- pointer to med_mesh structure associated with the mesh.
 *   writer        <-- pointer to associated writer.
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_l(const fvmc_nodal_t     *mesh,
                        fvmc_to_med_mesh_t     *med_mesh,
                        fvmc_to_med_writer_t   *writer)
{
  int  i, i_dim, name_size;
  fvmc_lnum_t  i_vtx;
  fvmc_datatype_t  datatype;

  fvmc_lnum_t  idx = 0;
  fvmc_lnum_t  n_extra_vertices = 0;

  char  *unit_name = NULL;
  char  *coord_name = NULL;
  med_float  *med_coords = NULL;
  med_float  *all_med_coords = NULL;
  fvmc_coord_t  *extra_vertex_coords = NULL;

  const int  dim = mesh->dim;
  const fvmc_lnum_t  n_vertices = mesh->n_vertices;
  const char         *const XYZ[3] = {"X", "Y", "Z"};
  const fvmc_lnum_t   *parent_vertex_num = mesh->parent_vertex_num;
  const fvmc_coord_t  *vertex_coords = mesh->vertex_coords;

  med_err retval = 0;


  assert(writer->is_open == true);
  assert(med_mesh != NULL);

  /* MED does not accept an empty mesh */

  if (n_vertices == 0)
    bftc_error(__FILE__, __LINE__, 0,
              "MED does not allow to export an empty mesh,\n"
              "Mesh: \"%s\" has no vertex.\n"
              "Associated file: \"%s\".",
              mesh->name, writer->filename);

  /* Get fvm datatype */

#if (FVMC_COORD == FVMC_DOUBLE)
  datatype = FVMC_DOUBLE;
#else
  datatype = FVMC_FLOAT;
#endif

  /* Coordinates name and unit names */

  name_size = MED_TAILLE_PNOM * dim;
  BFTC_MALLOC(coord_name, name_size + 1, char);
  BFTC_MALLOC(unit_name, name_size + 1, char);

  coord_name[0] = '\0';
  unit_name[0] = '\0';

  for (i = 1; i < name_size; i++) {
    coord_name[i] = ' ';
    unit_name[i] = ' ';
  }

  for (i_dim = 0; i_dim < dim; i_dim++)
    sprintf(&coord_name[i_dim * MED_TAILLE_PNOM],"Coord%s",XYZ[i_dim]);

  coord_name[name_size] = '\0';
  unit_name[name_size] = '\0';

  /* Compute extra vertex coordinates if present */

  _count_extra_vertices(writer,
                        mesh,
                        NULL,
                        &n_extra_vertices);

  extra_vertex_coords = _extra_vertex_coords(writer,
                                             mesh);


  /* Vertex coordinates export */
  /*---------------------------*/

  BFTC_MALLOC(med_coords, n_vertices * dim, med_float);

  if (parent_vertex_num != NULL) {

    for (i_vtx = 0; i_vtx < n_vertices; i_vtx++) {
      for (i_dim = 0; i_dim < dim; i_dim++)
        med_coords[idx++] = (med_float)
          vertex_coords[(parent_vertex_num[i_vtx]-1) * dim + i_dim];
    }

  }
  else
    _convert_float_fvmc_to_med(vertex_coords,
                              datatype,
                              med_coords,
                              n_vertices * dim);


  /* Extra vertices coordinates */
  /*----------------------------*/

  if (n_extra_vertices > 0) {

    BFTC_REALLOC(med_coords,
                (n_vertices + n_extra_vertices)*dim,
                med_float);

    /* Convert fvmc_coord_t to med_float */

    _convert_float_fvmc_to_med(extra_vertex_coords,
                              datatype,
                              med_coords + n_vertices * dim,
                              n_extra_vertices * dim);

  } /* End if n_extra_vertices > 0 */

  /* Write all coordinates */

  if (med_coords != NULL)
    retval = MEDcoordEcr(writer->fid,
                         med_mesh->name,
                         /* According to MED documentation */
                         med_mesh->entity_dim,
                         med_coords,
                         MED_FULL_INTERLACE,
                         (med_int)(n_vertices + n_extra_vertices),
                         MED_CART,
                         coord_name,
                         unit_name);

  if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDcoordEcr() failed to write coords:\"%s\"\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                coord_name, writer->name, med_mesh->name);

  /* Free buffers */

  BFTC_FREE(med_coords);

  if (all_med_coords != NULL)
    BFTC_FREE(all_med_coords);

  BFTC_FREE(coord_name);
  BFTC_FREE(unit_name);

  if (extra_vertex_coords != NULL)
    BFTC_FREE(extra_vertex_coords);

  return;
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write strided elements connectivity to a MED file in parallel mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   med_family            <-- med family number associated.
 *   export_connect        <-- buffer to export connectivity.
 *
 * returns:
 *  pointer to next MED section structure in list
 *----------------------------------------------------------------------------*/

static const fvmc_writer_section_t *
_export_connect_g(const fvmc_writer_section_t  *export_sections,
                  fvmc_to_med_writer_t         *writer,
                  const fvmc_nodal_t           *mesh,
                  fvmc_to_med_mesh_t           *med_mesh,
                  med_int                     *med_family,
                  char                        *export_connect)
{
  int vertex_order[FVMC_MED_MAX_N_NODES];
  fvmc_lnum_t  i_vtx, i_elt;
  med_geometrie_element  med_section_type;

  fvmc_lnum_t  i_num = 0;
  fvmc_gnum_t  global_num_start = 0, global_num_end = 0;
  fvmc_gnum_t  n_export_elements = 0;

  fvmc_gnum_t *_fvmc_export_connect = (fvmc_gnum_t *)export_connect;
  med_int    *_med_export_connect = (med_int *)export_connect;
  fvmc_gather_slice_t *elements_slice = NULL;

  const int stride = fvmc_nodal_n_vertices_element[export_sections->type];
  const fvmc_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  /* Get MED element vertex order */

  med_section_type = _get_med_elt_type(current_section->type);

  _get_vertex_order(med_section_type,
                    vertex_order);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvmc_nodal_section_t *section = current_section->section;
    const fvmc_gnum_t  n_g_elements_section
      = fvmc_nodal_section_n_g_elements(section);

    if (section->type == current_section->type) {

      /* Ordinary section */
      /*------------------*/

      fvmc_lnum_t *_vertex_num = NULL;

      const fvmc_lnum_t *vertex_num = section->vertex_num;

      /* Convert FVM connectivity to be congruent with MED standard */

      BFTC_MALLOC(_vertex_num, stride * section->n_elements, med_int);

      i_num = 0;
      for (i_elt = 0; i_elt < section->n_elements; i_elt++) {
        for (i_vtx = 0; i_vtx < stride; i_vtx++)
          _vertex_num[i_num++] =
            vertex_num[i_elt * stride + vertex_order[i_vtx] - 1];
      }

      elements_slice = fvmc_gather_slice_create(section->global_element_num,
                                               n_g_elements_section,
                                               writer->comm);

      /* Gather slices from other ranks to rank 0 */

      while (fvmc_gather_slice_advance(elements_slice,
                                      &global_num_start,
                                      &global_num_end) == 0) {

        fvmc_gather_strided_connect(_vertex_num,
                                   _fvmc_export_connect
                                   + n_export_elements * stride,
                                   stride,
                                   mesh->global_vertex_num,
                                   section->global_element_num,
                                   writer->comm,
                                   elements_slice);

      } /* End of slice advance */

      if (writer->rank == 0)
        n_export_elements +=  n_g_elements_section;

      BFTC_FREE(_vertex_num);
      fvmc_gather_slice_destroy(elements_slice);

    }
    else {

      /* Tesselated section */
      /*--------------------*/

      size_t i_tmp;
      fvmc_lnum_t  n_elts_slice = 0;
      fvmc_lnum_t  start_id = 0, end_id = 0;
      fvmc_gnum_t  buffer_size = 0;
      fvmc_gnum_t  buffer_size_prev = 0;
      fvmc_lnum_t  n_sub_elements_max = 0;
      fvmc_gnum_t  n_g_sub_elements = 0;

      fvmc_lnum_t  *local_idx = NULL;
      fvmc_gnum_t  *global_idx = NULL;
      fvmc_gnum_t  *sub_elt_vertex_num = NULL;
      fvmc_gnum_t  tmp_connect[5];

      const fvmc_lnum_t *sub_elt_index
        = fvmc_tesselation_sub_elt_index(section->tesselation,
                                        current_section->type);

      fvmc_tesselation_get_global_size(section->tesselation,
                                      current_section->type,
                                      &n_g_sub_elements,
                                      &n_sub_elements_max);

      BFTC_MALLOC(local_idx, section->n_elements + 1, fvmc_lnum_t);
      BFTC_MALLOC(global_idx, n_g_elements_section + 1, fvmc_gnum_t);

      buffer_size = FVMC_MAX((fvmc_gnum_t)(10 * n_sub_elements_max),
                            (fvmc_gnum_t)(n_g_sub_elements / writer->n_ranks));
      buffer_size *= stride;
      buffer_size_prev = buffer_size;

      BFTC_MALLOC(sub_elt_vertex_num, buffer_size, fvmc_gnum_t);

      elements_slice = fvmc_gather_slice_create(section->global_element_num,
                                               n_g_sub_elements,
                                               writer->comm);

      while (fvmc_gather_slice_advance(elements_slice,
                                      &global_num_start,
                                      &global_num_end) == 0) {

        /* Build element->vertices index */

        end_id
          = fvmc_tesselation_range_index_g(section->tesselation,
                                          current_section->type,
                                          stride,
                                          start_id,
                                          buffer_size,
                                          &global_num_end,
                                          local_idx,
                                          writer->comm);

        /* Check if the maximum id returned on some ranks leads to a
           lower global_num_end than initially required (due to the
           local buffer being too small) and adjust slice if necessary */

        fvmc_gather_slice_limit(elements_slice, &global_num_end);

        /* Gather element->vertices index */

        fvmc_gather_slice_index(local_idx,
                               global_idx,
                               section->global_element_num,
                               writer->comm,
                               elements_slice);

        /* Recompute maximum value of global_num_end for this slice */

        fvmc_gather_resize_indexed_slice(10,
                                        &global_num_end,
                                        &buffer_size,
                                        writer->comm,
                                        global_idx,
                                        elements_slice);

        /* If the buffer already allocated is too small, reallocate it */

        if (buffer_size_prev < buffer_size) {
          BFTC_REALLOC(sub_elt_vertex_num, buffer_size, fvmc_gnum_t);
          buffer_size_prev = buffer_size;
        }

        /* Now decode tesselation */

        end_id = fvmc_tesselation_decode_g(section->tesselation,
                                          current_section->type,
                                          start_id,
                                          buffer_size,
                                          &global_num_end,
                                          mesh->global_vertex_num,
                                          current_section->extra_vertex_base,
                                          sub_elt_vertex_num,
                                          writer->comm);

        /* Convert FVM connectivity to be congruent with MED standard */

        i_num = 0;
        n_elts_slice = sub_elt_index[end_id] - sub_elt_index[start_id];
        for (i_elt = 0; i_elt < n_elts_slice; i_elt++) {
          for (i_tmp = 0, i_vtx = 0; i_vtx < stride; i_vtx++)
            tmp_connect[i_tmp++]
              = sub_elt_vertex_num[  (i_elt * stride)
                                   + (vertex_order[i_vtx] - 1)];
          for (i_tmp = 0, i_vtx = 0; i_vtx < stride; i_vtx++)
            sub_elt_vertex_num[i_num++] = tmp_connect[i_tmp++];
        }

        /* No need to check if the maximum id returned on some ranks
           leads to a lower global_num_end than initially required
           (due to local buffer being full), as this was already done
           above for the local index */

        /* Now gather decoded element->vertices connectivity */

        fvmc_gather_indexed(sub_elt_vertex_num,
                           _fvmc_export_connect
                           + n_export_elements * stride,
                           FVMC_MPI_GNUM,
                           local_idx,
                           section->global_element_num,
                           writer->comm,
                           global_idx,
                           elements_slice);

        if (writer->rank == 0)
          n_export_elements
            += ((global_idx[global_num_end - global_num_start]) / stride);

        start_id = end_id;

      } /* End of slice advance */

      BFTC_FREE(local_idx);
      BFTC_FREE(global_idx);
      BFTC_FREE(sub_elt_vertex_num);

      fvmc_gather_slice_destroy(elements_slice);

    } /* End of tesselated section */

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  /* Write buffers into MED file */
  /*-----------------------------*/

  if (writer->rank == 0) {

    _convert_fvmc_gnum_to_med_int(_fvmc_export_connect,
                                 _med_export_connect,
                                 n_export_elements * stride);

    /* Write connectivity */

    retval  = MEDconnEcr(writer->fid,
                         med_mesh->name,
                         med_mesh->entity_dim,
                         _med_export_connect,
                         MED_FULL_INTERLACE,
                         (med_int)n_export_elements,
                         MED_MAILLE,
                         med_section_type,
                         MED_NOD);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDconnEcr() failed to write strided connectivity:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh_name: \"%s\"\n"
                  "Associated MED geometrical element: \"%i\"\n"),
                writer->name, med_mesh->name, med_section_type);

    /* Write family number */

    retval  = MEDfamEcr(writer->fid,
                        med_mesh->name,
                        med_family,
                        (med_int)n_export_elements,
                        MED_MAILLE,
                        med_section_type);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDfamEcr() failed to write family numbers:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh_name: \"%s\"\n"
                  "Associated MED geometrical element: \"%i\"\n"),
                writer->name, med_mesh->name, med_section_type);

  } /* If rank == 0 */

  return current_section;
}

#endif /* FVMC_HAVE_MPI */

/*----------------------------------------------------------------------------
 * Write strided elements connectivity to a MED file in serial mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   med_family            <-- med family number associated.
 *   export_connect        <-- buffer to export connectivity.
 *
 * returns:
 *  pointer to next MED section structure in list
 *----------------------------------------------------------------------------*/

static const fvmc_writer_section_t *
_export_connect_l(const fvmc_writer_section_t  *export_sections,
                  fvmc_to_med_writer_t         *writer,
                  fvmc_to_med_mesh_t           *med_mesh,
                  med_int                     *med_family,
                  med_int                     *med_export_connect)
{
  int vertex_order[FVMC_MED_MAX_N_NODES];
  fvmc_lnum_t  i_vtx, i_elt, i_num;
  med_geometrie_element  med_section_type;

  med_int  n_export_elements = 0;

  const int stride = fvmc_nodal_n_vertices_element[export_sections->type];
  const fvmc_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  med_section_type = _get_med_elt_type(current_section->type);

  /* Get MED element type vertex order */

  _get_vertex_order(med_section_type,
                    vertex_order);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvmc_nodal_section_t *section = current_section->section;

    if (section->type == current_section->type) {

      /* Ordinary section */
      /*------------------*/

      const fvmc_lnum_t *vertex_num = section->vertex_num;

      /* Convert FVM connectivity to be congruent with MED standard */

      i_num = n_export_elements * stride;
      for (i_elt = 0; i_elt < section->n_elements; i_elt++) {
        for (i_vtx = 0; i_vtx < stride; i_vtx++)
          med_export_connect[i_num++] =
            (med_int)vertex_num[i_elt * stride + vertex_order[i_vtx] - 1];
      }

      n_export_elements +=  section->n_elements;

    }
    else {

      /* Tesselated section */
      /*--------------------*/

      size_t  buffer_size = 0;
      fvmc_lnum_t  start_id = 0, end_id = 0;
      fvmc_lnum_t  n_sub_elts = 0, n_sub_loc = 0, n_sub_elts_max = 0;

      fvmc_lnum_t  *sub_elt_vertex_num = NULL;

      const fvmc_lnum_t *sub_elt_index
        = fvmc_tesselation_sub_elt_index(section->tesselation,
                                        current_section->type);

      n_sub_elts = fvmc_tesselation_n_sub_elements(section->tesselation,
                                                  current_section->type);

      fvmc_tesselation_get_global_size(section->tesselation,
                                      current_section->type,
                                      NULL,
                                      &n_sub_elts_max);

      buffer_size = FVMC_MAX(n_sub_elts_max, n_sub_elts/4 + 1);
      buffer_size = FVMC_MAX(256, buffer_size);

      BFTC_MALLOC(sub_elt_vertex_num, buffer_size * stride, fvmc_lnum_t);

      do {

        end_id
          = fvmc_tesselation_decode(section->tesselation,
                                   current_section->type,
                                   start_id,
                                   buffer_size,
                                   current_section->extra_vertex_base,
                                   sub_elt_vertex_num);

        /* Convert FVM connectivity to be congruent with MED standard */

        n_sub_loc = sub_elt_index[end_id] - sub_elt_index[start_id];

        i_num = n_export_elements * stride;
        for (i_elt = 0; i_elt < n_sub_loc; i_elt++) {
          for (i_vtx = 0; i_vtx < stride; i_vtx++)
            med_export_connect[i_num++]
              = (med_int)sub_elt_vertex_num[  (i_elt * stride)
                                            + (vertex_order[i_vtx] - 1)];
        }

        n_export_elements += n_sub_loc;
        start_id = end_id;

      } while (end_id < section->n_elements);

      BFTC_FREE(sub_elt_vertex_num);

    } /* End of tesselated section */

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  /* Write buffers into MED file */
  /*-----------------------------*/

  retval  = MEDconnEcr(writer->fid,
                       med_mesh->name,
                       med_mesh->entity_dim,
                       med_export_connect,
                       MED_FULL_INTERLACE,
                       n_export_elements,
                       MED_MAILLE,
                       med_section_type,
                       MED_NOD);

  if (retval < 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("MEDconnEcr() failed to write strided connectivity:\n"
                "Associated writer: \"%s\"\n"
                "Associated med_mesh_name: \"%s\"\n"
                "Associated MED geometrical element: \"%i\"\n"),
              writer->name, med_mesh->name, med_section_type);

  /* Write family numbers */

  retval  = MEDfamEcr(writer->fid,
                      med_mesh->name,
                      med_family,
                      n_export_elements,
                      MED_MAILLE,
                      med_section_type);

  if (retval < 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("MEDfamEcr() failed to write family numbers:\n"
                "Associated writer: \"%s\"\n"
                "Associated med_mesh_name: \"%s\"\n"
                "Associated MED geometrical element: \"%i\"\n"),
              writer->name, med_mesh->name, med_section_type);

  return current_section;
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polygonal elements connectivity to a MED file in parallel mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   med_family            <-- med family number associated.
 *   export_connect        <-- buffer to export connectivity.
 *----------------------------------------------------------------------------*/

static const fvmc_writer_section_t *
_export_nodal_polygons_g(const fvmc_writer_section_t  *export_sections,
                         fvmc_to_med_writer_t         *writer,
                         const fvmc_nodal_t           *mesh,
                         fvmc_to_med_mesh_t           *med_mesh,
                         med_int                     *med_family,
                         char                        *export_connect)
{
  fvmc_gnum_t  i;
  med_geometrie_element  med_section_type;

  int   n_passes = 0;
  fvmc_gnum_t  _n_connect_size = 0, n_g_connect_size = 0;
  fvmc_gnum_t  global_num_start = 0, global_num_end = 0;
  fvmc_gnum_t  n_export_connect = 0;
  fvmc_gnum_t  n_export_elements = 0;

  char  *global_vtx_idx_buffer = NULL;
  fvmc_gnum_t *fvmc_global_vtx_idx = NULL;
  med_int    *med_global_vtx_idx = NULL;

  fvmc_gather_slice_t *polygons_slice = NULL;

  fvmc_gnum_t *_fvmc_export_connect = (fvmc_gnum_t *)export_connect;
  med_int    *_med_export_connect = (med_int *)export_connect;

  const size_t export_datasize = FVMC_MAX(sizeof(fvmc_gnum_t),sizeof(med_int));
  const fvmc_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  /* Get MED element type */

  med_section_type = _get_med_elt_type(current_section->type);

  assert(med_section_type == MED_POLYGONE);
  assert(writer->discard_polygons == false);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvmc_nodal_section_t *section = current_section->section;
    const fvmc_gnum_t  n_g_element_section
      = fvmc_nodal_section_n_g_elements(section);

    n_passes++;

    _n_connect_size = (fvmc_gnum_t)section->connectivity_size;
    MPI_Allreduce(&_n_connect_size,
                  &n_g_connect_size,
                  1,
                  FVMC_MPI_GNUM,
                  MPI_SUM,
                  writer->comm);

    /* Allocate global buffers */

    BFTC_REALLOC(global_vtx_idx_buffer,
                (n_export_elements + n_g_element_section + 1)*export_datasize,
                char);

    fvmc_global_vtx_idx = (fvmc_gnum_t *)global_vtx_idx_buffer;

    /* Compute global buffers */

    polygons_slice =
      fvmc_gather_slice_create(section->global_element_num,
                              n_g_element_section,
                              writer->comm);

    while (fvmc_gather_slice_advance(polygons_slice,
                                    &global_num_start,
                                    &global_num_end) == 0) {

      /*
        Build global vertex connectivity. First, we have to create a global
        index. Then, we gather local vertex connectivity into
        _fvmc_export_connect. We use an offset of n_export_elements as there
        can be several sections of the same type (MED allows only one
        section per type and per mesh).
      */

      fvmc_gather_slice_index(section->vertex_index,
                             fvmc_global_vtx_idx + n_export_elements,
                             section->global_element_num,
                             writer->comm,
                             polygons_slice);

      fvmc_gather_indexed_numbers(section->vertex_index,
                                 section->vertex_num,
                                 _fvmc_export_connect + n_export_connect,
                                 mesh->global_vertex_num,
                                 section->global_element_num,
                                 writer->comm,
                                 fvmc_global_vtx_idx + n_export_elements,
                                 polygons_slice);

    }

    fvmc_gather_slice_destroy(polygons_slice);

    n_export_connect += n_g_connect_size;
    n_export_elements += n_g_element_section + 1;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  assert(n_passes >= 1);

  /* Write buffers into MED file */
  /*-----------------------------*/

  if (writer->rank == 0) {

    int i_count = 0;
    fvmc_gnum_t offset = 1;

    /* Create faces -> vertices index for MED from fvmc_global_vtx_idx */

    n_export_elements = n_export_elements + 1 - n_passes;
    fvmc_global_vtx_idx[0] = 1;

    for (i = 1; i < n_export_elements; i++) {

      if (fvmc_global_vtx_idx[i + i_count] == 0) {
        i_count++;
        offset = fvmc_global_vtx_idx[i - 1];
      }
      fvmc_global_vtx_idx[i] = fvmc_global_vtx_idx[i + i_count] + offset;
    }

    assert(n_passes == i_count + 1);

    med_global_vtx_idx = (med_int *)global_vtx_idx_buffer;
    _convert_fvmc_gnum_to_med_int(fvmc_global_vtx_idx,
                                 med_global_vtx_idx,
                                 n_export_elements);

    /* FVM connectivity to MED connectivity */

    _convert_fvmc_gnum_to_med_int(_fvmc_export_connect,
                                 _med_export_connect,
                                 n_export_connect);

    /* Write polygonal connectivity into MED file */

    retval = MEDpolygoneConnEcr(writer->fid,
                                med_mesh->name,
                                med_global_vtx_idx,
                                (med_int)n_export_elements,
                                _med_export_connect,
                                MED_MAILLE,
                                MED_NOD);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDpolygoneConnEcr() failed to write connectivity:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh_name: \"%s\"\n"),
                writer->name, med_mesh->name, med_section_type);

    /* Write family numbers into MED writer */

    retval  = MEDfamEcr(writer->fid,
                        med_mesh->name,
                        med_family,
                        (med_int)n_export_elements,
                        MED_MAILLE,
                        MED_POLYGONE);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDfamEcr() failed to write family numbers:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh_name: \"%s\"\n"
                  "Associated MED geometrical element: \"%i\"\n"),
                writer->name, med_mesh->name, med_section_type);

  } /* rank == 0 */

  BFTC_FREE(global_vtx_idx_buffer);

  return current_section;
}

#endif /* FVMC_HAVE_MPI */

/*----------------------------------------------------------------------------
 * Write polygonal element connectivity to a MED file in serial mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   med_family            <-- med family number associated.
 *   med_export_connect    <-- buffer to export connectivity.
 *----------------------------------------------------------------------------*/

static const fvmc_writer_section_t *
_export_nodal_polygons_l(const fvmc_writer_section_t  *export_sections,
                         fvmc_to_med_writer_t         *writer,
                         fvmc_to_med_mesh_t           *med_mesh,
                         med_int                     *med_family,
                         med_int                     *med_export_connect)
{
  int i_count;
  fvmc_gnum_t  i;
  med_geometrie_element  med_section_type;

  int   n_passes = 0;
  fvmc_gnum_t  n_export_connect = 0;
  fvmc_gnum_t  n_export_elements = 0;
  med_int   offset = 1;
  med_int  *med_global_vtx_idx = NULL;

  const fvmc_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  /* Get MED element type */

  med_section_type = _get_med_elt_type(current_section->type);

  assert(med_section_type == MED_POLYGONE);
  assert(writer->discard_polygons == false);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvmc_nodal_section_t *section = current_section->section;
    const fvmc_gnum_t n_connect_size = (fvmc_gnum_t)section->connectivity_size;
    const fvmc_gnum_t n_elements_section
      = fvmc_nodal_section_n_g_elements(section);

      n_passes++;

      /* Allocate global buffers */

      BFTC_REALLOC(med_global_vtx_idx,
                  n_export_elements + n_elements_section + 1,
                  med_int);

      /*
        Build global vertex connectivity. First, we have to create a global
        index: global_vertex_index. Then, we gather local vertex connectivity
        into med_export_connect. We have to use an offset of n_export_elements
        as there can be several sections of the same type used in fvmc_nodal_t.
      */

      for (i = 0; i < n_elements_section + 1; i++)
        med_global_vtx_idx[i + n_export_elements] =
          (med_int)section->vertex_index[i];

      for (i = 0; i < n_connect_size; i++)
        med_export_connect[i + n_export_connect] =
          (med_int)section->vertex_num[i];

      n_export_connect += n_connect_size;
      n_export_elements += n_elements_section + 1;

      current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  assert(n_passes >= 1);

  /* Concatenate med_global_vtx_idx */

  i_count = 0;
  n_export_elements = n_export_elements + 1 - n_passes;
  med_global_vtx_idx[0] = 1;

  for (i = 1; i < n_export_elements; i++) {

    if (med_global_vtx_idx[i + i_count] == 0) {
      i_count++;
      offset = med_global_vtx_idx[i - 1];
    }
    med_global_vtx_idx[i] = med_global_vtx_idx[i + i_count] + offset;
  }

  assert(n_passes == i_count + 1);

  /* Write buffers into MED file */
  /*-----------------------------*/

  /* Write polygonal connectivity into MED file */

  retval = MEDpolygoneConnEcr(writer->fid,
                              med_mesh->name,
                              med_global_vtx_idx,
                              (med_int)n_export_elements,
                              med_export_connect,
                              MED_MAILLE,
                              MED_NOD);
  if (retval < 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("MEDpolygoneConnEcr() failed to write connectivity:\n"
                "Associated writer: \"%s\"\n"
                "Associated med_mesh_name: \"%s\"\n"),
              writer->name, med_mesh->name, med_section_type);

  /* Write family numbers into MED writer */

  retval  = MEDfamEcr(writer->fid,
                      med_mesh->name,
                      med_family,
                      (med_int)n_export_elements,
                      MED_MAILLE,
                      MED_POLYGONE);

  if (retval < 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("MEDfamEcr() failed to write family numbers:\n"
                "Associated writer: \"%s\"\n"
                "Associated med_mesh_name: \"%s\"\n"
                "Associated MED geometrical element: \"%i\"\n"),
              writer->name, med_mesh->name, med_section_type);

  BFTC_FREE(med_global_vtx_idx);

  return current_section;
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polyhedral elements connectivity to a MED file in parallel mode
 *
 * parameters:
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   writer                <-- pointer to associated writer.
 *   med_elt_type          <-- MED geometrical element type to export.
 *   n_g_element_section   <-- global number of elements per section.
 *   med_family            <-- med family number associated.
 *   export_connect        <-- buffer to export connectivity.
 *----------------------------------------------------------------------------*/

static const fvmc_writer_section_t *
_export_nodal_polyhedra_g(const fvmc_writer_section_t  *export_sections,
                          fvmc_to_med_writer_t         *writer,
                          const fvmc_nodal_t           *mesh,
                          fvmc_to_med_mesh_t           *med_mesh,
                          med_int                     *med_family,
                          char                        *export_connect)
{
  fvmc_lnum_t  i_elt, vtx_id, face_id;
  fvmc_lnum_t  i_face_idx, cell_vtx_length;
  fvmc_gnum_t  i;
  med_geometrie_element  med_section_type;

  fvmc_lnum_t  i_face = 0, i_connect = 0;
  fvmc_gnum_t  global_num_start = 0, global_num_end = 0;
  fvmc_gnum_t  n_export_connect = 0;
  fvmc_gnum_t  n_export_elements = 0;
  fvmc_gnum_t  n_export_faces = 0;

  char       *global_face_lengths_buffer = NULL;
  fvmc_gnum_t *fvmc_global_face_lengths = NULL;
  med_int    *med_export_faces_index = NULL;

  char       *global_cell_lengths_buffer = NULL;
  fvmc_gnum_t *fvmc_global_cell_lengths = NULL;
  med_int    *med_export_vertices_index = NULL;

  fvmc_gnum_t *fvmc_export_connect = (fvmc_gnum_t *)export_connect;
  med_int    *med_export_connect = (med_int *)export_connect;

  fvmc_gather_slice_t *polyhedra_slice = NULL;

  const size_t export_datasize = FVMC_MAX(sizeof(fvmc_gnum_t),sizeof(med_int));
  const fvmc_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  /* Get MED element type */

  med_section_type = _get_med_elt_type(current_section->type);

  assert(med_section_type == MED_POLYEDRE);
  assert(writer->discard_polyhedra == false);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    fvmc_lnum_t  n_l_faces_section = 0, n_l_connect_section = 0;
    fvmc_gnum_t  n_g_faces_section = 0, n_g_connect_section = 0;

    fvmc_lnum_t *_face_lengths = NULL;
    fvmc_lnum_t *_cell_vtx_idx = NULL;
    fvmc_lnum_t *_cell_connect = NULL;
    fvmc_gnum_t *_cell_lengths = NULL;
    fvmc_gnum_t *global_cell_vtx_idx = NULL;
    fvmc_gnum_t *global_cell_face_idx = NULL;

    const fvmc_nodal_section_t *section = current_section->section;
    const fvmc_gnum_t  n_g_elements_section
      = fvmc_nodal_section_n_g_elements(section);

    n_l_faces_section = section->face_index[section->n_elements];
    MPI_Allreduce(&n_l_faces_section,
                  &n_g_faces_section,
                  1,
                  FVMC_MPI_LNUM,
                  MPI_SUM,
                  writer->comm);

    /*
      Build locally:
       - face_lengths (number of vertex per face),
       - cell_lengths (number of faces per cell) and
       - cells -> vertices index for each polyhedral section
    */

    BFTC_MALLOC(_face_lengths, n_l_faces_section + 1, fvmc_lnum_t);
    BFTC_MALLOC(_cell_lengths, section->n_elements, fvmc_gnum_t);
    BFTC_MALLOC(_cell_vtx_idx, section->n_elements + 1, fvmc_lnum_t);

    _cell_vtx_idx[0] = 0;

    for (i_elt = 0; i_elt < section->n_elements; i_elt++) {

      _cell_lengths[i_elt] = (fvmc_gnum_t)(  section->face_index[i_elt + 1]
                                          - section->face_index[i_elt]);
      cell_vtx_length = 0;

      for (i_face_idx = section->face_index[i_elt];
           i_face_idx < section->face_index[i_elt + 1];
           i_face_idx++) {

        face_id = FVMC_ABS(section->face_num[i_face_idx]) - 1;
        _face_lengths[i_face] = (  section->vertex_index[face_id + 1]
                                 - section->vertex_index[face_id]);
        cell_vtx_length += _face_lengths[i_face];
        i_face++;
      }

      _cell_vtx_idx[i_elt + 1] = _cell_vtx_idx[i_elt] + cell_vtx_length;
    }

    /* Build locally _cell_connect */

    n_l_connect_section = _cell_vtx_idx[section->n_elements];
    BFTC_MALLOC(_cell_connect, n_l_connect_section, fvmc_lnum_t);

    i_connect = 0;

    for (i_elt = 0; i_elt < section->n_elements; i_elt++) {
      for (i_face_idx = section->face_index[i_elt];
           i_face_idx < section->face_index[i_elt + 1];
           i_face_idx++) {

        if (section->face_num[i_face_idx] > 0) {
          face_id = section->face_num[i_face_idx] - 1;
          for (vtx_id = section->vertex_index[face_id];
               vtx_id < section->vertex_index[face_id + 1];
               vtx_id++)
            _cell_connect[i_connect++] = section->vertex_num[vtx_id];
        }
        else { /* face_num < 0 */
          face_id = -section->face_num[i_face_idx] - 1;
          vtx_id = section->vertex_index[face_id];
          _cell_connect[i_connect++] = section->vertex_num[vtx_id];
          for (vtx_id = section->vertex_index[face_id + 1] - 1;
               vtx_id > section->vertex_index[face_id];
               vtx_id--)
            _cell_connect[i_connect++] = section->vertex_num[vtx_id];
        }
      }
    }

    /* Compute global connectivity size */

    MPI_Allreduce(&n_l_connect_section,
                  &n_g_connect_section,
                  1,
                  FVMC_MPI_LNUM,
                  MPI_SUM,
                  writer->comm);

    /* Allocate global buffers to export MED connectivity */

    BFTC_REALLOC(global_cell_lengths_buffer,
                (n_export_elements + n_g_elements_section + 1)*export_datasize,
                char);

    BFTC_REALLOC(global_face_lengths_buffer,
                (n_export_faces + 1 + n_g_faces_section)*export_datasize,
                char);

    fvmc_global_cell_lengths = (fvmc_gnum_t *)global_cell_lengths_buffer;
    fvmc_global_face_lengths = (fvmc_gnum_t *)global_face_lengths_buffer;

    /* Allocate global buffers used in "fvmc_gather" operations */

    BFTC_MALLOC(global_cell_face_idx, n_g_elements_section + 1, fvmc_gnum_t);
    BFTC_MALLOC(global_cell_vtx_idx, n_g_elements_section + 1, fvmc_gnum_t);

    /* Compute global buffers */

    polyhedra_slice =
      fvmc_gather_slice_create(section->global_element_num,
                              n_g_elements_section,
                              writer->comm);

    while (fvmc_gather_slice_advance(polyhedra_slice,
                                    &global_num_start,
                                    &global_num_end) == 0) {

      /*
        Build global cell lengths: offset of n_export_elements because
        there can be several sections implied in the export.
        Offset of 1 because we have one more element when we transform
        cell_lengths into index.
      */

      fvmc_gather_array(_cell_lengths,
                       fvmc_global_cell_lengths + 1 + n_export_elements,
                       FVMC_MPI_GNUM,
                       1,
                       section->global_element_num,
                       writer->comm,
                       polyhedra_slice);

      /* Build global cells -> faces index used in gather_indexed_numbers */

      fvmc_gather_slice_index(section->face_index,
                             global_cell_face_idx,
                             section->global_element_num,
                             writer->comm,
                             polyhedra_slice);

      /*
        Build global face lengths: offset of n_export_faces because several
        sections can be implied in the export and offset of 1 because we
        have one more element when we transform face_lengths into an index.
      */

      fvmc_gather_indexed_numbers(section->face_index,
                                 _face_lengths,
                                 fvmc_global_face_lengths + 1 + n_export_faces,
                                 NULL,
                                 section->global_element_num,
                                 writer->comm,
                                 global_cell_face_idx,
                                 polyhedra_slice);

      /*
        Build global vertex connectivity. First, we have to create a global
        index. Then, we gather local vertex connectivity into
        fvmc_export_connect
      */

      fvmc_gather_slice_index(_cell_vtx_idx,
                             global_cell_vtx_idx,
                             section->global_element_num,
                             writer->comm,
                             polyhedra_slice);

      fvmc_gather_indexed_numbers(_cell_vtx_idx,
                                 _cell_connect,
                                 fvmc_export_connect + n_export_connect,
                                 mesh->global_vertex_num,
                                 section->global_element_num,
                                 writer->comm,
                                 global_cell_vtx_idx,
                                 polyhedra_slice);

    }

    BFTC_FREE(_face_lengths);
    BFTC_FREE(_cell_lengths);
    BFTC_FREE(_cell_connect);
    BFTC_FREE(_cell_vtx_idx);

    BFTC_FREE(global_cell_face_idx);
    BFTC_FREE(global_cell_vtx_idx);

    fvmc_gather_slice_destroy(polyhedra_slice);

    n_export_connect += n_g_connect_section;
    n_export_elements += n_g_elements_section;
    n_export_faces += n_g_faces_section;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  /* Write buffers into MED file */
  /*-----------------------------*/

  if (writer->rank == 0) {

    /* Create MED cells -> faces index from fvmc_global_cell_lengths. */

    fvmc_global_cell_lengths[0] = 1;
    for (i = 1; i < n_export_elements + 1; i++)
      fvmc_global_cell_lengths[i] =
        fvmc_global_cell_lengths[i] + fvmc_global_cell_lengths[i-1];

    med_export_faces_index = (med_int *)fvmc_global_cell_lengths;
    _convert_fvmc_gnum_to_med_int(fvmc_global_cell_lengths,
                                 med_export_faces_index,
                                 n_export_elements);

    /* Create MED faces -> vertices index from fvmc_global_face_lengths. */

    fvmc_global_face_lengths[0] = 1;
    for (i = 1; i < n_export_faces + 1; i++)
      fvmc_global_face_lengths[i] =
        fvmc_global_face_lengths[i] + fvmc_global_face_lengths[i-1];

    med_export_vertices_index = (med_int *)fvmc_global_face_lengths;
    _convert_fvmc_gnum_to_med_int(fvmc_global_face_lengths,
                                 med_export_vertices_index,
                                 n_export_faces + 1);

    /* FVM connectivity to MED connectivity */

    _convert_fvmc_gnum_to_med_int(fvmc_export_connect,
                                 med_export_connect,
                                 n_export_connect);

    /* Write polyhedral connectivity into MED file */

    retval = MEDpolyedreConnEcr(writer->fid,
                                med_mesh->name,
                                med_export_faces_index,
                                (med_int)n_export_elements + 1,
                                med_export_vertices_index,
                                (med_int)n_export_faces + 1,
                                med_export_connect,
                                MED_NOD);
    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDpolyedreConnEcr() failed to write connectivity:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh_name: \"%s\"\n"),
                writer->name, med_mesh->name, med_section_type);

    /* Write family numbers into MED writer */

    retval  = MEDfamEcr(writer->fid,
                        med_mesh->name,
                        med_family,
                        (med_int)n_export_elements,
                        MED_MAILLE,
                        MED_POLYEDRE);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDfamEcr() failed to write family numbers:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh_name: \"%s\"\n"
                  "Associated MED geometrical element: \"%i\"\n"),
                writer->name, med_mesh->name, med_section_type);

  } /* rank == 0 */

  BFTC_FREE(global_cell_lengths_buffer);
  BFTC_FREE(global_face_lengths_buffer);

  return current_section;
}

#endif /* FVMC_HAVE_MPI */

/*----------------------------------------------------------------------------
 * Write polyhedral elements connectivity to a MED file in serial mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   med_elt_type          <-- MED geometrical element type to export.
 *   n_g_element_section   <-- global number of elements per section.
 *   med_family            <-- med family number associated.
 *   export_connect        <-- buffer to export connectivity.
 *----------------------------------------------------------------------------*/

static const fvmc_writer_section_t *
_export_nodal_polyhedra_l(const fvmc_writer_section_t  *export_sections,
                          fvmc_to_med_writer_t         *writer,
                          fvmc_to_med_mesh_t           *med_mesh,
                          med_int                     *med_family,
                          med_int                     *med_export_connect)
{
  fvmc_lnum_t  i_face_idx, face_id, vtx_id, i_elt;
  fvmc_gnum_t  i;
  med_geometrie_element  med_section_type;

  fvmc_gnum_t  i_cell = 1, i_face = 1, i_connect = 0;
  fvmc_gnum_t  n_export_elements = 0;
  fvmc_gnum_t  n_export_faces = 0;

  med_int  *med_face_lengths = NULL;
  med_int  *med_cell_lengths = NULL;

  const fvmc_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  /* Get MED element type */

  med_section_type = _get_med_elt_type(current_section->type);

  assert(med_section_type == MED_POLYEDRE);
  assert(writer->discard_polyhedra == false);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvmc_nodal_section_t *section = current_section->section;

    fvmc_gnum_t  n_faces_section = 0;

    n_faces_section = section->face_index[section->n_elements];

    BFTC_REALLOC(med_face_lengths,
                n_faces_section + 1 + n_export_faces,
                med_int);

    BFTC_REALLOC(med_cell_lengths,
                section->n_elements + 1 + n_export_elements,
                med_int);

    /*
      Build locally:
      - med_face_lengths (number of vertex per face),
      - med_cell_lengths (number of faces per cell) and
      - cells -> vertices connectivity for each polyhedral section
    */

    for (i_elt = 0; i_elt < section->n_elements; i_elt++) {

      med_cell_lengths[i_cell++] = (  (med_int)section->face_index[i_elt + 1]
                                    - (med_int)section->face_index[i_elt]);

      for (i_face_idx = section->face_index[i_elt];
           i_face_idx < section->face_index[i_elt + 1];
           i_face_idx++) {

        face_id = FVMC_ABS(section->face_num[i_face_idx]) - 1;

        med_face_lengths[i_face] = (  (med_int)section->vertex_index[face_id + 1]
                                    - (med_int)section->vertex_index[face_id]);
        i_face++;

        if (section->face_num[i_face_idx] > 0) {
          for (vtx_id = section->vertex_index[face_id];
               vtx_id < section->vertex_index[face_id + 1];
               vtx_id++)
            med_export_connect[i_connect++] =
              (med_int)section->vertex_num[vtx_id];
          }
          else { /* face_num < 0 */
            vtx_id = section->vertex_index[face_id];
            med_export_connect[i_connect++ ] =
              (med_int)section->vertex_num[vtx_id];

            for (vtx_id = section->vertex_index[face_id + 1] - 1;
                 vtx_id > section->vertex_index[face_id];
                 vtx_id--)
              med_export_connect[i_connect++] =
                (med_int)section->vertex_num[vtx_id];
          }

      } /* End of loop on faces */

    } /* End of loop on cells */

    n_export_elements += section->n_elements;
    n_export_faces += n_faces_section;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  /* Write buffers into MED file */
  /*-----------------------------*/

  /* Create MED cells -> faces index from cell_lengths. */

  med_cell_lengths[0] = 1;

  for (i = 1; i < n_export_elements + 1; i++)
    med_cell_lengths[i] = med_cell_lengths[i] + med_cell_lengths[i-1];

  /* Create MED faces -> vertices index to export from face_lengths */

  med_face_lengths[0] = 1;

  for (i = 1; i < n_export_faces + 1; i++)
    med_face_lengths[i] = med_face_lengths[i] + med_face_lengths[i-1];

  /* Write polyhedral connectivity into MED file */

  retval = MEDpolyedreConnEcr(writer->fid,
                              med_mesh->name,
                              med_cell_lengths,
                              (med_int)n_export_elements + 1,
                              med_face_lengths,
                              (med_int)n_export_faces + 1,
                              med_export_connect,
                              MED_NOD);

  if (retval < 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("MEDpolyedreConnEcr() failed to write connectivity:\n"
                "Associated writer: \"%s\"\n"
                "Associated med_mesh_name: \"%s\"\n"),
              writer->name, med_mesh->name, med_section_type);

  /* Write family numbers into MED writer */

  retval  = MEDfamEcr(writer->fid,
                      med_mesh->name,
                      med_family,
                      (med_int)n_export_elements,
                      MED_MAILLE,
                      MED_POLYEDRE);

  if (retval < 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("MEDfamEcr() failed to write family numbers:\n"
                "Associated writer: \"%s\"\n"
                "Associated med_mesh_name: \"%s\"\n"
                "Associated MED geometrical element: \"%i\"\n"),
              writer->name, med_mesh->name, med_section_type);

  BFTC_FREE(med_cell_lengths);
  BFTC_FREE(med_face_lengths);

  return current_section;
}

/*----------------------------------------------------------------------------
 * Write per element field values into MED writer.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   export_section     <-- pointer to section helper structure
 *   helper             <-- pointer to general writer helper structure
 *   writer             <-- pointer to associated writer.
 *   input_dim          <-- input field dimension
 *   output_dim         <-- output field dimension
 *   interlace          <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   med_mesh_name      <-- MED name of the mesh on which the field is defined
 *   med_field_name     <-- MED name of the field to export.
 *   time_step          <-- number of the current time step
 *   time_value         <-- associated time value
 *   output_buffer_size <-- minimum output buffer size
 *   output_buffer      --- output buffer (size output_buffer_size on ranks
 *                          > 0, full output size on rank 0)
 *
 * returns:
 *  pointer to next section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvmc_writer_section_t *
_export_field_values_e(const fvmc_writer_section_t      *export_section,
                       fvmc_writer_field_helper_t       *helper,
                       fvmc_to_med_writer_t             *writer,
                       int                              input_dim,
                       int                              output_dim,
                       fvmc_interlace_t                  interlace,
                       int                              n_parent_lists,
                       const fvmc_lnum_t                 parent_num_shift[],
                       fvmc_datatype_t                   datatype,
                       const void                *const field_values[],
                       char                            *med_mesh_name,
                       char                            *med_field_name,
                       int                              time_step,
                       double                           time_value,
                       size_t                           output_buffer_size,
                       unsigned char                    output_buffer[])
{
  med_geometrie_element  med_section_type;
  med_err  retval = 0;

  size_t  datatype_size = 0;
  size_t  output_size = 0;
  size_t  output_size_tot = 0;
  unsigned char  *output_buffer_slice = output_buffer;

  const fvmc_writer_section_t  *current_section = export_section;

  _Bool loop_on_sections = true;

  med_section_type = _get_med_elt_type(current_section->type);
  datatype_size
    = fvmc_datatype_size[fvmc_writer_field_helper_datatype(helper)];

  /* Loop on sections of same type to fill output buffer */
  /*-----------------------------------------------------*/

  while (loop_on_sections == true) {

    while (fvmc_writer_field_helper_step_e(helper,
                                          current_section,
                                          input_dim,
                                          0,
                                          interlace,
                                          n_parent_lists,
                                          parent_num_shift,
                                          datatype,
                                          field_values,
                                          output_buffer_slice,
                                          output_buffer_size,
                                          &output_size) == 0) {

      output_size_tot += output_size;
      output_buffer_slice =   output_buffer + (output_size_tot * datatype_size);

    }

    current_section = current_section->next;

    if (   current_section == NULL
        || current_section->continues_previous == false)
      loop_on_sections = false;

  } /* while (loop on sections) */

  /* Write buffer to MED file */
  /*--------------------------*/

  if (writer->rank == 0) {

    char med_unit_name[] = MED_PNOM_BLANC;
    char med_nopfl[] = MED_NOPFL;
    char med_nogauss[] = MED_NOGAUSS;

    retval = MEDchampEcr(writer->fid,
                         med_mesh_name,
                         med_field_name,
                         output_buffer,
                         MED_FULL_INTERLACE,
                         output_size_tot / output_dim,
                         med_nogauss,
                         MED_ALL,
                         med_nopfl,
                         MED_NO_PFLMOD,
                         MED_MAILLE,
                         med_section_type,
                         time_step,
                         med_unit_name,
                         time_value,
                         MED_NONOR);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                "_export_field_e() failed to write per element field values\n"
                "Associated fieldname: \"%s\"\n"
                "Associated med mesh: \"%s\"\n"
                "Associated writer name: \"%s\"\n",
                med_field_name, med_mesh_name, writer->name);

  }

  return current_section;
}

/*----------------------------------------------------------------------------
 * Write per node field values into MED writer.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   export_section     <-- pointer to section helper structure
 *   helper             <-- pointer to general writer helper structure
 *   writer             <-- pointer to associated writer.
 *   input_dim          <-- input field dimension
 *   output_dim         <-- output field dimension
 *   interlace          <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   med_mesh_name      <-- MED name of the mesh on which the field is defined
 *   med_field_name     <-- MED name of the field to export.
 *   time_step          <-- number of the current time step
 *   time_value         <-- associated time value
 *   output_buffer_size <-- minimum output buffer size
 *   output_buffer      --- output buffer (size output_buffer_size on ranks
 *                          > 0, full output size on rank 0)
 *
 * returns:
 *  pointer to next section helper structure in list
 *----------------------------------------------------------------------------*/

static void
_export_field_values_n(const fvmc_nodal_t               *mesh,
                       fvmc_writer_field_helper_t       *helper,
                       fvmc_to_med_writer_t             *writer,
                       int                              input_dim,
                       int                              output_dim,
                       fvmc_interlace_t                  interlace,
                       int                              n_parent_lists,
                       const fvmc_lnum_t                 parent_num_shift[],
                       fvmc_datatype_t                   datatype,
                       const void                *const field_values[],
                       char                            *med_mesh_name,
                       char                            *med_field_name,
                       int                              time_step,
                       double                           time_value,
                       size_t                           output_buffer_size,
                       unsigned char                    output_buffer[])
{
  med_err  retval = 0;

  size_t  output_size = 0;
  size_t  output_size_tot = 0;
  unsigned char  *output_buffer_slice = output_buffer;

  size_t datatype_size
    = fvmc_datatype_size[fvmc_writer_field_helper_datatype(helper)];

  /* Fill output buffer */
  /*--------------------*/

  while (fvmc_writer_field_helper_step_n(helper,
                                        mesh,
                                        input_dim,
                                        0,
                                        interlace,
                                        n_parent_lists,
                                        parent_num_shift,
                                        datatype,
                                        field_values,
                                        output_buffer_slice,
                                        output_buffer_size,
                                        &output_size) == 0) {

    output_size_tot += output_size;
    output_buffer_slice =   output_buffer + (output_size_tot * datatype_size);

  }

  /* Write buffer to MED file */
  /*--------------------------*/

  if (writer->rank == 0) {

    char med_unit_name[] = MED_PNOM_BLANC;
    char med_nopfl[] = MED_NOPFL;
    char med_nogauss[] = MED_NOGAUSS;

    retval = MEDchampEcr(writer->fid,
                         med_mesh_name,
                         med_field_name,
                         output_buffer,
                         MED_FULL_INTERLACE,
                         output_size_tot / output_dim,
                         med_nogauss,
                         MED_ALL,
                         med_nopfl,
                         MED_NO_PFLMOD,
                         MED_NOEUD,
                         MED_POINT1,
                         time_step,
                         med_unit_name,
                         time_value,
                         MED_NONOR);

    if (retval < 0)
      bftc_error(__FILE__, __LINE__, 0,
                "_export_field_e() failed to write per node field values\n"
                "Associated fieldname: \"%s\"\n"
                "Associated med mesh: \"%s\"\n"
                "Associated writer name: \"%s\"\n",
                med_field_name, med_mesh_name, writer->name);

  }

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to MED file writer.
 *
 * Options are:
 *   discard_polygons    do not output polygons or related values
 *   discard_polyhedra   do not output polyhedra or related values
 *   divide_polygons     tesselate polygons with triangles
 *   divide_polyhedra    tesselate polyhedra with tetrahedra and pyramids
 *                       (adding a vertex near each polyhedron's center)
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque MED writer structure.
 *----------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)
void *
fvmc_to_med_init_writer(const char                   *name,
                       const char                   *path,
                       const char                   *options,
                       const fvmc_writer_time_dep_t   time_dependency,
                       const MPI_Comm                comm)
#else
void *
fvmc_to_med_init_writer(const char                   *name,
                       const char                   *path,
                       const char                   *options,
                       const fvmc_writer_time_dep_t   time_dependency)
#endif
{
  fvmc_to_med_writer_t  *writer = NULL;

  int  i, fid = 0;
  int  filename_length, name_length, path_length;

  /* Initialize writer */

  BFTC_MALLOC(writer, 1, fvmc_to_med_writer_t);

  writer->n_med_meshes = 0;
  writer->n_fields  = 0;
  writer->med_meshes   = NULL;
  writer->fields = NULL;

  writer->n_time_steps   = 0;
  writer->time_steps     = NULL;
  writer->time_values    = NULL;
  writer->time_dependency = time_dependency;

  /* Parallel parameters */

  writer->rank = 0;
  writer->n_ranks = 1;

#if defined(FVMC_HAVE_MPI)
  {
    int mpi_flag, rank, n_ranks;
    MPI_Initialized(&mpi_flag);

    if (mpi_flag && comm != MPI_COMM_NULL) {
      writer->comm = comm;
      MPI_Comm_rank(writer->comm, &rank);
      MPI_Comm_size(writer->comm, &n_ranks);
      writer->rank = rank;
      writer->n_ranks = n_ranks;
    }
    else
      writer->comm = MPI_COMM_NULL;
  }
#endif /* defined(FVMC_HAVE_MPI) */

  /* Writer name and file name */

  name_length = strlen(name);
  if (name_length == 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("No MED filename: \"%s\"\n"),
              *name);

  BFTC_MALLOC(writer->name, name_length + 1, char);
  strcpy(writer->name, name);

  for (i = 0; i < name_length; i++) {
    if (writer->name[i] == ' ' || writer->name[i] == '\t')
      writer->name[i] = '_';
  }

  if (path != NULL)
    path_length = strlen(path);
  else
    path_length = 0;
  filename_length = path_length + name_length + strlen(".med");
  BFTC_MALLOC(writer->filename, filename_length + 1, char);

  if (path != NULL)
    strcpy(writer->filename, path);
  else
    writer->filename[0] = '\0';

  strcat(writer->filename, writer->name);
  strcat(writer->filename, ".med");

  writer->filename[filename_length] = '\0';
  writer->name[name_length] = '\0';

  /* Reading options */

  writer->discard_polygons = false;
  writer->discard_polyhedra = false;
  writer->divide_polygons = false;
  writer->divide_polyhedra = false;

  if (options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1 ; i2 < l_tot && options[i2] != ' ' ; i2++);
      l_opt = i2 - i1;

      if (        (l_opt == 16)
               && (strncmp(options + i1, "discard_polygons", l_opt) == 0))
        writer->discard_polygons = true;

      else if (   (l_opt == 17)
               && (strncmp(options + i1, "discard_polyhedra", l_opt) == 0))
        writer->discard_polyhedra = true;

      else if (   (l_opt == 15)
               && (strncmp(options + i1, "divide_polygons", l_opt) == 0))
        writer->divide_polygons = true;

      else if (   (l_opt == 16)
               && (strncmp(options + i1, "divide_polyhedra", l_opt) == 0))
        writer->divide_polyhedra = true;

      for (i1 = i2 + 1 ; i1 < l_tot && options[i1] == ' ' ; i1++);

    }

  }

  /* Open MED file */

  writer->is_open = false;
  writer->fid = 0;

  if (writer->rank == 0) {

    fid = MEDouvrir(writer->filename, MED_CREATION);
    if (fid < 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDouvrir() failed to open file: %s\n"
                  "Associated writer: \"%s\"\n"),
                writer->filename, writer->name);

  }

  writer->is_open = true;

#if defined(FVMC_HAVE_MPI)
  if (writer->n_ranks > 1)
    MPI_Bcast(&fid, 1, MPI_INT, 0, writer->comm);
#endif

  writer->fid = fid;
  return writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to MED file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque MED writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvmc_to_med_finalize_writer(void  *this_writer_p)
{
  int i;

  fvmc_to_med_writer_t  *writer
                        = (fvmc_to_med_writer_t *)this_writer_p;

  const int rank = writer->rank;

  assert(writer != NULL);

  if (rank == 0) {

    /* Close MED File */
    /*----------------*/

    if (writer->is_open == true) {

      if (MEDfermer(writer->fid) != 0)
        bftc_error(__FILE__, __LINE__, 0,
                  _("MEDfermer() failed to close file \"%s\"\n"),
                  writer->filename);

      writer->fid = 0;

    }

  } /* End if rank = 0 */

  writer->is_open = false;

  /* Free structures */
  /*-----------------*/

  BFTC_FREE(writer->name);
  BFTC_FREE(writer->filename);
  BFTC_FREE(writer->time_values);
  BFTC_FREE(writer->time_steps);

  /* Free fvmc_to_med_mesh_t structure */

  for (i = 0; i < writer->n_med_meshes; i++)
    writer->med_meshes[i] = BFTC_FREE(writer->med_meshes[i]);
  BFTC_FREE(writer->med_meshes);

  for (i = 0; i < writer->n_fields; i++)
    BFTC_FREE(writer->fields[i]);

  BFTC_FREE(writer->fields);

  /* Free fvmc_to_med_writer_t structure */

  BFTC_FREE(writer);
  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with a MED mesh.
 *
 * parameters:
 *   this_writer <-- pointer to associated writer
 *   time_step   <-- time step number
 *   time_value  <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvmc_to_med_set_mesh_time(void          *this_writer,
                         int            time_step,
                         double         time_value)
{
  int n_vals;
  fvmc_to_med_writer_t  *writer = (fvmc_to_med_writer_t *)this_writer;

  const char time_value_err_string[] =
    N_("The time value associated with time step <%d> equals <%g>,\n"
       "but time value <%g> has already been associated with this time step.\n");

  assert(writer != NULL);

  /* First verification on time step */

  if (time_step < 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("The given time step value should be >= 0, and not %d\n"),
              time_step);

  if (   writer->time_steps != NULL
      && writer->time_values != NULL) {

    n_vals = writer->n_time_steps;
    if (time_step < writer->time_steps[n_vals - 1])
      bftc_error(__FILE__, __LINE__, 0,
                _("The given time step value should be >= %d, and not %d\n"),
                writer->time_steps[n_vals - 1], time_step);

    /* Verifications on time value */

    else if (time_step == writer->time_steps[n_vals - 1]) {
      if (   time_value < writer->time_values[n_vals - 1] - 1.e-16
          || time_value > writer->time_values[n_vals - 1] + 1.e-16)
        bftc_error(__FILE__, __LINE__, 0,
                  _(time_value_err_string), time_step,
                  time_value, writer->time_values[n_vals - 1]);
    }
    else { /* Add a new time step and time value */
      writer->n_time_steps += 1;
      n_vals = writer->n_time_steps;

      BFTC_REALLOC(writer->time_values, n_vals, double);
      BFTC_REALLOC(writer->time_steps, n_vals, int);

      writer->time_values[n_vals - 1] = time_value;
      writer->time_steps[n_vals - 1] = time_step;
    }
  }
  else { /* Setting of the first time step and time value */
    writer->n_time_steps += 1;
    n_vals = writer->n_time_steps;

    BFTC_REALLOC(writer->time_values, n_vals, double);
    BFTC_REALLOC(writer->time_steps, n_vals, int);

    writer->time_values[n_vals - 1] = time_value;
    writer->time_steps[n_vals - 1] = time_step;
  }

  return;
}

/*----------------------------------------------------------------------------
 * Indicate if a elements of a given type in a mesh associated to a given
 * MED file writer need to be tesselated.
 *
 * parameters:
 *   this_writer  <-- pointer to associated writer
 *   mesh         <-- pointer to nodal mesh structure that should be written
 *   element_type <-- element type we are interested in
 *
 * returns:
 *   1 if tesselation of the given element type is needed, 0 otherwise
 *----------------------------------------------------------------------------*/

int
fvmc_to_med_needs_tesselation(fvmc_writer_t       *this_writer,
                             const fvmc_nodal_t  *mesh,
                             fvmc_element_t       element_type)
{
  int  i;
  int  retval = 0;
  fvmc_to_med_writer_t  *writer = (fvmc_to_med_writer_t *)this_writer;

  const int  export_dim = fvmc_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  if (   (   element_type == FVMC_FACE_POLY
          && writer->divide_polygons == true)
      || (   element_type == FVMC_CELL_POLY
          && writer->divide_polyhedra == true)) {

    for (i = 0 ; i < mesh->n_sections ; i++) {

      const fvmc_nodal_section_t  *section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (section->entity_dim == export_dim) {
        if (section->type == element_type)
          retval = 1;
      }

    }

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a MED file
 *
 * parameters:
 *   this_writer  <-- pointer to associated writer.
 *   mesh         <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvmc_to_med_export_nodal(void               *this_writer,
                        const fvmc_nodal_t  *mesh)
{
  int   i_char, n_chars, med_mesh_num;
  fvmc_gnum_t i;
  size_t connect_type_size;
  med_geometrie_element med_type;

  fvmc_gnum_t  global_connect_slice_size = 0;

  char  med_mesh_name[MED_TAILLE_NOM + 1];

  char  *export_connect_buffer = NULL;
  med_int  *med_export_connect = NULL;
  med_int  *med_family = NULL;
  fvmc_to_med_mesh_t  *med_mesh = NULL;
  fvmc_writer_section_t  *export_list = NULL;

  fvmc_to_med_writer_t  *writer = (fvmc_to_med_writer_t *)this_writer;

  const fvmc_writer_section_t  *export_sections = NULL;
  const int  n_ranks = writer->n_ranks;
  const int  rank = writer->rank;

  /* Re-open MED file */

  if (writer->is_open == false) {

    if (writer->rank == 0) {
      writer->fid = MEDouvrir(writer->filename, MED_LECTURE_ECRITURE);
      if (writer->fid < 0)
        bftc_error(__FILE__, __LINE__, 0,
                  _("MEDouvrir() failed to re-open file: %s\n"
                    "Associated writer: \"%s\"\n"),
                  writer->filename, writer->name);
    }
    writer->is_open = true;

  }

  /* Get med_mesh structure */
  /*------------------------*/

  /* Clean mesh->name */

  if (mesh->name == NULL)
    bftc_error(__FILE__, __LINE__, 0,
              _("Mesh name required to continue.\n"));

  strncpy(med_mesh_name, mesh->name, MED_TAILLE_NOM);
  n_chars = strlen(med_mesh_name);

  for (i_char = n_chars + 1; i_char < MED_TAILLE_NOM; i_char++)
    med_mesh_name[i_char] = ' ';
  med_mesh_name[MED_TAILLE_NOM] = '\0';

  /* Get med_mesh_num */

  med_mesh_num = _get_med_mesh_num(writer,
                                   med_mesh_name);

  if (med_mesh_num == 0 )
    med_mesh_num = _add_med_mesh(writer,
                                 med_mesh_name,
                                 mesh);

  med_mesh = writer->med_meshes[med_mesh_num - 1];

  /*---------------------------*/
  /* Export vertex coordinates */
  /*---------------------------*/

#if defined(FVMC_HAVE_MPI)
  if (n_ranks > 1)
    _export_vertex_coords_g(mesh,
                            med_mesh,
                            writer);
#endif /* FVMC_HAVE_MPI */

  if (n_ranks == 1)
    _export_vertex_coords_l(mesh,
                            med_mesh,
                            writer);

  /*---------------------*/
  /* Export connectivity */
  /*---------------------*/

  /* Build list of sections that are used here, in order of output */

  export_list = fvmc_writer_export_list(mesh,
                                       true,
                                       writer->discard_polygons,
                                       writer->discard_polyhedra,
                                       writer->divide_polygons,
                                       writer->divide_polyhedra);

  if (export_list == NULL)
    bftc_error(__FILE__, __LINE__, 0,
              _("MED must have a connectivity.\n"
                "Mesh: \"%s\"\n"
                "Writer: \"%s\"\n"),
              med_mesh->name, writer->name);

  /* Compute connectivity buffer size */

  export_sections = export_list;

  global_connect_slice_size = _get_connect_buffer_size(writer,
                                                       export_sections);

  /* Allocate connectivity buffer */

#if defined(FVMC_HAVE_MPI)
  if (n_ranks > 1) {

    connect_type_size = FVMC_MAX(sizeof(fvmc_gnum_t), sizeof(med_int));
    BFTC_MALLOC(export_connect_buffer,
               global_connect_slice_size * connect_type_size,
               char);

  }
#endif /* FVMC_HAVE_MPI */

  if (n_ranks == 1)
    BFTC_MALLOC(med_export_connect, global_connect_slice_size, med_int);

  /* MED Family numbers treatment */
  /*------------------------------*/

  if (rank == 0) {

    BFTC_MALLOC(med_family, global_connect_slice_size, med_int);
    for (i = 0; i < global_connect_slice_size; i++)
      med_family[i] = 0;

  }

  /* Loop on MED element types */
  /*---------------------------*/

  while (export_sections != NULL) {

    med_type = _get_med_elt_type(export_sections->type);

#if defined(FVMC_HAVE_MPI)
    if (n_ranks > 1) { /* Parallel treatment */

      if (med_type == MED_POLYGONE) {

        if (writer->discard_polygons == false)
          export_sections =
            _export_nodal_polygons_g(export_sections,
                                     writer,
                                     mesh,
                                     med_mesh,
                                     med_family,
                                     export_connect_buffer);
        else
          /* discard_polygons == true */
          break;

      }
      else if (med_type == MED_POLYEDRE) {

        if (writer->discard_polyhedra == false)
          export_sections =
            _export_nodal_polyhedra_g(export_sections,
                                      writer,
                                      mesh,
                                      med_mesh,
                                      med_family,
                                      export_connect_buffer);
        else
          /* discard_polyhedra == true */
          break;

      }
      else
        export_sections =
          _export_connect_g(export_sections,
                            writer,
                            mesh,
                            med_mesh,
                            med_family,
                            export_connect_buffer);

    }
#endif /* FVMC_HAVE_MPI */

    if (n_ranks == 1) { /* Serial treatment */

      if (med_type == MED_POLYGONE) {

        if (writer->discard_polygons == false)
          export_sections =
            _export_nodal_polygons_l(export_sections,
                                     writer,
                                     med_mesh,
                                     med_family,
                                     med_export_connect);
        else
          /* discard_polygons == true */
          break;

      }
      else if (med_type == MED_POLYEDRE) {

        if (writer->discard_polyhedra == false)
          export_sections =
            _export_nodal_polyhedra_l(export_sections,
                                      writer,
                                      med_mesh,
                                      med_family,
                                      med_export_connect);
        else
          /* discard_polyhedra == true */
          break;

      }
      else
        export_sections =
          _export_connect_l(export_sections,
                            writer,
                            med_mesh,
                            med_family,
                            med_export_connect);

    } /* end of serial treatment */

  } /* End of loop on med_types */

  /* Free buffers */

  if (rank == 0)
    BFTC_FREE(med_family);

  if (med_export_connect != NULL)
    BFTC_FREE(med_export_connect);

  if (export_connect_buffer != NULL)
    BFTC_FREE(export_connect_buffer);

  if (export_list != NULL)
    BFTC_FREE(export_list);

  /* Close MED file (to force its update) */

  if (rank == 0) {
    if (MEDfermer(writer->fid) != 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDfermer() failed to close file \"%s\"\n"),
                writer->filename);
    writer->fid = 0;
  }

  writer->is_open = false;
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a MED file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer      <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- variable definition location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvmc_to_med_export_field(void                            *this_writer,
                        const fvmc_nodal_t               *mesh,
                        const char                      *name,
                        const fvmc_writer_var_loc_t       location,
                        const int                        dimension,
                        const fvmc_interlace_t            interlace,
                        const int                        n_parent_lists,
                        const fvmc_lnum_t                 parent_num_shift[],
                        const fvmc_datatype_t             datatype,
                        const int                        time_step,
                        const double                     time_value,
                        const void                *const field_values[])
{
  int   i_char, n_chars, data_sizeof;
  int   med_mesh_num;
  int   output_dim;
  char  med_mesh_name[MED_TAILLE_NOM + 1];
  char  med_fieldname[MED_TAILLE_NOM + 1];

  fvmc_datatype_t  datatype_convert = FVMC_DATATYPE_NULL;
  med_type_champ  datatype_med = MED_FLOAT64;

  size_t   input_size = 0, output_size = 0, min_var_buffer_size = 0;
  size_t   max_grouped_elements_out = 0;
  size_t   alloc_size = 0;
  size_t   var_buffer_size = 0;

  unsigned char *var_buffer = NULL;
  fvmc_to_med_mesh_t *med_mesh = NULL;
  fvmc_writer_section_t  *export_list = NULL;
  fvmc_to_med_writer_t  *writer = (fvmc_to_med_writer_t *)this_writer;

  const fvmc_writer_section_t  *export_sections = NULL;
  fvmc_writer_field_helper_t  *helper = NULL;
  const int  n_ranks = writer->n_ranks;

  /* Re-open MED file */

  if (writer->is_open == false) {

    if (writer->rank == 0) {
      writer->fid = MEDouvrir(writer->filename, MED_LECTURE_ECRITURE);
      if (writer->fid < 0)
        bftc_error(__FILE__, __LINE__, 0,
                  _("MEDouvrir() failed to re-open file: %s\n"
                    "Associated writer: \"%s\"\n"),
                  writer->filename, writer->name);
    }
    writer->is_open = true;

  }

  /* Adapt dimension */

  output_dim = dimension;

  assert(output_dim > 0);

  if (dimension == 2)
    output_dim = 3;
  else if (dimension > 3 && dimension != 6 && dimension != 9)
    bftc_error(__FILE__, __LINE__, 0,
              _("Data of dimension %d not handled"), dimension);

  /* Get med_mesh_name */

  strncpy(med_mesh_name, mesh->name, MED_TAILLE_NOM);
  n_chars = strlen(med_mesh_name);

  for (i_char = n_chars + 1; i_char < MED_TAILLE_NOM; i_char++)
    med_mesh_name[i_char] = ' ';

  med_mesh_name[MED_TAILLE_NOM] = '\0';

  /* Get MED mesh structure */
  /*------------------------*/

  med_mesh_num = _get_med_mesh_num(writer,
                                   med_mesh_name);

  if (med_mesh_num == 0 )
    med_mesh_num = _add_med_mesh(writer,
                                 med_mesh_name,
                                 mesh);

  med_mesh = writer->med_meshes[med_mesh_num - 1];

  /* Adapt fvmc_datatype and find corresponding MED datatype */

  _get_datatypes(datatype,
                 &datatype_convert,
                 &datatype_med,
                 &data_sizeof);

  /* Create MED field if necessary. Always return MED fieldname */

  _get_med_fieldname(writer,
                     name,
                     datatype_med,
                     output_dim,
                     med_fieldname);

  /* Build list of sections that are used here, in order of output */
  /*---------------------------------------------------------------*/

  export_list = fvmc_writer_export_list(mesh,
                                       true,
                                       writer->discard_polygons,
                                       writer->discard_polyhedra,
                                       writer->divide_polygons,
                                       writer->divide_polyhedra);

  /* Build writer helper */
  /*---------------------*/

  helper = fvmc_writer_field_helper_create(mesh,
                                          export_list,
                                          output_dim,
                                          FVMC_INTERLACE,
                                          datatype_convert,
                                          location);

#if defined(FVMC_HAVE_MPI)

  fvmc_writer_field_helper_init_g(helper,
                                 export_list,
                                 mesh,
                                 writer->comm);

#endif

  /* Buffer size computation and allocation */
  /*----------------------------------------*/

  fvmc_writer_field_helper_get_size(helper,
                                   &input_size,
                                   &output_size,
                                   &max_grouped_elements_out,
                                   &min_var_buffer_size);

  /* Slicing allows for arbitrary buffer size, but should be small enough
     to add little additional memory requirement (in proportion), large
     enough to limit number of write and gather calls. No slicing is
     possible on rank 0, as MED does not provide partial writes */

  input_size *= output_dim;
  output_size *= output_dim;
  max_grouped_elements_out *= output_dim;

  var_buffer_size = input_size / n_ranks;

  var_buffer_size = FVMC_MAX(var_buffer_size, min_var_buffer_size);
  var_buffer_size = FVMC_MAX(var_buffer_size, 128*(size_t)output_dim);

  alloc_size = var_buffer_size;

  if (location == FVMC_WRITER_PER_NODE) {
    var_buffer_size = FVMC_MIN(var_buffer_size, output_size);
    if (writer->rank == 0)
      alloc_size = output_size;
  }
  else if (location == FVMC_WRITER_PER_ELEMENT) {
    var_buffer_size = FVMC_MIN(var_buffer_size, max_grouped_elements_out);
    if (writer->rank == 0)
      alloc_size = max_grouped_elements_out;
  }

  alloc_size *= fvmc_datatype_size[datatype_convert];

  BFTC_MALLOC(var_buffer, alloc_size, unsigned char);

  /* Export field */
  /*--------------*/

  if (location == FVMC_WRITER_PER_NODE) {

    _export_field_values_n(mesh,
                           helper,
                           writer,
                           dimension,
                           output_dim,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values,
                           med_mesh->name,
                           med_fieldname,
                           time_step,
                           time_value,
                           var_buffer_size,
                           var_buffer);

  }

  else if (location == FVMC_WRITER_PER_ELEMENT) {

    if (export_list == NULL)
      bftc_error(__FILE__, __LINE__, 0,
                _("MED must have entities.\n"
                  "Mesh: \"%s\"\n"
                  "Writer: \"%s\"\n"),
                med_mesh->name, writer->name);

    /* Loop on MED element types */

    export_sections = export_list;

    while (export_sections != NULL) {

      export_sections = _export_field_values_e(export_sections,
                                               helper,
                                               writer,
                                               dimension,
                                               output_dim,
                                               interlace,
                                               n_parent_lists,
                                               parent_num_shift,
                                               datatype,
                                               field_values,
                                               med_mesh->name,
                                               med_fieldname,
                                               time_step,
                                               time_value,
                                               var_buffer_size,
                                               var_buffer);

    } /* End of loop on MED element types */

  }
  else
    bftc_error(__FILE__, __LINE__, 0,
              "fvmc_to_med_export_field(): field location not managed.\n"
              "Associated writer: \"%s\"\n"
              "Associated med_mesh: \"%s\"\n"
              "Associated fieldname: \"%s\"\n"
              "Associated location: %i\n",
              writer->name, med_mesh_name, med_fieldname, location);

  /* Free buffers and helper structures */
  /*------------------------------------*/

  BFTC_FREE(var_buffer);

  helper = fvmc_writer_field_helper_destroy(helper);

  BFTC_FREE(export_list);

  /* Close MED file (to force its update) */
  /*--------------------------------------*/

  if (writer->rank == 0) {
    if (MEDfermer(writer->fid) != 0)
      bftc_error(__FILE__, __LINE__, 0,
                _("MEDfermer() failed to close file \"%s\"\n"),
                writer->filename);
    writer->fid = 0;
  }

  writer->is_open = false;
}

/*----------------------------------------------------------------------------*/



#ifdef __cplusplus
}
#endif /* __cplusplus */

#else
typedef struct {

  int              index;      /* CGNS base index */

} fvmc_to_med_fake_t;


#endif
