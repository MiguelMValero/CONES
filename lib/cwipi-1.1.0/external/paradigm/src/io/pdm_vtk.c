
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_vtk.h"
#include "pdm_error.h"

#include "pdm_logging.h"

#include "pdm_mesh_nodal.h"
#include "pdm_distrib.h"
#include "pdm_block_to_block.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_array.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

static char *
_fortran_to_c_string (
  const char *application_name_f,
  const int l_application_name_f
)
{
  char *application_name_c;
  int imin = 0;
  int imax = 0;

  while (imin < l_application_name_f && application_name_f[imin] == ' ') {
    imin++;
  }
  while (imax < l_application_name_f && application_name_f[l_application_name_f - imax - 1] == ' ') {
    imax++;
  }

  imax = l_application_name_f - imax - 1;

  assert(imax >= imin);

  if ((imax == l_application_name_f) || (imin == l_application_name_f)) {
    application_name_c =  (char *) malloc (sizeof(char) * (1));
    // application_name_c = (char *) malloc(sizeof(char) * 1);
    application_name_c[0] = '\0';
  }
  else {
    int size = imax - imin + 2;
    application_name_c =  (char *) malloc (sizeof(char) * (size));
    // application_name_c = (char *) malloc(sizeof(char) * size);
    int index = 0;
    for (int k = imin ; k <= imax ; k++) {
      application_name_c[index++] = application_name_f[k];
    }
    application_name_c[index] = '\0';
  }

  return application_name_c;
}

/**
 * \brief Export a polygonal mesh to ASCII VTK format (polydata)
 *
 * \param [in]  filename      Output file name
 * \param [in]  l_filename    Length of filename
 * \param [in]  n_vtx         Number of vertices
 * \param [in]  vtx_coord     Coordinates of the vertices (size = 3 * \ref n_vtx)
 *                            (x0, y0, z0, x1, ...)
 * \param [in]  vtx_g_num     Global ids of the vertices (or NULL)
 * \param [in]  n_face        Number of faces
 * \param [in]  face_vtx_idx  Index of the face-vertex connectivity (size = \ref n_face + 1)
 * \param [in]  face_vtx      Face-vertex connectivity (size = \ref face_vtx_idx[\ref n_face])
 * \param [in]  face_g_num    Global ids of the faces (or NULL)
 * \param [in]  face_color    Integer color of the faces (or NULL)
 *
 */

void
PDM_vtk_write_polydata_cf
(
 const char         *filename,
 const int           l_filename,
 const int           n_vtx,
 const double*       vtx_coord,
 const PDM_g_num_t*  vtx_g_num,
 const int           n_face,
 const int*          face_vtx_idx,
 const int*          face_vtx,
 const PDM_g_num_t * face_g_num,
 const int*          face_color
)
{
  char *c_filename = NULL;

  c_filename = _fortran_to_c_string(filename, l_filename);

  PDM_vtk_write_polydata(c_filename,
                         n_vtx,
                         vtx_coord,
                         vtx_g_num,
                         n_face,
                         face_vtx_idx,
                         face_vtx,
                         face_g_num,
                         face_color);

  free(c_filename);
}

/**
 * \brief Export a point cloud to ASCII VTK format (unstructured grid of points)
 *
 * \param [in]  filename      Output file name
 * \param [in]  l_filename    Length of filename
 * \param [in]  n_vtx         Number of points
 * \param [in]  vtx_coord     Coordinates of the points (size = 3 * \ref n_vtx)
 *                            (x0, y0, z0, x1, ...)
 * \param [in]  vtx_g_num     Global ids of the points (or NULL)
 * \param [in]  color         Integer color of the points (or NULL)
 *
 */

void
PDM_vtk_write_point_cloud_cf
(
 const char        *filename,
 const int           l_filename,
 const int          n_vtx,
 const double*      vtx_coord,
 const PDM_g_num_t* vtx_g_num,
 const int*         color
)
{
  char *c_filename = NULL;

  c_filename = _fortran_to_c_string(filename, l_filename);

  PDM_vtk_write_point_cloud(c_filename,
                            n_vtx,
                            vtx_coord,
                            vtx_g_num,
                            color);

  free(c_filename);
}

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static int _vtk_elt_type
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order
 )
{
  int vtk_elt_type = -1;

  switch (elt_type) {
    case PDM_MESH_NODAL_POINT:
      vtk_elt_type = 1;
      break;
    case PDM_MESH_NODAL_BAR2:
      vtk_elt_type = 3;
      break;
    case PDM_MESH_NODAL_TRIA3:
      vtk_elt_type = 5;
      break;
    case PDM_MESH_NODAL_QUAD4:
      vtk_elt_type = 9;
      break;
    case PDM_MESH_NODAL_TETRA4:
      vtk_elt_type = 10;
      break;
    case PDM_MESH_NODAL_PYRAMID5:
      vtk_elt_type = 14;
      break;
    case PDM_MESH_NODAL_PRISM6:
      vtk_elt_type = 13;
      break;
    case PDM_MESH_NODAL_HEXA8:
      vtk_elt_type = 12;
      break;

    case PDM_MESH_NODAL_BARHO:
      vtk_elt_type = 68;
      break;
    case PDM_MESH_NODAL_TRIAHO:
      vtk_elt_type = 69;
      break;
    case PDM_MESH_NODAL_QUADHO:
      vtk_elt_type = 70;
      break;
    case PDM_MESH_NODAL_TETRAHO:
      vtk_elt_type = 71;
      break;
    case PDM_MESH_NODAL_PYRAMIDHO:
      if (order == 2) {
        vtk_elt_type = 27;
      } else {
        vtk_elt_type = 66;//74;//??
      }
      break;
    case PDM_MESH_NODAL_PRISMHO:
      vtk_elt_type = 73;
      break;
    case PDM_MESH_NODAL_HEXAHO:
      vtk_elt_type = 72;
      break;

    case PDM_MESH_NODAL_BARHO_BEZIER:
      vtk_elt_type = 75;
      break;
    case PDM_MESH_NODAL_TRIAHO_BEZIER:
      vtk_elt_type = 76;
      break;

    default:
      PDM_error(__FILE__, __LINE__, 0, "type %d is not a valid std elt type\n", elt_type);
  }

  return vtk_elt_type;
}



static void
_ijk_to_vtk
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order,
 int                        idx[]
)
{
  switch (elt_type) {

  case PDM_MESH_NODAL_POINT:
  case PDM_MESH_NODAL_BAR2:
  case PDM_MESH_NODAL_TRIA3:
  case PDM_MESH_NODAL_QUAD4:
  case PDM_MESH_NODAL_TETRA4:
  case PDM_MESH_NODAL_PYRAMID5:
  case PDM_MESH_NODAL_PRISM6:
  case PDM_MESH_NODAL_HEXA8:
  {
    int n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    for (int i = 0; i < n_vtx; i++) {
      idx[i] = i;
    }
    break;
  }
  case PDM_MESH_NODAL_BARHO:
    idx[0] = 0;
    idx[1] = order;
    for (int i = 1; i < order; i++) {
      idx[i+1] = i;
    }
    break;

  case PDM_MESH_NODAL_TRIAHO:
    if (order == 1) {
      idx[2] = 2;
      idx[0] = 0; idx[1] = 1;
    } else if (order == 2) {
      idx[2] = 5;
      idx[5] = 3; idx[4] = 4;
      idx[0] = 0; idx[3] = 1; idx[1] = 2;
    } else if (order == 3) {
      idx[2] = 9;
      idx[7] = 7; idx[6] = 8;
      idx[8] = 4; idx[9] = 5; idx[5] = 6;
      idx[0] = 0; idx[3] = 1; idx[4] = 2; idx[1] = 3;
    } else if (order == 4) {
      idx[ 2] = 14;
      idx[ 9] = 12; idx[ 8] = 13;
      idx[10] =  9; idx[14] = 10; idx[ 7] = 11;
      idx[11] =  5; idx[12] =  6; idx[13] =  7; idx[ 6] =  8;
      idx[ 0] =  0; idx[ 3] =  1; idx[ 4] =  2; idx[ 5] =  3; idx[ 1] =  4;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "TRIA VTK ordering not implemented for order %d\n", order);
    }
    break;

  case PDM_MESH_NODAL_QUADHO:
    if (order == 1) {
      idx[3] = 2; idx[2] = 3;
      idx[0] = 0; idx[1] = 1;
    } else if (order == 2) {
      idx[3] = 6; idx[6] = 7; idx[2] = 8;
      idx[7] = 3; idx[8] = 4; idx[5] = 5;
      idx[0] = 0; idx[4] = 1; idx[1] = 2;
    } else if (order == 3) {
      idx[ 3] = 12; idx[ 8] = 13; idx[ 9] = 14; idx[ 2] = 15;
      idx[11] =  8; idx[14] =  9; idx[15] = 10; idx[ 7] = 11;
      idx[10] =  4; idx[12] =  5; idx[13] =  6; idx[ 6] =  7;
      idx[ 0] =  0; idx[ 4] =  1; idx[ 5] =  2; idx[ 1] =  3;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "QUAD VTK ordering not implemented for order %d\n", order);
    }
    break;

  case PDM_MESH_NODAL_TETRAHO:
    if (order == 1) {
      idx[3] = 3;

      idx[2] = 2;
      idx[0] = 0; idx[1] = 1;
    } else if (order == 2) {
      idx[3] = 9;

      idx[9] = 8;
      idx[7] = 6; idx[8] = 7;

      idx[2] = 5;
      idx[6] = 3; idx[5] = 4;
      idx[0] = 0; idx[4] = 1; idx[1] = 2;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "TETRA VTK ordering not implemented for order %d\n", order);
    }
    break;

  case PDM_MESH_NODAL_PRISMHO:
    if (order == 1) {
      idx[5] = 4;
      idx[3] = 3; idx[4] = 5;

      idx[2] = 1;
      idx[0] = 0; idx[1] = 2;
    } else if (order == 2) {
      idx[ 5] = 14;
      idx[11] = 13; idx[10] = 16;
      idx[ 3] = 12; idx[ 9] = 15; idx[ 4] = 17;

      idx[14] =  8;
      idx[17] =  7; idx[16] = 10;
      idx[12] =  6; idx[15] =  9; idx[13] = 11;

      idx[ 2] =  2;
      idx[ 8] =  1; idx[ 7] =  4;
      idx[ 0] =  0; idx[ 6] =  3; idx[ 1] =  5;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "PRISM VTK ordering not implemented for order %d\n", order);
    }
    break;

  case PDM_MESH_NODAL_HEXAHO:
    if (order == 1) {
      idx[7] = 6; idx[6] = 7;
      idx[4] = 4; idx[5] = 5;

      idx[3] = 2; idx[2] = 3;
      idx[0] = 0; idx[1] = 1;
    } else if (order == 2) {
      idx[ 7] = 24; idx[14] = 25; idx[ 6] = 26;
      idx[15] = 21; idx[21] = 22; idx[13] = 23;
      idx[ 4] = 18; idx[12] = 19; idx[ 5] = 20;

      idx[19] = 15; idx[24] = 16; idx[18] = 17;
      idx[25] = 12; idx[26] = 13; idx[23] = 14;
      idx[16] =  9; idx[22] = 10; idx[17] = 11;

      idx[ 3] =  6; idx[10] =  7; idx[ 2] =  8;
      idx[11] =  3; idx[20] =  4; idx[ 9] =  5;
      idx[ 0] =  0; idx[ 8] =  1; idx[ 1] =  2;
    } else {
      PDM_error(__FILE__, __LINE__, 0, "HEXA VTK ordering not implemented for order %d\n", order);
    }
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "VTK ordering not implemented for element type %d at order %d\n", (int) elt_type, order);
    break;
  }
}



static int *_vtk_lagrange_bar_to_ijk (const int order) {

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_BARHO, order);
  int *ijk = malloc (sizeof(int) * n_nodes);

  int idx = 0;
  ijk[idx++] = 0;

  ijk[idx++] = order;

  for (int i = 1; i < order; i++) {
    ijk[idx++] = i;
  }


  return ijk;
}


static int *_vtk_lagrange_tria_to_ijk (const int order) {

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_TRIAHO, order);
  int *ijk = malloc (sizeof(int) * n_nodes * 2);

  int idx = 0;

  // vertices
  ijk[idx++] = 0;
  ijk[idx++] = 0;

  ijk[idx++] = order;
  ijk[idx++] = 0;

  ijk[idx++] = 0;
  ijk[idx++] = order;

  // edges
  for (int i = 1; i < order; i++) {
    ijk[idx++] = i;
    ijk[idx++] = 0;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = order-i;
    ijk[idx++] = i;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = 0;
    ijk[idx++] = order-i;
  }

  // face
  if (order == 3) {
    ijk[idx++] = 1;
    ijk[idx++] = 1;
  }
  else if (order > 3) {
    int n_sub = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_TRIAHO, order-3);
    int *ijk_sub = _vtk_lagrange_tria_to_ijk (order-3);
    for (int i = 0; i < n_sub; i++) {
      ijk[idx++] = ijk_sub[2*i  ] + 1;
      ijk[idx++] = ijk_sub[2*i+1] + 1;
    }
    free (ijk_sub);
  }

  return ijk;
}


static int *_vtk_lagrange_quad_to_ijk (const int order) {

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_QUADHO, order);
  int *ijk = malloc (sizeof(int) * n_nodes * 2);

  int idx = 0;

  // vertices
  ijk[idx++] = 0;
  ijk[idx++] = 0;

  ijk[idx++] = order;
  ijk[idx++] = 0;

  ijk[idx++] = order;
  ijk[idx++] = order;

  ijk[idx++] = 0;
  ijk[idx++] = order;

  // edges
  for (int i = 1; i < order; i++) {
    ijk[idx++] = i;
    ijk[idx++] = 0;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = order;
    ijk[idx++] = i;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = i;
    ijk[idx++] = order;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = 0;
    ijk[idx++] = i;
  }

  // face
  for (int j = 1; j < order; j++) {
    for (int i = 1; i < order; i++) {
      ijk[idx++] = i;
      ijk[idx++] = j;
    }
  }

  return ijk;
}


static int *_vtk_lagrange_tetra_to_ijk (const int order) {

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_TETRAHO, order);
  int *ijk = malloc (sizeof(int) * n_nodes * 3);

  int idx = 0;

  // vertices
  ijk[idx++] = 0;
  ijk[idx++] = 0;
  ijk[idx++] = 0;

  ijk[idx++] = order;
  ijk[idx++] = 0;
  ijk[idx++] = 0;

  ijk[idx++] = 0;
  ijk[idx++] = order;
  ijk[idx++] = 0;

  ijk[idx++] = 0;
  ijk[idx++] = 0;
  ijk[idx++] = order;

  // edges
  for (int i = 1; i < order; i++) {
    ijk[idx++] = i;
    ijk[idx++] = 0;
    ijk[idx++] = 0;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = order-i;
    ijk[idx++] = i;
    ijk[idx++] = 0;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = 0;
    ijk[idx++] = order-i;
    ijk[idx++] = 0;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = 0;
    ijk[idx++] = 0;
    ijk[idx++] = i;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = order-i;
    ijk[idx++] = 0;
    ijk[idx++] = i;
  }

  for (int i = 1; i < order; i++) {
    ijk[idx++] = 0;
    ijk[idx++] = order-i;
    ijk[idx++] = i;
  }


  if (order == 3) {
    // face v=0
    ijk[idx++] = 1;
    ijk[idx++] = 0;
    ijk[idx++] = 1;

    // face u+v+w=order
    ijk[idx++] = 1;
    ijk[idx++] = 1;
    ijk[idx++] = 1;

    // face u=0
    ijk[idx++] = 0;
    ijk[idx++] = 1;
    ijk[idx++] = 1;

    // face w=0
    ijk[idx++] = 1;
    ijk[idx++] = 1;
    ijk[idx++] = 0;
  }
  else if (order > 3) {
    int n_sub = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_TRIAHO, order-3);
    int *ijk_sub = _vtk_lagrange_tria_to_ijk (order-3);

    // face v=0
    for (int i = 0; i < n_sub; i++) {
      ijk[idx++] = ijk_sub[2*i  ] + 1;
      ijk[idx++] = 0;
      ijk[idx++] = ijk_sub[2*i+1] + 1;
    }

    // face u+v+w=order
    for (int i = 0; i < n_sub; i++) {
      ijk[idx++] = ijk_sub[2*i+1] + 1;
      ijk[idx++] = order - (ijk_sub[2*i] + ijk_sub[2*i+1] + 2);
      ijk[idx++] = ijk_sub[2*i  ] + 1;
    }

    // face u=0
    for (int i = 0; i < n_sub; i++) {
      ijk[idx++] = 0;
      ijk[idx++] = ijk_sub[2*i+1] + 1;
      ijk[idx++] = ijk_sub[2*i  ] + 1;
    }

    // face w=0
    for (int i = 0; i < n_sub; i++) {
      ijk[idx++] = ijk_sub[2*i+1] + 1;
      ijk[idx++] = ijk_sub[2*i  ] + 1;
      ijk[idx++] = 0;
    }
    free (ijk_sub);

    // volume
    if (order == 4) {
      ijk[idx++] = 1;
      ijk[idx++] = 1;
      ijk[idx++] = 1;
    }
    else {
      n_sub = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_TETRAHO, order-4);
      ijk_sub = _vtk_lagrange_tetra_to_ijk (order-4);
      for (int i = 0; i < 3*n_sub; i++) {
        ijk[idx++] = ijk_sub[i] + 1;
      }
      free (ijk_sub);
    }
  }

  return ijk;
}


static int *_vtk_lagrange_prism_to_ijk (const int order) {

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_PRISMHO, order);
  int *ijk = malloc (sizeof(int) * n_nodes * 3);

  int idx = 0;

  // vertices
  for (int k = 0; k < 2; k++) {
    ijk[idx++] = 0;
    ijk[idx++] = 0;
    ijk[idx++] = order*k;

    ijk[idx++] = order;
    ijk[idx++] = 0;
    ijk[idx++] = order*k;

    ijk[idx++] = 0;
    ijk[idx++] = order;
    ijk[idx++] = order*k;
  }

  // edges
  for (int k = 0; k < 2; k++) {
    for (int i = 1; i < order; i++) {
      ijk[idx++] = i;
      ijk[idx++] = 0;
      ijk[idx++] = order*k;
    }

    for (int i = 1; i < order; i++) {
      ijk[idx++] = order-i;
      ijk[idx++] = i;
      // ijk[idx++] = i;
      // ijk[idx++] = order-i;
      ijk[idx++] = order*k;
    }

    for (int i = 1; i < order; i++) {
      ijk[idx++] = 0;
      ijk[idx++] = order-i;
      ijk[idx++] = order*k;
    }
  }

  for (int k = 1; k < order; k++) {
    ijk[idx++] = 0;
    ijk[idx++] = 0;
    ijk[idx++] = k;
  }

  for (int k = 1; k < order; k++) {
    ijk[idx++] = order;
    ijk[idx++] = 0;
    ijk[idx++] = k;
  }

  for (int k = 1; k < order; k++) {
    ijk[idx++] = 0;
    ijk[idx++] = order;
    ijk[idx++] = k;
  }

  // triangular faces
  for (int k = 0; k < 2; k++) {
    for (int j = 1; j < order; j++) {
      for (int i = 1; i < order-j; i++) {
        ijk[idx++] = i;
        ijk[idx++] = j;
        ijk[idx++] = order*k;
      }
    }
  }

  // quadrilateral faces
  for (int k = 1; k < order; k++) {
    for (int i = 1; i < order; i++) {
      ijk[idx++] = i;
      ijk[idx++] = 0;
      ijk[idx++] = k;
    }
  }

  for (int k = 1; k < order; k++) {
    for (int i = 1; i < order; i++) {
      ijk[idx++] = order-i;
      ijk[idx++] = i;
      ijk[idx++] = k;
    }
  }

  for (int k = 1; k < order; k++) {
    for (int j = 1; j < order; j++) {
      ijk[idx++] = 0;
      // ijk[idx++] = order-j;
      ijk[idx++] = j;
      ijk[idx++] = k;
    }
  }

  // volume
  for (int k = 1; k < order; k++) {
    for (int j = 1; j < order; j++) {
      for (int i = 1; i < order-j; i++) {
        ijk[idx++] = i;
        ijk[idx++] = j;
        ijk[idx++] = k;
      }
    }
  }

  return ijk;
}


static int *_vtk_lagrange_hexa_to_ijk (const int order) {

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_HEXAHO, order);
  int *ijk = malloc (sizeof(int) * n_nodes * 3);

  int idx = 0;

  // vertices
  for (int k = 0; k < 2; k++) {
    ijk[idx++] = 0;
    ijk[idx++] = 0;
    ijk[idx++] = order*k;

    ijk[idx++] = order;
    ijk[idx++] = 0;
    ijk[idx++] = order*k;

    ijk[idx++] = order;
    ijk[idx++] = order;
    ijk[idx++] = order*k;

    ijk[idx++] = 0;
    ijk[idx++] = order;
    ijk[idx++] = order*k;
  }

  // horizontal edges
  for (int k = 0; k < 2; k++) {
    for (int i = 1; i < order; i++) {
      ijk[idx++] = i;
      ijk[idx++] = 0;
      ijk[idx++] = order*k;
    }

    for (int i = 1; i < order; i++) {
      ijk[idx++] = order;
      ijk[idx++] = i;
      ijk[idx++] = order*k;
    }

    for (int i = 1; i < order; i++) {
      ijk[idx++] = i;
      ijk[idx++] = order;
      ijk[idx++] = order*k;
    }

    for (int i = 1; i < order; i++) {
      ijk[idx++] = 0;
      ijk[idx++] = i;
      ijk[idx++] = order*k;
    }
  }

  // vertical edges
  for (int j = 0; j < 2; j++) {
    for (int i = 0; i < 2; i++) {
      for (int k = 1; k < order; k++) {
        ijk[idx++] = order*i;
        ijk[idx++] = order*j;
        ijk[idx++] = k;
      }
    }
  }

  // faces normal to X
  for (int i = 0; i < 2; i++) {
    for (int k = 1; k < order; k++) {
      for (int j = 1; j < order; j++) {
        ijk[idx++] = order*i;
        ijk[idx++] = j;
        ijk[idx++] = k;
      }
    }
  }

  // faces normal to Y
  for (int j = 0; j < 2; j++) {
    for (int k = 1; k < order; k++) {
      for (int i = 1; i < order; i++) {
        ijk[idx++] = i;
        ijk[idx++] = order*j;
        ijk[idx++] = k;
      }
    }
  }

  // faces normal to Z
  for (int k = 0; k < 2; k++) {
    for (int j = 1; j < order; j++) {
      for (int i = 1; i < order; i++) {
        ijk[idx++] = i;
        ijk[idx++] = j;
        ijk[idx++] = order*k;
      }
    }
  }

  // volume
  for (int k = 1; k < order; k++) {
    for (int j = 1; j < order; j++) {
      for (int i = 1; i < order; i++) {
        ijk[idx++] = i;
        ijk[idx++] = j;
        ijk[idx++] = k;
      }
    }
  }

  return ijk;
}




static PDM_Mesh_nodal_elt_t _vtk_to_pdm_elt_type
(
 const int vtk_elt_type
 )
{
  switch (vtk_elt_type) {
  case 1:
    return PDM_MESH_NODAL_POINT;
  case 3:
    return PDM_MESH_NODAL_BAR2;
  case 5:
    return PDM_MESH_NODAL_TRIA3;
  case 9:
    return PDM_MESH_NODAL_QUAD4;
  case 10:
    return PDM_MESH_NODAL_TETRA4;
  case 14:
    return PDM_MESH_NODAL_PYRAMID5;
  case 13:
    return PDM_MESH_NODAL_PRISM6;
  case 12:
    return PDM_MESH_NODAL_HEXA8;
    default:
      PDM_error(__FILE__, __LINE__, 0, "VTK type %d is not supported\n", vtk_elt_type);
  }

  return PDM_MESH_NODAL_N_ELEMENT_TYPES;
}


typedef struct _prepa_vtk_field_t {

  int           n_field;
  char        **field_name;
  PDM_data_t   *field_type;
  int          *field_stride;
  void        **field_value;

} _prepa_vtk_field_t;


typedef struct _prepa_vtk_t {

  PDM_g_num_t            n_elt;
  PDM_g_num_t            n_vtx;
  double                *vtx_coord;
  int                   *elt_vtx_idx;
  PDM_g_num_t           *elt_vtx;
  PDM_Mesh_nodal_elt_t  *elt_type;

  _prepa_vtk_field_t     field[2];

} _prepa_vtk_t;


static void
_vtk_read_points
 (
  FILE         *f,
  PDM_g_num_t  *n_vtx,
  double      **vtx_coord
  )
{
  char word[999];

  fscanf(f, PDM_FMT_G_NUM, n_vtx);
  fscanf(f, "%s", word); // coord type

  *vtx_coord = malloc(sizeof(double) * (*n_vtx) * 3);

  for (int i = 0; i < 3*(*n_vtx); i++) {
    fscanf(f, "%le", &(*vtx_coord)[i]);
  }
}

static void
_vtk_read_unstructured_grid
 (
  FILE         *f,
  _prepa_vtk_t *prepa
 )
{
  char word[999];
  while (1) {
    int stat = fscanf(f, "%s", word);

    if (stat == EOF) {
      // End of file
      break;
    }

    if (strstr(word, "POINTS") != NULL) {
      _vtk_read_points(f,
                       &prepa->n_vtx,
                       &prepa->vtx_coord);
    }

    if (strstr(word, "CELLS") != NULL) {
      PDM_g_num_t gn_elt    = 0;
      PDM_g_num_t s_elt_vtx = 0;
      fscanf(f, PDM_FMT_G_NUM" "PDM_FMT_G_NUM, &gn_elt, &s_elt_vtx);
      if (prepa->n_elt < 0) {
        prepa->n_elt = gn_elt;
      }
      else if (prepa->n_elt != gn_elt) {
        PDM_error(__FILE__, __LINE__, 0, "Incoherent number of cells (block CELLS)\n");
      }
      s_elt_vtx -= prepa->n_elt;

      prepa->elt_vtx_idx = malloc(sizeof(int) * (prepa->n_elt + 1));
      prepa->elt_vtx_idx[0] = 0;
      prepa->elt_vtx = malloc(sizeof(PDM_g_num_t) * s_elt_vtx);

      for (int i = 0; i < prepa->n_elt; i++) {
        prepa->elt_vtx_idx[i+1] = prepa->elt_vtx_idx[i];
        int elt_vtx_n = 0;
        fscanf(f, "%d", &elt_vtx_n);
        for (int j = 0; j < elt_vtx_n; j++) {
          PDM_g_num_t vtx_id = 0;
          fscanf(f, PDM_FMT_G_NUM, &vtx_id);
          prepa->elt_vtx[prepa->elt_vtx_idx[i+1]++] = vtx_id + 1;
        }
      }
    }

    if (strstr(word, "CELL_TYPES") != NULL) {
      PDM_g_num_t gn_elt    = 0;
      fscanf(f, PDM_FMT_G_NUM, &gn_elt);
      if (prepa->n_elt < 0) {
        prepa->n_elt = gn_elt;
      }
      else if (prepa->n_elt != gn_elt) {
        PDM_error(__FILE__, __LINE__, 0, "Incoherent number of cells (block CELL_TYPES)\n");
      }

      prepa->elt_type = malloc(sizeof(PDM_Mesh_nodal_elt_t) * prepa->n_elt);
      for (int i = 0; i < prepa->n_elt; i++) {
        int vtk_elt_type = -1;
        fscanf(f, "%d", &vtk_elt_type);
        prepa->elt_type[i] = _vtk_to_pdm_elt_type(vtk_elt_type);
      }
      break;
    }
  }
}



static void
_vtk_read_polydata
 (
  FILE         *f,
  _prepa_vtk_t *prepa
 )
{
  char word[999];
  while (1) {
    int stat = fscanf(f, "%s", word);

    if (stat == EOF) {
      // End of file
      break;
    }

    if (strstr(word, "POINTS") != NULL) {
      _vtk_read_points(f,
                       &prepa->n_vtx,
                       &prepa->vtx_coord);
    }

    if (strstr(word, "POLYGONS") != NULL) {
      PDM_g_num_t gn_elt    = 0;
      PDM_g_num_t s_elt_vtx = 0;
      fscanf(f, PDM_FMT_G_NUM" "PDM_FMT_G_NUM, &gn_elt, &s_elt_vtx);
      if (prepa->n_elt < 0) {
        prepa->n_elt = gn_elt;
      }
      else if (prepa->n_elt != gn_elt) {
        PDM_error(__FILE__, __LINE__, 0, "Incoherent number of cells (block CELLS)\n");
      }
      s_elt_vtx -= prepa->n_elt;

      prepa->elt_vtx_idx = malloc(sizeof(int) * (prepa->n_elt + 1));
      prepa->elt_vtx_idx[0] = 0;
      prepa->elt_vtx = malloc(sizeof(PDM_g_num_t) * s_elt_vtx);
      prepa->elt_type = malloc(sizeof(PDM_Mesh_nodal_elt_t) * prepa->n_elt);

      for (int i = 0; i < prepa->n_elt; i++) {
        prepa->elt_type[i] = PDM_MESH_NODAL_POLY_2D;
        prepa->elt_vtx_idx[i+1] = prepa->elt_vtx_idx[i];
        int elt_vtx_n = 0;
        fscanf(f, "%d", &elt_vtx_n);
        for (int j = 0; j < elt_vtx_n; j++) {
          PDM_g_num_t vtx_id = 0;
          fscanf(f, PDM_FMT_G_NUM, &vtx_id);
          prepa->elt_vtx[prepa->elt_vtx_idx[i+1]++] = vtx_id + 1;
        }
      }
      break;
    }

  }
}


typedef enum {
  POINT_DATA,
  CELL_DATA
} _field_location_t;


static PDM_data_t
_vtk_str_to_field_type
(
const char *s
 )
{
  if (strstr(s, "int")  != NULL ||
      strstr(s, "long") != NULL) {
    return PDM_INT;
  }
  else if (strstr(s, "double")  != NULL ||
           strstr(s, "float")   != NULL) {
    return PDM_DOUBLE;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Field data type '%s' is not supported\n");
  }

  return PDM_WRONG_DATA;
}


static size_t
_pdm_data_size
(
 const PDM_data_t data_type
 )
{
  switch (data_type) {
  case PDM_INT:
    return sizeof(int);
  case PDM_DOUBLE:
    return sizeof(double);
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid data type %d\n", (int) data_type);
  }

  return 0;
}

static void
_vtk_read_field_values
(
       FILE         *f,
 const PDM_g_num_t   n,
 const int           stride,
 const PDM_data_t    data_type,
       void        **values
 )
{
  size_t s_data = _pdm_data_size(data_type);
  *values = malloc(s_data * n * stride);

  if (data_type == PDM_INT) {
    int *_field_value = (int *) *values;
    for (PDM_g_num_t j = 0; j < n; j++) {
      for (int k = 0; k < stride; k++) {
        fscanf(f, "%d", &_field_value[stride*j + k]);
      }
    }
  }
  else if (data_type == PDM_DOUBLE) {
    double *_field_value = (double *) *values;
    for (PDM_g_num_t j = 0; j < n; j++) {
      for (int k = 0; k < stride; k++) {
        fscanf(f, "%le", &_field_value[stride*j + k]);
      }
    }
  }
}



static void
_vtk_read_fields
 (
  FILE              *f,
  _prepa_vtk_t      *prepa,
  _field_location_t  location
  )
{
  char word[999];

  int n_field = 0;

  char       **field_name   = NULL;
  PDM_data_t  *field_type   = NULL;
  int         *field_stride = NULL;
  void       **field_value  = NULL;

  PDM_g_num_t n = 0;
  fscanf(f, PDM_FMT_G_NUM, &n);

  while (1) {
    int stat = fscanf(f, "%s", word);

    if (stat == EOF) {
      // End of file
      break;
    }

    if (strstr(word, "FIELD") != NULL) {
      fscanf(f, "%s", word); // FieldData

      fscanf(f, "%d", &n_field);

      field_name   = malloc(sizeof(char     *) * n_field);
      field_type   = malloc(sizeof(PDM_data_t) * n_field);
      field_stride = malloc(sizeof(int       ) * n_field);
      field_value  = malloc(sizeof(void     *) * n_field);

      for (int i = 0; i < n_field; i++) {
        fscanf(f, "%s", word);

        field_name[i] = malloc(sizeof(char) * (strlen(word) + 1));
        strcpy(field_name[i], word);

        fscanf(f, "%d", &field_stride[i]);

        PDM_g_num_t n2 = 0;
        fscanf(f, PDM_FMT_G_NUM, &n2);

        assert(n2 == n);

        fscanf(f, "%s", word);
        field_type[i] = _vtk_str_to_field_type(word);

        _vtk_read_field_values(f,
                               n,
                               field_stride[i],
                               field_type  [i],
                               &field_value[i]);
      }

      break; // allow FIELD + SCALARS + VECTORS...?
    }

    else if (strstr(word, "SCALARS") != NULL ||
             strstr(word, "VECTORS") != NULL ||
             strstr(word, "TENSORS") != NULL) {

      n_field = 1;

      field_name   = malloc(sizeof(char     *) * n_field);
      field_type   = malloc(sizeof(PDM_data_t) * n_field);
      field_stride = malloc(sizeof(int       ) * n_field);
      field_value  = malloc(sizeof(void     *) * n_field);

      if (strstr(word, "SCALARS") != NULL) {
       field_stride[0] = 1;
      }
      else if (strstr(word, "VECTORS") != NULL) {
       field_stride[0] = 3;
      }
      else if (strstr(word, "TENSORS") != NULL) {
       field_stride[0] = 9;
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "Invalid field type '%s'\n", word);
      }

      fscanf(f, "%s", word);
      field_name[0] = malloc(sizeof(char) * (strlen(word) + 1));
      strcpy(field_name[0], word);

      fscanf(f, "%s", word);
      field_type[0] = _vtk_str_to_field_type(word);

      while (1) {
        stat = fscanf(f, "%s", word);

        if (stat == EOF) {
          // End of file
          break;
        }

        if (strstr(word, "LOOKUP_TABLE") != NULL) {
          break;
        }
      }
      fscanf(f, "%s", word);
      assert(strstr(word, "default") != NULL);



      _vtk_read_field_values(f,
                             n,
                             field_stride[0],
                             field_type  [0],
                             &field_value[0]);

      break;
    }
  }

  prepa->field[location].n_field      = n_field;
  prepa->field[location].field_name   = field_name;
  prepa->field[location].field_type   = field_type;
  prepa->field[location].field_stride = field_stride;
  prepa->field[location].field_value  = field_value;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Export a set of boxes to ASCII VTK format (unstructured grid of hexahedra)
 *
 * \param [in]  filename     Output file name
 * \param [in]  n_box        Number of boxes
 * \param [in]  box_extents  Extents of the boxes (size = 6 * \ref n_box)
 *                           (xmin0, ymin0, zmin0, xmax0, ymax0, zmax0, xmin1, ...)
 * \param [in]  box_g_num    Global ids of the boxes (or NULL)
 */

void
PDM_vtk_write_boxes
(
 const char        *filename,
 const int          n_box,
 const double      *box_extents,
 const PDM_g_num_t *box_g_num
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nboxes\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_box);
  for (int i = 0; i < n_box; i++) {
    const double *e = box_extents + 6*i;
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[0], e[1], e[2]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[3], e[1], e[2]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[3], e[4], e[2]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[0], e[4], e[2]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[0], e[1], e[5]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[3], e[1], e[5]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[3], e[4], e[5]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[0], e[4], e[5]);
  }

  fprintf(f, "CELLS %d %d\n", n_box, 9*n_box);
  for (int i = 0; i < n_box; i++) {
    int j = 8*i;
    fprintf(f, "8 %d %d %d %d %d %d %d %d\n", j, j+1, j+2, j+3, j+4, j+5, j+6, j+7);
  }

  fprintf(f, "CELL_TYPES %d\n", n_box);
  for (int i = 0; i < n_box; i++) {
    fprintf(f, "12\n");
  }

  if (box_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_box);
    fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_box; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", box_g_num[i]);
    }
  }

  fclose(f);
}

void
PDM_vtk_write_boxes_with_field
(
 const char        *filename,
 const int          n_box,
 const double      *box_extents,
 const PDM_g_num_t *box_g_num,
 const int          n_box_field,
 const char        *box_field_name[],
 const double      *box_field[]
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nboxes\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_box);
  for (int i = 0; i < n_box; i++) {
    const double *e = box_extents + 6*i;
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[0], e[1], e[2]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[3], e[1], e[2]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[3], e[4], e[2]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[0], e[4], e[2]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[0], e[1], e[5]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[3], e[1], e[5]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[3], e[4], e[5]);
    fprintf(f, "%20.16f %20.16f %20.16f\n", e[0], e[4], e[5]);
  }

  fprintf(f, "CELLS %d %d\n", n_box, 9*n_box);
  for (int i = 0; i < n_box; i++) {
    int j = 8*i;
    fprintf(f, "8 %d %d %d %d %d %d %d %d\n", j, j+1, j+2, j+3, j+4, j+5, j+6, j+7);
  }

  fprintf(f, "CELL_TYPES %d\n", n_box);
  for (int i = 0; i < n_box; i++) {
    fprintf(f, "12\n");
  }

  if (box_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_box);
    fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_box; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", box_g_num[i]);
    }
  }

  if (n_box_field > 0) {
    assert (box_field != NULL);

    if (box_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_box);
    }

    fprintf(f, "FIELD box_field %d\n", n_box_field);
    for (int i = 0; i < n_box_field; i++) {
      assert (box_field_name[i] != NULL);

      fprintf(f, "%s 1 %d double\n", box_field_name[i], n_box);
      for (int j = 0; j < n_box; j++) {
        fprintf(f, "%lf ", box_field[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}



/**
 * \brief Export a set of circles to ASCII VTK format (unstructured grid of line segments)
 *
 * \param [in]  filename     Output file name
 * \param [in]  n_circles    Number of circles
 * \param [in]  center       Centers of the circles (size = 3 * \ref n_circles)
 *                           (x0, y0, z0, x1, ...)
 * \param [in]  radius       Radii of the circles (size = \ref n_circles)
 * \param [in]  g_num        Global ids of the circles (or NULL)
 * \param [in]  color        Integer color of the circles (or NULL)
 * \param [in]  resolution   Number of segments on each circle
 *
 */

void
PDM_vtk_write_circles
(
 const char        *filename,
 const int          n_circles,
 const double      *center,
 const double      *radius,
 const PDM_g_num_t *g_num,
 const int         *color,
 const int          resolution
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "circles\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  double step = 2. * PDM_PI / (double) (resolution - 1);

  fprintf(f, "POINTS %d double\n", n_circles * resolution);
  for (int i = 0; i < n_circles; i++) {
    for (int j = 0; j < resolution; j++) {
      fprintf(f, "%20.16f %20.16f %20.16f\n",
              center[3*i]     + radius[i]*cos(j*step),
              center[3*i + 1] + radius[i]*sin(j*step),
              center[3*i + 2]);
    }
  }

  fprintf(f, "CELLS %d %d\n", n_circles, n_circles * (resolution + 1));
  for (int i = 0; i < n_circles; i++) {
    fprintf(f, "%d \n", resolution);
    for (int j = 0; j < resolution-1; j++) {
      fprintf(f, "%d ", resolution*i + j);
    }
    fprintf(f, "%d\n", resolution*i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_circles);
  for (int i = 0; i < n_circles; i++) {
    fprintf(f, "4\n");
  }

  if (g_num != NULL && color == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_circles);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_circles; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", g_num[i]);
     }
  } else if (color != NULL && g_num == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_circles);
    fprintf(f, "SCALARS color int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_circles; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  } else if (g_num != NULL && color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_circles);
    fprintf(f, "FIELD field 2\n");
    fprintf(f, "gnum 1 %d long\n", n_circles);
    for (int i = 0; i < n_circles; i++) {
      fprintf(f, PDM_FMT_G_NUM" ", g_num[i]);
    }
    fprintf(f, "\ncolor 1 %d int\n", n_circles);
    for (int i = 0; i < n_circles; i++) {
      fprintf(f, "%d ", color[i]);
    }
  }

  fclose(f);
}



/**
 * \brief Export a polygonal mesh to ASCII VTK format (polydata)
 *
 * \param [in]  filename      Output file name
 * \param [in]  n_vtx         Number of vertices
 * \param [in]  vtx_coord     Coordinates of the vertices (size = 3 * \ref n_vtx)
 *                            (x0, y0, z0, x1, ...)
 * \param [in]  vtx_g_num     Global ids of the vertices (or NULL)
 * \param [in]  n_face        Number of faces
 * \param [in]  face_vtx_idx  Index of the face-vertex connectivity (size = \ref n_face + 1)
 * \param [in]  face_vtx      Face-vertex connectivity (size = \ref face_vtx_idx[\ref n_face])
 * \param [in]  face_g_num    Global ids of the faces (or NULL)
 * \param [in]  face_color    Integer color of the faces (or NULL)
 *
 */

void
PDM_vtk_write_polydata
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const int          face_vtx[],
 const PDM_g_num_t  face_g_num[],
 const int          face_color[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int s_face_vtx = 0;
  if (n_face > 0) {
    assert(face_vtx_idx != NULL);
    s_face_vtx = face_vtx_idx[n_face];
  }
  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + s_face_vtx);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, " %d", face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (face_g_num != NULL && face_color == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[i]);
     }
  } else if (face_color != NULL && face_g_num == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_color int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%d\n", face_color[i]);
    }
  } else if (face_g_num != NULL && face_color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "FIELD field 2\n");
    fprintf(f, "face_gnum 1 %d long\n", n_face);
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM" ", face_g_num[i]);
    }
    fprintf(f, "\nface_color 1 %d int\n", n_face);
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%d ", face_color[i]);
    }
  }




  fclose(f);
}

void
PDM_vtk_write_polydata_with_field
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const int          face_vtx[],
 const PDM_g_num_t  face_g_num[],
 const int          face_color[],
 const int          n_elt_ifield,
 const char        *elt_ifield_name[],
 const int         *elt_ifield[]
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int s_face_vtx = 0;
  if (n_face > 0) {
    assert(face_vtx_idx != NULL);
    s_face_vtx = face_vtx_idx[n_face];
  }
  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + s_face_vtx);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, " %d", face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (face_g_num != NULL && face_color == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[i]);
     }
  } else if (face_color != NULL && face_g_num == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_color int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%d\n", face_color[i]);
    }
  } else if (face_g_num != NULL && face_color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "FIELD field 2\n");
    fprintf(f, "face_gnum 1 %d long\n", n_face);
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM" ", face_g_num[i]);
    }
    fprintf(f, "\nface_color 1 %d int\n", n_face);
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%d ", face_color[i]);
    }
  }

  if (n_elt_ifield > 0) {
    assert (elt_ifield != NULL);

    if (face_g_num == NULL && face_color == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_face);
    }

    fprintf(f, "FIELD elt_field %d\n", n_elt_ifield);
    for (int i = 0; i < n_elt_ifield; i++) {
      // assert (elt_ifield[i] != NULL);
      assert (elt_ifield_name[i] != NULL);

      fprintf(f, "%s 1 %d int\n", elt_ifield_name[i], n_face);
      for (int j = 0; j < n_face; j++) {
        fprintf(f, "%d ", elt_ifield[i][j]);
      }
      fprintf(f, "\n");
    }

  }

  fclose(f);
}



void
PDM_vtk_write_polydata_field
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const int          face_vtx[],
 const PDM_g_num_t  face_g_num[],
 const char         face_field_name[],
 const double       face_field[],
 const char         vtx_field_name[],
 const double       vtx_field[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int s_face_vtx = 0;
  if (n_face > 0) {
    assert(face_vtx_idx != NULL);
    s_face_vtx = face_vtx_idx[n_face];
  }
  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + s_face_vtx);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, " %d", face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }


  if (vtx_g_num != NULL && vtx_field == NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
     }
  } else if (vtx_field != NULL && vtx_g_num == NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS %s double 1\n", vtx_field_name);
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, "%20.16f\n", vtx_field[i]);
    }
  } else if (vtx_g_num != NULL && vtx_field != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "FIELD field 2\n");
    fprintf(f, "vtx_gnum 1 %d long\n", n_vtx);
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM" ", vtx_g_num[i]);
    }
    fprintf(f, "\n%s 1 %d double\n", vtx_field_name, n_vtx);
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, "%20.16f ", vtx_field[i]);
    }
  }

  if (face_g_num != NULL && face_field == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[i]);
     }
  } else if (face_field != NULL && face_g_num == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS %s double 1\n", face_field_name);
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%20.16f\n", face_field[i]);
    }
  } else if (face_g_num != NULL && face_field != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "FIELD field 2\n");
    fprintf(f, "face_gnum 1 %d long\n", n_face);
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM" ", face_g_num[i]);
    }
    fprintf(f, "\n%s 1 %d double\n", face_field_name, n_face);
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%20.16f ", face_field[i]);
    }
  }




  fclose(f);
}


/**
 * \brief Export a point cloud to ASCII VTK format (unstructured grid of points)
 *
 * \param [in]  filename      Output file name
 * \param [in]  n_vtx         Number of points
 * \param [in]  vtx_coord     Coordinates of the points (size = 3 * \ref n_vtx)
 *                            (x0, y0, z0, x1, ...)
 * \param [in]  vtx_g_num     Global ids of the points (or NULL)
 * \param [in]  color         Integer color of the points (or NULL)
 *
 */

void
PDM_vtk_write_point_cloud
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          color[]
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "point cloud\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_vtx, 2*n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    fprintf(f, "1\n");
  }

  if (vtx_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  } else if (color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS color int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  }

  fclose(f);
}



/**
 * \brief Export a set of lines to ASCII VTK format (unstructured grid of line segments)
 *
 * \param [in]  filename      Output file name
 * \param [in]  n_line        Number of lines
 * \param [in]  coord         Coordinates of the vertices (size = 6 * \ref n_line)
 *                            (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [in]  g_num         Global ids of the lines (or NULL)
 * \param [in]  color         Integer color of the lines (or NULL)
 *
 */

void
PDM_vtk_write_lines
(
 const char        *filename,
 const int          n_line,
 const double      *coord,
 const PDM_g_num_t *g_num,
 const int         *color
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "lines\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_line * 2);
  for (int i = 0; i < 2*n_line; i++) {
    fprintf(f, "%20.16f %20.16f %20.16f\n",
            coord[3*i    ],
            coord[3*i + 1],
            coord[3*i + 2]);
  }

  fprintf(f, "CELLS %d %d\n", n_line, n_line * 3);
  for (int i = 0; i < n_line; i++) {
    fprintf(f, "%d %d %d\n", 2, 2*i, 2*i+1);
  }

  fprintf(f, "CELL_TYPES %d\n", n_line);
  for (int i = 0; i < n_line; i++) {
    fprintf(f, "3\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_line);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_line; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  if (color != NULL) {
    if (g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_line);
    }
    fprintf(f, "FIELD line_field 1\n");
    fprintf(f, "color 1 %d int\n", n_line);
    for (int i = 0; i < n_line; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  }

  fclose(f);
}



/**
 * \brief Export a block of standard elements to ASCII VTK format (unstructured grid)
 *
 * Export a block of elements with multiple cell-based, integer-valued fields.
 *
 * \param [in]  filename        Output file name
 * \param [in]  n_vtx           Number of vertices
 * \param [in]  vtx_coord       Coordinates of the vertices (size = 3 * \ref n_vtx)
 *                              (x0, y0, z0, x1, ...)
 * \param [in]  vtx_g_num       Global ids of the vertices (or NULL)
 * \param [in]  elt_type        Type of elements
 * \param [in]  n_elt           Number of elements
 * \param [in]  elt_vtx         Element-vertex connectivity (size = \ref n_elt * n_vtx_per_elt)
 * \param [in]  elt_g_num       Global ids of the elements (or NULL)
 * \param [in]  n_elt_ifield    Number of fields
 * \param [in]  elt_ifield_name Name of the fields (or NULL)
 * \param [in]  elt_ifield      Fields (or NULL)
 *
 */

void
PDM_vtk_write_std_elements
(
 const char                 *filename,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_ifield,
 const char                 *elt_ifield_name[],
 const int                  *elt_ifield[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get (elt_type, 1);

  fprintf(f, "CELLS %d %d\n", n_elt, n_elt * (1 + n_vtx_elt));
  if (elt_type == PDM_MESH_NODAL_POINT) {
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, "1 %d\n", i);
    }
  }
  else {
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, "%d", n_vtx_elt);
      for (int j = 0; j < n_vtx_elt; j++) {
        fprintf(f, " %d", elt_vtx[n_vtx_elt*i + j] - 1);
      }
      fprintf(f, "\n");
    }
  }

  int vtk_elt_type = _vtk_elt_type (elt_type, 1);
  fprintf(f, "CELL_TYPES %d\n", n_elt);
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d\n", vtk_elt_type);
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (elt_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_elt);
    fprintf(f, "SCALARS elt_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", elt_g_num[i]);
    }
  }

  if (n_elt_ifield > 0) {
    assert (elt_ifield != NULL);

    if (elt_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_elt);
    }

    fprintf(f, "FIELD elt_field %d\n", n_elt_ifield);
    for (int i = 0; i < n_elt_ifield; i++) {
      // assert (elt_ifield[i] != NULL);
      assert (elt_ifield_name[i] != NULL);

      fprintf(f, "%s 1 %d int\n", elt_ifield_name[i], n_elt);
      for (int j = 0; j < n_elt; j++) {
        fprintf(f, "%d ", elt_ifield[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}



/**
 * \brief Export a block of standard elements to ASCII VTK format (unstructured grid)
 *
 * Export a block of elements with multiple cell-based, real-valued fields.
 *
 * \param [in]  filename        Output file name
 * \param [in]  n_vtx           Number of vertices
 * \param [in]  vtx_coord       Coordinates of the vertices (size = 3 * \ref n_vtx)
 *                              (x0, y0, z0, x1, ...)
 * \param [in]  vtx_g_num       Global ids of the vertices (or NULL)
 * \param [in]  elt_type        Type of elements
 * \param [in]  n_elt           Number of elements
 * \param [in]  elt_vtx         Element-vertex connectivity (size = \ref n_elt * n_vtx_per_elt)
 * \param [in]  elt_g_num       Global ids of the elements (or NULL)
 * \param [in]  n_elt_field     Number of fields
 * \param [in]  elt_field_name  Name of the fields (or NULL)
 * \param [in]  elt_field       Fields (or NULL)
 *
 */

void
PDM_vtk_write_std_elements_double
(
 const char                 *filename,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_field,
 const char                 *elt_field_name[],
 const double               *elt_field[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get (elt_type, 1);

  fprintf(f, "CELLS %d %d\n", n_elt, n_elt * (1 + n_vtx_elt));
  if (elt_type == PDM_MESH_NODAL_POINT) {
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, "1 %d\n", i);
    }
  }
  else {
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, "%d", n_vtx_elt);
      for (int j = 0; j < n_vtx_elt; j++) {
        fprintf(f, " %d", elt_vtx[n_vtx_elt*i + j] - 1);
      }
      fprintf(f, "\n");
    }
  }


  int vtk_elt_type = _vtk_elt_type (elt_type, 1);
  fprintf(f, "CELL_TYPES %d\n", n_elt);
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d\n", vtk_elt_type);
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (elt_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_elt);
    fprintf(f, "SCALARS elt_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", elt_g_num[i]);
     }
  }

  if (n_elt_field > 0) {
    assert (elt_field != NULL);

    if (elt_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_elt);
    }

    fprintf(f, "FIELD elt_field %d\n", n_elt_field);
    for (int i = 0; i < n_elt_field; i++) {
      // assert (elt_field[i] != NULL);
      assert (elt_field_name[i] != NULL);

      fprintf(f, "%s 1 %d double\n", elt_field_name[i], n_elt);
      for (int j = 0; j < n_elt; j++) {
        fprintf(f, "%lf ", elt_field[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}

/**
 * \brief Export a point cloud to ASCII VTK format (unstructured grid of points)
 *
 * \param [in]  filename              Output file name
 * \param [in]  n_vtx                 Number of points
 * \param [in]  vtx_coord             Coordinates of the points (size = 3 * \ref n_vtx)
 *                                    (x0, y0, z0, x1, ...)
 * \param [in]  vtx_g_num             Global ids of the points (or NULL)
 * \param [in]  color                 Integer color of the points (or NULL)
 * \param [in]  n_vtx_field           Number of vertex fields
 * \param [in]  vtx_field_name        Name of those vertex fields
 * \param [in]  vtx_field             Vertex fields
 * \param [in]  n_vtx_vector_field    Number of vertex vector fields
 * \param [in]  vtx_vector_field_name Name of those vertex vector fields
 * \param [in]  vtx_vector_field      Vertex vector fields
 * \param [in]  n_vtx_normal_field    Number of vertex normal fields
 * \param [in]  vtx_normal_field_name Name of those vertex normal fields
 * \param [in]  vtx_normal_field      Vertex normal fields
 *
 */

void
PDM_vtk_write_point_cloud_with_field
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          color[],
 const int          n_vtx_field,
 const char        *vtx_field_name[],
 const double      *vtx_field[],
 const int          n_vtx_vector_field,
 const char        *vtx_vector_field_name[],
 const double      *vtx_vector_field[],
 const int          n_vtx_normal_field,
 const char        *vtx_normal_field_name[],
 const double      *vtx_normal_field[]
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "point cloud\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_vtx, 2*n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    fprintf(f, "1\n");
  }

  if (n_vtx_vector_field > 0 || n_vtx_normal_field > 0) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
  }

  if (n_vtx_vector_field > 0) {
    assert (vtx_vector_field != NULL);
    for (int h = 0; h < n_vtx_vector_field; h++) {
      assert (vtx_vector_field_name[h] != NULL);
      fprintf(f, "VECTORS %s double\n", vtx_vector_field_name[h]);
      for (int i = 0; i < n_vtx; i++) {
        for (int j = 0; j < 3; j++) {
          fprintf(f, "%.20lf ", vtx_vector_field[h][3*i+j]);
        }
        fprintf(f, "\n");
      }
    }
  }

  if (n_vtx_normal_field > 0) {
    assert (vtx_normal_field != NULL);
    for (int h = 0; h < n_vtx_normal_field; h++) {
      assert (vtx_normal_field_name[h] != NULL);
      fprintf(f, "NORMALS %s double\n", vtx_normal_field_name[h]);
      for (int i = 0; i < n_vtx; i++) {
        for (int j = 0; j < 3; j++) {
          fprintf(f, "%.20lf ", vtx_normal_field[h][3*i+j]);
        }
        fprintf(f, "\n");
      }
    }
  }

  if (vtx_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  } else if (color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS color int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  }

  if (n_vtx_field > 0) {
    assert (vtx_field != NULL);

    if (vtx_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_vtx);
    }

    fprintf(f, "FIELD vtx_field %d\n", n_vtx_field);
    for (int i = 0; i < n_vtx_field; i++) {
      // assert (vtx_field[i] != NULL);
      assert (vtx_field_name[i] != NULL);

      fprintf(f, "%s 1 %d double\n", vtx_field_name[i], n_vtx);
      for (int j = 0; j < n_vtx; j++) {
        fprintf(f, "%lf ", vtx_field[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}


/**
 * \brief Export a block of elements of arbitray order to ASCII VTK format (unstructured grid)
 *
 * Export a block of elements of arbitray order with multiple cell-based, real-valued fields.
 *
 * \param [in]  filename        Output file name
 * \param [in]  order           Geometric order of the elements
 * \param [in]  n_vtx           Number of vertices
 * \param [in]  vtx_coord       Coordinates of the vertices (size = 3 * \ref n_vtx)
 *                              (x0, y0, z0, x1, ...)
 * \param [in]  vtx_g_num       Global ids of the vertices (or NULL)
 * \param [in]  elt_type        Type of elements
 * \param [in]  n_elt           Number of elements
 * \param [in]  elt_vtx         Element-vertex connectivity (size = \ref n_elt * n_vtx_per_elt)
 * \param [in]  elt_g_num       Global ids of the elements (or NULL)
 * \param [in]  n_elt_field     Number of fields
 * \param [in]  elt_field_name  Name of the fields (or NULL)
 * \param [in]  elt_field       Fields (or NULL)
 *
 */

void
PDM_vtk_write_std_elements_ho
(
 const char                 *filename,
 const int                   order,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_field,
 const char                 *elt_field_name[],
 const double               *elt_field[]
 )
{
  //assert (order < 3);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get (elt_type, order);
  if(0 == 1) {
    int *vtk_idx = malloc (sizeof(int) * n_vtx_elt);
    _ijk_to_vtk (elt_type, order, vtk_idx);
  }

  fprintf(f, "CELLS %d %d\n", n_elt, n_elt * (1 + n_vtx_elt));
  if (elt_type == PDM_MESH_NODAL_POINT) {
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, "1 %d\n", i);
    }
  }
  else {
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, "%d", n_vtx_elt);
      for (int j = 0; j < n_vtx_elt; j++) {
      //fprintf(f, " %d", elt_vtx[n_vtx_elt*i + vtk_idx[j]] - 1);
        fprintf(f, " %d", elt_vtx[n_vtx_elt*i + j] - 1);
      }
      fprintf(f, "\n");
    }
  }


  int vtk_elt_type = _vtk_elt_type (elt_type, order);
  fprintf(f, "CELL_TYPES %d\n", n_elt);
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d\n", vtk_elt_type);
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (elt_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_elt);
    fprintf(f, "SCALARS elt_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", elt_g_num[i]);
     }
  }

  if (n_elt_field > 0) {
    assert (elt_field != NULL);

    if (elt_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_elt);
    }

    fprintf(f, "FIELD elt_field %d\n", n_elt_field);
    for (int i = 0; i < n_elt_field; i++) {
      // assert (elt_field[i] != NULL);
      assert (elt_field_name[i] != NULL);

      fprintf(f, "%s 1 %d double\n", elt_field_name[i], n_elt);
      for (int j = 0; j < n_elt; j++) {
        fprintf(f, "%lf ", elt_field[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}



void
PDM_vtk_write_std_elements_ho_with_vtx_field
(
 const char                 *filename,
 const int                   order,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_field,
 const char                 *elt_field_name[],
 const double               *elt_field[],
 const int                   n_vtx_field,
 const char                 *vtx_field_name[],
 const double               *vtx_field[]
 )
{
  //assert (order < 3);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get (elt_type, order);
  if(0 == 1) {
    int *vtk_idx = malloc (sizeof(int) * n_vtx_elt);
    _ijk_to_vtk (elt_type, order, vtk_idx);
  }

  fprintf(f, "CELLS %d %d\n", n_elt, n_elt * (1 + n_vtx_elt));
  if (elt_type == PDM_MESH_NODAL_POINT) {
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, "1 %d\n", i);
    }
  }
  else {
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, "%d", n_vtx_elt);
      for (int j = 0; j < n_vtx_elt; j++) {
      //fprintf(f, " %d", elt_vtx[n_vtx_elt*i + vtk_idx[j]] - 1);
        fprintf(f, " %d", elt_vtx[n_vtx_elt*i + j] - 1);
      }
      fprintf(f, "\n");
    }
  }


  int vtk_elt_type = _vtk_elt_type (elt_type, order);
  fprintf(f, "CELL_TYPES %d\n", n_elt);
  for (int i = 0; i < n_elt; i++) {
    fprintf(f, "%d\n", vtk_elt_type);
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (n_vtx_field > 0) {
    assert (vtx_field != NULL);

    if (vtx_g_num == NULL) {
      fprintf(f, "POINT_DATA %d\n", n_vtx);
    }

    fprintf(f, "FIELD vtx_field %d\n", n_vtx_field);
    for (int i = 0; i < n_vtx_field; i++) {
      // assert (vtx_field[i] != NULL);
      assert (vtx_field_name[i] != NULL);

      fprintf(f, "%s 1 %d double\n", vtx_field_name[i], n_vtx);
      for (int j = 0; j < n_vtx; j++) {
        fprintf(f, "%lf ", vtx_field[i][j]);
      }
      fprintf(f, "\n");
    }
  }


  if (elt_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_elt);
    fprintf(f, "SCALARS elt_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_elt; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", elt_g_num[i]);
     }
  }

  if (n_elt_field > 0) {
    assert (elt_field != NULL);

    if (elt_g_num == NULL) {
      fprintf(f, "CELL_DATA %d\n", n_elt);
    }

    fprintf(f, "FIELD elt_field %d\n", n_elt_field);
    for (int i = 0; i < n_elt_field; i++) {
      // assert (elt_field[i] != NULL);
      assert (elt_field_name[i] != NULL);

      fprintf(f, "%s 1 %d double\n", elt_field_name[i], n_elt);
      for (int j = 0; j < n_elt; j++) {
        fprintf(f, "%lf ", elt_field[i][j]);
      }
      fprintf(f, "\n");
    }
  }

  fclose(f);
}



/**
 * \brief Export a set of ellipses to ASCII VTK format (unstructured grid of line segments)
 *
 * \param [in]  filename     Output file name
 * \param [in]  n_ellipse    Number of ellipses
 * \param [in]  center       Centers of the ellipses (size = 3 * \ref n_circles)
 *                           (x0, y0, z0, x1, ...)
 * \param [in]  axes         Axes of the ellipses (size = 6 * \ref n_circles)
 *                           (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [in]  radii       Radii of the ellipses (size = 2 * \ref n_ellipse)
 * \param [in]  g_num        Global ids of the ellipses (or NULL)
 * \param [in]  color        Integer color of the ellipses (or NULL)
 * \param [in]  resolution   Number of segments on each ellipse
 *
 */

void
PDM_vtk_write_ellipses
(
 const char        *filename,
 const int          n_ellipse,
 const double      *center,
 const double      *axes,
 const double      *radii,
 const PDM_g_num_t *g_num,
 const int         *color,
 const int          resolution
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "ellipses\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  double step = 2. * PDM_PI / (double) (resolution - 1);

  fprintf(f, "POINTS %d double\n", n_ellipse * resolution);
  for (int i = 0; i < n_ellipse; i++) {
    for (int j = 0; j < resolution; j++) {
      double x = radii[2*i  ] * cos(j*step);
      double y = radii[2*i+1] * sin(j*step);
      fprintf(f, "%20.16f %20.16f %20.16f\n",
              center[3*i    ] + x*axes[6*i  ] + y*axes[6*i+3],
              center[3*i + 1] + x*axes[6*i+1] + y*axes[6*i+4],
              center[3*i + 2] + x*axes[6*i+2] + y*axes[6*i+5]);
    }
  }

  fprintf(f, "CELLS %d %d\n", n_ellipse, n_ellipse * (resolution + 1));
  for (int i = 0; i < n_ellipse; i++) {
    fprintf(f, "%d \n", resolution);
    for (int j = 0; j < resolution-1; j++) {
      fprintf(f, "%d ", resolution*i + j);
    }
    fprintf(f, "%d\n", resolution*i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_ellipse);
  for (int i = 0; i < n_ellipse; i++) {
    fprintf(f, "4\n");
  }

  if (g_num != NULL && color == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_ellipse);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_ellipse; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", g_num[i]);
     }
  } else if (color != NULL && g_num == NULL) {
    fprintf(f, "CELL_DATA %d\n", n_ellipse);
    fprintf(f, "SCALARS color int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_ellipse; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  } else if (g_num != NULL && color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_ellipse);
    fprintf(f, "FIELD field 2\n");
    fprintf(f, "gnum 1 %d long\n", n_ellipse);
    for (int i = 0; i < n_ellipse; i++) {
      fprintf(f, PDM_FMT_G_NUM" ", g_num[i]);
    }
    fprintf(f, "\ncolor 1 %d int\n", n_ellipse);
    for (int i = 0; i < n_ellipse; i++) {
      fprintf(f, "%d ", color[i]);
    }
  }

  fclose(f);
}



/**
 * \brief Get the ijk-coordinates of nodes in a VTK Lagrange high-order element
 *
 * \param [in]  elt_type  Type of element
 * \param [in]  order     Geometric order of the element
 *
 * \return                Array of ijk-coordinates of the nodes
 *                        (size = n_nodes * dim_elt)
 */

int *
PDM_vtk_lagrange_to_ijk
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order
 )
{
  switch (elt_type) {
  case PDM_MESH_NODAL_POINT: {
    int *_ijk = malloc (sizeof(int));
    _ijk[0] = 0;
    return _ijk;
  }

  case PDM_MESH_NODAL_BARHO:
    return _vtk_lagrange_bar_to_ijk(order);

  case PDM_MESH_NODAL_TRIAHO:
    return _vtk_lagrange_tria_to_ijk(order);

  case PDM_MESH_NODAL_QUADHO:
    return _vtk_lagrange_quad_to_ijk(order);

  case PDM_MESH_NODAL_TETRAHO:
    return _vtk_lagrange_tetra_to_ijk(order);

  case PDM_MESH_NODAL_PRISMHO:
    return _vtk_lagrange_prism_to_ijk(order);

  case PDM_MESH_NODAL_HEXAHO:
    return _vtk_lagrange_hexa_to_ijk(order);

  default:
    PDM_error(__FILE__, __LINE__, 0, "VTK lagrange ordering not implemented for elt type %d\n", (int) elt_type);
  }

  return NULL;
}


/**
 *
 * \brief Create a dmesh nodal from a file in ASCII VTK mesh format
 *
 * \param[in]  comm                MPI communicator
 * \param[in]  filename            Filename
 *
 * \return Pointer to PDM_dmesh_nodal object
 *
 */


PDM_dmesh_nodal_t *
PDM_vtk_read_to_dmesh_nodal
(
 const PDM_MPI_Comm    comm,
 const char           *filename,
       int            *n_vtx_field,
       char         ***vtx_field_name,
       PDM_data_t    **vtx_field_type,
       int           **vtx_field_stride,
       void         ***vtx_field_value,
       int            *n_elt_field,
       char         ***elt_field_name,
       PDM_data_t    **elt_field_type,
       int           **elt_field_stride,
       void         ***elt_field_value
)
{
  PDM_UNUSED(elt_field_value);
  int dbg_enabled = 0;

  PDM_UNUSED(elt_field_value);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  _prepa_vtk_t prepa;
  prepa.n_elt = -1;
  prepa.n_vtx =  0;
  prepa.vtx_coord   = NULL;
  prepa.elt_vtx_idx = NULL;
  prepa.elt_vtx     = NULL;
  prepa.elt_type    = NULL;

  for (int i = 0; i < 2; i++) {
    prepa.field[i].n_field = 0;
    prepa.field[i].field_name   = NULL;
    prepa.field[i].field_type   = NULL;
    prepa.field[i].field_stride = NULL;
    prepa.field[i].field_value  = NULL;
  }

  *n_vtx_field = 0;
  *n_elt_field = 0;

  if (i_rank == 0) {

    FILE *f = fopen(filename, "r");

    if (f == NULL) {
      PDM_error(__FILE__, __LINE__, 0, "Failed to open file '%s'\n", filename);
    }

    char word[999];
    while (1) {

      int stat = fscanf(f, "%s", word);

      if (stat == EOF) {
        // End of file
        break;
      }


      if (strstr(word, "DATASET") != NULL) {
        // Get dataset type
        stat = fscanf(f, "%s", word);

        // UNSTRUCTURED_GRID
        // POLYDATA
        // STRUCTURED_GRID --> unstructured hexahedra (TO DO)

        if (strstr(word, "UNSTRUCTURED_GRID") != NULL) {
          _vtk_read_unstructured_grid(f, &prepa);
        }
        else if (strstr(word, "POLYDATA") != NULL) {
          _vtk_read_polydata(f, &prepa);
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "Dataset '%s' not supported\n", word);
        }

      }

      else if (strstr(word, "CELL_DATA") != NULL) {
        _vtk_read_fields(f,
                         &prepa,
                         CELL_DATA);
      }

      else if (strstr(word, "POINT_DATA") != NULL) {
        _vtk_read_fields(f,
                         &prepa,
                         POINT_DATA);
      }

    }

    fclose(f);

  }



  // Redistribute and assemble dmesh nodal
  int          dn_elt[PDM_MESH_NODAL_N_ELEMENT_TYPES] = {0};
  PDM_g_num_t  gn_elt[PDM_MESH_NODAL_N_ELEMENT_TYPES] = {0};
  int          dn_vtx = 0;
  PDM_g_num_t  gn_vtx = 0;
  double      *dvtx_coord = NULL;
  double      *gvtx_coord = NULL;
  PDM_g_num_t *delt_vtx[PDM_MESH_NODAL_N_ELEMENT_TYPES] = {NULL};
  PDM_g_num_t *gelt_vtx[PDM_MESH_NODAL_N_ELEMENT_TYPES] = {NULL};
  int         *gpoly2d_vtx_n = NULL;

  if (i_rank == 0) {
    int idx[PDM_MESH_NODAL_N_ELEMENT_TYPES] = {0};
    for (PDM_g_num_t i = 0; i < prepa.n_elt; i++) {
      dn_elt[prepa.elt_type[i]]++;
      gn_elt[prepa.elt_type[i]]++;
      idx[prepa.elt_type[i]] += prepa.elt_vtx_idx[i+1] - prepa.elt_vtx_idx[i];
    }

    for (PDM_Mesh_nodal_elt_t t = (PDM_Mesh_nodal_elt_t) 0; t < PDM_MESH_NODAL_N_ELEMENT_TYPES; t++) {
      if (gn_elt[t] > 0) {
        gelt_vtx[t] = malloc(sizeof(PDM_g_num_t) * idx[t]);
        if (t == PDM_MESH_NODAL_POLY_2D) {
          gpoly2d_vtx_n = malloc(sizeof(int) * gn_elt[t]);
        }
      }
      idx[t] = 0;
    }

    int ipoly2d = 0;
    for (PDM_g_num_t i = 0; i < prepa.n_elt; i++) {
      for (int j = prepa.elt_vtx_idx[i]; j < prepa.elt_vtx_idx[i+1]; j++) {
        gelt_vtx[prepa.elt_type[i]][idx[prepa.elt_type[i]]++] = prepa.elt_vtx[j];
      }
      if (prepa.elt_type[i] == PDM_MESH_NODAL_POLY_2D) {
        gpoly2d_vtx_n[ipoly2d++] = prepa.elt_vtx_idx[i+1] - prepa.elt_vtx_idx[i];
      }
    }

    dn_vtx = (int) prepa.n_vtx;
    gn_vtx = prepa.n_vtx;
    gvtx_coord = prepa.vtx_coord;

    if (prepa.elt_vtx_idx != NULL) {
      free(prepa.elt_vtx_idx);
    }

    if (prepa.elt_vtx != NULL) {
      free(prepa.elt_vtx);
    }

    if (prepa.elt_type != NULL) {
      free(prepa.elt_type);
    }
  }


  PDM_MPI_Bcast(&gn_vtx, 1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(gn_elt, PDM_MESH_NODAL_N_ELEMENT_TYPES, PDM__PDM_MPI_G_NUM, 0, comm);
  if (dbg_enabled) {
    PDM_log_trace_array_long(gn_elt, PDM_MESH_NODAL_N_ELEMENT_TYPES, "gn_elt :");
  }


  int n_field[2] = {0};
  if (i_rank == 0) {
    n_field[0] = prepa.field[0].n_field;
    n_field[1] = prepa.field[1].n_field;
  }

  PDM_MPI_Bcast(n_field, 2, PDM_MPI_INT, 0, comm);

  if (dbg_enabled) {
    PDM_log_trace_array_int(n_field, 2, "n_field : ");
  }

  for (int i = 0; i < 2; i++) {

    if (n_field[i] > 0) {
      // Bcast field name, stride and type
      int *l_name = malloc(sizeof(int) * n_field[i]);
      if (i_rank == 0) {
        for (int j = 0; j < n_field[i]; j++) {
          l_name[j] = strlen(prepa.field[i].field_name[j]) + 1;
        }
      }
      else {
        prepa.field[i].n_field = n_field[i];
        prepa.field[i].field_name   = malloc(sizeof(char     *) * n_field[i]);
        prepa.field[i].field_type   = malloc(sizeof(PDM_data_t) * n_field[i]);
        prepa.field[i].field_stride = malloc(sizeof(int       ) * n_field[i]);
        prepa.field[i].field_value  = malloc(sizeof(void     *) * n_field[i]);
      }

      PDM_MPI_Bcast(l_name, n_field[i], PDM_MPI_INT, 0, comm);

      int l_char_buf = 0;
      for (int j = 0; j < n_field[i]; j++) {
        l_char_buf += l_name[j];
      }
      char *char_buf = malloc(sizeof(char *) * l_char_buf);
      if (i_rank == 0) {
        int idx = 0;
        for (int j = 0; j < n_field[i]; j++) {
          for (int k = 0; k < l_name[j]-1; k++) {
            char_buf[idx++] = prepa.field[i].field_name[j][k];
          }
          char_buf[idx++] = '\0';
        }
      }

      PDM_MPI_Bcast((void *) char_buf, l_char_buf, PDM_MPI_CHAR, 0, comm);

      if (i_rank != 0) {
        int idx = 0;
        for (int j = 0; j < n_field[i]; j++) {
          prepa.field[i].field_name[j] = malloc(sizeof(char) * l_name[j]);
          for (int k = 0; k < l_name[j]; k++) {
            prepa.field[i].field_name[j][k] = char_buf[idx++];
          }
        }
      }
      free(char_buf);
      free(l_name);

      PDM_MPI_Bcast((void *) prepa.field[i].field_type  , n_field[i], PDM_MPI_INT, 0, comm);
      PDM_MPI_Bcast((void *) prepa.field[i].field_stride, n_field[i], PDM_MPI_INT, 0, comm);



      if (dbg_enabled) {
        log_trace("i = %d\n", i);
        log_trace("n_field = %d\n", prepa.field[i].n_field);
        for (int j = 0; j < prepa.field[i].n_field; j++) {
          log_trace("  %s : type %d, stride %d\n",
                    prepa.field[i].field_name[j],
                    prepa.field[i].field_type[j],
                    prepa.field[i].field_stride[j]);
        }
      }

      if (i == 0) {
        *n_vtx_field      = n_field[i];
        *vtx_field_name   = prepa.field[i].field_name;
        *vtx_field_type   = prepa.field[i].field_type;
        *vtx_field_stride = prepa.field[i].field_stride;
      }
      else {
        *n_elt_field      = n_field[i];
        *elt_field_name   = prepa.field[i].field_name;
        *elt_field_type   = prepa.field[i].field_type;
        *elt_field_stride = prepa.field[i].field_stride;
      }

    }

  }





  int mesh_dimension = -1;
  PDM_g_num_t gn_elt_dim[3] = {0};
  for (PDM_Mesh_nodal_elt_t t = (PDM_Mesh_nodal_elt_t) 0; t < PDM_MESH_NODAL_N_ELEMENT_TYPES; t++) {
    if (gn_elt[t] > 0) {
      int dim = PDM_Mesh_nodal_elt_dim_get(t);
      mesh_dimension = PDM_MAX(mesh_dimension, dim);
      gn_elt_dim[dim-1]++;
    }
  }

  if (dbg_enabled) {
    log_trace("mesh_dimension = %d\n", mesh_dimension);
  }


  /* Vertices*/
  PDM_g_num_t *init_distrib_vtx = PDM_compute_entity_distribution        (comm, dn_vtx);
  PDM_g_num_t *distrib_vtx      = PDM_compute_uniform_entity_distribution(comm, gn_vtx);
  PDM_block_to_block_t *btb_vtx = PDM_block_to_block_create(init_distrib_vtx,
                                                            distrib_vtx,
                                                            comm);
  PDM_block_to_block_exch(btb_vtx,
                          3*sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void  *) gvtx_coord,
                          NULL,
                (void **) &dvtx_coord);

  if (gvtx_coord != NULL) {
    free(gvtx_coord);
  }

  if (*n_vtx_field > 0) {
    *vtx_field_value = malloc(sizeof(void *) * (*n_vtx_field));

    for (int i = 0; i < *n_vtx_field; i++) {
      PDM_block_to_block_exch(btb_vtx,
                              _pdm_data_size(prepa.field[0].field_type[i]),
                              PDM_STRIDE_CST_INTERLACED,
                              prepa.field[0].field_stride[i],
                              NULL,
                    (void  *) prepa.field[0].field_value[i],
                              NULL,
                    (void **) &(*vtx_field_value)[i]);

      if (i_rank == 0) {
        free(prepa.field[0].field_value[i]);
      }
    }
    free(prepa.field[0].field_value);
  }


  PDM_block_to_block_free(btb_vtx);
  free(init_distrib_vtx);



  /* Elements */
  PDM_g_num_t *distrib_elt[PDM_MESH_NODAL_N_ELEMENT_TYPES] = {NULL};
  int *dpoly2d_vtx_n = NULL;

  for (PDM_Mesh_nodal_elt_t t = (PDM_Mesh_nodal_elt_t) 0; t < PDM_MESH_NODAL_N_ELEMENT_TYPES; t++) {
    if (gn_elt[t] > 0) {
      PDM_g_num_t *init_distrib = PDM_compute_entity_distribution        (comm, dn_elt[t]);
      distrib_elt[t]            = PDM_compute_uniform_entity_distribution(comm, gn_elt[t]);

      PDM_block_to_block_t *btb = PDM_block_to_block_create(init_distrib,
                                                            distrib_elt[t],
                                                            comm);

      if (t == PDM_MESH_NODAL_POLY_2D) {
        dpoly2d_vtx_n = malloc(sizeof(int) * (distrib_elt[t][i_rank+1] - distrib_elt[t][i_rank]));
        PDM_block_to_block_exch(btb,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR_INTERLACED,
                                0,
                                gpoly2d_vtx_n,
                      (void  *) gelt_vtx[t],
                                dpoly2d_vtx_n,
                      (void **) &delt_vtx[t]);
        if (gpoly2d_vtx_n != NULL) {
          free(gpoly2d_vtx_n);
        }
      }
      else if (t == PDM_MESH_NODAL_POLY_3D) {
        PDM_error(__FILE__, __LINE__, 0, "Poly3d are not supported\n");
      }
      else {
        int stride = PDM_Mesh_nodal_n_vtx_elt_get(t, 1); // high-order??

        PDM_block_to_block_exch(btb,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_CST_INTERLACED,
                                stride,
                                NULL,
                      (void  *) gelt_vtx[t],
                                NULL,
                      (void **) &delt_vtx[t]);
      }

      if (gelt_vtx[t] != NULL) {
        free(gelt_vtx[t]);
      }
      PDM_block_to_block_free(btb);
      free(init_distrib);
    }
  }


  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  mesh_dimension,
                                                  gn_vtx,
                                                  gn_elt_dim[2],
                                                  gn_elt_dim[1],
                                                  gn_elt_dim[0]);


  /* Vertices */
  dn_vtx = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);
  free(distrib_vtx);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  /* Sections */
  for (PDM_Mesh_nodal_elt_t t = (PDM_Mesh_nodal_elt_t) 0; t < PDM_MESH_NODAL_N_ELEMENT_TYPES; t++) {
    if (gn_elt[t] > 0) {
      dn_elt[t] = (int) (distrib_elt[t][i_rank+1] - distrib_elt[t][i_rank]);

      PDM_dmesh_nodal_elmts_t *dmne = NULL;
      switch (PDM_Mesh_nodal_elt_dim_get(t)) {
      case 0:
        dmne = dmn->corner;
        break;
      case 1:
        dmne = dmn->ridge;
        break;
      case 2:
        dmne = dmn->surfacic;
        break;
      case 3:
        dmne = dmn->volumic;
        break;
      }

      int id_section = PDM_DMesh_nodal_elmts_section_add(dmne,
                                                         t);

      if (t == PDM_MESH_NODAL_POLY_2D) {
        int *dpoly2d_vtx_idx = PDM_array_new_idx_from_sizes_int(dpoly2d_vtx_n, dn_elt[t]);
        free(dpoly2d_vtx_n);

        PDM_DMesh_nodal_elmts_section_poly2d_set(dmne,
                                                 id_section,
                                                 dn_elt[t],
                                                 dpoly2d_vtx_idx,
                                                 delt_vtx[t],
                                                 PDM_OWNERSHIP_KEEP);
        if (dbg_enabled) {
          log_trace("elt type %d\n", t);
          PDM_log_trace_connectivity_long(dpoly2d_vtx_idx,
                                          delt_vtx[t],
                                          dn_elt[t],
                                          "delt_vtx : ");
        }
      }
      else {
        // high-order??
        PDM_DMesh_nodal_elmts_section_std_set(dmne,
                                              id_section,
                                              dn_elt[t],
                                              delt_vtx[t],
                                              PDM_OWNERSHIP_KEEP);
        if (dbg_enabled) {
          log_trace("elt type %d\n", t);
          int *connec_idx = PDM_array_new_idx_from_const_stride_int(PDM_Mesh_nodal_n_vtx_elt_get(t, 1),
                                                                    dn_elt[t]);
          PDM_log_trace_connectivity_long(connec_idx,
                                          delt_vtx[t],
                                          dn_elt[t],
                                          "delt_vtx : ");
          free(connec_idx);
        }
      }


      free(distrib_elt[t]);
    }
  }

  if (*n_elt_field > 0) {
    for (int i = 0; i < *n_elt_field; i++) {
      if (i_rank == 0) {
        free(prepa.field[1].field_value[i]);
      }
    }
    free(prepa.field[1].field_value);
  }


  return dmn;
}
