/*
 * \file
 */

#ifndef __PDM_VTK_H__
#define __PDM_VTK_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/
#include <string.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"

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
 * Fortran function header
 *============================================================================*/

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
);

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
);

/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
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
);

/**
 * \brief Export a set of boxes to ASCII VTK format (unstructured grid of hexahedra)
 *
 * \param [in]  filename     Output file name
 * \param [in]  n_box        Number of boxes
 * \param [in]  box_extents  Extents of the boxes (size = 6 * \ref n_box)
 *                           (xmin0, ymin0, zmin0, xmax0, ymax0, zmax0, xmin1, ...)
 * \param [in]  box_g_num    Global ids of the boxes (or NULL)
 * \param [in]  n_box_field           Number of box fields
 * \param [in]  box_field_name        Name of those box fields
 * \param [in]  box_field             Box fields
 */

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
);


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
);


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
 );


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
 );

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
);


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
 );


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
 );


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
);

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
 );

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
 );


/**
 * \brief Export a set of ellipsed to ASCII VTK format (unstructured grid of line segments)
 *
 * \param [in]  filename     Output file name
 * \param [in]  n_ellipse    Number of ellipses
 * \param [in]  center       Centers of the ellipses (size = 3 * \ref n_circles)
 *                           (x0, y0, z0, x1, ...)
 * \param [in]  axes         Axes of the ellipses (size = 6 * \ref n_circles)
 *                           (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [in]  radii        Radii of the ellipses (size = 2 * \ref n_ellipse)
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
 );


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
 );


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
 );


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
 );


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
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_VTK_H__ */
