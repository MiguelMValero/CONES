#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_multipart.h"
#include "pdm_dcube_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_part_to_block.h"

#include "pdm_part_extension.h"
#include "pdm_vtk.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_triangulate.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -t      <level>  Bounding boxes tolerance.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -octree          Use octree-based method.\n\n"
     "  -dbbree          Use dbbtree-based method.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scotch       Call PT-Scotch.\n\n"
     "  -hilbert         Call Hilbert.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                          argc,
           char                       **argv,
           PDM_g_num_t                 *n_vtx_seg,
           double                      *length,
           double                      *separation_x,
           double                      *separation_y,
           double                      *separation_z,
           int                         *deform,
           double                      *tolerance,
           double                      *marge,
           int                         *n_part,
           PDM_g_num_t                 *n_pts,
           int                         *post,
           int                         *part_method,
           PDM_mesh_location_method_t  *loc_method,
           int                         *disable_uvw,
           int                         *use_tgt_nodes,
           int                         *extension_depth_tgt,
           int                         *extension_depth_src,
           PDM_Mesh_nodal_elt_t        *elt_type)
{
  int i = 1;

  PDM_UNUSED (post);

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sep") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_y = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_z = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-def") == 0) {
      *deform = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *marge = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_pts = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_PART_SPLIT_HILBERT;
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = PDM_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DBBTREE;
    }
    else if (strcmp(argv[i], "-doctree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DOCTREE;
    }
    else if (strcmp(argv[i], "-locate_all_tgt") == 0) {
      *loc_method = PDM_MESH_LOCATION_LOCATE_ALL_TGT;
    }
    else if (strcmp(argv[i], "-no_uvw") == 0) {
      *disable_uvw = 1;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-nodes") == 0) {
      *use_tgt_nodes = 1;
    }
    else if (strcmp(argv[i], "-ext_depth") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_tgt = atoi(argv[i]);
        *extension_depth_src = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ext_depth_tgt") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_tgt = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ext_depth_src") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_src = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static double R[3][3] =
{
  {-0.14547275709949994,  0.8415293589391187 , -0.5202557207618055 },
  { 0.9893622576902102 ,  0.12373586628506748, -0.07649678720582984},
  { 0.                 , -0.5258495730132333 , -0.8505775840931856 }
};

static void _rotate (const int  n_pts,
                     double    *coord)
{

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i];
    double y = coord[3*i+1];
    double z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
    }
  }
}

static void
_unrotate(const int n_pts, double *coord) {
  for (int i = 0 ; i < n_pts ; i++) {
    double x = coord[3 * i];
    double y = coord[3 * i + 1];
    double z = coord[3 * i + 2];

    for (int j = 0 ; j < 3 ; j++) {
      coord[3 * i + j] = R[0][j] * x + R[1][j] * y + R[2][j] * z;
    }
  }
}


static void
_cube_mesh2
(
 const PDM_MPI_Comm            comm,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 const PDM_g_num_t             n_vtx_seg,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  length,
 const int                     deform,
 const int                     part_extension_depth,
 const PDM_Mesh_nodal_elt_t    elt_type,
 int                         **pn_cell,
 int                         **pn_face,
 int                         **pn_vtx,
 int                        ***pcell_face_idx,
 int                        ***pcell_face,
 int                        ***pface_vtx_idx,
 int                        ***pface_vtx,
 double                     ***pvtx_coord,
 PDM_g_num_t                ***pcell_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn,
 PDM_g_num_t                ***pvtx_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_domain = 1;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);


  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  *pface_vtx_idx = malloc(sizeof(int *) * n_part);
  *pface_vtx     = malloc(sizeof(int *) * n_part);

  if (elt_type < PDM_MESH_NODAL_POLY_3D) {
    PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           length,
                                                           xmin,
                                                           ymin,
                                                           zmin,
                                                           elt_type,
                                                           1,
                                                           PDM_OWNERSHIP_KEEP);
    PDM_dcube_nodal_gen_build (dcube);


    PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);

    /*
     * Split mesh
     */
    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
    PDM_multipart_compute(mpart);

    PDM_dcube_nodal_gen_free(dcube);
  }

  else {
    // Polyhedral mesh
    PDM_g_num_t  ng_cell = 0;
    PDM_g_num_t  ng_face = 0;
    PDM_g_num_t  ng_vtx  = 0;
    int          dn_cell = 0;
    int          dn_face = 0;
    int          dn_edge = 0;
    int          dn_vtx  = 0;
    double      *dvtx_coord     = NULL;
    int         *dcell_face_idx = NULL;
    PDM_g_num_t *dcell_face     = NULL;
    int          n_face_group    = 0;
    PDM_g_num_t *dface_cell      = NULL;
    int         *dface_vtx_idx   = NULL;
    PDM_g_num_t *dface_vtx       = NULL;
    int         *dface_group_idx = NULL;
    PDM_g_num_t *dface_group     = NULL;

    PDM_poly_vol_gen (comm,
                      xmin,
                      ymin,
                      zmin,
                      length,
                      length,
                      length,
                      n_vtx_seg,
                      n_vtx_seg,
                      n_vtx_seg,
                      1, // randomize
                      0,
                      &ng_cell,
                      &ng_face,
                      &ng_vtx,
                      &n_face_group,
                      &dn_cell,
                      &dn_face,
                      &dn_vtx,
                      &dcell_face_idx,
                      &dcell_face,
                      &dface_cell,
                      &dface_vtx_idx,
                      &dface_vtx,
                      &dvtx_coord,
                      &dface_group_idx,
                      &dface_group);

    /* Generate dmesh */
    PDM_dmesh_t *dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
                                           dn_cell,
                                           dn_face,
                                           dn_edge,
                                           dn_vtx,
                                           comm);

    PDM_dmesh_vtx_coord_set(dmesh,
                            dvtx_coord,
                            PDM_OWNERSHIP_USER);


    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               dface_vtx,
                               dface_vtx_idx,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               dface_cell,
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_FACE,
                        n_face_group,
                        dface_group,
                        dface_group_idx,
                        PDM_OWNERSHIP_USER);


    PDM_multipart_dmesh_set (mpart, 0, dmesh);


    PDM_multipart_compute(mpart);

    PDM_dmesh_free(dmesh);
    free(dvtx_coord);
    free(dcell_face_idx);
    free(dcell_face);
    free(dface_cell);
    free(dface_vtx_idx);
    free(dface_vtx);
    free(dface_group_idx);
    free(dface_group);
  }



  if (deform) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      double *_vtx_coord = NULL;
      int _n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                    0,
                                                    i_part,
                                                    &_vtx_coord,
                                                    PDM_OWNERSHIP_KEEP);

      _rotate(_n_vtx, _vtx_coord);
    }
  }

  for (int i_part = 0; i_part < n_part; i_part++) {

    int *face_vtx_idx = NULL;
    int *face_vtx     = NULL;
    int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                     &face_vtx_idx,
                                                     &face_vtx,
                                                     PDM_OWNERSHIP_KEEP);

    int *_face_vtx_idx = face_vtx_idx;
    int *_face_vtx     = face_vtx;
    if (face_vtx == NULL) {
      int *_face_edge_idx = NULL;
      int *_face_edge     = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          &_face_edge_idx,
                                          &_face_edge,
                                          PDM_OWNERSHIP_KEEP);

      int *_edge_vtx_idx = NULL;
      int *_edge_vtx     = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &_edge_vtx_idx,
                                          &_edge_vtx,
                                          PDM_OWNERSHIP_KEEP);

      _face_vtx_idx = _face_edge_idx;
      PDM_compute_face_vtx_from_face_and_edge(n_face,
                                              _face_edge_idx,
                                              _face_edge,
                                              _edge_vtx,
                                              &_face_vtx);
    }
    (*pface_vtx_idx)[i_part] = malloc(sizeof(int) * (n_face + 1));
    memcpy((*pface_vtx_idx)[i_part], _face_vtx_idx, sizeof(int) * (n_face + 1));

    (*pface_vtx)[i_part] = malloc(sizeof(int) * _face_vtx_idx[n_face]);
    memcpy((*pface_vtx)[i_part], _face_vtx, sizeof(int) * _face_vtx_idx[n_face]);

    if (_face_vtx != face_vtx) {
      free(_face_vtx);
    }
  }


  PDM_part_extension_t *part_ext = NULL;
  if (part_extension_depth > 0) {
    part_ext = PDM_part_extension_create(1,
                                         &n_part,
                                         PDM_EXTEND_FROM_FACE,
                                         part_extension_depth,
                                         comm,
                                         PDM_OWNERSHIP_KEEP);
    for (int i_part = 0; i_part < n_part; i_part++) {

      PDM_g_num_t* cell_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      int *cell_face     = NULL;
      int *cell_face_idx = NULL;
      int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       PDM_OWNERSHIP_KEEP);


      int *face_cell     = NULL;
      int *face_cell_idx = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                          &face_cell_idx,
                                          &face_cell,
                                          PDM_OWNERSHIP_KEEP);
      assert(face_cell_idx == NULL);

      int *face_vtx     = NULL;
      int *face_vtx_idx = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &face_vtx_idx,
                                          &face_vtx,
                                          PDM_OWNERSHIP_KEEP);

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                       &face_edge_idx,
                                                       &face_edge,
                                                       PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* face_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx_idx,
                                                       &edge_vtx,
                                                       PDM_OWNERSHIP_KEEP);
      assert(edge_vtx_idx == NULL);
      PDM_g_num_t* edge_ln_to_gn = NULL;
      int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart,
                                                    0,
                                                    i_part,
                                                    PDM_MESH_ENTITY_EDGE,
                                                    &edge_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);
      assert(n_edge2 == n_edge);

      PDM_g_num_t* vtx_ln_to_gn = NULL;
      int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                  0,
                                                  i_part,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);
      int          pn_face_group        = 0;
      int         *face_group_idx      = NULL;
      int         *face_group          = NULL;
      PDM_g_num_t *face_group_ln_to_gn = NULL;
      PDM_multipart_group_get(mpart,
                              0,
                              i_part,
                              PDM_MESH_ENTITY_FACE,
                              &pn_face_group,
                              &face_group_idx,
                              &face_group,
                              &face_group_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);

      int *vtx_part_bound_proc_idx = NULL;
      int *vtx_part_bound_part_idx = NULL;
      int *vtx_part_bound          = NULL;
      PDM_multipart_part_graph_comm_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &vtx_part_bound_proc_idx,
                                        &vtx_part_bound_part_idx,
                                        &vtx_part_bound,
                                        PDM_OWNERSHIP_KEEP);

      int *face_part_bound_proc_idx = NULL;
      int *face_part_bound_part_idx = NULL;
      int *face_part_bound          = NULL;
      PDM_multipart_part_graph_comm_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &face_part_bound_proc_idx,
                                        &face_part_bound_part_idx,
                                        &face_part_bound,
                                        PDM_OWNERSHIP_KEEP);
      double *vtx = NULL;
      PDM_multipart_part_vtx_coord_get(mpart,
                                       0,
                                       i_part,
                                       &vtx,
                                       PDM_OWNERSHIP_KEEP);

      PDM_part_extension_set_part(part_ext,
                                  0,
                                  i_part,
                                  n_cell,
                                  n_face,
                                  0, //  USELESS n_face_part_bound,
                                  pn_face_group,
                                  0,   // n_edge
                                  n_vtx,
                                  cell_face_idx,
                                  cell_face,
                                  face_cell,
                                  NULL, // face_edge_idx
                                  NULL, // face_edge
                                  (*pface_vtx_idx)[i_part],
                                  (*pface_vtx    )[i_part],
                                  NULL, //edge_vtx
                                  face_group_idx,
                                  face_group,
                                  NULL, // face_join_idx
                                  NULL, // face_join
                                  face_part_bound_proc_idx,
                                  face_part_bound_part_idx,
                                  face_part_bound,
                                  NULL, // vtx_part_bound_proc_idx
                                  NULL, // vtx_part_bound_part_idx
                                  NULL, // vtx_part_bound
                                  cell_ln_to_gn,
                                  face_ln_to_gn,
                                  NULL, // edge_ln_to_gn
                                  vtx_ln_to_gn,
                                  face_group_ln_to_gn,
                                  vtx);
    }

    PDM_part_extension_compute(part_ext);
  }


  *pn_cell        = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_face        = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_vtx         = (int *)          malloc(sizeof(int *)          * n_part);
  *pcell_face_idx = (int **)         malloc(sizeof(int **)         * n_part);
  *pcell_face     = (int **)         malloc(sizeof(int **)         * n_part);
  *pcell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pface_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_coord     = (double **)      malloc(sizeof(double **)      * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_g_num_t* cell_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    &cell_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    int *cell_face     = NULL;
    int *cell_face_idx = NULL;
    int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                     &cell_face_idx,
                                                     &cell_face,
                                                     PDM_OWNERSHIP_KEEP);


    int *face_cell     = NULL;
    int *face_cell_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                        &face_cell_idx,
                                        &face_cell,
                                        PDM_OWNERSHIP_KEEP);
    assert(face_cell_idx == NULL);

    int *face_vtx     = NULL;
    int *face_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &face_vtx_idx,
                                        &face_vtx,
                                        PDM_OWNERSHIP_KEEP);

    int *face_edge     = NULL;
    int *face_edge_idx = NULL;
    int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                     &face_edge_idx,
                                                     &face_edge,
                                                     PDM_OWNERSHIP_KEEP);

    PDM_g_num_t* face_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    &face_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    int *edge_vtx     = NULL;
    int *edge_vtx_idx = NULL;
    int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                     &edge_vtx_idx,
                                                     &edge_vtx,
                                                     PDM_OWNERSHIP_KEEP);
    assert(edge_vtx_idx == NULL);
    PDM_g_num_t* edge_ln_to_gn = NULL;
    int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart,
                                                  0,
                                                  i_part,
                                                  PDM_MESH_ENTITY_EDGE,
                                                  &edge_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);
    assert(n_edge2 == n_edge);

    PDM_g_num_t* vtx_ln_to_gn = NULL;
    int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                0,
                                                i_part,
                                                PDM_MESH_ENTITY_VTX,
                                                &vtx_ln_to_gn,
                                                PDM_OWNERSHIP_KEEP);
    int          pn_face_group        = 0;
    int         *face_group_idx      = NULL;
    int         *face_group          = NULL;
    PDM_g_num_t *face_group_ln_to_gn = NULL;
    PDM_multipart_group_get(mpart,
                            0,
                            i_part,
                            PDM_MESH_ENTITY_FACE,
                            &pn_face_group,
                            &face_group_idx,
                            &face_group,
                            &face_group_ln_to_gn,
                            PDM_OWNERSHIP_KEEP);

    int *vtx_part_bound_proc_idx = NULL;
    int *vtx_part_bound_part_idx = NULL;
    int *vtx_part_bound          = NULL;
    PDM_multipart_part_graph_comm_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_part_bound_proc_idx,
                                      &vtx_part_bound_part_idx,
                                      &vtx_part_bound,
                                      PDM_OWNERSHIP_KEEP);

    int *face_part_bound_proc_idx = NULL;
    int *face_part_bound_part_idx = NULL;
    int *face_part_bound          = NULL;
    PDM_multipart_part_graph_comm_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_part_bound_proc_idx,
                                      &face_part_bound_part_idx,
                                      &face_part_bound,
                                      PDM_OWNERSHIP_KEEP);
    double *vtx = NULL;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     0,
                                     i_part,
                                     &vtx,
                                     PDM_OWNERSHIP_KEEP);

    int n_ext_vtx  = 0;
    int n_ext_face = 0;
    int n_ext_cell = 0;
    double      *ext_vtx_coord     = NULL;
    PDM_g_num_t *ext_vtx_ln_to_gn  = NULL;
    int         *ext_cell_face     = NULL;
    int         *ext_cell_face_idx = NULL;
    PDM_g_num_t *ext_cell_ln_to_gn = NULL;
    int         *ext_face_vtx      = NULL;
    int         *ext_face_vtx_idx  = NULL;
    PDM_g_num_t *ext_face_ln_to_gn = NULL;

    if (part_extension_depth > 0) {
      /* Vertices */
      n_ext_vtx = PDM_part_extension_vtx_coord_get(part_ext,
                                                   0,
                                                   i_part,
                                                   &ext_vtx_coord);

      PDM_part_extension_ln_to_gn_get(part_ext,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &ext_vtx_ln_to_gn);


      /* Cells */
      n_ext_cell = PDM_part_extension_connectivity_get(part_ext,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &ext_cell_face_idx,
                                                       &ext_cell_face);

      PDM_part_extension_ln_to_gn_get(part_ext,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &ext_cell_ln_to_gn);


      /* Faces */
      n_ext_face = PDM_part_extension_connectivity_get(part_ext,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                       &ext_face_vtx_idx,
                                                       &ext_face_vtx);

      if (ext_face_vtx == NULL) {
        abort();
      //   int *_face_edge_idx = NULL;
      //   int *_face_edge     = NULL;
      //   PDM_part_extension_connectivity_get(part_ext,
      //                                       0,
      //                                       i_part,
      //                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
      //                                       &_face_edge,
      //                                       &_face_edge_idx);

      //   int *_edge_vtx_idx = NULL;
      //   int *_edge_vtx     = NULL;
      //   PDM_part_extension_connectivity_get(part_ext,
      //                                       0,
      //                                       i_part,
      //                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
      //                                       &_edge_vtx,
      //                                       &_edge_vtx_idx);

      //   ext_face_vtx_idx = _face_edge_idx
      //   PDM_compute_face_vtx_from_face_and_edge(n_ext_face,
      //                                           _face_edge_idx,
      //                                           _face_edge,
      //                                           _edge_vtx,
      //                                           &ext_face_vtx);
      }

      PDM_part_extension_ln_to_gn_get(part_ext,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &ext_face_ln_to_gn);
    }

    (*pn_cell)[i_part] = n_cell + n_ext_cell;
    (*pn_face)[i_part] = n_face + n_ext_face;
    (*pn_vtx)[i_part]  = n_vtx  + n_ext_vtx;

    /* Vertices */
    (*pvtx_coord)[i_part] = (double *) malloc(sizeof(double) * 3 * (n_vtx + n_ext_vtx));
    memcpy((*pvtx_coord)[i_part], vtx, sizeof(double) * 3 * n_vtx);

    (*pvtx_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_vtx + n_ext_vtx));
    memcpy((*pvtx_ln_to_gn)[i_part], vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);


    /* Cells */
    (*pcell_face_idx)[i_part] = (int *) malloc(sizeof(int) * (n_cell + n_ext_cell + 1));
    memcpy((*pcell_face_idx)[i_part], cell_face_idx, sizeof(int) * (n_cell + 1));

    int s_cell_face = cell_face_idx[n_cell];
    if (part_extension_depth > 0) {
      s_cell_face += ext_cell_face_idx[n_ext_cell];
    }
    (*pcell_face)[i_part] = (int *) malloc(sizeof(int) * s_cell_face);
    memcpy((*pcell_face)[i_part], cell_face, sizeof(int) * cell_face_idx[n_cell]);

    (*pcell_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_cell + n_ext_cell));
    memcpy((*pcell_ln_to_gn)[i_part], cell_ln_to_gn, sizeof(PDM_g_num_t) * n_cell);


    /* Faces */
    // (*pface_vtx_idx)[i_part] = (int *) malloc(sizeof(int) * (n_face + n_ext_face + 1));
    // memcpy((*pface_vtx_idx)[i_part], face_vtx_idx, sizeof(int) * (n_face + 1));

    int s_face_vtx = (*pface_vtx_idx)[i_part][n_face];
    if (part_extension_depth > 0) {
      s_face_vtx += ext_face_vtx_idx[n_ext_face];
    }
    // (*pface_vtx)[i_part] = (int *) malloc(sizeof(int) * s_face_vtx);
    // memcpy((*pface_vtx)[i_part], face_vtx, sizeof(int) * face_vtx_idx[n_face]);


    (*pface_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_face + n_ext_face));
    memcpy((*pface_ln_to_gn)[i_part], face_ln_to_gn, sizeof(PDM_g_num_t) * n_face);


    if (part_extension_depth > 0) {
      /* Vertices */
      memcpy((*pvtx_coord)[i_part] + 3*n_vtx, ext_vtx_coord, sizeof(double) * 3 * n_ext_vtx);
      memcpy((*pvtx_ln_to_gn)[i_part] + n_vtx, ext_vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_vtx);

      /* Cells */
      for (int i = 1; i <= n_ext_cell; i++) {
        (*pcell_face_idx)[i_part][n_cell + i] = cell_face_idx[n_cell] + ext_cell_face_idx[i];
      }

      memcpy((*pcell_face)[i_part] + cell_face_idx[n_cell],
             ext_cell_face,
             sizeof(int) * ext_cell_face_idx[n_ext_cell]);

      memcpy((*pcell_ln_to_gn)[i_part] + n_cell, ext_cell_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_cell);

      /* Faces */
      (*pface_vtx_idx)[i_part] = realloc((*pface_vtx_idx)[i_part], sizeof(int) * (n_face + n_ext_face + 1));
      for (int i = 1; i <= n_ext_face; i++) {
        (*pface_vtx_idx)[i_part][n_face + i] = (*pface_vtx_idx)[i_part][n_face] + ext_face_vtx_idx[i];
      }

      (*pface_vtx)[i_part] = realloc((*pface_vtx)[i_part], sizeof(int) * s_face_vtx);
      memcpy((*pface_vtx)[i_part] + (*pface_vtx_idx)[i_part][n_face],
             ext_face_vtx,
             sizeof(int) * ext_face_vtx_idx[n_ext_face]);

      memcpy((*pface_ln_to_gn)[i_part] + n_face, ext_face_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_face);


      if (1) {
        char filename[999];
        sprintf(filename, "faces_%d.vtk", n_part*i_rank + i_part);

        int *is_extended = malloc(sizeof(int) * (n_face+n_ext_face));
        for (int i = 0; i < n_face; i++) {
          is_extended[i] = 0;
        }
        for (int i = 0; i < n_ext_face; i++) {
          is_extended[n_face+i] = 1;
        }

        PDM_vtk_write_polydata(filename,
                               (*pn_vtx)        [i_part],
                               (*pvtx_coord)    [i_part],
                               (*pvtx_ln_to_gn) [i_part],
                               (*pn_face)       [i_part],
                               (*pface_vtx_idx) [i_part],
                               (*pface_vtx)     [i_part],
                               (*pface_ln_to_gn)[i_part],
                               is_extended);

        free(is_extended);
      }
    }

    // if(_face_vtx != face_vtx) {
    //   free(_face_vtx);
    // }
  }

  PDM_multipart_free (mpart);
  PDM_part_extension_free (part_ext);
}

static void _export_point_cloud
(
 char         *filename,
 int           n_part,
 int          *n_pts,
 double      **coord,
 PDM_g_num_t **g_num,
 PDM_g_num_t **location
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\npoints\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  int n_pts_t = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    n_pts_t += n_pts[ipart];
  }

  fprintf(f, "POINTS %d double\n", n_pts_t);
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < n_pts[ipart]; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(f, "%f ", coord[ipart][3*i + j]);
      }
      fprintf(f, "\n");
    }
  }

  fprintf(f, "CELLS %d %d\n", n_pts_t, 2*n_pts_t);
  for (int i = 0; i < n_pts_t; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts_t);
  for (int i = 0; i < n_pts_t; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_pts_t);
  fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < n_pts[ipart]; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[ipart][i]);
    }
  }

  if (location != NULL) {
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "location 1 %d int\n", n_pts_t);
    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int i = 0; i < n_pts[ipart]; i++) {
        fprintf(f, ""PDM_FMT_G_NUM"\n", location[ipart][i]);
      }
    }
  }

  fclose(f);
}


static inline double
_eval_field
(
 double *xyz
 )
{
  return 1 + 2*xyz[0] + 3*xyz[1] + 4*xyz[2];
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t n_vtx_seg           = 10;
  double      length              = 1.;
  double      separation_x        = 2.;
  double      separation_y        = 0.;
  double      separation_z        = 0.;
  int         deform              = 0;
  double      tolerance           = 1e-6;
  double      marge               = 0.;
  int         n_part              = 1;
  int         post                = 0;
  int         extension_depth_tgt = 0;
  int         extension_depth_src = 0;
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  PDM_g_num_t n_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;
  int disable_uvw = 0;
  int use_tgt_nodes = 0;

  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_HEXA8;
  //  2 -> tria
  //  3 -> quad
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  // 11 -> tria_ho
  // 12 -> quad_ho
  // 13 -> tetra_ho
  // 14 -> pyramid_ho
  // 15 -> prism_ho
  // 16 -> hexa_ho

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length,
              &separation_x,
              &separation_y,
              &separation_z,
              &deform,
              &tolerance,
              &marge,
              &n_part,
              &n_pts,
              &post,
              (int *) &part_method,
              &loc_method,
              &disable_uvw,
              &use_tgt_nodes,
              &extension_depth_tgt,
              &extension_depth_src,
              &elt_type);


  printf("separation_x : %12.5e / %12.5e / %12.5e \n", separation_x, separation_y, separation_z);

  // assert(//elt_type == PDM_MESH_NODAL_BAR2     ||
  //        //elt_type == PDM_MESH_NODAL_TRIA3    ||
  //        //elt_type == PDM_MESH_NODAL_QUAD4    ||
  //        elt_type == PDM_MESH_NODAL_TETRA4   ||
  //        elt_type == PDM_MESH_NODAL_PYRAMID5 ||
  //        elt_type == PDM_MESH_NODAL_PRISM6   ||
  //        elt_type == PDM_MESH_NODAL_HEXA8);
  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 3);


  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  const double xmin = -0.5*length;//0;
  const double ymin = -0.5*length;//0;
  const double zmin = -0.5*length;//0;

  double xyz_min[3] = {xmin, ymin, zmin};
  double xyz_max[3] = {xmin + length, ymin + length, zmin + length};

  /*
   *  Source cube
   */
  int          *src_n_cell        = NULL;
  int          *src_n_face        = NULL;
  int          *src_n_vtx         = NULL;
  int         **src_cell_face_idx = NULL;
  int         **src_cell_face     = NULL;
  int         **src_face_vtx_idx  = NULL;
  int         **src_face_vtx      = NULL;
  double      **src_vtx_coord     = NULL;
  PDM_g_num_t **src_cell_ln_to_gn = NULL;
  PDM_g_num_t **src_face_ln_to_gn = NULL;
  PDM_g_num_t **src_vtx_ln_to_gn  = NULL;

  _cube_mesh2(comm,
              n_part,
              part_method,
              n_vtx_seg,
              xmin,
              ymin,
              zmin,
              length,
              deform,
              extension_depth_src,
              elt_type,
              &src_n_cell,
              &src_n_face,
              &src_n_vtx,
              &src_cell_face_idx,
              &src_cell_face,
              &src_face_vtx_idx,
              &src_face_vtx,
              &src_vtx_coord,
              &src_cell_ln_to_gn,
              &src_face_ln_to_gn,
              &src_vtx_ln_to_gn);

  /*
   *  Target cube
   */
  int          *tgt_n_cell        = NULL;
  int          *tgt_n_face        = NULL;
  int          *tgt_n_vtx         = NULL;
  int         **tgt_cell_face_idx = NULL;
  int         **tgt_cell_face     = NULL;
  int         **tgt_face_vtx_idx  = NULL;
  int         **tgt_face_vtx      = NULL;
  double      **tgt_vtx_coord     = NULL;
  PDM_g_num_t **tgt_cell_ln_to_gn = NULL;
  PDM_g_num_t **tgt_face_ln_to_gn = NULL;
  PDM_g_num_t **tgt_vtx_ln_to_gn  = NULL;

  _cube_mesh2(comm,
              n_part,
              part_method,
              2 * n_vtx_seg,
              xmin + separation_x,//*length,
              ymin + separation_y,//*length,
              zmin + separation_z,//*length,
              length,
              deform,
              extension_depth_tgt,
              elt_type, //PDM_MESH_NODAL_HEXA8,
              &tgt_n_cell,
              &tgt_n_face,
              &tgt_n_vtx,
              &tgt_cell_face_idx,
              &tgt_cell_face,
              &tgt_face_vtx_idx,
              &tgt_face_vtx,
              &tgt_vtx_coord,
              &tgt_cell_ln_to_gn,
              &tgt_face_ln_to_gn,
              &tgt_vtx_ln_to_gn);

  /*
   *  Mesh location structure initialization
   */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create (1,
                                                            comm,
                                                            PDM_OWNERSHIP_KEEP);

  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,
                                      n_part);

  int          *n_tgt     = NULL;
  PDM_g_num_t **tgt_g_num = NULL;
  double      **tgt_coord = NULL;
  if (use_tgt_nodes) {
    n_tgt     = tgt_n_vtx;
    tgt_g_num = tgt_vtx_ln_to_gn;
    tgt_coord = tgt_vtx_coord;
  }
  else {
    n_tgt     = tgt_n_cell;
    tgt_g_num = tgt_cell_ln_to_gn;
    tgt_coord = (double **) malloc (sizeof(double *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      const int is_oriented = 1;
      double *cell_volume = (double *) malloc(sizeof(double) * n_tgt[ipart]);
      tgt_coord[ipart] = (double *) malloc(sizeof(double) * 3 * n_tgt[ipart]);

      PDM_geom_elem_polyhedra_properties (is_oriented,
                                          tgt_n_cell[ipart],
                                          tgt_n_face[ipart],
                                          tgt_face_vtx_idx[ipart],
                                          tgt_face_vtx[ipart],
                                          tgt_cell_face_idx[ipart],
                                          tgt_cell_face[ipart],
                                          tgt_n_vtx[ipart],
                                          tgt_vtx_coord[ipart],
                                          cell_volume,
                                          tgt_coord[ipart],
                                          NULL,
                                          NULL);

      free (cell_volume);
    }
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_mesh_location_cloud_set (mesh_loc,
                                 0,
                                 ipart,
                                 n_tgt[ipart],
                                 tgt_coord[ipart],
                                 tgt_g_num[ipart]);
  }


  /* Set source mesh */
  PDM_mesh_location_mesh_n_part_set (mesh_loc,
                                          n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_mesh_location_part_set (mesh_loc,
                                ipart,
                                src_n_cell[ipart],
                                src_cell_face_idx[ipart],
                                src_cell_face[ipart],
                                src_cell_ln_to_gn[ipart],
                                src_n_face[ipart],
                                src_face_vtx_idx[ipart],
                                src_face_vtx[ipart],
                                src_face_ln_to_gn[ipart],
                                src_n_vtx[ipart],
                                src_vtx_coord[ipart],
                                src_vtx_ln_to_gn[ipart]);
  }


  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc,
                                   tolerance);

  PDM_mesh_location_method_set (mesh_loc,
                                loc_method);

  /*
   *  Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  PDM_mesh_location_compute(mesh_loc);

  PDM_mesh_location_dump_times(mesh_loc);




  /*
   *  Check result from target PoV
   */
  PDM_g_num_t **tgt_location   = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **tgt_proj_coord = malloc (sizeof(double *)      * n_part);

  int n_wrong = 0;
  const PDM_g_num_t n_cell_seg = n_vtx_seg - 1;
  const double cell_side = length / ((double) n_cell_seg);

  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_located = PDM_mesh_location_n_located_get (mesh_loc,
                                                     0,//i_point_cloud,
                                                     ipart);

    int *located = PDM_mesh_location_located_get (mesh_loc,
                                                  0,//i_point_cloud,
                                                  ipart);

    int n_unlocated = PDM_mesh_location_n_unlocated_get (mesh_loc,
                                                         0,//i_point_cloud,
                                                         ipart);

    int *unlocated = PDM_mesh_location_unlocated_get (mesh_loc,
                                                      0,//i_point_cloud,
                                                      ipart);

    printf("[%d] part %d, n_located = %d, n_unlocated = %d\n", i_rank, ipart, n_located, n_unlocated);
    assert(n_located + n_unlocated == n_tgt[ipart]);

    tgt_location[ipart] = PDM_array_const_gnum (n_tgt[ipart], 0);
    tgt_proj_coord[ipart] = malloc (sizeof(double) * n_tgt[ipart] * 3);

    PDM_g_num_t *p_location    = NULL;
    double      *p_dist2  = NULL;
    double      *p_proj_coord  = NULL;
    PDM_mesh_location_point_location_get (mesh_loc,
                                          0,//i_point_cloud,
                                          ipart,
                                          &p_location,
                                          &p_dist2,
                                          &p_proj_coord);

    for (int j = 0; j < n_located; j++) {
      int i = located[j] - 1;
      tgt_location[ipart][i] = p_location[j];
      for (int k = 0; k < 3; k++) {
        tgt_proj_coord[ipart][3*i+k] = p_proj_coord[3*j+k];
      }
    }

    for (int j = 0; j < n_unlocated; j++) {
      int i = unlocated[j] - 1;
      tgt_location[ipart][i] = -1;
      for (int k = 0; k < 3; k++) {
        tgt_proj_coord[ipart][3*i+k] = tgt_coord[ipart][3*i+k];
      }
    }


    if (elt_type == PDM_MESH_NODAL_HEXA8) {//!deform) {

      for (int k1 = 0; k1 < n_located; k1++) {

        int is_wrong = 0;

        int ipt = located[k1] - 1;
        double *p = tgt_coord[ipart] + 3*ipt;

        double coord[3];
        memcpy(coord, p, sizeof(double)*3);
        if (deform) {
          _unrotate(1, coord);
        }
          for (int i = 0; i < 3; i++) {
            coord[i] = PDM_MIN(PDM_MAX(coord[i], xyz_min[i]), xyz_max[i]);
          }
        if (deform) {
          _rotate(1, coord);
        }

        double err = 0;
        for (int i = 0; i < 3; i++) {
          double delta = coord[i] - p_proj_coord[3*k1+i];
          err += delta*delta;
        }
        err = sqrt(err);
        if (err > tolerance) {
          log_trace("pt %ld (%f %f %f) (%f %f %f) err = %e\n",
                    tgt_g_num[ipart][ipt],
                    p[0], p[1], p[2],
                    coord[0], coord[1], coord[2],
                    err);
          is_wrong = 1;
        }

        int i = (int) floor (coord[0] / cell_side);
        int j = (int) floor (coord[1] / cell_side);
        int k = (int) floor (coord[2] / cell_side);

        PDM_g_num_t box_gnum = 1 + i + n_cell_seg*(j + n_cell_seg*k);

        if (coord[0] < -tolerance || coord[0] > length + tolerance ||
            coord[1] < -tolerance || coord[1] > length + tolerance ||
            coord[2] < -tolerance || coord[2] > length + tolerance) {
          box_gnum = -1;
        }

        if (p_location[k1] != box_gnum) {
          double cell_min[3] = {cell_side * i,     cell_side * j,     cell_side * k};
          double cell_max[3] = {cell_side * (i+1), cell_side * (j+1), cell_side * (k+1)};

          // double dist = HUGE_VAL;
          // for (int idim = 0; idim < 3; idim++) {
          //   double _dist1 = PDM_ABS (coord[idim] - cell_min[idim]);
          //   double _dist2 = PDM_ABS (coord[idim] - cell_max[idim]);
          //   double _dist = PDM_MIN (_dist1, _dist2);
          //   dist = PDM_MIN (dist, _dist);
          // }
          double dist = HUGE_VAL;
          for (int idim = 0; idim < 3; idim++) {
            double _dist1 = 0;
            if (coord[idim] < cell_min[idim]) {
              _dist1 = cell_min[idim] - coord[idim];
            }
            double _dist2 = 0;
            if (coord[idim] > cell_max[idim]) {
              _dist2 = PDM_ABS (coord[idim] - cell_max[idim]);
            }
            double _dist = PDM_MIN (_dist1, _dist2);
            dist = PDM_MIN (dist, _dist);
          }

          if (dist > tolerance) {
            log_trace("pt %ld (%f %f %f) dist = %e\n",
                      tgt_g_num[ipart][ipt],
                      p[0], p[1], p[2],
                      dist);
            is_wrong = 1;
          }
        }

        n_wrong += is_wrong;
      }



      for (int k1 = 0; k1 < n_unlocated; k1++) {
        int ipt = unlocated[k1] - 1;

        double coord[3];
        memcpy(coord, &tgt_coord[ipart][3*ipt], sizeof(double)*3);

        if (deform) {
          _unrotate(1, coord);
        }

        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        if (x >= xmin && x <= xmin + length &&
            y >= ymin && y <= ymin + length &&
            z >= zmin && z <= zmin + length) {
          log_trace("pt %ld should have been located\n", tgt_g_num[ipart][ipt]);
          n_wrong++;
        }
      }

    }
  }

  int g_n_wrong;
  PDM_MPI_Allreduce (&n_wrong, &g_n_wrong, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  if ((i_rank == 0) && (elt_type == PDM_MESH_NODAL_HEXA8)) {
    printf("Viewed from target: g_n_wrong = %d / "PDM_FMT_G_NUM"\n", g_n_wrong, (n_vtx_seg-1)*(n_vtx_seg-1)*(n_vtx_seg-1));
  }









  /*
   *  Check result from source PoV
   */
  if (elt_type == PDM_MESH_NODAL_HEXA8) {//!deform) {

    n_wrong = 0;

    for (int ipart = 0; ipart < n_part; ipart++) {
      int *cell_vtx_idx;
      int *cell_vtx;
      PDM_mesh_location_cell_vertex_get (mesh_loc,
                                         ipart,
                                         &cell_vtx_idx,
                                         &cell_vtx);

      int         *elt_pts_inside_idx;
      PDM_g_num_t *points_gnum;
      double      *points_coords;
      double      *points_uvw;
      int         *points_weights_idx;
      double      *points_weights;
      double      *points_dist2;
      double      *points_projected_coords;

      PDM_mesh_location_points_in_elt_get(mesh_loc,
                                          0,//i_point_cloud,
                                          ipart,
                                          &elt_pts_inside_idx,
                                          &points_gnum,
                                          &points_coords,
                                          &points_uvw,
                                          &points_weights_idx,
                                          &points_weights,
                                          &points_dist2,
                                          &points_projected_coords);
                                           
      // PDM_log_trace_connectivity_long(elt_pts_inside_idx,
      //                                 points_gnum,
      //                                 src_n_cell[ipart],
      //                                 "pie->gnum : ");

      for (int i = 0; i < src_n_cell[ipart]; i++) {

        PDM_g_num_t ck = (src_cell_ln_to_gn[ipart][i] - 1) / (n_cell_seg * n_cell_seg);
        PDM_g_num_t ci = (src_cell_ln_to_gn[ipart][i] - 1) % n_cell_seg;
        PDM_g_num_t cj = (src_cell_ln_to_gn[ipart][i] - 1 - ck*n_cell_seg*n_cell_seg) / n_cell_seg;

        for (int j = elt_pts_inside_idx[i]; j < elt_pts_inside_idx[i+1]; j++) {
          // double *p = points_coords + 3*j;
          double *p = points_projected_coords + 3*j;

          double coord[3];
          memcpy(coord, p, sizeof(double)*3);

          if (deform) {
            _unrotate(1, coord);
          }

          PDM_g_num_t pi = (PDM_g_num_t) floor (coord[0] / cell_side);
          PDM_g_num_t pj = (PDM_g_num_t) floor (coord[1] / cell_side);
          PDM_g_num_t pk = (PDM_g_num_t) floor (coord[2] / cell_side);

          if (ci != pi || cj != pj || ck != pk) {

            double cell_min[3] = {cell_side * ci,     cell_side * cj,     cell_side * ck};
            double cell_max[3] = {cell_side * (ci+1), cell_side * (cj+1), cell_side * (ck+1)};

            double dist = HUGE_VAL;
            for (int idim = 0; idim < 3; idim++) {
              double _dist1 = 0;
              if (coord[idim] < cell_min[idim]) {
                _dist1 = cell_min[idim] - coord[idim];
              }
              double _dist2 = 0;
              if (coord[idim] > cell_max[idim]) {
                _dist2 = PDM_ABS (coord[idim] - cell_max[idim]);
              }
              double _dist = PDM_MIN (_dist1, _dist2);
              dist = PDM_MIN (dist, _dist);
            }

            if (dist > tolerance) {
              // log_trace("!!! part %d, from source cell "PDM_FMT_G_NUM", point "PDM_FMT_G_NUM"\n", ipart, src_g_num[ipart][i], points_gnum[j]);
              log_trace("!!! part %d, from source cell %d, point "PDM_FMT_G_NUM"\n", ipart, i, points_gnum[j]);
              n_wrong++;
            }
          }
        }
      }
    }


    PDM_MPI_Allreduce (&n_wrong, &g_n_wrong, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    if (i_rank == 0) {
      printf("Viewed from source: g_n_wrong = %d / "PDM_FMT_G_NUM"\n", g_n_wrong, (n_vtx_seg-1)*(n_vtx_seg-1)*(n_vtx_seg-1));
    }
  }




  if (post) {
    char filename[999];

    sprintf(filename, "tgt_location_%3.3d.vtk", i_rank);
    _export_point_cloud (filename,
                         n_part,
                         n_tgt,
                         tgt_coord,
                         tgt_g_num,
                         tgt_location);

    sprintf(filename, "tgt_proj_coord_%3.3d.vtk", i_rank);

    _export_point_cloud (filename,
                         n_part,
                         n_tgt,
                         tgt_proj_coord,
                         tgt_g_num,
                         tgt_location);
  }





  /*
   *  Check location (interpolation of an affine field)
   */
  double **src_field = malloc(sizeof(double *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    src_field[ipart] = malloc(sizeof(double) * src_n_vtx[ipart]);
    for (int i = 0; i < src_n_vtx[ipart]; i++) {
      double coord[3];
      memcpy(coord, &src_vtx_coord[ipart][3*i], sizeof(double)*3);
      if (deform) {
        _unrotate(1, coord);
      }
      src_field[ipart][i] = _eval_field(coord);
    }
  }

  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);

  double **send_field = malloc(sizeof(double *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    int         *elt_pts_idx        = NULL;
    PDM_g_num_t *elt_pts_gnum       = NULL;
    double      *elt_pts_coord      = NULL;
    double      *elt_pts_uvw        = NULL;
    int         *elt_pts_weight_idx = NULL;
    double      *elt_pts_weight     = NULL;
    double      *elt_pts_dist2      = NULL;
    double      *elt_pts_proj_coord = NULL;
    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        0, // i_point_cloud,
                                        ipart,
                                        &elt_pts_idx,
                                        &elt_pts_gnum,
                                        &elt_pts_coord,
                                        &elt_pts_uvw,
                                        &elt_pts_weight_idx,
                                        &elt_pts_weight,
                                        &elt_pts_dist2,
                                        &elt_pts_proj_coord);

    int *cell_vtx_idx = NULL;
    int *cell_vtx     = NULL;
    PDM_mesh_location_cell_vertex_get(mesh_loc,
                                      ipart,
                                      &cell_vtx_idx,
                                      &cell_vtx);

    send_field[ipart] = malloc(sizeof(double) * elt_pts_idx[src_n_cell[ipart]]);
    for (int ielt = 0; ielt < src_n_cell[ipart]; ielt++) {

      int *cv = cell_vtx + cell_vtx_idx[ielt];

      for (int idx_pt = elt_pts_idx[ielt]; idx_pt < elt_pts_idx[ielt+1]; idx_pt++) {
        send_field[ipart][idx_pt] = 0.;
        int idx_vtx = 0;
        // double e[3] = {
        //   elt_pts_proj_coord[3*idx_pt  ],
        //   elt_pts_proj_coord[3*idx_pt+1],
        //   elt_pts_proj_coord[3*idx_pt+2]
        // };
        double e[3] = {
          elt_pts_coord[3*idx_pt  ],
          elt_pts_coord[3*idx_pt+1],
          elt_pts_coord[3*idx_pt+2]
        };

        for (int idx_w = elt_pts_weight_idx[idx_pt]; idx_w < elt_pts_weight_idx[idx_pt+1]; idx_w++) {
          int vtx_id = cv[idx_vtx++] - 1;
          send_field[ipart][idx_pt] += elt_pts_weight[idx_w] * src_field[ipart][vtx_id];
          for (int j = 0; j < 3; j++) {
            e[j] -= elt_pts_weight[idx_w] * src_vtx_coord[ipart][3*vtx_id+j];
          }
        }

        // log_trace("pt "PDM_FMT_G_NUM" (%f %f %f), in elt "PDM_FMT_G_NUM" : dist = %e\n",
        //           elt_pts_gnum[idx_pt],
        //           elt_pts_coord[3*idx_pt], elt_pts_coord[3*idx_pt+1], elt_pts_coord[3*idx_pt+2],
        //           src_cell_ln_to_gn[ipart][ielt],
        //           PDM_MODULE(e));
      }

    }
  }


  double **recv_field = NULL;
  int request = -1;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) send_field,
                         NULL,
        (      void ***) &recv_field,
                         &request);

  PDM_part_to_part_iexch_wait(ptp, request);

  double lmax_err = 0.;
  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_located = PDM_mesh_location_n_located_get (mesh_loc,
                                                     0,//i_point_cloud,
                                                     ipart);

    int *located = PDM_mesh_location_located_get (mesh_loc,
                                                  0,//i_point_cloud,
                                                  ipart);

    for (int i = 0; i < n_located; i++) {
      int pt_id = located[i] - 1;

      double coord[3];
      memcpy(coord, &tgt_coord[ipart][3*pt_id], sizeof(double)*3);
      if (deform) {
        _unrotate(1, coord);
      }

      for (int j = 0; j < 3; j++) {
        coord[j] = PDM_MAX(coord[j], xyz_min[j]);
        coord[j] = PDM_MIN(coord[j], xyz_max[j]);
      }

      double f = _eval_field(coord);

      double err = PDM_ABS(recv_field[ipart][i] - f);
      lmax_err = PDM_MAX(lmax_err, err);

      if (err > 1.e-12) {
        log_trace("point "PDM_FMT_G_NUM" (%f %f %f) located in elt "PDM_FMT_G_NUM" : error = %e (%20.16f / %20.16f)\n",
                  tgt_g_num[ipart][pt_id],
                  tgt_coord[ipart][3*pt_id], tgt_coord[ipart][3*pt_id+1], tgt_coord[ipart][3*pt_id+2],
                  tgt_location[ipart][pt_id], err, recv_field[ipart][i], f);
      }
    }
    free(recv_field[ipart]);
  }
  free(recv_field);


  double gmax_err;
  PDM_MPI_Allreduce(&lmax_err, &gmax_err, 1, PDM_MPI_DOUBLE,
                    PDM_MPI_MAX, PDM_MPI_COMM_WORLD);


  if (i_rank == 0) {
    printf("global max interpolation error = %e\n", gmax_err);
  }
  if (gmax_err > tolerance) { // scaling??
    PDM_error(__FILE__, __LINE__, 0, "Large interpolation error!\n");
  }




  /*
   *  Free memory
   */
  for (int ipart = 0; ipart < n_part; ipart++) {
    free(src_cell_face_idx[ipart]);
    free(src_cell_face    [ipart]);
    free(src_face_vtx_idx [ipart]);
    free(src_face_vtx     [ipart]);
    free(src_vtx_coord    [ipart]);
    free(src_cell_ln_to_gn[ipart]);
    free(src_face_ln_to_gn[ipart]);
    free(src_vtx_ln_to_gn [ipart]);

    free(tgt_cell_face_idx[ipart]);
    free(tgt_cell_face    [ipart]);
    free(tgt_face_vtx_idx [ipart]);
    free(tgt_face_vtx     [ipart]);
    free(tgt_vtx_coord    [ipart]);
    free(tgt_cell_ln_to_gn[ipart]);
    free(tgt_face_ln_to_gn[ipart]);
    free(tgt_vtx_ln_to_gn [ipart]);

    free (tgt_location[ipart]);
    free (tgt_proj_coord[ipart]);

    if (!use_tgt_nodes) {
      free(tgt_coord[ipart]);
    }

    free(src_field [ipart]);
    free(send_field[ipart]);
  }

  free(src_n_cell       );
  free(src_n_face       );
  free(src_n_vtx        );
  free(src_cell_face_idx);
  free(src_cell_face    );
  free(src_face_vtx_idx );
  free(src_face_vtx     );
  free(src_vtx_coord    );
  free(src_cell_ln_to_gn);
  free(src_face_ln_to_gn);
  free(src_vtx_ln_to_gn );

  free(tgt_n_cell       );
  free(tgt_n_face       );
  free(tgt_n_vtx        );
  free(tgt_cell_face_idx);
  free(tgt_cell_face    );
  free(tgt_face_vtx_idx );
  free(tgt_face_vtx     );
  free(tgt_vtx_coord    );
  free(tgt_cell_ln_to_gn);
  free(tgt_face_ln_to_gn);
  free(tgt_vtx_ln_to_gn );

  free (tgt_location);
  free (tgt_proj_coord);

  if (!use_tgt_nodes) {
    free(tgt_coord);
  }

  free(src_field);
  free(send_field);

  PDM_mesh_location_free(mesh_loc);
  PDM_part_to_part_free (ptp);


  double min_elaps_create;
  double max_elaps_create;
  double min_cpu_create;
  double max_cpu_create;
  double min_elaps_create2;
  double max_elaps_create2;
  double min_cpu_create2;
  double max_cpu_create2;
  double min_elaps_exch;
  double max_elaps_exch;
  double min_cpu_exch;
  double max_cpu_exch;

  PDM_part_to_block_global_timer_get (PDM_MPI_COMM_WORLD,
                                      &min_elaps_create,
                                      &max_elaps_create,
                                      &min_cpu_create,
                                      &max_cpu_create,
                                      &min_elaps_create2,
                                      &max_elaps_create2,
                                      &min_cpu_create2,
                                      &max_cpu_create2,
                                      &min_elaps_exch,
                                      &max_elaps_exch,
                                      &min_cpu_exch,
                                      &max_cpu_exch);

  if (i_rank == 0) {
    printf("Global time in PDM_part_to_block : \n");
    printf("   - min max elaps create  : %12.5e %12.5e\n", min_elaps_create, max_elaps_create);
    printf("   - min max elaps create2 : %12.5e %12.5e\n", min_elaps_create2, max_elaps_create2);
    printf("   - min max elaps exch    : %12.5e %12.5e\n", min_elaps_exch, max_elaps_exch);
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return g_n_wrong || gmax_err > 1e-12;
}

