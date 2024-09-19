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
#include "pdm_reader_gamma.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_part_to_part.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"

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
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                         argc,
           char                      **argv,
           char                      **filename,
           int                        *post,
           int                        *n_part,
           PDM_g_num_t                *gn_pts,
           PDM_split_dual_t           *part_method,
           PDM_mesh_location_method_t *loc_method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }

    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
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
        long n = atol(argv[i]);
        *gn_pts = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
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

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static void
_dump_mesh
(
 const char     *name,
 PDM_MPI_Comm    comm,
 int             n_part,
 int            *pn_vtx,
 double        **pvtx_coord,
 PDM_g_num_t   **pvtx_ln_to_gn,
 int            *pn_face,
 int           **pface_vtx_idx,
 int           **pface_vtx,
 int            *pn_cell,
 int           **pcell_face_idx,
 int           **pcell_face,
 PDM_g_num_t   **pcell_ln_to_gn,
 const int       n_elt_field,
 const char    **elt_field_name,
 double       ***elt_field_value,
 const int       n_vtx_field,
 const char    **vtx_field_name,
 double       ***vtx_field_value
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        name,
                                        name,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       n_part);

  int id_var_num_part = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var_elt[n_elt_field];
  for (int i = 0; i < n_elt_field; i++) {
    id_var_elt[i] = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          elt_field_name[i]);
  }

  int id_var_vtx[n_vtx_field];
  for (int i = 0; i < n_vtx_field; i++) {
    id_var_vtx[i] = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_VERTICES,
                                          vtx_field_name[i]);
  }


  PDM_writer_step_beg(wrt, 0.);


  PDM_real_t **val_num_part = malloc(sizeof(PDM_real_t  *) * n_part);
  int **pface_vtx_n  = malloc(sizeof(int *) * n_part);
  int **pcell_face_n = malloc(sizeof(int *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              ipart,
                              pn_vtx[ipart],
                              pvtx_coord[ipart],
                              pvtx_ln_to_gn[ipart],
                              PDM_OWNERSHIP_USER);

    pface_vtx_n [ipart] = (int *) malloc(sizeof(int) * pn_face[ipart]);
    pcell_face_n[ipart] = (int *) malloc(sizeof(int) * pn_cell[ipart]);

    for (int i = 0; i < pn_cell[ipart]; i++) {
      pcell_face_n[ipart][i] = pcell_face_idx[ipart][i+1] - pcell_face_idx[ipart][i];
    }

    for (int i = 0; i < pn_face[ipart]; i++) {
      pface_vtx_n[ipart][i] = pface_vtx_idx[ipart][i+1] - pface_vtx_idx[ipart][i];
    }

    PDM_writer_geom_cell3d_cellface_add(wrt,
                                        id_geom,
                                        ipart,
                                        pn_cell       [ipart],
                                        pn_face       [ipart],
                                        pface_vtx_idx [ipart],
                                        pface_vtx_n   [ipart],
                                        pface_vtx     [ipart],
                                        pcell_face_idx[ipart],
                                        pcell_face_n  [ipart],
                                        pcell_face    [ipart],
                                        pcell_ln_to_gn[ipart]);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    val_num_part[ipart] = malloc(sizeof(PDM_real_t) * pn_cell[ipart]);
    for (int i = 0; i < pn_cell[ipart]; i++) {
      val_num_part[ipart][i] = n_part*i_rank + ipart;
    }

    PDM_writer_var_set(wrt,
                       id_var_num_part,
                       id_geom,
                       ipart,
                       val_num_part[ipart]);

    for (int i = 0; i < n_elt_field; i++) {
      PDM_writer_var_set(wrt,
                         id_var_elt[i],
                         id_geom,
                         ipart,
                         elt_field_value[i][ipart]);
    }

    for (int i = 0; i < n_vtx_field; i++) {
      PDM_writer_var_set(wrt,
                         id_var_vtx[i],
                         id_geom,
                         ipart,
                         vtx_field_value[i][ipart]);
    }
  }

  PDM_writer_geom_write(wrt,
                        id_geom);

  PDM_writer_var_write(wrt,
                       id_var_num_part);
  PDM_writer_var_free(wrt,
                      id_var_num_part);

  for (int i = 0; i < n_elt_field; i++) {
    PDM_writer_var_write(wrt,
                         id_var_elt[i]);
    PDM_writer_var_free(wrt,
                        id_var_elt[i]);
  }

  for (int i = 0; i < n_vtx_field; i++) {
    PDM_writer_var_write(wrt,
                         id_var_vtx[i]);
    PDM_writer_var_free(wrt,
                        id_var_vtx[i]);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(val_num_part[ipart]);
    free(pface_vtx_n [ipart]);
    free(pcell_face_n[ipart]);
  }
  free(val_num_part);
  free(pface_vtx_n );
  free(pcell_face_n);

  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);
}




static void
_dump_point_cloud
(
 const char     *name,
 PDM_MPI_Comm    comm,
 int             n_part,
 int            *pn_pts,
 double        **ppts_coord,
 PDM_g_num_t   **ppts_ln_to_gn,
 const int       n_field,
 const char    **field_name,
 double       ***field_value
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        name,
                                        name,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       n_part);

  int id_var_num_part = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var[n_field];
  for (int i = 0; i < n_field; i++) {
    id_var[i] = PDM_writer_var_create(wrt,
                                      PDM_WRITER_OFF,
                                      PDM_WRITER_VAR_SCALAR,
                                      PDM_WRITER_VAR_ELEMENTS,
                                      field_name[i]);
  }

  PDM_writer_step_beg(wrt, 0.);

  int id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom,
                                          PDM_WRITER_POINT,
                                          PDM_OWNERSHIP_USER);

  PDM_real_t **val_num_part = malloc(sizeof(PDM_real_t  *) * n_part);
  int **connec = malloc(sizeof(int *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              ipart,
                              pn_pts[ipart],
                              ppts_coord[ipart],
                              ppts_ln_to_gn[ipart],
                              PDM_OWNERSHIP_USER);

    connec[ipart] = malloc(sizeof(int) * pn_pts[ipart]);
    for (int i = 0; i < pn_pts[ipart]; i++) {
      connec[ipart][i] = i+1;
    }
    PDM_writer_geom_bloc_std_set(wrt,
                                 id_geom,
                                 id_block,
                                 ipart,
                                 pn_pts[ipart],
                                 connec[ipart],
                                 ppts_ln_to_gn[ipart]);

    val_num_part[ipart] = malloc(sizeof(PDM_real_t) * pn_pts[ipart]);
    for (int i = 0; i < pn_pts[ipart]; i++) {
      val_num_part[ipart][i] = n_part*i_rank + ipart;
    }

    PDM_writer_var_set(wrt,
                       id_var_num_part,
                       id_geom,
                       ipart,
                       val_num_part[ipart]);

    for (int i = 0; i < n_field; i++) {
      PDM_writer_var_set(wrt,
                         id_var[i],
                         id_geom,
                         ipart,
                         field_value[i][ipart]);
    }
  }

  PDM_writer_geom_write(wrt,
                        id_geom);

  PDM_writer_var_write(wrt,
                       id_var_num_part);
  PDM_writer_var_free(wrt,
                      id_var_num_part);

  for (int i = 0; i < n_field; i++) {
    PDM_writer_var_write(wrt,
                         id_var[i]);
    PDM_writer_var_free(wrt,
                        id_var[i]);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(val_num_part[ipart]);
    free(connec[ipart]);
  }
  free(val_num_part);
  free(connec);

  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);
}




static double
_eval_field
(
 const double x,
 const double y,
 const double z
 )
{
  return x + 2*y + 3*z;
}



int main(int argc, char *argv[])
{
  /*
   *  Init
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /*
   *  Read args
   */
  char       *filename = NULL;
  int         post     = 0;
  int         n_part   = 1;
  PDM_g_num_t gn_pts = 10;

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;

  _read_args(argc,
             argv,
             &filename,
             &post,
             &n_part,
             &gn_pts,
             &part_method,
             &loc_method);

  if (filename == NULL) {
    filename = (char *) PDM_MESH_DIR"box.mesh";
  }

  PDM_dmesh_nodal_t *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                                        filename,
                                                        0,
                                                        0);



  double l_extents[6] = {
    HUGE_VAL, HUGE_VAL, HUGE_VAL,
    -HUGE_VAL, -HUGE_VAL, -HUGE_VAL
  };

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  for (int i = 0; i < dn_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      double x = dvtx_coord[3*i + j];
      l_extents[j  ] = PDM_MIN(l_extents[j  ], x);
      l_extents[3+j] = PDM_MAX(l_extents[3+j], x);
    }
  }

  double g_extents[6];
  PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);



  int n_domain = 1;
  int n_part_mesh = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_mesh,
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

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);
  PDM_DMesh_nodal_free(dmn);







  /*
   *  Generate point cloud
   */
  int n_part_cloud = 1;
  int          *pn_pts        = malloc(sizeof(int          ) * n_part_cloud);
  double      **ppts_coord    = malloc(sizeof(double      *) * n_part_cloud);
  PDM_g_num_t **ppts_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part_cloud);
  PDM_point_cloud_gen_random(comm,
                             0, // seed
                             0, // geometric_g_num
                             gn_pts,
                             g_extents[0],
                             g_extents[1],
                             g_extents[2],
                             g_extents[3],
                             g_extents[4],
                             g_extents[5],
                             &pn_pts[0],
                             &ppts_coord[0],
                             &ppts_ln_to_gn[0]);

  // PDM_point_cloud_gen_cartesian(comm,
  //                               nx,
  //                               ny,
  //                               nz,
  //                               g_extents[0],
  //                               g_extents[1],
  //                               g_extents[2],
  //                               g_extents[3],
  //                               g_extents[4],
  //                               g_extents[5],
  //                               &n_pts,
  //                               &pts_coord,
  //                               &pts_ln_to_gn);


  /*
   *  Mesh location
   */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create(1,
                                                           comm,
                                                           PDM_OWNERSHIP_KEEP);

  /* Set point cloud */
  PDM_mesh_location_n_part_cloud_set(mesh_loc,
                                     0,//i_point_cloud,
                                     n_part_cloud);

  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    PDM_mesh_location_cloud_set(mesh_loc,
                                0,//i_point_cloud,
                                ipart,
                                pn_pts[ipart],
                                ppts_coord[ipart],
                                ppts_ln_to_gn[ipart]);
  }

  /* Set mesh */
  PDM_mesh_location_mesh_n_part_set(mesh_loc,
                                         n_part_mesh);

  int **pcell_face_idx = malloc(sizeof(int *) * n_part_mesh);
  int **pcell_face     = malloc(sizeof(int *) * n_part_mesh);
  int **pface_vtx_idx  = malloc(sizeof(int *) * n_part_mesh);
  int **pface_vtx      = malloc(sizeof(int *) * n_part_mesh);

  int          *pn_face       = malloc(sizeof(int          ) * n_part_mesh);
  int          *pn_elt        = malloc(sizeof(int          ) * n_part_mesh);
  PDM_g_num_t **pelt_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part_mesh);
  int          *pn_vtx        = malloc(sizeof(int          ) * n_part_mesh);
  PDM_g_num_t **pvtx_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part_mesh);
  double      **pvtx_coord    = malloc(sizeof(double      *) * n_part_mesh);

  int use_edge = 0;

  for (int ipart = 0; ipart < n_part_mesh; ipart++) {

    int *cell_face_idx;
    int *cell_face;
    int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     ipart,
                                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                     &cell_face_idx,
                                                     &cell_face,
                                                     PDM_OWNERSHIP_KEEP);

    int *face_vtx_idx;
    int *face_vtx;
    int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     ipart,
                                                     PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                     &face_vtx_idx,
                                                     &face_vtx,
                                                     PDM_OWNERSHIP_KEEP);



    if (face_vtx == NULL) {
      use_edge = 1;

      int *face_edge;
      int *face_edge_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          &face_edge_idx,
                                          &face_edge,
                                          PDM_OWNERSHIP_KEEP);

      int *edge_vtx;
      int *edge_vtx_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &edge_vtx_idx,
                                          &edge_vtx,
                                          PDM_OWNERSHIP_KEEP);

      pface_vtx_idx[ipart] = face_edge_idx;
      PDM_compute_face_vtx_from_face_and_edge(n_face,
                                              face_edge_idx,
                                              face_edge,
                                              edge_vtx,
                                              &pface_vtx[ipart]);
    }
    else {
      pface_vtx_idx[ipart] = face_vtx_idx;
      pface_vtx    [ipart] = face_vtx;
    }

    double *vtx_coord;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                 0,
                                                 ipart,
                                                 &vtx_coord,
                                                 PDM_OWNERSHIP_KEEP);
    PDM_g_num_t *cell_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    ipart,
                                    PDM_MESH_ENTITY_CELL,
                                    &cell_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *face_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    ipart,
                                    PDM_MESH_ENTITY_FACE,
                                    &face_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *vtx_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    ipart,
                                    PDM_MESH_ENTITY_VTX,
                                    &vtx_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);


    double *volume = malloc(sizeof(double) * n_cell);
    double *center = malloc(sizeof(double) * n_cell * 3);
    PDM_geom_elem_polyhedra_properties(1,
                                       n_cell,
                                       n_face,
                                       pface_vtx_idx[ipart],
                                       pface_vtx[ipart],
                                       cell_face_idx,
                                       cell_face,
                                       n_vtx,
                                       vtx_coord,
                                       volume,
                                       center,
                                       NULL,
                                       NULL);
    free(center);

    for (int i = 0; i < n_cell; i++) {
      if (volume[i] < 0) {
        // log_trace("flip cell %d (%ld)\n", i, cell_ln_to_gn[i]);
        for (int idx = cell_face_idx[i]; idx < cell_face_idx[i+1]; idx++) {
          cell_face[idx] *= -1;
        }
      }
    }
    free(volume);



    PDM_mesh_location_part_set(mesh_loc,
                               ipart,
                               n_cell,
                               cell_face_idx,
                               cell_face,
                               cell_ln_to_gn,
                               n_face,
                               pface_vtx_idx[ipart],
                               pface_vtx[ipart],
                               face_ln_to_gn,
                               n_vtx,
                               vtx_coord,
                               vtx_ln_to_gn);

    pn_elt       [ipart] = n_cell;
    pelt_ln_to_gn[ipart] = cell_ln_to_gn;

    pn_vtx       [ipart] = n_vtx;
    pvtx_ln_to_gn[ipart] = vtx_ln_to_gn;
    pvtx_coord   [ipart] = vtx_coord;

    pcell_face_idx[ipart] = cell_face_idx;
    pcell_face    [ipart] = cell_face;

    pn_face       [ipart] = n_face;
  }


  /* Set location parameters */
  PDM_mesh_location_tolerance_set(mesh_loc,
                                  1e-6);

  PDM_mesh_location_method_set(mesh_loc,
                               loc_method);

  /*
   * Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }


  PDM_mesh_location_compute(mesh_loc);


  /*
   *  Get and exchange
   */
  /*
   *  Exchange data between point cloud and mesh
   */
  int                 **pelt_pts_idx        = malloc(sizeof(int         *) * n_part_mesh);
  PDM_g_num_t         **pelt_pts_gnum       = malloc(sizeof(PDM_g_num_t *) * n_part_mesh);
  double              **pelt_pts_coord      = malloc(sizeof(double      *) * n_part_mesh);
  double              **pelt_pts_uvw        = malloc(sizeof(double      *) * n_part_mesh);
  int                 **pelt_pts_weight_idx = malloc(sizeof(int         *) * n_part_mesh);
  double              **pelt_pts_weight     = malloc(sizeof(double      *) * n_part_mesh);
  double              **pelt_pts_dist2      = malloc(sizeof(double      *) * n_part_mesh);
  double              **pelt_pts_proj_coord = malloc(sizeof(double      *) * n_part_mesh);

  int                 **pelt_vtx_idx        = malloc(sizeof(int         *) * n_part_mesh);
  int                 **pelt_vtx            = malloc(sizeof(int         *) * n_part_mesh);

  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        ipart,
                                        0, // i_point_cloud,
                                        &pelt_pts_idx       [ipart],
                                        &pelt_pts_gnum      [ipart],
                                        &pelt_pts_coord     [ipart],
                                        &pelt_pts_uvw       [ipart],
                                        &pelt_pts_weight_idx[ipart],
                                        &pelt_pts_weight    [ipart],
                                        &pelt_pts_dist2     [ipart],
                                        &pelt_pts_proj_coord[ipart]);

    PDM_mesh_location_cell_vertex_get(mesh_loc,
                                      ipart,
                                      &pelt_vtx_idx[ipart],
                                      &pelt_vtx    [ipart]);

    if (0) {
      log_trace("\n\n\n\n");
      // PDM_log_trace_array_long(pelt_ln_to_gn[ipart], pn_elt[ipart], "pelt_ln_to_gn : ");
      // PDM_log_trace_connectivity_long(pelt_pts_idx[ipart],
      //                                 pelt_pts_gnum[ipart],
      //                                 pn_elt[ipart],
      //                                 "pelt_pts_gnum : ");
      for (int i = 0; i < pn_elt[ipart]; i++) {
        if (pelt_pts_idx[ipart][i+1] > pelt_pts_idx[ipart][i]) {
          log_trace("elt %ld: pts ", pelt_ln_to_gn[ipart][i]);
          PDM_log_trace_array_long(pelt_pts_gnum[ipart] + pelt_pts_idx[ipart][i],
                                   pelt_pts_idx[ipart][i+1] - pelt_pts_idx[ipart][i],
                                   "");
        }
      }
      PDM_log_trace_array_int(pelt_pts_weight_idx[ipart], pelt_pts_idx[ipart][pn_elt[ipart]], "pelt_pts_weight_idx : ");
    }
  }

  if (0) {
    for (int ipart = 0; ipart < n_part_cloud; ipart++) {
      int n_located = PDM_mesh_location_n_located_get(mesh_loc,
                                                      0,
                                                      ipart);

      int *located = PDM_mesh_location_located_get(mesh_loc,
                                                   0,
                                                   ipart);

      PDM_g_num_t *location;
      double      *dist2;
      double      *projected_coord;
      PDM_mesh_location_point_location_get(mesh_loc,
                                           0,
                                           ipart,
                                           &location,
                                           &dist2,
                                           &projected_coord);
      for (int i = 0; i < n_located; i++) {
        int j = located[i]-1;
        log_trace("point %ld located in elt %ld\n", ppts_ln_to_gn[ipart][j], location[i]);
      }
    }
  }


  /* Create exchange protocol */
  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);

  if (ptp == NULL) {
    ptp = PDM_part_to_part_create((const PDM_g_num_t **) pelt_ln_to_gn,
                                  (const int          *) pn_elt,
                                  n_part_mesh,
                                  (const PDM_g_num_t **) ppts_ln_to_gn,
                                  (const int          *) pn_pts,
                                  n_part_cloud,
                                  (const int         **) pelt_pts_idx,
                                  (const PDM_g_num_t **) pelt_pts_gnum,
                                  comm);
  }


  /*
   *  Exchange mesh (cell-centered) -> point cloud
   */
  double **pelt_field1 = malloc(sizeof(double *) * n_part_mesh);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    pelt_field1[ipart] = malloc(sizeof(double) * pn_elt[ipart]);

    for (int i = 0; i < pn_elt[ipart]; i++) {
      double y = 0.;

      for (int idx = pelt_vtx_idx[ipart][i]; idx < pelt_vtx_idx[ipart][i+1]; idx++) {
        int vtx_id = pelt_vtx[ipart][idx] - 1;
        y += pvtx_coord[ipart][3*vtx_id+1];
      }

      if (pelt_vtx_idx[ipart][i+1] > pelt_vtx_idx[ipart][i]) {
        y /= (double) (pelt_vtx_idx[ipart][i+1] - pelt_vtx_idx[ipart][i]);
      }

      pelt_field1[ipart][i] = sin(5*y);
    }
  }

  int request_elt_pts;
  double **tmp_ppts_field = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) pelt_field1,
                         NULL,
        (      void ***) &tmp_ppts_field,
                         &request_elt_pts);

  PDM_part_to_part_iexch_wait(ptp, request_elt_pts);


  int  *n_located;
  int **located;
  PDM_part_to_part_ref_lnum2_get(ptp,
                                 &n_located,
                                 &located);

  double **ppts_field1 = malloc(sizeof(double *) * n_part_cloud);
  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    ppts_field1[ipart] = malloc(sizeof(double) * pn_pts[ipart]);

    for (int i = 0; i < pn_pts[ipart]; i++) {
      ppts_field1[ipart][i] = 0.;
    }

    for (int i = 0; i < n_located[ipart]; i++) {
      int pt_id = located[ipart][i] - 1;

      ppts_field1[ipart][pt_id] = tmp_ppts_field[ipart][i];
    }

    free(tmp_ppts_field[ipart]);
  }
  free(tmp_ppts_field);



  /*
   *  Exchange mesh <- point cloud
   */
  double **ppts_field2 = malloc(sizeof(double *) * n_part_cloud);
  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    ppts_field2[ipart] = malloc(sizeof(double) * pn_pts[ipart]);

    for (int i = 0; i < pn_pts[ipart]; i++) {
      ppts_field2[ipart][i] = cos(3*ppts_coord[ipart][3*i]);
    }
  }

  int request_pts_elt;
  double **tmp_pelt_field = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(double),
                                 NULL,
                (const void  **) ppts_field2,
                                 NULL,
                (      void ***) &tmp_pelt_field,
                                 &request_pts_elt);

  PDM_part_to_part_reverse_iexch_wait(ptp, request_pts_elt);

  double **pelt_field2 = malloc(sizeof(double *) * n_part_mesh);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    pelt_field2[ipart] = malloc(sizeof(double) * pn_elt[ipart]);

    for (int i = 0; i < pn_elt[ipart]; i++) {
      pelt_field2[ipart][i] = 0.;

      for (int idx = pelt_pts_idx[ipart][i]; idx < pelt_pts_idx[ipart][i+1]; idx++) {
        pelt_field2[ipart][i] += tmp_pelt_field[ipart][idx];
      }

      if (pelt_pts_idx[ipart][i+1] > pelt_pts_idx[ipart][i]) {
        pelt_field2[ipart][i] /= (double) (pelt_pts_idx[ipart][i+1] - pelt_pts_idx[ipart][i]);
      }
    }

    free(tmp_pelt_field[ipart]);
  }
  free(tmp_pelt_field);



  /*
   *  Exchange mesh (node-centered) -> point cloud
   */

  double **pvtx_field3      = malloc(sizeof(double *) * n_part_mesh);
  double **pelt_pts_field3 = malloc(sizeof(double *) * n_part_mesh);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    pvtx_field3     [ipart] = malloc(sizeof(double) * pn_vtx[ipart]);
    for (int i = 0; i < pn_vtx[ipart]; i++) {
      pvtx_field3[ipart][i] = _eval_field(pvtx_coord[ipart][3*i  ],
                                          pvtx_coord[ipart][3*i+1],
                                          pvtx_coord[ipart][3*i+2]);//cos(6*pvtx_coord[ipart][3*i+2]);
    }

    pelt_pts_field3[ipart] = malloc(sizeof(double) * pelt_pts_idx[ipart][pn_elt[ipart]]);
    for (int i = 0; i < pn_elt[ipart]; i++) {
      int *ev = pelt_vtx[ipart] + pelt_vtx_idx[ipart][i];

      for (int idx_pt = pelt_pts_idx[ipart][i]; idx_pt < pelt_pts_idx[ipart][i+1]; idx_pt++) {
        pelt_pts_field3[ipart][idx_pt] = 0.;

        int idx_vtx = 0;
        for (int idx_w = pelt_pts_weight_idx[ipart][idx_pt]; idx_w < pelt_pts_weight_idx[ipart][idx_pt+1]; idx_w++) {
          pelt_pts_field3[ipart][idx_pt] += pelt_pts_weight[ipart][idx_w] * pvtx_field3[ipart][ev[idx_vtx++]-1];
        }
      }
    }
  }

  int request_elt_pts3;
  double **tmp_ppts_field3 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) pelt_pts_field3,
                         NULL,
        (      void ***) &tmp_ppts_field3,
                         &request_elt_pts3);

  PDM_part_to_part_iexch_wait(ptp, request_elt_pts3);


  double max_err = 0.;

  double **ppts_field3 = malloc(sizeof(double *) * n_part_cloud);
  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    ppts_field3[ipart] = malloc(sizeof(double) * pn_pts[ipart]);

    for (int i = 0; i < pn_pts[ipart]; i++) {
      ppts_field3[ipart][i] = 0.;
    }

    for (int i = 0; i < n_located[ipart]; i++) {
      int pt_id = located[ipart][i] - 1;

      ppts_field3[ipart][pt_id] = tmp_ppts_field3[ipart][i];

      double exact = _eval_field(ppts_coord[ipart][3*i  ],
                                 ppts_coord[ipart][3*i+1],
                                 ppts_coord[ipart][3*i+2]);

      double err = PDM_ABS(exact - ppts_field3[ipart][pt_id]);
      max_err = PDM_MAX(max_err, err);
    }

    free(tmp_ppts_field3[ipart]);
  }
  free(tmp_ppts_field3);

  printf("[%d] max_err = %e\n", i_rank, max_err);

  /* Visu */
  if (post) {
    const char    *elt_field_name[] = {"field1", "field2"};
    double **elt_field_value[2] = {pelt_field1, pelt_field2};
    const char    *vtx_field_name[] = {"field3"};
    double **vtx_field_value[1] = {pvtx_field3};
    PDM_log_trace_array_int(pn_vtx, n_part_mesh, "pn_vtx : ");
    PDM_log_trace_array_int(pn_face, n_part_mesh, "pn_face : ");
    PDM_log_trace_array_int(pn_elt, n_part_mesh, "pn_elt : ");
    _dump_mesh("mesh_location_from_file__mesh",
               comm,
               n_part_mesh,
               pn_vtx,
               pvtx_coord,
               pvtx_ln_to_gn,
               pn_face,
               pface_vtx_idx,
               pface_vtx,
               pn_elt,
               pcell_face_idx,
               pcell_face,
               pelt_ln_to_gn,
               2,
               elt_field_name,
               elt_field_value,
               1,
               vtx_field_name,
               vtx_field_value);

    const char    *field_name[] = {"field1", "field2", "field3"};
    double **pts_field_value[3] = {ppts_field1, ppts_field2, ppts_field3};
    _dump_point_cloud("mesh_location_from_file__cloud",
                      comm,
                      n_part_cloud,
                      pn_pts,
                      ppts_coord,
                      ppts_ln_to_gn,
                      3,
                      field_name,
                      pts_field_value);
  }

  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    free(pelt_field1    [ipart]);
    free(pelt_field2    [ipart]);
    free(pvtx_field3     [ipart]);
    free(pelt_pts_field3[ipart]);
  }
  free(pelt_field1);
  free(pelt_field2);
  free(pvtx_field3);
  free(pelt_pts_field3);

  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    free(ppts_field1[ipart]);
    free(ppts_field2[ipart]);
    free(ppts_field3[ipart]);
  }
  free(ppts_field1);
  free(ppts_field2);
  free(ppts_field3);

  free(pelt_pts_idx       );
  free(pelt_pts_gnum      );
  free(pelt_pts_coord     );
  free(pelt_pts_uvw       );
  free(pelt_pts_weight_idx);
  free(pelt_pts_weight    );
  free(pelt_pts_dist2     );
  free(pelt_pts_proj_coord);
  free(pelt_vtx_idx       );
  free(pelt_vtx           );


  /* Free memory */
  PDM_part_to_part_free(ptp);
  PDM_mesh_location_free(mesh_loc);


  /* Free memory */
  if (use_edge) {
    for (int ipart = 0; ipart < n_part_mesh; ipart++) {
      free(pface_vtx[ipart]);
    }
  }
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pcell_face_idx);
  free(pcell_face);
  free(pelt_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pvtx_coord);
  free(pn_elt);
  free(pn_face);
  free(pn_vtx);

  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    free(ppts_coord   [ipart]);
    free(ppts_ln_to_gn[ipart]);
  }
  free(ppts_coord);
  free(ppts_ln_to_gn);
  free(pn_pts);

  PDM_multipart_free(mpart);


  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;
}
