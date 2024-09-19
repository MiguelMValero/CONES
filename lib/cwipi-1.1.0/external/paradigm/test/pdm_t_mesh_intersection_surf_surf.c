#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_mesh_intersection.h"
#include "pdm_multipart.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_part_to_part.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -nA     <level>  Number vtx in side of mesh A (default : 10).\n\n"
     "  -nA     <level>  Number vtx in side of mesh B (default : 10).\n\n"
     "  -n_part <level>  Number vtx in side of mesh B (default : 10).\n\n"
     "  -t               Element kind .\n\n"
     "  -h               This message.\n\n");
  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int                    argc,
 char                 **argv,
 PDM_g_num_t           *n_vtx_a,
 PDM_g_num_t           *n_vtx_b,
 int                   *n_part,
 PDM_Mesh_nodal_elt_t  *elt_type_a,
 PDM_Mesh_nodal_elt_t  *elt_type_b,
 int                   *nodal_a,
 int                   *nodal_b,
 int                   *verbose
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nA") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_vtx_a = atol(argv[i]);
        *n_vtx_a = (PDM_g_num_t) _n_vtx_a;
      }
    }
    else if (strcmp(argv[i], "-nB") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_vtx_b = atol(argv[i]);
        *n_vtx_b = (PDM_g_num_t) _n_vtx_b;
      }
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tA") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type_a = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-tB") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type_b = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-nodalA") == 0) {
      *nodal_a = 1;
    }
    else if (strcmp(argv[i], "-nodalB") == 0) {
      *nodal_b = 1;
    }
    else if (strcmp(argv[i], "-verbose") == 0) {
      *verbose = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


static void _rotate_coord(const double angle,
                                double *coord)
{
  double Rz[3][3] = {{cos(angle), -sin(angle), 0},
                     {sin(angle),  cos(angle), 0},
                     {0            ,  0            , 1}};

  double x = coord[0];
  double y = coord[1];
  double z = coord[2];

  for (int j = 0; j < 3; j++) {
    coord[j] = Rz[j][0]*x + Rz[j][1]*y + Rz[j][2]*z;
  }
}



static void
_add_depth
(
 double *coord
 )
{
  // double x = coord[0];
  double y = coord[1];
  double z = coord[2];

  coord[1] = 0.8*y - 0.6*z;
  coord[2] = 0.6*y + 0.8*z;

  // // angular sector
  // double t = PDM_PI * (0.5 + (x - 0.5 + 0.3*cos(3*y)) / 6.);
  // double r = 0.3 + 0.6 * y;
  // double scale = 0.2;
  // // x = r * cos(t);
  // // y = r * sin(t);
  // coord[0] = 0.07*y*(1-y)*cos(2*PDM_PI*x) + 0.05*sin(5*y);//;0.5 * scale * (cos(3*(x+y) + .2) + sin(5*y + .1));
  // coord[1] = r * cos(t);
  // coord[2] = r * sin(t);

}


static
void
_generate_surface_mesh
(
 const PDM_MPI_Comm           comm,
 const PDM_g_num_t            n_vtx_seg,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const int                    rotate,
 const double                 xmin,
 const double                 ymin,
 const double                 zmin,
 const double                 length,
 const PDM_split_dual_t       part_method,
 const int                    n_part,
       PDM_dmesh_nodal_t    **_dmn,
       PDM_multipart_t      **_mpart
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_dmesh_nodal_t *dmn = NULL;
  if (0) {
    PDM_sphere_surf_icosphere_gen_nodal(comm,
                                        n_vtx_seg,
                                        0, 0, 0,
                                        1,
                                        &dmn);
  }
  else {
    PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           length,
                                                           xmin,
                                                           ymin,
                                                           zmin,
                                                           elt_type,
                                                           1,
                                                           PDM_OWNERSHIP_USER);
    PDM_dcube_nodal_gen_build (dcube);
    dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);
    PDM_dcube_nodal_gen_free(dcube);
  }

  PDM_g_num_t* distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
  double* vtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

  // randomize
  if (0) {
    double noise = 0.2*length/(double) (n_vtx_seg - 1);
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      if (PDM_ABS(vtx_coord[3*i_vtx  ] - xmin         ) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx  ] - xmin - length) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx+1] - ymin         ) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx+1] - ymin - length) > 1.e-9) {
        srand(distrib_vtx[i_rank] + i_vtx);
        for (int i = 0; i < 2; i++) {
          vtx_coord[3*i_vtx+i] += noise*0.5*(2*rand()/(double) RAND_MAX - 1);
        }
      }
    }
  }

  if(rotate) {
    // Do something
    double pi = 4 * atan(1.);
    double angle = pi/5.;
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      _rotate_coord(angle, &vtx_coord[3*i_vtx]);
    }
  }

  if (1) {
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      _add_depth(&vtx_coord[3*i_vtx]);
    }
  }

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "sphere_surf_");
  }

  int n_domain = 1;
  // int n_part_domains = {n_part};
  int *n_part_domains = (int *) malloc(sizeof(int) * n_domain);
  n_part_domains[0] = n_part;

  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                n_part_domains,
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

  free(n_part_domains);


  *_mpart = mpart;
  *_dmn   = dmn;

}

static
void
_set_mesh
(
 PDM_mesh_intersection_t *mi,
 int                      i_mesh,
 PDM_multipart_t         *mpart,
 int                      n_part
)
{
  int i_rank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);

  PDM_mesh_intersection_n_part_set(mi, i_mesh, n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    double *vtx_coord  = NULL;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                 0,
                                                 i_part,
                                                 &vtx_coord,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *face_ln_to_gn  = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    &face_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *edge_ln_to_gn  = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_EDGE,
                                    &edge_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *vtx_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    &vtx_ln_to_gn,
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

    int *edge_vtx     = NULL;
    int *edge_vtx_idx = NULL;
    int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                     &edge_vtx_idx,
                                                     &edge_vtx,
                                                     PDM_OWNERSHIP_KEEP);

    // if (i_mesh == 1) {
    //   for (int i = 0; i < n_face; i++) {
    //     log_trace(PDM_FMT_G_NUM" : %d %d %d\n",
    //               face_ln_to_gn[i],
    //               i_rank,
    //               i_part,
    //               i);
    //   }
    // }

    PDM_mesh_intersection_part_set(mi,
                                   i_mesh,
                                   i_part,
                                   0, //n_cell
                                   n_face,
                                   n_edge,
                                   n_vtx,
                                   NULL,//cell_face_idx,
                                   NULL,//cell_face,
                                   face_edge_idx,
                                   face_edge,
                                   edge_vtx,
                                   NULL,//face_vtx_idx,
                                   NULL,//face_vtx,
                                   NULL,//cell_ln_to_gn,
                                   face_ln_to_gn,
                                   edge_ln_to_gn,
                                   vtx_ln_to_gn,
                                   vtx_coord);
  }

}

static
PDM_part_mesh_nodal_t *
_set_mesh_nodal
(
 PDM_mesh_intersection_t *mi,
 int                      i_mesh,
 PDM_multipart_t         *mpart
)
{
  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmn, PDM_OWNERSHIP_KEEP);

  if (0 == 1) {
    char filename[999];
    sprintf(filename, "check_pmn_%d", i_mesh);
    PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_SURFACIC, filename);
  }

  PDM_mesh_intersection_mesh_nodal_set(mi, i_mesh, pmn);

  return pmn;
}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t n_vtx_a   = 10;
  PDM_g_num_t n_vtx_b   = 10;
  int         nodal_a   = 0;
  int         nodal_b   = 0;
  PDM_Mesh_nodal_elt_t elt_type_a = PDM_MESH_NODAL_TRIA3;
  PDM_Mesh_nodal_elt_t elt_type_b = PDM_MESH_NODAL_TRIA3;

  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_HILBERT;

  int n_part = 1;
  int verbose = 0;

  _read_args(argc,
             argv,
             &n_vtx_a,
             &n_vtx_b,
             &n_part,
             &elt_type_a,
             &elt_type_b,
             &nodal_a,
             &nodal_b,
             &verbose);

  /*
   * Generate meshA
   */
  double length_a = 1.;
  int rotate_a = 0;
  PDM_dmesh_nodal_t     *dmn_surf_a   = NULL;
  PDM_multipart_t       *mpart_surf_a = NULL;
  _generate_surface_mesh (comm,
                          n_vtx_a,
                          elt_type_a,
                          rotate_a,
                          0.,
                          0.,
                          0.,
                          length_a,
                          part_method,
                          n_part,
                          &dmn_surf_a,
                          &mpart_surf_a);


  double length_b = 1.;
  int rotate_b = 1;
  PDM_dmesh_nodal_t     *dmn_surf_b   = NULL;
  PDM_multipart_t       *mpart_surf_b = NULL;
  _generate_surface_mesh (comm,
                          n_vtx_b,
                          elt_type_b,
                          rotate_b,
                          0.5,
                          0.5,
                          0.,
                          length_b,
                          part_method,
                          n_part,
                          &dmn_surf_b,
                          &mpart_surf_b);

  if(verbose) {
    PDM_dmesh_nodal_dump_vtk(dmn_surf_a,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "dmn_surf_a_");
    PDM_dmesh_nodal_dump_vtk(dmn_surf_b,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "dmn_surf_b_");
  }

  /*
   * Mesh_intersection
   */
  int dim_mesh_a = 2;
  int dim_mesh_b = 2;
  PDM_mesh_intersection_t* mi = PDM_mesh_intersection_create(PDM_MESH_INTERSECTION_KIND_WEIGHT,
                                                             dim_mesh_a,
                                                             dim_mesh_b,
                                                             1e-6,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);

  /*
   * Set mesh_a and mesh_b
   */
  PDM_part_mesh_nodal_t *pmn_a = NULL;
  PDM_part_mesh_nodal_t *pmn_b = NULL;
  if (nodal_a) {
    pmn_a = _set_mesh_nodal(mi, 0, mpart_surf_a);
  }
  else {
    _set_mesh(mi, 0, mpart_surf_a, n_part);
  }
  if (nodal_b) {
    pmn_b = _set_mesh_nodal(mi, 1, mpart_surf_b);
  }
  else {
    _set_mesh(mi, 1, mpart_surf_b, n_part);
  }

  PDM_mesh_intersection_compute(mi);


  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_intersection_part_to_part_get(mi,
                                         &ptp,
                                         PDM_OWNERSHIP_USER);

  if (ptp != NULL) {

    // Check ptp
    int  *n_ref_b = NULL;
    int **ref_b   = NULL;
    PDM_part_to_part_ref_lnum2_get(ptp,
                                   &n_ref_b,
                                   &ref_b);

    int         **pelt_b_elt_a_idx = NULL;
    PDM_g_num_t **pelt_b_elt_a     = NULL;
    PDM_part_to_part_gnum1_come_from_get(ptp,
                                         &pelt_b_elt_a_idx,
                                         &pelt_b_elt_a);

    if (verbose) {
      log_trace("FROM A USER POV\n");
    }
    int    **pelt_a_elt_b_n      = malloc(sizeof(int    *) * n_part);
    double **pelt_a_elt_b_volume = malloc(sizeof(double *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {
      int         *elt_a_elt_b_idx = NULL;
      PDM_g_num_t *elt_a_elt_b     = NULL;
      PDM_mesh_intersection_result_from_a_get(mi,
                                              ipart,
                                              &elt_a_elt_b_idx,
                                              &elt_a_elt_b,
                                              &pelt_a_elt_b_volume[ipart]);

      PDM_g_num_t *elt_a_ln_to_gn = NULL;
      int n_elt_a = PDM_multipart_part_ln_to_gn_get(mpart_surf_a,
                                                    0,
                                                    ipart,
                                                    PDM_MESH_ENTITY_FACE,
                                                    &elt_a_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);

      pelt_a_elt_b_n[ipart] = malloc(sizeof(int) * n_elt_a);
      for (int i = 0; i < n_elt_a; i++) {
        pelt_a_elt_b_n[ipart][i] = elt_a_elt_b_idx[i+1] - elt_a_elt_b_idx[i];

        if (verbose) {
          log_trace("elt_a "PDM_FMT_G_NUM" : ", elt_a_ln_to_gn[i]);
          for (int j = elt_a_elt_b_idx[i]; j < elt_a_elt_b_idx[i+1]; j++) {
            log_trace("("PDM_FMT_G_NUM", %f)  ", elt_a_elt_b[j], pelt_a_elt_b_volume[ipart][j]);
          }
          log_trace("\n");
        }
      }
    }

    double **pelt_b_elt_a_volume = NULL;
    int request = -1;
    PDM_part_to_part_iexch(ptp,
                           PDM_MPI_COMM_KIND_P2P,
                           PDM_STRIDE_CST_INTERLACED,
                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                           1,
                           sizeof(double),
                           NULL,
                           (const void  **) pelt_a_elt_b_volume,
                           NULL,
                           (      void ***) &pelt_b_elt_a_volume,
                           &request);
    PDM_part_to_part_iexch_wait(ptp, request);

    if (verbose) {
      log_trace("FROM B USER POV\n");
    }
    for (int ipart = 0; ipart < n_part; ipart++) {

      PDM_g_num_t *elt_b_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart_surf_b,
                                      0,
                                      ipart,
                                      PDM_MESH_ENTITY_FACE,
                                      &elt_b_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      double *volume = NULL;
      if (!nodal_b) { // en attendant le calcul des volumes en nodal...
        PDM_mesh_intersection_elt_volume_get(mi,
                                             1,
                                             ipart,
                                             &volume);
      }

      double *elt_b_elt_a_volume = NULL;
      PDM_mesh_intersection_result_from_b_get(mi,
                                              ipart,
                                              &elt_b_elt_a_volume);


      for (int i = 0; i < n_ref_b[ipart]; i++) {
        int faceB_id = ref_b[ipart][i] - 1;
        if (verbose) {
          log_trace("elt_b "PDM_FMT_G_NUM" : ", elt_b_ln_to_gn[faceB_id]);
          double sum_w = 0;
          for (int j = pelt_b_elt_a_idx[ipart][i]; j < pelt_b_elt_a_idx[ipart][i+1]; j++) {
            log_trace("("PDM_FMT_G_NUM", %f (%f))  ",
                      pelt_b_elt_a[ipart][j],
                      pelt_b_elt_a_volume[ipart][j],
                      elt_b_elt_a_volume[j]);
            sum_w += elt_b_elt_a_volume[j];
          }
          log_trace("\n");
          if (volume != NULL) {
            log_trace("  sum_w = %f/%f (%.1f%%)\n", sum_w, volume[faceB_id], 100*sum_w/volume[faceB_id]);
          }
        }
      }
    }

    PDM_MPI_Barrier(comm);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free(pelt_a_elt_b_n[ipart]);
    }
    free(pelt_a_elt_b_n     );
    free(pelt_a_elt_b_volume);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free(pelt_b_elt_a_volume[ipart]);
    }
    free(pelt_b_elt_a_volume);


    PDM_part_to_part_free(ptp);
  }
  PDM_mesh_intersection_free(mi);

  if (nodal_a) {
    PDM_part_mesh_nodal_free(pmn_a);
  }
  if (nodal_b) {
    PDM_part_mesh_nodal_free(pmn_b);
  }

  PDM_DMesh_nodal_free(dmn_surf_b);
  PDM_multipart_free(mpart_surf_b);

  PDM_DMesh_nodal_free(dmn_surf_a);
  PDM_multipart_free(mpart_surf_a);

  PDM_MPI_Barrier(comm);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize ();

  return 0;

}
