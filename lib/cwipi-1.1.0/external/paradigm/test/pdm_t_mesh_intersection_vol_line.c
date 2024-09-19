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
 int                   *post,
 PDM_Mesh_nodal_elt_t  *elt_type
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
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
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


static
void
_generate_volume_mesh
(
 const PDM_MPI_Comm           comm,
 const PDM_g_num_t            n_vtx_seg,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const int                    rotate,
 const double                 xmin,
 const double                 ymin,
 const double                 zmin,
 const double                 lenght,
 const PDM_split_dual_t       part_method,
 const int                    n_part,
       PDM_dmesh_nodal_t    **_dmn,
       PDM_multipart_t      **_mpart
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         lenght,
                                                         xmin,
                                                         ymin,
                                                         zmin,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_USER);
  PDM_dcube_nodal_gen_build (dcube);
  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);
  PDM_dcube_nodal_gen_free(dcube);

  if(rotate) {
    // Do something
    double pi = 4 * atan(1.);
    double angle = pi/5.;
    PDM_g_num_t* distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
    double* vtx_coord = PDM_DMesh_nodal_vtx_get(dmn);
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      _rotate_coord(angle, &vtx_coord[3*i_vtx]);
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
  PDM_mesh_intersection_n_part_set(mi, i_mesh, n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int *face_edge_idx = NULL;
    int *face_edge     = NULL;
    int *edge_vtx_idx  = NULL;
    int *edge_vtx      = NULL;

    double *vtx_coord = NULL;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                 0,
                                                 i_part,
                                                 &vtx_coord,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *cell_ln_to_gn = NULL;
    int n_cell = PDM_multipart_part_ln_to_gn_get(mpart,
                                                 0,
                                                 i_part,
                                                 PDM_MESH_ENTITY_CELL,
                                                 &cell_ln_to_gn,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    &face_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *edge_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_EDGE,
                                    &edge_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);
    PDM_g_num_t *vtx_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    &vtx_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    int *cell_face_idx = NULL;
    int *cell_face     = NULL;
    n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                 0,
                                                 i_part,
                                                 PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                 &cell_face_idx,
                                                 &cell_face,
                                                 PDM_OWNERSHIP_KEEP);
    int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                     &face_edge_idx,
                                                     &face_edge,
                                                     PDM_OWNERSHIP_KEEP);
    int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                     &edge_vtx_idx,
                                                     &edge_vtx,
                                                     PDM_OWNERSHIP_KEEP);

    PDM_mesh_intersection_part_set(mi,
                                   i_mesh,
                                   i_part,
                                   n_cell,
                                   n_face,
                                   n_edge,
                                   n_vtx,
                                   cell_face_idx,
                                   cell_face,
                                   face_edge_idx,
                                   face_edge,
                                   edge_vtx,
                                   NULL, // face_vtx_idx,
                                   NULL, // face_vtx,
                                   cell_ln_to_gn,
                                   face_ln_to_gn,
                                   edge_ln_to_gn,
                                   vtx_ln_to_gn,
                                   vtx_coord);
  }

}

static
void
_generate_ray
(
       PDM_MPI_Comm           comm,
 const PDM_g_num_t            gn_ray,
       double                 min_size,
       double                 max_size,
       double                 x_min,
       double                 y_min,
       double                 z_min,
       double                 x_max,
       double                 y_max,
       double                 z_max,
       double               **vtx_coord_out,
       PDM_g_num_t          **vtx_ln_to_gn_out,
       PDM_g_num_t          **edge_ln_to_gn_out,
       int                  **edge_vtx_out,
       int                   *pn_ray
)
{
  PDM_UNUSED(min_size);
  PDM_UNUSED(max_size);

  double origin[3] = {x_min, y_min, z_min};
  double length[3] = {x_max - x_min, y_max - y_min, z_max - z_min};

  int seed = 0;
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_g_num_t *distrib_ray = PDM_compute_uniform_entity_distribution (comm,
                                                                      gn_ray);

  int dn_ray = distrib_ray[i_rank+1] - distrib_ray[i_rank];

  int dn_vtx = 2 * dn_ray;

  PDM_g_num_t *vtx_ln_to_gn  = malloc(    dn_vtx * sizeof(PDM_g_num_t));
  PDM_g_num_t *edge_ln_to_gn = malloc(    dn_ray * sizeof(PDM_g_num_t));
  double      *vtx_coord     = malloc(3 * dn_vtx * sizeof(double     ));
  int         *edge_vtx      = malloc(3 * dn_ray * sizeof(int        ));

  int i_vtx = 0;
  for(int i = 0; i < dn_ray; ++i) {
    unsigned int _seed = (unsigned int) (distrib_ray[i_rank] + i) + 1;
    srand(_seed + seed);

    for (int j = 0; j < 3; j++) {
      vtx_coord[6*i + j    ] = origin[j] + length[j] * ((double) rand() / (double) RAND_MAX);
      vtx_coord[6*i + j + 3] = origin[j] + length[j] * ((double) rand() / (double) RAND_MAX);
    }

    edge_ln_to_gn[  i  ] =     (distrib_ray[i_rank] + i) + 1;
    vtx_ln_to_gn [2*i  ] = 2 * (distrib_ray[i_rank] + i) + 1;
    vtx_ln_to_gn [2*i+1] = 2 * (distrib_ray[i_rank] + i) + 2;

    edge_vtx[2*i  ] = i_vtx+1;
    edge_vtx[2*i+1] = i_vtx+2;
    i_vtx += 2;

  }

  *vtx_coord_out     = vtx_coord;
  *vtx_ln_to_gn_out  = vtx_ln_to_gn;
  *edge_ln_to_gn_out = edge_ln_to_gn;
  *edge_vtx_out      = edge_vtx;
  *pn_ray            = dn_ray;

  free(distrib_ray);
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
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t n_vtx_a   = 10;
  PDM_g_num_t n_g_ray   = 10;
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_HEXA8;

  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;

  int n_part = 1;
  int post   = 0;

  _read_args(argc,
             argv,
             &n_vtx_a,
             &n_g_ray,
             &n_part,
             &post,
             &elt_type);

  /*
   * Generate meshA
   */
  double lenght_a = 1.;
  int rotate_a = 0;
  PDM_dmesh_nodal_t     *dmn_vol_a   = NULL;
  PDM_multipart_t       *mpart_vol_a = NULL;
  _generate_volume_mesh (comm,
                         n_vtx_a,
                         elt_type,
                         rotate_a,
                         0.,
                         0.,
                         0.,
                         lenght_a,
                         part_method,
                         n_part,
                         &dmn_vol_a,
                         &mpart_vol_a);

  // Generate ray
  int pn_ray = 0;
  double      *pvtx_coord     = NULL;
  int         *pedge_vtx      = NULL;
  PDM_g_num_t *pedge_ln_to_gn = NULL;
  PDM_g_num_t *pvtx_ln_to_gn = NULL;
  int min_size = 0.25;
  int max_size = 0.5;
  double x_min = 0.;
  double y_min = 0.;
  double z_min = 0.;
  double x_max = 1.;
  double y_max = 1.;
  double z_max = 1.;
  _generate_ray(comm,
                n_g_ray,
                min_size,
                max_size,
                x_min,
                y_min,
                z_min,
                x_max,
                y_max,
                z_max,
                &pvtx_coord,
                &pvtx_ln_to_gn,
                &pedge_ln_to_gn,
                &pedge_vtx,
                &pn_ray);

  if(post) {
    char filename[999];
    sprintf(filename, "out_ray_%i.vtk", i_rank);
    PDM_vtk_write_std_elements(filename,
                               2 * pn_ray,
                               pvtx_coord,
                               pvtx_ln_to_gn,
                               PDM_MESH_NODAL_BAR2,
                               pn_ray,
                               pedge_vtx,
                               pedge_ln_to_gn,
                               0,
                               NULL,
                               NULL);

    PDM_dmesh_nodal_dump_vtk(dmn_vol_a,
                             PDM_GEOMETRY_KIND_VOLUMIC,
                             "dmn_vol_a_");
  }

  /*
   * Mesh_intersection
   */
  int dim_mesh_a = 3;
  int dim_mesh_b = 1;
  PDM_mesh_intersection_t* mi = PDM_mesh_intersection_create(PDM_MESH_INTERSECTION_KIND_WEIGHT,
                                                             dim_mesh_a,
                                                             dim_mesh_b,
                                                             1e-6,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);

  /*
   * Set mesh_a and mesh_b
   */
  _set_mesh(mi, 0, mpart_vol_a , n_part);

  /*
   * Set line mesh
   */
  PDM_mesh_intersection_n_part_set(mi, 1, 1);

  PDM_mesh_intersection_part_set(mi,
                                 1, // i_mesh
                                 0,
                                 0,
                                 0,
                                 pn_ray,
                                 2 * pn_ray,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 pedge_vtx,
                                 NULL, // face_vtx_idx,
                                 NULL, // face_vtx,
                                 NULL,
                                 NULL,
                                 pedge_ln_to_gn,
                                 pvtx_ln_to_gn,
                                 pvtx_coord);


  // PDM_log_trace_array_long(pvtx_ln_to_gn, 2 * pn_ray, "pvtx_ln_to_gn :");
  // PDM_log_trace_array_long(pedge_ln_to_gn, pn_ray, "pedge_ln_to_gn :");

  PDM_mesh_intersection_compute(mi);

  PDM_mesh_intersection_free(mi);


  PDM_DMesh_nodal_free(dmn_vol_a);
  PDM_multipart_free(mpart_vol_a);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  free(pvtx_coord);
  free(pedge_vtx);
  free(pedge_ln_to_gn);
  free(pvtx_ln_to_gn);
  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize ();

  return 0;

}
