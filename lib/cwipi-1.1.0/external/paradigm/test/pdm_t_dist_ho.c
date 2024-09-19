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
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"

#include "pdm_vtk.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"



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
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *n_vtx_seg,
           PDM_g_num_t           *n_g_pts_clouds,
           double                *length,
           int                   *n_part,
           int                   *post,
           int                   *order,
           PDM_Mesh_nodal_elt_t  *elt_type,
           int                   *method)
{
  int i = 1;

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
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_g_pts_clouds = atol(argv[i]);
        *n_g_pts_clouds = (PDM_g_num_t) _n_g_pts_clouds;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
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
    else if (strcmp(argv[i], "-order") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_deformation
(
 const double       length,
 const PDM_g_num_t  n_vtx_seg,
 const int          n_vtx,
 double            *vtx_coord
 )
{
  PDM_UNUSED(n_vtx_seg);
  PDM_UNUSED(length);
  //double amplitude = 0.1;//0.07;
  //double frequency = 4.;

  for (int i = 0; i < n_vtx; i++) {
    // double x = (vtx_coord[3*i    ] - 0.5) / length;
    // double y = (vtx_coord[3*i + 1] - 0.5) / length;
    // double z = (vtx_coord[3*i + 2] - 0.5) / length;

    // vtx_coord[3*i    ] += amplitude*length*cos(frequency*y);
    // vtx_coord[3*i + 1] += amplitude*length*cos(frequency*z);
    // vtx_coord[3*i + 2] += amplitude*length*cos(frequency*x);
    vtx_coord[3*i+2] = vtx_coord[3*i]*vtx_coord[3*i];
  }

}


static
void
_generate_surface_mesh
(
 const PDM_MPI_Comm           comm,
 const int                    order,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const PDM_g_num_t            n_vtx_seg,
 const double                 zero_x,
 const double                 zero_y,
 const double                 zero_z,
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

  /* First: generate a dcube nodal */
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        zero_x,
                                                        zero_y,
                                                        zero_z,
                                                        elt_type,
                                                        order,
                                                        PDM_OWNERSHIP_USER);

  PDM_dcube_nodal_gen_ordering_set(dcube, "PDM_HO_ORDERING_CGNS");

  PDM_dcube_nodal_gen_build(dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  if (i_rank == 0) {
    fflush(stdout);
  }


  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  // double amplitude = 0.1;//0.07;
  // double frequence = 4.;
  _deformation(length,
               n_vtx_seg,
               dn_vtx,
               dvtx_coord);


  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "surface_mesh_");
  }

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

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);

  PDM_dcube_nodal_gen_free(dcube);

  *_mpart = mpart;
  *_dmn   = dmn;
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

  PDM_g_num_t           n_vtx_seg      = 10;
  PDM_g_num_t           n_g_pts_clouds = 10;
  double                length         = 1.;
  int                   n_part         = 1;
  int                   post           = 0;
  int                   order          = 1;
  PDM_Mesh_nodal_elt_t  elt_type       = PDM_MESH_NODAL_TRIA3;

  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;

  setenv("PDM_DIST_CLOUD_SURF_OPTIM", "1", 1);

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &n_g_pts_clouds,
             &length,
             &n_part,
             &post,
             &order,
             &elt_type,
     (int *) &part_method);

  //double radius = length;//2*length;



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
   *  Generate cloud
   */
  if (i_rank == 0) {
    printf("-- Generate point cloud\n");
    fflush(stdout);
  }
  int          n_pts_clouds;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  PDM_point_cloud_gen_random(comm,
                             0, // seed
                             0, // geometric_g_num
                             n_g_pts_clouds,
                             -0.5*length, -0.5*length, -0.5*length,// -0.1*radius, -0.1*radius, -0.5*radius,
                             0.5*length, 0.5*length, 0.5*length,// 1.1*radius, 1.1*radius, 0.5*radius,
                             &n_pts_clouds,
                             &pts_coord,
                             &pts_g_num);


  /*
   *  Generate surface mesh
   */
  if (i_rank == 0) {
    printf("-- Generate surface mesh\n");
    fflush(stdout);
  }

  PDM_dmesh_nodal_t *dmn_surf   = NULL;
  PDM_multipart_t   *mpart_surf = NULL;
  _generate_surface_mesh(comm,
                         order,
                         elt_type,
                         n_vtx_seg,
                         -0.5*length,
                         -0.5*length,
                         -0.5*length,
                         length,
                         part_method,
                         n_part,
                         &dmn_surf,
                         &mpart_surf);


  /*
   *  Create dist_cloud_surf object
   */

  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t *dist = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_NODAL_SHARED,
                                                           n_point_cloud,
                                                           comm,
                                                           PDM_OWNERSHIP_KEEP);

  /* Set point cloud */
  int n_part_cloud  = 1;
  int i_point_cloud = 0;
  int i_part_cloud  = 0;
  PDM_dist_cloud_surf_n_part_cloud_set(dist, i_part_cloud, n_part_cloud);

  PDM_dist_cloud_surf_cloud_set(dist,
                                i_point_cloud,
                                i_part_cloud,
                                n_pts_clouds,
                                pts_coord,
                                pts_g_num);

  /* Set surface mesh */
  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart_surf, 0, &pmn, PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_nodal_mesh_set(dist, pmn);

  if (post) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "dist_ho_mesh_");
  }


  /* Compute distance */
  PDM_dist_cloud_surf_compute(dist);

  /* Pretty visu */
  if (post) {
    char filename[999];
    sprintf(filename, "dist_ho_cloud_%d.vtk", i_rank);

    double      *distance         = NULL;
    double      *projected        = NULL;
    PDM_g_num_t *closest_elt_gnum = NULL;
    PDM_dist_cloud_surf_get(dist,
                            0,
                            0,
                            &distance,
                            &projected,
                            &closest_elt_gnum);

    double *coord = malloc(sizeof(double) * n_pts_clouds * 2 * 3);
    memcpy(coord,                  pts_coord, sizeof(double) * n_pts_clouds * 3);
    memcpy(coord + 3*n_pts_clouds, projected, sizeof(double) * n_pts_clouds * 3);

    int *connec = malloc(sizeof(int) * n_pts_clouds * 2);
    for (int i = 0; i < n_pts_clouds; i++) {
      connec[2*i  ] = i+1;
      connec[2*i+1] = i+1 + n_pts_clouds;
    }

    PDM_g_num_t *g_num = malloc(sizeof(PDM_g_num_t) * n_pts_clouds * 2);
    memcpy(g_num,                pts_g_num,        sizeof(PDM_g_num_t) * n_pts_clouds);
    memcpy(g_num + n_pts_clouds, closest_elt_gnum, sizeof(PDM_g_num_t) * n_pts_clouds);

    double *_closest_elt_gnum = malloc(sizeof(double) * n_pts_clouds);
    for (int i = 0; i < n_pts_clouds; i++) {
      _closest_elt_gnum[i] = (double) closest_elt_gnum[i];
    }

    const char   *field_name []  = {"distance2", "closest_face"};
    const double *field_value[2] = {distance, _closest_elt_gnum};

    PDM_vtk_write_std_elements_double(filename,
                                      n_pts_clouds * 2,
                                      coord,
                                      g_num,
                                      PDM_MESH_NODAL_BAR2,
                                      n_pts_clouds,
                                      connec,
                                      pts_g_num,
                                      2,
                                      field_name,
                                      field_value);
    free(_closest_elt_gnum);
    free(coord);
    free(g_num);
    free(connec);
  }


  /*
   *  Free memory
   */
  PDM_dist_cloud_surf_free(dist);
  PDM_part_mesh_nodal_free(pmn);

  PDM_multipart_free(mpart_surf);
  PDM_DMesh_nodal_free(dmn_surf);

  free(pts_coord);
  free(pts_g_num);


  PDM_MPI_Barrier(comm);
  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
