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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_priv.h"
#include "pdm_multipart.h"

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
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
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
_read_args
(
 int                    argc,
 char                 **argv,
 PDM_g_num_t           *n_vtx_seg,
 double                *length,
 int                   *n_part,
 int                   *post,
 PDM_Mesh_nodal_elt_t  *elt_type,
 int                   *order,
 int                   *method
 )
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
    else if (strcmp(argv[i], "-type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-order") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
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
  double amplitude = 0.07;
  double frequency = 8.;

  for (int i = 0; i < n_vtx; i++) {
    double x = (vtx_coord[3*i    ] - 0.5) / length;
    double y = (vtx_coord[3*i + 1] - 0.5) / length;
    double z = (vtx_coord[3*i + 2] - 0.5) / length;

    vtx_coord[3*i    ] += amplitude*length*cos(frequency*y);
    vtx_coord[3*i + 1] += amplitude*length*cos(frequency*z);
    vtx_coord[3*i + 2] += amplitude*length*cos(frequency*x);
  }

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

  PDM_g_num_t          n_vtx_seg   = 10;
  double               length      = 1.;
  int                  n_part      = 1;
  int                  post        = 0;
  PDM_Mesh_nodal_elt_t elt_type    = PDM_MESH_NODAL_HEXA8;
  int                  order       = 1;
  PDM_split_dual_t     part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
             &elt_type,
             &order,
             (int *) &part_method);

  /*
   *  Init
   */

  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        elt_type,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);
  // PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
  //                                                        n_vtx_seg,
  //                                                        n_vtx_seg,
  //                                                        n_vtx_seg,
  //                                                        length,
  //                                                        0.,
  //                                                        0.,
  //                                                        0.,
  //                                                        PDM_MESH_NODAL_TRIAHO,//HEXAHO,
  //                                                        2,
  //                                                        PDM_OWNERSHIP_KEEP);
  // // PDM_dcube_nodal_gen_ordering_set(dcube, "PDM_HO_ORDERING_VTK");
  // PDM_dcube_nodal_gen_ordering_set(dcube, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube);



  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
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

  PDM_part_mesh_nodal_t* pmsh_nodal = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmsh_nodal, PDM_OWNERSHIP_KEEP); // Ownership keep is mandatory in C

  PDM_multipart_free(mpart);


  if (post) {
    int elt_dim = PDM_Mesh_nodal_elt_dim_get(elt_type);
    if (elt_dim > 2) {
      PDM_part_mesh_nodal_dump_vtk(pmsh_nodal, PDM_GEOMETRY_KIND_VOLUMIC , "volumic_ho_" );
    }
    if (elt_dim > 1) {
      PDM_part_mesh_nodal_dump_vtk(pmsh_nodal, PDM_GEOMETRY_KIND_SURFACIC, "surfacic_ho_");
    }
    PDM_part_mesh_nodal_dump_vtk(pmsh_nodal, PDM_GEOMETRY_KIND_RIDGE   , "ridge_ho_"   );
  }
  PDM_part_mesh_nodal_free(pmsh_nodal);

  gettimeofday(&t_elaps_debut, NULL);
  PDM_dcube_nodal_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}


