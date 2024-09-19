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
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_para_octree.h"

#include "pdm_multipart.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_vtk.h"
#include "pdm_reader_gamma.h"
#include "pdm_reader_stl.h"
#include "pdm_writer.h"

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
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
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
 int            argc,
 char         **argv,
 PDM_g_num_t   *nPts,
 double        *radius,
 int           *local,
 int           *rand,
 char         **filename,
 int           *visu,
 int           *points_in_leaf_max
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nPts = atol(argv[i]);
        *nPts = (PDM_g_num_t) _nPts;
      }
    }

    else if (strcmp(argv[i], "-radius") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *radius = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else if (strcmp(argv[i], "-rand") == 0) {
      *rand = 1;
    }

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *filename = argv[i];
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else if (strcmp(argv[i], "-points_in_leaf_max") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *points_in_leaf_max = atoi(argv[i]);
      }
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}



//https://stackoverflow.com/questions/5309471/getting-file-extension-in-c
static const char *
_get_filename_extension
(
 const char *filename
 )
{
  const char *dot = strrchr(filename, '.');
  if(!dot || dot == filename) return "";
  return dot + 1;
}

static void
_read_cloud_from_mesh
(
 const PDM_MPI_Comm   comm,
 const char          *filename,
       int           *n_pts,
       double       **pts_coord,
       PDM_g_num_t  **pts_ln_to_gn,
 const int            visu
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;


  const char *file_extension = _get_filename_extension(filename);

  PDM_dmesh_nodal_t *dmn = NULL;

  if (strcmp(file_extension, "vtk") == 0) {
    int            n_vtx_field      = 0;
    char         **vtx_field_name   = NULL;
    PDM_data_t    *vtx_field_type   = NULL;
    int           *vtx_field_stride = NULL;
    void         **vtx_field_value  = NULL;
    int            n_elt_field      = 0;
    char         **elt_field_name   = NULL;
    PDM_data_t    *elt_field_type   = NULL;
    int           *elt_field_stride = NULL;
    void         **elt_field_value  = NULL;
    dmn = PDM_vtk_read_to_dmesh_nodal(comm,
                                      filename,
                                      &n_vtx_field,
                                      &vtx_field_name,
                                      &vtx_field_type,
                                      &vtx_field_stride,
                                      &vtx_field_value,
                                      &n_elt_field,
                                      &elt_field_name,
                                      &elt_field_type,
                                      &elt_field_stride,
                                      &elt_field_value);

    if (n_vtx_field > 0) {
      for (int i = 0; i < n_vtx_field; i++) {
        free(vtx_field_name [i]);
        free(vtx_field_value[i]);
      }
      free(vtx_field_name );
      free(vtx_field_type );
      free(vtx_field_stride);
      free(vtx_field_value);
    }

    if (n_elt_field > 0) {
      for (int i = 0; i < n_elt_field; i++) {
        free(elt_field_name [i]);
        // free(elt_field_value[i]);
      }
      free(elt_field_name );
      free(elt_field_type );
      free(elt_field_stride);
      free(elt_field_value);
    }
  }
  else if (strcmp(file_extension, "mesh") == 0) {
    dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                       filename,
                                       0,
                                       0);
  }
  else if (strcmp(file_extension, "mesh") == 0) {
    dmn = PDM_reader_stl_dmesh_nodal(comm,
                                     filename);
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "File format '%s' is not supported\n", file_extension);
  }


  int n_domain = 1;
  int n_part = 1;
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


  double      *vtx_coord    = NULL;
  PDM_g_num_t *vtx_ln_to_gn = NULL;
  int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                               0,
                                               0,
                                               &vtx_coord,
                                               PDM_OWNERSHIP_KEEP); // USER is broken

  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  0,
                                  PDM_MESH_ENTITY_VTX,
                                  &vtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP); // USER is broken

  *n_pts        = n_vtx;
  // *pts_coord    = vtx_coord;
  // *pts_ln_to_gn = vtx_ln_to_gn;

  *pts_coord = malloc(sizeof(double) * n_vtx * 3);
  memcpy(*pts_coord, vtx_coord, sizeof(double) * n_vtx * 3);

  *pts_ln_to_gn = malloc(sizeof(PDM_g_num_t) * n_vtx);
  memcpy(*pts_ln_to_gn, vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);

  if (visu) {
    PDM_part_mesh_nodal_t *pmn = NULL;
    PDM_multipart_get_part_mesh_nodal(mpart,
                                      0,
                                      &pmn,
                                      PDM_OWNERSHIP_KEEP);

    PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                          PDM_WRITER_FMT_BIN,
                                          PDM_WRITER_TOPO_CST,
                                          PDM_WRITER_OFF,
                                          "para_octree_mesh",
                                          "para_octree_mesh",
                                          comm,
                                          PDM_IO_KIND_MPI_SIMPLE,
                                          1.,
                                          NULL);

    int id_geom = PDM_writer_geom_create_from_mesh_nodal(wrt,
                                                         "para_octree_mesh",
                                                         pmn);

    int id_var_rank = PDM_writer_var_create(wrt,
                                            PDM_WRITER_OFF,
                                            PDM_WRITER_VAR_SCALAR,
                                            PDM_WRITER_VAR_ELEMENTS,
                                            "i_rank");

    PDM_writer_step_beg(wrt, 0.);

    PDM_writer_geom_write(wrt,
                          id_geom);

    PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_SURFACIC; // !!!
    int n_elt = PDM_part_mesh_nodal_n_elmts_get(pmn,
                                                geom_kind,
                                                0);

    PDM_real_t *val_rank = malloc(sizeof(PDM_real_t) * n_elt);
    for (int i = 0; i < n_elt; i++) {
      val_rank[i] = (PDM_real_t) i_rank;
    }

    PDM_writer_var_set(wrt,
                       id_var_rank,
                       id_geom,
                       0,
                       val_rank);

    PDM_writer_var_write(wrt,
                         id_var_rank);
    PDM_writer_var_free(wrt,
                        id_var_rank);

    PDM_writer_step_end(wrt);

    PDM_writer_free(wrt);

    free(val_rank);


    PDM_part_mesh_nodal_free(pmn);
  }


  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);
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
int   argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int numProcs;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  PDM_g_num_t  nPts               = 10;
  double       radius             = 10.;
  int          local              = 0;
  int          rand               = 0;
  char        *filename           = NULL;
  int          visu               = 0;
  int          points_in_leaf_max = 1;

  _read_args(argc,
             argv,
             &nPts,
             &radius,
             &local,
             &rand,
             &filename,
             &visu,
             &points_in_leaf_max);

  /* Initialize random */

  if (rand) {
    srand(time(NULL));
  }
  else {
    srand(i_rank);
  }

  int _n_pts_l;
  double      *coords = NULL;
  PDM_g_num_t *gnum   = NULL;
  if (filename == NULL) {
    /* Random point cloud */
    PDM_point_cloud_gen_random (PDM_MPI_COMM_WORLD,
                                0, // seed
                                0, // geometric_g_num
                                nPts,
                                -radius, -radius, -radius,
                                radius, radius, radius,
                                &_n_pts_l,
                                &coords,
                                &gnum);
  }
  else {
    _read_cloud_from_mesh(PDM_MPI_COMM_WORLD,
                          filename,
                          &_n_pts_l,
                          &coords,
                          &gnum,
                          visu);
  }

  /* Parallel octree */

  const int n_point_cloud = 1;
  const int depth_max = 31;
  // const int points_in_leaf_max = 1;

  const int build_leaf_neighbours = 1;
  PDM_para_octree_t *octree = PDM_para_octree_create (n_point_cloud,
                                                      depth_max,
                                                      points_in_leaf_max,
                                                      build_leaf_neighbours,
                                                      PDM_MPI_COMM_WORLD);

  PDM_para_octree_point_cloud_set (octree, 0, _n_pts_l, coords, gnum);

  PDM_para_octree_build (octree, NULL);

  // PDM_para_octree_dump (octree);

  PDM_para_octree_dump_times (octree);

  if (visu) {
    PDM_para_octree_export_vtk(octree, "para_octree");
  }

  PDM_para_octree_free (octree);

  /* Free */

  free (coords);
  free (gnum);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
