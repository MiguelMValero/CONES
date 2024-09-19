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
#include "pdm_sphere_vol_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_priv.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_ho_ordering.h"



/*============================================================================
 * Type definitions
 *============================================================================*/

#define _MIN(a,b) ((a) < (b) ? (a) : (b))
#define _MAX(a,b) ((a) > (b) ? (a) : (b))


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
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *nx,
           PDM_g_num_t           *ny,
           PDM_g_num_t           *nz,
           double                *center_x,
           double                *center_y,
           double                *center_z,
           int                   *order,
           PDM_Mesh_nodal_elt_t  *t_elt,
           double                *radius)
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
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-o") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-r") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *radius = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-cx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *center_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-cy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *center_y = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-cz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *center_z = atof(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_dmesh_nodal_dump_vtk
(
       PDM_dmesh_nodal_t   *dmn,
       int                  order,
       PDM_geometry_kind_t  geom_kind,
 const char                *filename_patter
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  int* sections_id = PDM_DMesh_nodal_sections_id_get(dmn, geom_kind);
  int n_section    = PDM_DMesh_nodal_n_section_get(dmn, geom_kind);


  const char *field_name = "group";
  int n_field = 0;
  _pdm_dmesh_nodal_elts_t *dmne = NULL;
  int *delt_group_idx = NULL;
  int *delt_group     = NULL;
  double **field = NULL;
  if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
    dmne = dmn->ridge;
  } else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {// && dmn->mesh_dimension == 3) {
    dmne = dmn->surfacic;
  }

  PDM_g_num_t *distrib_elt = NULL;
  if (dmne != NULL) {
    distrib_elt = PDM_compute_uniform_entity_distribution(dmn->comm,
                                                          dmne->n_g_elmts);
    // PDM_log_trace_array_long(distrib_elt, n_rank+1, "distrib_elt : ");

    PDM_dgroup_entity_transpose(dmne->n_group_elmt,
                                dmne->dgroup_elmt_idx,
                                dmne->dgroup_elmt,
                (PDM_g_num_t *) distrib_elt,
                                &delt_group_idx,
                                &delt_group,
                                dmn->comm);
  }

  PDM_g_num_t shift = 0;
  for(int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = sections_id[i_section];
    const PDM_g_num_t    *delmt_distribution = PDM_DMesh_nodal_distrib_section_get(dmn, geom_kind, id_section);
    int                   n_elt              = PDM_DMesh_nodal_section_n_elt_get  (dmn, geom_kind, id_section);
    PDM_g_num_t          *dconnec            = PDM_DMesh_nodal_section_std_get    (dmn, geom_kind, id_section);
    PDM_Mesh_nodal_elt_t  t_elt              = PDM_DMesh_nodal_section_type_get   (dmn, geom_kind, id_section);

    int         *dconnec_idx    = (int         * ) malloc( (n_elt+1) * sizeof(int        ));
    PDM_g_num_t *delmt_ln_to_gn = (PDM_g_num_t * ) malloc( (n_elt  ) * sizeof(PDM_g_num_t));

    int strid = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    dconnec_idx[0] = 0;
    for(int i = 0; i < n_elt; ++i) {
      dconnec_idx[i+1] = dconnec_idx[i] + strid;
      delmt_ln_to_gn[i] = delmt_distribution[i_rank] + i + 1;
    }

    // log_trace("section %d (%d) :\n", i_section, id_section);
    // PDM_log_trace_array_long(delmt_ln_to_gn, n_elt, "  delmt_ln_to_gn : ");

    PDM_g_num_t *pvtx_ln_to_gn;
    int         *pcell_vtx_idx;
    int         *pcell_vtx;
    int          pn_vtx;
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(dmn->comm,
                                                             delmt_distribution,
                                                             dconnec_idx,
                                                             dconnec,
                                                             n_elt,
                                    (const PDM_g_num_t *)    delmt_ln_to_gn,
                                                            &pn_vtx,
                                                            &pvtx_ln_to_gn,
                                                            &pcell_vtx_idx,
                                                            &pcell_vtx);

    /*
     * Coordinates
     */
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    // int          dn_vtx   = PDM_DMesh_nodal_n_vtx_get(dln->dmesh_nodal_in);
    // assert(dn_vtx == (vtx_distrib[i_rank+1]-vtx_distrib[i_rank]));
    double** tmp_pvtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(dmn->comm,
                                          1,
                                          vtx_distrib,
                                          dvtx_coord,
                                          &pn_vtx,
                   (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                          &tmp_pvtx_coord);

    double* pvtx_coord_out = tmp_pvtx_coord[0];

    /*
     * Groups
     */
    if (dmne != NULL) {
      for(int i = 0; i < n_elt; ++i) {
        delmt_ln_to_gn[i] += shift;
      }

      int **tmp_elt_group_idx = NULL;
      int **tmp_elt_group     = NULL;
      PDM_part_dentity_group_to_pentity_group(dmn->comm,
                                              1,
                                              distrib_elt,
                                              delt_group_idx,
                                              delt_group,
                                              &n_elt,
                      (const PDM_g_num_t **)  &delmt_ln_to_gn,
                                              &tmp_elt_group_idx,
                                              &tmp_elt_group);
      int *pelt_group_idx = tmp_elt_group_idx[0];
      int *pelt_group     = tmp_elt_group    [0];
      // PDM_log_trace_connectivity_int(pelt_group_idx, pelt_group, n_elt, "pelt_group : ");
      free (tmp_elt_group_idx);
      free (tmp_elt_group);
      // PDM_log_trace_array_long(delmt_ln_to_gn, n_elt, "  delmt_ln_to_gn (shifted) : ");

      n_field = 1;
      field = malloc (sizeof(double *) * n_field);
      field[0] = malloc (sizeof(double) * n_elt);
      for (int i = 0; i < n_elt; i++) {
        assert (pelt_group_idx[i+1] == pelt_group_idx[i] + 1);
        field[0][i] = (double) pelt_group[i];
      }
      free (pelt_group);
      free (pelt_group_idx);
    }

    /*
     *  Dump
     */
    char filename[999];
    sprintf(filename, "%s_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    if (order == 1) {
      PDM_vtk_write_std_elements_double(filename,
                                        pn_vtx,
                                        pvtx_coord_out,
                                        pvtx_ln_to_gn,
                                        t_elt,
                                        n_elt,
                                        pcell_vtx,
                                        delmt_ln_to_gn,
                                        n_field,
                                        (const char   **) &field_name,
                                        (const double **) field);
    } else {
      PDM_vtk_write_std_elements_ho(filename,
                                    order,
                                    pn_vtx,
                                    pvtx_coord_out,
                                    pvtx_ln_to_gn,
                                    t_elt,
                                    n_elt,
                                    pcell_vtx,
                                    delmt_ln_to_gn,
                                    n_field,
                                    (const char   **) &field_name,
                                    (const double **) field);
    }
    free(tmp_pvtx_coord);
    free(pvtx_ln_to_gn);
    free(pcell_vtx_idx);
    free(pcell_vtx);

    free(dconnec_idx);
    free(delmt_ln_to_gn);

    free(pvtx_coord_out);

    shift += delmt_distribution[n_rank];

    if (dmne != NULL) {
      free (field[0]);
      free (field);
    }
  }

  if (dmne != NULL) {
    free (delt_group_idx);
    free (delt_group);
    free (distrib_elt);
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
  PDM_g_num_t          nx       = 10;
  PDM_g_num_t          ny       = 10;
  PDM_g_num_t          nz       = 10;
  double               center_x = 0.;
  double               center_y = 0.;
  double               center_z = 0.;
  int                  order    = 1;
  double               radius   = 1.;
  PDM_Mesh_nodal_elt_t t_elt    = PDM_MESH_NODAL_TETRA4;
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  // 13 -> tetra_ho
  // 14 -> pyramid_ho
  // 15 -> prism_ho
  // 16 -> hexa_ho

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nx,
             &ny,
             &nz,
             &center_x,
             &center_y,
             &center_z,
             &order,
             &t_elt,
             &radius);


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
   *  Create distributed spheres
   */
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_sphere_vol_gen_nodal(comm,
                           nx,
                           ny,
                           nz,
                           radius,
                           center_x,
                           center_y,
                           center_z,
                           t_elt,
                           order,
                           &dmn);

  PDM_dmesh_nodal_t *dmn2 = NULL;
  PDM_sphere_surf_icosphere_gen_nodal(comm,
                                      nx,
                                      center_x,
                                      center_y,
                                      center_z,
                                      radius,
                                      &dmn2);


  /* Reorder HO for visu */
  if (t_elt > PDM_MESH_NODAL_HEXA8) {
    if (order > 3) {
      int *ijk = NULL;

      for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_BARHO;
           type <= PDM_MESH_NODAL_HEXAHO;
           type++) {

        if (type == PDM_MESH_NODAL_PYRAMIDHO) continue;

        ijk = PDM_vtk_lagrange_to_ijk(type, order);
        PDM_ho_ordering_user_to_ijk_add ("PDM_HO_ORDERING_VTK",
                                         type,
                                         order,
                                         PDM_Mesh_nodal_n_vtx_elt_get(type, order),
                                         ijk);
        free (ijk);
      }
    }

    PDM_dmesh_nodal_reorder(dmn,
                            "PDM_HO_ORDERING_VTK");
  }

  if(0 == 1) {
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_SURFACIC, "sphere_vol_surface");
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_VOLUMIC,  "sphere_vol_volume");

    _dmesh_nodal_dump_vtk(dmn2, 1, PDM_GEOMETRY_KIND_SURFACIC, "sphere_surf_surface");
  }

  // Free memory
  PDM_DMesh_nodal_free(dmn);
  PDM_DMesh_nodal_free(dmn2);

  PDM_MPI_Finalize();

  return 0;
}
