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

#include "pdm_lagrange_to_bezier.h"
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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *nx,
           PDM_g_num_t   *ny,
           PDM_g_num_t   *nz,
           int           *order,
           int           *t_elt,
           double        *length,
           int           *n_part,
           int           *post,
           int           *method)
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
  } else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC && dmn->mesh_dimension == 3) {
    dmne = dmn->surfacic;
  } else if (geom_kind == PDM_GEOMETRY_KIND_CORNER) {
    dmne = dmn->corner;
    log_trace("Corners\n");
  }

  PDM_g_num_t *distrib_elt = NULL;
  if (dmne != NULL) {
    distrib_elt = PDM_compute_uniform_entity_distribution(dmn->comm,
                                                          dmne->n_g_elmts);
    PDM_log_trace_array_long(distrib_elt, n_rank+1, "distrib_elt : ");

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

    log_trace("section %d (%d) :\n", i_section, id_section);
    PDM_log_trace_array_long(delmt_ln_to_gn, n_elt, "  delmt_ln_to_gn : ");

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
      PDM_log_trace_connectivity_int(pelt_group_idx, pelt_group, n_elt, "pelt_group : ");
      free (tmp_elt_group_idx);
      free (tmp_elt_group);
      PDM_log_trace_array_long(delmt_ln_to_gn, n_elt, "  delmt_ln_to_gn (shifted) : ");

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



static void
_bezier_bounding_boxes
(
 PDM_dmesh_nodal_t   *dmn,
 int                  order,
 PDM_geometry_kind_t  geom_kind,
 const char          *filename_patter
 )
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, geom_kind);
  int  n_section   = PDM_DMesh_nodal_n_section_get  (dmn, geom_kind);

  PDM_g_num_t shift = 0;
  for (int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = sections_id[i_section];
    const PDM_g_num_t    *delmt_distribution = PDM_DMesh_nodal_distrib_section_get(dmn, geom_kind, id_section);
    int                   n_elt              = PDM_DMesh_nodal_section_n_elt_get  (dmn, geom_kind, id_section);
    PDM_g_num_t          *dconnec            = PDM_DMesh_nodal_section_std_get    (dmn, geom_kind, id_section);
    PDM_Mesh_nodal_elt_t  t_elt              = PDM_DMesh_nodal_section_type_get   (dmn, geom_kind, id_section);

    if (t_elt != PDM_MESH_NODAL_BAR2      &&
        t_elt != PDM_MESH_NODAL_TRIA3     &&
        t_elt != PDM_MESH_NODAL_QUAD4     &&
        t_elt != PDM_MESH_NODAL_BARHO     &&
        t_elt != PDM_MESH_NODAL_TRIAHO    &&
        t_elt != PDM_MESH_NODAL_QUADHO    &&
        t_elt != PDM_MESH_NODAL_TETRA4    &&
        t_elt != PDM_MESH_NODAL_TETRAHO   &&
        t_elt != PDM_MESH_NODAL_PYRAMID5  &&
        t_elt != PDM_MESH_NODAL_PYRAMIDHO &&
        t_elt != PDM_MESH_NODAL_PRISM6    &&
        t_elt != PDM_MESH_NODAL_PRISMHO   &&
        t_elt != PDM_MESH_NODAL_HEXA8     &&
        t_elt != PDM_MESH_NODAL_HEXAHO) continue;

    int         *dconnec_idx    = (int         * ) malloc( (n_elt+1) * sizeof(int        ));
    PDM_g_num_t *delmt_ln_to_gn = (PDM_g_num_t * ) malloc( (n_elt  ) * sizeof(PDM_g_num_t));

    int strid = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    dconnec_idx[0] = 0;
    for(int i = 0; i < n_elt; ++i) {
      dconnec_idx[i+1] = dconnec_idx[i] + strid;
      delmt_ln_to_gn[i] = delmt_distribution[i_rank] + i + 1;
    }

    int *ijk_to_vtk = PDM_ho_ordering_ijk_to_user_get("PDM_HO_ORDERING_VTK",
                                                      t_elt,
                                                      order);


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

    int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    double *lagrange_coord = malloc (sizeof(double) * n_nodes * 3);
    double *bezier_coord   = malloc (sizeof(double) * n_nodes * 3);
    double *elt_coord      = malloc (sizeof(double) * n_elt * n_nodes * 3);
    int    *elt_vtx        = malloc (sizeof(int)    * n_elt * n_nodes);

    double *matrix = NULL;
    if (order > 3) {
      int n_nodes_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);
      matrix = malloc(sizeof(double) * n_nodes_quad * n_nodes_quad);
    }

    double *extents = malloc (sizeof(double) * n_elt * 6);
    int idx2 = 0;
    for (int i = 0; i < n_elt; i++) {
      double *_min = extents + 6*i;
      double *_max = _min + 3;

      for (int j = 0; j < 3; j++) {
        _min[j] =  1e30;
        _max[j] = -1e30;
      }

      int idx = 0;
      for (int k = pcell_vtx_idx[i]; k < pcell_vtx_idx[i+1]; k++) {
        int ivtx = pcell_vtx[k] - 1;
        for (int j = 0; j < 3; j++) {
          lagrange_coord[idx++] = pvtx_coord_out[3*ivtx + j];
        }
      }

      if (t_elt == PDM_MESH_NODAL_BAR2 ||
          t_elt == PDM_MESH_NODAL_BARHO) {
        PDM_lagrange_to_bezier_bar(order, lagrange_coord, bezier_coord, matrix);
      }
      else if (t_elt == PDM_MESH_NODAL_TRIA3 ||
               t_elt == PDM_MESH_NODAL_TRIAHO) {
        PDM_lagrange_to_bezier_tria(order, lagrange_coord, bezier_coord, matrix);
      }
      else if (t_elt == PDM_MESH_NODAL_QUAD4 ||
               t_elt == PDM_MESH_NODAL_QUADHO) {
        PDM_lagrange_to_bezier_quad(order, lagrange_coord, bezier_coord, matrix);
      }
      else if (t_elt == PDM_MESH_NODAL_TETRA4 ||
               t_elt == PDM_MESH_NODAL_TETRAHO) {
        PDM_lagrange_to_bezier_tetra(order, lagrange_coord, bezier_coord, matrix);
      }
      else if (t_elt == PDM_MESH_NODAL_PYRAMID5 ||
               t_elt == PDM_MESH_NODAL_PYRAMIDHO) {
        PDM_lagrange_to_bezier_pyramid(order, lagrange_coord, bezier_coord, matrix);
      }
      else if (t_elt == PDM_MESH_NODAL_PRISM6 ||
               t_elt == PDM_MESH_NODAL_PRISMHO) {
        PDM_lagrange_to_bezier_prism(order, lagrange_coord, bezier_coord, matrix);
      }
      else if (t_elt == PDM_MESH_NODAL_HEXA8 ||
               t_elt == PDM_MESH_NODAL_HEXAHO) {
        PDM_lagrange_to_bezier_hexa(order, lagrange_coord, bezier_coord, matrix);
      }

      for (int k = 0; k < n_nodes; k++) {
        for (int j = 0; j < 3; j++) {
          elt_coord[3*idx2 + j] = bezier_coord[3*k + j];//lagrange_coord[3*k + j];//
          _min[j] = _MIN(_min[j], bezier_coord[3*k + j]);
          _max[j] = _MAX(_max[j], bezier_coord[3*k + j]);
        }
        elt_vtx[n_nodes*i + ijk_to_vtk[k]] = ++idx2;
      }
    }

    if (matrix != NULL) {
      free(matrix);
    }

    /*
     *  Dump
     */
    for(int i = 0; i < n_elt; ++i) {
      delmt_ln_to_gn[i] += shift;
    }

    char filename[999];
    sprintf(filename, "%s_bezier_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    PDM_vtk_write_std_elements_ho(filename,
                                  order,
                                  n_elt * n_nodes,
                                  elt_coord,
                                  NULL,
                                  t_elt,
                                  n_elt,
                                  elt_vtx,
                                  delmt_ln_to_gn,
                                  0,
                                  NULL,
                                  NULL);

    // sprintf(filename, "%s_bezier_coord_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    // PDM_vtk_write_point_cloud(filename,
    //                           n_elt * n_nodes,
    //                           elt_coord,
    //                           NULL,
    //                           NULL);
    free(elt_vtx);
    free(elt_coord);

    sprintf(filename, "%s_boxes_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    PDM_vtk_write_boxes(filename,
                        n_elt,
                        extents,
                        delmt_ln_to_gn);
    free (extents);
    free (bezier_coord);
    free (lagrange_coord);

    free(tmp_pvtx_coord);
    free(pvtx_ln_to_gn);
    free(pcell_vtx_idx);
    free(pcell_vtx);

    free(dconnec_idx);
    free(delmt_ln_to_gn);

    free(pvtx_coord_out);

    shift += delmt_distribution[n_rank];
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
  double amplitude = 0.1;//0.07;
  double frequency = 4.;

  for (int i = 0; i < n_vtx; i++) {
    double x = (vtx_coord[3*i    ] - 0.5) / length;
    double y = (vtx_coord[3*i + 1] - 0.5) / length;
    double z = (vtx_coord[3*i + 2] - 0.5) / length;

    vtx_coord[3*i    ] += amplitude*length*cos(frequency*y);
    vtx_coord[3*i + 1] += amplitude*length*cos(frequency*z);
    vtx_coord[3*i + 2] += amplitude*length*cos(frequency*x);
  }

  // double amplitude = 0.25 * length / (double) (n_vtx_seg - 1);
  // double frequency = PDM_PI * (n_vtx_seg - 1) / length;

  // for (int i = 0; i < n_vtx; i++) {
  //   double d[3];
  //   double dmax = -1;
  //   int jmax = -1;
  //   for (int j = 0; j < 3; j++) {
  //     d[j] = vtx_coord[3*i + j] - 0.5*length;
  //     if (PDM_ABS(d[j]) > dmax) {
  //       dmax = PDM_ABS(d[j]);
  //       jmax = j;
  //     }
  //   }

  //   double md = PDM_MODULE(d);
  //   if (md > 1e-15) {
  //     md = 1. / md;

  //     double x1 = vtx_coord[3*i + (jmax+1)%3];
  //     double x2 = vtx_coord[3*i + (jmax+2)%3];

  //     double s = amplitude * PDM_ABS( sin(frequency*x1) * sin(frequency*x2) );

  //     // double s = -amplitude * (sin(frequency*x1) + sin(frequency*x2));

  //     // for (int j = 0; j < 3; j++) {
  //     //   vtx_coord[3*i + j] += s * d[j] * md;
  //     // }
  //     vtx_coord[3*i + jmax] += s * PDM_SIGN(d[jmax]) * md;
  //   }
  // }
}


/**
 *
 * \brief  Main
 *
 */
// @@@param[n_proc] : 1,2,3,4
// @@@param[n] : 10,20,30,40
// @@@param[l] : 1.
// @@@args[part_kind] : -parmetis, -pt-scotch
int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */
  PDM_g_num_t          nx     = 10;
  PDM_g_num_t          ny     = 10;
  PDM_g_num_t          nz     = 10;
  int                  order  = 1;
  double               length = 1.;
  int                  n_part = 1;
  int                  post   = 0;
  PDM_Mesh_nodal_elt_t t_elt  = PDM_MESH_NODAL_TRIA3;
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

  PDM_part_split_t method  = PDM_PART_SPLIT_HILBERT;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nx,
             &ny,
             &nz,
             &order,
             (int *) &t_elt,
             &length,
             &n_part,
             &post,
             (int *) &method);

  if (t_elt == PDM_MESH_NODAL_TRIA3    ||
      t_elt == PDM_MESH_NODAL_QUAD4    ||
      t_elt == PDM_MESH_NODAL_TETRA4   ||
      t_elt == PDM_MESH_NODAL_PYRAMID5 ||
      t_elt == PDM_MESH_NODAL_PRISM6   ||
      t_elt == PDM_MESH_NODAL_HEXA8) {
    if (order != 1) {
      PDM_error(__FILE__, __LINE__, 0, "Invalid order %d for linear element type %d\n", order, (int) t_elt);
    }
  }

  /*
   *  Init
   */
  struct timeval t_elaps_debut;

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dim = 3;
  if (t_elt == PDM_MESH_NODAL_TRIA3  ||
      t_elt == PDM_MESH_NODAL_QUAD4  ||
      t_elt == PDM_MESH_NODAL_TRIAHO ||
      t_elt == PDM_MESH_NODAL_QUADHO) {
    dim = 2;
  }

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

  /*
   *  Create distributed cube
   */
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        nx,
                                                        ny,
                                                        nz,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        t_elt,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);

  /*PDM_dcube_nodal_gen_ordering_set (dcube,
    "PDM_HO_ORDERING_VTK");*/

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);


  /* Deform */
  // double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
  //                   {0.3129918,  0.9447025, -0.0978434},
  //                   {-0.1593451,  0.1537920,  0.9751703}};

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  // double amplitude = 0.1;//0.07;
  // double frequence = 4.;

  // if (1) {
  //   for (int i = 0; i < dn_vtx; i++) {
  //     double x = (dvtx_coord[3*i    ] - 0.5) / length;
  //     double y = (dvtx_coord[3*i + 1] - 0.5) / length;
  //     double z = (dvtx_coord[3*i + 2] - 0.5) / length;

  //     //double scale = length * pow(2, order-1);

  //     if (dim == 2) {
  //       //dvtx_coord[3*i + 2] = scale * (pow(x, order) + pow(y, order));
  //       dvtx_coord[3*i + 2] = length * (x*x + y*y);
  //     } else {
  //       dvtx_coord[3*i    ] += amplitude*length*cos(frequence*y);
  //       dvtx_coord[3*i + 1] += amplitude*length*cos(frequence*z);
  //       dvtx_coord[3*i + 2] += amplitude*length*cos(frequence*x);
  //     }
  //   }

  //   if (1) {
  //     for (int i = 0; i < dn_vtx; i++) {
  //       double x = dvtx_coord[3*i  ];
  //       double y = dvtx_coord[3*i+1];
  //       double z = dvtx_coord[3*i+2];

  //       for (int j = 0; j < 3; j++) {
  //         dvtx_coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
  //       }
  //     }
  //   }
  // }
  _deformation(length,
               nx,
               dn_vtx,
               dvtx_coord);

  if (post) {
    if (t_elt > PDM_MESH_NODAL_HEXA8) {
      /* Bounding boxes */
      if (dim == 3) {
        _bezier_bounding_boxes(dmn, order, PDM_GEOMETRY_KIND_VOLUMIC, "out_volumic");
      }
      _bezier_bounding_boxes(dmn, order, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
      _bezier_bounding_boxes(dmn, order, PDM_GEOMETRY_KIND_RIDGE,    "out_ridge");

      /* Reorder */
      PDM_dmesh_nodal_reorder (dmn,
                               "PDM_HO_ORDERING_VTK");
    }


    if (dim == 3) {
      _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_VOLUMIC, "out_volumic");
    }
    //PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    //PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_RIDGE,    "out_ridge");
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_CORNER,   "out_corner");
  }


  //PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  gettimeofday(&t_elaps_debut, NULL);
  PDM_dcube_nodal_gen_free(dcube);


  PDM_MPI_Finalize();

  return 0;
}
