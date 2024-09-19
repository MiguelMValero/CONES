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
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_polygon.h"
#include "pdm_array.h"
#include "pdm_triangulate.h"
#include "pdm_triangle.h"
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
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *n_vtx_seg,
           double                *length,
           int                   *n_part,
           int                   *post,
           int                   *part_method,
           PDM_Mesh_nodal_elt_t  *elt_type,
           int                   *triangulate)
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
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-tri") == 0) {
      *triangulate = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


// static void
// _compute_cell_vtx_hexa
// (
//  int  n_cell,
//  int *cell_face,
//  int *face_vtx_idx,
//  int *face_vtx,
//  int *cell_vtx
//  )
// {
//   int dbg = 0;


//   for (int icell = 0; icell < n_cell; icell++) {

//     int *cf = cell_face + 6*icell;
//     int *cv = cell_vtx  + 8*icell;
//     for (int i = 0; i < 8; i++) {
//       cv[i] = 0;
//     }

//     if (dbg) {
//       log_trace("\nCell %d\n", icell);
//       for (int i = 0; i < 6; i++) {
//         int iface = PDM_ABS(cf[i]) - 1;
//         int *fv = face_vtx + face_vtx_idx[iface];
//         log_trace(" Face %6d : %6d %6d %6d %6d\n", cf[i], fv[0], fv[1], fv[2], fv[3]);
//       }
//     }


//     // first face
//     int iface = cf[0];

//     if (iface < 0) {
//       iface = -iface - 1;
//       int *fv = face_vtx + face_vtx_idx[iface];
//       assert(face_vtx_idx[iface+1] - face_vtx_idx[iface] == 4);
//       for (int i = 0; i < 4; i++) {
//         cv[i] = fv[i];
//       }
//     }

//     else {
//       iface = iface - 1;
//       int *fv = face_vtx + face_vtx_idx[iface];
//       assert(face_vtx_idx[iface+1] - face_vtx_idx[iface] == 4);
//       for (int i = 0; i < 4; i++) {
//         cv[i] = fv[3-i];
//       }
//     }

//     if (dbg) {
//       log_trace("first 4 vtx : %d %d %d %d\n", cv[0], cv[1], cv[2], cv[3]);
//     }

//     // opposite face
//     int count = 0;

//     for (int idx_face = 1; idx_face < 6; idx_face++) {
//       int jface = PDM_ABS(cf[idx_face]) - 1;
//       int sgn   = PDM_SIGN(cf[idx_face]);

//       if (dbg) {
//         log_trace("  face %6d\n", cf[idx_face]);
//       }

//       int *fv = face_vtx + face_vtx_idx[jface];
//       assert(face_vtx_idx[jface+1] - face_vtx_idx[jface] == 4);

//       for (int i = 0; i < 4; i++) {
//         int ivtx1, ivtx2;
//         if (sgn > 0) {
//           ivtx1 = fv[i];
//           ivtx2 = fv[(i+1)%4];
//         } else {
//           ivtx2 = fv[i];
//           ivtx1 = fv[(i+1)%4];
//         }

//         if (dbg) {
//           log_trace("    edge %6d %6d\n", ivtx1, ivtx2);
//         }

//         if (ivtx1 == cv[0] && ivtx2 == cv[1]) {
//           if (sgn < 0) {
//             cv[4] = fv[(i+2)%4];
//             cv[5] = fv[(i+3)%4];
//           } else {
//             cv[4] = fv[(i+3)%4];
//             cv[5] = fv[(i+2)%4];
//           }
//           count++;
//         }

//         else if (ivtx1 == cv[2] && ivtx2 == cv[3]) {
//           if (sgn < 0) {
//             cv[6] = fv[(i+2)%4];
//             cv[7] = fv[(i+3)%4];
//           } else {
//             cv[6] = fv[(i+3)%4];
//             cv[7] = fv[(i+2)%4];
//           }
//           count++;
//         }
//       }

//       if (count == 2) {
//         break;
//       }

//     }

//     if (dbg) {
//       log_trace("count = %d\n", count);
//       log_trace("cv = %d %d %d %d %d %d %d %d\n",
//                 cv[0], cv[1], cv[2], cv[3],
//                 cv[4], cv[5], cv[6], cv[7]);
//     }
//     assert(count == 2);
//   }
// }



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

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TETRA4;
  int                  triangulate = 0;
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  //  9 -> poly3d

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method,
             &elt_type,
             &triangulate);
  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         0.,
                                                         0.,
                                                         0.,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  if (1) {
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

    double noise = 0.2 / (double) (n_vtx_seg - 1);
    for (int i = 0; i < 3*dn_vtx; i++) {
      dvtx_coord[i] += noise * (rand() / (double) RAND_MAX - 0.5);
    }
  }

  if(post) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  /*
   * Partitionnement
   */
  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
               &n_part_domains,
               PDM_FALSE,
               part_method,
               PDM_PART_SIZE_HOMOGENEOUS,
               NULL,
               comm,
               PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart, -1, "PDM_PART_RENUM_CELL_NONE",
                                                     NULL,
                                                     "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);

  for (int i_domain = 0; i_domain < n_domain; i_domain++){
    for (int i_part = 0; i_part < n_part; i_part++){

      int *cell_face     = NULL;
      int *cell_face_idx = NULL;
      int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                       i_domain,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       PDM_OWNERSHIP_KEEP);

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       i_domain,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                       &face_edge_idx,
                                                       &face_edge,
                                                       PDM_OWNERSHIP_KEEP);
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                       i_domain,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx_idx,
                                                       &edge_vtx,
                                                       PDM_OWNERSHIP_KEEP);
      assert(edge_vtx_idx == NULL);
      PDM_g_num_t* cell_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
      PDM_g_num_t* face_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* edge_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &edge_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* vtx_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);


      double *vtx = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx,
                                                   PDM_OWNERSHIP_KEEP);

      /*
       *  Compute additionnal connectivity
       */
      int *face_vtx_idx = NULL;
      int *face_vtx     = NULL;
      PDM_compute_face_vtx_from_face_and_edge(n_face,
                                              face_edge_idx,
                                              face_edge,
                                              edge_vtx,
                                              &face_vtx);
      face_vtx_idx = (int *) malloc(sizeof(int) * (n_face + 1));
      memcpy(face_vtx_idx, face_edge_idx, sizeof(int) * (n_face + 1));

      // PDM_log_trace_connectivity_int(face_vtx_idx, face_vtx, n_face, "face_vtx : ");


      int* cell_vtx_idx = NULL;
      int* cell_vtx     = NULL;
      PDM_combine_connectivity(n_cell, cell_face_idx, cell_face, face_vtx_idx, face_vtx, &cell_vtx_idx, &cell_vtx);

      int* vtx_cell_idx = NULL;
      int* vtx_cell     = NULL;
      PDM_connectivity_transpose(n_cell, n_vtx, cell_vtx_idx, cell_vtx, &vtx_cell_idx, &vtx_cell);



      /* Flag boundary edges (for error-checking) */
      int *_face_cell_idx = NULL;
      int *_face_cell     = NULL;
      PDM_connectivity_transpose(n_cell, n_face, cell_face_idx, cell_face, &_face_cell_idx, &_face_cell);

      int *is_bdr_face = PDM_array_zeros_int(n_face);
      int *is_bdr_vtx  = PDM_array_zeros_int(n_vtx);
      for (int iface = 0; iface < n_face; iface++) {
        if (_face_cell_idx[iface+1] == _face_cell_idx[iface]+1) {
          is_bdr_face[iface] = 1;
          for (int idx_vtx = face_vtx_idx[iface]; idx_vtx < face_vtx_idx[iface+1]; idx_vtx++) {
            int vtx_id = PDM_ABS(face_vtx[idx_vtx]) - 1;
            is_bdr_vtx[vtx_id] = 1;
          }
        }
      }
      free(_face_cell_idx);
      free(_face_cell);

      int *is_bdr_edge = PDM_array_zeros_int(n_edge);
      for (int iedge = 0; iedge < n_edge; iedge++) {
        for (int idx_vtx = 2*iedge; idx_vtx < 2*(iedge+1); idx_vtx++) {
          int vtx_id = edge_vtx[idx_vtx] - 1;
          if (is_bdr_vtx[vtx_id]) {
            is_bdr_edge[iedge] = 1;
            break;
          }
        }
      }

      // PDM_vtk_write_polydata("bdr_faces.vtk",
      //                        n_vtx,
      //                        vtx,
      //                        vtx_ln_to_gn,
      //                        n_face,
      //                        face_vtx_idx,
      //                        face_vtx,
      //                        face_ln_to_gn,
      //                        is_bdr_face);

      // const char *field_name[] = {"is_bdr_edge"};
      // PDM_vtk_write_std_elements("bdr_edges.vtk",
      //                            n_vtx,
      //                            vtx,
      //                            vtx_ln_to_gn,
      //                            PDM_MESH_NODAL_BAR2,
      //                            n_edge,
      //                            edge_vtx,
      //                            edge_ln_to_gn,
      //                            1,
      //                            field_name,
      //                            (const int **) &is_bdr_edge);

      int    *upwind_cell_out    = NULL;
      int    *downwind_cell_out  = NULL;
      int    *upwind_face_out    = NULL;
      int    *downwind_face_out  = NULL;
      double *upwind_point_out   = NULL;
      double *downwind_point_out = NULL;

      double *face_center = NULL;
      double *face_normal = NULL;
      if (!triangulate) {
        face_center = malloc(sizeof(double) * n_face * 3);
        face_normal = malloc(sizeof(double) * n_face * 3);

        PDM_geom_elem_polygon_properties(n_face,
                                         face_vtx_idx,
                                         face_vtx,
                                         vtx,
                                         face_normal,
                                         face_center,
                                         NULL,
                                         NULL);
      }

      PDM_geom_elem_edge_upwind_and_downwind(n_face,
                                             n_edge,
                                             cell_ln_to_gn,
                                             cell_face_idx,
                                             cell_face,
                                             face_vtx_idx,
                                             face_vtx,
                                             vtx_cell_idx,
                                             vtx_cell,
                                             edge_vtx,
                                             vtx,
                                             face_center,
                                             face_normal,
                                             &upwind_cell_out,
                                             &downwind_cell_out,
                                             &upwind_face_out,
                                             &downwind_face_out,
                                             &upwind_point_out,
                                             &downwind_point_out);

      if (!triangulate) {
        free(face_center);
        free(face_normal);
      }

      /* Check for errors */
      for (int i = 0; i < n_edge; i++) {
        if (is_bdr_edge[i] == 0) {
          if (upwind_face_out[i] < 0 || downwind_face_out[i] < 0) {
            printf("error edge %d ("PDM_FMT_G_NUM"): %d %d\n",
                   i, edge_ln_to_gn[i],
                   upwind_face_out[i], downwind_face_out[i]);
          }
        }
      }

      free(is_bdr_face);
      free(is_bdr_vtx);
      free(is_bdr_edge);


      /* Visualisation */
      if (post) {
        PDM_vtk_write_polydata("check_faces.vtk",
                               n_vtx,
                               vtx,
                               vtx_ln_to_gn,
                               n_face,
                               face_vtx_idx,
                               face_vtx,
                               face_ln_to_gn,
                               NULL);

        PDM_vtk_write_std_elements("check_edges.vtk",
                                   n_vtx,
                                   vtx,
                                   vtx_ln_to_gn,
                                   PDM_MESH_NODAL_BAR2,
                                   n_edge,
                                   edge_vtx,
                                   edge_ln_to_gn,
                                   0,
                                   NULL,
                                   NULL);

        double *updown = malloc(sizeof(double) * n_edge * 6);
        for (int i = 0; i < n_edge; i++) {
          int ivtx1 = edge_vtx[2*i  ]-1;
          int ivtx2 = edge_vtx[2*i+1]-1;

          if (upwind_face_out[i] < 0) {
            for (int j = 0; j < 3; j++) {
              updown[6*i+j] = 0.5*(vtx[3*ivtx1+j] + vtx[3*ivtx2+j]);
            }
          }
          else {
            memcpy(updown + 6*i, upwind_point_out + 3*i, sizeof(double)*3);
          }

          if (downwind_face_out[i] < 0) {
            for (int j = 0; j < 3; j++) {
              updown[6*i+3+j] = 0.5*(vtx[3*ivtx1+j] + vtx[3*ivtx2+j]);
            }
          }
          else {
            memcpy(updown + 6*i+3, downwind_point_out + 3*i, sizeof(double)*3);
          }
        }

        PDM_vtk_write_lines("updown.vtk",
                            n_edge,
                            updown,
                            edge_ln_to_gn,
                            NULL);




        int _n_face = 0;
        int    *_face_vtx_idx = malloc(sizeof(int) * (2*n_edge + 1));
        int    *_face_vtx     = malloc(sizeof(int) * 2*face_vtx_idx[n_face]);
        double *_face_edge    = malloc(sizeof(double) * 2*n_edge);
        double *_face_updown  = malloc(sizeof(double) * 2*n_edge);
        PDM_g_num_t *_face_ln_to_gn = malloc(sizeof(PDM_g_num_t) * 2*n_edge);
        _face_vtx_idx[0] = 0;



        int _n_cellface = 0;
        int _s_cellface = 6*2*n_edge;
        int _s_cellface_vtx = _s_cellface * 8;
        int *_cellface_vtx_idx = malloc(sizeof(int) * (_s_cellface + 1));
        int *_cellface_vtx     = malloc(sizeof(int) * _s_cellface_vtx);
        _cellface_vtx_idx[0] = 0;

        double      *_cellface_updown   = malloc(sizeof(double     ) * _s_cellface);
        PDM_g_num_t *_cellface_ln_to_gn = malloc(sizeof(PDM_g_num_t) * _s_cellface);

        double *_updown_pts = malloc(sizeof(double) * 2*n_edge*3);
        for (int iedge = 0; iedge < n_edge; iedge++) {

          if (upwind_face_out[iedge] >= 0) {
            _face_vtx_idx[_n_face+1] = _face_vtx_idx[_n_face];
            int face_id = upwind_face_out[iedge];
            for (int idx_vtx = face_vtx_idx[face_id]; idx_vtx < face_vtx_idx[face_id+1]; idx_vtx++) {
              _face_vtx[_face_vtx_idx[_n_face+1]++] = face_vtx[idx_vtx];
            }

            _face_edge[_n_face]   = iedge;
            _face_updown[_n_face] = 1;
            _face_ln_to_gn[_n_face] = edge_ln_to_gn[iedge];//face_ln_to_gn[face_id];
            memcpy(_updown_pts + 3*_n_face, upwind_point_out + 3*iedge, sizeof(double)*3);
            _n_face++;


            // realloc cellface?

            int cell_id = upwind_cell_out[iedge];
            for (int idx_face = cell_face_idx[cell_id]; idx_face < cell_face_idx[cell_id+1]; idx_face++) {
              int _face_id = PDM_ABS(cell_face[idx_face]) - 1;
              _cellface_vtx_idx[_n_cellface+1] = _cellface_vtx_idx[_n_cellface];
              for (int idx_vtx = face_vtx_idx[_face_id]; idx_vtx < face_vtx_idx[_face_id+1]; idx_vtx++) {
                _cellface_vtx[_cellface_vtx_idx[_n_cellface+1]++] = face_vtx[idx_vtx];
              }

              _cellface_updown  [_n_cellface] = 1;
              _cellface_ln_to_gn[_n_cellface] = edge_ln_to_gn[iedge];
              _n_cellface++;
            }
          }

          if (downwind_face_out[iedge] >= 0) {
            _face_vtx_idx[_n_face+1] = _face_vtx_idx[_n_face];
            int face_id = downwind_face_out[iedge];
            for (int idx_vtx = face_vtx_idx[face_id]; idx_vtx < face_vtx_idx[face_id+1]; idx_vtx++) {
              _face_vtx[_face_vtx_idx[_n_face+1]++] = face_vtx[idx_vtx];
            }

            _face_edge[_n_face]   = iedge;
            _face_updown[_n_face] = -1;
            _face_ln_to_gn[_n_face] = edge_ln_to_gn[iedge];//face_ln_to_gn[face_id];
            memcpy(_updown_pts + 3*_n_face, downwind_point_out + 3*iedge, sizeof(double)*3);
            _n_face++;


            // realloc cellface?

            int cell_id = downwind_cell_out[iedge];
            for (int idx_face = cell_face_idx[cell_id]; idx_face < cell_face_idx[cell_id+1]; idx_face++) {
              int _face_id = PDM_ABS(cell_face[idx_face]) - 1;
              _cellface_vtx_idx[_n_cellface+1] = _cellface_vtx_idx[_n_cellface];
              for (int idx_vtx = face_vtx_idx[_face_id]; idx_vtx < face_vtx_idx[_face_id+1]; idx_vtx++) {
                _cellface_vtx[_cellface_vtx_idx[_n_cellface+1]++] = face_vtx[idx_vtx];
              }

              _cellface_updown  [_n_cellface] = -1;
              _cellface_ln_to_gn[_n_cellface] = edge_ln_to_gn[iedge];
              _n_cellface++;
            }
          }

        }


        PDM_vtk_write_polydata_field("updown_faces.vtk",
                                     n_vtx,
                                     vtx,
                                     vtx_ln_to_gn,
                                     _n_face,
                                     _face_vtx_idx,
                                     _face_vtx,
                                     _face_ln_to_gn,
                                     "updown",
                                     _face_updown,
                                     NULL,
                                     NULL);


        PDM_vtk_write_polydata_field("updown_cells.vtk",
                                     n_vtx,
                                     vtx,
                                     vtx_ln_to_gn,
                                     _n_cellface,
                                     _cellface_vtx_idx,
                                     _cellface_vtx,
                                     _cellface_ln_to_gn,
                                     "updown",
                                     _cellface_updown,
                                     NULL,
                                     NULL);

        const char   *field_name[1]  = {"updown"};
        const double *field_value[1] = {_face_updown};
        PDM_vtk_write_std_elements_double("updown_pts.vtk",
                                          _n_face,
                                          _updown_pts,
                                          _face_ln_to_gn,
                                          PDM_MESH_NODAL_POINT,
                                          _n_face,
                                          NULL,
                                          _face_ln_to_gn,
                                          1,
                                          field_name,
                                          field_value);

        free(_face_vtx_idx);
        free(_face_vtx);
        free(_face_edge);
        free(_face_updown);
        free(_face_ln_to_gn);
        free(_cellface_vtx_idx);
        free(_cellface_vtx);
        free(_cellface_updown);
        free(_cellface_ln_to_gn);
      }




      free(upwind_cell_out);
      free(downwind_cell_out);
      free(upwind_face_out);
      free(downwind_face_out);
      free(upwind_point_out);
      free(downwind_point_out);


      free(vtx_cell_idx);
      free(vtx_cell);
      free(cell_vtx_idx);
      free(cell_vtx);
      free(face_vtx_idx);
      free(face_vtx);

    }
  }


  PDM_multipart_free(mpart);


  PDM_dcube_nodal_gen_free(dcube);


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;

}
