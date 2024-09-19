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
#include "pdm_dist_cloud_surf.h"
#include "pdm_gnum.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dbbtree.h"
#include "pdm_geom_elem.h"
#include "pdm_mesh_intersection.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

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
     "  -n      <level>  Number vtx in side of mesh A (default : 10).\n\n"
     "  -n_part <level>  Number of partition                        .\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
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
 int                   *n_part,
 int                   *post,
 PDM_Mesh_nodal_elt_t  *elt_type,
 double                *tolerance
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
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-tol") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
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


  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  double pi = 4 * atan(1.);
  for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {

    double x = dvtx_coord[3*i_vtx  ];
    double y = dvtx_coord[3*i_vtx+1];
    double z = dvtx_coord[3*i_vtx+2];

    dvtx_coord[3*i_vtx+1] = tanh( y * y * y);
    y = dvtx_coord[3*i_vtx+1];

    double angle = -pi/4;
    double Rz[3][3] = {{cos(angle), -sin(angle), 0},
                       {sin(angle),  cos(angle), 0},
                       {0         ,  0         , 1}};

    // dvtx_coord[3*i_vtx  ] = x +
    // dvtx_coord[3*i_vtx+1] = cos(y) - sin(x);

    if( x > 0.) {
      dvtx_coord[3*i_vtx  ] *= 4;
    }
    x = dvtx_coord[3*i_vtx];

    for (int j = 0; j < 1; j++) {
      dvtx_coord[3*i_vtx+j] = Rz[j][0]*x + Rz[j][1]*y + Rz[j][2]*z;
    }


  }

  PDM_UNUSED(rotate);
  // if(rotate) {
  //   // Do something
  //   double pi = 4 * atan(1.);
  //   double angle = pi/5.;
  //   PDM_g_num_t* distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  //   int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
  //   double* vtx_coord = PDM_DMesh_nodal_vtx_get(dmn);
  //   for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
  //     _rotate_coord(angle, &vtx_coord[3*i_vtx]);
  //   }
  // }

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
_cell_center_3d
(
  int      pn_cell,
  int     *pcell_face_idx,
  int     *pcell_face,
  int     *pface_edge_idx,
  int     *pface_edge,
  int     *pface_vtx_idx,
  int     *pface_vtx,
  int     *pedge_vtx,
  double  *pvtx_coord,
  double **cell_center
)
{
  int from_edge = 0;
  int from_face = 0;
  if(pface_edge     != NULL) {
    from_edge = 1;
  }
  if(pface_vtx     != NULL) {
    from_face = 1;
  }
  assert(pvtx_coord     != NULL);

  double* entity_center = malloc(3 * pn_cell * sizeof(double ));

  if(from_face == 1) {
    for(int i_cell = 0; i_cell < pn_cell; ++i_cell) {

      entity_center[3*i_cell  ] = 0.;
      entity_center[3*i_cell+1] = 0.;
      entity_center[3*i_cell+2] = 0.;

      double inv = 1./((double) pcell_face_idx[i_cell+1] - pcell_face_idx[i_cell]);

      for(int idx_face = pcell_face_idx[i_cell]; idx_face < pcell_face_idx[i_cell+1]; ++idx_face) {
        int i_face = PDM_ABS(pcell_face[idx_face])-1;

        double inv2 = 1./((double)  pface_vtx_idx[i_face+1] - pface_vtx_idx[i_face]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_vtx = pface_vtx_idx[i_face]; idx_vtx < pface_vtx_idx[i_face+1]; ++idx_vtx) {
          int i_vtx = pface_vtx[idx_vtx]-1;
          fcx += pvtx_coord[3*i_vtx  ];
          fcy += pvtx_coord[3*i_vtx+1];
          fcz += pvtx_coord[3*i_vtx+2];
        }
        fcx = fcx * inv2;
        fcy = fcy * inv2;
        fcz = fcz * inv2;

        entity_center[3*i_cell  ] += fcx;
        entity_center[3*i_cell+1] += fcy;
        entity_center[3*i_cell+2] += fcz;
      }

      entity_center[3*i_cell  ] = entity_center[3*i_cell  ] * inv;
      entity_center[3*i_cell+1] = entity_center[3*i_cell+1] * inv;
      entity_center[3*i_cell+2] = entity_center[3*i_cell+2] * inv;
    } /* End cell */
  } else if( from_edge == 1) {
    for(int i_cell = 0; i_cell < pn_cell; ++i_cell) {

      entity_center[3*i_cell  ] = 0.;
      entity_center[3*i_cell+1] = 0.;
      entity_center[3*i_cell+2] = 0.;

      double inv = 1./((double)  pcell_face_idx[i_cell+1] - pcell_face_idx[i_cell]);

      double fcx = 0;
      double fcy = 0;
      double fcz = 0;
      for(int idx_face = pcell_face_idx[i_cell]; idx_face < pcell_face_idx[i_cell+1]; ++idx_face) {
        int i_face = PDM_ABS(pcell_face[idx_face])-1;

        double inv2 = 1./((double)  pface_edge_idx[i_face+1] - pface_edge_idx[i_face]);

        for(int idx_edge = pface_edge_idx[i_face]; idx_edge < pface_edge_idx[i_face+1]; ++idx_edge) {
          int i_edge = PDM_ABS(pface_edge[idx_edge])-1;
          int i_vtx1 = pedge_vtx[2*i_edge  ] - 1;
          int i_vtx2 = pedge_vtx[2*i_edge+1] - 1;
          fcx += 0.5 * (pvtx_coord[3*i_vtx1  ] + pvtx_coord[3*i_vtx2  ]);
          fcy += 0.5 * (pvtx_coord[3*i_vtx1+1] + pvtx_coord[3*i_vtx2+1]);
          fcz += 0.5 * (pvtx_coord[3*i_vtx1+2] + pvtx_coord[3*i_vtx2+2]);
        }
        fcx = fcx * inv2;
        fcy = fcy * inv2;
        fcz = fcz * inv2;

        entity_center[3*i_cell  ] += fcx;
        entity_center[3*i_cell+1] += fcy;
        entity_center[3*i_cell+2] += fcz;
      }

      entity_center[3*i_cell  ] = entity_center[3*i_cell  ] * inv;
      entity_center[3*i_cell+1] = entity_center[3*i_cell+1] * inv;
      entity_center[3*i_cell+2] = entity_center[3*i_cell+2] * inv;
    } /* End cell */
  }

  *cell_center = entity_center;
}

static
void
_create_wall_surf
(
 const PDM_MPI_Comm           comm,
 const int                    n_part,
       PDM_multipart_t       *mpart,
       int                  **n_surf_vtx_out,
       int                  **n_surf_face_out,
       double              ***psurf_vtx_coord_out,
       int                 ***psurf_face_vtx_idx_out,
       int                 ***psurf_face_vtx_out,
       PDM_g_num_t         ***psurf_face_ln_to_gn_out,
       PDM_g_num_t         ***psurf_vtx_ln_to_gn_out,
       int                  **pn_cell_out,
       PDM_g_num_t         ***pcell_ln_to_gn_out,
       double              ***cell_center_out,
       int                   *pequi_surf_nface_out,
       int                   *pequi_surf_nvtx_out,
       int                  **pequi_surf_face_vtx_idx_out,
       int                  **pequi_surf_face_vtx_out,
       PDM_g_num_t          **pequi_surf_face_ln_to_gn_out,
       PDM_g_num_t          **pequi_surf_parent_face_ln_to_gn_out,
       PDM_g_num_t          **pequi_surf_vtx_ln_to_gn_out,
       double               **pequi_surf_vtx_coord_out
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  int          *n_surf_vtx           = malloc(n_part * sizeof(int          ));
  int          *n_surf_face          = malloc(n_part * sizeof(int          ));

  double      **cell_center          = malloc(n_part * sizeof(double      *));
  double      **psurf_vtx_coord      = malloc(n_part * sizeof(double      *));
  int         **psurf_face_vtx_idx   = malloc(n_part * sizeof(int         *));
  int         **psurf_face_vtx_n     = malloc(n_part * sizeof(int         *));
  int         **psurf_face_vtx       = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **psurf_face_ln_to_gn  = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **psurf_vtx_ln_to_gn   = malloc(n_part * sizeof(PDM_g_num_t *));

  PDM_g_num_t **psurf_face_vtx_g_num = malloc(n_part * sizeof(PDM_g_num_t *));

  int          *pn_cell        = malloc(n_part * sizeof(int          ));
  PDM_g_num_t **pcell_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));

  /* Compute gnum for vtx and faces */
  PDM_gen_gnum_t* gnum_face = PDM_gnum_create(3,
                                              n_part,
                                              PDM_FALSE,
                                              1.e-6,
                                              comm,
                                              PDM_OWNERSHIP_USER);

  PDM_gen_gnum_t* gnum_vtx = PDM_gnum_create(3,
                                             n_part,
                                             PDM_FALSE,
                                             1.e-6,
                                             comm,
                                             PDM_OWNERSHIP_USER);

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

    PDM_g_num_t *face_ln_to_gn = NULL;
    int  n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                  0,
                                                  i_part,
                                                  PDM_MESH_ENTITY_FACE,
                                                  &face_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *vtx_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    &vtx_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    &pcell_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);


    int *pcell_face     = NULL;
    int *pcell_face_idx = NULL;
    int n_cell =   PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &pcell_face_idx,
                                                       &pcell_face,
                                                       PDM_OWNERSHIP_KEEP);
    // PDM_log_trace_array_long(pcell_ln_to_gn[i_part], n_cell, "pcell_ln_to_gn : ");

    n_face = PDM_multipart_part_connectivity_get(mpart,
                                                 0,
                                                 i_part,
                                                 PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                 &face_edge_idx,
                                                 &face_edge,
                                                 PDM_OWNERSHIP_KEEP);
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx_idx,
                                        &edge_vtx,
                                        PDM_OWNERSHIP_KEEP);

    int  n_bound = 0;
    int* group_face_idx      = NULL;
    int* group_face          = NULL;
    PDM_g_num_t* face_group_ln_to_gn = NULL;

    PDM_multipart_group_get(mpart,
                            0,
                            i_part,
                            PDM_MESH_ENTITY_FACE,
                            &n_bound,
                            &group_face_idx,
                            &group_face,
                            &face_group_ln_to_gn,
                            PDM_OWNERSHIP_KEEP);

    int *face_vtx = NULL;
    int *face_vtx_idx = face_edge_idx;
    PDM_compute_face_vtx_from_face_and_edge(n_face, face_edge_idx, face_edge, edge_vtx, &face_vtx);

    pn_cell[i_part] = n_cell;
    _cell_center_3d(n_cell,
                    pcell_face_idx,
                    pcell_face,
                    NULL,
                    NULL,
                    face_vtx_idx,
                    face_vtx,
                    NULL,
                    vtx_coord,
                    &cell_center[i_part]);
    /*
     * Nez de la plaque en x = 0
     */
    int i_group = 4;
    n_surf_face[i_part] = 0;
    int n_surf_face_vtx = 0;
    int n_face_in_group = group_face_idx[i_group+1] - group_face_idx[i_group];
    int* face_bnd = malloc(n_face_in_group * sizeof(int));
    for(int idx_face = group_face_idx[i_group]; idx_face < group_face_idx[i_group+1]; ++idx_face) {
      int i_face = group_face[idx_face]-1;

      int is_in_plate = 1;
      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = face_vtx[idx_vtx]-1;
        if(vtx_coord[3*i_vtx] < 0.) {
          is_in_plate = 0;
        }
      }

      if(is_in_plate == 0) {
        continue;
      }

      face_bnd[n_surf_face[i_part]++] = i_face;

      n_surf_face_vtx += face_vtx_idx[i_face+1] - face_vtx_idx[i_face];

    }


    psurf_vtx_coord     [i_part] = malloc(3 * n_surf_face_vtx          * sizeof(double     ));
    psurf_face_vtx_idx  [i_part] = malloc(   ( n_surf_face[i_part] +1) * sizeof(int        ));
    psurf_face_vtx_n    [i_part] = malloc(   ( n_surf_face[i_part]   ) * sizeof(int        ));
    psurf_face_vtx      [i_part] = malloc(    n_surf_face_vtx          * sizeof(int        ));
    psurf_face_ln_to_gn [i_part] = malloc(    n_surf_face[i_part]      * sizeof(PDM_g_num_t));
    psurf_vtx_ln_to_gn  [i_part] = malloc(    n_surf_face_vtx          * sizeof(PDM_g_num_t));
    psurf_face_vtx_g_num[i_part] = malloc(    n_surf_face_vtx          * sizeof(PDM_g_num_t));

    double      *_psurf_vtx_coord      = psurf_vtx_coord     [i_part];
    int         *_psurf_face_vtx_idx   = psurf_face_vtx_idx  [i_part];
    int         *_psurf_face_vtx_n     = psurf_face_vtx_n    [i_part];
    int         *_psurf_face_vtx       = psurf_face_vtx      [i_part];
    PDM_g_num_t *_psurf_face_ln_to_gn  = psurf_face_ln_to_gn [i_part];
    PDM_g_num_t *_psurf_vtx_ln_to_gn   = psurf_vtx_ln_to_gn  [i_part];
    PDM_g_num_t *_psurf_face_vtx_g_num = psurf_face_vtx_g_num[i_part];

    int *vtx_flags = malloc(n_vtx * sizeof(int));
    for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
      vtx_flags[i_vtx] = -100;
    }
    n_surf_vtx[i_part] = 0;

    _psurf_face_vtx_idx[0] = 0;
    for(int idx_face = 0; idx_face < n_surf_face[i_part]; ++idx_face) {
      int i_face = face_bnd[idx_face];

      _psurf_face_vtx_idx[idx_face+1] = _psurf_face_vtx_idx[idx_face];
      _psurf_face_ln_to_gn[idx_face] = face_ln_to_gn[i_face];

      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = face_vtx[idx_vtx]-1;

        if(vtx_flags[i_vtx] == -100) {
          vtx_flags[i_vtx] = n_surf_vtx[i_part];
          _psurf_vtx_coord[3*n_surf_vtx[i_part]  ] = vtx_coord[3*i_vtx  ];
          _psurf_vtx_coord[3*n_surf_vtx[i_part]+1] = vtx_coord[3*i_vtx+1];
          _psurf_vtx_coord[3*n_surf_vtx[i_part]+2] = vtx_coord[3*i_vtx+2];

          _psurf_vtx_ln_to_gn[n_surf_vtx[i_part]] = vtx_ln_to_gn[i_vtx];

          n_surf_vtx[i_part]++;
        }

        _psurf_face_vtx      [_psurf_face_vtx_idx[idx_face+1]] = vtx_flags[i_vtx]+1;
        _psurf_face_vtx_g_num[_psurf_face_vtx_idx[idx_face+1]] = vtx_ln_to_gn[i_vtx];
        _psurf_face_vtx_idx[idx_face+1]++;
      }

      _psurf_face_vtx_n[idx_face] = face_vtx_idx[i_face+1] - face_vtx_idx[i_face];

    }

    // printf("n_surf_face = %i \n", n_surf_face[i_part]);
    // printf("n_surf_vtx  = %i \n", n_surf_vtx [i_part]);

    psurf_vtx_coord    [i_part] = realloc(psurf_vtx_coord    [i_part], 3 * n_surf_vtx[i_part] * sizeof(double     ));
    psurf_vtx_ln_to_gn [i_part] = realloc(psurf_vtx_ln_to_gn [i_part],     n_surf_vtx[i_part] * sizeof(PDM_g_num_t));

    PDM_gnum_set_from_parents(gnum_face, i_part, n_surf_face[i_part], _psurf_face_ln_to_gn);
    PDM_gnum_set_from_parents(gnum_vtx , i_part, n_surf_vtx [i_part], psurf_vtx_ln_to_gn [i_part] );

    free(vtx_flags);
    free(face_bnd);
    free(face_vtx);

  }

  PDM_gnum_compute(gnum_face);
  PDM_gnum_compute(gnum_vtx );



  /* Vtk en légende */
  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      char filename[999];
      sprintf(filename, "face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                             n_surf_vtx [i_part],
                             psurf_vtx_coord    [i_part],
                             psurf_vtx_ln_to_gn [i_part],
                             n_surf_face        [i_part],
                             psurf_face_vtx_idx [i_part],
                             psurf_face_vtx     [i_part],
                             psurf_face_ln_to_gn[i_part],
                             NULL);
    }

  }

  /*
   * Equilibrate wall
   */
  PDM_part_to_block_t *ptb_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                           psurf_vtx_ln_to_gn,
                                                           NULL,
                                                           n_surf_vtx,
                                                           n_part,
                                                           comm);

  double *dvtx_coords = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         3,
                         NULL,
              (void **)  psurf_vtx_coord,
                         NULL,
              (void **)  &dvtx_coords);

  int dn_vtx = PDM_part_to_block_n_elt_block_get(ptb_vtx);
  PDM_g_num_t *blk_vtx_gnum = PDM_part_to_block_block_gnum_get(ptb_vtx);
  // PDM_log_trace_array_long(blk_vtx_gnum, dn_vtx, "blk_vtx_gnum :");

  /* A laisser après le ptb vu qu'on ecrase psurf_vtx_ln_to_gn */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(psurf_face_ln_to_gn[i_part]);
    free(psurf_vtx_ln_to_gn [i_part]);
    psurf_face_ln_to_gn[i_part] = PDM_gnum_get(gnum_face, 0);
    psurf_vtx_ln_to_gn [i_part] = PDM_gnum_get(gnum_vtx , 0);
  }

  PDM_gnum_free(gnum_face);
  PDM_gnum_free(gnum_vtx);

  PDM_part_to_block_t *ptb_face = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                           psurf_face_ln_to_gn,
                                                           NULL,
                                                           n_surf_face,
                                                           n_part,
                                                           comm);

  PDM_g_num_t *distrib_face  = PDM_part_to_block_distrib_index_get(ptb_face);
  PDM_g_num_t *blk_face_gnum = PDM_part_to_block_block_gnum_get   (ptb_face);
  int dn_face = PDM_part_to_block_n_elt_block_get(ptb_face);

  int         *dface_vtx_n = NULL;
  PDM_g_num_t *dface_vtx   = NULL;
  PDM_part_to_block_exch(ptb_face,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                         psurf_face_vtx_n,
          (void **)      psurf_face_vtx_g_num,
                         &dface_vtx_n,
          (void **)      &dface_vtx);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(psurf_face_vtx_g_num[i_part]);
    free(psurf_face_vtx_n    [i_part]);
  }
  free(psurf_face_vtx_g_num);
  free(psurf_face_vtx_n);

  PDM_g_num_t *distrib_face_child = PDM_compute_entity_distribution(comm, dn_face);

  PDM_g_num_t* dface_ln_to_gn        = malloc(dn_face * sizeof(PDM_g_num_t));
  PDM_g_num_t* pface_parent_ln_to_gn = malloc(dn_face * sizeof(PDM_g_num_t));
  for(int i_face = 0; i_face < dn_face; ++i_face) {
    dface_ln_to_gn       [i_face] = distrib_face_child[i_rank] + i_face + 1;
    pface_parent_ln_to_gn[i_face] = blk_face_gnum[i_face];
  }

  int* dface_vtx_idx = PDM_array_new_idx_from_sizes_int(dface_vtx_n, dn_face);

  int          *tmp_pequi_surf_nvtx         = NULL;
  PDM_g_num_t **tmp_pequi_surf_vtx_ln_to_gn = NULL;
  int         **tmp_pequi_surf_face_vtx_idx = NULL;
  int         **tmp_pequi_surf_face_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               distrib_face,
                                               dface_vtx_idx,
                                               dface_vtx,
                                               1,
                                               &dn_face,
              (const PDM_g_num_t **)           &dface_ln_to_gn,
                                               &tmp_pequi_surf_nvtx,
                                               &tmp_pequi_surf_vtx_ln_to_gn,
                                               &tmp_pequi_surf_face_vtx_idx,
                                               &tmp_pequi_surf_face_vtx);

  int          pequi_surf_nvtx         = tmp_pequi_surf_nvtx        [0];
  PDM_g_num_t *pequi_surf_vtx_ln_to_gn = tmp_pequi_surf_vtx_ln_to_gn[0];
  int         *pequi_surf_face_vtx_idx = tmp_pequi_surf_face_vtx_idx[0];
  int         *pequi_surf_face_vtx     = tmp_pequi_surf_face_vtx    [0];
  free(tmp_pequi_surf_nvtx        );
  free(tmp_pequi_surf_vtx_ln_to_gn);
  free(tmp_pequi_surf_face_vtx_idx);
  free(tmp_pequi_surf_face_vtx    );
  free(dface_vtx_idx    );

  PDM_block_to_part_t* btp_vtx = PDM_block_to_part_create_from_sparse_block(blk_vtx_gnum,
                                                                            dn_vtx,
                                                     (const PDM_g_num_t **) &pequi_surf_vtx_ln_to_gn,
                                                                            &pequi_surf_nvtx,
                                                                            1,
                                                                            comm);

  int stride_one = 1;
  double **tmp_pvtx_coord = NULL;
  PDM_block_to_part_exch(btp_vtx,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         dvtx_coords,
                         NULL,
            (void ***)   &tmp_pvtx_coord);
  double* pequi_surf_vtx_coord = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);
  free(dvtx_coords);

  PDM_block_to_part_free(btp_vtx);

  PDM_part_to_block_free(ptb_face);
  PDM_part_to_block_free(ptb_vtx);
  free(distrib_face_child);

  free(dface_vtx_n);
  free(dface_vtx);


  *pn_cell_out        = pn_cell;
  *pcell_ln_to_gn_out = pcell_ln_to_gn;
  *cell_center_out    = cell_center;
  *n_surf_vtx_out     = n_surf_vtx;
  *n_surf_face_out    = n_surf_face;

  *psurf_vtx_coord_out      = psurf_vtx_coord;
  *psurf_face_vtx_idx_out   = psurf_face_vtx_idx;
  *psurf_face_vtx_out       = psurf_face_vtx;
  *psurf_face_ln_to_gn_out  = psurf_face_ln_to_gn;
  *psurf_vtx_ln_to_gn_out   = psurf_vtx_ln_to_gn;

  /* Vtk en légende */
  if(0 == 1) {
    char filename[999];
    sprintf(filename, "equi_face_vtx_coord_%3.3d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           pequi_surf_nvtx,
                           pequi_surf_vtx_coord,
                           pequi_surf_vtx_ln_to_gn,
                           dn_face,
                           pequi_surf_face_vtx_idx,
                           pequi_surf_face_vtx,
                           dface_ln_to_gn,
                           NULL);
  }


  /* Equi */
  *pequi_surf_nface_out                = dn_face;
  *pequi_surf_nvtx_out                 = pequi_surf_nvtx;
  *pequi_surf_face_vtx_idx_out         = pequi_surf_face_vtx_idx;
  *pequi_surf_face_vtx_out             = pequi_surf_face_vtx;
  *pequi_surf_face_ln_to_gn_out        = dface_ln_to_gn;
  *pequi_surf_parent_face_ln_to_gn_out = pface_parent_ln_to_gn;
  *pequi_surf_vtx_ln_to_gn_out         = pequi_surf_vtx_ln_to_gn;
  *pequi_surf_vtx_coord_out            = pequi_surf_vtx_coord;


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

    int *face_edge_idx;
    int *face_edge;
    int *edge_vtx_idx;
    int *edge_vtx;

    double *vtx_coord;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                 0,
                                                 i_part,
                                                 &vtx_coord,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *cell_ln_to_gn;
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
    PDM_multipart_part_connectivity_get(mpart,
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
_create_wall_ray
(
 const PDM_MPI_Comm     comm,
 const int              n_part,
       int            *n_surf_vtx,
       int            *n_surf_face,
       double        **psurf_vtx_coord,
       int           **psurf_face_vtx_idx,
       int           **psurf_face_vtx,
       PDM_g_num_t   **psurf_face_ln_to_gn,
       PDM_g_num_t   **psurf_vtx_ln_to_gn,
       int            *pn_ray_out,
       PDM_g_num_t   **pvtx_ln_to_gn_out,
       PDM_g_num_t   **pray_ln_to_gn_out,
       int           **pray_vtx_out,
       double        **pvtx_coord_out,
       double        **psurf_face_normal_out,
       double        **psurf_face_center_out
)
{
  PDM_UNUSED(n_surf_vtx);
  PDM_UNUSED(psurf_vtx_ln_to_gn);
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  int pn_ray = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_ray += n_surf_face[i_part];
  }

  PDM_g_num_t *pray_ln_to_gn = malloc(    pn_ray * sizeof(PDM_g_num_t));
  PDM_g_num_t *pvtx_ln_to_gn = malloc(2 * pn_ray * sizeof(PDM_g_num_t));
  double      *pray_coord    = malloc(6 * pn_ray * sizeof(double     ));
  double      *pface_normal  = malloc(3 * pn_ray * sizeof(double     ));
  double      *pface_center  = malloc(3 * pn_ray * sizeof(double     ));

  PDM_g_num_t *distrib_vtx = PDM_compute_entity_distribution (comm, 2 * pn_ray);

  int *pray_vtx = malloc(2 * pn_ray * sizeof(int));

  pn_ray = 0;
  int pn_vtx = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {

    double      *face_normal    = &pface_normal[pn_ray];
    double      *face_center    = &pface_center[pn_ray];
    PDM_geom_elem_polygon_properties(n_surf_face[i_part],
                                     psurf_face_vtx_idx[i_part],
                                     psurf_face_vtx    [i_part],
                                     psurf_vtx_coord   [i_part],
                                     face_normal,
                                     face_center,
                                     NULL,
                                     NULL);

    double dmax = 0.2;
    for(int i_face = 0; i_face < n_surf_face[i_part]; ++i_face) {

      pray_coord[6*pn_ray  ] = face_center[3*i_face  ];
      pray_coord[6*pn_ray+1] = face_center[3*i_face+1];
      pray_coord[6*pn_ray+2] = face_center[3*i_face+2];

      double nx = face_normal[3*i_face  ];
      double ny = face_normal[3*i_face+1];
      double nz = face_normal[3*i_face+2];

      double sn = sqrt(nx * nx + ny * ny + nz * nz);
      nx = -nx / sn;
      ny = -ny / sn;
      nz = -nz / sn;

      double dnx = nx * dmax;
      double dny = ny * dmax;
      double dnz = nz * dmax;

      double xb = pray_coord[6*pn_ray  ] + dnx;
      double yb = pray_coord[6*pn_ray+1] + dny;
      double zb = pray_coord[6*pn_ray+2] + dnz;

      pray_coord[6*pn_ray+3] = xb;
      pray_coord[6*pn_ray+4] = yb;
      pray_coord[6*pn_ray+5] = zb;

      pray_vtx[2*pn_ray  ] = pn_vtx+1;
      pray_vtx[2*pn_ray+1] = pn_vtx+2;

      pray_ln_to_gn[pn_ray++] = psurf_face_ln_to_gn[i_part][i_face];

      pvtx_ln_to_gn[pn_vtx  ] = distrib_vtx[i_rank] + pn_vtx + 1;
      pvtx_ln_to_gn[pn_vtx+1] = distrib_vtx[i_rank] + pn_vtx + 2;
      pn_vtx += 2;

    }
  }

  free(distrib_vtx);

  if(0 == 1) {
    char filename[999];
    sprintf(filename, "ray_%i.vtk", i_rank);
    PDM_vtk_write_lines(filename,
                        pn_ray,
                        pray_coord,
                        pray_ln_to_gn,
                        NULL);
  }


  *pn_ray_out            = pn_ray;
  *pvtx_ln_to_gn_out     = pvtx_ln_to_gn;
  *pray_ln_to_gn_out     = pray_ln_to_gn;
  *pvtx_coord_out        = pray_coord;
  *pray_vtx_out          = pray_vtx;
  *psurf_face_normal_out = pface_normal;
  *psurf_face_center_out = pface_center;

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
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_HEXA8;

  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;

  int n_part = 1;
  int post   = 0;
  double tolerance = 1e-6;

  _read_args(argc,
             argv,
             &n_vtx_a,
             &n_part,
             &post,
             &elt_type,
             &tolerance);

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
                         -0.2,
                         0.,
                         0.,
                         lenght_a,
                         part_method,
                         n_part,
                         &dmn_vol_a,
                         &mpart_vol_a);

  if(post) {
    PDM_dmesh_nodal_dump_vtk(dmn_vol_a,
                             PDM_GEOMETRY_KIND_VOLUMIC,
                             "dmn_vol_a_");
    PDM_dmesh_nodal_dump_vtk(dmn_vol_a,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "dmn_surf_a_");

  }

  /*
   * Extract boundary wall
   */
  int          *pn_cell              = NULL;
  PDM_g_num_t **pcell_ln_to_gn       = NULL;
  int          *psurf_vtx            = NULL;
  int          *psurf_face           = NULL;
  double      **psurf_vtx_coord      = NULL;
  int         **psurf_face_vtx_idx   = NULL;
  int         **psurf_face_vtx       = NULL;
  PDM_g_num_t **psurf_face_ln_to_gn  = NULL;
  PDM_g_num_t **psurf_vtx_ln_to_gn   = NULL;
  double      **cell_center          = NULL;

  int          pequi_surf_nface                = 0;
  int          pequi_surf_nvtx                 = 0;
  int         *pequi_surf_face_vtx_idx         = NULL;
  int         *pequi_surf_face_vtx             = NULL;
  PDM_g_num_t *pequi_surf_face_ln_to_gn        = NULL;
  PDM_g_num_t *pequi_surf_parent_face_ln_to_gn = NULL;
  PDM_g_num_t *pequi_surf_vtx_ln_to_gn         = NULL;
  double      *pequi_surf_vtx_coord            = NULL;
  _create_wall_surf(comm,
                    n_part,
                    mpart_vol_a,
                    &psurf_vtx,
                    &psurf_face,
                    &psurf_vtx_coord,
                    &psurf_face_vtx_idx,
                    &psurf_face_vtx,
                    &psurf_face_ln_to_gn,
                    &psurf_vtx_ln_to_gn,
                    &pn_cell,
                    &pcell_ln_to_gn,
                    &cell_center,
                    &pequi_surf_nface,
                    &pequi_surf_nvtx,
                    &pequi_surf_face_vtx_idx,
                    &pequi_surf_face_vtx,
                    &pequi_surf_face_ln_to_gn,
                    &pequi_surf_parent_face_ln_to_gn,
                    &pequi_surf_vtx_ln_to_gn,
                    &pequi_surf_vtx_coord);

  /* Wall distance */
  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t* dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            n_point_cloud,
                                                            comm,
                                                            PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (dist,
                                                 n_part);

  PDM_dist_cloud_surf_n_part_cloud_set (dist, 0, n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_dist_cloud_surf_surf_mesh_part_set (dist,
                                            i_part,
                                            psurf_face         [i_part],
                                            psurf_face_vtx_idx [i_part],
                                            psurf_face_vtx     [i_part],
                                            psurf_face_ln_to_gn[i_part],
                                            psurf_vtx          [i_part],
                                            psurf_vtx_coord    [i_part],
                                            psurf_vtx_ln_to_gn [i_part]);

    PDM_dist_cloud_surf_cloud_set (dist,
                                   0,
                                   i_part,
                                   pn_cell       [i_part],
                                   cell_center   [i_part],
                                   pcell_ln_to_gn[i_part]);
  }

  PDM_dist_cloud_surf_compute(dist);

  PDM_g_num_t **closest_elt_gnum = malloc(n_part * sizeof(PDM_g_num_t *));
  for (int i_part = 0; i_part < n_part; i_part++) {

    double      *distance;
    double      *projected;

    PDM_dist_cloud_surf_get (dist,
                             0,
                             i_part,
                             &distance,
                             &projected,
                             &closest_elt_gnum[i_part]);
    if(post) {
      char filename[999];
      sprintf(filename, "distance_%3.3d_%3.3d.vtk", i_part, i_rank);

      const char   *vector_field_name[1] = {"distance"};
      const double *vector_field     [1] = {distance};
      PDM_vtk_write_point_cloud_with_field(filename,
                                           pn_cell       [i_part],
                                           cell_center   [i_part],
                                           pcell_ln_to_gn[i_part],
                                           NULL,
                                           1,
                                           vector_field_name,
                                           vector_field,
                                           0,
                                           NULL,
                                           NULL,
                                           0,
                                           NULL,
                                           NULL );
    }
  }

  /* Create field of speed */
  double **velocity = malloc(n_part * sizeof(double));
  for (int i_part = 0; i_part < n_part; i_part++) {

    double      *distance;
    double      *projected;

    PDM_dist_cloud_surf_get (dist,
                             0,
                             i_part,
                             &distance,
                             &projected,
                             &closest_elt_gnum[i_part]);

    velocity[i_part] = malloc(pn_cell[i_part] * sizeof(double));

    for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {
      if(cell_center   [i_part][3*i_cell+1] < 0.1 * cell_center   [i_part][3*i_cell]) {
        velocity[i_part][i_cell] = tanh(6*distance[i_cell]);
      } else {
        velocity[i_part][i_cell] = 1.;
      }
    }


    if(post) {
      char filename[999];
      sprintf(filename, "velocity_%3.3d_%3.3d.vtk", i_part, i_rank);

      const char   *vector_field_name[1] = {"velocity"};
      const double *vector_field     [1] = {velocity[i_part]};
      PDM_vtk_write_point_cloud_with_field(filename,
                                           pn_cell       [i_part],
                                           cell_center   [i_part],
                                           pcell_ln_to_gn[i_part],
                                           NULL,
                                           1,
                                           vector_field_name,
                                           vector_field,
                                           0,
                                           NULL,
                                           NULL,
                                           0,
                                           NULL,
                                           NULL );
    }
  }

  int          n_lines           = 0;
  double      *ray_coord         = NULL;
  int         *pray_vtx          = NULL;
  PDM_g_num_t *pvtx_ln_to_gn     = NULL;
  PDM_g_num_t *pray_ln_to_gn     = NULL;
  double      *psurf_face_normal = NULL;
  double      *psurf_face_center = NULL;
  // _create_wall_ray(comm,
  //                  n_part,
  //                  psurf_vtx,
  //                  psurf_face,
  //                  psurf_vtx_coord,
  //                  psurf_face_vtx_idx,
  //                  psurf_face_vtx,
  //                  psurf_face_ln_to_gn,
  //                  psurf_vtx_ln_to_gn,
  //                  &n_lines,
  //                  &pvtx_ln_to_gn,
  //                  &pray_ln_to_gn,
  //                  &pray_vtx,
  //                  &ray_coord,
  //                  &psurf_face_normal,
  //                  &psurf_face_center);

  _create_wall_ray(comm,
                   1,
                   &pequi_surf_nvtx,
                   &pequi_surf_nface,
                   &pequi_surf_vtx_coord,
                   &pequi_surf_face_vtx_idx,
                   &pequi_surf_face_vtx,
                   &pequi_surf_face_ln_to_gn,
                   &pequi_surf_vtx_ln_to_gn,
                   &n_lines,
                   &pvtx_ln_to_gn,
                   &pray_ln_to_gn,
                   &pray_vtx,
                   &ray_coord,
                   &psurf_face_normal,
                   &psurf_face_center);

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
  _set_mesh(mi, 0, mpart_vol_a, n_part);

  /*
   * Set line mesh
   */
  PDM_mesh_intersection_n_part_set(mi, 1, 1);

  PDM_mesh_intersection_part_set(mi,
                                 1, // i_mesh
                                 0,
                                 0,
                                 0,
                                 n_lines,
                                 2 * n_lines,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 pray_vtx,
                                 NULL, // face_vtx_idx,
                                 NULL, // face_vtx,
                                 NULL,
                                 NULL,
                                 pray_ln_to_gn,
                                 pvtx_ln_to_gn,
                                 ray_coord);

  PDM_mesh_intersection_compute(mi);


  /*
   * Use of meshintersection for bl post-treatment
   */
  PDM_part_to_part_t* ptp = NULL;
  PDM_mesh_intersection_part_to_part_get(mi,
                                         &ptp,
                                         PDM_OWNERSHIP_KEEP);

  int n_part_cell = 0;
  int n_part_line = 0;
  PDM_part_to_part_n_part_get(ptp, &n_part_cell, &n_part_line);

  int  *n_ref_b = NULL;
  int **ref_b   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp,
                                 &n_ref_b,
                                 &ref_b);

  for(int i_part = 0; i_part < n_part_line; ++i_part) {
    printf("n_ref_b[%i] = %i \n", i_part, n_ref_b[i_part]);
  }
  assert(n_part_line == 1);


  /*
   * Pour Guillaume/Lucas :
   *     - On doit garder le lien avec la num absolu des faces
   *     - part_to_block avec les numéros de faces de wall
   *     - Donc il faut un block_to_part_from_sparse_block
   *     - A partir de ce block -> on refait un block de normal
   *     - A l'issu des l'intersection de maillage :
   *          - part1 = La maillage user classique
   *          - part2 = Le "singleton" paroi
   *     - le equi surface suit des nums asbolue croissant : on peut le prendre comme un bloc
   */
  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(pequi_surf_parent_face_ln_to_gn,
                                                                        pequi_surf_nface,
                                              (const PDM_g_num_t **)    closest_elt_gnum,
                                                                        pn_cell,
                                                                        n_part,
                                                                        comm);


  /*
   * Envoi pour chaque ligne le centre cellule
   */

  double **pline_to_cell_center = NULL;
  int request = -1;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         3 * sizeof(double),
                         NULL,
        (const void  **) cell_center,
                         NULL,
        (      void ***) &pline_to_cell_center,
                         &request);
  PDM_part_to_part_iexch_wait(ptp, request);

  double **pline_to_cell_velocity = NULL;
  request = -1;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) velocity,
                         NULL,
        (      void ***) &pline_to_cell_velocity,
                         &request);
  PDM_part_to_part_iexch_wait(ptp, request);



  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp,
                                       &gnum1_come_from_idx,
                                       &gnum1_come_from);
  int         *_gnum1_come_from_idx    = gnum1_come_from_idx   [0];
  PDM_g_num_t *_gnum1_come_from        = gnum1_come_from       [0];
  double      *_pline_to_cell_center   = pline_to_cell_center  [0];
  double      *_pline_to_cell_velocity = pline_to_cell_velocity[0];

  // PDM_log_trace_connectivity_long(_gnum1_come_from_idx, _gnum1_come_from, n_ref_b[0], "_gnum1_come_from ::");
  assert(n_lines == n_ref_b[0]);

  double *pseudo_distance = malloc(    _gnum1_come_from_idx[n_lines] * sizeof(double));
  double *pseudo_coords   = malloc(3 * _gnum1_come_from_idx[n_lines] * sizeof(double));
  int    *order_by_dist   = malloc(3 * _gnum1_come_from_idx[n_lines] * sizeof(int   ));

  double *dline_data = malloc(n_lines * sizeof(double));

  for(int idx_line = 0; idx_line < n_ref_b[0]; ++idx_line) {
    int i_line = ref_b[0][idx_line] - 1;

    int beg = _gnum1_come_from_idx[idx_line];
    int n_cell_connect = _gnum1_come_from_idx[idx_line+1] - beg;
    int         *_order_by_dist     = &order_by_dist          [    beg];
    PDM_g_num_t *_cell_g_num        = &_gnum1_come_from       [    beg];
    double      *_cell_center_coord = &_pline_to_cell_center  [3 * beg];
    double      *_pseudo_coords     = &pseudo_coords          [3 * beg];
    double      *_pseudo_distance   = &pseudo_distance        [    beg];
    double      *_velocity          = &_pline_to_cell_velocity[    beg];

    /* Compute distance */
    for(int j = 0; j < n_cell_connect; ++j ) {
      _order_by_dist[j] = j;

      double nx = psurf_face_normal[3*i_line  ];
      double ny = psurf_face_normal[3*i_line+1];
      double nz = psurf_face_normal[3*i_line+2];

      double sn = sqrt(nx * nx + ny * ny + nz * nz);
      nx = nx/sn;
      ny = ny/sn;
      nz = nz/sn;
      double proj = nx *  _cell_center_coord[3*j  ]
                  + ny *  _cell_center_coord[3*j+1]
                  + nz *  _cell_center_coord[3*j+2];

      _pseudo_coords[3*j  ] = psurf_face_center[3*i_line  ] + proj * nx;
      _pseudo_coords[3*j+1] = psurf_face_center[3*i_line+1] + proj * ny;
      _pseudo_coords[3*j+2] = psurf_face_center[3*i_line+2] + proj * nz;

      double dx = psurf_face_center[3*i_line  ] - _pseudo_coords[3*j  ];
      double dy = psurf_face_center[3*i_line+1] - _pseudo_coords[3*j+1];
      double dz = psurf_face_center[3*i_line+2] - _pseudo_coords[3*j+2];

      _pseudo_distance[j] = sqrt(dx * dx + dy * dy + dz * dz);

    }

    /*
     * Indirect sort by distance
     */
    PDM_sort_double(_pseudo_distance, _order_by_dist, n_cell_connect);

    dline_data[i_line] = 0.; //pequi_surf_parent_face_ln_to_gn[i_line];

    /* Fake integral */
    for(int j = 0; j < n_cell_connect; ++j ) {
      if(_velocity[j] < 1.){
        dline_data[i_line] += _velocity[j];
      }
    }

    if(post) {
      char filename[999];
      sprintf(filename, "line_to_cell_vtx_coords_%i_%i.vtk", i_rank, i_line);
      PDM_vtk_write_point_cloud(filename,
                                n_cell_connect,
                                _pseudo_coords,
                                _cell_g_num,
                                NULL);
    }

  }


  /*
   * Renvoi vers les super cellules pour faire un super modele de plus
   */
  int stride_one = 1;
  double **pline_data = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         dline_data,
                         NULL,
          (void ***)     &pline_data);
  free(dline_data);


  for(int i_part = 0; i_part < n_part_line; ++i_part) {
    if(post) {
      char filename[999];
      sprintf(filename, "pline_data_%3.3d_%3.3d.vtk", i_part, i_rank);

      const char   *vector_field_name[1] = {"pline_data"};
      const double *vector_field     [1] = {pline_data[i_part]};
      PDM_vtk_write_point_cloud_with_field(filename,
                                           pn_cell       [i_part],
                                           cell_center   [i_part],
                                           pcell_ln_to_gn[i_part],
                                           NULL,
                                           1,
                                           vector_field_name,
                                           vector_field,
                                           0,
                                           NULL,
                                           NULL,
                                           0,
                                           NULL,
                                           NULL );
    }
  }



  PDM_block_to_part_free(btp);

  for(int i_part = 0; i_part < n_part_line; ++i_part) {
    free(pline_to_cell_center  [i_part]);
    free(pline_to_cell_velocity[i_part]);
    free(pline_data[i_part]);
  }
  free(pline_to_cell_center);
  free(pline_to_cell_velocity);
  free(pseudo_distance);
  free(pseudo_coords  );
  free(order_by_dist  );
  free(closest_elt_gnum);
  free(pline_data);


  PDM_mesh_intersection_free(mi);

  PDM_dist_cloud_surf_free(dist);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(psurf_vtx_coord    [i_part]);
    free(psurf_face_vtx_idx [i_part]);
    free(psurf_face_vtx     [i_part]);
    free(psurf_face_ln_to_gn[i_part]);
    free(psurf_vtx_ln_to_gn [i_part]);
    free(cell_center        [i_part]);
    free(velocity           [i_part]);
  }
  free(psurf_vtx );
  free(psurf_face);
  free(psurf_vtx_coord);
  free(psurf_face_vtx_idx);
  free(psurf_face_vtx);
  free(psurf_face_ln_to_gn);
  free(psurf_vtx_ln_to_gn);
  free(cell_center);
  free(pcell_ln_to_gn);
  free(pn_cell);
  free(velocity);
  free(pray_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(psurf_face_normal);
  free(psurf_face_center);
  free(ray_coord);
  free(pray_vtx);

  free(pequi_surf_face_vtx_idx        );
  free(pequi_surf_face_vtx            );
  free(pequi_surf_face_ln_to_gn       );
  free(pequi_surf_parent_face_ln_to_gn);
  free(pequi_surf_vtx_ln_to_gn        );
  free(pequi_surf_vtx_coord           );


  PDM_DMesh_nodal_free(dmn_vol_a);
  PDM_multipart_free(mpart_vol_a);
  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize ();

  return 0;

}
