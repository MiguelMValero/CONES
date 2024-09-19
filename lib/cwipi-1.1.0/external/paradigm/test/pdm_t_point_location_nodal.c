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
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_part_geom.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_array.h"
#include "pdm_point_location.h"
#include "pdm_dmesh.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_distrib.h"
#include "pdm_reader_gamma.h"

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
           char                 **filename)
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
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_gen_mesh
(
 const PDM_MPI_Comm                   comm,
 const int                            n_part,
 const PDM_Mesh_nodal_elt_t           t_elt,
 const PDM_g_num_t                    n_vtx_seg,
 const double                         length,
 const PDM_split_dual_t               part_method,
       PDM_part_mesh_nodal_elmts_t  **pmne,
       int                          **pn_elt,
       PDM_g_num_t                 ***pelt_ln_to_gn,
       int                          **pn_vtx,
       double                      ***pvtx_coord
 )
{
  const int randomize = 0;

  *pn_elt        = malloc(sizeof(int          ) * n_part);
  *pelt_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);
  *pn_vtx        = malloc(sizeof(int          ) * n_part);
  *pvtx_coord    = malloc(sizeof(double      *) * n_part);

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

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  if (t_elt == PDM_MESH_NODAL_POLY_3D) {
    PDM_g_num_t  gn_cell = 0;
    PDM_g_num_t  gn_face = 0;
    PDM_g_num_t  gn_vtx  = 0;
    int          dn_cell = 0;
    int          dn_face = 0;
    int          dn_edge = 0;
    int          dn_vtx  = 0;
    double      *dvtx_coord     = NULL;
    int         *dcell_face_idx = NULL;
    PDM_g_num_t *dcell_face     = NULL;
    int          n_face_group    = 0;
    PDM_g_num_t *dface_cell      = NULL;
    int         *dface_vtx_idx   = NULL;
    PDM_g_num_t *dface_vtx       = NULL;
    int         *dface_group_idx = NULL;
    PDM_g_num_t *dface_group     = NULL;

    PDM_poly_vol_gen(comm,
                     0.,
                     0.,
                     0.,
                     length,
                     length,
                     length,
                     n_vtx_seg,
                     n_vtx_seg,
                     n_vtx_seg,
                     randomize,
                     0,
                     &gn_cell,
                     &gn_face,
                     &gn_vtx,
                     &n_face_group,
                     &dn_cell,
                     &dn_face,
                     &dn_vtx,
                     &dcell_face_idx,
                     &dcell_face,
                     &dface_cell,
                     &dface_vtx_idx,
                     &dface_vtx,
                     &dvtx_coord,
                     &dface_group_idx,
                     &dface_group);

    /* Generate dmesh */
    PDM_dmesh_t *dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
                                           dn_cell,
                                           dn_face,
                                           dn_edge,
                                           dn_vtx,
                                           comm);

    PDM_dmesh_vtx_coord_set(dmesh,
                            dvtx_coord,
                            PDM_OWNERSHIP_USER);


    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               dface_vtx,
                               dface_vtx_idx,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               dface_cell,
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_FACE,
                        n_face_group,
                        dface_group,
                        dface_group_idx,
                        PDM_OWNERSHIP_USER);

    PDM_multipart_dmesh_set (mpart, 0, dmesh);

    PDM_multipart_compute(mpart);

    PDM_dmesh_free(dmesh);
    free(dvtx_coord);
    free(dcell_face_idx);
    free(dcell_face);
    free(dface_cell);
    free(dface_vtx_idx);
    free(dface_vtx);
    free(dface_group_idx);
    free(dface_group);


    /* Get parts */
    int         **pcell_face_idx = malloc(sizeof(int         *) * n_part);
    int         **pcell_face     = malloc(sizeof(int         *) * n_part);
    int         **pface_vtx_idx  = malloc(sizeof(int         *) * n_part);
    int         **pface_vtx      = malloc(sizeof(int         *) * n_part);
    PDM_g_num_t **pface_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);

    int  *pn_face        = malloc(sizeof(int) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {


      double *_vtx_coord;
      (*pn_vtx)[ipart] = PDM_multipart_part_vtx_coord_get(mpart,
                                                          0,
                                                          ipart,
                                                          &_vtx_coord,
                                                          PDM_OWNERSHIP_USER);
      (*pvtx_coord)[ipart] = _vtx_coord;

      PDM_g_num_t *_elt_ln_to_gn;
      (*pn_elt)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                         0,
                                                         ipart,
                                                         PDM_MESH_ENTITY_CELL,
                                                         &_elt_ln_to_gn,
                                                         PDM_OWNERSHIP_USER);
      (*pelt_ln_to_gn)[ipart] = _elt_ln_to_gn;

      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      ipart,
                                      PDM_MESH_ENTITY_FACE,
                                      &pface_ln_to_gn[ipart],
                                      PDM_OWNERSHIP_KEEP);

      int *_face_vtx;
      int *_face_vtx_idx;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                       &_face_vtx_idx,
                                                       &_face_vtx,
                                                       PDM_OWNERSHIP_KEEP);
      pn_face[ipart] = n_face;
      int *_cell_face;
      int *_cell_face_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &_cell_face_idx,
                                          &_cell_face,
                                          PDM_OWNERSHIP_KEEP);

      pcell_face_idx[ipart] = _cell_face_idx;
      pcell_face    [ipart] = _cell_face;
      pface_vtx_idx [ipart] = _face_vtx_idx;
      pface_vtx     [ipart] = _face_vtx;
    }

    *pmne = PDM_part_mesh_nodal_elmts_create_from_part3d(n_part,
                             (const int               *) (*pn_elt),
                             (const int               *) pn_face,
                             (const int              **) pface_vtx_idx,
                             (const int              **) pface_vtx,
                             (const PDM_g_num_t      **) pface_ln_to_gn,
                             (const int              **) pcell_face_idx,
                             (const int              **) pcell_face,
                             (const double           **) (*pvtx_coord),
                             (const PDM_g_num_t      **) (*pelt_ln_to_gn),
                                                         comm);
    free(pn_face       );
    free(pface_vtx_idx );
    free(pface_vtx     );
    free(pcell_face_idx);
    free(pcell_face    );
    free(pface_ln_to_gn);
  }

  else if (t_elt == PDM_MESH_NODAL_POLY_2D) {

    PDM_g_num_t     gn_face         = 0;
    PDM_g_num_t     gn_vtx          = 0;
    PDM_g_num_t     gn_edge         = 0;
    int             dn_vtx          = 0;
    int             dn_face         = 0;
    int             dn_edge         = 0;
    int             n_edge_group    = 0;
    double         *dvtx_coord      = NULL;
    int            *dface_edge_idx  = NULL;
    PDM_g_num_t    *dface_vtx       = NULL;
    PDM_g_num_t    *dface_edge      = NULL;
    int            *dedge_vtx_idx   = NULL;
    PDM_g_num_t    *dedge_vtx       = NULL;
    PDM_g_num_t    *dedge_face      = NULL;
    int            *dedge_group_idx = NULL;
    PDM_g_num_t    *dedge_group     = NULL;

    PDM_poly_surf_gen(comm,
                      0,
                      length,
                      0,
                      length,
                      randomize,
                      0,
                      n_vtx_seg,
                      n_vtx_seg,
                      &gn_face,
                      &gn_vtx,
                      &gn_edge,
                      &dn_vtx,
                      &dvtx_coord,
                      &dn_face,
                      &dface_edge_idx,
                      &dface_vtx,
                      &dface_edge,
                      &dn_edge,
                      &dedge_vtx,
                      &dedge_face,
                      &n_edge_group,
                      &dedge_group_idx,
                      &dedge_group);

    dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

    int n_bound = n_edge_group;

    int *dedge_bnd_idx  = (int *) malloc((n_bound + 1) * sizeof(int));
    dedge_bnd_idx[0] = 0;

    // First pass to count and allocate
    int i_bnd = 1;
    for (int igroup = 0; igroup < n_edge_group; igroup++) {
      int group_size = dedge_group_idx[igroup+1] - dedge_group_idx[igroup];
      dedge_bnd_idx[i_bnd++] = group_size;
    }
    for (int i = 0; i < n_bound; i++) {
      dedge_bnd_idx[i+1] = dedge_bnd_idx[i+1] + dedge_bnd_idx[i];
    }

    // Second pass to copy
    PDM_g_num_t *dedge_bnd  = (PDM_g_num_t *) malloc(dedge_bnd_idx[n_bound] * sizeof(PDM_g_num_t));

    i_bnd = 0;
    for (int igroup = 0; igroup < n_edge_group; igroup++) {
      for (int i = dedge_group_idx[igroup]; i < dedge_group_idx[igroup+1]; i++) {
        dedge_bnd[i_bnd++] = dedge_group[i];
      }
    }
    PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                          0,
                                          dn_face,
                                          dn_edge,
                                          dn_vtx,
                                          comm);

    PDM_dmesh_vtx_coord_set(dmesh,
                            dvtx_coord,
                            PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               dedge_vtx,
                               dedge_vtx_idx,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                               dedge_face,
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_EDGE,
                        n_edge_group,
                        dedge_bnd,
                        dedge_bnd_idx,
                        PDM_OWNERSHIP_USER);

    PDM_multipart_dmesh_set(mpart, 0, dmesh);

    /* Run */
    PDM_multipart_compute(mpart);

    PDM_dmesh_free(dmesh);
    free(dedge_bnd_idx   );
    free(dedge_bnd       );
    free(dvtx_coord      );
    free(dface_edge_idx  );
    free(dface_vtx       );
    free(dface_edge      );
    free(dedge_vtx_idx   );
    free(dedge_vtx       );
    free(dedge_face      );
    free(dedge_group_idx );
    free(dedge_group     );

    /* Get parts */
    int  *pn_edge        = malloc(sizeof(int  ) * n_part);
    int **pface_edge_idx = malloc(sizeof(int *) * n_part);
    int **pface_edge     = malloc(sizeof(int *) * n_part);
    int **pedge_vtx_idx  = malloc(sizeof(int *) * n_part);
    int **pedge_vtx      = malloc(sizeof(int *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      double *_vtx_coord;
      (*pn_vtx)[ipart] = PDM_multipart_part_vtx_coord_get(mpart,
                                                          0,
                                                          ipart,
                                                          &_vtx_coord,
                                                          PDM_OWNERSHIP_USER);
      (*pvtx_coord)[ipart] = _vtx_coord;

      PDM_g_num_t *_elt_ln_to_gn;
      (*pn_elt)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                         0,
                                                         ipart,
                                                         PDM_MESH_ENTITY_CELL,
                                                         &_elt_ln_to_gn,
                                                         PDM_OWNERSHIP_USER);
      (*pelt_ln_to_gn)[ipart] = _elt_ln_to_gn;

      int *_face_vtx;
      int *_face_vtx_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &_face_vtx_idx,
                                          &_face_vtx,
                                          PDM_OWNERSHIP_KEEP);

      int *_face_edge;
      int *_face_edge_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &_face_edge_idx,
                                          &_face_edge,
                                          PDM_OWNERSHIP_KEEP);

      pface_edge_idx[ipart] = _face_edge_idx;
      pface_edge    [ipart] = _face_edge;

      if (_face_edge != NULL) {

        int *_edge_vtx;
        int *_edge_vtx_idx;
        int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         ipart,
                                                         PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                         &_edge_vtx_idx,
                                                         &_edge_vtx,
                                                         PDM_OWNERSHIP_KEEP);
        pn_edge[ipart] = n_edge;
        pedge_vtx_idx[ipart] = _edge_vtx_idx;
        pedge_vtx    [ipart] = _edge_vtx;
      }
      else {
        abort();
        assert(_face_vtx != NULL);
      }
    }
     *pmne = PDM_part_mesh_nodal_elmts_create_from_part2d(n_part,
                                   (const int          *) (*pn_elt),
                                   (const int          *) pn_edge,
                                   (const int          *) (*pn_vtx),
                                   (const int         **) pedge_vtx_idx,
                                   (const int         **) pedge_vtx,
                                   (const int         **) pface_edge_idx,
                                   (const int         **) pface_edge,
                                   (const PDM_g_num_t **) (*pelt_ln_to_gn),
                                                          comm);
     free(pn_edge       );
     free(pedge_vtx_idx );
     free(pedge_vtx     );
     free(pface_edge_idx);
     free(pface_edge    );

  }

  else {

    PDM_g_num_t n_vtx_seg_z = n_vtx_seg;
    PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
    PDM_mesh_entities_t entity_type = PDM_MESH_ENTITY_CELL;
    if (PDM_Mesh_nodal_elt_dim_get(t_elt) == 2) {
      n_vtx_seg_z = 1;
      geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
      entity_type = PDM_MESH_ENTITY_FACE;
    }

    PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                          n_vtx_seg,
                                                          n_vtx_seg,
                                                          n_vtx_seg_z,
                                                          length,
                                                          0.,
                                                          0.,
                                                          0.,
                                                          t_elt,
                                                          1,
                                                          PDM_OWNERSHIP_KEEP);

    PDM_dcube_nodal_gen_build(dcube);

    PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);


    /* Add noise */
    if (1) {
      int i_rank;
      PDM_MPI_Comm_rank(comm, &i_rank);

      PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
      double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
      int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];
      double noise = 0.2 * length / (double) (n_vtx_seg - 1);
      for (int i = 0; i < dn_vtx; i++) {
        for (int j = 0; j < 3; j++) {
          dvtx_coord[3*i+j] += noise * (rand() / (double) RAND_MAX - 0.5);
        }
      }
    }


    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

    PDM_multipart_compute(mpart);

    PDM_g_num_t **pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      PDM_g_num_t *_elt_ln_to_gn;
      (*pn_elt)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                         0,
                                                         ipart,
                                                         entity_type,
                                                         &_elt_ln_to_gn,
                                                         PDM_OWNERSHIP_USER);
      (*pelt_ln_to_gn)[ipart] = _elt_ln_to_gn;

      (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                         0,
                                                         ipart,
                                                         PDM_MESH_ENTITY_VTX,
                                                         &pvtx_ln_to_gn[ipart],
                                                         PDM_OWNERSHIP_KEEP);

      double *_vtx_coord;
      PDM_multipart_part_vtx_coord_get(mpart,
                                       0,
                                       ipart,
                                       &_vtx_coord,
                                       PDM_OWNERSHIP_USER);
      (*pvtx_coord)[ipart] = _vtx_coord;
    }


    *pmne = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn,
                                                     geom_kind,
                                                     n_part,
                                                     *pn_vtx,
                                                     pvtx_ln_to_gn,
                                                     *pn_elt,
                                                     *pelt_ln_to_gn,
                                                     NULL);
    free(pvtx_ln_to_gn);

    PDM_dcube_nodal_gen_free(dcube);

  }
  PDM_multipart_free(mpart);
}



static void
_mesh_from_file
(
 const PDM_MPI_Comm                   comm,
 const int                            n_part,
 const char                          *filename,
 const PDM_split_dual_t               part_method,
       PDM_part_mesh_nodal_elmts_t  **pmne,
       int                          **pn_elt,
       PDM_g_num_t                 ***pelt_ln_to_gn,
       int                          **pn_vtx,
       double                      ***pvtx_coord
 )
{
  PDM_dmesh_nodal_t *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                                        filename,
                                                        0,
                                                        0);

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

  *pn_elt        = malloc(sizeof(int          ) * n_part);
  *pelt_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);
  *pn_vtx        = malloc(sizeof(int          ) * n_part);
  *pvtx_coord    = malloc(sizeof(double      *) * n_part);

  PDM_g_num_t **pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_g_num_t *_elt_ln_to_gn;
    (*pn_elt)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_CELL,
                                                       &_elt_ln_to_gn,
                                                       PDM_OWNERSHIP_USER);
    (*pelt_ln_to_gn)[ipart] = malloc(sizeof(PDM_g_num_t) * (*pn_elt)[ipart]);
    memcpy((*pelt_ln_to_gn)[ipart], _elt_ln_to_gn, sizeof(PDM_g_num_t) * (*pn_elt)[ipart]);

    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VTX,
                                                       &pvtx_ln_to_gn[ipart],
                                                       PDM_OWNERSHIP_USER);

    double *_vtx_coord;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     0,
                                     ipart,
                                     &_vtx_coord,
                                     PDM_OWNERSHIP_USER);
    (*pvtx_coord)[ipart] = malloc(sizeof(double) * (*pn_vtx)[ipart] * 3);
    memcpy((*pvtx_coord)[ipart], _vtx_coord, sizeof(double) * (*pn_vtx)[ipart] * 3);
  }


  *pmne = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn,
                                                   PDM_GEOMETRY_KIND_VOLUMIC,
                                                   n_part,
                                                   *pn_vtx,
                                                   pvtx_ln_to_gn,
                                                   *pn_elt,
                                                   *pelt_ln_to_gn,
                                                   NULL);
  free(pvtx_ln_to_gn);

  PDM_multipart_free(mpart);
  PDM_DMesh_nodal_free(dmn);
}


static void
_compute_cell_centers
(
 PDM_part_mesh_nodal_elmts_t   *pmne,
 int                            n_part,
 double                       **pvtx_coord,
 int                          **pn_pts,
 double                      ***pts_coord
 )
{
  *pn_pts    = malloc(sizeof(int     ) * n_part);
  *pts_coord = malloc(sizeof(double *) * n_part);

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int ipart = 0; ipart < n_part; ipart++) {

    int pn_elt = 0;
    for (int isection = 0; isection < n_section; isection++) {
      int id_section = sections_id[isection];
      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            id_section,
                                                            ipart);

      pn_elt += n_elt;
    }

    (*pn_pts)   [ipart] = pn_elt;
    (*pts_coord)[ipart] = malloc(sizeof(double) * pn_elt * 3);

    for (int isection = 0; isection < n_section; isection++) {

      int id_section = sections_id[isection];

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            id_section,
                                                            ipart);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                            id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 ipart,
                                                                 PDM_OWNERSHIP_KEEP);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {

        /* Polygonal section */
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec_idx,
                                                   &connec,
                                                   PDM_OWNERSHIP_KEEP);
        // PDM_log_trace_connectivity_int(connec_idx,
        //                                connec,
        //                                n_elt,
        //                                "connec : ");

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int n_vtx = connec_idx[ielt+1] - connec_idx[ielt];

          double *pc = (*pts_coord)[ipart] + 3*icell;
          pc[0] = pc[1] = pc[2] = 0.;

          for (int idx = connec_idx[ielt]; idx < connec_idx[ielt+1]; idx++) {
            int vtx_id = connec[idx] - 1;
            for (int j = 0; j < 3; j++) {
              pc[j] += pvtx_coord[ipart][3*vtx_id + j];
            }
          }

          double in = 1./(double) n_vtx;
          for (int j = 0; j < 3; j++) {
            pc[j] *= in;
          }

        } // End of loop on polygons

      } // End if PDM_MESH_NODAL_POLY_2D

      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {

        /* Polyhedral section */
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get(pmne,
                                                                    id_section,
                                                                    ipart,
                                                                    &connec_idx,
                                                                    &connec,
                                                                    PDM_OWNERSHIP_KEEP);


        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int n_vtx = connec_idx[ielt+1] - connec_idx[ielt];

          double *pc = (*pts_coord)[ipart] + 3*icell;
          pc[0] = pc[1] = pc[2] = 0.;

          for (int idx = connec_idx[ielt]; idx < connec_idx[ielt+1]; idx++) {
            int vtx_id = connec[idx] - 1;
            for (int j = 0; j < 3; j++) {
              pc[j] += pvtx_coord[ipart][3*vtx_id + j];
            }
          }

          double in = 1./(double) n_vtx;
          for (int j = 0; j < 3; j++) {
            pc[j] *= in;
          }

        } // End of loop on polyhedra

      } // End if PDM_MESH_NODAL_POLY_3D


      else {

        /* Standard section */
              int         *connec              = NULL;
              PDM_g_num_t *numabs              = NULL;
              int         *_parent_num         = NULL;
              PDM_g_num_t *parent_entity_g_num = NULL;
              int          order               = 0;
        const char        *ho_ordering         = NULL;
        PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec,
                                                   &numabs,
                                                   &_parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_KEEP);

        int n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                 order);
        // double in = 1./(double) n_vtx;

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          double *pc = (*pts_coord)[ipart] + 3*icell;
          pc[0] = pc[1] = pc[2] = 0.;
          double sum_w = 0.;

          for (int idx = n_vtx*ielt; idx < n_vtx*(ielt+1); idx++) {
            int vtx_id = connec[idx] - 1;
            double w = (double) rand() / (double) RAND_MAX;
            // double w = 1.;
            sum_w += w;
            for (int j = 0; j < 3; j++) {
              pc[j] += w * pvtx_coord[ipart][3*vtx_id + j];
            }
          }

          for (int j = 0; j < 3; j++) {
            pc[j] /= sum_w;
          }

        } // End of loop on elts

      } // End if std elt

    } // End of loop on sections
  } // End of loop on parts
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

  PDM_g_num_t          n_vtx_seg   = 5;
  double               length      = 1.;
  int                  n_part      = 1;
  int                  post        = 0;
  PDM_Mesh_nodal_elt_t t_elt       = PDM_MESH_NODAL_HEXA8;
  double               tolerance   = 1e-6;
  PDM_split_dual_t     part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  char       *filename_mesh = NULL;

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
             &t_elt,
             &filename_mesh);

  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /* Generate a part_mesh_nodal */
  PDM_part_mesh_nodal_elmts_t *pmne;
  int          *pn_elt        = NULL;
  PDM_g_num_t **pelt_ln_to_gn = NULL;
  int          *pn_vtx        = NULL;
  double      **pvtx_coord    = NULL;
  if (filename_mesh == NULL) {
    _gen_mesh(comm,
              n_part,
              t_elt,
              n_vtx_seg,
              length,
              part_method,
              &pmne,
              &pn_elt,
              &pelt_ln_to_gn,
              &pn_vtx,
              &pvtx_coord);
  }
  else {
    _mesh_from_file(comm,
                    n_part,
                    filename_mesh,
                    part_method,
                    &pmne,
                    &pn_elt,
                    &pelt_ln_to_gn,
                    &pn_vtx,
                    &pvtx_coord);
  }

  if (post) {
    char filename[999];

    int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
    int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

    for (int ipart = 0; ipart < n_part; ipart++) {

      for (int isection = 0; isection < n_section; isection++) {

        sprintf(filename, "pmne_part_%d_section_%d_%3.3d.vtk",
                ipart, isection, i_rank);

        int id_section = sections_id[isection];
        PDM_Mesh_nodal_elt_t _t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                               id_section);

        int _n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                               id_section,
                                                               ipart);

        int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                   id_section,
                                                                   ipart,
                                                                   PDM_OWNERSHIP_KEEP);


        PDM_g_num_t *gnum = malloc(sizeof(PDM_g_num_t) * _n_elt);
        for (int i = 0; i < _n_elt; i++) {
          int icell = i;
          if (parent_num != NULL) {
            icell = parent_num[i];
          }
          gnum[i] = pelt_ln_to_gn[ipart][icell];
        }


        if (_t_elt == PDM_MESH_NODAL_POLY_2D) {

          int *connec_idx;
          int *connec;
          PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                     id_section,
                                                     ipart,
                                                     &connec_idx,
                                                     &connec,
                                                     PDM_OWNERSHIP_KEEP);

          PDM_vtk_write_polydata(filename,
                                 pn_vtx[ipart],
                                 pvtx_coord[ipart],
                                 NULL,
                                 _n_elt,
                                 connec_idx,
                                 connec,
                                 gnum,
                                 NULL);

        }

        else if (_t_elt == PDM_MESH_NODAL_POLY_3D) {

          int  n_face;
          int *face_vtx_idx;
          int *face_vtx;
          int *cell_face_idx;
          int *cell_face;
          PDM_g_num_t *face_ln_to_gn            = NULL;
          int         *_parent_num              = NULL;
          PDM_g_num_t *numabs                   = NULL;
          PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
          PDM_part_mesh_nodal_elmts_section_poly3d_get(pmne,
                                                     id_section,
                                                     ipart,
                                                     &n_face,
                                                     &face_ln_to_gn,
                                                     &face_vtx_idx,
                                                     &face_vtx,
                                                     &numabs,
                                                     &cell_face_idx,
                                                     &cell_face,
                                                     &_parent_num,
                                                     &parent_entitity_ln_to_gn,
                                                     PDM_OWNERSHIP_KEEP);

          PDM_vtk_write_polydata(filename,
                                 pn_vtx[ipart],
                                 pvtx_coord[ipart],
                                 NULL,
                                 n_face,
                                 face_vtx_idx,
                                 face_vtx,
                                 face_ln_to_gn,
                                 NULL);

        }

        else {
          int         *elmt_vtx                 = NULL;
          int         *_parent_num              = NULL;
          PDM_g_num_t *numabs                   = NULL;
          PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
          PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                  id_section,
                                                  ipart,
                                                  &elmt_vtx,
                                                  &numabs,
                                                  &_parent_num,
                                                  &parent_entitity_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);


          PDM_vtk_write_std_elements(filename,
                                     pn_vtx[ipart],
                                     pvtx_coord[ipart],
                                     NULL,
                                     _t_elt,
                                     _n_elt,
                                     elmt_vtx,
                                     gnum,//parent_entitity_ln_to_gn,
                                     0,
                                     NULL,
                                     NULL);
        }

        free(gnum);
      }
    }
  }

  /* Generate point cloud (cell centers) */
  int     *pn_pts    = NULL;
  double **pts_coord = NULL;
  _compute_cell_centers(pmne,
                        n_part,
                        pvtx_coord,
                        &pn_pts,
                        &pts_coord);

  int **pts_idx = malloc(sizeof(int *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    pts_idx[ipart] = PDM_array_new_idx_from_const_stride_int(1,
                                                             pn_pts[ipart]);
  }


  if (post) {
    char filename[999];

    for (int ipart = 0; ipart < n_part; ipart++) {
      sprintf(filename, "pts_coord_%d_%3.3d.vtk", ipart, i_rank);

      assert(pn_elt[ipart] == pn_pts[ipart]);
      PDM_vtk_write_point_cloud(filename,
                                pts_idx[ipart][pn_elt[ipart]],
                                pts_coord[ipart],
                                pelt_ln_to_gn[ipart],
                                NULL);
    }
  }




  /*
   *  Location
   */
  double **distance        = NULL;
  double **projected_coord = NULL;
  int    **bar_coord_idx   = NULL;
  double **bar_coord       = NULL;
  double **uvw             = NULL;
  PDM_point_location_nodal(pmne,
                            n_part,
          (const double **) pvtx_coord,
          (const int    **) pts_idx,
          (const double **) pts_coord,
                            tolerance,
                            &distance,
                            &projected_coord,
                            &bar_coord_idx,
                            &bar_coord,
                            &uvw);

  if (post) {
    char filename[999];

    for (int ipart = 0; ipart < n_part; ipart++) {
      sprintf(filename, "proj_coord_%d_%3.3d.vtk", ipart, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                pts_idx[ipart][pn_elt[ipart]],
                                projected_coord[ipart],
                                pelt_ln_to_gn[ipart],
                                NULL);
    }
  }

  // for (int ipart = 0; ipart < n_part; ipart++) {
  //   log_trace("--- part %d ---\n", ipart);
  //   for (int ielt = 0; ielt < pn_elt[ipart]; ielt++) {
  //     log_trace("  elt %d\n", ielt);
  //     log_trace("    distance2 = %e\n", distance[ipart][ielt]);
  //     log_trace("    proj_coord = %f %f %f\n",
  //               projected_coord[ipart][3*ielt  ],
  //               projected_coord[ipart][3*ielt+1],
  //               projected_coord[ipart][3*ielt+2]);
  //     PDM_log_trace_array_double(bar_coord[ipart] + bar_coord_idx[ipart][ielt],
  //                                bar_coord_idx[ipart][ielt+1] - bar_coord_idx[ipart][ielt],
  //                                "    bar_coord : ");
  //   }
  // }



  /* Free memory */
  for (int ipart = 0; ipart < n_part; ipart++) {
    free(pts_idx        [ipart]);
    free(pts_coord      [ipart]);
    free(distance       [ipart]);
    free(projected_coord[ipart]);
    free(bar_coord_idx  [ipart]);
    free(bar_coord      [ipart]);
    free(uvw            [ipart]);
    free(pelt_ln_to_gn  [ipart]);
    free(pvtx_coord     [ipart]);
  }
  free(pts_idx);
  free(pts_coord);
  free(distance);
  free(projected_coord);
  free(bar_coord_idx);
  free(bar_coord);
  free(uvw);
  free(pelt_ln_to_gn);
  free(pvtx_coord);

  free(pn_pts);
  free(pn_vtx);
  free(pn_elt);

  PDM_part_mesh_nodal_elmts_free(pmne);


  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
