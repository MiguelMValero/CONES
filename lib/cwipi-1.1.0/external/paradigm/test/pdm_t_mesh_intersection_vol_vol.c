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
#include "pdm_mesh_intersection_priv.h"
#include "pdm_multipart.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_reader_gamma.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_connectivity_transform.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum
{
  TETRA_POINT,
  TETRA_CENTER,
  USER
} point_t;

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
     "  -nB     <level>  Number vtx in side of mesh B (default : 10).\n\n"
     "  -n_part <level>  Number of partitions         (default : 1).\n\n"
     "  -t               Element type.\n\n"
     "  -tA              Element type of mesh A (default : HEXA8).\n\n"
     "  -tB              Element type of mesh B (default : HEXA8).\n\n"
     "  -tP              Tetraisation point type.\n\n"
     "  -coordP          Tetraisation point coordinates if type is user.\n\n"
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
 double                *noise_a,
 double                *noise_b,
 point_t               *tetraisation_pt_type,
 double               **tetraisation_pt_coord,
 double                *shift_b,
 int                   *verbose,
 char                  *filenames[]
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
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type_a = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
        *elt_type_b = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tA") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type_a = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tB") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type_b = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nodalA") == 0) {
      *nodal_a = 1;
    }
    else if (strcmp(argv[i], "-nodalB") == 0) {
      *nodal_b = 1;
    }
    else if (strcmp(argv[i], "-noiseA") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *noise_a = (double) atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-noiseB") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *noise_b = (double) atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tP") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *tetraisation_pt_type = (point_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-coordP") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        (*tetraisation_pt_coord)[0] = (double) atof(argv[i]);
        i++;
        (*tetraisation_pt_coord)[1] = (double) atof(argv[i]);
        i++;
        (*tetraisation_pt_coord)[2] = (double) atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-shiftB") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        shift_b[0] = (double) atof(argv[i]);
        i++;
        shift_b[1] = (double) atof(argv[i]);
        i++;
        shift_b[2] = (double) atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-verbose") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-fA") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        filenames[0] = argv[i];
      }
    }
    else if (strcmp(argv[i], "-fB") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        filenames[1] = argv[i];
      }
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
 const double                 noise_scale,
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


  PDM_dmesh_nodal_t *dmn = NULL;
  int dn_vtx = 0;
  double      *dvtx_coord       = NULL;
  int         *dcell_face_idx   = NULL;
  PDM_g_num_t *dcell_face       = NULL;
  PDM_g_num_t *dface_cell       = NULL;
  int         *dface_vtx_idx    = NULL;
  PDM_g_num_t *dface_vtx        = NULL;
  int         *dface_group_idx  = NULL;
  PDM_g_num_t *dface_group      = NULL;
  PDM_dmesh_t *dmesh            = NULL;

  if (elt_type < PDM_MESH_NODAL_POLY_3D) {

    PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
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
    PDM_dcube_nodal_gen_build(dcube);
    dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);
    PDM_dcube_nodal_gen_free(dcube);

    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

    PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);
    dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

    if (noise_scale > 0) {
      double eps = 1e-12;
      double step = length / (double) (n_vtx_seg - 1);
      double origin[3] = {xmin, ymin, zmin};

      for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
        for (int j = 0; j < 3; j++) {
          double x = dvtx_coord[3*i_vtx+j];
          if (x > origin[j]+eps && x < origin[j]+length-eps) {
            dvtx_coord[3*i_vtx+j] += noise_scale * step * (2*(double) rand()/(double) RAND_MAX - 1.);
          }
        }
      }
    }
  }

  else {
    // Polyhedral mesh
    PDM_g_num_t ng_cell      = 0;
    PDM_g_num_t ng_face      = 0;
    PDM_g_num_t ng_vtx       = 0;
    int         dn_cell      = 0;
    int         dn_face      = 0;
    int         dn_edge      = 0;
    int         n_face_group = 0;

    PDM_poly_vol_gen(comm,
                     xmin,
                     ymin,
                     zmin,
                     length,
                     length,
                     length,
                     n_vtx_seg,
                     n_vtx_seg,
                     n_vtx_seg,
                     (noise_scale > 0),
                     0,
                     &ng_cell,
                     &ng_face,
                     &ng_vtx,
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
    dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
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

    PDM_multipart_dmesh_set(mpart, 0, dmesh);

  }

  if(rotate) {
    // Do something
    double pi = 4 * atan(1.);
    double angle = pi/5.;
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      _rotate_coord(angle, &dvtx_coord[3*i_vtx]);
    }
  }

  PDM_multipart_compute(mpart);

  if (elt_type == PDM_MESH_NODAL_POLY_3D) {
    PDM_dmesh_free(dmesh);
    free(dvtx_coord);
    free(dcell_face_idx);
    free(dcell_face);
    free(dface_cell);
    free(dface_vtx_idx);
    free(dface_vtx);
    free(dface_group_idx);
    free(dface_group);
  }

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

    int *face_vtx_idx = NULL;
    int *face_vtx     = NULL;

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
                                   face_vtx_idx,
                                   face_vtx,
                                   cell_ln_to_gn,
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

  if (pmn == NULL) {
    int n_part = 1;//!!!!!!

    int          *n_cell        = malloc(sizeof(int          ) * n_part);
    int          *n_face        = malloc(sizeof(int          ) * n_part);
    int         **face_vtx_idx  = malloc(sizeof(int         *) * n_part);
    int         **face_vtx      = malloc(sizeof(int         *) * n_part);
    PDM_g_num_t **face_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);
    int         **cell_face_idx = malloc(sizeof(int         *) * n_part);
    int         **cell_face     = malloc(sizeof(int         *) * n_part);
    double      **vtx_coord     = malloc(sizeof(double      *) * n_part);
    PDM_g_num_t **cell_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);

    pmn = PDM_part_mesh_nodal_create(3, n_part, mi->comm);

    // int i_section = PDM_part_mesh_nodal_section_add(pmn, PDM_MESH_NODAL_POLY_3D);

    for (int ipart = 0; ipart < n_part; ipart++) {

      n_cell[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                      0,
                                                      ipart,
                                                      PDM_MESH_ENTITY_CELL,
                                                      &cell_ln_to_gn[ipart],
                                                      PDM_OWNERSHIP_KEEP);

      n_face[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                      0,
                                                      ipart,
                                                      PDM_MESH_ENTITY_FACE,
                                                      &face_ln_to_gn[ipart],
                                                      PDM_OWNERSHIP_KEEP);
      PDM_g_num_t *vtx_ln_to_gn = NULL;
      int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                  0,
                                                  ipart,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_vtx_coord_get(mpart,
                                       0,
                                       ipart,
                                       &vtx_coord[ipart],
                                       PDM_OWNERSHIP_KEEP);

      PDM_part_mesh_nodal_coord_set(pmn,
                                    ipart,
                                    n_vtx,
                                    vtx_coord[ipart],
                                    vtx_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &cell_face_idx[ipart],
                                          &cell_face    [ipart],
                                          PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &face_vtx_idx[ipart],
                                          &face_vtx    [ipart],
                                          PDM_OWNERSHIP_KEEP);

      if (face_vtx[ipart] == NULL) {
        int *edge_vtx_idx = NULL;
        int *edge_vtx     = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            ipart,
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &edge_vtx_idx,
                                            &edge_vtx,
                                            PDM_OWNERSHIP_KEEP);

        int *face_edge = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            ipart,
                                            PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                            &face_vtx_idx[ipart],
                                            &face_edge,
                                            PDM_OWNERSHIP_KEEP);

        PDM_compute_face_vtx_from_face_and_edge(n_face[ipart],
                                                face_vtx_idx[ipart],
                                                face_edge,
                                                edge_vtx,
                                                &face_vtx_idx[ipart]);
      }

      PDM_part_mesh_nodal_cell3d_cellface_add(pmn,
                                              ipart,
                                              n_cell       [ipart],
                                              n_face       [ipart],
                                              face_vtx_idx [ipart],
                                              face_vtx     [ipart],
                                              face_ln_to_gn[ipart],
                                              cell_face_idx[ipart],
                                              cell_face    [ipart],
                                              cell_ln_to_gn[ipart],
                                              PDM_OWNERSHIP_KEEP);
      // PDM_part_mesh_nodal_section_poly3d_set(pmn,
      //                                        i_section,
      //                                        ipart,
      //                                        n_cell       [ipart],
      //                                        n_face       [ipart],
      //                                        face_vtx_idx [ipart],
      //                                        face_vtx     [ipart],
      //                                        face_ln_to_gn[ipart],
      //                                        cell_face_idx[ipart],
      //                                        cell_face    [ipart],
      //                                        cell_ln_to_gn[ipart],
      //                                        NULL,
      //                                        NULL,
      //                                        PDM_OWNERSHIP_KEEP);
    }

    // PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_elmts_create_from_part3d(n_part,
    //                                                           (const int          *) n_cell,
    //                                                           (const int          *) n_face,
    //                                                           (const int         **) face_vtx_idx,
    //                                                           (const int         **) face_vtx,
    //                                                           (const PDM_g_num_t **) face_ln_to_gn,
    //                                                           (const int         **) cell_face_idx,
    //                                                           (const int         **) cell_face,
    //                                                           (const double      **) vtx_coord,
    //                                                           (const PDM_g_num_t **) cell_ln_to_gn,
    //                                                                                  mi->comm);

    // PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmne, PDM_OWNERSHIP_KEEP);

    free(n_cell       );
    free(n_face       );
    free(face_vtx_idx );
    free(face_vtx     );
    free(face_ln_to_gn);
    free(cell_face_idx);
    free(cell_face    );
    free(vtx_coord    );
    free(cell_ln_to_gn);
  }

  // if (1) {
  //   /* Fix orientation if necessary */
  //   PDM_part_mesh_nodal_elmts_t *pmne = pmn->volumic;

  //   int n_part = PDM_part_mesh_nodal_n_part_get(pmn);

  //   int n_section = PDM_part_mesh_nodal_elmts_n_section_get(pmne);

  //   int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  //   for (int ipart = 0; ipart < n_part; ipart++) {

  //     double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, ipart);

  //     for (int isection = 0; isection < n_section; isection++) {
  //       int id_section = sections_id[isection];

  //       PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
  //                                                                               id_section);

  //       int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
  //                                                               id_section,
  //                                                               ipart);

  //       double *volume = malloc(sizeof(double) * n_elt);
  //       double *center = malloc(sizeof(double) * n_elt * 3);

  //       switch (t_elt) {

  //         case PDM_MESH_NODAL_PYRAMID5: {
  //           int         *connec;
  //           PDM_g_num_t *numabs;
  //           int         *parent_num;
  //           PDM_g_num_t *parent_entity_g_num;
  //           PDM_part_mesh_nodal_elmts_section_std_get(pmne,
  //                                                     id_section,
  //                                                     ipart,
  //                                                     &connec,
  //                                                     &numabs,
  //                                                     &parent_num,
  //                                                     &parent_entity_g_num);

  //           PDM_geom_e

  //           break;
  //         }

  //         default:
  //           break;
  //       }

  //     }
  //   }
  // }



  if (1 == 1) {
    char filename[999];
    sprintf(filename, "check_pmn_%d", i_mesh);
    PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_VOLUMIC, filename);
  }

  PDM_mesh_intersection_mesh_nodal_set(mi, i_mesh, pmn);

  return pmn;
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

static
void
_read_volume_mesh
(
 const PDM_MPI_Comm           comm,
 const char                  *filename,
 const int                    rotate,
 const PDM_split_dual_t       part_method,
 const int                    n_part,
       PDM_dmesh_nodal_t    **dmn,
       PDM_multipart_t      **mpart
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);


  const char *file_extension = _get_filename_extension(filename);

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
    *dmn = PDM_vtk_read_to_dmesh_nodal(comm,
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
    *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                        filename,
                                        0,
                                        0);
  }

  if (rotate) {
    double pi = 4 * atan(1.);
    double angle = pi/5.;

    PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(*dmn);
    double *dvtx_coord = PDM_DMesh_nodal_vtx_get(*dmn);
    int     dn_vtx     = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      _rotate_coord(angle, &dvtx_coord[3*i_vtx]);
    }
  }

  int n_domain = 1;
  *mpart = PDM_multipart_create(n_domain,
                                &n_part,
                                PDM_FALSE,
                                part_method,
                                PDM_PART_SIZE_HOMOGENEOUS,
                                NULL,
                                comm,
                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(*mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(*mpart, 0, *dmn);
  PDM_multipart_compute(*mpart);
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
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t          n_vtx_a               = 10;
  PDM_g_num_t          n_vtx_b               = 10;
  PDM_Mesh_nodal_elt_t elt_type_a            = PDM_MESH_NODAL_HEXA8;
  PDM_Mesh_nodal_elt_t elt_type_b            = PDM_MESH_NODAL_HEXA8;
  int                  nodal_a               = 0;
  int                  nodal_b               = 0;
  double               noise_a               = 0;
  double               noise_b               = 0;
  point_t              tetraisation_pt_type  = TETRA_POINT;
  double              *tetraisation_pt_coord = malloc(sizeof(double) * 3);
  double               shift_b[3]            = {0.5, 0.5, 0.5};
  char                *filenames[2]          = {NULL, NULL};

  PDM_split_dual_t     part_method           = PDM_SPLIT_DUAL_WITH_HILBERT;

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
             &noise_a,
             &noise_b,
             &tetraisation_pt_type,
             &tetraisation_pt_coord,
             shift_b,
             &verbose,
             filenames);

  /*
   * Generate meshA
   */
  double length_a = 1.;
  int rotate_a = 0;
  PDM_dmesh_nodal_t *dmn_a   = NULL;
  PDM_multipart_t   *mpart_a = NULL;
  if (filenames[0] != NULL) {
    nodal_a = 1;
    _read_volume_mesh(comm,
                      filenames[0],
                      rotate_a,
                      part_method,
                      n_part,
                      &dmn_a,
                      &mpart_a);
  }
  else {
    _generate_volume_mesh(comm,
                          n_vtx_a,
                          elt_type_a,
                          rotate_a,
                          noise_a,
                          0.,
                          0.,
                          0.,
                          length_a,
                          part_method,
                          n_part,
                          &dmn_a,
                          &mpart_a);
  }


  double length_b = 1.;
  int rotate_b = 0;
  PDM_dmesh_nodal_t *dmn_b   = NULL;
  PDM_multipart_t   *mpart_b = NULL;
  if (filenames[1] != NULL) {
    nodal_b = 1;
    _read_volume_mesh(comm,
                      filenames[1],
                      rotate_a,
                      part_method,
                      n_part,
                      &dmn_a,
                      &mpart_a);
  }
  else {
    _generate_volume_mesh(comm,
                          n_vtx_b,
                          elt_type_b,
                          rotate_b,
                          noise_b,
                          shift_b[0],
                          shift_b[1],
                          shift_b[2],
                          length_b,
                          part_method,
                          n_part,
                          &dmn_b,
                          &mpart_b);
  }

  if(0 == 1) {
    if (dmn_a != NULL) {
      PDM_dmesh_nodal_dump_vtk(dmn_a,
                               PDM_GEOMETRY_KIND_VOLUMIC,
                               "dmn_a_");
    }
    if (dmn_b != NULL) {
      PDM_dmesh_nodal_dump_vtk(dmn_b,
                               PDM_GEOMETRY_KIND_VOLUMIC,
                               "dmn_b_");
    }
  }

  /*
   * Mesh_intersection
   */
  int dim_mesh_a = 3;
  int dim_mesh_b = 3;
  PDM_mesh_intersection_t *mi = PDM_mesh_intersection_create(PDM_MESH_INTERSECTION_KIND_WEIGHT,
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
    pmn_a = _set_mesh_nodal(mi, 0, mpart_a);
  }
  else {
    _set_mesh(mi, 0, mpart_a, n_part);
  }
  if (nodal_b) {
    pmn_b = _set_mesh_nodal(mi, 1, mpart_b);
  }
  else {
    _set_mesh(mi, 1, mpart_b, n_part);
  }

  // tetraisation point
  PDM_mesh_intersection_tetraisation_pt_set(mi,
                                            tetraisation_pt_type,
                                            tetraisation_pt_coord);

  // compute
  PDM_mesh_intersection_compute(mi);



  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_intersection_part_to_part_get(mi,
                                         &ptp,
                                         PDM_OWNERSHIP_USER);

  double l_volume_AB = 0.;
  // Check ptp
  if (ptp != NULL) {
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


    for (int ipart = 0; ipart < n_part; ipart++) {
      double *elt_b_elt_a_volume = NULL;
      PDM_mesh_intersection_result_from_b_get(mi,
                                              ipart,
                                              &elt_b_elt_a_volume);
      for (int i = 0; i < pelt_b_elt_a_idx[ipart][n_ref_b[ipart]]; i++) {
        l_volume_AB += elt_b_elt_a_volume[i];
      }
    }

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
      int n_elt_a = PDM_multipart_part_ln_to_gn_get(mpart_a,
                                                    0,
                                                    ipart,
                                                    PDM_MESH_ENTITY_CELL,
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
      PDM_multipart_part_ln_to_gn_get(mpart_b,
                                      0,
                                      ipart,
                                      PDM_MESH_ENTITY_CELL,
                                      &elt_b_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);


      for (int i = 0; i < n_ref_b[ipart]; i++) {
        int faceB_id = ref_b[ipart][i] - 1;
        if (verbose) {
          log_trace("elt_b "PDM_FMT_G_NUM" : ", elt_b_ln_to_gn[faceB_id]);
          for (int j = pelt_b_elt_a_idx[ipart][i]; j < pelt_b_elt_a_idx[ipart][i+1]; j++) {
            log_trace("("PDM_FMT_G_NUM", %f)  ", pelt_b_elt_a[ipart][j], pelt_b_elt_a_volume[ipart][j]);
          }
          log_trace("\n");
        }
      }
    }

    PDM_part_to_part_free(ptp);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free(pelt_a_elt_b_n[ipart]);
    }
    free(pelt_a_elt_b_n     );
    free(pelt_a_elt_b_volume);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free(pelt_b_elt_a_volume[ipart]);
    }
    free(pelt_b_elt_a_volume);
  }


  double g_volume_AB;
  PDM_MPI_Allreduce(&l_volume_AB, &g_volume_AB, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  if (i_rank == 0) {
    printf("total volume of A inter B : %20.16f\n", g_volume_AB);
    double exact = 1;
    for (int i = 0; i < 3; i++) {
      if (shift_b[i] < 0) {
        exact *= PDM_MAX(0, PDM_MIN(length_a, length_b + shift_b[i]));
      }
      else {
        exact *= PDM_MAX(0, length_a - shift_b[i]);
      }
    }
    printf("error : absolute = %e, relative = %e\n",
           PDM_ABS(g_volume_AB - exact),
           PDM_ABS(g_volume_AB - exact)/exact);
  }



  /*
   * Free
   */

  PDM_mesh_intersection_free(mi);

  if (nodal_a) {
    PDM_part_mesh_nodal_free(pmn_a);
  }
  if (nodal_b) {
    PDM_part_mesh_nodal_free(pmn_b);
  }

  PDM_DMesh_nodal_free(dmn_b);
  PDM_multipart_free(mpart_b);

  PDM_DMesh_nodal_free(dmn_a);
  PDM_multipart_free(mpart_a);

  free(tetraisation_pt_coord);

  PDM_MPI_Barrier(comm);

  // if (i_rank == 0) {
  //   PDM_printf ("-- End\n");
  //   fflush(stdout);
  // }
  PDM_MPI_Finalize ();

  return 0;

}
