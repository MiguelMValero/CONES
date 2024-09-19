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
#include "pdm_mesh_location.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_part_to_block.h"

#include "pdm_vtk.h"

#include "pdm_generate_mesh.h"
#include "pdm_part_mesh_nodal.h"


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
     "  -t      <level>  Bounding boxes tolerance.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -octree          Use octree-based method.\n\n"
     "  -dbbree          Use dbbtree-based method.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scotch       Call PT-Scotch.\n\n"
     "  -hilbert         Call Hilbert.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                          argc,
           char                       **argv,
           PDM_g_num_t                 *n_vtx_seg,
           double                      *length,
           double                      *separation_x,
           double                      *separation_y,
           double                      *separation_z,
           int                         *n_part,
           int                         *n_block,
           int                         *post,
           int                         *part_method,
           PDM_mesh_location_method_t  *loc_method,
           PDM_Mesh_nodal_elt_t        *elt_type)
{
  int i = 1;

  PDM_UNUSED (post);

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
    else if (strcmp(argv[i], "-sep") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_y = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_z = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_block") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_block = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_PART_SPLIT_HILBERT;
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = PDM_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DBBTREE;
    }
    else if (strcmp(argv[i], "-doctree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DOCTREE;
    }
    else if (strcmp(argv[i], "-locate_all_tgt") == 0) {
      *loc_method = PDM_MESH_LOCATION_LOCATE_ALL_TGT;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
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
 const int                   active_rank,
 const PDM_MPI_Comm          comm,
 const int                   n_part,
 const int                   n_block,
 const PDM_split_dual_t      part_method,
 const PDM_g_num_t           n,
 const double                xmin,
 const double                ymin,
 const double                zmin,
 const double                length,
       int                ***n_cell,
       int               ****cell_vtx,
       PDM_g_num_t       ****cell_g_num,
       int               ****cell_parent_num,
       int                 **n_vtx,
       double             ***vtx_coord,
       PDM_g_num_t        ***vtx_g_num
 )
{
  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_HEXA8;
  int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, 1);

  *n_vtx     = malloc(sizeof(int          ) * n_part);
  *vtx_coord = malloc(sizeof(double      *) * n_part);
  *vtx_g_num = malloc(sizeof(PDM_g_num_t *) * n_part);


  *n_cell          = malloc(sizeof(int          *) * n_block);
  *cell_vtx        = malloc(sizeof(int         **) * n_block);
  *cell_g_num      = malloc(sizeof(PDM_g_num_t **) * n_block);
  *cell_parent_num = malloc(sizeof(int         **) * n_block);

  if (active_rank) {
    PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_parallelepiped(comm,
                                                                  elt_type,
                                                                  1,
                                                                  NULL,
                                                                  xmin,
                                                                  ymin,
                                                                  zmin,
                                                                  length,
                                                                  length,
                                                                  length,
                                                                  n,
                                                                  n,
                                                                  n,
                                                                  n_part,
                                                                  part_method);


    for (int iblock = 0; iblock < n_block; iblock++) {
      (*n_cell)         [iblock] = PDM_array_zeros_int(n_part);
      (*cell_vtx)       [iblock] = malloc(sizeof(int         *) * n_part);
      (*cell_g_num)     [iblock] = malloc(sizeof(PDM_g_num_t *) * n_part);
      (*cell_parent_num)[iblock] = malloc(sizeof(int         *) * n_part);
    }

    for (int ipart = 0; ipart < n_part; ipart++) {

      int part_n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn, 0, ipart);

      int         *connec              = NULL;
      PDM_g_num_t *numabs              = NULL;
      int         *parent_num          = NULL;
      PDM_g_num_t *parent_entity_g_num = NULL;
      PDM_part_mesh_nodal_section_std_get(pmn,
                                          0,
                                          ipart,
                                          &connec,
                                          &numabs,
                                          &parent_num,
                                          &parent_entity_g_num,
                                          PDM_OWNERSHIP_KEEP);

      int block_n_elt = part_n_elt / n_block;

      for (int iblock = 0; iblock < n_block-1; iblock++) {
        (*n_cell)[iblock][ipart] = block_n_elt;
      }

      (*n_cell)[n_block-1][ipart] = part_n_elt - (n_block-1)*block_n_elt;


      int idx = 0;
      for (int iblock = 0; iblock < n_block; iblock++) {
        (*cell_vtx)[iblock][ipart] = malloc(sizeof(int) * elt_vtx_n * (*n_cell)[iblock][ipart]);

        memcpy((*cell_vtx)[iblock][ipart],
               connec + elt_vtx_n * idx,
               sizeof(int) * elt_vtx_n * (*n_cell)[iblock][ipart]);

        (*cell_g_num)[iblock][ipart] = malloc(sizeof(PDM_g_num_t) * (*n_cell)[iblock][ipart]);
        memcpy((*cell_g_num)[iblock][ipart],
               numabs + idx,
               sizeof(PDM_g_num_t) * (*n_cell)[iblock][ipart]);

        (*cell_parent_num)[iblock][ipart] = malloc(sizeof(int) * (*n_cell)[iblock][ipart]);
        for (int i = 0; i < (*n_cell)[iblock][ipart]; i++) {
          (*cell_parent_num)[iblock][ipart][i] = idx + i;
        }

        idx += (*n_cell)[iblock][ipart];
      }

      (*n_vtx)[ipart]         = PDM_part_mesh_nodal_n_vtx_get(pmn, ipart);
      double      *_vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, ipart);
      PDM_g_num_t *_vtx_g_num = PDM_part_mesh_nodal_vtx_g_num_get(pmn, ipart);

      (*vtx_coord)[ipart] = malloc(sizeof(double) * (*n_vtx)[ipart] * 3);
      memcpy((*vtx_coord)[ipart],
             _vtx_coord,
             sizeof(double) * (*n_vtx)[ipart] * 3);

      (*vtx_g_num)[ipart] = malloc(sizeof(PDM_g_num_t) * (*n_vtx)[ipart]);
      memcpy((*vtx_g_num)[ipart],
             _vtx_g_num,
             sizeof(PDM_g_num_t) * (*n_vtx)[ipart]);

    }

    PDM_part_mesh_nodal_free(pmn);
  }
  else {

    for (int iblock = 0; iblock < n_block; iblock++) {
      (*n_cell)         [iblock] = PDM_array_zeros_int(n_part);
      (*cell_vtx)       [iblock] = malloc(sizeof(int         *) * n_part);
      (*cell_g_num)     [iblock] = malloc(sizeof(PDM_g_num_t *) * n_part);
      (*cell_parent_num)[iblock] = malloc(sizeof(int         *) * n_part);
    }

    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int iblock = 0; iblock < n_block; iblock++) {
        (*cell_vtx)       [iblock][ipart] = malloc(sizeof(int)         * elt_vtx_n * (*n_cell)[iblock][ipart]);
        (*cell_g_num)     [iblock][ipart] = malloc(sizeof(PDM_g_num_t) * elt_vtx_n);
        (*cell_parent_num)[iblock][ipart] = malloc(sizeof(int)         * elt_vtx_n);
      }

      (*n_vtx)    [ipart] = 0;
      (*vtx_coord)[ipart] = malloc(sizeof(double     ) * (*n_vtx)[ipart] * 3);
      (*vtx_g_num)[ipart] = malloc(sizeof(PDM_g_num_t) * (*n_vtx)[ipart]);
    }

  }
}




static PDM_part_mesh_nodal_t *
_multisection_pmn
(
 const PDM_MPI_Comm   comm,
 const int            n_part,
 const int            n_section,
       int          **n_cell,
       int         ***cell_vtx,
       PDM_g_num_t ***cell_g_num,
       int         ***cell_parent_num,
       int           *n_vtx,
       double       **vtx_coord,
       PDM_g_num_t  **vtx_g_num
 )
{
  PDM_part_mesh_nodal_t *pmn = PDM_part_mesh_nodal_create(3,
                                                          n_part,
                                                          comm);


  for (int i_section = 0; i_section < n_section; i_section++) {
    PDM_part_mesh_nodal_section_add(pmn,
                                    PDM_MESH_NODAL_HEXA8);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_nodal_coord_set(pmn,
                                  i_part,
                                  n_vtx    [i_part],
                                  vtx_coord[i_part],
                                  vtx_g_num[i_part],
                                  PDM_OWNERSHIP_KEEP);

    for (int i_section = 0; i_section < n_section; i_section++) {
      PDM_part_mesh_nodal_section_std_set(pmn,
                                          i_section,
                                          i_part,
                                          n_cell         [i_section][i_part],
                                          cell_vtx       [i_section][i_part],
                                          cell_g_num     [i_section][i_part],
                                          cell_parent_num[i_section][i_part],
                                          NULL, //?cell_g_num[i_section][i_part],
                                          PDM_OWNERSHIP_KEEP);
    }
  }

  return pmn;
}


static PDM_part_mesh_nodal_t *
_gen_pmn
(
 const int                   current_rank_has_mesh,
 const PDM_MPI_Comm          comm,
 const int                   n_part,
 const int                   n_section,
 const PDM_split_dual_t      part_method,
 const PDM_g_num_t           n_vtx_seg,
 const double                xmin,
 const double                ymin,
 const double                zmin,
 const double                length
)
{
  int          **n_cell          = NULL;
  int         ***cell_vtx        = NULL;
  PDM_g_num_t ***cell_g_num      = NULL;
  int         ***cell_parent_num = NULL;
  int           *n_vtx           = NULL;
  double       **vtx_coord       = NULL;
  PDM_g_num_t  **vtx_g_num       = NULL;
  _gen_mesh(current_rank_has_mesh,
            comm,
            n_part,
            n_section,
            part_method,
            n_vtx_seg,
            xmin,
            ymin,
            zmin,
            length,
            &n_cell,
            &cell_vtx,
            &cell_g_num,
            &cell_parent_num,
            &n_vtx,
            &vtx_coord,
            &vtx_g_num);

  PDM_part_mesh_nodal_t *pmn = _multisection_pmn(comm,
                                                 n_part,
                                                 n_section,
                                                 n_cell,
                                                 cell_vtx,
                                                 cell_g_num,
                                                 cell_parent_num,
                                                 n_vtx,
                                                 vtx_coord,
                                                 vtx_g_num);

  for (int i_section = 0; i_section < n_section; i_section++) {
    free(n_cell         [i_section]);
    free(cell_vtx       [i_section]);
    free(cell_g_num     [i_section]);
    free(cell_parent_num[i_section]);
  }

  free(n_cell         );
  free(cell_vtx       );
  free(cell_g_num     );
  free(cell_parent_num);
  free(n_vtx          );
  free(vtx_coord      );
  free(vtx_g_num      );

  return pmn;
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

  PDM_g_num_t                n_vtx_seg    = 10;
  double                     length       = 1.;
  double                     separation_x = 2.;
  double                     separation_y = 0.;
  double                     separation_z = 0.;
  double                     tolerance    = 1e-6;
  int                        n_part       = 1;
  int                        n_block      = 1;
  int                        post         = 0;
  PDM_split_dual_t           part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_mesh_location_method_t loc_method   = PDM_MESH_LOCATION_OCTREE;
  PDM_Mesh_nodal_elt_t       elt_type     = PDM_MESH_NODAL_HEXA8;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &separation_x,
             &separation_y,
             &separation_z,
             &n_part,
             &n_block,
             &post,
     (int *) &part_method,
             &loc_method,
             &elt_type);


  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 3);


  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  const double xmin = -0.5*length;
  const double ymin = -0.5*length;
  const double zmin = -0.5*length;

  int current_rank_has_mesh = 1;


  /*
   *  Source mesh
   */
  PDM_part_mesh_nodal_t *src_pmn = _gen_pmn(current_rank_has_mesh,
                                            comm,
                                            n_part,
                                            n_block,
                                            part_method,
                                            n_vtx_seg,
                                            xmin,
                                            ymin,
                                            zmin,
                                            length);

  /*
   *  Target mesh
   */
  PDM_part_mesh_nodal_t *tgt_pmn = _gen_pmn(current_rank_has_mesh,
                                            comm,
                                            n_part,
                                            1,//n_block,
                                            part_method,
                                            n_vtx_seg,
                                            xmin + separation_x,
                                            ymin + separation_y,
                                            zmin + separation_z,
                                            length);


  /*
   *  Mesh location structure initialization
   */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create(1,
                                                           comm,
                                                           PDM_OWNERSHIP_KEEP);

  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set(mesh_loc,
                                     0,
                                     n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int          n_tgt     = PDM_part_mesh_nodal_n_vtx_get    (tgt_pmn, i_part);
    double      *tgt_coord = PDM_part_mesh_nodal_vtx_coord_get(tgt_pmn, i_part);
    PDM_g_num_t *tgt_g_num = PDM_part_mesh_nodal_vtx_g_num_get(tgt_pmn, i_part);

    PDM_mesh_location_cloud_set(mesh_loc,
                                0,
                                i_part,
                                n_tgt,
                                tgt_coord,
                                tgt_g_num);
  }


  /* Set source mesh */
  PDM_mesh_location_shared_nodal_mesh_set(mesh_loc, src_pmn);


  /* Set location parameters */
  PDM_mesh_location_tolerance_set(mesh_loc,
                                  tolerance);

  PDM_mesh_location_method_set(mesh_loc,
                               loc_method);

  /*
   *  Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  PDM_mesh_location_compute(mesh_loc);

  PDM_mesh_location_dump_times(mesh_loc);



  /* Check */
  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_KEEP);

  double **send_field = malloc(sizeof(double *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {

    double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(src_pmn, i_part);

    int n_cell = PDM_part_mesh_nodal_n_elmts_get(src_pmn,
                                                 PDM_GEOMETRY_KIND_VOLUMIC,
                                                 i_part);

    int         *elt_pts_idx        = NULL;
    PDM_g_num_t *elt_pts_gnum       = NULL;
    double      *elt_pts_coord      = NULL;
    double      *elt_pts_uvw        = NULL;
    int         *elt_pts_weight_idx = NULL;
    double      *elt_pts_weight     = NULL;
    double      *elt_pts_dist2      = NULL;
    double      *elt_pts_proj_coord = NULL;
    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        i_part,
                                        0, // i_point_cloud,
                                        &elt_pts_idx,
                                        &elt_pts_gnum,
                                        &elt_pts_coord,
                                        &elt_pts_uvw,
                                        &elt_pts_weight_idx,
                                        &elt_pts_weight,
                                        &elt_pts_dist2,
                                        &elt_pts_proj_coord);

    int *cell_vtx_idx = NULL;
    int *cell_vtx     = NULL;
    PDM_mesh_location_cell_vertex_get(mesh_loc,
                                      i_part,
                                      &cell_vtx_idx,
                                      &cell_vtx);

    if (0 && post) {
      int          n_vtx     = PDM_part_mesh_nodal_n_vtx_get    (src_pmn, i_part);
      PDM_g_num_t *vtx_g_num = PDM_part_mesh_nodal_vtx_g_num_get(src_pmn, i_part);

      char filename[999];
      sprintf(filename, "mesh_location_multisection_src_%d_%d.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements(filename,
                                 n_vtx,
                                 vtx_coord,
                                 vtx_g_num,
                                 PDM_MESH_NODAL_HEXA8,
                                 n_cell,
                                 cell_vtx,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL);
    }


    send_field[i_part] = malloc(sizeof(double) * elt_pts_idx[n_cell]);
    for (int icell = 0; icell < n_cell; icell++) {

      int *cv = cell_vtx + cell_vtx_idx[icell];

      for (int idx_pt = elt_pts_idx[icell]; idx_pt < elt_pts_idx[icell+1]; idx_pt++) {
        send_field[i_part][idx_pt] = 0.;
        int idx_vtx = 0;

        for (int idx_w = elt_pts_weight_idx[idx_pt]; idx_w < elt_pts_weight_idx[idx_pt+1]; idx_w++) {
          int vtx_id = cv[idx_vtx++] - 1;
          for (int i = 0; i < 3; i++) {
            send_field[i_part][idx_pt] += elt_pts_weight[idx_w] * vtx_coord[3*vtx_id+i];
          }
        }
      }

    }


  }


  double **recv_field = NULL;
  int request = -1;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) send_field,
                         NULL,
        (      void ***) &recv_field,
                         &request);


  PDM_part_to_part_iexch_wait(ptp, request);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(send_field[i_part]);
  }
  free(send_field);


  double lmax_err = 0.;
  for (int i_part = 0; i_part < n_part; i_part++) {
    double      *tgt_coord = PDM_part_mesh_nodal_vtx_coord_get(tgt_pmn, i_part);
    PDM_g_num_t *tgt_g_num = PDM_part_mesh_nodal_vtx_g_num_get(tgt_pmn, i_part);


    int n_located = PDM_mesh_location_n_located_get(mesh_loc,
                                                    0,//i_point_cloud,
                                                    i_part);

    int *located = PDM_mesh_location_located_get(mesh_loc,
                                                 0,//i_point_cloud,
                                                 i_part);

    PDM_g_num_t *p_location    = NULL;
    double      *p_dist2  = NULL;
    double      *p_proj_coord  = NULL;
    PDM_mesh_location_point_location_get(mesh_loc,
                                         0,//i_point_cloud,
                                         i_part,
                                         &p_location,
                                         &p_dist2,
                                         &p_proj_coord);

    for (int i = 0; i < n_located; i++) {
      int pt_id = located[i] - 1;

      double f = 0.;
      for (int j = 0; j < 3; j++) {
        f += tgt_coord[3*pt_id + j];
      }

      double err = PDM_ABS(recv_field[i_part][i] - f);
      lmax_err = PDM_MAX(lmax_err, err);

      if (err > 1.e-12) {
        log_trace("point "PDM_FMT_G_NUM" (%f %f %f) located in elt "PDM_FMT_G_NUM" at dist %e: error = %e (%20.16f / %20.16f)\n",
                  tgt_g_num[pt_id],
                  tgt_coord[3*pt_id], tgt_coord[3*pt_id+1], tgt_coord[3*pt_id+2],
                  p_location[i],
                  PDM_SIGN(p_dist2[i])*sqrt(PDM_ABS(p_dist2[i])),
                  err, recv_field[i_part][i], f);
      }
    }

    if (post) {
      int n_tgt = PDM_part_mesh_nodal_n_vtx_get(tgt_pmn, i_part);
      char filename[999];

      double *exact    = malloc(sizeof(double) * n_tgt);
      double *interp   = malloc(sizeof(double) * n_tgt);
      double *location = malloc(sizeof(double) * n_tgt);
      double *proj     = malloc(sizeof(double) * n_tgt * 3);
      for (int i = 0; i < n_tgt; i++) {
        exact   [i] = tgt_coord[3*i] + tgt_coord[3*i+1] + tgt_coord[3*i+2];
        interp  [i] = -1;
        location[i] = -1;
        proj[3*i] = proj[3*i+1] = proj[3*i+2] = -1;
      }
      for (int i = 0; i < n_located; i++) {
        int pt_id = located[i] - 1;
        interp  [pt_id] = recv_field[i_part][i];
        location[pt_id] = (double) p_location[i];
        memcpy(proj + 3*pt_id, p_proj_coord + 3*i, sizeof(double) * 3);
      }

      const char   *field_name [3] = {"exact", "interp", "location"};
      const double *field_value[3] = {exact, interp, location};

      sprintf(filename, "mesh_location_multisection_tgt_%d_%d.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements_double(filename,
                                        n_tgt,
                                        tgt_coord,
                                        tgt_g_num,
                                        PDM_MESH_NODAL_POINT,
                                        n_tgt,
                                        NULL,
                                        tgt_g_num,
                                        3,
                                        field_name,
                                        field_value);

      sprintf(filename, "mesh_location_multisection_tgt_proj_%d_%d.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements_double(filename,
                                        n_tgt,
                                        proj,
                                        tgt_g_num,
                                        PDM_MESH_NODAL_POINT,
                                        n_tgt,
                                        NULL,
                                        tgt_g_num,
                                        3,
                                        field_name,
                                        field_value);
      free(exact);
      free(interp);
      free(location);
      free(proj);
    }

    free(recv_field[i_part]);
  }
  free(recv_field);

  double gmax_err;
  PDM_MPI_Allreduce(&lmax_err, &gmax_err, 1, PDM_MPI_DOUBLE,
                    PDM_MPI_MAX, PDM_MPI_COMM_WORLD);


  if (i_rank == 0) {
    printf("global max interpolation error = %e\n", gmax_err);
  }
  if (gmax_err > tolerance) { // scaling??
    PDM_error(__FILE__, __LINE__, 0, "Large interpolation error!\n");
  }


  /* Free memory */
  PDM_mesh_location_free(mesh_loc);
  PDM_part_mesh_nodal_free(src_pmn);
  PDM_part_mesh_nodal_free(tgt_pmn);


  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return EXIT_SUCCESS;
}




