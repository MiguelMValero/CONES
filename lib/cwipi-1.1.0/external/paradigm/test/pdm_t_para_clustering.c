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
#include "pdm_printf.h"
#include "pdm_gnum.h"
#include "pdm_dcube_gen.h"
#include "pdm_part.h"
#include "pdm_sort.h"
#include "pdm_part_to_part.h"
#include "pdm_array.h"
#include "pdm_closest_points.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_para_octree.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"

#include "pdm_mpi.h"
#include "pdm_config.h"


/**
 *
 * \brief  Compute gnum
 *
 */

static void
_compute_gnum
(
  const int           n_part,
  const int          *n_entity,
        PDM_g_num_t **entity_ln_to_gn,
        PDM_MPI_Comm  comm
) {

  PDM_gen_gnum_t *gnum = PDM_gnum_create(3,
                                         n_part,
                                         PDM_FALSE,
                                         1.,
                                         comm,
                                         PDM_OWNERSHIP_USER);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_parents(gnum,
                              i_part,
                              n_entity[i_part],
                              entity_ln_to_gn[i_part]);
  }

  PDM_gnum_compute(gnum);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(entity_ln_to_gn[i_part]);
    entity_ln_to_gn[i_part] = PDM_gnum_get(gnum, i_part);
  }

  PDM_gnum_free(gnum);

}

static
void
_mean_data_per_leaf
(
  const int      n_explicit_nodes,
  const int     *n_points,
  const int     *range,
  const int     *leaf_id,
  const int      stride,
  const double  *pts_data_octree,
        int     *n_leaf,
        double **leaf_data_octree
)
{

  int _n_leaf = 0;

  int *extract_leaf_node_id = malloc(n_explicit_nodes * sizeof(int));
  for(int i_node = 0; i_node < n_explicit_nodes; i_node++) {
    if(leaf_id[i_node] != -1){
      extract_leaf_node_id[_n_leaf++] = i_node;
    }
  }
  extract_leaf_node_id = realloc(extract_leaf_node_id, _n_leaf * sizeof(int));

  double *_leaf_data_octree = malloc(stride * _n_leaf * sizeof(double));

  for (int i_leaf = 0; i_leaf < _n_leaf; i_leaf++) {
    _leaf_data_octree[3*i_leaf    ] = 0.0;
    _leaf_data_octree[3*i_leaf + 1] = 0.0;
    _leaf_data_octree[3*i_leaf + 2] = 0.0;
  }

  for (int i_leaf = 0; i_leaf < _n_leaf; i_leaf++) {
    for (int i_pts = range[extract_leaf_node_id[i_leaf]]; i_pts < range[extract_leaf_node_id[i_leaf]] + n_points[extract_leaf_node_id[i_leaf]]; i_pts++) {
      for (int i = 0; i < stride; i++) {
        _leaf_data_octree[stride*i_leaf+i] = _leaf_data_octree[stride*i_leaf+i] + pts_data_octree[stride*i_pts+i];
      }
    }
    for (int i = 0; i < stride; i++) {
      _leaf_data_octree[stride*i_leaf+i] /= n_points[extract_leaf_node_id[i_leaf]];
    }
  }

  free(extract_leaf_node_id);

  *n_leaf            = _n_leaf;
  *leaf_data_octree = _leaf_data_octree;

}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init(&argc, &argv);

  /*
   * Parameters
   */

  PDM_bool_t  show_mesh        = PDM_FALSE;
  PDM_g_num_t n_vtx_seg        = 30;
  double      length           = 1e-2;
  double      zero_x           = -length/2;
  double      zero_y           = -length/2;
  double      zero_z           = -length/2;
  int         n_part           = 2;
  double      min_dist         = pow(length/20, 2);
  int         n_layer          = 4;
  int        *n_leaf_per_layer = malloc (n_layer*sizeof(int));
  int         n_vtx_min_dist   = 10;
  int         depth_max        = 50;
  int         n_var            = 4;

  n_leaf_per_layer[0] = 31;
  n_leaf_per_layer[1] = 62;
  n_leaf_per_layer[2] = 125;
  n_leaf_per_layer[3] = 250;

  /*
   * Generate a cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int          dn_face_group;
  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int          dsface_vtx;
  int          dsface_group;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  if (i_rank == 0) {
    printf("-- Generate cube\n");
    fflush(stdout);
  }

  PDM_dcube_t *dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          zero_x,
                                          zero_y,
                                          zero_z,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                       &dn_face_group,
                       &dn_cell,
                       &dn_face,
                       &dn_vtx,
                       &dsface_vtx,
                       &dsface_group);

  PDM_dcube_gen_data_get(dcube,
                        &dface_cell,
                        &dface_vtx_idx,
                        &dface_vtx,
                        &dvtx_coord,
                        &dface_group_idx,
                        &dface_group);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  if (i_rank == 0) {
    printf("-- Part cube\n");
    fflush(stdout);
  }

  PDM_part_t* ppart = PDM_part_create(comm,
                                      PDM_PART_SPLIT_HILBERT,
                                      "PDM_PART_RENUM_CELL_NONE",
                                      "PDM_PART_RENUM_FACE_NONE",
                                      n_property_cell,
                                      renum_properties_cell,
                                      n_property_face,
                                      renum_properties_face,
                                      n_part,
                                      dn_cell,
                                      dn_face,
                                      dn_vtx,
                                      dn_face_group,
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL,
                                      have_dcell_part,
                                      dcell_part,
                                      dface_cell,
                                      dface_vtx_idx,
                                      dface_vtx,
                                      NULL,
                                      dvtx_coord,
                                      NULL,
                                      dface_group_idx,
                                      dface_group);

  free(dcell_part);

  PDM_g_num_t **gnum_bnd_vtx = malloc (sizeof(PDM_g_num_t *) * (n_part+1));
  int         **lnum_bnd_vtx = malloc (sizeof(int         *) *  n_part   );
  double      **bnd_vtx      = malloc (sizeof(double      *) * (n_part+1));
  double      **data_bnd_vtx = malloc (sizeof(double      *) * (n_part+1));
  int          *n_bnd_vtx    = malloc (sizeof(int          ) * (n_part+1));

  PDM_g_num_t **gnum_all_vtx = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **all_vtx      = malloc (sizeof(double      *) * n_part);
  int          *n_all_vtx    = malloc (sizeof(int          ) * n_part);

  PDM_g_num_t **bnd_to_all_vtx     = malloc (sizeof(PDM_g_num_t *) * n_part);
  int         **bnd_to_all_vtx_idx = malloc (sizeof(int         *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nface_group;

    PDM_part_part_dim_get(ppart,
                          i_part,
                         &n_cell,
                         &n_face,
                         &n_face_part_bound,
                         &n_vtx,
                         &n_proc,
                         &n_total_part,
                         &scell_face,
                         &sface_vtx,
                         &sface_group,
                         &nface_group);

    int         *cell_tag                 = NULL;
    int         *cell_face_idx            = NULL;
    int         *cell_face                = NULL;
    PDM_g_num_t *cell_ln_to_gn            = NULL;
    int         *face_tag                 = NULL;
    int         *face_cell                = NULL;
    int         *face_vtx_idx             = NULL;
    int         *face_vtx                 = NULL;
    PDM_g_num_t *face_ln_to_gn            = NULL;
    int         *face_part_bound_proc_idx = NULL;
    int         *face_part_bound_part_idx = NULL;
    int         *face_part_bound          = NULL;
    int         *vtx_tag                  = NULL;
    double      *vtx                      = NULL;
    PDM_g_num_t *vtx_ln_to_gn             = NULL;
    int         *face_group_idx           = NULL;
    int         *face_group               = NULL;
    PDM_g_num_t *face_group_ln_to_gn      = NULL;

    PDM_part_part_val_get(ppart,
                          i_part,
                         &cell_tag,
                         &cell_face_idx,
                         &cell_face,
                         &cell_ln_to_gn,
                         &face_tag,
                         &face_cell,
                         &face_vtx_idx,
                         &face_vtx,
                         &face_ln_to_gn,
                         &face_part_bound_proc_idx,
                         &face_part_bound_part_idx,
                         &face_part_bound,
                         &vtx_tag,
                         &vtx,
                         &vtx_ln_to_gn,
                         &face_group_idx,
                         &face_group,
                         &face_group_ln_to_gn);

    if (show_mesh) {
      char filename[999];
      sprintf(filename, "ini_face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                             n_vtx,
                             vtx,
                             vtx_ln_to_gn,
                             n_face,
                             face_vtx_idx,
                             face_vtx,
                             face_ln_to_gn,
                             NULL);
    }

    n_all_vtx[i_part] = n_vtx;
    gnum_all_vtx[i_part] = malloc (sizeof(PDM_g_num_t)     * n_all_vtx[i_part]);
    all_vtx     [i_part] = malloc (sizeof(double     ) * 3 * n_all_vtx[i_part]);

    for (int i_vtx = 0; i_vtx < n_all_vtx[i_part]; i_vtx++) {
      gnum_all_vtx[i_part][  i_vtx] = vtx_ln_to_gn[i_vtx];
      all_vtx     [i_part][3*i_vtx    ] = vtx[3*i_vtx    ];
      all_vtx     [i_part][3*i_vtx + 1] = vtx[3*i_vtx + 1];
      all_vtx     [i_part][3*i_vtx + 2] = vtx[3*i_vtx + 2];
    }

   /*
    * Get boundary vertices
    */

    n_bnd_vtx[i_part] = 0;

    for (int i_group = 0; i_group < nface_group; i_group++) {
      for (int i_face_group = face_group_idx[i_group]; i_face_group < face_group_idx[i_group+1]; i_face_group++) {
        int i_face = face_group[i_face_group]-1;
        n_bnd_vtx[i_part] += face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
      }
    }

    lnum_bnd_vtx[i_part] = malloc (sizeof(int) * n_bnd_vtx[i_part]);

    n_bnd_vtx[i_part] = 0;

    for (int i_group = 0; i_group < nface_group; i_group++) {
      for (int i_face_group = face_group_idx[i_group]; i_face_group < face_group_idx[i_group+1]; i_face_group++) {
        int i_face = face_group[i_face_group]-1;
        for (int i_vtx = face_vtx_idx[i_face]; i_vtx < face_vtx_idx[i_face+1]; i_vtx++) {
          lnum_bnd_vtx[i_part][n_bnd_vtx[i_part]++] = face_vtx[i_vtx]-1;
        }
      }
    }

    PDM_sort_int(lnum_bnd_vtx[i_part], NULL, n_bnd_vtx[i_part]);

    int n_vtx_tmp = n_bnd_vtx[i_part];
    n_bnd_vtx[i_part] = 1;

    for (int i_vtx = 1; i_vtx < n_vtx_tmp; i_vtx++) {
      if (lnum_bnd_vtx[i_part][i_vtx] > lnum_bnd_vtx[i_part][i_vtx-1]){
        lnum_bnd_vtx[i_part][n_bnd_vtx[i_part]++] = lnum_bnd_vtx[i_part][i_vtx];
      }
    }

    lnum_bnd_vtx[i_part] = realloc (lnum_bnd_vtx[i_part], sizeof(int) * n_bnd_vtx[i_part]);

    gnum_bnd_vtx[i_part] = malloc (sizeof(PDM_g_num_t)         * n_bnd_vtx[i_part]);
    bnd_vtx     [i_part] = malloc (sizeof(double     ) * 3     * n_bnd_vtx[i_part]);
    data_bnd_vtx[i_part] = malloc (sizeof(double     ) * n_var * n_bnd_vtx[i_part]);

    bnd_to_all_vtx    [i_part] = malloc (sizeof(PDM_g_num_t) * n_bnd_vtx[i_part]);
    bnd_to_all_vtx_idx[i_part] = PDM_array_new_idx_from_const_stride_int(1, n_bnd_vtx[i_part]);

    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      gnum_bnd_vtx  [i_part][  i_vtx    ] = vtx_ln_to_gn[lnum_bnd_vtx[i_part][i_vtx]];
      bnd_to_all_vtx[i_part][  i_vtx    ] = vtx_ln_to_gn[lnum_bnd_vtx[i_part][i_vtx]];
      bnd_vtx       [i_part][3*i_vtx    ] = vtx[3*lnum_bnd_vtx[i_part][i_vtx]    ];
      bnd_vtx       [i_part][3*i_vtx + 1] = vtx[3*lnum_bnd_vtx[i_part][i_vtx] + 1];
      bnd_vtx       [i_part][3*i_vtx + 2] = vtx[3*lnum_bnd_vtx[i_part][i_vtx] + 2];
      for (int i_var = 0; i_var < n_var; i_var++) {
        data_bnd_vtx[i_part][n_var*i_vtx+i_var] = 0.0;
      }
      data_bnd_vtx[i_part][n_var*i_vtx    ] = 1.0;
      data_bnd_vtx[i_part][n_var*i_vtx + 1] = 2e-3*(1.0-(bnd_vtx[i_part][3*i_vtx]+length/2)/length);
    }

  }

  if (i_rank == 0) {
    printf("-- Gnum boundary vertices\n");
    fflush(stdout);
  }

  _compute_gnum(n_part,
                n_bnd_vtx,
                gnum_bnd_vtx,
                comm);

  PDM_part_to_part_t *ptp = PDM_part_to_part_create((const PDM_g_num_t **) gnum_bnd_vtx,
                                                    (const int          *) n_bnd_vtx,
                                                                           n_part,
                                                    (const PDM_g_num_t **) gnum_all_vtx,
                                                    (const int          *) n_all_vtx,
                                                                           n_part,
                                                    (const int         **) bnd_to_all_vtx_idx,
                                                    (const PDM_g_num_t **) bnd_to_all_vtx,
                                                                           comm);

  int  *n_unref_vtx = NULL;
  int **unref_vtx   = NULL;

  PDM_part_to_part_unref_lnum2_get(ptp,
                                  &n_unref_vtx,
                                  &unref_vtx);

  PDM_g_num_t **gnum_int_vtx          = malloc (sizeof(PDM_g_num_t *) * n_part);
  int         **lnum_int_vtx          = malloc (sizeof(int         *) * n_part);
  double      **int_vtx               = malloc (sizeof(double      *) * n_part);
  int          *n_int_vtx             = malloc (sizeof(int          ) * n_part);
  PDM_g_num_t **gnum_bnd_vtx_min_dist = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **bnd_vtx_min_dist      = malloc (sizeof(double      *) * n_part);
  int          *n_bnd_vtx_min_dist    = malloc (sizeof(int          ) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

   /*
    * Get interior vertices
    */

    n_int_vtx[i_part] = n_unref_vtx[i_part];

    lnum_int_vtx[i_part] = malloc (sizeof(int) * n_int_vtx[i_part]);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      lnum_int_vtx[i_part][i_vtx] = unref_vtx[i_part][i_vtx]-1;
    }

    int         *cell_tag                 = NULL;
    int         *cell_face_idx            = NULL;
    int         *cell_face                = NULL;
    PDM_g_num_t *cell_ln_to_gn            = NULL;
    int         *face_tag                 = NULL;
    int         *face_cell                = NULL;
    int         *face_vtx_idx             = NULL;
    int         *face_vtx                 = NULL;
    PDM_g_num_t *face_ln_to_gn            = NULL;
    int         *face_part_bound_proc_idx = NULL;
    int         *face_part_bound_part_idx = NULL;
    int         *face_part_bound          = NULL;
    int         *vtx_tag                  = NULL;
    double      *vtx                      = NULL;
    PDM_g_num_t *vtx_ln_to_gn             = NULL;
    int         *face_group_idx           = NULL;
    int         *face_group               = NULL;
    PDM_g_num_t *face_group_ln_to_gn      = NULL;

    PDM_part_part_val_get(ppart,
                          i_part,
                         &cell_tag,
                         &cell_face_idx,
                         &cell_face,
                         &cell_ln_to_gn,
                         &face_tag,
                         &face_cell,
                         &face_vtx_idx,
                         &face_vtx,
                         &face_ln_to_gn,
                         &face_part_bound_proc_idx,
                         &face_part_bound_part_idx,
                         &face_part_bound,
                         &vtx_tag,
                         &vtx,
                         &vtx_ln_to_gn,
                         &face_group_idx,
                         &face_group,
                         &face_group_ln_to_gn);

    gnum_int_vtx[i_part] = malloc (sizeof(PDM_g_num_t)     * n_int_vtx[i_part]);
    int_vtx     [i_part] = malloc (sizeof(double     ) * 3 * n_int_vtx[i_part]);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      gnum_int_vtx[i_part][i_vtx] = vtx_ln_to_gn[lnum_int_vtx[i_part][i_vtx]];
      int_vtx     [i_part][3*i_vtx    ] = vtx[3*lnum_int_vtx[i_part][i_vtx]    ];
      int_vtx     [i_part][3*i_vtx + 1] = vtx[3*lnum_int_vtx[i_part][i_vtx] + 1];
      int_vtx     [i_part][3*i_vtx + 2] = vtx[3*lnum_int_vtx[i_part][i_vtx] + 2];
    }

  }

  if (i_rank == 0) {
    printf("-- Gnum interior vertices\n");
    fflush(stdout);
  }

  _compute_gnum(n_part,
                n_int_vtx,
                gnum_int_vtx,
                comm);

  PDM_part_to_part_free(ptp);
  PDM_dcube_gen_free(dcube);

  /*
   * Build layers of parallel octree
   */

  if (i_rank == 0) {
    printf("-- Build layers of parallel octree\n");
    fflush(stdout);
  }

  PDM_g_num_t _max_g_num = 0;
  PDM_g_num_t  max_g_num = 0;

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      _max_g_num = PDM_MAX(_max_g_num, gnum_bnd_vtx[i_part][i_vtx]);
    }
  }

  PDM_MPI_Allreduce(&_max_g_num, &max_g_num, 1, PDM_MPI_LONG, PDM_MPI_MAX, comm);

  _max_g_num = max_g_num + 1;

  PDM_g_num_t **gnum_layer_vtx   = malloc (sizeof(PDM_g_num_t *) *  n_layer   );
  PDM_g_num_t **distri_layer_vtx = malloc (sizeof(int         *) *  n_layer   );
  double      **layer_vtx        = malloc (sizeof(double      *) *  n_layer   );
  double      **data_layer_vtx   = malloc (sizeof(double      *) *  n_layer   );
  int          *n_layer_vtx      = malloc (sizeof(int          ) *  n_layer   );
  int          *m_layer_vtx      = malloc (sizeof(int          ) *  n_layer   );
  int          *layer_vtx_idx    = malloc (sizeof(int          ) * (n_layer+1));
  layer_vtx_idx[0] = 0;

  for (int i_layer = 0; i_layer < n_layer; i_layer++) {

    if (i_rank == 0) {
      printf("-- Build layers of parallel octree: layer %i/%i\n", i_layer+1, n_layer);
      fflush(stdout);
    }

    PDM_para_octree_t *octree = PDM_para_octree_create(1,
                                                       depth_max,
                                                       n_leaf_per_layer[i_layer],
                                                       0,
                                                       comm);

    for (int i_part = 0; i_part < n_part; i_part++) {

      PDM_para_octree_point_cloud_set(octree,
                                      i_part,
                                      n_bnd_vtx[i_part],
                                      bnd_vtx[i_part],
                                      gnum_bnd_vtx[i_part]);

    }

    PDM_para_octree_build(octree, NULL);

    int  n_explicit_nodes = 0;
    int *n_points         = NULL;
    int *range            = NULL;
    int *leaf_id          = NULL;
    int *children_id      = NULL;
    int *ancestor_id      = NULL;
    int  n_child          = 0;
    int  stack_size       = 0;
    PDM_para_octree_explicit_node_get(octree,
                                     &n_explicit_nodes,
                                     &n_points,
                                     &range,
                                     &leaf_id,
                                     &children_id,
                                     &ancestor_id,
                                     &n_child,
                                     &stack_size);

    int          n_pts_octree     = 0;
    double      *pts_coord_octree = NULL;
    PDM_g_num_t *pts_gnum_octree  = NULL;
    PDM_para_octree_points_get(octree,
                              &n_pts_octree,
                              &pts_coord_octree,
                              &pts_gnum_octree);

    /*
     * Extend boundary vertices gnum by leaves barycenters and data at leaves barycenters
     */

    int     n_leaf            = 0;
    double *leaf_coord_octree = NULL;
    _mean_data_per_leaf(n_explicit_nodes,
                        n_points,
                        range,
                        leaf_id,
                        3,
                        pts_coord_octree,
                       &n_leaf,
                       &leaf_coord_octree);

    n_layer_vtx[i_layer] = n_leaf;
    layer_vtx  [i_layer] = leaf_coord_octree;

    layer_vtx_idx[i_layer+1] = layer_vtx_idx[i_layer] + n_layer_vtx[i_layer];

    int *part_octree_to_bnd_vtx_idx = PDM_array_new_idx_from_const_stride_int(1, n_pts_octree);

    PDM_part_to_part_t *ptp_octree = PDM_part_to_part_create((const PDM_g_num_t **) &pts_gnum_octree,
                                                             (const int          *) &n_pts_octree,
                                                                                     1,
                                                             (const PDM_g_num_t **)  gnum_bnd_vtx,
                                                             (const int          *)  n_bnd_vtx,
                                                                                     n_part,
                                                             (const int         **) &part_octree_to_bnd_vtx_idx,
                                                             (const PDM_g_num_t **) &pts_gnum_octree,
                                                                                     comm);

    int      request         = 0;
    double **tmp_data_octree = NULL;

    PDM_part_to_part_reverse_iexch(ptp_octree,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   n_var,
                                   sizeof(double),
                                   NULL,
                 (const void **)   data_bnd_vtx,
                                   NULL,
                       (void ***) &tmp_data_octree,
                                  &request);

    PDM_part_to_part_reverse_iexch_wait(ptp_octree, request);

    double *leaf_data_octree = NULL;
    _mean_data_per_leaf(n_explicit_nodes,
                        n_points,
                        range,
                        leaf_id,
                        n_var,
                        tmp_data_octree[0],
                       &n_leaf,
                       &leaf_data_octree);

    data_layer_vtx[i_layer] = leaf_data_octree;

    distri_layer_vtx[i_layer] = malloc (sizeof(PDM_g_num_t) * (n_rank+1));
    PDM_distrib_compute(n_leaf, distri_layer_vtx[i_layer], -1, comm);

    m_layer_vtx   [i_layer] = distri_layer_vtx[i_layer][n_rank];
    gnum_layer_vtx[i_layer] = malloc (sizeof(PDM_g_num_t) * m_layer_vtx[i_layer]);

    for (int i_vtx = 0; i_vtx < m_layer_vtx[i_layer]; i_vtx++) {
      gnum_layer_vtx[i_layer][i_vtx] = _max_g_num + i_vtx;
    }

    _max_g_num = _max_g_num + m_layer_vtx[i_layer];

    free(tmp_data_octree[0]);
    free(tmp_data_octree);
    free(part_octree_to_bnd_vtx_idx);
    PDM_part_to_part_free(ptp_octree);
    PDM_para_octree_free(octree);

  }

  /*
   * Extend boundary parts by layers
   */

  n_bnd_vtx[n_part] = layer_vtx_idx[n_layer];

  gnum_bnd_vtx[n_part] = malloc (sizeof(PDM_g_num_t)         * n_bnd_vtx[n_part]);
  bnd_vtx     [n_part] = malloc (sizeof(double     ) * 3     * n_bnd_vtx[n_part]);
  data_bnd_vtx[n_part] = malloc (sizeof(double     ) * n_var * n_bnd_vtx[n_part]);

  int idx_write = 0;

  for (int i_layer = 0; i_layer < n_layer; i_layer++) {
    for (int i_vtx = 0; i_vtx < n_layer_vtx[i_layer]; i_vtx++) {
      gnum_bnd_vtx[n_part][  idx_write    ] = gnum_layer_vtx[i_layer][  i_vtx + distri_layer_vtx[i_layer][i_rank]];
      bnd_vtx     [n_part][3*idx_write    ] = layer_vtx     [i_layer][3*i_vtx                                    ];
      bnd_vtx     [n_part][3*idx_write + 1] = layer_vtx     [i_layer][3*i_vtx + 1                                ];
      bnd_vtx     [n_part][3*idx_write + 2] = layer_vtx     [i_layer][3*i_vtx + 2                                ];
      for (int i_var = 0; i_var < n_var; i_var++){
        data_bnd_vtx[n_part][n_var*idx_write + i_var] = data_layer_vtx[i_layer][n_var*i_vtx + i_var];
      }
      idx_write++;
    }
  }

  /*
   * Get distance from interior to closest boundary vertex (true boundary vertices)
   */

  if (i_rank == 0) {
    printf("-- Dist from interior to closest boundary vertex\n");
    fflush(stdout);
  }

  int         **n_int_to_bnd_vtx   = malloc (sizeof(int         *) * n_part );
  int         **id_int_to_bnd_vtx  = malloc (sizeof(int         *) * n_part );
  int         **int_to_bnd_vtx_idx = malloc (sizeof(int         *) * n_part );
  PDM_g_num_t **int_to_bnd_vtx     = malloc (sizeof(PDM_g_num_t *) * n_part );

  PDM_closest_point_t* clsp = PDM_closest_points_create(comm,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set(clsp,
                                      n_part,
                                      n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_closest_points_src_cloud_set(clsp,
                                     i_part,
                                     n_bnd_vtx   [i_part],
                                     bnd_vtx     [i_part],
                                     gnum_bnd_vtx[i_part]);

    PDM_closest_points_tgt_cloud_set(clsp,
                                     i_part,
                                     n_int_vtx   [i_part],
                                     int_vtx     [i_part],
                                     gnum_int_vtx[i_part]);

  }

  PDM_closest_points_compute(clsp);

  double _max_dist = 0.0;
  double  max_dist = 0.0;

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_g_num_t *closest_src_gnum = NULL;
    double      *closest_src_dist = NULL;

    PDM_closest_points_get(clsp,
                           i_part,
                          &closest_src_gnum,
                          &closest_src_dist);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      _max_dist = PDM_MAX(_max_dist, closest_src_dist[i_vtx]);
    }

  }

  PDM_MPI_Allreduce(&_max_dist, &max_dist, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  for (int i_part = 0; i_part < n_part; i_part++) {

    gnum_bnd_vtx_min_dist[i_part] = malloc (sizeof(PDM_g_num_t)     * n_int_vtx[i_part]);
    bnd_vtx_min_dist     [i_part] = malloc (sizeof(double     ) * 3 * n_int_vtx[i_part]);
    n_bnd_vtx_min_dist   [i_part] = 0;

    n_int_to_bnd_vtx  [i_part]    = malloc (sizeof(int) *  n_int_vtx[i_part]   );
    id_int_to_bnd_vtx [i_part]    = malloc (sizeof(int) *  n_int_vtx[i_part]   );
    int_to_bnd_vtx_idx[i_part]    = malloc (sizeof(int) * (n_int_vtx[i_part]+1));
    int_to_bnd_vtx_idx[i_part][0] = 0;

    PDM_g_num_t *sorted_gnum = malloc (sizeof(PDM_g_num_t) * n_int_vtx[i_part]);
    int         *order       = malloc (sizeof(int        ) * n_int_vtx[i_part]);

    PDM_g_num_t *closest_src_gnum = NULL;
    double      *closest_src_dist = NULL;

    PDM_closest_points_get(clsp,
                           i_part,
                          &closest_src_gnum,
                          &closest_src_dist);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {

      id_int_to_bnd_vtx[i_part][i_vtx] = 0;

      if (closest_src_dist[i_vtx] <= min_dist) {
        for (int j_vtx = 0; j_vtx < n_bnd_vtx_min_dist[i_part]; j_vtx++) {
          sorted_gnum[j_vtx] = gnum_bnd_vtx_min_dist[i_part][j_vtx];
          order[j_vtx]       = j_vtx;
        }
        PDM_sort_long(sorted_gnum, order, n_bnd_vtx_min_dist[i_part]);
        int i_gnum = -1;
        if (n_bnd_vtx_min_dist[i_part] > 0) {
          i_gnum = PDM_binary_search_long(closest_src_gnum[i_vtx], sorted_gnum, n_bnd_vtx_min_dist[i_part]);
        }
        n_int_to_bnd_vtx [i_part][i_vtx] = n_vtx_min_dist;
        if (i_gnum < 0) {
          gnum_bnd_vtx_min_dist[i_part][  n_bnd_vtx_min_dist[i_part]    ] = closest_src_gnum[i_vtx];
          bnd_vtx_min_dist     [i_part][3*n_bnd_vtx_min_dist[i_part]    ] = int_vtx[i_part][3*i_vtx    ];
          bnd_vtx_min_dist     [i_part][3*n_bnd_vtx_min_dist[i_part] + 1] = int_vtx[i_part][3*i_vtx + 1];
          bnd_vtx_min_dist     [i_part][3*n_bnd_vtx_min_dist[i_part] + 2] = int_vtx[i_part][3*i_vtx + 2];
          id_int_to_bnd_vtx    [i_part][i_vtx] = n_bnd_vtx_min_dist[i_part] + 1;
          n_bnd_vtx_min_dist[i_part]++;
        } else {
          id_int_to_bnd_vtx    [i_part][i_vtx] = order[i_gnum] + 1;
        }
      } else {
        for (int i_layer = 0; i_layer < n_layer; i_layer++) {
          double a_dist = min_dist + (max_dist - min_dist)/n_layer* i_layer;
          double b_dist = min_dist + (max_dist - min_dist)/n_layer*(i_layer+1);
          if (closest_src_dist[i_vtx] > a_dist && closest_src_dist[i_vtx] <= b_dist) {
            n_int_to_bnd_vtx [i_part][i_vtx] = m_layer_vtx[i_layer];
            id_int_to_bnd_vtx[i_part][i_vtx] = -(i_layer + 1);
          }
        }
      }

      assert(n_int_to_bnd_vtx [i_part][i_vtx] != 0);
      assert(id_int_to_bnd_vtx[i_part][i_vtx] != 0);

    }

    gnum_bnd_vtx_min_dist[i_part] = realloc (gnum_bnd_vtx_min_dist[i_part], sizeof(PDM_g_num_t)     * n_bnd_vtx_min_dist[i_part]);
    bnd_vtx_min_dist     [i_part] = realloc (bnd_vtx_min_dist     [i_part], sizeof(double     ) * 3 * n_bnd_vtx_min_dist[i_part]);

    free(sorted_gnum);
    free(order);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      int_to_bnd_vtx_idx[i_part][i_vtx+1] = int_to_bnd_vtx_idx[i_part][i_vtx] + n_int_to_bnd_vtx[i_part][i_vtx];
    }

    int_to_bnd_vtx[i_part] = malloc (sizeof(PDM_g_num_t) * int_to_bnd_vtx_idx[i_part][n_int_vtx[i_part]]);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {

      int idx_write2 = int_to_bnd_vtx_idx[i_part][i_vtx];

      if (id_int_to_bnd_vtx[i_part][i_vtx] < 0) {
        int i_layer = -id_int_to_bnd_vtx[i_part][i_vtx] - 1;
        for (int j_vtx = 0; j_vtx < m_layer_vtx[i_layer]; j_vtx++) {
          int_to_bnd_vtx[i_part][idx_write2++] = gnum_layer_vtx[i_layer][j_vtx];
        }
      }

    }

  }

  PDM_closest_points_free(clsp);

  /*
   * Get groups of vertices for selected boundary vertices
   */

  if (i_rank == 0) {
    printf("-- Get groups of vertices for selected boundary vertices\n");
    fflush(stdout);
  }

  _compute_gnum(n_part,
                n_bnd_vtx_min_dist,
                gnum_bnd_vtx_min_dist,
                comm);

  PDM_closest_point_t* clsp_min_dist = PDM_closest_points_create(comm,
                                                                 n_vtx_min_dist,
                                                                 PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set(clsp_min_dist,
                                      n_part,
                                      n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_closest_points_src_cloud_set(clsp_min_dist,
                                     i_part,
                                     n_bnd_vtx   [i_part],
                                     bnd_vtx     [i_part],
                                     gnum_bnd_vtx[i_part]);

    PDM_closest_points_tgt_cloud_set(clsp_min_dist,
                                     i_part,
                                     n_bnd_vtx_min_dist   [i_part],
                                     bnd_vtx_min_dist     [i_part],
                                     gnum_bnd_vtx_min_dist[i_part]);

  }

  PDM_closest_points_compute(clsp_min_dist);

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_g_num_t *closest_src_gnum = NULL;
    double      *closest_src_dist = NULL;

    PDM_closest_points_get(clsp_min_dist,
                           i_part,
                          &closest_src_gnum,
                          &closest_src_dist);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {

      if (id_int_to_bnd_vtx[i_part][i_vtx] > 0) {
        int idx_write3 = int_to_bnd_vtx_idx[i_part][i_vtx];
        int i_bnd_vtx_min_dist = id_int_to_bnd_vtx[i_part][i_vtx] - 1;
        for (int j_vtx = 0; j_vtx < n_vtx_min_dist; j_vtx++) {
          int_to_bnd_vtx[i_part][idx_write3++] = closest_src_gnum[n_vtx_min_dist*i_bnd_vtx_min_dist+j_vtx];
        }
      }

    }

  }

  PDM_closest_points_free(clsp_min_dist);

  /*
   * Reequilibrate interior vertices
   */

  if (i_rank == 0) {
    printf("-- Reequilibrate interior vertices\n");
    fflush(stdout);
  }

  double **weight_int_vtx = malloc (sizeof(double*) * n_part );
  int    **stride_3       = malloc (sizeof(int   *) * n_part );

  for (int i_part = 0; i_part < n_part; i_part++) {

    weight_int_vtx[i_part] = malloc (sizeof(double) * n_int_vtx[i_part]);
    stride_3      [i_part] = malloc (sizeof(int   ) * n_int_vtx[i_part]);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      weight_int_vtx[i_part][i_vtx] = (double) n_int_to_bnd_vtx[i_part][i_vtx];
      stride_3      [i_part][i_vtx] = 3;
    }

  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      gnum_int_vtx,
                                                      weight_int_vtx,
                                                      n_int_vtx,
                                                      n_part,
                                                      comm);

  int          blk_n_int_vtx        = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *blk_gnum_int_vtx     = PDM_part_to_block_block_gnum_get(ptb);
  int         *blk_n_int_to_bnd_vtx = NULL;
  PDM_g_num_t *blk_int_to_bnd_vtx   = NULL;
  double      *blk_int_vtx          = NULL;
  int         *blk_stride_3         = NULL;

  int blk_size = PDM_part_to_block_exch(ptb,
                                        sizeof(PDM_g_num_t),
                                        PDM_STRIDE_VAR_INTERLACED,
                                        1,
                                        n_int_to_bnd_vtx,
                             (void **)  int_to_bnd_vtx,
                                       &blk_n_int_to_bnd_vtx,
                             (void **) &blk_int_to_bnd_vtx);

  PDM_UNUSED(blk_size);

  int *blk_int_to_bnd_vtx_idx = malloc (sizeof(int) * (blk_n_int_vtx+1));

  blk_int_to_bnd_vtx_idx[0] = 0;
  for (int i_vtx = 0; i_vtx < blk_n_int_vtx; i_vtx++) {
    assert(blk_n_int_to_bnd_vtx[i_vtx] != 0);
    blk_int_to_bnd_vtx_idx[i_vtx+1] = blk_int_to_bnd_vtx_idx[i_vtx] + blk_n_int_to_bnd_vtx[i_vtx];
  }

  free(blk_n_int_to_bnd_vtx);

  blk_size = PDM_part_to_block_exch(ptb,
                                    sizeof(double),
                                    PDM_STRIDE_VAR_INTERLACED,
                                    1,
                                    stride_3,
                         (void **)  int_vtx,
                                   &blk_stride_3,
                         (void **) &blk_int_vtx);

  for (int i_vtx = 0; i_vtx < blk_n_int_vtx; i_vtx++) {
    assert(blk_stride_3[i_vtx] == 3);
  }

  /*
   * Exchange data from boundary to interior vertices
   */

  if (i_rank == 0) {
    printf("-- Exchange data from boundary to interior vertices\n");
    fflush(stdout);
  }

  PDM_part_to_part_t *ptp_int_to_bnd = PDM_part_to_part_create((const PDM_g_num_t **) &blk_gnum_int_vtx,
                                                               (const int          *) &blk_n_int_vtx,
                                                                                       1,
                                                               (const PDM_g_num_t **)  gnum_bnd_vtx,
                                                               (const int          *)  n_bnd_vtx,
                                                                                       n_part+1,
                                                               (const int         **) &blk_int_to_bnd_vtx_idx,
                                                               (const PDM_g_num_t **) &blk_int_to_bnd_vtx,
                                                                                       comm);

  int      request            = 0;
  double **data_from_bnd_vtx  = NULL;
  double **coord_from_bnd_vtx = NULL;

  PDM_part_to_part_reverse_iexch(ptp_int_to_bnd,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 3,
                                 sizeof(double),
                                 NULL,
               (const void **)   bnd_vtx,
                                 NULL,
                     (void ***) &coord_from_bnd_vtx,
                                &request);

  PDM_part_to_part_reverse_iexch_wait(ptp_int_to_bnd, request);

  PDM_part_to_part_reverse_iexch(ptp_int_to_bnd,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 n_var,
                                 sizeof(double),
                                 NULL,
               (const void **)   data_bnd_vtx,
                                 NULL,
                     (void ***) &data_from_bnd_vtx,
                                &request);

  PDM_part_to_part_reverse_iexch_wait(ptp_int_to_bnd, request);

  PDM_part_to_part_free(ptp_int_to_bnd);

  /*
   * Move vertices
   */

  if (i_rank == 0) {
    printf("-- Move vertices\n");
    fflush(stdout);
  }

  double *dx = malloc (sizeof(double) * 3);
  double *dr = malloc (sizeof(double) * 3);
  double  sdist;
  double  dist;

  PDM_g_num_t  m_vtx = 0;
  PDM_g_num_t _m_vtx = 0;
  double       aire  = 0.0;
  double      _aire  = 0.0;
  double       l1    = 0.0;
  double      _l1    = 0.0;
  double       l2    = 0.0;
  double      _l2    = 0.0;

  dx[0] = 0.0;
  dx[1] = 0.0;
  dx[2] = 0.0;

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      m_vtx++;
      dx[0] = dx[0] + all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx]    ];
      dx[1] = dx[1] + all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx] + 1];
      dx[2] = dx[2] + all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx] + 2];
    }
  }

  PDM_MPI_Allreduce (&m_vtx, &_m_vtx, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce (dx, dr, 3, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  dx[0] = dr[0]/_m_vtx;
  dx[1] = dr[1]/_m_vtx;
  dx[2] = dr[2]/_m_vtx;

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      dr[0] = dx[0] - all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx]    ];
      dr[1] = dx[1] - all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx] + 1];
      dr[2] = dx[2] - all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx] + 2];
      l1 = PDM_MAX(l1, pow(dr[0],2)+pow(dr[1],2)+pow(dr[2],2));
    }
  }

  PDM_MPI_Allreduce (&l1, &_l1, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
  l1 = sqrt(_l1);

  dx[0] = 0.0;
  dx[1] = 0.0;
  dx[2] = 0.0;
  aire  = 0;

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      aire  = aire  + data_bnd_vtx[i_part][n_var*i_vtx];
      dx[0] = dx[0] + data_bnd_vtx[i_part][n_var*i_vtx]*data_bnd_vtx[i_part][n_var*i_vtx + 1];
      dx[1] = dx[1] + data_bnd_vtx[i_part][n_var*i_vtx]*data_bnd_vtx[i_part][n_var*i_vtx + 2];
      dx[2] = dx[2] + data_bnd_vtx[i_part][n_var*i_vtx]*data_bnd_vtx[i_part][n_var*i_vtx + 3];
    }
  }

  PDM_MPI_Allreduce (dx, dr, 3, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce (&aire, &_aire, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  dx[0] = dr[0]/_aire;
  dx[1] = dr[1]/_aire;
  dx[2] = dr[2]/_aire;

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      dr[0] = dx[0] - data_bnd_vtx[i_part][n_var*i_vtx + 1];
      dr[1] = dx[1] - data_bnd_vtx[i_part][n_var*i_vtx + 2];
      dr[2] = dx[2] - data_bnd_vtx[i_part][n_var*i_vtx + 3];
      l2 = PDM_MAX(l2, pow(dr[0],2)+pow(dr[1],2)+pow(dr[2],2));
    }
  }

  PDM_MPI_Allreduce (&l2, &_l2, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
  l2 = PDM_MAX(0.1, 5.0*sqrt(_l2)/l1)*l1;

  for (int i_vtx = 0; i_vtx < blk_n_int_vtx; i_vtx++) {
    dx[0] = 0.0;
    dx[1] = 0.0;
    dx[2] = 0.0;
    sdist = 0.0;
    for (int j_vtx = blk_int_to_bnd_vtx_idx[i_vtx]; j_vtx < blk_int_to_bnd_vtx_idx[i_vtx+1]; j_vtx++) {
      dr[0] = coord_from_bnd_vtx[0][3*j_vtx    ] - blk_int_vtx[3*i_vtx    ];
      dr[1] = coord_from_bnd_vtx[0][3*j_vtx + 1] - blk_int_vtx[3*i_vtx + 1];
      dr[2] = coord_from_bnd_vtx[0][3*j_vtx + 2] - blk_int_vtx[3*i_vtx + 2];
      dist  = sqrt(pow(dr[0],2)+pow(dr[1],2)+pow(dr[2],2));
      dist  = data_from_bnd_vtx[0][n_var*j_vtx]*(pow(l1/dist,3) + pow(l2/dist,5));
      sdist = sdist + dist;
      dx[0] = dx[0] + data_from_bnd_vtx[0][n_var*j_vtx + 1]*dist;
      dx[1] = dx[1] + data_from_bnd_vtx[0][n_var*j_vtx + 2]*dist;
      dx[2] = dx[2] + data_from_bnd_vtx[0][n_var*j_vtx + 3]*dist;
    }
    blk_int_vtx[3*i_vtx    ] = blk_int_vtx[3*i_vtx    ] + dx[0]/sdist;
    blk_int_vtx[3*i_vtx + 1] = blk_int_vtx[3*i_vtx + 1] + dx[1]/sdist;
    blk_int_vtx[3*i_vtx + 2] = blk_int_vtx[3*i_vtx + 2] + dx[2]/sdist;
  }

  free(dx);
  free(dr);

  /*
   * Get interior vertices coordinates
   */

  if (i_rank == 0) {
    printf("-- Get interior vertices coordinates\n");
    fflush(stdout);
  }

  double **int_vtx_from_blk = NULL;

  PDM_part_to_block_reverse_exch(ptb,
                                 sizeof(double),
                                 PDM_STRIDE_CST_INTERLACED,
                                 3,
                                 NULL,
                     (void   *)  blk_int_vtx,
                                 NULL,
                     (void ***) &int_vtx_from_blk);

  PDM_part_to_block_free(ptb);

  for (int i_part = 0; i_part < n_part; i_part++) {

    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx]    ] = all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx]    ] + data_bnd_vtx[i_part][n_var*i_vtx + 1];
      all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx] + 1] = all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx] + 1] + data_bnd_vtx[i_part][n_var*i_vtx + 2];
      all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx] + 2] = all_vtx[i_part][3*lnum_bnd_vtx[i_part][i_vtx] + 2] + data_bnd_vtx[i_part][n_var*i_vtx + 3];
    }

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      all_vtx[i_part][3*lnum_int_vtx[i_part][i_vtx]    ] = int_vtx_from_blk[i_part][3*i_vtx    ];
      all_vtx[i_part][3*lnum_int_vtx[i_part][i_vtx] + 1] = int_vtx_from_blk[i_part][3*i_vtx + 1];
      all_vtx[i_part][3*lnum_int_vtx[i_part][i_vtx] + 2] = int_vtx_from_blk[i_part][3*i_vtx + 2];
    }

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nface_group;

    PDM_part_part_dim_get(ppart,
                          i_part,
                         &n_cell,
                         &n_face,
                         &n_face_part_bound,
                         &n_vtx,
                         &n_proc,
                         &n_total_part,
                         &scell_face,
                         &sface_vtx,
                         &sface_group,
                         &nface_group);

    int         *cell_tag                 = NULL;
    int         *cell_face_idx            = NULL;
    int         *cell_face                = NULL;
    PDM_g_num_t *cell_ln_to_gn            = NULL;
    int         *face_tag                 = NULL;
    int         *face_cell                = NULL;
    int         *face_vtx_idx             = NULL;
    int         *face_vtx                 = NULL;
    PDM_g_num_t *face_ln_to_gn            = NULL;
    int         *face_part_bound_proc_idx = NULL;
    int         *face_part_bound_part_idx = NULL;
    int         *face_part_bound          = NULL;
    int         *vtx_tag                  = NULL;
    double      *vtx                      = NULL;
    PDM_g_num_t *vtx_ln_to_gn             = NULL;
    int         *face_group_idx           = NULL;
    int         *face_group               = NULL;
    PDM_g_num_t *face_group_ln_to_gn      = NULL;

    PDM_part_part_val_get(ppart,
                          i_part,
                         &cell_tag,
                         &cell_face_idx,
                         &cell_face,
                         &cell_ln_to_gn,
                         &face_tag,
                         &face_cell,
                         &face_vtx_idx,
                         &face_vtx,
                         &face_ln_to_gn,
                         &face_part_bound_proc_idx,
                         &face_part_bound_part_idx,
                         &face_part_bound,
                         &vtx_tag,
                         &vtx,
                         &vtx_ln_to_gn,
                         &face_group_idx,
                         &face_group,
                         &face_group_ln_to_gn);

    if (show_mesh) {
      char filename[999];
      sprintf(filename, "end_face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                            n_vtx,
                            all_vtx[i_part],
                            vtx_ln_to_gn,
                            n_face,
                            face_vtx_idx,
                            face_vtx,
                            face_ln_to_gn,
                            NULL);
    }

  }

  PDM_part_free(ppart);

  /*
   * Free memory
   */

  if (i_rank == 0) {
    printf("-- Free memory\n");
    fflush(stdout);
  }

  for (int i_part = 0; i_part < n_part+1; i_part++) {
    free(gnum_bnd_vtx[i_part]);
    free(bnd_vtx     [i_part]);
    free(data_bnd_vtx[i_part]);
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    free(lnum_bnd_vtx         [i_part]);
    free(gnum_int_vtx         [i_part]);
    free(lnum_int_vtx         [i_part]);
    free(weight_int_vtx       [i_part]);
    free(int_vtx              [i_part]);
    free(int_vtx_from_blk     [i_part]);
    free(gnum_all_vtx         [i_part]);
    free(all_vtx              [i_part]);
    free(bnd_to_all_vtx       [i_part]);
    free(bnd_to_all_vtx_idx   [i_part]);
    free(gnum_bnd_vtx_min_dist[i_part]);
    free(bnd_vtx_min_dist     [i_part]);
    free(n_int_to_bnd_vtx     [i_part]);
    free(id_int_to_bnd_vtx    [i_part]);
    free(int_to_bnd_vtx_idx   [i_part]);
    free(int_to_bnd_vtx       [i_part]);
    free(stride_3             [i_part]);
  }
  for (int i_layer = 0; i_layer < n_layer; i_layer++) {
    free(gnum_layer_vtx  [i_layer]);
    free(distri_layer_vtx[i_layer]);
    free(layer_vtx       [i_layer]);
    free(data_layer_vtx  [i_layer]);
  }
  free(gnum_bnd_vtx          );
  free(lnum_bnd_vtx          );
  free(bnd_vtx               );
  free(data_bnd_vtx          );
  free(n_bnd_vtx             );
  free(gnum_int_vtx          );
  free(lnum_int_vtx          );
  free(weight_int_vtx        );
  free(int_vtx               );
  free(int_vtx_from_blk      );
  free(n_int_vtx             );
  free(gnum_all_vtx          );
  free(all_vtx               );
  free(n_all_vtx             );
  free(bnd_to_all_vtx        );
  free(bnd_to_all_vtx_idx    );
  free(gnum_bnd_vtx_min_dist );
  free(bnd_vtx_min_dist      );
  free(n_bnd_vtx_min_dist    );
  free(gnum_layer_vtx        );
  free(distri_layer_vtx      );
  free(layer_vtx             );
  free(data_layer_vtx        );
  free(n_layer_vtx           );
  free(m_layer_vtx           );
  free(layer_vtx_idx         );
  free(n_int_to_bnd_vtx      );
  free(id_int_to_bnd_vtx     );
  free(int_to_bnd_vtx_idx    );
  free(int_to_bnd_vtx        );
  free(blk_int_to_bnd_vtx_idx);
  free(blk_int_to_bnd_vtx    );
  free(blk_int_vtx           );
  free(blk_stride_3          );
  free(stride_3              );
  free(data_from_bnd_vtx [0] );
  free(coord_from_bnd_vtx[0] );
  free(data_from_bnd_vtx     );
  free(coord_from_bnd_vtx    );
  free(n_leaf_per_layer      );

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
