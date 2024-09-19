/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_overlay.h"
#include "pdm_overlay_priv.h"
#include "pdm_dbbtree.h"
#include "pdm_sort.h"
#include "pdm_box_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_hash_tab.h"
#include "pdm_dhash_tab.h"
#include "pdm_poly_clipp.h"
#include "pdm_surf_mesh_priv.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"

#include "pdm_part_to_block_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define _DOT_PRODUCT(v1, v2) \
  (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define _CROSS_PRODUCT_3D(cross_v1_v2, v1, v2) ( \
 cross_v1_v2[0] = v1[1]*v2[2] - v1[2]*v2[1],   \
 cross_v1_v2[1] = v1[2]*v2[0] - v1[0]*v2[2],   \
 cross_v1_v2[2] = v1[0]*v2[1] - v1[1]*v2[0]  )

#define _MODULE(v) \
  sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum PDM_data_t
 * \brief Type of data
 *
 */

typedef enum {

  XMIN    = 0,  /*!< XMIN */
  YMIN    = 1,  /*!< YMIN */
  ZMIN    = 2,  /*!< ZMIN */
  XMAX    = 3,  /*!< XMAX */
  YMAX    = 4,  /*!< YMAX */
  ZMAX    = 5,  /*!< ZMAX */

} _extents_t;


/**
 * \struct _sub_edge_t
 * \brief  sub edge definition to store in a hash table structure
 *
 */

typedef struct {

  PDM_g_num_t vtx1;       /*!< First Vertex */
  double     coords1[3];/*!< Vertex1 Coordinates */

  PDM_g_num_t vtx2;      /*!< second Vertex */
  double     coords2[3];/*!< Vertex 2 Coordinates */

} _sub_edge_t;


/**
 * \struct _sub_edge_t
 * \brief  sub edge definition to store in a hash table structure
 *
 */

typedef struct {

  PDM_g_num_t  vtx1;                 /*!< First Vertex */
  PDM_g_num_t  vtx2;                 /*!< second Vertex */
  int          sIntEdge;             /*!< Size of arrays */
  int          nIntEdge;             /*!< Number of sub-vertices */

  PDM_g_num_t *vtxIntEdge;
  double      *uIntEdge;
  double      *coordsIntEdge;

} _sub_vertices_origin_edge_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/**
 * \brief Dump elapsed an CPU time
 *
 *
 * \param [in]  id                PDM_ol identifier.
 *
 */

static void
_ol_dump_times
(
 PDM_ol_t *ol
)
{

  _ol_timer_step_t tfin = OL_SORT_A_B_CONNECT_GRAPH;

  ol->times_elapsed[END] = ol->times_elapsed[tfin];
  ol->times_cpu[END]     = ol->times_cpu[tfin];
  ol->times_cpu_u[END]   = ol->times_cpu_u[tfin];
  ol->times_cpu_s[END]   = ol->times_cpu_s[tfin];

  double t1;
  double t2;

  t1 = ol->times_elapsed[END] - ol->times_elapsed[BEGIN];
  t2 = ol->times_cpu[END] - ol->times_cpu[BEGIN];

  PDM_printf( "ol times ALL : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[INIT_DEF_DATA] - ol->times_elapsed[BEGIN];
  t2 = ol->times_cpu[INIT_DEF_DATA] - ol->times_cpu[BEGIN];

  // PDM_printf( "ol times INIT_DEF_DATA : %12.5e %12.5e\n", t1, t2);

  t1 += ol->times_elapsed[INIT_BUILD_EDGES] - ol->times_elapsed[INIT_DEF_DATA];
  t2 += ol->times_cpu[INIT_BUILD_EDGES] - ol->times_cpu[INIT_DEF_DATA];

  // PDM_printf( "ol times INIT_BUILD_EDGES : %12.5e %12.5e\n", t1, t2);

  t1 += ol->times_elapsed[INIT_BUILD_EXCH_GRAPH] - ol->times_elapsed[INIT_BUILD_EDGES];
  t2 += ol->times_cpu[INIT_BUILD_EXCH_GRAPH] - ol->times_cpu[INIT_BUILD_EDGES];

  // PDM_printf( "ol times INIT_BUILD_EXCH_GRAPH : %12.5e %12.5e\n", t1, t2);

  t1 += ol->times_elapsed[INIT_COMPUTE_CAR_LGTH] - ol->times_elapsed[INIT_BUILD_EXCH_GRAPH];
  t2 += ol->times_cpu[INIT_COMPUTE_CAR_LGTH] - ol->times_cpu[INIT_BUILD_EXCH_GRAPH];

  //PDM_printf( "ol times INIT_COMPUTE_CAR_LGTH : %12.5e %12.5e\n", t1, t2);
  PDM_printf( "ol times - Add elements to the meshe data structures (edges, characteristic length,...) : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_FACE_EXTENTS] - ol->times_elapsed[INIT_COMPUTE_CAR_LGTH];
  t2 = ol->times_cpu[OL_FACE_EXTENTS] - ol->times_cpu[INIT_COMPUTE_CAR_LGTH];

  // PDM_printf( "ol times OL_FACE_EXTENTS : %12.5e %12.5e\n", t1, t2);

  t1 += ol->times_elapsed[OL_IS_PLANE] - ol->times_elapsed[OL_FACE_EXTENTS];
  t2 += ol->times_cpu[OL_IS_PLANE] - ol->times_cpu[OL_FACE_EXTENTS];

  // PDM_printf( "ol times OL_IS_PLANE : %12.5e %12.5e\n", t1, t2);

  t1 += ol->times_elapsed[OL_BUILD_DBBTREE] - ol->times_elapsed[OL_IS_PLANE];
  t2 += ol->times_cpu[OL_BUILD_DBBTREE] - ol->times_cpu[OL_IS_PLANE];

  PDM_printf( "ol times - Build parallel bounding box tree with meshA faces : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_BOXES_INTERSECT] - ol->times_elapsed[OL_BUILD_DBBTREE];
  t2 = ol->times_cpu[OL_BOXES_INTERSECT] - ol->times_cpu[OL_BUILD_DBBTREE];

  PDM_printf( "ol times - Look for intersection candidates by the intersection of bounding box tree with meshB face boxes : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_DISTRIB_BOXESA_BLOCK] - ol->times_elapsed[OL_BOXES_INTERSECT];
  t2 = ol->times_cpu[OL_DISTRIB_BOXESA_BLOCK] - ol->times_cpu[OL_BOXES_INTERSECT];

  PDM_printf( "ol times - Distribute intersections from the absolute number of meshA faces to ensure a good load balancing  : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_EDGES_INTERSECTION] - ol->times_elapsed[OL_DISTRIB_BOXESA_BLOCK];
  t2 = ol->times_cpu[OL_EDGES_INTERSECTION] - ol->times_cpu[OL_DISTRIB_BOXESA_BLOCK];

  PDM_printf( "ol times - Clipp step 1 - Edge intersections : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_EDGES_SYNCHRO] - ol->times_elapsed[OL_EDGES_INTERSECTION];
  t2 = ol->times_cpu[OL_EDGES_SYNCHRO] - ol->times_cpu[OL_EDGES_INTERSECTION];

  PDM_printf( "ol times - Clipp step 2 : Synchronisation of the edge intersection to ensure robustnsess  : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_CLIPPING] - ol->times_elapsed[OL_EDGES_SYNCHRO];
  t2 = ol->times_cpu[OL_CLIPPING] - ol->times_cpu[OL_EDGES_SYNCHRO];

  PDM_printf( "ol times - Clipp step 3 : Build sub-faces : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESA] - ol->times_elapsed[OL_CLIPPING];
  t2 = ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESA] - ol->times_cpu[OL_CLIPPING];

  PDM_printf( "ol times - Build additional sub-faces of the meshA for the faces partielly covered by meshB : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_DISTRIB_BOXESB_BLOCK] - ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESA];
  t2 = ol->times_cpu[OL_DISTRIB_BOXESB_BLOCK] - ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESA];

  PDM_printf( "ol times - Distribute clipping results from the absolute number of meshB : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESB] - ol->times_elapsed[OL_DISTRIB_BOXESB_BLOCK];
  t2 = ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESB] - ol->times_cpu[OL_DISTRIB_BOXESB_BLOCK];

  PDM_printf( "ol times - Build additional sub-faces of the meshB for the faces partielly covered by meshA : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTA] - ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESB];
  t2 = ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTA] - ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESB];

  PDM_printf( "ol times - Distribute results to the initial partitionning of meshA : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTB] - ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTA];
  t2 = ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTB] - ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTA];

  PDM_printf( "ol times - Distribute results to the initial partitionning of meshB : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTA] - ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTB];
  t2 = ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTA] - ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTB];

  PDM_printf( "ol times - Compute local connectivity of new meshA taking into account sub-faces : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTB] - ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTA];
  t2 = ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTB] - ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTA];

  PDM_printf( "ol times - Compute local connectivity of new meshB taking into account sub-faces : %12.5e %12.5e\n", t1, t2);

  t1 = ol->times_elapsed[OL_UPDATE_A_B_CONNECT_GRAPH] - ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTB];
  t2 = ol->times_cpu[OL_UPDATE_A_B_CONNECT_GRAPH] - ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTB];

  //PDM_printf( "ol times OL_UPDATE_A_B_CONNECT_GRAPH : %12.5e %12.5e\n", t1, t2);

  t1 += ol->times_elapsed[OL_SORT_A_B_CONNECT_GRAPH] - ol->times_elapsed[OL_UPDATE_A_B_CONNECT_GRAPH];
  t2 += ol->times_cpu[OL_SORT_A_B_CONNECT_GRAPH] - ol->times_cpu[OL_UPDATE_A_B_CONNECT_GRAPH];

  PDM_printf( "ol times - Finalize parallel graph communication between shared sub-faces of new meshes : %12.5e %12.5e\n", t1, t2);
}

/**
 * \brief Compute edge entities
 *
 * This function defines edges of initial meshes and computes their conectivities
 *
 * \param [in]  ol      overlay objet
 *
 */

static void
_build_edges
(
 PDM_ol_t *ol
)
{

  /*
   * Build edges on each part
   */

  PDM_surf_mesh_build_edges (ol->meshA);
  PDM_surf_mesh_build_edges (ol->meshB);

  /*
   * Define global numbering
   */

  PDM_surf_mesh_build_edges_gn_and_edge_part_bound (ol->meshA);
  PDM_surf_mesh_build_edges_gn_and_edge_part_bound (ol->meshB);

}


/**
 * \brief Build communication graph between internal partitions of each initial mesh
 *
 * This function builds the communication graph between internal partitions
 * of each initial mesh
 *
 * \param [in]  ol       overlay object
 *
 */

static void
_build_exchange_graph
(
 PDM_ol_t *ol
)
{

  PDM_surf_mesh_build_exchange_graph (ol->meshA);

  PDM_surf_mesh_build_exchange_graph (ol->meshB);

}


/**
 * \brief Build ghost faces and ghost edges
 *
 * This function builds ghost faces and ghost edges
 * of each initial mesh
 *
 * \param [in]  ol       overlay object
 *
 */

static void
_compute_carLgthVtx
(
 PDM_ol_t *ol
)
{

  PDM_surf_mesh_compute_carLgthVtx (ol->meshA);

  PDM_surf_mesh_compute_carLgthVtx (ol->meshB);

}


/**
 * \brief Compute face extents
 *
 * This function computes face extents
 *
 * \param [in]  ol       overlay object
 *
 */

static void
_compute_faceExtents
(
 PDM_ol_t *ol
)
{

  PDM_surf_mesh_compute_faceExtentsMesh(ol->meshA,
                                        ol->extentsTol);

  PDM_surf_mesh_compute_faceExtentsMesh(ol->meshB,
                                        ol->extentsTol);
}


/**
 * \brief Chek if the two meshes are plane surfaces
 *
 * This function cheks the two meshes are plane surfaces
 * and returns the distance between the two meshes
 *
 * \param [in]   ol     overlay object
 * \param [out]  dist   Distance between the two planes
 *
 */

static int
_is_same_plane
(
 PDM_ol_t *ol,
 double    *dist
)
{
  double planeEquationA[4] = {0, 0, 0, 0};
  double planeEquationB[4] = {0, 0, 0, 0};
  double barycenterA[3] = {0, 0, 0};
  double barycenterB[3] = {0, 0, 0};

  int isPlaneMeshA = PDM_surf_mesh_is_plane_surface (ol->meshA,
                                                     ol->samePlaneTol,
                                                     planeEquationA,
                                                     barycenterA);

  int isPlaneMeshB = PDM_surf_mesh_is_plane_surface (ol->meshB,
                                                     ol->samePlaneTol,
                                                     planeEquationB,
                                                     barycenterB);

  int isSamePlane = isPlaneMeshA && isPlaneMeshB;

  if (isSamePlane) {
    if (   (1 - fabs (_DOT_PRODUCT (planeEquationA, planeEquationB))
         > ol->samePlaneTol)) {
      isSamePlane = 0;
    }
  }

  if (isSamePlane) {
    double _v1[3];
    double _v2[3];

    for (int i = 0; i < 3; i++){
      _v1[i] = barycenterB[i] - barycenterA[i];
      _v2[i] = -_v1[i];
    }
    double dist1 = fabs (_DOT_PRODUCT (planeEquationA, _v1));
    double dist2 = fabs (_DOT_PRODUCT (planeEquationB, _v2));
    *dist = _MAX (dist1, dist2);


  }

  return isSamePlane;
}

/**
 * \brief Compute overlay planes
 *
 * This function computes overlay mesh of two plane meshes
 *
 * \param [in]  ol       overlay object
 *
 */

static void
_compute_overlay_planes
(
 PDM_ol_t *ol
)
{
  const int vb = 0;

  //TODO: Cette fonction est très grosse. Elle sera decoupee en petites fonctions
  //      lors de l'optimisation de l'algorithme.

  int i_rank;
  PDM_MPI_Comm_rank (ol->comm, &i_rank);
  int lComm;
  PDM_MPI_Comm_size (ol->comm, &lComm);


  /*****************************************************************************
   *                                                                           *
   * Create a DBBtree structure (dbbtreeA) to store mesh A elements boxes      *
   * (boxes B)                                                                 *
   *                                                                           *
   ****************************************************************************/

  /* for (int i = 0; i < ol->meshA->n_part; i++) { */
	/*   PDM_surf_part_dump(ol->meshA->part[i]); */
  /* } */

  /* for (int i = 0; i < ol->meshA->n_part; i++) { */
	/*   PDM_surf_part_dump(ol->meshA->part[i]); */
  /* } */

  double global_extents[6] = {HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                             -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

  const int dim = 3;

  PDM_surf_mesh_t *meshA = ol->meshA;
  const int n_partA = PDM_surf_mesh_n_part_get (meshA);

  int *nEltsA = (int *) malloc (sizeof(int) * n_partA);
  const PDM_g_num_t **gNumA =
                  (const PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_partA);
  const double **extentsA = (const double **) malloc (sizeof(double *) * n_partA);

  for (int i = 0; i < n_partA; i++) {
    nEltsA[i] = PDM_surf_mesh_part_n_face_get (meshA, i);
    extentsA[i] = PDM_surf_mesh_part_extents_get (meshA, i);
    for (int j = 0; j < nEltsA[i]; j++) {
      for (int k = 0; k < dim; k++) {
        global_extents[k]       = PDM_MIN(extentsA[i][2*dim*j + k], global_extents[k]);
        global_extents[dim + k] = PDM_MAX(extentsA[i][(2*j+1) * dim + k], global_extents[dim+k]);
      }
    }
    gNumA[i] = PDM_surf_mesh_part_face_g_num_get (meshA, i);
  }

  PDM_surf_mesh_t *meshB = ol->meshB;
  const int n_partB = PDM_surf_mesh_n_part_get (meshB);

  int *nEltsB = (int *) malloc (sizeof(int) * n_partB);
  const PDM_g_num_t **gNumB = (const PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_partB);
  const double **extentsB = (const double **) malloc (sizeof(double *) * n_partB);

  for (int i = 0; i < n_partB; i++) {
    nEltsB[i] = PDM_surf_mesh_part_n_face_get (meshB, i);
    extentsB[i] = PDM_surf_mesh_part_extents_get (meshB, i);
    for (int j = 0; j < nEltsB[i]; j++) {
      for (int k = 0; k < dim; k++) {
        global_extents[k]       = PDM_MIN(extentsB[i][2*dim*j + k], global_extents[k]);
        global_extents[dim + k] = PDM_MAX(extentsB[i][(2*j+1) * dim + k], global_extents[dim+k]);
      }
    }
    gNumB[i] = PDM_surf_mesh_part_face_g_num_get (meshB, i);
  }

  double g_global_extents[6];
  PDM_MPI_Allreduce(global_extents, g_global_extents, dim,
                    PDM_MPI_DOUBLE, PDM_MPI_MIN, ol->comm);
  PDM_MPI_Allreduce(global_extents+dim, g_global_extents+dim, dim,
                    PDM_MPI_DOUBLE, PDM_MPI_MAX, ol->comm);
  double max_range = -HUGE_VAL;
  double min_range = HUGE_VAL;

  for (int k = 0; k < dim; k++) {
    max_range = PDM_MAX(max_range, (g_global_extents[dim+k] - g_global_extents[k]));
    min_range = PDM_MIN(min_range, (g_global_extents[dim+k] - g_global_extents[k]));
  }

  for (int k = 0; k < dim; k++) {
    g_global_extents[k]     += -max_range * 1.1e-3; // On casse la symetrie !
    g_global_extents[dim+k] +=  max_range * 1e-3;
  }

  PDM_dbbtree_t *dbbtreeA = PDM_dbbtree_create (ol->comm, dim, g_global_extents);

  PDM_box_set_t  *boxesA = PDM_dbbtree_boxes_set (dbbtreeA,
                                                  n_partA,
                                                  nEltsA,
                                                  extentsA,
                                                  gNumA);

  free (nEltsA);
  free (gNumA);
  free (extentsA);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_BUILD_DBBTREE] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_BUILD_DBBTREE]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_BUILD_DBBTREE]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_BUILD_DBBTREE]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  if (i_rank == 0 && vb == 1) {
    printf ("!!!! build dbbtree !!!!\n");
    fflush(stdout);
  }

  /*****************************************************************************
   *                                                                           *
   * Look for intersections between each mesh A boxes and mesh B boxes         *
   * the result is stored inside boxesB_intersection_index (size = n_eltA + 1);*
   *                             boxesB_intersection_index_l_num               *
   *                             (local number of B boxes)                     *
   *                                                                           *
   ****************************************************************************/

  int *boxesB_intersection_index;
  int *boxesB_intersection_l_num;
  PDM_box_set_t  *boxesB = PDM_dbbtree_intersect_boxes_set (dbbtreeA,
                                                            n_partB,
                                                            nEltsB,
                                                            extentsB,
                                                            gNumB,
                                                            &boxesB_intersection_index,
                                                            &boxesB_intersection_l_num);


  PDM_dbbtree_free (dbbtreeA);
  free (nEltsB);
  free (gNumB);
  free (extentsB);

  /*
   * Check boxesA and boxesB entities
   */

  const int *originA = PDM_box_set_origin_get (boxesA);

  const int *originB = PDM_box_set_origin_get (boxesB);

  int              n_eltA    = PDM_box_set_get_size (boxesA);
  int              n_eltB    = PDM_box_set_get_size (boxesB);

  PDM_g_num_t *gnum_eltA = (PDM_g_num_t *) PDM_box_set_get_g_num (boxesA);
  PDM_g_num_t *gnum_eltB = (PDM_g_num_t *) PDM_box_set_get_g_num (boxesB);

  if (1 == 0) {
    printf ("intersections A\n");
    for (int i = 0; i < n_eltA; i++) {
      printf(PDM_FMT_G_NUM" :", gnum_eltA[i]);
      for (int j = boxesB_intersection_index[i]; j < boxesB_intersection_index[i+1]; j++) {
        printf(" "PDM_FMT_G_NUM, gnum_eltB[boxesB_intersection_l_num[j]]);
      }
      printf("\n");
    }
  }

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_BOXES_INTERSECT] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_BOXES_INTERSECT]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_BOXES_INTERSECT]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_BOXES_INTERSECT]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  if (i_rank == 0  && vb == 1) {
    printf ("!!!! boxes intersect !!!!\n");
    fflush(stdout);
  }

  /*****************************************************************************
   *                                                                           *
   *  Transfer intersection information from partitions to blocks              *
   * with PDM_part_to_block_exch function                                       *
   *                                                                           *
   *  Results :                                                                *
   *      - blockA_boxesB_idx                                                  *
   *      - blockA_boxesB_gnum_data                                            *
   *                                                                           *
   ****************************************************************************/


  PDM_part_to_block_t *ptb_boxesA = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_MERGE,
                                                              1.,
                                                              (PDM_g_num_t **) &gnum_eltA,
                                                              NULL,
                                                              &n_eltA,
                                                              1,
                                                              ol->comm);

  int n_elt_blockA = PDM_part_to_block_n_elt_block_get (ptb_boxesA);

  PDM_g_num_t *block_gnumA = PDM_part_to_block_block_gnum_get (ptb_boxesA);

  int *part_strideA = (int *) malloc (sizeof(int) * n_eltA);

  for (int i = 0; i < n_eltA; i++) {
    part_strideA[i] =
            boxesB_intersection_index[i+1] - boxesB_intersection_index[i];
  }


  PDM_g_num_t *boxesB_intersection_g_num = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) *
                                                               boxesB_intersection_index[n_eltA]);

  for (int k = 0; k < boxesB_intersection_index[n_eltA]; k++) {
    boxesB_intersection_g_num[k] =  gnum_eltB[boxesB_intersection_l_num[k]];
  }

  int       *blockA_boxesB_stride;
  PDM_g_num_t *blockA_boxesB_gnum_data;

  PDM_part_to_block_exch (ptb_boxesA,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &part_strideA,
                         (void **) &boxesB_intersection_g_num,
                         &blockA_boxesB_stride,
                         (void **) &blockA_boxesB_gnum_data);

  int *blockA_boxesB_idx = PDM_array_zeros_int(n_elt_blockA + 1);

  for (int i = 0; i < n_elt_blockA; i++) {
    blockA_boxesB_idx[i+1] = blockA_boxesB_idx[i] + blockA_boxesB_stride[i];
  }

  free (blockA_boxesB_stride);

  int idx = 0;
  int *blockA_boxesB_idx_new = PDM_array_zeros_int(n_elt_blockA + 1);

  for (int i = 0; i < n_elt_blockA; i++) {

    PDM_g_num_t *ideb = blockA_boxesB_gnum_data + blockA_boxesB_idx[i];
    int length = blockA_boxesB_idx[i+1] - blockA_boxesB_idx[i];

    PDM_sort_long (ideb, NULL, length);

    PDM_g_num_t pre = -1;
    for (int k = 0; k < length; k++) {
      if (pre != ideb[k]) {
        blockA_boxesB_gnum_data[idx++] = ideb[k];
        blockA_boxesB_idx_new[i+1] += 1;
        pre = ideb[k];
      }
    }
  }

  for (int i = 0; i < n_elt_blockA; i++) {
    blockA_boxesB_idx_new[i+1] += blockA_boxesB_idx_new[i];
  }

  for (int i = 0; i < n_elt_blockA + 1; i++) {
    blockA_boxesB_idx[i] = blockA_boxesB_idx_new[i];
  }

  free (blockA_boxesB_idx_new);

  /*****************************************************************************
   *                                                                           *
   * Redistribute boxesA intersections to ensure a good load balacing          *
   * in ol->comm MPI communicator                                              *
   * This step removes intersections found many times on different ranks       *
   *                                                                           *
   * After this step, data are stored in a block with n_elt_blockA, block_gnumA,
   * part_strideA                                                              *
   *                                                                           *
   ****************************************************************************/

  /*
   * Redistribute data boxesA from blockB distribution with a PDM_box_distrib_t
   * structure
   * TODO: - Build a new PDM_box_distrib_t structure more simple
   *       - Hide PDM_box_distrib_t attributes
   */

  PDM_l_num_t *destination = PDM_part_to_block_destination_get (ptb_boxesA);

  PDM_g_num_t n_g_eltA = PDM_box_set_get_global_size(boxesA);

  PDM_box_distrib_t *distribA =
          PDM_box_distrib_create (n_eltA,
                                  n_g_eltA,
                                  1, // Don't use in this case
                                  ol->comm);


  PDM_g_num_t n_g_eltB = PDM_box_set_get_global_size(boxesB);

  PDM_box_distrib_t *distribB = PDM_box_distrib_create (n_eltB,
                                                        n_g_eltB,
                                                        1, // Don't use in this case
                                                        ol->comm);


  int *countEltsA = (int *) malloc (sizeof(int) * lComm);
  int *countEltsB = (int *) malloc (sizeof(int) * lComm);

  for (int i = 0; i < lComm + 1; i++) {
    distribA->index[i] = 0;
    distribB->index[i] = 0;
  }

  for (int i = 0; i < lComm; i++) {
    countEltsA[i] = 0;
    countEltsB[i] = 0;
  }

  for (int i = 0; i < n_eltA; i++) {
    int iProc = destination[i] + 1;
    distribA->index[iProc]++;
    distribB->index[iProc] += part_strideA[i];
  }

  for (int i = 0; i < lComm; i++) {
    distribA->index[i+1] += distribA->index[i];
    distribB->index[i+1] += distribB->index[i];
  }

  distribA->list = (int *) malloc (sizeof(int) * distribA->index[lComm]);
  distribB->list = (int *) malloc (sizeof(int) * distribB->index[lComm]);

  for (int i = 0; i < n_eltA; i++) {
    int iProc = destination[i]; // EQU + 1; mais ce n est pas necessaire
    int idxA = distribA->index[iProc] + (countEltsA[iProc]++);
    distribA->list[idxA] = i;
    int idxB = distribB->index[iProc] + countEltsB[iProc];
    countEltsB[iProc] += part_strideA[i];
    int k=0;
    for (int j = boxesB_intersection_index[i]; j < boxesB_intersection_index[i+1]; j++) {
      distribB->list[idxB+k++] = boxesB_intersection_l_num[j];
    }
  }

  PDM_g_num_t *_block_gnumA = malloc (sizeof(PDM_g_num_t) * n_elt_blockA);
  memcpy(_block_gnumA, block_gnumA, sizeof(PDM_g_num_t) * n_elt_blockA);

  PDM_part_to_block_free (ptb_boxesA);

  block_gnumA = _block_gnumA;

  free (boxesB_intersection_g_num);
  free (boxesB_intersection_l_num);
  free (boxesB_intersection_index);
  free (part_strideA);
  free (countEltsA);
  free (countEltsB);

  PDM_box_distrib_clean (distribA);
  PDM_box_distrib_clean (distribB);

  PDM_box_set_redistribute (distribA, boxesA);
  PDM_box_set_redistribute (distribB, boxesB);

  PDM_box_distrib_destroy (&distribA);
  PDM_box_distrib_destroy (&distribB);

  PDM_box_set_remove_duplicate (boxesA);
  PDM_box_set_remove_duplicate (boxesB);

  /*****************************************************************************
   *                                                                           *
   *  Transfer from origin data to compute polygon intersections               *
   * (function PDM_box_set_recv_data_from_origin_distrib)                      *
   *      - faceA->edgesA connectivities (absolute number) (mapping)           *
   *      - faceB->edgesB connectivities (absolute number) (mapping)           *
   *      - faceA->vtxA connectivities (absolute number) (mapping)             *
   *      - faceB->vtxB connectivities (absolute number) (mapping)             *
   *      - faceA->coordsVtxA (To be allocate)                                 *
   *      - faceB->coordsVtxB (To be allocate)                                 *
   *      - faceA->epsVtxA (To be allocate)                                    *
   *      - faceB->epsVtxB (To be allocate)                                    *
   *                                                                           *
   ****************************************************************************/

  int        *faceStrideCurrent[2];
  int        *faceStrideCurrent3[2];
  PDM_g_num_t *faceToEdgeCurrent[2];
  PDM_g_num_t *faceToVtxCurrent[2];
  double     *face_vtxCooCurrent[2];
  double     *face_vtxEpsCurrent[2];

  PDM_surf_mesh_t *meshes[2] = {meshA, meshB};
  PDM_box_set_t  *boxes[2] = {boxesA, boxesB};
  //PDM_g_num_t nGF[2] = {ol->meshA->nGFace, ol->meshB->nGFace};

  for (int imesh = 0; imesh < 2; imesh++) {

    size_t data_size_l = sizeof(PDM_g_num_t);
    size_t data_size_r = sizeof(double);
    PDM_surf_mesh_t *mesh = meshes[imesh];

    int n_part = PDM_surf_mesh_n_part_get (mesh);

    int         **faceStrideOrigin = (int **) malloc (sizeof(int *) * n_part);
    PDM_g_num_t  **faceToEdgeOrigin = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);
    PDM_g_num_t  **faceToVtxOrigin  = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);
    double      **face_vtxCooOrigin = (double **) malloc (sizeof(double *) * n_part);
    double      **face_vtxEpsOrigin = (double **) malloc (sizeof(double *) * n_part);

    for (int i = 0; i < n_part; i++) {
      int               nEltPart        = PDM_surf_mesh_part_n_face_get (mesh, i);
      const int        *partFaceEdgeIdx = PDM_surf_mesh_part_face_edge_idx_get (mesh, i);
      const int        *partFaceEdge    = PDM_surf_mesh_part_face_edge_get (mesh, i);
      const int        *partface_vtx     = PDM_surf_mesh_part_face_vtx_get (mesh, i);
      const PDM_g_num_t *partEdgeGNum    = PDM_surf_mesh_part_edge_g_num_get (mesh, i);
      const PDM_g_num_t *partVtxGNum     = PDM_surf_mesh_part_vtx_g_num_get (mesh, i);
      //const PDM_g_num_t *partFaceGNum     = PDM_surf_mesh_part_face_g_num_get (mesh, i);
      const double     *partVtxCoord    = PDM_surf_mesh_part_vtx_get (mesh, i);
      const double     *partVtxEps      = PDM_surf_mesh_part_carLgthVtx_get (mesh, i);

      faceStrideOrigin[i] = (int *) malloc (sizeof(int) * nEltPart);
      // numero de face local et numero d aretes global
      faceToEdgeOrigin[i] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * partFaceEdgeIdx[nEltPart]);
      faceToVtxOrigin[i]  = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) *partFaceEdgeIdx[nEltPart]);
      face_vtxCooOrigin[i] = (double *) malloc (sizeof(double) * 3 * partFaceEdgeIdx[nEltPart]);
      face_vtxEpsOrigin[i] = (double *) malloc (sizeof(double) * partFaceEdgeIdx[nEltPart]);

      int        *_faceStrideOrigin = faceStrideOrigin[i];
      PDM_g_num_t *_faceToEdgeOrigin = faceToEdgeOrigin[i];
      PDM_g_num_t *_faceToVtxOrigin  = faceToVtxOrigin[i];
      double     *_face_vtxCooOrigin = face_vtxCooOrigin[i];
      double     *_face_vtxEpsOrigin = face_vtxEpsOrigin[i];

      for (int j = 0; j < nEltPart; j++) {
        _faceStrideOrigin[j] = partFaceEdgeIdx[j+1] - partFaceEdgeIdx[j];
      }

      for (int j = 0; j < partFaceEdgeIdx[nEltPart]; j++) {
        int edge1 = PDM_ABS(partFaceEdge[j]) - 1;
        int vtx1 = partface_vtx[j] - 1;
        _faceToEdgeOrigin[j] = partFaceEdge[j]/PDM_ABS(partFaceEdge[j]) * partEdgeGNum[edge1];
        _faceToVtxOrigin[j] = partVtxGNum[vtx1];
        _face_vtxEpsOrigin[j] = partVtxEps[vtx1];
        for (int k = 0; k < 3; k++) {
          _face_vtxCooOrigin[3*j+k] = partVtxCoord[3*vtx1+k];
        }
      }
    }

    faceStrideCurrent[imesh]  = NULL;
    faceStrideCurrent3[imesh] = NULL;
    faceToEdgeCurrent[imesh]  = NULL;
    faceToVtxCurrent[imesh]   = NULL;
    face_vtxCooCurrent[imesh]  = NULL;
    face_vtxEpsCurrent[imesh]  = NULL;

    PDM_box_set_recv_data_from_origin_distrib (boxes[imesh],
                                               PDM_STRIDE_VAR_INTERLACED,
                                               1,
                                               data_size_l,
                                               faceStrideOrigin,
                                               (void **) faceToEdgeOrigin,
                                               &(faceStrideCurrent[imesh]),
                                               (void **) &faceToEdgeCurrent[imesh]);

    PDM_box_set_recv_data_from_origin_distrib (boxes[imesh],
                                               PDM_STRIDE_VAR_INTERLACED,
                                               1,
                                               data_size_l,
                                               faceStrideOrigin,
                                               (void **) faceToVtxOrigin,
                                               &(faceStrideCurrent[imesh]),
                                               (void **) &faceToVtxCurrent[imesh]);

    PDM_box_set_recv_data_from_origin_distrib (boxes[imesh],
                                               PDM_STRIDE_VAR_INTERLACED,
                                               1,
                                               data_size_r,
                                               faceStrideOrigin,
                                               (void **) face_vtxEpsOrigin,
                                               &(faceStrideCurrent[imesh]),
                                               (void **) &face_vtxEpsCurrent[imesh]);

    for (int i = 0; i < n_part; i++) {
      int nEltPart        = PDM_surf_mesh_part_n_face_get (mesh, i);
      for (int j = 0; j < nEltPart; j++) {
        faceStrideOrigin[i][j] = 3 * faceStrideOrigin[i][j];
      }
    }

    PDM_box_set_recv_data_from_origin_distrib (boxes[imesh],
                                               PDM_STRIDE_VAR_INTERLACED,
                                               1,
                                               data_size_r,
                                               faceStrideOrigin,
                                               (void **) face_vtxCooOrigin,
                                               &(faceStrideCurrent3[imesh]),
                                               (void **) &face_vtxCooCurrent[imesh]);

    for (int i = 0; i < n_part; i++) {
      free (faceStrideOrigin[i]);
      free (faceToEdgeOrigin[i]);
      free (faceToVtxOrigin[i]);
      free (face_vtxCooOrigin[i]);
      free (face_vtxEpsOrigin[i]);
    }

    free (faceStrideOrigin);
    free (faceToEdgeOrigin);
    free (faceToVtxOrigin);
    free (face_vtxCooOrigin);
    free (face_vtxEpsOrigin);

  }

  /*****************************************************************************
   *                                                                           *
   *   Find BoxesB local number to build : blockA_boxesB_lnum_data             *
   *                                                                           *
   *  Steps :                                                                  *
   *      - Sort Boxes B gnum                                                  *
   *      - Binary search in sorted boxes B gnum array                         *
   *                                                                           *
   ****************************************************************************/

  int *blockA_boxesB_lnum_data =
    (int *) malloc (sizeof(int) * blockA_boxesB_idx[n_elt_blockA]);

  gnum_eltB = (PDM_g_num_t *) PDM_box_set_get_g_num (boxesB);

  n_eltB    = PDM_box_set_get_size (boxesB);

  int *lnum = (int *) malloc (sizeof(int) * n_eltB);
  PDM_g_num_t *gnum_eltB_cp = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_eltB);

  for (int i = 0; i < n_eltB; i++) {
    lnum[i] = i + 1;
    gnum_eltB_cp[i] = gnum_eltB[i];
  }

  PDM_sort_long (gnum_eltB_cp, lnum, n_eltB);

  for (int i = 0; i <  blockA_boxesB_idx[n_elt_blockA]; i++) {
    PDM_g_num_t box_gnum = blockA_boxesB_gnum_data[i];
    int idx1 = PDM_binary_search_long (box_gnum, gnum_eltB_cp, n_eltB);

    blockA_boxesB_lnum_data[i] = lnum[idx1];
  }

  free (lnum);
  free (gnum_eltB_cp);

  gnum_eltA = (PDM_g_num_t *) PDM_box_set_get_g_num (boxesA);

  n_eltA    = PDM_box_set_get_size (boxesA);

  int *blockA_lnum_data =
    (int *) malloc (sizeof(int) * n_elt_blockA);

  PDM_g_num_t *gnum_eltA_cp = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_eltA);

  for (int i = 0; i < n_eltA; i++) {
    blockA_lnum_data[i] = i + 1;
    gnum_eltA_cp[i] = gnum_eltA[i];
  }

  PDM_sort_long (gnum_eltA_cp, blockA_lnum_data, n_eltA);


  free (gnum_eltA_cp);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_DISTRIB_BOXESA_BLOCK] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_DISTRIB_BOXESA_BLOCK]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_DISTRIB_BOXESA_BLOCK]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_DISTRIB_BOXESA_BLOCK]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  if (i_rank == 0  && vb == 1) {
    printf ("!!!! ol_distrib_boxesa_block !!!!\n");
    fflush(stdout);
  }

  /*****************************************************************************
   *                                                                           *
   *   Perform polygon clippings                                               *
   *      - Use :                                                              *
   *              * faceStrideCurrent                                          *
   *              * faceStrideCurrent3                                         *
   *              * faceToEdgeCurrent                                          *
   *              * faceToVtxCurrent                                           *
   *              * face_vtxCooCurrent                                         *
   *              * face_vtxEpsCurrent                                         *
   *                                                                           *
   ****************************************************************************/

  /* Stocker les intersection dans une structure */

  int *faceIdxCurrent[2];
  faceIdxCurrent[0] = PDM_array_new_idx_from_sizes_int(faceStrideCurrent[0], n_elt_blockA);
  faceIdxCurrent[1] = PDM_array_new_idx_from_sizes_int(faceStrideCurrent[1], n_eltB);

  /*****************************************************************************
   *                                                                           *
   *                  Pre-clipping : Perform intersections                     *
   *                                                                           *
   ****************************************************************************/

  PDM_g_num_t maxGNEdgeA = PDM_surf_mesh_n_g_edge_get (meshA);
  PDM_g_num_t maxGNEdgeB = PDM_surf_mesh_n_g_edge_get (meshB);

  PDM_edges_intersect_t *intersect = PDM_edges_intersect_create (maxGNEdgeA,
                                                                 maxGNEdgeB,
                                                                 ol->vtxCarLengthTol,
                                                                 ol->comm);

  for (int i = 0; i < n_elt_blockA; i++) {
    int lnum_boxA = blockA_lnum_data[i] - 1;
    int n_vtxA = faceIdxCurrent[0][lnum_boxA  + 1] - faceIdxCurrent[0][lnum_boxA ];
    PDM_g_num_t *_faceToEdgeA = faceToEdgeCurrent[0] + faceIdxCurrent[0][lnum_boxA ];
    PDM_g_num_t *_faceToVtxA  = faceToVtxCurrent[0] + faceIdxCurrent[0][lnum_boxA ];
    double     *_face_vtxCooA = face_vtxCooCurrent[0] + (3 * faceIdxCurrent[0][lnum_boxA]);
    double     *_face_vtxEpsA = face_vtxEpsCurrent[0] + faceIdxCurrent[0][lnum_boxA ];
    //PDM_g_num_t  gnum_boxA = block_gnumA[i];

    /*
     * Perform each edge-edge intersection :
     *      - Each elementA edge intersects each elementB edge
     *      - Storage in a hash table the result of each intersection
     *        (key = sum of global number of vertices)
     *      - Build sub face element
     */

    for (int j = blockA_boxesB_idx[i]; j < blockA_boxesB_idx[i+1]; j++) {
      int lnum_boxB = blockA_boxesB_lnum_data[j]-1;
      //PDM_g_num_t gnum_boxB = blockA_boxesB_gnum_data[j];

      int n_vtxB = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
      PDM_g_num_t *_faceToEdgeB = faceToEdgeCurrent[1] + faceIdxCurrent[1][lnum_boxB];
      PDM_g_num_t *_faceToVtxB  = faceToVtxCurrent[1] + faceIdxCurrent[1][lnum_boxB];
  	  double     *_face_vtxCooB = face_vtxCooCurrent[1] + (3 * faceIdxCurrent[1][lnum_boxB]);
      double     *_face_vtxEpsB = face_vtxEpsCurrent[1] + faceIdxCurrent[1][lnum_boxB];

      /*
       * Build polygon structure according to clipping polygon structure
       */

      PDM_edges_intersect_poly_add (intersect,
                                    n_vtxA,
                                    _faceToEdgeA,
                                    _faceToVtxA,
                                    _face_vtxCooA,
                                    _face_vtxEpsA,
                                    n_vtxB,
                                    _faceToEdgeB,
                                    _faceToVtxB,
                                    _face_vtxCooB,
                                    _face_vtxEpsB);

    }
  }

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_EDGES_INTERSECTION] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_EDGES_INTERSECTION]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_EDGES_INTERSECTION]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_EDGES_INTERSECTION]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  if (i_rank == 0  && vb == 1) {
    printf ("!!!! ol_edges_intersection !!!!\n");
    fflush(stdout);
  }

  /*****************************************************************************
   *                                                                           *
   *  Synchronize intersections                                                *
   *                                                                           *
   *      Synchronize function Defines global number for :                     *
   *      - new A vertices in A overlay mesh                                   *
   *      - new B vertices in B overlay mesh                                   *
   *                                                                           *
   ****************************************************************************/

  PDM_g_num_t n_g_vtxA = PDM_surf_mesh_n_g_vtx_get (meshA);
  PDM_g_num_t n_g_vtxB = PDM_surf_mesh_n_g_vtx_get (meshB);

  PDM_g_num_t n_g_newVtxA = 0;
  PDM_g_num_t n_g_newVtxB = 0;

  PDM_edges_intersect_synchronize (intersect, n_g_vtxA, n_g_vtxB,
                                   &n_g_newVtxA, &n_g_newVtxB);



  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_EDGES_SYNCHRO] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_EDGES_SYNCHRO]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_EDGES_SYNCHRO]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_EDGES_SYNCHRO]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_edges_synchro !!!!\n");
    fflush(stdout);
  }

  /*****************************************************************************
   *                                                                           *
   *                             Clipping                                      *
   *                                                                           *
   ****************************************************************************/

  int iclipp = 0;

  int s_subFacesToFaces = 5 * n_elt_blockA;
  PDM_g_num_t *subFacesToFaces = malloc (sizeof(PDM_g_num_t) * s_subFacesToFaces);

  int s_subFacesConnecIdx = (1 + n_elt_blockA);
  int *subFacesConnecIdx = malloc (sizeof(int) * s_subFacesConnecIdx);
  subFacesConnecIdx[0] = 0;

  int s_subFacesConnecA = 4 * n_elt_blockA;
  int s_subFacesConnecB = 4 * n_elt_blockA;
  PDM_g_num_t *subFacesConnecA = malloc (sizeof(PDM_g_num_t) * s_subFacesConnecA);
  int s_gNumSubFacesA = n_elt_blockA;
  PDM_g_num_t *gNumSubFacesA   = malloc (sizeof(PDM_g_num_t) * s_gNumSubFacesA);
  PDM_g_num_t *subFacesConnecB = malloc (sizeof(PDM_g_num_t) * s_subFacesConnecB);

  int s_subFacesCoordsA = 3 * s_subFacesConnecA;
  double *subFacesCoordsA = malloc (sizeof(double) * s_subFacesCoordsA);

  int *facesToSubFacesAIdx = malloc (sizeof(int) * (1 + n_elt_blockA));
  int *facesToSubFacesBIdx = malloc (sizeof(int) * (1 + blockA_boxesB_idx[n_elt_blockA]));

  facesToSubFacesAIdx[0] = 0;
  facesToSubFacesBIdx[0] = 0;

  originB = PDM_box_set_origin_get (boxesB);

  for (int i = 0; i < n_elt_blockA; i++) {

    int lnum_boxA = blockA_lnum_data[i] - 1;

    PDM_g_num_t gnum_boxA = block_gnumA[i];

    int n_vtxA = faceIdxCurrent[0][lnum_boxA + 1] - faceIdxCurrent[0][lnum_boxA];
    PDM_g_num_t *_faceToEdgeA = faceToEdgeCurrent[0] + faceIdxCurrent[0][lnum_boxA];
    PDM_g_num_t *_faceToVtxA  = faceToVtxCurrent[0] + faceIdxCurrent[0][lnum_boxA];
    double     *_face_vtxCooA = face_vtxCooCurrent[0] + 3 * faceIdxCurrent[0][lnum_boxA];

    facesToSubFacesAIdx[i+1] = facesToSubFacesAIdx[i];

    for (int j = blockA_boxesB_idx[i]; j < blockA_boxesB_idx[i+1]; j++) {

      PDM_g_num_t gnum_boxB = blockA_boxesB_gnum_data[j];

      int lnum_boxB = blockA_boxesB_lnum_data[j] - 1;
      int n_vtxB = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
      PDM_g_num_t *_faceToEdgeB = faceToEdgeCurrent[1] + faceIdxCurrent[1][lnum_boxB];
      PDM_g_num_t *_faceToVtxB  = faceToVtxCurrent[1] + faceIdxCurrent[1][lnum_boxB];
      double     *_face_vtxCooB = face_vtxCooCurrent[1] + 3 * faceIdxCurrent[1][lnum_boxB];

      int         nPolyClippA      = 0;
      int        *polyClippIdxA    = NULL;
      PDM_g_num_t *polyClippConnecA = NULL;
      double     *polyClippCoordsA = NULL;

      int         nPolyClippB      = 0;
      int        *polyClippIdxB    = NULL;
      PDM_g_num_t *polyClippConnecB = NULL;
      double     *polyClippCoordsB = NULL;

      if (0 && (gnum_boxA == 1) && (gnum_boxB == 504101)) {

        FILE *f = fopen("poly_unit.case", "w");

        fprintf(f,
                "FORMAT\n"
                "type: ensight gold\n");

        /* Output geometry */

        fprintf(f,
                "\n"
                "GEOMETRY\n");

        fprintf(f, "model: poly_unit.geo\n");

        fclose(f);

        f = fopen("poly_unit.geo", "w");
        fprintf(f,"poly_unit\n");
        fprintf(f,"Sortie par CEDRE (V ?.?.?.?)\n");
        fprintf(f,"node id given\n");
        fprintf(f,"element id given\n");
        fprintf(f,"part\n");
        fprintf(f,"%10d\n", 1);
        fprintf(f,"poly1\n");
        fprintf(f,"coordinates\n");
        fprintf(f,"%10d\n", n_vtxA);
        for (int k = 0; k < n_vtxA; k++) {
          fprintf(f,"%10d\n", (int)_faceToVtxA[k]);
        }
        for (int k = 0; k < n_vtxA; k++) {
          fprintf(f,"%12.5e\n", _face_vtxCooA[3*k]);
        }
        for (int k = 0; k < n_vtxA; k++) {
          fprintf(f,"%12.5e\n", _face_vtxCooA[3*k+1]);
        }
        for (int k = 0; k < n_vtxA; k++) {
          fprintf(f,"%12.5e\n", _face_vtxCooA[3*k+2]);
        }
        fprintf(f,"nsided\n");
        fprintf(f,"%10d\n", 1);
        fprintf(f,"%10d\n", (int) gnum_boxA);
        fprintf(f,"%10d\n", n_vtxA);
        for (int k = 0; k < n_vtxA; k++) {
          //          fprintf(f,"%10d", (int) _faceToVtxA[k]);
          fprintf(f,"%10d", (int) (k+1));
        }
        fprintf(f,"\n");

        fprintf(f,"part\n");
        fprintf(f,"%10d\n", 2);
        fprintf(f,"poly2\n");
        fprintf(f,"coordinates\n");
        fprintf(f,"%10d\n", n_vtxB);
        for (int k = 0; k < n_vtxB; k++) {
          fprintf(f,"%10d\n", (int) _faceToVtxB[k]);
        }
        for (int k = 0; k < n_vtxB; k++) {
          fprintf(f,"%12.5e\n", _face_vtxCooB[3*k]);
        }
        for (int k = 0; k < n_vtxB; k++) {
          fprintf(f,"%12.5e\n", _face_vtxCooB[3*k+1]);
        }
        for (int k = 0; k < n_vtxB; k++) {
          fprintf(f,"%12.5e\n", _face_vtxCooB[3*k+2]);
        }
        fprintf(f,"nsided\n");
        fprintf(f,"%10d\n", 1);
        fprintf(f,"%10d\n", (int) gnum_boxB);
        fprintf(f,"%10d\n", n_vtxB);
        for (int k = 0; k < n_vtxB; k++) {
          fprintf(f,"%10d", (int) (k+1));
          //          fprintf(f,"%10d", (int) _faceToVtxB[k]);
        }
        fprintf(f,"\n");
        fclose(f);

        for (int k1 = 0; k1 < n_vtxA; k1++) {
          for (int k2 = 0; k2 < n_vtxB; k2++) {
            int n_intersect;
            PDM_edges_intersect_res_t **eir1 = PDM_edges_intersect_get (intersect,
                                                                        PDM_EDGES_GET_FROM_AB,
                                                                        PDM_ABS(_faceToEdgeA[k1]),
                                                                        PDM_ABS(_faceToEdgeB[k2]),
                                                                        &n_intersect);
            if (eir1 != NULL) {
              PDM_edges_intersect_res_dump(*eir1);
            }
          }
        }
      }

      PDM_poly_clipp (intersect,
                      gnum_boxA,
                      gnum_boxB,
                      n_vtxA,
                      _faceToEdgeA,
                      _faceToVtxA,
                      _face_vtxCooA,
                      n_vtxB,
                      _faceToEdgeB,
                      _faceToVtxB,
                      _face_vtxCooB,
                      PDM_POLY_CLIPP_CLIP,
                      &nPolyClippA,
                      &polyClippIdxA,
                      &polyClippConnecA,
                      &polyClippCoordsA,
                      &nPolyClippB,
                      &polyClippIdxB,
                      &polyClippConnecB,
                      &polyClippCoordsB);

      if (vb) {
        printf(" ---> poly_clipp gnumEltA gnumEltB : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" fin\n", gnum_boxA, gnum_boxB);
      }
      if (nPolyClippA != nPolyClippB) {
         printf("error poly_clipp gnumEltA gnumEltB : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n", gnum_boxA, gnum_boxB);
      }
      assert (nPolyClippA == nPolyClippB);

      facesToSubFacesAIdx[i+1] += nPolyClippA;

      facesToSubFacesBIdx[j+1] = nPolyClippB;

      int newSize = iclipp + nPolyClippA + 1;
      if (newSize > s_subFacesConnecIdx) {
        while (newSize > s_subFacesConnecIdx) {
          s_subFacesConnecIdx += PDM_MAX(1, s_subFacesConnecIdx/3);
        }
        subFacesConnecIdx = realloc (subFacesConnecIdx,
                                  sizeof(int) * s_subFacesConnecIdx);
      }

      newSize = iclipp + nPolyClippA;
      if (newSize > s_gNumSubFacesA) {
        while (newSize > s_gNumSubFacesA) {
          s_gNumSubFacesA += PDM_MAX(1, s_gNumSubFacesA/3);
        }
        gNumSubFacesA  = realloc (gNumSubFacesA,
                                  sizeof(PDM_g_num_t) * s_gNumSubFacesA);
      }

      int ibeg = subFacesConnecIdx[iclipp];
      int ibeg2 = 3 *ibeg;
      newSize = ibeg + polyClippIdxA[nPolyClippA];
      if (newSize > s_subFacesConnecA) {
        while (newSize > s_subFacesConnecA) {
          s_subFacesConnecA += PDM_MAX(1, s_subFacesConnecA/3);
        }
        subFacesConnecA = realloc (subFacesConnecA,
                                  sizeof(PDM_g_num_t) * s_subFacesConnecA);
      }

      if (newSize > s_subFacesConnecB) {
        while (newSize > s_subFacesConnecB) {
          s_subFacesConnecB += PDM_MAX(1, s_subFacesConnecB/3);
        }
        subFacesConnecB = realloc (subFacesConnecB,
                                  sizeof(PDM_g_num_t) * s_subFacesConnecB);
      }

      newSize = 3 * (ibeg + polyClippIdxA[nPolyClippA]);
      if (newSize > s_subFacesCoordsA) {
        while (newSize > s_subFacesCoordsA) {
          s_subFacesCoordsA += PDM_MAX(1, s_subFacesCoordsA/3);
        }
        subFacesCoordsA = realloc (subFacesCoordsA,
                                  sizeof(double) * s_subFacesCoordsA);
      }

      newSize = 5*(iclipp + nPolyClippA);
      if (newSize > s_subFacesToFaces) {
        while (newSize > s_subFacesToFaces) {
          s_subFacesToFaces += PDM_MAX(1, s_subFacesToFaces/3);
        }
        subFacesToFaces = realloc (subFacesToFaces,
                                  sizeof(PDM_g_num_t) * s_subFacesToFaces);
      }

      for (int k = 0; k < nPolyClippA; k++) {
        subFacesConnecIdx[iclipp+1] =  subFacesConnecIdx[iclipp]
                                     + (polyClippIdxA[k+1] - polyClippIdxA[k]);

        for (int k1 = polyClippIdxA[k]; k1 < polyClippIdxA[k+1]; k1++) {
          subFacesConnecA[ibeg] = polyClippConnecA[k1];
          subFacesConnecB[ibeg] = polyClippConnecB[k1];
          ibeg += 1;
        }

        for (int k1 = 3*polyClippIdxA[k]; k1 < 3*polyClippIdxA[k+1]; k1++) {
          subFacesCoordsA[ibeg2++] = polyClippCoordsA[k1];
        }

        subFacesToFaces[5*iclipp ]  = i;         // A local number
        subFacesToFaces[5*iclipp+1] = gnum_boxA; // A global number
        subFacesToFaces[5*iclipp+2] = originB[3*lnum_boxB];         // B local number
        subFacesToFaces[5*iclipp+3] = originB[3*lnum_boxB+1];         // B local number
        subFacesToFaces[5*iclipp+4] = gnum_boxB; // B global number

        //printf("link "PDM_FMT_G_NUM" : %d %d "PDM_FMT_G_NUM"\n", gnum_boxA, originB[3*lnum_boxB], originB[3*lnum_boxB+1], gnum_boxB);

        iclipp += 1;
      }

      if (polyClippIdxB != polyClippIdxA) {
        free (polyClippIdxB);
        polyClippIdxB = NULL;
      }

      if (polyClippConnecB != NULL) {
        free (polyClippConnecB);
        polyClippConnecB = NULL;
      }

      if (polyClippIdxA != NULL) {
        free (polyClippIdxA);
      }

      if (polyClippConnecA != NULL) {
        free (polyClippConnecA);
      }
      if (polyClippCoordsA != NULL) {
        free (polyClippCoordsA);
      }

    }
  }

  const int nSharedSubFaces = iclipp;

  /*
   * Global numbering Definition of sub-faces
   * This numbering is shared between A and B
   * Partial split faces are stored behind sub-faces
   * Not split faces are stored behind partial split faces
   *
   */

  PDM_g_num_t nSharedSubFaces_l = nSharedSubFaces;
  PDM_g_num_t beg_nSharedSubFaces = -1;

  PDM_MPI_Scan(&nSharedSubFaces_l, &beg_nSharedSubFaces,
           1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  beg_nSharedSubFaces += -nSharedSubFaces;

  for (int i = 0; i < nSharedSubFaces; i++) {
    gNumSubFacesA[i] = beg_nSharedSubFaces + i + 1;
  }

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_CLIPPING] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_CLIPPING]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_CLIPPING]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_CLIPPING]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  if (i_rank == 0 && vb == 1) {
    printf ("!!!! clipping !!!!\n");
    fflush(stdout);
  }

  /*****************************************************************************
   *                                                                           *
   * Compute global number for points of B into A and point of A into B        *
   *                                                                           *
   ****************************************************************************/

  if (1 == 0) {
    printf("A sub-faces before gnum:\n");
    for (int i = 0; i < nSharedSubFaces; i++) {
      printf ("%d :", i);
      for (int j = subFacesConnecIdx[i]; j < subFacesConnecIdx[i+1]; j++) {
      printf (" "PDM_FMT_G_NUM, subFacesConnecA[j]);
      }
      printf("\n");
    }

    printf("B sub-faces before gnum:\n");
    for (int i = 0; i < nSharedSubFaces; i++) {
      printf ("%d :", i);
      for (int j = subFacesConnecIdx[i]; j < subFacesConnecIdx[i+1]; j++) {
        printf (" "PDM_FMT_G_NUM, subFacesConnecB[j]);
      }
      printf("\n");
    }
  }

  const int          n_part = 1;
  const PDM_bool_t   merge = PDM_FALSE;
  const double       tolerance = 1e-8;

  PDM_gen_gnum_t *gnum_B_into_A = PDM_gnum_create (dim, n_part, merge, tolerance, ol->comm, PDM_OWNERSHIP_KEEP);

  int n_elt_B_into_A = 0;
  for (int i = 0; i < subFacesConnecIdx[nSharedSubFaces]; i++) {
    if (subFacesConnecA[i] < 0) {
      n_elt_B_into_A += 1;
    }
  }

  PDM_g_num_t *elt_B_into_A_g_num = malloc (sizeof(PDM_g_num_t) * n_elt_B_into_A);
  n_elt_B_into_A = 0;
  for (int i = 0; i < subFacesConnecIdx[nSharedSubFaces]; i++) {
    if (subFacesConnecA[i] < 0) {
      elt_B_into_A_g_num[n_elt_B_into_A++] = -subFacesConnecA[i];
    }
  }

  PDM_gnum_set_from_parents (gnum_B_into_A, 0, n_elt_B_into_A, elt_B_into_A_g_num);

  PDM_gnum_compute (gnum_B_into_A);

  PDM_g_num_t *elt_B_into_A_g_num_compress = PDM_gnum_get (gnum_B_into_A, 0);

  n_elt_B_into_A = 0;
  PDM_g_num_t _max = n_g_newVtxA;
  for (int i = 0; i < subFacesConnecIdx[nSharedSubFaces]; i++) {
    if (subFacesConnecA[i] < 0) {
      subFacesConnecA[i] = n_g_newVtxA + elt_B_into_A_g_num_compress[n_elt_B_into_A++];
      _max = PDM_MAX (_max, subFacesConnecA[i]);
    }
  }

  PDM_MPI_Allreduce(&_max, &n_g_newVtxA, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ol->comm);

  PDM_gnum_free (gnum_B_into_A);
  free (elt_B_into_A_g_num);

  PDM_gen_gnum_t *gnum_A_into_B = PDM_gnum_create (dim, n_part, merge, tolerance, ol->comm, PDM_OWNERSHIP_KEEP);

  int n_elt_A_into_B = 0;
  for (int i = 0; i < subFacesConnecIdx[nSharedSubFaces]; i++) {
    if (subFacesConnecB[i] < 0) {
      n_elt_A_into_B += 1;
    }
  }

  PDM_g_num_t *elt_A_into_B_g_num = malloc (sizeof(PDM_g_num_t) * n_elt_A_into_B);
  n_elt_A_into_B = 0;
  for (int i = 0; i < subFacesConnecIdx[nSharedSubFaces]; i++) {
    if (subFacesConnecB[i] < 0) {
      elt_A_into_B_g_num[n_elt_A_into_B++] = -subFacesConnecB[i];
    }
  }

  PDM_gnum_set_from_parents (gnum_A_into_B, 0, n_elt_A_into_B, elt_A_into_B_g_num);

  PDM_gnum_compute (gnum_A_into_B);

  PDM_g_num_t *elt_A_into_B_g_num_compress = PDM_gnum_get (gnum_A_into_B, 0);

  n_elt_A_into_B = 0;
  _max =  n_g_newVtxB;
  for (int i = 0; i < subFacesConnecIdx[nSharedSubFaces]; i++) {
    if (subFacesConnecB[i] < 0) {
      subFacesConnecB[i] = n_g_newVtxB + elt_A_into_B_g_num_compress[n_elt_A_into_B++];
      _max = PDM_MAX (_max, subFacesConnecB[i]);
    }
  }

  PDM_MPI_Allreduce(&_max, &n_g_newVtxB, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ol->comm);

  PDM_gnum_free (gnum_A_into_B);
  free (elt_A_into_B_g_num);

  if (1 == 0) {
    printf("\n-- Liste des sous-facettes des elements de A (intersection)\n");
    for (int i = 0; i < n_elt_blockA; i++) {

      PDM_g_num_t gnum_boxA = block_gnumA[i];

      printf ("face origin A : "PDM_FMT_G_NUM"\n", gnum_boxA);

      for (int j = facesToSubFacesAIdx[i]; j < facesToSubFacesAIdx[i+1]; j++) {
        printf (" - sub face %d:", j);
        for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
          printf(" "PDM_FMT_G_NUM, subFacesConnecA[k]);
        }
        printf("\n");
      }
    }
  }

  /* printf("\n-- Liste des sous-facettes des elements de B\n"); */
  /* for (int i = 0; i < n_elt_blockB; i++) { */

  /*   PDM_g_num_t gnum_boxB = block_gnumB[i]; */

  /*   printf ("face origin B : "PDM_FMT_G_NUM"\n", gnum_boxB); */

  /*   for (int j = facesToSubFacesBIdx[i]; j < facesToSubFacesBIdx[i+1]; j++) { */
  /*     printf (" - sub face %d:", j); */
  /*     for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) { */
  /*       printf(" "PDM_FMT_G_NUM"", subFacesConnecB[k]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /* } */

  /* printf("A sub-faces after gnum:\n"); */
  /* for (int i = 0; i < nSharedSubFaces; i++) { */
  /*   printf ("%d :", i); */
  /*   for (int j = subFacesConnecIdx[i]; j < subFacesConnecIdx[i+1]; j++) { */
  /*     printf (" "PDM_FMT_G_NUM, subFacesConnecA[j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* printf("B sub-faces after gnum:\n"); */
  /* for (int i = 0; i < nSharedSubFaces; i++) { */
  /*   printf ("%d :", i); */
  /*   for (int j = subFacesConnecIdx[i]; j < subFacesConnecIdx[i+1]; j++) { */
  /*     printf (" "PDM_FMT_G_NUM, subFacesConnecB[j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  //
  // TODO: Optimisation : Commencer les communications du part_to_block
  //                      non-bloquantes ici
  //

  /*****************************************************************************
   *                                                                           *
   * Compute sub-faces coming from A faces with partial B cover.               *
   *                                                                           *
   * Check that each sub-edge is referenced twice                              *
   *    - Take into sub-edges of sub-faces                                     *
   *    - Take into sub-edges of origin elements                               *
   *                                                                           *
   * If sub-edges are not twice referenced, the original face is partial       *
   * covered by B. Additional sub-faces have to compute.                       *
   *                                                                           *
   ****************************************************************************/

  int *facesToAddSubFacesAIdx = malloc (sizeof(int) * (1 + n_elt_blockA));
  facesToAddSubFacesAIdx[0] = facesToSubFacesAIdx[n_elt_blockA];

  int *faceIniVtxIdxA = malloc (sizeof(int) * (1 + n_elt_blockA));
  faceIniVtxIdxA[0] = 0;

  int s_faceIniVtxA = 4 * n_elt_blockA;
  PDM_g_num_t *faceIniVtxA = malloc (sizeof(PDM_g_num_t) * s_faceIniVtxA);

  idx = 0;
  int idxFaceIni = 0;
  PDM_g_num_t _max_key = (2*n_g_newVtxA) / lComm + 1;
  int max_key = (int) _max_key;
  PDM_hash_tab_t *htSubEdgeA = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                                    &max_key);
  int s_tab_g = 0;

  double *u_inter = NULL;
  double *coords_inter  = NULL;
  PDM_g_num_t *nG_inter = NULL;
  int *order = NULL;
  PDM_g_num_t *nG_inter_tmp = NULL;
  double *coords_inter_tmp = NULL;

  for (int i = 0; i < n_elt_blockA; i++) {

    PDM_g_num_t gnum_boxA = block_gnumA[i];

    int lnum_boxA = blockA_lnum_data[i] - 1;

    /* Take into account sub-faces */

    for (int j = facesToSubFacesAIdx[i]; j < facesToSubFacesAIdx[i+1]; j++) {
      int nElt = subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
      int iBeg = subFacesConnecIdx[j];
      int iEnd = subFacesConnecIdx[j+1];
      for (int k = iBeg; k < iEnd; k++) {
        _sub_edge_t *se = malloc(sizeof(_sub_edge_t));
        int next = iBeg + ((k + 1 - iBeg) % nElt);

        int imin = (subFacesConnecA[k] < subFacesConnecA[next]) ?  k : next;
        int imax = (k == imin) ? next : k;

        se->vtx1 = subFacesConnecA[imin];
        se->vtx2 = subFacesConnecA[imax];

        for (int k1 = 0; k1 < 3; k1++) {
          se->coords1[k1] = subFacesCoordsA[3*imin+k1];
          se->coords2[k1] = subFacesCoordsA[3*imax+k1];
        }

        PDM_g_num_t _key = (se->vtx1 + se->vtx2) / lComm;
        int key = (int) _key;
				PDM_hash_tab_data_add (htSubEdgeA, &key, (void *) se);
      }
    }

    int n_vtxA = faceIdxCurrent[0][lnum_boxA + 1] - faceIdxCurrent[0][lnum_boxA];
    PDM_g_num_t *_faceToEdgeA = faceToEdgeCurrent[0] + faceIdxCurrent[0][lnum_boxA];
    PDM_g_num_t *_faceToVtxA  = faceToVtxCurrent[0] + faceIdxCurrent[0][lnum_boxA];
    double     *_face_vtxCooA = face_vtxCooCurrent[0] + 3 * faceIdxCurrent[0][lnum_boxA];

    /* Take into account origin edges */

    for (int j = 0; j < n_vtxA; j++) {

      int n_intersect = -1;
      PDM_edges_intersect_res_t **eirA = PDM_edges_intersect_get (intersect,
                                                                  PDM_EDGES_GET_FROM_A,
                                                                  PDM_ABS(_faceToEdgeA[j]),
                                                                  0,
                                                                  &n_intersect);

      int vtx_intersect = 0;
      for (int k = 0; k < n_intersect; k++) {

        PDM_line_intersect_t         _tIntersect;

        PDM_g_num_t                   _nGEdgeA;
        PDM_g_num_t                   _originEdgeA;
        PDM_g_num_t                   _endEdgeA;
        int                          _nNewPointsA;
        PDM_edges_intersect_point_t *_oNewPointsA;
        PDM_g_num_t                  *_linkA;
        PDM_g_num_t                  *_gNumVtxA;
        double                      *_coordsA;
        double                      *_uA;

        /*
         * Get intersection properties
         */

        PDM_edges_intersect_res_data_get (eirA[k],
                                          PDM_EDGES_INTERSECT_MESHA,
                                          &_nGEdgeA,
                                          &_originEdgeA,
                                          &_endEdgeA,
                                          &_tIntersect,
                                          &_nNewPointsA,
                                          &_oNewPointsA,
                                          &_linkA,
                                          &_gNumVtxA,
                                          &_coordsA,
                                          &_uA);

        for (int k1 = 0; k1 < _nNewPointsA; k1++) {
          if (_oNewPointsA[k1] != PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
            vtx_intersect += 1;
          }
        }
      }

      int s_tab = vtx_intersect + 2;
      if (s_tab > s_tab_g) {
        s_tab_g = s_tab;
        u_inter = realloc(u_inter, sizeof(double) * s_tab);
        coords_inter = realloc(coords_inter, sizeof(double) * 3 * s_tab);
        nG_inter = realloc(nG_inter, sizeof(PDM_g_num_t) * s_tab);

        order = realloc(order, sizeof(double) * s_tab);

        nG_inter_tmp = realloc(nG_inter_tmp, sizeof(PDM_g_num_t) * s_tab);
        coords_inter_tmp = realloc(coords_inter_tmp, sizeof(double) * 3 * s_tab);
      }

      int newSize = idxFaceIni + s_tab - 1;
      if (newSize > s_faceIniVtxA) {
        while (newSize > s_faceIniVtxA) {
          s_faceIniVtxA += PDM_MAX(1, s_faceIniVtxA);
        }
        faceIniVtxA = realloc (faceIniVtxA,
                               sizeof(PDM_g_num_t) * s_faceIniVtxA);
      }

      int next = (j+1) % n_vtxA;

      u_inter[0]       = 0;
      u_inter[s_tab-1] = 1.;
      nG_inter[0]      = _faceToVtxA[j];
      nG_inter[s_tab-1]  = _faceToVtxA[next];

      for (int k = 0; k < 3; k++) {
        coords_inter[0+k] = _face_vtxCooA[3*j+k];
        coords_inter[3*(s_tab-1)+k] = _face_vtxCooA[3*next+k];
      }

      vtx_intersect = 1;
      for (int k = 0; k < n_intersect; k++) {

        PDM_line_intersect_t         _tIntersect;

        PDM_g_num_t                   _nGEdgeA;
        PDM_g_num_t                   _originEdgeA;
        PDM_g_num_t                   _endEdgeA;
        int                          _nNewPointsA;
        PDM_edges_intersect_point_t *_oNewPointsA;
        PDM_g_num_t                  *_linkA;
        PDM_g_num_t                  *_gNumVtxA;
        double                      *_coordsA;
        double                      *_uA;

        /*
         * Get intersection properties
         */

        PDM_edges_intersect_res_data_get (eirA[k],
                                          PDM_EDGES_INTERSECT_MESHA,
                                          &_nGEdgeA,
                                          &_originEdgeA,
                                          &_endEdgeA,
                                          &_tIntersect,
                                          &_nNewPointsA,
                                          &_oNewPointsA,
                                          &_linkA,
                                          &_gNumVtxA,
                                          &_coordsA,
                                          &_uA);

        bool reverse = false;
        if (_originEdgeA != _faceToVtxA[j]) {
          assert (_faceToVtxA[next] == _originEdgeA);
          reverse = true;
        }

        for (int k1 = 0; k1 < _nNewPointsA; k1++) {
          if (_oNewPointsA[k1] != PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {

            u_inter[vtx_intersect]  = reverse ? (1 - _uA[k1]) : _uA[k1];
            nG_inter[vtx_intersect] = _gNumVtxA[k1];

            for (int k2 = 0; k2 < 3; k2++) {
              coords_inter[3*vtx_intersect + k2] = _coordsA[3*k1+k2];
            }
            vtx_intersect += 1;
          }
        }

      }

      for (int k = 0; k < s_tab; k++) {
        order[k] = k;
      }

      PDM_sort_double (u_inter, order, s_tab);

      for (int k = 0; k < s_tab; k++) {
        nG_inter_tmp[k] = nG_inter[order[k]];
        for (int k2 = 0; k2 < 3; k2++) {
          coords_inter_tmp[3*k+k2] = coords_inter[3*order[k] + k2];
        }
      }

      PDM_g_num_t *nG_inter_tmp2 = nG_inter;
      double *coords_inter_tmp2 = coords_inter;

      nG_inter = nG_inter_tmp;
      coords_inter = coords_inter_tmp;

      nG_inter_tmp = nG_inter_tmp2;
      coords_inter_tmp = coords_inter_tmp2;

      int k = 0;
      int idx1 = 0;
      while (s_tab > 0) {
        PDM_g_num_t val = nG_inter[k];
        nG_inter[idx1] = nG_inter[k];
        u_inter[idx1] = u_inter[k];
        for (int k2 = 0; k2 < 3; k2++) {
          coords_inter[3*idx1+k2] = coords_inter[3*k + k2];
        }
        idx1 += 1;
        while (nG_inter[k] == val) {
          k += 1;
          if (k >= s_tab) {
            break;
          }
        }
        if (k >= s_tab) {
          break;
        }
      }

      s_tab = idx1;

      for (int k1 = 0; k1 < s_tab; k1++) {
        int next1 = (k1 + 1) % s_tab;

        int imin = (nG_inter[k1] < nG_inter[next1]) ?  k1 : next1;
        int imax = (k1 == imin) ? next1 : k1;

        if (next1 != 0) {
          _sub_edge_t *se = malloc(sizeof(_sub_edge_t));
          se->vtx1 = nG_inter[imin];
          se->vtx2 = nG_inter[imax];

          faceIniVtxA[idxFaceIni++] = nG_inter[k1];

          for (int k2 = 0; k2 < 3; k2++) {
            se->coords1[k2] = coords_inter[3*imin+k2];
            se->coords2[k2] = coords_inter[3*imax+k2];
          }

          PDM_g_num_t _key = (se->vtx1 + se->vtx2) / lComm;
          int key = (int) _key;
          PDM_hash_tab_data_add (htSubEdgeA, &key, (void *) se);
        }
      }

      free (eirA);

    }

    /* Count the number of reference of each sub-edge */

    int n_used_keys = PDM_hash_tab_n_used_keys_get (htSubEdgeA);
    PDM_g_num_t *used_keys = PDM_hash_tab_used_keys_get (htSubEdgeA);

    int t_n_data = 0;
    for (int j = 0; j < n_used_keys; j++) {
      int _key = used_keys[j];
      int n_data = PDM_hash_tab_n_data_get (htSubEdgeA, &_key);
      t_n_data += n_data;
    }

    PDM_g_num_t *oneRef = malloc(sizeof(PDM_g_num_t) * 2 * t_n_data);
    double     *coordsOneRef = malloc(sizeof(double) * 6 * t_n_data);
    int nOneRef = 0;

    for (int j = 0; j < n_used_keys; j++) {
      int _key = used_keys[j];
      int n_data = PDM_hash_tab_n_data_get (htSubEdgeA, &_key);
      int tag[n_data];
      PDM_array_reset_int(tag, n_data, 0);

      _sub_edge_t **se = (_sub_edge_t **) PDM_hash_tab_data_get (htSubEdgeA, &_key);
      for (int k1 = 0; k1 < n_data; k1++) {
        if (tag[k1] == 1) continue;
        _sub_edge_t *_se1 = se[k1];
        for (int k2 = k1+1; k2 < n_data; k2++) {
          if (tag[k2] == 1) continue;
          _sub_edge_t *_se2 = se[k2];
          if (((_se1->vtx1 == _se2->vtx1) && (_se1->vtx2 == _se2->vtx2)) ||
              ((_se1->vtx1 == _se2->vtx2) && (_se1->vtx2 == _se2->vtx1))) {
            tag[k1] = 1;
            tag[k2] = 1;
          }
        }
      }
      for (int k1 = 0; k1 < n_data; k1++) {
        if (tag[k1] == 0) {
          _sub_edge_t *_se1 = se[k1];
          oneRef[2 * nOneRef    ] = _se1->vtx1;
          oneRef[2 * nOneRef + 1] = _se1->vtx2;
          for (int k2 = 0; k2 < 3; k2++) {
            coordsOneRef[6 * nOneRef + k2] = _se1->coords1[k2];
            coordsOneRef[6 * nOneRef + 3 + k2] = _se1->coords2[k2];
          }
          nOneRef += 1;
        }
      }
    }

    PDM_g_num_t *tag = malloc(sizeof(PDM_g_num_t) * nOneRef);
    int nAddSubFace = 0;
    int *addSubFaceIdx = malloc(sizeof(int) * (nOneRef + 1));
    PDM_g_num_t *addSubFace = malloc(sizeof(PDM_g_num_t) * nOneRef);
    double *addSubFaceCoords = malloc(sizeof(double) * 3 * nOneRef);

    addSubFaceIdx[0] = 0;
    for (int k1 = 0; k1 < nOneRef; k1++) {
      tag[k1] = 0;
      addSubFaceIdx[k1+1] = 0;
    }

    /* Add additional sub-faces */

    PDM_g_num_t iniVal = -1;
    PDM_g_num_t nextVal = -1;
    idx = 0;
    int idx2 = 0;

    for (int k1 = 0; k1 < nOneRef; k1++) {
      if (tag[k1] == 1) continue;

      iniVal = oneRef[2*k1];
      nextVal = oneRef[2*k1 + 1];

      addSubFace[idx++] = iniVal;
      addSubFace[idx++] = oneRef[2*k1 + 1];

      for (int k2 = 0; k2 < 3; k2++) {
        addSubFaceCoords[idx2++] = coordsOneRef[6*k1+k2];
      }
      for (int k2 = 0; k2 < 3; k2++) {
        addSubFaceCoords[idx2++] = coordsOneRef[6*k1+3+k2];
      }

      nAddSubFace += 1;
      addSubFaceIdx[nAddSubFace]++;
      addSubFaceIdx[nAddSubFace]++;

      int k_prev = k1;
      while (iniVal != nextVal) {
        for (int k2 = 0; k2 < nOneRef; k2++) {
          if ((tag[k2] == 1) || (k_prev == k2)) continue;
          if (nextVal == oneRef[2*k2]) {
            nextVal = oneRef[2*k2 + 1];
            if (nextVal != iniVal) {
              addSubFace[idx++] = nextVal;
              for (int k3 = 0; k3 < 3; k3++) {
                addSubFaceCoords[idx2++] = coordsOneRef[6*k2+3+k3];
              }
              addSubFaceIdx[nAddSubFace]++;
            }
            tag[k2] = 1;
            k_prev = k2;
            break;
          }
          else if (nextVal == oneRef[2*k2 + 1]) {
            nextVal = oneRef[2*k2];
            if (nextVal != iniVal) {
              addSubFace[idx++] = nextVal;
              for (int k3 = 0; k3 < 3; k3++) {
                addSubFaceCoords[idx2++] = coordsOneRef[6*k2+k3];
              }
              addSubFaceIdx[nAddSubFace]++;
            }
            tag[k2] = 1;
            k_prev = k2;
            break;
          }
        }
      }
    }

    free (oneRef);
    free (coordsOneRef);

    for (int k1 = 0; k1 < nAddSubFace; k1++) {
      addSubFaceIdx[k1+1] += addSubFaceIdx[k1];
    }

    int newSize = iclipp + nAddSubFace + 1;
    if (newSize > s_subFacesConnecIdx) {
      while (newSize > s_subFacesConnecIdx) {
        s_subFacesConnecIdx += PDM_MAX (1, s_subFacesConnecIdx/3);
      }
      subFacesConnecIdx = realloc(subFacesConnecIdx,
                                sizeof(int) * s_subFacesConnecIdx);
    }

    int ibeg = subFacesConnecIdx[iclipp];
    newSize = ibeg + addSubFaceIdx[nAddSubFace];

    if (newSize > s_subFacesConnecA) {
      while (newSize > s_subFacesConnecA) {
        s_subFacesConnecA += PDM_MAX (1, s_subFacesConnecA/3);
      }
      subFacesConnecA = realloc(subFacesConnecA,
                                sizeof(PDM_g_num_t) * s_subFacesConnecA);
    }

    int ibeg2 = 3 * ibeg;
    newSize = 3 * (ibeg + addSubFaceIdx[nAddSubFace]);

    if (newSize > s_subFacesCoordsA) {
      while (newSize > s_subFacesCoordsA) {
        s_subFacesCoordsA += PDM_MAX (1, s_subFacesCoordsA/3);
      }
      subFacesCoordsA = realloc (subFacesCoordsA,
                                sizeof(double) * s_subFacesCoordsA);
    }

    newSize = 5*(iclipp + nAddSubFace);
    if (newSize > s_subFacesToFaces) {
      while (newSize > s_subFacesToFaces) {
        s_subFacesToFaces += PDM_MAX (1, s_subFacesToFaces/3);
      }
      subFacesToFaces = realloc(subFacesToFaces,
                                sizeof(PDM_g_num_t) * s_subFacesToFaces);
    }

    for (int k1 = 0; k1 < addSubFaceIdx[nAddSubFace]; k1++) {
      subFacesConnecA[ibeg++] = addSubFace[k1];
    }

    for (int k1 = 0; k1 < 3*addSubFaceIdx[nAddSubFace]; k1++) {
      subFacesCoordsA[ibeg2++] = addSubFaceCoords[k1];
    }
    free (addSubFaceCoords);

    for (int k1 = 0; k1 < nAddSubFace; k1++) {
      subFacesConnecIdx[iclipp+1] = subFacesConnecIdx[iclipp] +
                                    (addSubFaceIdx[k1+1] - addSubFaceIdx[k1]);
      subFacesToFaces[5*iclipp  ] = i;          // A local number
      subFacesToFaces[5*iclipp+1] = gnum_boxA;  // A global number
      subFacesToFaces[5*iclipp+2] = -1;         // B local number
      subFacesToFaces[5*iclipp+3] = -1;         // B global number
      subFacesToFaces[5*iclipp+4] = -1;         // B global number

      iclipp += 1 ;
    }

    //iclipp += nAddSubFace;

    facesToAddSubFacesAIdx[i+1] = facesToAddSubFacesAIdx[i] + nAddSubFace;

    /* Free */

    free (tag);
    free (addSubFaceIdx);
    free (addSubFace);

    faceIniVtxIdxA[i+1] = idxFaceIni;

    PDM_hash_tab_purge (htSubEdgeA, PDM_TRUE);

  }

  free (u_inter);
  free (coords_inter);
  free (nG_inter);
  free (coords_inter_tmp);
  free (nG_inter_tmp);
  free (order);

  PDM_hash_tab_purge (htSubEdgeA, PDM_TRUE);
  PDM_hash_tab_free (htSubEdgeA);

  const int nSubFaces = iclipp;

  /*
   * Update Memory
   */

  faceIniVtxA = realloc (faceIniVtxA,
                          sizeof(PDM_g_num_t) * faceIniVtxIdxA[n_elt_blockA]);

  subFacesConnecA = realloc (subFacesConnecA,
                             sizeof(PDM_g_num_t) * subFacesConnecIdx[nSubFaces]);

  subFacesCoordsA = realloc (subFacesCoordsA,
                             3 * sizeof(double) * subFacesConnecIdx[nSubFaces]);

  subFacesConnecB = realloc (subFacesConnecB,
                             sizeof(PDM_g_num_t) * subFacesConnecIdx[nSharedSubFaces]);

  subFacesConnecIdx = realloc (subFacesConnecIdx,
                               sizeof(int) * (nSubFaces + 1));

  subFacesToFaces = realloc (subFacesToFaces,
                             sizeof(PDM_g_num_t) *  5 * nSubFaces);


  if (1 == 0) {
    printf("\n-- Liste des sous-facettes des elements de A (intersection + complement intersection)\n");
    for (int i = 0; i < n_elt_blockA; i++) {

      PDM_g_num_t gnum_boxA = block_gnumA[i];

      //int lnum_boxA = blockA_lnum_data[i] - 1;

      printf ("face origin A : "PDM_FMT_G_NUM"\n", gnum_boxA);

      for (int j = facesToSubFacesAIdx[i]; j < facesToSubFacesAIdx[i+1]; j++) {
        printf (" - sub face %d : ", j);
        for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
          printf(" "PDM_FMT_G_NUM, subFacesConnecA[k]);
        }
        printf("\n");
      }

      for (int j = facesToAddSubFacesAIdx[i]; j < facesToAddSubFacesAIdx[i+1]; j++) {
        printf (" - add sub face %d : ", j);
        for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
          printf(" "PDM_FMT_G_NUM, subFacesConnecA[k]);
          for (int k1 = 0; k1 < 3; k1++) {
            printf(" %12.5e", subFacesCoordsA[3*k+k1]);
          }
          printf("/");
        }
        printf("\n");
        /* for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) { */
        /*   for (int k1 = 0; k1 < 3; k1++) { */
        /*     printf(" %12.5e", subFacesCoordsA[3*k+k1]); */
        /*   } */
        /*   printf("/"); */
        /* } */
        /* printf("\n"); */
      }
    }

    printf("\n-- Faces d'origine de A avec nouveaux points\n");

    for (int i = 0; i < n_elt_blockA; i++) {
      PDM_g_num_t gnum_boxA = block_gnumA[i];

      printf ("face origin A : "PDM_FMT_G_NUM"\n", gnum_boxA);

      for (int j = faceIniVtxIdxA[i]; j < faceIniVtxIdxA[i+1]; j++) {
        printf(" "PDM_FMT_G_NUM, faceIniVtxA[j]);
      }
      printf("\n");
    }
  }

  /*
   * Global numbering definition of sub-faces
   * This numbering is shared between A and B
   * Partial split faces are stored behind sub-faces
   * Not split faces are stored behind partial split faces
   *
   */

  PDM_g_num_t nAddSubFaces_l = nSubFaces - nSharedSubFaces;
  PDM_g_num_t beg_nAddSubFaces;
  PDM_g_num_t end_subFaces;

  PDM_g_num_t _max_gnum_loc = -1;
  if (nSharedSubFaces > 0) {
    _max_gnum_loc = gNumSubFacesA[nSharedSubFaces-1];
  }
  PDM_MPI_Allreduce(&_max_gnum_loc, &end_subFaces, 1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ol->comm);

  PDM_MPI_Scan (&nAddSubFaces_l, &beg_nAddSubFaces,
            1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  beg_nAddSubFaces += -nAddSubFaces_l;


  gNumSubFacesA = realloc(gNumSubFacesA, sizeof(PDM_g_num_t) * nSubFaces);
  for (int i = nSharedSubFaces; i < nSubFaces; i++) {
    gNumSubFacesA[i] = end_subFaces + beg_nAddSubFaces + i - nSharedSubFaces + 1;
  }

  _max_gnum_loc = -1;
  if (nSharedSubFaces > 0) {
    _max_gnum_loc = gNumSubFacesA[nSharedSubFaces-1];
  }

  PDM_g_num_t nTSubFacesA;
  PDM_MPI_Allreduce(&_max_gnum_loc, &nTSubFacesA, 1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ol->comm);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESA] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESA]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_COMPUTE_ADD_SUB_FACESA]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_COMPUTE_ADD_SUB_FACESA]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_compute_add_sub_facesa !!!!\n");
    fflush(stdout);
  }

  if (1 == 0) {
    printf("\n-- Numero global des sous-facettes de A :");
    for (int i = 0; i < nSubFaces; i++) {
      printf(" "PDM_FMT_G_NUM, gNumSubFacesA[i]);
    }
    printf("\n");
  }


  /* exit(0); */

  /*****************************************************************************
   *                                                                           *
   * Compute sub-faces coming from B faces with partial A cover.               *
   *                                                                           *
   * Check that each sub-edge is referenced twice                              *
   *    - Take into sub-edges of sub-faces                                     *
   *    - Take into sub-edges of origin elements                               *
   *                                                                           *
   * If sub-edges are not twice referenced, the origin face is partial covered *
   * by A. Additional sub-faces have to compute.                               *
   *                                                                           *
   * Two steps :                                                               *
   *    - Redistribute B cells (cs_part_to_block) to merge results             *
   *    - Merge Data                                                           *
   *                                                                           *
   * TODO: Deplacer cette partie dans une fonction permettant d'enchainer les  *
   *  communications en mode non-bloquant en plusieurs temps pendant le calcul *
   *  des faces complementaires de A                                           *
   *                                                                           *
   * TODO: Commencer a faire le transfert du resultat du maillage A pendant    *
   * le calcul des sous-facettes complementaires de B                          *
   *                                                                           *
   *****************************************************************************/


  /* Send Number of vtx and global Number for each sub-face */


  PDM_g_num_t *blockA_boxesB_gnum_data_cp = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * blockA_boxesB_idx[n_elt_blockA]);
  int *blockB_lnum_data = (int *) malloc (sizeof(int) * blockA_boxesB_idx[n_elt_blockA]);

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    blockB_lnum_data[i] = i;
    blockA_boxesB_gnum_data_cp[i] = blockA_boxesB_gnum_data[i];
  }

  PDM_sort_long (blockA_boxesB_gnum_data_cp, blockB_lnum_data, blockA_boxesB_idx[n_elt_blockA]);

  int *idx_without_dupl = malloc(sizeof(int) * blockA_boxesB_idx[n_elt_blockA]);

  int _k1 = 0;
  int n_boxesB_without_dupl = 0;
  while (blockA_boxesB_idx[n_elt_blockA] > 0) {
    PDM_g_num_t val =  blockA_boxesB_gnum_data_cp[_k1];
    blockA_boxesB_gnum_data_cp[n_boxesB_without_dupl] =  blockA_boxesB_gnum_data_cp[_k1];
    while ( blockA_boxesB_gnum_data_cp[_k1] == val) {
      idx_without_dupl[_k1] = n_boxesB_without_dupl;
      _k1 += 1;
      if (_k1 >= blockA_boxesB_idx[n_elt_blockA]) {
        break;
      }
    }
    n_boxesB_without_dupl += 1;
    if (_k1 >= blockA_boxesB_idx[n_elt_blockA]) {
      break;
    }
  }

  blockA_boxesB_gnum_data_cp = realloc(blockA_boxesB_gnum_data_cp, sizeof(PDM_g_num_t) * n_boxesB_without_dupl);

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    facesToSubFacesBIdx[i+1] += facesToSubFacesBIdx[i];
  }

  double *weight = malloc (sizeof(double) * n_boxesB_without_dupl);
  int *_tmp_stride = malloc (sizeof(int) * n_boxesB_without_dupl);

  for (int i = 0; i < n_boxesB_without_dupl; i++) {
    weight[i] = 0;
    _tmp_stride[i] = 0;
  }

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int _nEltSubFaces = (facesToSubFacesBIdx[blockB_lnum_data[i]+1] - facesToSubFacesBIdx[blockB_lnum_data[i]]);
    _tmp_stride[idx_without_dupl[i]] += _nEltSubFaces;
    weight[idx_without_dupl[i]] += _nEltSubFaces;
  }

  int *_tmp_idx = malloc (sizeof(int) * (n_boxesB_without_dupl + 1));
  _tmp_idx[0] = 0;

  for (int i = 0; i < n_boxesB_without_dupl; i++) {
    _tmp_idx[i+1] = _tmp_idx[i] + _tmp_stride[i];
    _tmp_stride[i] = 0;
  }

  PDM_part_to_block_t *ptb_boxesB = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_MERGE,
                                                              1.,
                                                              (PDM_g_num_t **) &blockA_boxesB_gnum_data_cp,
                                                              &weight,
                                                              &n_boxesB_without_dupl,
                                                              1,
                                                              ol->comm);
  free (weight);

  int *recv_stride1 = NULL;

  PDM_g_num_t *subFacesConnecN =
          (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 4* nSharedSubFaces);

  originA = PDM_box_set_origin_get (boxesA);

  originB = PDM_box_set_origin_get (boxesB);

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int _idx = blockB_lnum_data[i];
    for (int j = facesToSubFacesBIdx[_idx]; j < facesToSubFacesBIdx[_idx+1]; j++) {
      int idx_subFacesConnecN = 4 * _tmp_idx[idx_without_dupl[i]] + _tmp_stride[idx_without_dupl[i]];
      int lParentFace  = (int) subFacesToFaces[5*j];
      PDM_g_num_t iProc = originA[3*(blockA_lnum_data[lParentFace] - 1)];
      PDM_g_num_t i_part = originA[3*(blockA_lnum_data[lParentFace] - 1)+1];
      subFacesConnecN[idx_subFacesConnecN++] = subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
      subFacesConnecN[idx_subFacesConnecN++] = iProc;
      subFacesConnecN[idx_subFacesConnecN++] = i_part;
      subFacesConnecN[idx_subFacesConnecN++] = gNumSubFacesA[j];
      _tmp_stride[idx_without_dupl[i]]+=4;
    }
  }

  PDM_g_num_t *recvSubFacesConnecBN;
  PDM_part_to_block_exch (ptb_boxesB,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &_tmp_stride,
                         (void **) &subFacesConnecN,
                         &recv_stride1,
                         (void **) &recvSubFacesConnecBN);

  free (subFacesConnecN);

  /* Send Connectivity for each sub-face */

  int *recv_stride2 = NULL;

  PDM_array_reset_int(_tmp_stride, n_boxesB_without_dupl, 0);

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int _idx = blockB_lnum_data[i];
    for (int j = facesToSubFacesBIdx[_idx]; j < facesToSubFacesBIdx[_idx+1]; j++) {
      _tmp_stride[idx_without_dupl[i]] += subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
    }
  }

  for (int i = 0; i < n_boxesB_without_dupl; i++) {
    _tmp_idx[i+1] = _tmp_idx[i] + _tmp_stride[i];
    _tmp_stride[i] = 0;
  }

  PDM_g_num_t *subFacesConnecB_ordered = malloc(sizeof(PDM_g_num_t) * subFacesConnecIdx[nSharedSubFaces]);
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int _idx = blockB_lnum_data[i];
    for (int j = facesToSubFacesBIdx[_idx]; j < facesToSubFacesBIdx[_idx+1]; j++) {
      int idx_subFacesConnecN = _tmp_idx[idx_without_dupl[i]] + _tmp_stride[idx_without_dupl[i]];
      for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
        subFacesConnecB_ordered[idx_subFacesConnecN++] = subFacesConnecB[k];
        _tmp_stride[idx_without_dupl[i]]++;
      }
    }
  }

  free (subFacesConnecB);
  PDM_g_num_t *recvSubFacesConnecB;
  PDM_part_to_block_exch (ptb_boxesB,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &_tmp_stride,
                         (void **) &subFacesConnecB_ordered,
                         &recv_stride2,
                         (void **) &recvSubFacesConnecB);

  free (subFacesConnecB_ordered);

  PDM_array_reset_int(_tmp_stride, n_boxesB_without_dupl, 0);

  double *subFacesCoordsA_ordered = malloc(3 * sizeof(double) * subFacesConnecIdx[nSharedSubFaces]);
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int _idx = blockB_lnum_data[i];
    for (int j = facesToSubFacesBIdx[_idx]; j < facesToSubFacesBIdx[_idx+1]; j++) {
      int idx_subFacesConnecN = 3 * _tmp_idx[idx_without_dupl[i]] + _tmp_stride[idx_without_dupl[i]];
      for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
        for (int k11 = 0; k11 < 3; k11++) {
          subFacesCoordsA_ordered[idx_subFacesConnecN++] = subFacesCoordsA[3*k+k11];
          _tmp_stride[idx_without_dupl[i]]++;
        }
      }
    }
  }

  int *recv_stride21;
  double *recvSubFacesCoordsB;
  PDM_part_to_block_exch (ptb_boxesB,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &_tmp_stride,
                         (void **) &subFacesCoordsA_ordered,
                         &recv_stride21,
                         (void **) &recvSubFacesCoordsB);

  free (subFacesCoordsA_ordered);

  /* Send result of intersection */

  /* Send edge and vertex connectivity */


  int idx1 = 0;

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    int _idx = blockB_lnum_data[i];
    int lnum_boxB = blockA_boxesB_lnum_data[_idx] - 1;
    _tmp_stride[idx_without_dupl[i]] = (faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB]);
  }

  for (int i = 0; i < n_boxesB_without_dupl; i++) {
    _tmp_idx[i+1] = _tmp_idx[i] + _tmp_stride[i];
    _tmp_stride[i] = 0;
  }

  /*
   * Send :
   *   - number of intersection points per edge per element
   *   - edge vertices
   *   - edge vertices coordinates (in an other exchange)
   */

  int *recv_stride4 = NULL;

  PDM_g_num_t *sendFaceToEdgeNPtInt = PDM_array_const_gnum(3*_tmp_idx[n_boxesB_without_dupl], 0);
  idx1 = 0;
  int s_properties = 0;

  int *_tmp_stride2 = PDM_array_zeros_int(n_boxesB_without_dupl);
  int *_tmp_idx2 = malloc (sizeof(int) * (n_boxesB_without_dupl + 1));

  _tmp_idx2[0] = 0;

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    if (_tmp_stride[idx_without_dupl[i]] == 0) {
      int _idx = blockB_lnum_data[i];
      int lnum_boxB = blockA_boxesB_lnum_data[_idx] - 1;
      int iBeg = faceIdxCurrent[1][lnum_boxB];
      int nEdge = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
      int _idx1 = 3 *_tmp_idx[idx_without_dupl[i]];
      assert (faceIdxCurrent[1][lnum_boxB] < faceIdxCurrent[1][lnum_boxB + 1]);
      assert (sendFaceToEdgeNPtInt[_idx1] == 0);
      for (int j = faceIdxCurrent[1][lnum_boxB]; j < faceIdxCurrent[1][lnum_boxB + 1]; j++) {
        int n_intersect = -1;

        int next = iBeg + ((j + 1 - iBeg) % nEdge);

        PDM_g_num_t vtx1 = faceToVtxCurrent[1][j];
        PDM_g_num_t vtx2 = faceToVtxCurrent[1][next];

        PDM_edges_intersect_res_t **eirB = PDM_edges_intersect_get (intersect,
                                                                    PDM_EDGES_GET_FROM_B,
                                                                    0,
                                                                    PDM_ABS(faceToEdgeCurrent[1][j]),
                                                                    &n_intersect);
        for (int k = 0; k < n_intersect; k++) {

          PDM_line_intersect_t         _tIntersect;

          PDM_g_num_t                   _nGEdgeB;
          PDM_g_num_t                   _originEdgeB;
          PDM_g_num_t                   _endEdgeB;
          int                          _nNewPointsB;
          PDM_edges_intersect_point_t *_oNewPointsB;
          PDM_g_num_t                  *_linkB;
          PDM_g_num_t                  *_gNumVtxB;
          double                      *_coordsB;
          double                      *_uB;

          PDM_edges_intersect_res_data_get (eirB[k],
                                            PDM_EDGES_INTERSECT_MESHB,
                                            &_nGEdgeB,
                                            &_originEdgeB,
                                            &_endEdgeB,
                                            &_tIntersect,
                                            &_nNewPointsB,
                                            &_oNewPointsB,
                                            &_linkB,
                                            &_gNumVtxB,
                                            &_coordsB,
                                            &_uB);

          sendFaceToEdgeNPtInt[_idx1] += _nNewPointsB;
        }
        free (eirB);
        _tmp_stride2[idx_without_dupl[i]] += (int) sendFaceToEdgeNPtInt[_idx1];
        s_properties += (int) sendFaceToEdgeNPtInt[_idx1++];
        sendFaceToEdgeNPtInt[_idx1++] = vtx1;
        sendFaceToEdgeNPtInt[_idx1++] = vtx2;
        _tmp_stride[idx_without_dupl[i]] += 3;
      }
    }
  }

  for (int i = 0; i < n_boxesB_without_dupl; i++) {
    _tmp_idx2[i+1] = _tmp_idx2[i] + _tmp_stride2[i];
    _tmp_stride2[i] = 0;
  }

  PDM_g_num_t *recvFaceToEdgeNPtInt;
  PDM_part_to_block_exch (ptb_boxesB,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &_tmp_stride,
                         (void **) &sendFaceToEdgeNPtInt,
                         &recv_stride4,
                         (void **) &recvFaceToEdgeNPtInt);

  free (sendFaceToEdgeNPtInt);

  PDM_array_reset_int(_tmp_stride, n_boxesB_without_dupl, 0);

  int *recv_stride41 = NULL;

   double *sendFaceToEdgeCoordsvtx = malloc (sizeof(double) * 6 * _tmp_idx[n_boxesB_without_dupl]);
  //double *recvFaceToEdgeCoordsvtx = malloc (sizeof(double) * idx1);
  double *recvFaceToEdgeCoordsvtx;

  idx1 = 0;
  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {
    if (_tmp_stride[idx_without_dupl[i]] == 0) {
      int _idx = blockB_lnum_data[i];
      int lnum_boxB = blockA_boxesB_lnum_data[_idx]-1;
      int iBeg = faceIdxCurrent[1][lnum_boxB];
      int nEdge = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
      int _idx1 = 6 * _tmp_idx[idx_without_dupl[i]];
      for (int j = faceIdxCurrent[1][lnum_boxB]; j < faceIdxCurrent[1][lnum_boxB + 1]; j++) {

        int next = iBeg + ((j + 1 - iBeg) % nEdge);

        for (int k11 = 0; k11 < 3; k11++) {
          sendFaceToEdgeCoordsvtx[_idx1++] = face_vtxCooCurrent[1][3*j+k11];
        }

        for (int k11 = 0; k11 < 3; k11++) {
          sendFaceToEdgeCoordsvtx[_idx1++] = face_vtxCooCurrent[1][3*next+k11];
        }
        _tmp_stride[idx_without_dupl[i]] += 6;
      }
    }
  }

  PDM_part_to_block_exch (ptb_boxesB,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &_tmp_stride,
                         (void **) &sendFaceToEdgeCoordsvtx,
                         &recv_stride41,
                         (void **) &recvFaceToEdgeCoordsvtx);

  free (sendFaceToEdgeCoordsvtx);

  /* Send origin, global number and U for each intersection point */

  int *recv_stride5 = NULL;
  int *recv_stride6 = NULL;

  for (int i = 0; i < n_boxesB_without_dupl; i++) {
    _tmp_stride[i] = 0;
    _tmp_stride2[i] = 0;
  }

  PDM_g_num_t *sendFaceToEdgeOrAndGnumPtInt = malloc (sizeof(PDM_g_num_t) * 3 * s_properties);

  double *sendFaceToEdgeUPtInt = malloc (sizeof(double) * 4 * s_properties);

  idx1 = 0;
  int idx2 = 0;
  int idx3 = 0;

  for (int i = 0; i < blockA_boxesB_idx[n_elt_blockA]; i++) {

    if (_tmp_stride[idx_without_dupl[i]] == 0) {
      int _idx = blockB_lnum_data[i];
      int lnum_boxB = blockA_boxesB_lnum_data[_idx] - 1;
      int _idx2 = 3 * _tmp_idx2[idx_without_dupl[i]];
      int _idx3 = 4 * _tmp_idx2[idx_without_dupl[i]];
      //int iBeg = faceIdxCurrent[1][lnum_boxB];
      //int nEdge = faceIdxCurrent[1][lnum_boxB + 1] - faceIdxCurrent[1][lnum_boxB];
      for (int j = faceIdxCurrent[1][lnum_boxB]; j < faceIdxCurrent[1][lnum_boxB + 1]; j++) {

        //int next = iBeg + ((j + 1 - iBeg) % nEdge);
        int n_intersect = -1;
        PDM_edges_intersect_res_t **eirB = PDM_edges_intersect_get (intersect,
                                                                    PDM_EDGES_GET_FROM_B,
                                                                    0,
                                                                    PDM_ABS(faceToEdgeCurrent[1][j]),
                                                                    &n_intersect);

        for (int k = 0; k < n_intersect; k++) {
          PDM_line_intersect_t         _tIntersect;

          PDM_g_num_t                   _nGEdgeB;
          PDM_g_num_t                   _originEdgeB;
          PDM_g_num_t                   _endEdgeB;
          int                          _nNewPointsB;
          PDM_edges_intersect_point_t *_oNewPointsB;
          PDM_g_num_t                  *_linkB;
          PDM_g_num_t                  *_gNumVtxB;
          double                      *_coordsB;
          double                      *_uB;

          PDM_edges_intersect_res_data_get (eirB[k],
                                            PDM_EDGES_INTERSECT_MESHB,
                                            &_nGEdgeB,
                                            &_originEdgeB,
                                            &_endEdgeB,
                                            &_tIntersect,
                                            &_nNewPointsB,
                                            &_oNewPointsB,
                                            &_linkB,
                                            &_gNumVtxB,
                                            &_coordsB,
                                            &_uB);

          for (int k11 = 0; k11 < _nNewPointsB; k11++) {
            _tmp_stride2[idx_without_dupl[i]] += 3;
            sendFaceToEdgeOrAndGnumPtInt[_idx2++] = _originEdgeB;
            sendFaceToEdgeOrAndGnumPtInt[_idx2++] = _gNumVtxB[k11];
            sendFaceToEdgeOrAndGnumPtInt[_idx2++] = _oNewPointsB[k11];
            sendFaceToEdgeUPtInt[_idx3++]         = _uB[k11];
            sendFaceToEdgeUPtInt[_idx3++]         = _coordsB[3*k11  ];
            sendFaceToEdgeUPtInt[_idx3++]         = _coordsB[3*k11+1];
            sendFaceToEdgeUPtInt[_idx3++]         = _coordsB[3*k11+2];
          }
        }
        free (eirB);
      }
      _tmp_stride[idx_without_dupl[i]]++;
    }
  }

  free (idx_without_dupl);

  PDM_edges_intersect_free (intersect);

  free (faceIdxCurrent[0]);
  free (faceIdxCurrent[1]);

  // int s_prop2 = 0;
  // for (int i = 0; i < n_boxesB_without_dupl; i++) {
  //   s_prop2 += _tmp_stride2[i];
  // }

  PDM_g_num_t *recvFaceToEdgeOrAndGnumPtInt;
  PDM_part_to_block_exch (ptb_boxesB,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &_tmp_stride2,
                         (void **) &sendFaceToEdgeOrAndGnumPtInt,
                         &recv_stride5,
                         (void **) &recvFaceToEdgeOrAndGnumPtInt);

  free (sendFaceToEdgeOrAndGnumPtInt);
  free (recv_stride5);

  for (int i = 0; i < n_boxesB_without_dupl; i++) {
    _tmp_stride2[i] = (_tmp_stride2[i]/3) * 4;
  }

  double *recvFaceToEdgeUPtInt;
  PDM_part_to_block_exch (ptb_boxesB,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &_tmp_stride2,
                         (void **) &sendFaceToEdgeUPtInt,
                         &recv_stride6,
                         (void **) &recvFaceToEdgeUPtInt);

  free (_tmp_stride);
  free (_tmp_stride2);
  free (_tmp_idx);
  free (_tmp_idx2);

  free (sendFaceToEdgeUPtInt);
  free (recv_stride6);

  /* Cleanup */

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_DISTRIB_BOXESB_BLOCK] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_DISTRIB_BOXESB_BLOCK]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_DISTRIB_BOXESB_BLOCK]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_DISTRIB_BOXESB_BLOCK]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_distrib_boxesb_block !!!!\n");
    fflush(stdout);
  }

  /* printf ("!!!! Arret temporaire apres distribution des faces de B en bloc !!!!\n"); */
  /* if (i_rank == 0) { */
  /*   _ol_dump_times (ol); */
  /* } */
  /* exit(0); */
  /* printf ("!!!! Arret temporaire apres distribution des faces de B en bloc !!!!\n"); */

  /*
   * Loop on faces
   */

  int n_elt_blockB = PDM_part_to_block_n_elt_block_get (ptb_boxesB);

  PDM_g_num_t *block_gnumB = PDM_part_to_block_block_gnum_get (ptb_boxesB);

  int nSubFacesB = 0;

  s_subFacesConnecB = 4 * n_elt_blockB;
  subFacesConnecB = malloc (sizeof(PDM_g_num_t) * s_subFacesConnecB);

  int s_subFacesCoordsB = 3 * 4 * n_elt_blockB;
  double *subFacesCoordsB = malloc (sizeof(double) * s_subFacesCoordsB);

  int s_subFacesConnecIdxB = 4 * n_elt_blockB;
  int *subFacesConnecIdxB = malloc (sizeof(int) * s_subFacesConnecIdxB);
  subFacesConnecIdxB[0] = 0;

  int s_subFacesToFaceB = 4 * n_elt_blockB;
  //FIXME : Le tableau subFacesToFaceB est probablement inutile (A verifier)
  PDM_g_num_t *subFacesToFaceB = malloc (sizeof(PDM_g_num_t) * s_subFacesToFaceB);
  PDM_g_num_t *gNumSubFacesB   = malloc (sizeof(PDM_g_num_t) * s_subFacesToFaceB);

  int s_subFacesToLinkA = 4 * n_elt_blockB;
  PDM_g_num_t *subFacesToLinkA = malloc (sizeof(PDM_g_num_t) * s_subFacesToLinkA);

      idx1 = 0;
      idx2 = 0;
      idx3 = 0;
  int idx4 = 0;
  int idx41 = 0;
  int idx5 = 0;
  int idx6 = 0;

  int sIntEdge = 4;
  PDM_g_num_t *vtxIntEdgeSorted   = malloc (sizeof(PDM_g_num_t) * sIntEdge);
  double     *coordsIntEdgeSorted = malloc (sizeof(double) * 3 * sIntEdge);
  order                           = malloc (sizeof(double) * sIntEdge);

  facesToSubFacesBIdx = realloc (facesToSubFacesBIdx, sizeof(int) * (n_elt_blockB + 1));
  facesToSubFacesBIdx[0] = 0;

  PDM_g_num_t n_t_nAddSubFace = 0;

  int *faceIniVtxIdxB = malloc (sizeof(int) * (1 + n_elt_blockB));
  faceIniVtxIdxB[0] = 0;

  int s_faceIniVtxB = 4 * n_elt_blockB;
  PDM_g_num_t *faceIniVtxB = malloc (sizeof(PDM_g_num_t) * s_faceIniVtxB);

  idx = 0;
  _max_key = (2*n_g_newVtxB) / lComm + 1;
  max_key = (int) _max_key;
  PDM_hash_tab_t *htSubEdgeB = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                                    &max_key);
  PDM_hash_tab_t *htEdgeB = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                                    &max_key);

  PDM_g_num_t *sum_vtx = NULL;
  int s_sum_vtx = 0;

  for (int i = 0; i < n_elt_blockB; i++) {

    faceIniVtxIdxB[i+1] = faceIniVtxIdxB[i];

    PDM_hash_tab_purge (htSubEdgeB, PDM_TRUE);

    int n_used_keys = PDM_hash_tab_n_used_keys_get (htEdgeB);
    PDM_g_num_t *used_keys = PDM_hash_tab_used_keys_get (htEdgeB);

    for (int j = 0; j < n_used_keys; j++) {
      int _key = (int) used_keys[j];
      int n_data = PDM_hash_tab_n_data_get (htEdgeB, &_key);
      _sub_vertices_origin_edge_t **data =
        (_sub_vertices_origin_edge_t **) PDM_hash_tab_data_get (htEdgeB, &_key);
      for (int j1 = 0; j1 < n_data; j1++) {
        free (data[j1]->vtxIntEdge);
        free (data[j1]->uIntEdge);
        free (data[j1]->coordsIntEdge);
      }
    }

    PDM_hash_tab_purge (htEdgeB, PDM_TRUE);

    PDM_g_num_t gnum_boxB = block_gnumB[i];

    int n_SubFaceFace = recv_stride1[i]/4;

    int newSize =  nSubFacesB + n_SubFaceFace;
    if (newSize > s_subFacesToFaceB) {
      while (newSize > s_subFacesToFaceB) {
        s_subFacesToFaceB += PDM_MAX (1, s_subFacesToFaceB/3);
      }
      subFacesToFaceB = realloc(subFacesToFaceB,
                                sizeof(PDM_g_num_t) * s_subFacesToFaceB);
      gNumSubFacesB = realloc(gNumSubFacesB,
                                sizeof(PDM_g_num_t) * s_subFacesToFaceB);
    }

    newSize =  3 * (nSubFacesB + n_SubFaceFace);
    if (newSize > s_subFacesToLinkA) {
      while (newSize > s_subFacesToLinkA) {
        s_subFacesToLinkA += PDM_MAX (1, s_subFacesToLinkA/3);
      }
      subFacesToLinkA = realloc(subFacesToLinkA,
                                sizeof(PDM_g_num_t) * s_subFacesToLinkA);
    }

    newSize =  nSubFacesB + n_SubFaceFace + 1;
    if (newSize > s_subFacesConnecIdxB) {
      while (newSize > s_subFacesConnecIdxB) {
        s_subFacesConnecIdxB += PDM_MAX (1, s_subFacesConnecIdxB/3);
      }
      subFacesConnecIdxB = realloc(subFacesConnecIdxB,
                                   sizeof(int) * s_subFacesConnecIdxB);
    }

    int ideb  = nSubFacesB;
    int ideb1 = 3 *nSubFacesB;
    int ideb2 = nSubFacesB + 1;

    /* Loop on sub-faces */

    for (int j = 0; j < n_SubFaceFace; j++) {

      int        n_vtx   = (int) recvSubFacesConnecBN[idx1++];
      int        iProcA = (int) recvSubFacesConnecBN[idx1++];
      int        iPartA = (int) recvSubFacesConnecBN[idx1++];
      PDM_g_num_t _gNumA  =       recvSubFacesConnecBN[idx1++];

      subFacesToFaceB[ideb]    = gnum_boxB;
      gNumSubFacesB[ideb++]    = _gNumA;

      subFacesToLinkA[ideb1++] = iProcA;
      subFacesToLinkA[ideb1++] = iPartA;
      subFacesToLinkA[ideb1++] = _gNumA;

      subFacesConnecIdxB[ideb2] = subFacesConnecIdxB[ideb2-1] + n_vtx;

      /* Loop on sub-edges */

      int iBeg = idx2;

      newSize =  subFacesConnecIdxB[ideb2];
      if (newSize > s_subFacesConnecB) {
        while (newSize > s_subFacesConnecB) {
          s_subFacesConnecB += PDM_MAX (1, s_subFacesConnecB/3);
        }
        subFacesConnecB = realloc(subFacesConnecB,
                                  sizeof(PDM_g_num_t) * s_subFacesConnecB);
      }

      newSize = 3 * subFacesConnecIdxB[ideb2];
      if (newSize > s_subFacesCoordsB) {
        while (newSize > s_subFacesCoordsB) {
          s_subFacesCoordsB += PDM_MAX (1, s_subFacesCoordsB/3);
        }
        subFacesCoordsB = realloc(subFacesCoordsB,
                                  sizeof(double) * s_subFacesCoordsB);
      }

      int iBeg1 = subFacesConnecIdxB[ideb2-1];
      int iBeg2 = 3 * subFacesConnecIdxB[ideb2-1];
      for (int k = 0; k < n_vtx; k++) {

        int next = iBeg + (k + 1) % n_vtx;
        _sub_edge_t *se = malloc(sizeof(_sub_edge_t));

        int imin = (recvSubFacesConnecB[idx2] < recvSubFacesConnecB[next]) ?  idx2 : next;
        int imax = (idx2 == imin) ? next : idx2;

        se->vtx1 = recvSubFacesConnecB[imin];
        se->vtx2 = recvSubFacesConnecB[imax];

        for (int k2 = 0; k2 < 3; k2++) {
          se->coords1[k2] = recvSubFacesCoordsB[3*imin+k2];
          se->coords2[k2] = recvSubFacesCoordsB[3*imax+k2];
        }

        PDM_g_num_t _key = (se->vtx1 + se->vtx2) / lComm;
        int key = (int) _key;

				PDM_hash_tab_data_add (htSubEdgeB, &key, (void *) se);

        subFacesConnecB[iBeg1++] = recvSubFacesConnecB[idx2];

        for (int k11 = 0; k11 < 3; k11++) {
          subFacesCoordsB[iBeg2++] = recvSubFacesCoordsB[3*idx2+k11];
        }

        idx2 += 1;

      }

      ideb2 += 1;

    }

    nSubFacesB += n_SubFaceFace;

    /************************************
     *                                  *
     * Add additional B sub-faces       *
     * Take into account initial edges  *
     *                                  *
     ************************************/

    int nEdge = recv_stride4[i]/3;

    if (nEdge > s_sum_vtx) {
      s_sum_vtx = nEdge;
      sum_vtx = realloc (sum_vtx, sizeof(PDM_g_num_t) * s_sum_vtx);
    }

    int n_vtx = 0;
    for (int j = 0; j < nEdge; j++) {


      int      nInt  = (int) recvFaceToEdgeNPtInt[idx4++];
      PDM_g_num_t n_vtx1 = recvFaceToEdgeNPtInt[idx4++];
      PDM_g_num_t n_vtx2 = recvFaceToEdgeNPtInt[idx4++];

      PDM_g_num_t s_vtx = n_vtx1 + n_vtx2;
      int _key_loc = s_vtx / lComm;

      int new_edge = 1;
      int n_data_key = PDM_hash_tab_n_data_get (htEdgeB, &_key_loc);

      _sub_vertices_origin_edge_t *curr_svoe = NULL;

      if (n_data_key > 0) {
        _sub_vertices_origin_edge_t **svoe_key =
          (_sub_vertices_origin_edge_t **) PDM_hash_tab_data_get (htEdgeB, &_key_loc);

        for (int j1 = 0; j1 < n_data_key; j1++) {
          curr_svoe = svoe_key[j1];
          if (((curr_svoe->vtx1 == n_vtx1) && (curr_svoe->vtx2 == n_vtx2)) ||
              ((curr_svoe->vtx2 == n_vtx1) && (curr_svoe->vtx1 == n_vtx2))) {
            new_edge = 0;
            break;
          }
        }
      }

      if (new_edge) {

        sum_vtx[n_vtx++]     = n_vtx1;

        _sub_vertices_origin_edge_t *svoe = malloc (sizeof(_sub_vertices_origin_edge_t));
        svoe->vtx1                = n_vtx1;
        svoe->vtx2                = n_vtx2;
        svoe->sIntEdge            = 4;
        svoe->nIntEdge            = 0;
        svoe->vtxIntEdge          = malloc (sizeof(PDM_g_num_t) * sIntEdge);
        svoe->uIntEdge            = malloc (sizeof(double) * sIntEdge);
        svoe->coordsIntEdge       = malloc (sizeof(double) * 3 * sIntEdge);

        PDM_hash_tab_data_add (htEdgeB, &_key_loc, svoe);

        curr_svoe = svoe;

        svoe->vtxIntEdge[svoe->nIntEdge] = n_vtx1;

        for (int k = 0; k < 3; k++) {
          svoe->coordsIntEdge[3 * svoe->nIntEdge + k] = recvFaceToEdgeCoordsvtx[idx41++];
        }

        svoe->uIntEdge[svoe->nIntEdge++] = 0.;

        svoe->vtxIntEdge[svoe->nIntEdge] = n_vtx2;

        for (int k = 0; k < 3; k++) {
          svoe->coordsIntEdge[3 * svoe->nIntEdge + k] = recvFaceToEdgeCoordsvtx[idx41++];
        }

        svoe->uIntEdge[svoe->nIntEdge++] = 1.;
      }

      else {
        for (int k = 0; k < 6; k++) {
          idx41++;
        }
      }

      for (int k = 0; k < nInt; k++) {

        PDM_g_num_t _originEdgeB = recvFaceToEdgeOrAndGnumPtInt[idx5++];
        PDM_g_num_t _gnumB       = recvFaceToEdgeOrAndGnumPtInt[idx5++];
        PDM_g_num_t _oNewPtB     = recvFaceToEdgeOrAndGnumPtInt[idx5++];
        double u                 = recvFaceToEdgeUPtInt[idx6++];
        double coords[3]         = {recvFaceToEdgeUPtInt[idx6],
                                    recvFaceToEdgeUPtInt[idx6+1],
                                    recvFaceToEdgeUPtInt[idx6+2]};
        idx6 += 3;

        if (n_vtx1 != _originEdgeB) {
          assert (n_vtx2 == _originEdgeB);
          u = 1 - u;
        }

        if (curr_svoe->nIntEdge >= curr_svoe->sIntEdge) {
          while (curr_svoe->nIntEdge >= curr_svoe->sIntEdge) {
            curr_svoe->sIntEdge += PDM_MAX (1, curr_svoe->sIntEdge/3);
          }

          curr_svoe->vtxIntEdge =
            realloc (curr_svoe->vtxIntEdge, sizeof(PDM_g_num_t) * curr_svoe->sIntEdge);
          curr_svoe->uIntEdge =
            realloc (curr_svoe->uIntEdge, sizeof(double) * curr_svoe->sIntEdge);
          curr_svoe->coordsIntEdge =
            realloc (curr_svoe->coordsIntEdge, sizeof(double) * 3 * curr_svoe->sIntEdge);
        }

        if (sIntEdge <= curr_svoe->sIntEdge) {
          sIntEdge = curr_svoe->sIntEdge;
          vtxIntEdgeSorted =
            realloc (vtxIntEdgeSorted, sizeof(PDM_g_num_t) * sIntEdge);
          coordsIntEdgeSorted =
            realloc (coordsIntEdgeSorted, sizeof(double) * 3 * sIntEdge);
          order =
            realloc (order, sizeof(int) * sIntEdge);
        }

        if (_oNewPtB != PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
          curr_svoe->vtxIntEdge[curr_svoe->nIntEdge] = _gnumB;
          curr_svoe->uIntEdge[curr_svoe->nIntEdge]   = u;

          for (int k11 = 0; k11 < 3; k11++) {
            curr_svoe->coordsIntEdge[3*curr_svoe->nIntEdge+k11] = coords[k11];
          }

          curr_svoe->nIntEdge += 1;

        }
      }
    }

    for (int k = 0; k < n_vtx; k++) {
      PDM_g_num_t n_vtx1 = sum_vtx[k];
      PDM_g_num_t n_vtx2 = sum_vtx[(k+1) % n_vtx];

      PDM_g_num_t s_vtx = n_vtx1 + n_vtx2;
      int _key_loc = s_vtx / lComm;

      int n_data_key = PDM_hash_tab_n_data_get (htEdgeB, &_key_loc);

      _sub_vertices_origin_edge_t *curr_svoe = NULL;

      int is_found = 0;
      if (n_data_key > 0) {
        _sub_vertices_origin_edge_t **svoe_key =
          (_sub_vertices_origin_edge_t **) PDM_hash_tab_data_get (htEdgeB, &_key_loc);

        for (int j1 = 0; j1 < n_data_key; j1++) {
          curr_svoe = svoe_key[j1];
          if (((curr_svoe->vtx1 == n_vtx1) && (curr_svoe->vtx2 == n_vtx2)) ||
              ((curr_svoe->vtx2 == n_vtx1) && (curr_svoe->vtx1 == n_vtx2))) {
            is_found = 1;
            break;
          }
        }
      }

      assert (is_found == 1);

      for (int k11 = 0; k11 < curr_svoe->nIntEdge; k11++) {
        order[k11] = k11;
      }

      PDM_sort_double (curr_svoe->uIntEdge, order, curr_svoe->nIntEdge);

      for (int k2 = 0; k2 < curr_svoe->nIntEdge; k2++) {
        vtxIntEdgeSorted[k2] = curr_svoe->vtxIntEdge[order[k2]];
        for (int k11 = 0; k11 < 3; k11++) {
          coordsIntEdgeSorted[3*k2+k11] = curr_svoe->coordsIntEdge[3*order[k2]+k11];
        }
      }

      /* Remove multiple points */

      int _nIntEdge = 0;
      for (int k2 = 0; k2 < curr_svoe->nIntEdge; k2++) {
        if (vtxIntEdgeSorted[k2] != vtxIntEdgeSorted[_nIntEdge]) {
          vtxIntEdgeSorted[++_nIntEdge] = vtxIntEdgeSorted[k2];
          for (int k11 = 0; k11 < 3; k11++) {
            coordsIntEdgeSorted[3*_nIntEdge+k11] = coordsIntEdgeSorted[3*k2+k11];
          }
          curr_svoe->uIntEdge[_nIntEdge] = curr_svoe->uIntEdge[k2];
        }
      }
      _nIntEdge++;

      faceIniVtxIdxB[i+1] += _nIntEdge - 1;

      int _newSize = faceIniVtxIdxB[i+1];
      if (_newSize > s_faceIniVtxB) {
        while (_newSize > s_faceIniVtxB) {
          s_faceIniVtxB += PDM_MAX (1, s_faceIniVtxB/3);
        }
        faceIniVtxB = realloc (faceIniVtxB,
                               sizeof(PDM_g_num_t) * s_faceIniVtxB);
      }

      for (int k3 = 0; k3 < _nIntEdge; k3++) {
        int next = (k3 + 1) % _nIntEdge;

        int imin = (vtxIntEdgeSorted[k3] < vtxIntEdgeSorted[next]) ?  k3 : next;
        int imax = (k3 == imin) ? next : k3;

        if (next != 0) {

          _sub_edge_t *se = malloc(sizeof(_sub_edge_t));
          faceIniVtxB[idx++] = vtxIntEdgeSorted[k3];
          se->vtx1 = vtxIntEdgeSorted[imin];
          se->vtx2 = vtxIntEdgeSorted[imax];

          for (int k2 = 0; k2 < 3; k2++) {
            se->coords1[k2] = coordsIntEdgeSorted[3*imin+k2];
            se->coords2[k2] = coordsIntEdgeSorted[3*imax+k2];
          }

          PDM_g_num_t _key = (se->vtx1 + se->vtx2) / lComm;
          int key = (int) _key;

          PDM_hash_tab_data_add (htSubEdgeB, &key, (void *) se);
        }
      }
    }

    /* Count the number of reference of each sub-edge */

    n_used_keys = PDM_hash_tab_n_used_keys_get (htSubEdgeB);
    used_keys = PDM_hash_tab_used_keys_get (htSubEdgeB);

    int t_n_data = 0;
    for (int j = 0; j < n_used_keys; j++) {
      int _key = (int) used_keys[j];
      int n_data = PDM_hash_tab_n_data_get (htSubEdgeB, &_key);
      t_n_data += n_data;
    }

    PDM_g_num_t *oneRef = malloc(sizeof(PDM_g_num_t) * 2 * t_n_data);
    double     *coordsOneRef = malloc(sizeof(double) * 6 * t_n_data);

    int nOneRef = 0;

    for (int j = 0; j < n_used_keys; j++) {
      int _key = (int) used_keys[j];
      int n_data = PDM_hash_tab_n_data_get (htSubEdgeB, &_key);
      int tag[n_data];
      PDM_array_reset_int(tag, n_data, 0);

      _sub_edge_t **se = (_sub_edge_t **) PDM_hash_tab_data_get (htSubEdgeB, &_key);
      for (int k11 = 0; k11 < n_data; k11++) {
        if (tag[k11] == 1) continue;
        _sub_edge_t *_se1 = se[k11];
        for (int k2 = k11+1; k2 < n_data; k2++) {
          if (tag[k2] == 1) continue;
          _sub_edge_t *_se2 = se[k2];
          if ((_se1->vtx1 == _se2->vtx1) &&
              (_se1->vtx2 == _se2->vtx2)) {
            tag[k11] = 1;
            tag[k2] = 1;
          }
        }
      }
      for (int k11 = 0; k11 < n_data; k11++) {
        if (tag[k11] == 0) {
          _sub_edge_t *_se1 = se[k11];
          oneRef[2 * nOneRef    ] = _se1->vtx1;
          oneRef[2 * nOneRef + 1] = _se1->vtx2;

          for (int k2 = 0; k2 < 3; k2++) {
            coordsOneRef[6 * nOneRef + k2] = _se1->coords1[k2];
          }

          for (int k2 = 0; k2 < 3; k2++) {
            coordsOneRef[6 * nOneRef + 3 + k2] = _se1->coords2[k2];
          }
          nOneRef += 1;
        }
      }
    }


    if (nOneRef < 2) {
      if (vb == 1 && nOneRef == 1) {
        printf("[%6d] Warning : nOneRef = 1\n", i_rank);
      }
      facesToSubFacesBIdx[i+1] = facesToSubFacesBIdx[i]
        + n_SubFaceFace;
      free(oneRef);
      free(coordsOneRef);
    }
    else {

      PDM_g_num_t *tag = malloc(sizeof(PDM_g_num_t) * nOneRef);
      int nAddSubFace = 0;
      int *addSubFaceIdx = malloc(sizeof(int) * (nOneRef + 1));
      PDM_g_num_t *addSubFace = malloc(sizeof(PDM_g_num_t) * nOneRef);
      double *addSubFaceCoords = malloc(sizeof(double) * 3 * nOneRef);
      for (int k11 = 0; k11 < nOneRef; k11++) {
        tag[k11] = 0;
        addSubFaceIdx[k11] = 0;
      }

      /* Add additional sub-faces */

      PDM_g_num_t iniVal = -1;
      PDM_g_num_t nextVal = -1;
      int idx_addSubFace1 = 0;
      int idx_addSubFace = 0;
      for (int k11 = 0; k11 < nOneRef; k11++) {
        if (tag[k11] == 1) continue;

        iniVal = oneRef[2*k11];
        nextVal = oneRef[2*k11 + 1];

        addSubFace[idx_addSubFace1++] = iniVal;
        addSubFace[idx_addSubFace1++] = oneRef[2*k11 + 1];

        for (int k3 = 0; k3 < 3; k3++) {
          addSubFaceCoords[idx_addSubFace++] = coordsOneRef[6*k11+k3];
        }
        for (int k3 = 0; k3 < 3; k3++) {
          addSubFaceCoords[idx_addSubFace++] = coordsOneRef[6*k11+3+k3];
        }

        nAddSubFace += 1;
        addSubFaceIdx[nAddSubFace] += 2;

        int k_prev = k11;
        int count = 0;
        while (iniVal != nextVal) {
          count++;
          if (count > 1000) {
            printf("[%6d] error nOneRef = %d, gnumB = "PDM_FMT_G_NUM", k11 = %d, iniVal = "PDM_FMT_G_NUM", nextVal = "PDM_FMT_G_NUM", oneRef[2*k11 + 1] = "PDM_FMT_G_NUM"\n", i_rank, nOneRef, block_gnumB[i], k11, iniVal, nextVal, oneRef[2*k11 + 1]);
            abort();
          }
          int found = 0;
          for (int k2 = 0; k2 < nOneRef; k2++) {
            if ((tag[k2] == 1) || (k_prev == k2)) continue;
            if (nextVal == oneRef[2*k2]) {
              nextVal = oneRef[2*k2 + 1];
              if (nextVal != iniVal) {
                addSubFace[idx_addSubFace1++] = nextVal;
                for (int k3 = 0; k3 < 3; k3++) {
                  addSubFaceCoords[idx_addSubFace++] = coordsOneRef[6*k2+3+k3];
                }
                addSubFaceIdx[nAddSubFace]++;
              }
              tag[k2] = 1;
              k_prev = k2;
              found = 1;
              break;
            }
            else if (nextVal == oneRef[2*k2 + 1]) {
              nextVal = oneRef[2*k2];
              if (nextVal != iniVal) {
                addSubFace[idx_addSubFace1++] = nextVal;
                for (int k3 = 0; k3 < 3; k3++) {
                  addSubFaceCoords[idx_addSubFace++] = coordsOneRef[6*k2+k3];
                }
                addSubFaceIdx[nAddSubFace]++;
              }
              tag[k2] = 1;
              k_prev = k2;
              found = 1;
              break;
            }
          }
          if (!found) {
            printf("[%6d] !found nOneRef = %d, gnumB = "PDM_FMT_G_NUM", k11 = %d, iniVal = "PDM_FMT_G_NUM", nextVal = "PDM_FMT_G_NUM", oneRef[2*k11 + 1] = "PDM_FMT_G_NUM"\n", i_rank, nOneRef, block_gnumB[i], k11, iniVal, nextVal, oneRef[2*k11 + 1]);
            abort();
          }
        }
      }

      free (tag);

      for (int k11 = 0; k11 < nAddSubFace; k11++) {
        addSubFaceIdx[k11+1] += addSubFaceIdx[k11];
      }

      newSize = nSubFacesB + nAddSubFace + 1;
      if (newSize > s_subFacesConnecIdxB) {
        while (newSize > s_subFacesConnecIdxB) {
          s_subFacesConnecIdxB += PDM_MAX (1, s_subFacesConnecIdxB/3);
        }
        subFacesConnecIdxB = realloc(subFacesConnecIdxB,
                                     sizeof(int) * s_subFacesConnecIdxB);
      }

      newSize =  nSubFacesB + nAddSubFace;
      if (newSize > s_subFacesToFaceB) {
        while (newSize > s_subFacesToFaceB) {
          s_subFacesToFaceB += PDM_MAX (1, s_subFacesToFaceB/3);
        }
        subFacesToFaceB = realloc(subFacesToFaceB,
                                  sizeof(PDM_g_num_t) * s_subFacesToFaceB);
        gNumSubFacesB = realloc(gNumSubFacesB,
                                sizeof(PDM_g_num_t) * s_subFacesToFaceB);
      }

      newSize =  3 * (nSubFacesB + nAddSubFace);
      if (newSize > s_subFacesToLinkA) {
        while (newSize > s_subFacesToLinkA) {
          s_subFacesToLinkA += PDM_MAX (1, s_subFacesToLinkA/3) ;
        }
        subFacesToLinkA = realloc(subFacesToLinkA,
                                  sizeof(PDM_g_num_t) * s_subFacesToLinkA);
      }

      n_t_nAddSubFace += nAddSubFace;

      for (int j = 0; j < nAddSubFace; j++) {

        int _n_vtx = addSubFaceIdx[j+1] - addSubFaceIdx[j];
        subFacesToFaceB[nSubFacesB] = gnum_boxB;
        gNumSubFacesB[nSubFacesB] = -1;

        subFacesToLinkA[3*nSubFacesB] = -1;
        subFacesToLinkA[3*nSubFacesB+1] = -1;
        subFacesToLinkA[3*nSubFacesB+2] = -1;

        subFacesConnecIdxB[nSubFacesB+1] =
          subFacesConnecIdxB[nSubFacesB] + _n_vtx;

        nSubFacesB += 1;

        newSize = subFacesConnecIdxB[nSubFacesB];
        if (newSize > s_subFacesConnecB) {
          while (newSize > s_subFacesConnecB) {
            s_subFacesConnecB += PDM_MAX(1, s_subFacesConnecB/3) ;
          }
          subFacesConnecB = realloc(subFacesConnecB,
                                    sizeof(PDM_g_num_t) * s_subFacesConnecB);
        }

        newSize = 3 * subFacesConnecIdxB[nSubFacesB];
        if (newSize > s_subFacesCoordsB) {
          while (newSize > s_subFacesCoordsB) {
            s_subFacesCoordsB += PDM_MAX(1, s_subFacesCoordsB/3);
          }
          subFacesCoordsB = realloc(subFacesCoordsB,
                                    sizeof(double) * s_subFacesCoordsB);
        }

        int _ideb2 =     subFacesConnecIdxB[nSubFacesB-1];
        int ideb3 = 3 * subFacesConnecIdxB[nSubFacesB-1];

        for (int k = addSubFaceIdx[j]; k < addSubFaceIdx[j+1]; k++) {
          subFacesConnecB[_ideb2++] = addSubFace[k];
          for (int k11 = 0; k11 < 3; k11++) {
            subFacesCoordsB[ideb3++] = addSubFaceCoords[3*k+k11];
          }
        }
      }

      facesToSubFacesBIdx[i+1] = facesToSubFacesBIdx[i]
        + n_SubFaceFace
        + nAddSubFace;
      /* Cleanup */

      free (oneRef);
      free (coordsOneRef);
      free (addSubFaceIdx);
      free (addSubFace);
      free (addSubFaceCoords);
    }

  }
  if (vb) {
    printf("[%6d] >> hash_tab purge, free\n", i_rank);
    fflush(stdout);
  }

  if (sum_vtx != NULL) {
    free (sum_vtx);
  }

  int n_used_keys = PDM_hash_tab_n_used_keys_get (htEdgeB);
  PDM_g_num_t *used_keys = PDM_hash_tab_used_keys_get (htEdgeB);

  for (int j = 0; j < n_used_keys; j++) {
    int _key = (int) used_keys[j];
    int n_data = PDM_hash_tab_n_data_get (htEdgeB, &_key);
    _sub_vertices_origin_edge_t **data =
      (_sub_vertices_origin_edge_t **) PDM_hash_tab_data_get (htEdgeB, &_key);
    for (int j1 = 0; j1 < n_data; j1++) {
      free (data[j1]->vtxIntEdge);
      free (data[j1]->uIntEdge);
      free (data[j1]->coordsIntEdge);
    }
  }
  PDM_hash_tab_purge (htEdgeB, PDM_TRUE);
  PDM_hash_tab_free (htEdgeB);

  PDM_hash_tab_purge (htSubEdgeB, PDM_TRUE);
  PDM_hash_tab_free (htSubEdgeB);

  /*
   * Update memory
   */
  if (vb) {
    printf("[%6d] >> Realloc\n", i_rank);
    fflush(stdout);
  }
  subFacesConnecB = realloc (subFacesConnecB, sizeof(PDM_g_num_t) * subFacesConnecIdxB[nSubFacesB]);
  subFacesCoordsB = realloc (subFacesCoordsB, sizeof(double) * 3 * subFacesConnecIdxB[nSubFacesB]);
  subFacesConnecIdxB = realloc (subFacesConnecIdxB, sizeof(int) * (nSubFacesB + 1));
  subFacesToFaceB = realloc (subFacesToFaceB, sizeof(PDM_g_num_t) * nSubFacesB);
  gNumSubFacesB = realloc (gNumSubFacesB, sizeof(PDM_g_num_t) * nSubFacesB);
  subFacesToLinkA = realloc (subFacesToLinkA, sizeof(PDM_g_num_t) * 3 * nSubFacesB);

  if (1 == 0) {
    printf("\n-- Liste des sous-facettes des elements de B (intersection + complement intersection)\n");
    for (int i = 0; i < n_elt_blockB; i++) {

      PDM_g_num_t gnum_boxB = block_gnumB[i];

      printf ("face origin B : "PDM_FMT_G_NUM"\n", gnum_boxB);

      for (int j = facesToSubFacesBIdx[i]; j < facesToSubFacesBIdx[i+1]; j++) {
        printf (" - sub face %d :", j);
        for (int k = subFacesConnecIdxB[j]; k < subFacesConnecIdxB[j+1]; k++) {
          printf(" "PDM_FMT_G_NUM, subFacesConnecB[k]);
          for (int k11 = 0; k11 < 3; k11++) {
            printf(" %12.5e", subFacesCoordsB[3*k+k11]);
          }
          printf("/");

        }
        printf("\n");
      }
    }

    printf("\n-- Faces d'origine avec nouveaux points (intersection + complement intersection)\n");
    for (int i = 0; i < n_elt_blockB; i++) {

      PDM_g_num_t gnum_boxB = block_gnumB[i];

      printf ("face origin B : "PDM_FMT_G_NUM"\n", gnum_boxB);
      for (int j = faceIniVtxIdxB[i]; j < faceIniVtxIdxB[i+1]; j++) {
        printf(" "PDM_FMT_G_NUM, faceIniVtxB[j]);
      }
      printf("\n");
    }
  }

  if (vb) {
    printf("[%6d] >> Scan n_t_nAddSubFace\n", i_rank);
    fflush(stdout);
  }
  PDM_MPI_Scan (&n_t_nAddSubFace, &beg_nAddSubFaces,
            1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  beg_nAddSubFaces += -n_t_nAddSubFace + 1;

  int j1 = 0;
  _max_gnum_loc = end_subFaces;
  for (int i = 0; i < nSubFacesB; i++) {
    if (gNumSubFacesB[i] == -1) {
      gNumSubFacesB[i] = end_subFaces + beg_nAddSubFaces + j1++ ;
      _max_gnum_loc = gNumSubFacesB[i];
    }
  }

  PDM_g_num_t nTSubFacesB;
  /* _max_gnum_loc = -1; */
  /* if (nSharedSubFaces > 0) { */
  /*   _max_gnum_loc = gNumSubFacesA[nSharedSubFaces-1]; */
  /* } */

  if (vb) {
    printf("[%6d] >> Allreduce _max_gnum_loc\n", i_rank);
    fflush(stdout);
  }
  PDM_MPI_Allreduce(&_max_gnum_loc, &nTSubFacesB, 1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ol->comm);

  /*
   * Cleanup
   */

  free (coordsIntEdgeSorted);
  free (vtxIntEdgeSorted);
  free (order);

  // sous-facettes calculees

  free (recv_stride1);
  free (recvSubFacesConnecBN);

  free (recv_stride2);
  free (recvSubFacesConnecB);

  free (recv_stride21);
  free (recvSubFacesCoordsB);

  // arete initiales du contour

  free (recv_stride4);
  free (recvFaceToEdgeNPtInt);

  free (recv_stride41);
  free (recvFaceToEdgeCoordsvtx);

  free (recvFaceToEdgeOrAndGnumPtInt);

  free (recvFaceToEdgeUPtInt);

  free (blockA_boxesB_idx);
  free (blockA_boxesB_gnum_data);
  free (blockA_boxesB_lnum_data);

  for (int i = 0; i < 2; i++) {
    free (faceStrideCurrent[i]);
    free (faceStrideCurrent3[i]);
    free (faceToEdgeCurrent[i]);
    free (faceToVtxCurrent[i]);
    free (face_vtxCooCurrent[i]);
    free (face_vtxEpsCurrent[i]);
  }

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_COMPUTE_ADD_SUB_FACESB] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_COMPUTE_ADD_SUB_FACESB]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_COMPUTE_ADD_SUB_FACESB]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_COMPUTE_ADD_SUB_FACESB]   = PDM_timer_cpu_sys(ol->timer);

  PDM_timer_resume(ol->timer);
  PDM_MPI_Barrier (ol->comm);
  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_compute_add_sub_facesb !!!!\n");
    fflush(stdout);
  }

  /* PDM_MPI_Barrier (ol->comm); */
  /* printf ("!!!! Arret temporaire apres calcul des sous-faces complementaires de B !!!!\n"); */
  /* if (i_rank == 0) { */
  /*   _ol_dump_times (ol); */
  /* } */
  /* exit(0); */

  /***********************************************************************************
   *                                                                                 *
   * Transfer results to origin distribution for mesh A :                            *
   * For each boxe :                                                                 *
   *     - Size of connectivity of each subface + numabs + grpah comm (iproc i_part) *
   *     - Connectivity of each subface (numabs)                                     *
   *     - Coordinates of vertices                                                   *
   *                                                                                 *
   * For each face, additional faces and sub-faces have to be merge                  *
   *                                                                                 *
   ***********************************************************************************/

  PDM_g_num_t *firstSend = malloc (sizeof(PDM_g_num_t) *
                          facesToAddSubFacesAIdx[n_elt_blockA] * 5);

  int        *firstSendStride = malloc (sizeof(int) * n_elt_blockA);

  int **firstRecvStrideA = malloc (sizeof(int *) * n_partA);

  for (int i = 0; i < n_partA; i++) {
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    firstRecvStrideA[i] = malloc(sizeof(int) * n_face);
  }

  PDM_g_num_t **firstRecvA = malloc (sizeof(PDM_g_num_t *) * n_partA);

  idx = 0;
  int n_T_vertex = 0;

  for (int i = 0; i < n_elt_blockA; i++) {
    blockA_lnum_data[i] = PDM_binary_search_long (gnum_eltA[i], block_gnumA, n_eltA);
  }

  for (int i = 0; i < n_elt_blockA; i++) {

    int lnum_boxA = blockA_lnum_data[i];

    int nAddSubFaces = facesToAddSubFacesAIdx[lnum_boxA+1] - facesToAddSubFacesAIdx[lnum_boxA];
    int _nSubFaces = facesToSubFacesAIdx[lnum_boxA+1] - facesToSubFacesAIdx[lnum_boxA];

    firstSendStride[i] =  5 * (_nSubFaces + nAddSubFaces);

    for (int j = facesToSubFacesAIdx[lnum_boxA]; j < facesToSubFacesAIdx[lnum_boxA+1]; j++) {

      int n_vtx = subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
      n_T_vertex += n_vtx;

      PDM_g_num_t iProcB = subFacesToFaces[5*j+2];
      PDM_g_num_t iPartB = subFacesToFaces[5*j+3];
      PDM_g_num_t numAbs = gNumSubFacesA[j];
      PDM_g_num_t oNumAbs = block_gnumA[lnum_boxA];

      firstSend[idx++] = n_vtx;
      firstSend[idx++] = numAbs;
      firstSend[idx++] = oNumAbs;
      firstSend[idx++] = iProcB;
      firstSend[idx++] = iPartB;

    }

    for (int j = facesToAddSubFacesAIdx[lnum_boxA]; j < facesToAddSubFacesAIdx[lnum_boxA+1]; j++) {

      int n_vtx = subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
      n_T_vertex += n_vtx;
      PDM_g_num_t iProcB = -1;
      PDM_g_num_t iPartB = -1;
      PDM_g_num_t numAbs = gNumSubFacesA[j];
      PDM_g_num_t oNumAbs = block_gnumA[lnum_boxA];

      firstSend[idx++] = n_vtx;
      firstSend[idx++] = numAbs;
      firstSend[idx++] = oNumAbs;
      firstSend[idx++] = iProcB;
      firstSend[idx++] = iPartB;

    }
  }

  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride,
                                           firstSend,
                                           firstRecvStrideA,
                                           (void ** )firstRecvA);

  PDM_g_num_t *secondSend = malloc (sizeof(PDM_g_num_t) * n_T_vertex);
  double     *thirdSend = malloc (3 * sizeof(double) * n_T_vertex);

  idx = 0;
  idx2 = 0;
  for (int i = 0; i < n_elt_blockA; i++) {

    int lnum_boxA = blockA_lnum_data[i];

    firstSendStride[i] = 0;

    for (int j = facesToSubFacesAIdx[lnum_boxA]; j < facesToSubFacesAIdx[lnum_boxA+1]; j++) {

      int n_vtx = subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
      firstSendStride[i] += n_vtx;

      for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
        secondSend[idx++] = subFacesConnecA[k];
        for (int k11 = 0; k11 < 3; k11++) {
          thirdSend[idx2++] = subFacesCoordsA[3*k+k11];
        }
      }

    }

    for (int j = facesToAddSubFacesAIdx[lnum_boxA]; j < facesToAddSubFacesAIdx[lnum_boxA+1]; j++) {

      int n_vtx = subFacesConnecIdx[j+1] - subFacesConnecIdx[j];
      firstSendStride[i] += n_vtx;

      for (int k = subFacesConnecIdx[j]; k < subFacesConnecIdx[j+1]; k++) {
        secondSend[idx++] = subFacesConnecA[k];
        for (int k11 = 0; k11 < 3; k11++) {
          thirdSend[idx2++] = subFacesCoordsA[3*k+k11];
        }
      }
    }
  }

  int **secondRecvStrideA = malloc (sizeof(int *) * n_partA);

  for (int i = 0; i < n_partA; i++) {
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    secondRecvStrideA[i] = malloc(sizeof(int) * n_face);
  }

  PDM_g_num_t **secondRecvA = malloc (sizeof(PDM_g_num_t *) * n_partA);

  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride,
                                           secondSend,
                                           secondRecvStrideA,
                                           (void ** ) secondRecvA);

  //

  for (int i = 0; i < n_elt_blockA; i++) {
    firstSendStride[i] *= 3;
  }

  int **thirdRecvStrideA = malloc (sizeof(int *) * n_partA);

  for (int i = 0; i < n_partA; i++) {
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    thirdRecvStrideA[i] = malloc(sizeof(int) * n_face);
  }

  double **thirdRecvA = malloc (sizeof(double *) * n_partA);

  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           sizeof(double),
                                           firstSendStride,
                                           thirdSend,
                                           thirdRecvStrideA,
                                           (void ** )thirdRecvA);

  //

  PDM_g_num_t *fourthSend = malloc (sizeof(PDM_g_num_t) * faceIniVtxIdxA[n_elt_blockA]);

  idx=0;
  for (int i = 0; i < n_elt_blockA; i++) {
    int lnum_boxA = blockA_lnum_data[i];

    firstSendStride[i] = faceIniVtxIdxA[lnum_boxA+1] - faceIniVtxIdxA[lnum_boxA];
    for (int j = faceIniVtxIdxA[lnum_boxA]; j < faceIniVtxIdxA[lnum_boxA+1]; j++) {
      fourthSend[idx++] = faceIniVtxA[j];
    }
  }

  free (faceIniVtxIdxA);
  free (faceIniVtxA);

  int **fourthRecvStrideA = malloc (sizeof(int *) * n_partA);
  for (int i = 0; i < n_partA; i++) {
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    fourthRecvStrideA[i] = malloc(sizeof(int) * n_face);
  }

  PDM_g_num_t **fourthRecvA = malloc (sizeof(PDM_g_num_t *) * n_partA);

  PDM_box_set_send_data_to_origin_distrib (boxesA,
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride,
                                           fourthSend,
                                           fourthRecvStrideA,
                                           (void ** )fourthRecvA);

  //

  PDM_box_set_destroy (&boxesA);

  free (firstSend);
  free (secondSend);
  free (thirdSend);
  free (fourthSend);
  free (block_gnumA);
  free (firstSendStride);
  free (blockA_lnum_data);

  free (subFacesConnecIdx);
  free (subFacesCoordsA);
  free (facesToSubFacesAIdx);
  free (facesToAddSubFacesAIdx);
  free (gNumSubFacesA);
  free (subFacesToFaces);
  free (subFacesConnecA);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTA] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTA]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_SEND_RESULTS_TO_INIT_PARTA]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_SEND_RESULTS_TO_INIT_PARTA]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  PDM_MPI_Barrier (ol->comm);
  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_send_results_to_init_parta !!!!\n");
    fflush(stdout);
  }

  /* PDM_MPI_Barrier (ol->comm); */
  /* printf ("!!!! Arret temporaire apres envoi des resultats au partitions d'origine de A !!!!\n"); */
  /* if (i_rank == 0) { */
  /*   _ol_dump_times (ol); */
  /* } */
  /* exit(0); */

  /*****************************************************************************
   *                                                                           *
   * Redistribute B Boxes to send results to B                                 *
   *                                                                           *
   ****************************************************************************/

  int box_size = PDM_box_set_get_size(boxesB); //***
  const PDM_g_num_t *boxes_gnum = PDM_box_set_get_g_num(boxesB);

  PDM_g_num_t *distrib_index = PDM_part_to_block_distrib_index_get(ptb_boxesB);

  block_gnumB = PDM_part_to_block_block_gnum_get (ptb_boxesB);

  n_g_eltB = PDM_box_set_get_global_size(boxesB);

  distribB = PDM_box_distrib_create (box_size,
                                     n_g_eltB,
                                     1, // Don't use in this case
                                     ol->comm);

  countEltsB = PDM_array_zeros_int(lComm);
  PDM_array_reset_int(distribB->index, lComm+1, 0);

  for (int i = 0; i < box_size; i++) {
    int iProc = PDM_binary_search_gap_long (boxes_gnum[i] - 1,
                                            distrib_index,
                                            lComm + 1) + 1;
    distribB->index[iProc]++;
  }

  for (int i = 0; i < lComm; i++) {
    distribB->index[i+1] += distribB->index[i];
  }

  distribB->list = (int *) malloc (sizeof(int) * distribB->index[lComm]);

  for (int i = 0; i < box_size; i++) {
    int iProc = PDM_binary_search_gap_long (boxes_gnum[i] - 1,
                                            distrib_index,
                                            lComm + 1);
    int idxB2 = distribB->index[iProc] + (countEltsB[iProc]++);
    distribB->list[idxB2] = i;
  }

  free (countEltsB);

  PDM_box_distrib_clean (distribB);
  PDM_box_set_redistribute (distribB, boxesB);

  /*
   * Cleanup
   */

  PDM_box_distrib_destroy (&distribB);

  PDM_box_set_remove_duplicate (boxesB);

  box_size = PDM_box_set_get_size(boxesB);
  boxes_gnum = PDM_box_set_get_g_num(boxesB);

  blockB_lnum_data = (int *) realloc (blockB_lnum_data, sizeof(int) * n_elt_blockB);

  for (int i = 0; i < n_elt_blockB; i++) {
    blockB_lnum_data[i] = PDM_binary_search_long (boxes_gnum[i], block_gnumB, n_elt_blockB);
  }

  /*****************************************************************************
   *                                                                           *
   * Transfer results to origin distribution for mesh B                        *
   *                                                                           *
   ****************************************************************************/

  firstSend = malloc (sizeof(PDM_g_num_t) *
                      facesToSubFacesBIdx[n_elt_blockB] * 5);

  firstSendStride = malloc (sizeof(int) * n_elt_blockB);

  int **firstRecvStrideB = malloc (sizeof(int *) * n_partB);

  for (int i = 0; i < n_partB; i++) {
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    firstRecvStrideB[i] = malloc(sizeof(int) * n_face);
  }

  PDM_g_num_t **firstRecvB = malloc (sizeof(PDM_g_num_t *) * n_partB);

  idx = 0;
  n_T_vertex = 0;

  for (int i = 0; i < n_elt_blockB; i++) {

    int lnum_boxB = blockB_lnum_data[i];
    int _nSubFaces = facesToSubFacesBIdx[lnum_boxB+1] - facesToSubFacesBIdx[lnum_boxB];

    firstSendStride[i] = 5 * _nSubFaces;

    for (int j = facesToSubFacesBIdx[lnum_boxB]; j < facesToSubFacesBIdx[lnum_boxB+1]; j++) {

      int n_vtx = subFacesConnecIdxB[j+1] - subFacesConnecIdxB[j];
      n_T_vertex += n_vtx;
      int iProcA = (int) subFacesToLinkA[3*j];
      int iPartA = (int) subFacesToLinkA[3*j + 1];
      PDM_g_num_t numAbs = gNumSubFacesB[j];
      PDM_g_num_t oNumAbs = block_gnumB[lnum_boxB];

      firstSend[idx++] = n_vtx;
      firstSend[idx++] = numAbs;
      firstSend[idx++] = oNumAbs;
      firstSend[idx++] = iProcA;
      firstSend[idx++] = iPartA;

    }
  }

  PDM_part_to_block_free (ptb_boxesB);
  free (blockA_boxesB_gnum_data_cp);

  PDM_box_set_send_data_to_origin_distrib (boxesB,
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride,
                                           firstSend,
                                           firstRecvStrideB,
                                           (void ** )firstRecvB);


  free (firstSend);

  secondSend = malloc (sizeof(PDM_g_num_t) * n_T_vertex);
  thirdSend = malloc (3 * sizeof(double) * n_T_vertex);

  idx = 0;
  idx2 = 0;

  for (int i = 0; i < n_elt_blockB; i++) {

    int lnum_boxB = blockB_lnum_data[i];
    firstSendStride[i] = 0;

    for (int j = facesToSubFacesBIdx[lnum_boxB]; j < facesToSubFacesBIdx[lnum_boxB+1]; j++) {

      int n_vtx = subFacesConnecIdxB[j+1] - subFacesConnecIdxB[j];

      firstSendStride[i] += n_vtx;

      for (int k = subFacesConnecIdxB[j]; k < subFacesConnecIdxB[j+1]; k++) {
        secondSend[idx++] = subFacesConnecB[k];
        for (int k1 = 0; k1 < 3; k1++) {
          thirdSend[idx2++] = subFacesCoordsB[3*k+k1];
        }
      }
    }
  }

  free (subFacesCoordsB);

  int **secondRecvStrideB = malloc (sizeof(int *) * n_partB);

  for (int i = 0; i < n_partB; i++) {
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    secondRecvStrideB[i] = malloc(sizeof(int) * n_face);
  }

  PDM_g_num_t **secondRecvB = malloc (sizeof(PDM_g_num_t *) * n_partB);

  PDM_box_set_send_data_to_origin_distrib (boxesB,
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride,
                                           secondSend,
                                           secondRecvStrideB,
                                           (void ** )secondRecvB);

  for (int i = 0; i < n_elt_blockB; i++) {
    firstSendStride[i] *= 3;
  }

  int **thirdRecvStrideB = malloc (sizeof(int *) * n_partB);

  for (int i = 0; i < n_partB; i++) {
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    thirdRecvStrideB[i] = malloc(sizeof(int) * n_face);
  }

  double **thirdRecvB = malloc (sizeof(double *) * n_partB);

  PDM_box_set_send_data_to_origin_distrib (boxesB,
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           sizeof(double),
                                           firstSendStride,
                                           thirdSend,
                                           thirdRecvStrideB,
                                           (void ** )thirdRecvB);

  fourthSend = malloc (sizeof(PDM_g_num_t) * faceIniVtxIdxB[n_elt_blockB]);

  idx=0;
  for (int i = 0; i < n_elt_blockB; i++) {
    int lnum_boxB = blockB_lnum_data[i];
    firstSendStride[i] = faceIniVtxIdxB[lnum_boxB+1] - faceIniVtxIdxB[lnum_boxB];
    for (int j = faceIniVtxIdxB[lnum_boxB]; j < faceIniVtxIdxB[lnum_boxB+1]; j++) {
      fourthSend[idx++] = faceIniVtxB[j];
    }
  }

  free (faceIniVtxIdxB);
  free (faceIniVtxB);

  int **fourthRecvStrideB = malloc (sizeof(int *) * n_partB);
  for (int i = 0; i < n_partB; i++) {
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    fourthRecvStrideB[i] = malloc(sizeof(int) * n_face);
  }

  PDM_g_num_t **fourthRecvB = malloc (sizeof(PDM_g_num_t *) * n_partB);

  PDM_box_set_send_data_to_origin_distrib (boxesB,
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           sizeof(PDM_g_num_t),
                                           firstSendStride,
                                           fourthSend,
                                           fourthRecvStrideB,
                                           (void ** )fourthRecvB);

  free (firstSendStride);
  free (secondSend);
  free (thirdSend);
  free (fourthSend);
  free (blockB_lnum_data);

  free (subFacesConnecB);
  free (subFacesConnecIdxB);

  free (subFacesToFaceB);
  free (gNumSubFacesB);
  free (subFacesToLinkA);

  free (facesToSubFacesBIdx);

  PDM_box_set_destroy (&boxesB);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_SEND_RESULTS_TO_INIT_PARTB] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_SEND_RESULTS_TO_INIT_PARTB]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_SEND_RESULTS_TO_INIT_PARTB]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_SEND_RESULTS_TO_INIT_PARTB]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  PDM_MPI_Barrier (ol->comm);
  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_send_results_to_init_partb !!!!\n");
    fflush(stdout);
  }

  /* PDM_MPI_Barrier (ol->comm); */
  /* printf ("!!!! Arret temporaire apres envoi des resultats au partitions d'origine de B !!!!\n"); */
  /* if (i_rank == 0) { */
  /*   _ol_dump_times (ol); */
  /* } */
  /* exit(0); */

  /*****************************************************************************
   *                                                                           *
   * Compute A local numbering for vertices and faces                          *
   *                                                                           *
   ****************************************************************************/

  PDM_g_num_t nUnChangedFaceA = 0;
  PDM_g_num_t nSubFaceA = 0;
  int *nUnChangedFacePartA = malloc(sizeof(int) * n_partA);
  int *nSubFacePartA = malloc(sizeof(int) * n_partA);
  int *s_olface_vtxA = PDM_array_zeros_int(n_partA);


  /*
   * Compute dimensions
   */

  for (int i = 0; i < n_partA; i++) {
    nUnChangedFacePartA[i] = 0;
    nSubFacePartA[i] = 0;
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    const PDM_g_num_t *gNumFaceA = PDM_surf_mesh_part_face_g_num_get (ol->meshA, i);
    const int *iniface_vtx_idxA =  PDM_surf_mesh_part_face_vtx_idx_get (ol->meshA, i);

    idx = 0;
    for (int j = 0; j < n_face; j++) {
      if (firstRecvStrideA[i][j] == 0) {
        nUnChangedFacePartA[i] += 1;
        s_olface_vtxA[i] += iniface_vtx_idxA[j+1] - iniface_vtx_idxA[j];
      }
      else {
        nSubFacePartA[i] += firstRecvStrideA[i][j]/5;
        for (int k = 0; k < firstRecvStrideA[i][j]/5; k++) {
          int n_vtx    = (int) firstRecvA[i][idx++];
          idx++;
          PDM_g_num_t oNumAbs = firstRecvA[i][idx++];
          idx += 2;

          s_olface_vtxA[i] += n_vtx;

          assert(oNumAbs == gNumFaceA[j]);
        }
      }
    }
    nUnChangedFaceA += nUnChangedFacePartA[i];
    nSubFaceA += nSubFacePartA[i];

  }

  PDM_g_num_t begUnChangedFaceA = 0;
  PDM_MPI_Scan (&nUnChangedFaceA, &begUnChangedFaceA,
            1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  begUnChangedFaceA += -nUnChangedFaceA + 1;

  PDM_g_num_t n_faceA = nUnChangedFaceA + nSubFaceA;
  PDM_g_num_t nGFaceA;
  PDM_MPI_Allreduce(&n_faceA, &nGFaceA, 1,
                PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  PDM_g_num_t nGVtxA = n_g_newVtxA;
  ol->olMeshA = _ol_mesh_create (nGFaceA, nGVtxA, n_partA);

  /*
   * Build local connectivity
   */

  int idxIpart = 0;

  PDM_g_num_t iniNGVtx = PDM_surf_mesh_n_g_vtx_get (ol->meshA);

  for (int i = 0; i < n_partA; i++) {

    _ol_part_t *olp = _ol_part_create (meshA->part[i]);
    ol->olMeshA->part[i] = olp;

    olp->sInitToOlFace   = 0;
    olp->initToOlFaceIdx = NULL;
    olp->initToOlFace    = NULL;
    olp->nLinkedFace     = 0;
    olp->linkedFaces     = NULL;

    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
    const int *iniface_vtx_idxA =  PDM_surf_mesh_part_face_vtx_idx_get (ol->meshA, i);
    const int *iniface_vtxA =  PDM_surf_mesh_part_face_vtx_get (ol->meshA, i);

    int n_facePart = nUnChangedFacePartA[i] + nSubFacePartA[i];
    olp->sInitToOlFace = n_facePart;
    olp->nLinkedFace = nSubFacePartA[i];

    olp->initToOlFaceIdx = malloc (sizeof(int) * (n_face + 1));
    olp->initToOlFaceIdx[0] = 0;

    olp->initToOlFace = malloc (sizeof(int) * n_facePart);
    olp->linkedFaces =  malloc (sizeof(int) * 4 * nSubFacePartA[i]);

    int *face_vtx_idxPart = malloc (sizeof(int) * (n_facePart + 1));
    face_vtx_idxPart[0] = 0;

    PDM_g_num_t *gface_vtxPart = malloc (sizeof(PDM_g_num_t) * s_olface_vtxA[i]);
    int *face_vtxPart = malloc(sizeof(int) * s_olface_vtxA[i]);
    PDM_g_num_t *face_ln_to_gn = malloc (sizeof(PDM_g_num_t) * n_facePart);

    int n_facePartStored = 0;

    idx1 = 0;
    int *initToOlTmp1 = malloc (sizeof(int) * n_facePart);
    int *initToOlTmp2 = malloc (sizeof(int) * n_facePart);
    int *initToOlTmpN = malloc (sizeof(int) * n_face);

    int itmp1 = 0;
    for (int j = 0; j < n_face; j++) {
      initToOlTmpN[j] = 0;
      if (firstRecvStrideA[i][j] == 0) {
        initToOlTmpN[j] = 1;
        face_vtx_idxPart[n_facePartStored+1] =
                face_vtx_idxPart[n_facePartStored] +
                (iniface_vtx_idxA[j+1] - iniface_vtx_idxA[j]);
        for (int k = iniface_vtx_idxA[j]; k < iniface_vtx_idxA[j+1]; k++) {
          face_vtxPart[idx1++] = iniface_vtxA[k];
        }
        initToOlTmp1[n_facePartStored] = n_facePartStored+1;
        face_ln_to_gn[n_facePartStored] = nTSubFacesA
                                    + begUnChangedFaceA
                                    + idxIpart
                                    + 1;
        idxIpart += 1;
        n_facePartStored += 1;
      }
    }

    const PDM_g_num_t *iniGNumVtxA = PDM_surf_mesh_part_vtx_g_num_get (ol->meshA, i);
    const int ini_n_vtx = PDM_surf_mesh_part_n_vtx_get (ol->meshA, i);

    PDM_g_num_t *cpIniGNumVtxA = malloc (sizeof(PDM_g_num_t) * ini_n_vtx);
    int *orderVtxA = malloc (sizeof(int) * ini_n_vtx);
    for (int j = 0; j < ini_n_vtx; j++) {
      cpIniGNumVtxA[j] = iniGNumVtxA[j];
      orderVtxA[j] = j+1;
    }

    PDM_sort_long (cpIniGNumVtxA, orderVtxA, ini_n_vtx);

    idx  = 0;
    idx2 = 0;
    idx3 = 0;
    idx4 = 0;
    idx5 = 0;
    int itmp2 = 0;
    int n_linked_faces = 0;

    for (int j = 0; j < n_face; j++) {

      if (firstRecvStrideA[i][j] > 0) {
        initToOlTmpN[j] = firstRecvStrideA[i][j]/5;
        for (int k = 0; k < firstRecvStrideA[i][j]/5; k++) {
          int n_vtx    = (int) firstRecvA[i][idx++];
          PDM_g_num_t numAbs  = firstRecvA[i][idx++];
          idx++;
          int iProcB  = (int) firstRecvA[i][idx++];
          int iPartB  = (int) firstRecvA[i][idx++];
          face_vtx_idxPart[n_facePartStored+1] =
                  face_vtx_idxPart[n_facePartStored] + n_vtx;

          face_ln_to_gn[n_facePartStored] = numAbs;

          if (iProcB >= 0) {
            olp->linkedFaces[4*n_linked_faces]   = n_facePartStored + 1;
            olp->linkedFaces[4*n_linked_faces+1] = iProcB;
            olp->linkedFaces[4*n_linked_faces+2] = iPartB;
            olp->linkedFaces[4*n_linked_faces+3] = -1; // Come from B in the next step
            n_linked_faces++;
          }

          initToOlTmp2[itmp2++] = n_facePartStored + 1;
          n_facePartStored += 1;

          for (int k1 = 0; k1 < n_vtx; k1++) {
            PDM_g_num_t currentVtx = secondRecvA[i][idx2++];
            if (currentVtx <= iniNGVtx) {
              int idxVtx = PDM_binary_search_long (currentVtx,
                                                   cpIniGNumVtxA,
                                                   ini_n_vtx);
              assert (idxVtx != -1);
              face_vtxPart[idx1++] = orderVtxA[idxVtx];
            }
            else {
              gface_vtxPart[idx3++] = currentVtx;
              face_vtxPart[idx1++] = -idx3;
            }

          }
        }
      }
    }

    free (orderVtxA);
    free (cpIniGNumVtxA);

    olp->nLinkedFace = n_linked_faces;
    olp->linkedFaces = realloc (olp->linkedFaces, sizeof(int) * 4 * n_linked_faces);

    for (int j = 0; j < n_face; j++) {
      olp->initToOlFaceIdx[j+1] = olp->initToOlFaceIdx[j] + initToOlTmpN[j];
      initToOlTmpN[j] = 0;
    }

    itmp1 = 0;
    itmp2 = 0;
    for (int j = 0; j < n_face; j++) {
      int ii = olp->initToOlFaceIdx[j] + initToOlTmpN[j];
      if (firstRecvStrideA[i][j] == 0) {
        olp->initToOlFace[ii] = initToOlTmp1[itmp1++];
      }
      else {
        for (int k = 0; k < firstRecvStrideA[i][j]/5; k++) {
                               olp->initToOlFace[ii++] = initToOlTmp2[itmp2++];
        }
      }
      initToOlTmpN[j]++;
    }

    free (initToOlTmp1);
    free (initToOlTmp2);
    free (initToOlTmpN);

    PDM_g_num_t *cpgface_vtxPart = malloc (sizeof(PDM_g_num_t) * idx3);
    int *ordergface_vtxPart = malloc (sizeof(int) * idx3);
    for (int j = 0; j < idx3; j++) {
      cpgface_vtxPart[j] = gface_vtxPart[j];
      ordergface_vtxPart [j] = j;
    }

    PDM_sort_long (cpgface_vtxPart, ordergface_vtxPart, idx3);

    int k1 = 0;
    int k2 = 0;

    while (idx3 > 0) {
      PDM_g_num_t val = cpgface_vtxPart[k1];
      cpgface_vtxPart[k2] = cpgface_vtxPart[k1];
      k2 += 1;
      while (cpgface_vtxPart[k1] == val) {
        k1 += 1;
        if (k1 >= idx3) {
          break;
        }
      }
      if (k1 >= idx3) {
        break;
      }
    }

    int inin_vtx = PDM_surf_mesh_part_n_vtx_get (ol->meshA, i);
    for (int j = 0; j < s_olface_vtxA[i]; j++) {
      int iniIdx = face_vtxPart[j];
      if (iniIdx < 0) {
        PDM_g_num_t iniVal = gface_vtxPart[-(iniIdx+1)];
        int idxVtx = PDM_binary_search_long (iniVal, cpgface_vtxPart, k2);
        assert(idxVtx != -1);
        face_vtxPart[j] = inin_vtx + idxVtx + 1;
      }
    }

    free (ordergface_vtxPart);

    int n_vtxPart = inin_vtx + k2;

    double *coordsPart = malloc (sizeof(double) * 3 * n_vtxPart);
    PDM_g_num_t *vtx_ln_to_gn_part = malloc (sizeof(PDM_g_num_t) * n_vtxPart);

    for (int j = 0; j < inin_vtx; j++) {
      vtx_ln_to_gn_part[j] = iniGNumVtxA[j];
    }

    for (int j = 0; j < k2; j++) {
      vtx_ln_to_gn_part[inin_vtx + j] = cpgface_vtxPart[j];
    }
    free (cpgface_vtxPart);
    free (gface_vtxPart);

    PDM_g_num_t *cpvtx_ln_to_gn_part = malloc (sizeof(PDM_g_num_t) * n_vtxPart);
    int *ordervtx_ln_to_gn_part = malloc (sizeof(int) * n_vtxPart);
    for (int j = 0; j < n_vtxPart; j++) {
      cpvtx_ln_to_gn_part[j] = vtx_ln_to_gn_part[j];
      ordervtx_ln_to_gn_part[j] = j;
    }

    PDM_sort_long (cpvtx_ln_to_gn_part, ordervtx_ln_to_gn_part, n_vtxPart);

    const double *iniVtx = PDM_surf_mesh_part_vtx_get (ol->meshA, i);

    for (int j = 0; j < inin_vtx; j++) {
      for (int k = 0; k < 3; k++) {
        coordsPart[3*j+k] = iniVtx[3*j+k];
      }
    }

    idx  = 0;
    idx2 = 0;
    idx3 = 0;
    idx4 = 0;
    idx5 = 0;
    for (int j = 0; j < n_face; j++) {
      if (firstRecvStrideA[i][j] > 0) {
        for (int k = 0; k < firstRecvStrideA[i][j]/5; k++) {
          PDM_g_num_t n_vtx    = firstRecvA[i][idx++];
          idx += 4;

          for (int k4 = 0; k4 < n_vtx; k4++) {
            PDM_g_num_t currentVtx = secondRecvA[i][idx2++];
            int idxVtx = PDM_binary_search_long (currentVtx,
                                                 cpvtx_ln_to_gn_part,
                                                 n_vtxPart);

            coordsPart[3*ordervtx_ln_to_gn_part[idxVtx]    ] = thirdRecvA[i][idx3++];
            coordsPart[3*ordervtx_ln_to_gn_part[idxVtx] + 1] = thirdRecvA[i][idx3++];
            coordsPart[3*ordervtx_ln_to_gn_part[idxVtx] + 2] = thirdRecvA[i][idx3++];
          }
        }
      }
    }

    olp->part = PDM_surf_part_create (n_facePart,
                                      face_vtx_idxPart,
                                      face_vtxPart,
                                      face_ln_to_gn,
                                      n_vtxPart,
                                      coordsPart,
                                      vtx_ln_to_gn_part);

    olp->faceIniVtxIdx = malloc (sizeof(int) * (n_face + 1));

    olp->faceIniVtxIdx[0] = 0;
    for (int j = 0; j < n_face; j++) {
      if (fourthRecvStrideA[i][j] > 0) {
        olp->faceIniVtxIdx[j+1] = olp->faceIniVtxIdx[j] + fourthRecvStrideA[i][j];
      }
      else {
        olp->faceIniVtxIdx[j+1] = olp->faceIniVtxIdx[j] +
                iniface_vtx_idxA[j+1] - iniface_vtx_idxA[j];

      }
    }

    olp->sFaceIniVtx   = olp->faceIniVtxIdx[n_face];
    olp->faceIniVtx    = malloc (sizeof(int) * olp->faceIniVtxIdx[n_face]);

    idx = 0;
    idx2 = 0;
    for (int j = 0; j < n_face; j++) {
      if (fourthRecvStrideA[i][j] > 0) {
        for (int k = 0; k < fourthRecvStrideA[i][j]; k++) {
          PDM_g_num_t currentVtx = fourthRecvA[i][idx2++];
          int idxVtx = PDM_binary_search_long (currentVtx,
                                               cpvtx_ln_to_gn_part,
                                               n_vtxPart);
          assert(idxVtx != -1);
          olp->faceIniVtx[idx++] = ordervtx_ln_to_gn_part[idxVtx];
        }
      }
      else {
        for (int k = iniface_vtx_idxA[j]; k < iniface_vtx_idxA[j+1]; k++) {
          olp->faceIniVtx[idx++] = iniface_vtxA[k];

        }
      }
    }

    free (ordervtx_ln_to_gn_part);
    free (cpvtx_ln_to_gn_part);
  }

  if (1 == 0) {
    printf("*** Maillage resultat A : \n");
    for (int i = 0; i < n_partA; i++) {
      printf("\n  ******  Partition : %d\n", i);
      PDM_surf_part_dump (ol->olMeshA->part[i]->part);
      int n_face = PDM_surf_mesh_part_n_face_get (ol->meshA, i);
      printf(" *** faceIniVtx :\n");
      for (int j = 0; j < n_face; j++) {
        for (int k = ol->olMeshA->part[i]->faceIniVtxIdx[j]; k < ol->olMeshA->part[i]->faceIniVtxIdx[j+1]; k++) {
          printf(" %d", ol->olMeshA->part[i]->faceIniVtx[k]);
        }
        printf("\n");
      }
      printf(" ***  initToOlFaceIdx:\n");
      for (int j = 0; j < n_face; j++) {
        for (int k = ol->olMeshA->part[i]->initToOlFaceIdx[j]; k < ol->olMeshA->part[i]->initToOlFaceIdx[j+1]; k++) {
          printf(" %d", ol->olMeshA->part[i]->initToOlFace[k]);
        }
        printf("\n");
      }

    }
  }


  for (int i = 0; i < n_partA; i++) {
    free (firstRecvA[i]);
    free (firstRecvStrideA[i]);
  }

  free (firstRecvA);
  free (firstRecvStrideA);

  for (int i = 0; i < n_partA; i++) {
    free (secondRecvA[i]);
    free (secondRecvStrideA[i]);
  }

  free (secondRecvA);
  free (secondRecvStrideA);

  for (int i = 0; i < n_partA; i++) {
    free (thirdRecvA[i]);
    free (thirdRecvStrideA[i]);
  }

  free (thirdRecvA);
  free (thirdRecvStrideA);

  for (int i = 0; i < n_partA; i++) {
    free (fourthRecvA[i]);
    free (fourthRecvStrideA[i]);
  }
  free (fourthRecvA);
  free (fourthRecvStrideA);

  free (nUnChangedFacePartA);
  free (nSubFacePartA);
  free (s_olface_vtxA);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTA] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTA]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_COMPUTE_LOCAL_CONNECTA]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_COMPUTE_LOCAL_CONNECTA]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);
  PDM_MPI_Barrier (ol->comm);
  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_compute_local_connecta !!!!\n");
    fflush(stdout);
  }

  /* PDM_MPI_Barrier (ol->comm); */
  /* printf ("!!!! Arret temporaire apres construction locale de la connectivite de A !!!!\n"); */
  /* if (i_rank == 0) { */
  /*   _ol_dump_times (ol); */
  /* } */
  /* exit(0); */

  /*****************************************************************************
   *                                                                           *
   * Compute B local numbering for vertices and faces                          *
   *                                                                           *
   ****************************************************************************/

  PDM_g_num_t nUnChangedFaceB = 0;
  PDM_g_num_t nSubFaceB = 0;
  int *nUnChangedFacePartB = malloc(sizeof(int) * n_partB);
  int *nSubFacePartB = malloc(sizeof(int) * n_partB);
  int *s_olface_vtxB = PDM_array_zeros_int(n_partB);

  /*
   * Compute dimensions
   */

  for (int i = 0; i < n_partB; i++) {
    nUnChangedFacePartB[i] = 0;
    nSubFacePartB[i] = 0;
    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    const PDM_g_num_t *gNumFaceB = PDM_surf_mesh_part_face_g_num_get (ol->meshB, i);
    const int *iniface_vtx_idxB =  PDM_surf_mesh_part_face_vtx_idx_get (ol->meshB, i);

    idx = 0;
    for (int j = 0; j < n_face; j++) {
      if (firstRecvStrideB[i][j] == 0) {
        nUnChangedFacePartB[i] += 1;
        s_olface_vtxB[i] += iniface_vtx_idxB[j+1] - iniface_vtx_idxB[j];
      }
      else {
        nSubFacePartB[i] += firstRecvStrideB[i][j]/5;
        for (int k = 0; k < firstRecvStrideB[i][j]/5; k++) {
          int n_vtx    = (int) firstRecvB[i][idx++];
          idx++;
          PDM_g_num_t oNumAbs = firstRecvB[i][idx++];
          idx += 2;

          s_olface_vtxB[i] += n_vtx;

          assert(oNumAbs == gNumFaceB[j]);
        }
      }
    }

    nUnChangedFaceB += nUnChangedFacePartB[i];
    nSubFaceB += nSubFacePartB[i];

  }

  PDM_g_num_t begUnChangedFaceB = 0;
  PDM_MPI_Scan (&nUnChangedFaceB, &begUnChangedFaceB,
            1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  begUnChangedFaceB += -nUnChangedFaceB + 1;

  PDM_g_num_t n_faceB = nUnChangedFaceB + nSubFaceB;
  PDM_g_num_t nGFaceB;
  PDM_MPI_Allreduce(&n_faceB, &nGFaceB, 1,
                PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ol->comm);

  PDM_g_num_t nGVtxB = n_g_newVtxB;
  ol->olMeshB = _ol_mesh_create (nGFaceB, nGVtxB, n_partB);

  /*
   * Build local connectivity
   */

  idxIpart = 0;

  iniNGVtx = PDM_surf_mesh_n_g_vtx_get (ol->meshB);

  for (int i = 0; i < n_partB; i++) {

    _ol_part_t *olp = _ol_part_create (meshB->part[i]);
    ol->olMeshB->part[i] = olp;

    olp->sInitToOlFace   = 0;
    olp->initToOlFaceIdx = NULL;
    olp->initToOlFace    = NULL;
    olp->nLinkedFace     = 0;
    olp->linkedFaces     = NULL;

    int n_face = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
    const int *iniface_vtx_idxB =  PDM_surf_mesh_part_face_vtx_idx_get (ol->meshB, i);
    const int *iniface_vtxB =  PDM_surf_mesh_part_face_vtx_get (ol->meshB, i);

    int n_facePart = nUnChangedFacePartB[i] + nSubFacePartB[i];
    olp->sInitToOlFace = n_facePart;
    olp->nLinkedFace = nSubFacePartB[i];

    olp->initToOlFaceIdx = malloc (sizeof(int) * (n_face + 1));
    olp->initToOlFaceIdx[0] = 0;

    olp->initToOlFace = malloc (sizeof(int) * n_facePart);
    olp->linkedFaces =  malloc (sizeof(int) * 4 * nSubFacePartB[i]);

    int *face_vtx_idxPart = malloc (sizeof(int) * (n_facePart + 1));
    face_vtx_idxPart[0] = 0;

    PDM_g_num_t *gface_vtxPart = malloc (sizeof(PDM_g_num_t) * s_olface_vtxB[i]);
    int *face_vtxPart = malloc(sizeof(int) * s_olface_vtxB[i]);
    PDM_g_num_t *face_ln_to_gn = malloc (sizeof(PDM_g_num_t) * n_facePart);

    int n_facePartStored = 0;

    idx1 = 0;
    int *initToOlTmp1 = malloc (sizeof(int) * n_facePart);
    int *initToOlTmp2 = malloc (sizeof(int) * n_facePart);
    int *initToOlTmpN = malloc (sizeof(int) * n_face);

    int itmp1 = 0;
    for (int j = 0; j < n_face; j++) {
      initToOlTmpN[j] = 0;
      if (firstRecvStrideB[i][j] == 0) {
        initToOlTmpN[j] = 1;
        face_vtx_idxPart[n_facePartStored+1] =
                face_vtx_idxPart[n_facePartStored] +
                (iniface_vtx_idxB[j+1] - iniface_vtx_idxB[j]);
        for (int k = iniface_vtx_idxB[j]; k < iniface_vtx_idxB[j+1]; k++) {
          face_vtxPart[idx1++] = iniface_vtxB[k];
        }
        initToOlTmp1[n_facePartStored] = n_facePartStored + 1;
        face_ln_to_gn[n_facePartStored] = nTSubFacesB
                                    + begUnChangedFaceB
                                    + idxIpart;

        idxIpart += 1;
        n_facePartStored += 1;
      }
    }

    const PDM_g_num_t *iniGNumVtxB = PDM_surf_mesh_part_vtx_g_num_get (ol->meshB, i);
    const int ini_n_vtx = PDM_surf_mesh_part_n_vtx_get (ol->meshB, i);

    PDM_g_num_t *cpIniGNumVtxB = malloc (sizeof(PDM_g_num_t) * ini_n_vtx);
    int *orderVtxB = malloc (sizeof(int) * ini_n_vtx);
    for (int j = 0; j < ini_n_vtx; j++) {
      cpIniGNumVtxB[j] = iniGNumVtxB[j];
      orderVtxB[j] = j + 1;
    }

    PDM_sort_long (cpIniGNumVtxB, orderVtxB, ini_n_vtx);

    idx  = 0;
    idx2 = 0;
    idx3 = 0;
    idx4 = 0;
    idx5 = 0;
    int itmp2 = 0;
    int n_linked_faces = 0;
    for (int j = 0; j < n_face; j++) {
      if (firstRecvStrideB[i][j] > 0) {
        initToOlTmpN[j] = firstRecvStrideB[i][j]/5;
        for (int k = 0; k < firstRecvStrideB[i][j]/5; k++) {
          int n_vtx    = (int) firstRecvB[i][idx++];
          PDM_g_num_t numAbs  = firstRecvB[i][idx++];
          idx++;
          int iProcA  = (int) firstRecvB[i][idx++];
          int iPartA  = (int) firstRecvB[i][idx++];
          face_vtx_idxPart[n_facePartStored+1] =
                  face_vtx_idxPart[n_facePartStored] + n_vtx;

          face_ln_to_gn[n_facePartStored] = numAbs;

          if (iProcA >= 0) {
            olp->linkedFaces[4*n_linked_faces]   = n_facePartStored + 1;
            olp->linkedFaces[4*n_linked_faces+1] = iProcA;
            olp->linkedFaces[4*n_linked_faces+2] = iPartA;
            olp->linkedFaces[4*n_linked_faces+3] = -1; // Come from A in the next step
            n_linked_faces++;
          }

          initToOlTmp2[itmp2++] = n_facePartStored + 1;
          n_facePartStored += 1;

          for (int k1 = 0; k1 < n_vtx; k1++) {
            PDM_g_num_t currentVtx = secondRecvB[i][idx2++];
            if (currentVtx <= iniNGVtx) {
              int idxVtx = PDM_binary_search_long (currentVtx,
                                                   cpIniGNumVtxB,
                                                   ini_n_vtx);
              assert (idxVtx != -1);
              face_vtxPart[idx1++] = orderVtxB[idxVtx];
            }
            else {
              gface_vtxPart[idx3++] = currentVtx;
              face_vtxPart[idx1++] = -idx3;
            }

          }
        }
      }
    }

    free (orderVtxB);
    free (cpIniGNumVtxB);

    olp->nLinkedFace = n_linked_faces;
    olp->linkedFaces = realloc (olp->linkedFaces, sizeof(int) * 4 *n_linked_faces);

    for (int j = 0; j < n_face; j++) {
      olp->initToOlFaceIdx[j+1] = olp->initToOlFaceIdx[j] + initToOlTmpN[j];
      initToOlTmpN[j] = 0;
    }

    itmp1 = 0;
    itmp2 = 0;
    for (int j = 0; j < n_face; j++) {
      int ii = olp->initToOlFaceIdx[j] + initToOlTmpN[j];
      if (firstRecvStrideB[i][j] == 0) {
        olp->initToOlFace[ii] = initToOlTmp1[itmp1++];
      }
      else {
        for (int k = 0; k < firstRecvStrideB[i][j]/5; k++) {
          olp->initToOlFace[ii++] = initToOlTmp2[itmp2++];
        }
      }
      initToOlTmpN[j]++;
    }

    free (initToOlTmp1);
    free (initToOlTmp2);
    free (initToOlTmpN);

    PDM_g_num_t *cpgface_vtxPart = malloc (sizeof(PDM_g_num_t) * idx3);
    int *ordergface_vtxPart = malloc (sizeof(int) * idx3);
    for (int j = 0; j < idx3; j++) {
      cpgface_vtxPart[j] = gface_vtxPart[j];
      ordergface_vtxPart [j] = j;
    }

    PDM_sort_long (cpgface_vtxPart, ordergface_vtxPart, idx3);

    int k1 = 0;
    int k2 = 0;
    while (idx3 > 0) {
      PDM_g_num_t val = cpgface_vtxPart[k1];
      cpgface_vtxPart[k2] = cpgface_vtxPart[k1];
      k2 += 1;
      while (cpgface_vtxPart[k1] == val) {
        k1 += 1;
        if (k1 >= idx3) {
          break;
        }
      }
      if (k1 >= idx3) {
        break;
      }
    }

    int inin_vtx = PDM_surf_mesh_part_n_vtx_get (ol->meshB, i);
    for (int j = 0; j < s_olface_vtxB[i]; j++) {
      int iniIdx = face_vtxPart[j];
      if (iniIdx < 0) {
        PDM_g_num_t iniVal = gface_vtxPart[-(iniIdx+1)];
        int idxVtx = PDM_binary_search_long (iniVal, cpgface_vtxPart, k2);
        assert(idxVtx != -1);
        face_vtxPart[j] = inin_vtx + idxVtx + 1;
      }
    }

    free (ordergface_vtxPart);

    int n_vtxPart = inin_vtx + k2;

    double *coordsPart = malloc (sizeof(double) * 3 * n_vtxPart);
    PDM_g_num_t *vtx_ln_to_gn_part = malloc (sizeof(PDM_g_num_t) * n_vtxPart);

    for (int j = 0; j < inin_vtx; j++) {
      vtx_ln_to_gn_part[j] = iniGNumVtxB[j];
    }

    for (int j = 0; j < k2; j++) {
      vtx_ln_to_gn_part[inin_vtx + j] = cpgface_vtxPart[j];
    }

    free (cpgface_vtxPart);
    free (gface_vtxPart);

    PDM_g_num_t *cpvtx_ln_to_gn_part = malloc (sizeof(PDM_g_num_t) * n_vtxPart);
    int *ordervtx_ln_to_gn_part = malloc (sizeof(int) * n_vtxPart);
    for (int j = 0; j < n_vtxPart; j++) {
      cpvtx_ln_to_gn_part[j] = vtx_ln_to_gn_part[j];
      ordervtx_ln_to_gn_part[j] = j;
    }

    PDM_sort_long (cpvtx_ln_to_gn_part, ordervtx_ln_to_gn_part, n_vtxPart);

    const double *iniVtx = PDM_surf_mesh_part_vtx_get (ol->meshB, i);

    for (int j = 0; j < inin_vtx; j++) {
      for (int k = 0; k < 3; k++) {
        coordsPart[3*j+k] = iniVtx[3*j+k];
      }
    }

    idx  = 0;
    idx2 = 0;
    idx3 = 0;
    idx4 = 0;
    idx5 = 0;
    for (int j = 0; j < n_face; j++) {
      if (firstRecvStrideB[i][j] > 0) {
        for (int k = 0; k < firstRecvStrideB[i][j]/5; k++) {
          PDM_g_num_t n_vtx    = firstRecvB[i][idx++];
          idx += 4;

          for (int k4 = 0; k4 < n_vtx; k4++) {
            PDM_g_num_t currentVtx = secondRecvB[i][idx2++];
            int idxVtx = PDM_binary_search_long (currentVtx,
                                                 cpvtx_ln_to_gn_part,
                                                 n_vtxPart);

            coordsPart[3*ordervtx_ln_to_gn_part[idxVtx]    ] = thirdRecvB[i][idx3++];
            coordsPart[3*ordervtx_ln_to_gn_part[idxVtx] + 1] = thirdRecvB[i][idx3++];
            coordsPart[3*ordervtx_ln_to_gn_part[idxVtx] + 2] = thirdRecvB[i][idx3++];

          }
        }
      }
    }

    olp->part = PDM_surf_part_create (n_facePart,
                                      face_vtx_idxPart,
                                      face_vtxPart,
                                      face_ln_to_gn,
                                      n_vtxPart,
                                      coordsPart,
                                      vtx_ln_to_gn_part);


    olp->faceIniVtxIdx = malloc (sizeof(int) * (n_face + 1));

    olp->faceIniVtxIdx[0] = 0;
    for (int j = 0; j < n_face; j++) {
      if (fourthRecvStrideB[i][j] > 0) {
        olp->faceIniVtxIdx[j+1] = olp->faceIniVtxIdx[j] + fourthRecvStrideB[i][j];
      }
      else {
        olp->faceIniVtxIdx[j+1] = olp->faceIniVtxIdx[j] +
                iniface_vtx_idxB[j+1] - iniface_vtx_idxB[j];

      }
    }

    olp->sFaceIniVtx   = olp->faceIniVtxIdx[n_face];
    olp->faceIniVtx    = malloc (sizeof(int) * olp->faceIniVtxIdx[n_face]);

    idx = 0;
    idx2 = 0;
    for (int j = 0; j < n_face; j++) {
      if (fourthRecvStrideB[i][j] > 0) {
        for (int k = 0; k < fourthRecvStrideB[i][j]; k++) {
          PDM_g_num_t currentVtx = fourthRecvB[i][idx2++];
          int idxVtx = PDM_binary_search_long (currentVtx,
                                               cpvtx_ln_to_gn_part,
                                               n_vtxPart);

          olp->faceIniVtx[idx++] = ordervtx_ln_to_gn_part[idxVtx];
        }
      }
      else {
        for (int k = iniface_vtx_idxB[j]; k < iniface_vtx_idxB[j+1]; k++) {
          olp->faceIniVtx[idx++] = iniface_vtxB[k];

        }
      }
    }

    free (ordervtx_ln_to_gn_part);
    free (cpvtx_ln_to_gn_part);

  }

  if (1 == 0) {
    printf("*** Maillage resultat B : \n");
    for (int i = 0; i < n_partB; i++) {
      printf("\n  ******  Partition : %d\n", i);
      PDM_surf_part_dump (ol->olMeshB->part[i]->part);
      int n_face = PDM_surf_mesh_part_n_face_get (ol->meshB, i);
      printf(" *** faceIniVtx :\n");
      for (int j = 0; j < n_face; j++) {
        for (int k = ol->olMeshB->part[i]->faceIniVtxIdx[j]; k < ol->olMeshB->part[i]->faceIniVtxIdx[j+1]; k++) {
          printf(" %d", ol->olMeshB->part[i]->faceIniVtx[k]);
        }
        printf("\n");
      }
      printf(" ***  initToOlFaceIdx:\n");
      for (int j = 0; j < n_face; j++) {
        for (int k = ol->olMeshB->part[i]->initToOlFaceIdx[j]; k < ol->olMeshB->part[i]->initToOlFaceIdx[j+1]; k++) {
          printf(" %d", ol->olMeshB->part[i]->initToOlFace[k]);
        }
        printf("\n");
      }
    }
  }

  for (int i = 0; i < n_partB; i++) {
    free (firstRecvB[i]);
    free (firstRecvStrideB[i]);
  }
  free (firstRecvB);
  free (firstRecvStrideB);

  for (int i = 0; i < n_partB; i++) {
    free (secondRecvB[i]);
    free (secondRecvStrideB[i]);
  }
  free (secondRecvB);
  free (secondRecvStrideB);

  for (int i = 0; i < n_partB; i++) {
    free (thirdRecvB[i]);
    free (thirdRecvStrideB[i]);
  }
  free (thirdRecvB);
  free (thirdRecvStrideB);

  for (int i = 0; i < n_partB; i++) {
    free (fourthRecvB[i]);
    free (fourthRecvStrideB[i]);
  }
  free (fourthRecvB);
  free (fourthRecvStrideB);

  free (nUnChangedFacePartB);
  free (nSubFacePartB);
  free (s_olface_vtxB);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_COMPUTE_LOCAL_CONNECTB] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_COMPUTE_LOCAL_CONNECTB]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_COMPUTE_LOCAL_CONNECTB]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_COMPUTE_LOCAL_CONNECTB]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  PDM_MPI_Barrier (ol->comm);
  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_compute_local_connectb !!!!\n");
    fflush(stdout);
  }
  /* PDM_MPI_Barrier (ol->comm); */
  /* printf ("!!!! Arret temporaire apres construction locale de la connectivite de B !!!!\n"); */
  /* if (i_rank == 0) { */
  /*   _ol_dump_times (ol); */
  /* } */
  /* exit(0); */

  /*****************************************************************************
   *                                                                           *
   * Update boundary graph with local face number (link)                       *
   *                                                                           *
   ****************************************************************************/

  int *sendIdx = malloc (sizeof(int)*lComm);
  int *sendN   = PDM_array_zeros_int(lComm);
  PDM_g_num_t *sendBuff = NULL;

  int *recvIdx = malloc (sizeof(int)*lComm);
  int *recvN   = malloc (sizeof(int)*lComm);
  PDM_g_num_t *recvBuff = NULL;


  int n_t_send = 0;
  for (int i = 0; i < n_partA; i++) {

    _ol_part_t *olp = ol->olMeshA->part[i];
    n_t_send += olp->nLinkedFace;

    for (int j = 0; j < olp->nLinkedFace; j++) {
      int iProc    = olp->linkedFaces[4*j + 1];
      sendN[iProc] += 1;
    }

  }

  int n_t_recv = 0;
  for (int i = 0; i < n_partB; i++) {

    _ol_part_t *olp = ol->olMeshB->part[i];
    n_t_recv += olp->nLinkedFace;
  }

  sendIdx[0] = 0;
  for (int i = 1; i < lComm; i++) {
    sendIdx[i] = sendIdx[i-1] + 3 * sendN[i-1];
    sendN[i-1] = 0;
  }
  sendN[lComm-1] = 0;

  sendBuff = malloc(sizeof(PDM_g_num_t) * n_t_send * 3);
  recvBuff = malloc(sizeof(PDM_g_num_t) * n_t_recv * 3);

  for (int i = 0; i < n_partA; i++) {

    _ol_part_t *olp = ol->olMeshA->part[i];

    const PDM_g_num_t *numAbs = PDM_surf_part_faceLnToGn_get (olp->part);

    for (int j = 0; j < olp->nLinkedFace; j++) {
      int iFacLoc  = olp->linkedFaces[4*j    ];
      int iProc    = olp->linkedFaces[4*j + 1];
      int i_part    = olp->linkedFaces[4*j + 2];

      int id = sendIdx[iProc] + sendN[iProc];
      sendN[iProc] += 3;

      sendBuff[id    ] = i_part;
      sendBuff[id + 1] = numAbs[iFacLoc-1];
      sendBuff[id + 2] = iFacLoc;

    }

  }

  PDM_MPI_Alltoall (sendN, 1, PDM_MPI_INT, recvN, 1, PDM_MPI_INT, ol->comm);

  recvIdx[0] = 0;
  for (int i = 1; i < lComm; i++) {
    recvIdx[i] = recvIdx[i-1] + recvN[i-1];
  }

  PDM_MPI_Alltoallv (sendBuff, sendN, sendIdx, PDM__PDM_MPI_G_NUM,
                     recvBuff, recvN, recvIdx, PDM__PDM_MPI_G_NUM, ol->comm);

  PDM_g_num_t **face_ln_to_gnBSorted = malloc (sizeof(PDM_g_num_t *) * n_partB);
  int **face_ln_to_gnBOrder = malloc (sizeof(int *) * n_partB);
  int **faceToLinked = malloc (sizeof(int *) * n_partB);

  for (int i = 0; i < n_partB; i++) {

    _ol_part_t *olp = ol->olMeshB->part[i];
    PDM_surf_part_t *_part = olp->part;

    face_ln_to_gnBSorted[i] = malloc (sizeof(PDM_g_num_t) * _part->n_face);
    face_ln_to_gnBOrder[i] = malloc (sizeof(int) * _part->n_face);
    faceToLinked[i] = malloc (sizeof(int) * _part->n_face);

    for (int j = 0; j < _part->n_face; j++) {
      face_ln_to_gnBSorted[i][j] = _part->face_ln_to_gn[j];
      face_ln_to_gnBOrder[i][j]  = j;
      faceToLinked[i][j] = -1;
    }

    for (int j = 0; j < olp->nLinkedFace; j++) {
      int iFacLoc  = olp->linkedFaces[4*j    ];

      faceToLinked[i][iFacLoc-1] = j;
    }

    PDM_sort_long (face_ln_to_gnBSorted[i],
                   face_ln_to_gnBOrder[i],
                   _part->n_face);

  }

  int k = 0;
  for (int i = 0; i < n_t_recv; i++) {

    int i_part    = (int) recvBuff[3*i    ];
    PDM_g_num_t numAbs   = recvBuff[3*i + 1];
    int iFacDist = (int) recvBuff[3*i + 2];

    _ol_part_t *olp = ol->olMeshB->part[i_part];
    PDM_surf_part_t *_part = olp->part;
    int id = PDM_binary_search_long (numAbs,
                                     face_ln_to_gnBSorted[i_part],
                                     _part->n_face);
    assert (id != -1);
    int iFacIdx = face_ln_to_gnBOrder[i_part][id];

    int iLinked = faceToLinked[i_part][iFacIdx];

    assert(iLinked != -1);
    olp->linkedFaces[4*iLinked + 3] = iFacDist;

    recvBuff[k++] = iFacIdx + 1;

  }

  for (int i = 0; i < n_partB; i++) {
    free (face_ln_to_gnBSorted[i]);
    free (face_ln_to_gnBOrder[i]);
    free (faceToLinked[i]);
  }

  free (face_ln_to_gnBSorted);
  free (face_ln_to_gnBOrder);
  free (faceToLinked);

  for (int i = 0; i < lComm; i++) {
    sendN[i]   = sendN[i]/3;
    sendIdx[i] = sendIdx[i]/3;
    recvN[i]   = recvN[i]/3;
    recvIdx[i] = recvIdx[i]/3;
  }

  PDM_MPI_Alltoallv (recvBuff, recvN, recvIdx, PDM__PDM_MPI_G_NUM,
                     sendBuff, sendN, sendIdx, PDM__PDM_MPI_G_NUM,
                     ol->comm);

  k = 0;

  PDM_array_reset_int(sendN, lComm, 0);

  for (int i = 0; i < n_partA; i++) {
    _ol_part_t *olp = ol->olMeshA->part[i];

    for (int j = 0; j < olp->nLinkedFace; j++) {
      int iProc    = olp->linkedFaces[4*j + 1];
      int id = sendIdx[iProc] + sendN[iProc]++;
      olp->linkedFaces[4*j + 3] = (int) sendBuff[id];
    }
  }

  free (sendN);
  free (sendIdx);
  free (sendBuff);

  free (recvN);
  free (recvIdx);
  free (recvBuff);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_UPDATE_A_B_CONNECT_GRAPH] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_UPDATE_A_B_CONNECT_GRAPH]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_UPDATE_A_B_CONNECT_GRAPH]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_UPDATE_A_B_CONNECT_GRAPH]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  PDM_MPI_Barrier (ol->comm);
  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_update_a_b_connect_graph !!!!\n");
    fflush(stdout);
  }

  /* PDM_MPI_Barrier (ol->comm); */
  /* printf ("!!!! Arret temporaire apres temporaire apres mis a jour du graphe de comm !!!!\n"); */
  /* if (i_rank == 0) { */
  /*   _ol_dump_times (ol); */
  /* } */
  /* exit(0); */


  /*****************************************************************************
   *                                                                           *
   *                     Sort boundary graph                                   *
   *                                                                           *
   ****************************************************************************/

  for (int i = 0; i < n_partA; i++) {

    _ol_part_t *olp = ol->olMeshA->part[i];

    int *sortGraph = malloc (sizeof(int) * olp->nLinkedFace);
    int *orderGraph = malloc (sizeof(int) * olp->nLinkedFace);

    olp->linkedFacesProcIdx = PDM_array_zeros_int(lComm + 1);

    for (int j = 0; j < olp->nLinkedFace; j++) {
      orderGraph[j] = j;
      sortGraph[j]  = olp->linkedFaces[4*j + 1];
      olp->linkedFacesProcIdx[sortGraph[j]+1]++;
    }

    for (int j = 0; j < lComm; j++) {
      olp->linkedFacesProcIdx[j+1] += olp->linkedFacesProcIdx[j];
    }

    PDM_sort_int (sortGraph, orderGraph, olp->nLinkedFace);

    int *tmpLinkedFaces = malloc (sizeof(int) * 4 * olp->nLinkedFace);

    for (int j = 0; j < olp->nLinkedFace; j++) {
      int newId = orderGraph[j];
      tmpLinkedFaces[4*j    ] = olp->linkedFaces[4*newId    ];
      tmpLinkedFaces[4*j + 1] = olp->linkedFaces[4*newId + 1];
      tmpLinkedFaces[4*j + 2] = olp->linkedFaces[4*newId + 2];
      tmpLinkedFaces[4*j + 3] = olp->linkedFaces[4*newId + 3];
    }

    free (olp->linkedFaces);
    olp->linkedFaces = tmpLinkedFaces;

    if (1 == 0) {
      printf ("*** Mesh A linked face\n");
      for (int j = 0; j < olp->nLinkedFace; j++) {
        printf ("%d : %d %d %d\n",
                olp->linkedFaces[4*j    ],
                olp->linkedFaces[4*j + 1],
                olp->linkedFaces[4*j + 2],
                olp->linkedFaces[4*j + 3]);
      }
    }
    free (sortGraph);
    free (orderGraph);

  }


  for (int i = 0; i < n_partB; i++) {

    _ol_part_t *olp = ol->olMeshB->part[i];


    int *sortGraph = malloc (sizeof(int) * olp->nLinkedFace);
    int *orderGraph = malloc (sizeof(int) * olp->nLinkedFace);

    olp->linkedFacesProcIdx = PDM_array_zeros_int(lComm + 1);

    for (int j = 0; j < olp->nLinkedFace; j++) {
      orderGraph[j] = j;
      sortGraph[j]  = olp->linkedFaces[4*j + 1];
      olp->linkedFacesProcIdx[sortGraph[j]+1]++;
    }

    for (int j = 0; j < lComm; j++) {
      olp->linkedFacesProcIdx[j+1] += olp->linkedFacesProcIdx[j];
    }

    PDM_sort_int (sortGraph, orderGraph, olp->nLinkedFace);

    int *tmpLinkedFaces = malloc (sizeof(int) * 4 * olp->nLinkedFace);

    if (1 == 0) {
      printf ("*** Mesh B linked face 0 \n");
      for (int j = 0; j < olp->nLinkedFace; j++) {
        printf ("%d : %d %d %d\n",
                olp->linkedFaces[4*j    ],
                olp->linkedFaces[4*j + 1],
                olp->linkedFaces[4*j + 2],
                olp->linkedFaces[4*j + 3]);
      }
    }

    for (int j = 0; j < olp->nLinkedFace; j++) {
      int newId = orderGraph[j];
      tmpLinkedFaces[4*j    ] = olp->linkedFaces[4*newId    ];
      tmpLinkedFaces[4*j + 1] = olp->linkedFaces[4*newId + 1];
      tmpLinkedFaces[4*j + 2] = olp->linkedFaces[4*newId + 2];
      tmpLinkedFaces[4*j + 3] = olp->linkedFaces[4*newId + 3];

    }

    free (olp->linkedFaces);
    olp->linkedFaces = tmpLinkedFaces;

    if (1 == 0) {
      printf ("*** Mesh B linked face\n");
      for (int j = 0; j < olp->nLinkedFace; j++) {
        printf ("%d : %d %d %d\n",
                olp->linkedFaces[4*j    ],
                olp->linkedFaces[4*j + 1],
                olp->linkedFaces[4*j + 2],
                olp->linkedFaces[4*j + 3]);
      }
    }


    free (sortGraph);
    free (orderGraph);

  }

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_SORT_A_B_CONNECT_GRAPH] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_SORT_A_B_CONNECT_GRAPH]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_SORT_A_B_CONNECT_GRAPH]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_SORT_A_B_CONNECT_GRAPH]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  PDM_MPI_Barrier (ol->comm);
  if (i_rank == 0 && vb == 1) {
    printf ("!!!! ol_sort_a_b_connect_graph !!!!\n");
    fflush(stdout);
  }
  /* PDM_MPI_Barrier (ol->comm); */
  /* printf ("!!!! Arret final !!!!\n"); */
  /* if (i_rank == 0) { */
  /*   _ol_dump_times (ol); */
  /* } */
  /* exit(0); */

  /*****************************************************************************
   *                                                                           *
   *                           Finish                                          *
   *                                                                           *
   ****************************************************************************/

}


/**
 * \brief Compute overlay general surfaces
 *
 * This function computes overlay mesh of two surface meshes
 *
 * \param [in]  ol       overlay object
 *
 */

static void
_compute_overlay_surfaces
(
 PDM_ol_t *ol
)
{
  PDM_UNUSED(ol);
  PDM_error(__FILE__, __LINE__, 0, "Error _compute_overlay_surfaces : Not yet implemented\n");
  exit(0);
}


/**
 * \brief Compute overlay mesh
 *
 * This function computes overlayh mesh after plane surface checking
 *
 * \param [in]  ol       overlay object
 *
 */

static void
_compute_overlay
(
 PDM_ol_t *ol
)
{
  /*
   * Compute faces extents
   */

  _compute_faceExtents(ol);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_FACE_EXTENTS] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_FACE_EXTENTS]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_FACE_EXTENTS]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_FACE_EXTENTS]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*
   * Check if initial meshes are plane sufaces ou general surfaces
   */

  double dist = 0.;
  int isPlaneSurfaces = _is_same_plane (ol, &dist);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[OL_IS_PLANE] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[OL_IS_PLANE]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[OL_IS_PLANE]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[OL_IS_PLANE]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  double gMinCarLgthVtxA = PDM_surf_mesh_gMinCarLgthVtx_get (ol->meshA);
  double gMinCarLgthVtxB = PDM_surf_mesh_gMinCarLgthVtx_get (ol->meshB);

  double carLgth = _MIN (gMinCarLgthVtxA, gMinCarLgthVtxB);

  if (dist > carLgth) {
    PDM_error(__FILE__, __LINE__, 0,
              "Warning _compute_overlay : distance (%12.5e) between the two planes may"
              " be too far to overlay its : \n"
              " translate one of them and recompute\n", dist);
  }

  if (isPlaneSurfaces) {
    _compute_overlay_planes(ol);
  }
  else {
    _compute_overlay_surfaces(ol);
  }

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Build and initialize an overlaying object
 *
 * This function builds an initializes an overlaying surface meshes object
 *
 * \param [in]  n_partMeshA   Number of local partitions of the meshA input
 * \param [in]  n_partMeshB   Number of local partitions of the meshB input
 * \param [in]  projectCoeff  Projection coefficient to define the overlay surface projection
 *                            If value == 0, the surface projection is MeshA
 *                            If value == 1, the surface projection is MeshB
 *                            If 0 < value < 1 , the projection surface is an intermediate surface
 * \param [in]  tolerance     Absolute tolerance used to define
 *                            local geometric tolerance (0 < tolerance).
 * \param [in]  comm          MPI communicator.
 *
 * \return      id            Overlay object identifier.
 *
 */

PDM_ol_t *
PDM_ol_create
(
 const int          n_partMeshA,
 const int          n_partMeshB,
 const double       projectCoeff,
 const PDM_MPI_Comm comm
)
{
  /*
   * Allocate structure
   */

  PDM_ol_t *ol = (PDM_ol_t *) malloc(sizeof(PDM_ol_t));

  /*
   * Initialization
   */

  ol->timer                = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    ol->times_elapsed[i] = 0.;
    ol->times_cpu[i] = 0.;
    ol->times_cpu_u[i] = 0.;
    ol->times_cpu_s[i] = 0.;
  }

  PDM_timer_resume(ol->timer);

  ol->projectCoeff         = projectCoeff;
  ol->vtxCarLengthTol      = 1e-3;
  ol->extentsTol           = 1e-3;
  ol->samePlaneTol         = 1e-6;
  ol->comm                 = comm;
  ol->meshA                = PDM_surf_mesh_create(n_partMeshA, comm);
  ol->meshB                = PDM_surf_mesh_create(n_partMeshB, comm);
  ol->olMeshA              = NULL;
  ol->olMeshB              = NULL;
  //ol->dbbtreeA             = NULL;

  PDM_timer_hang_on(ol->timer);

  return ol;
}



/**
 * \brief Set an overlay parameter
 *
 * This function sets en overlay parameter
 *
 * \param [in]  id          PDM_ol identifier
 * \param [in]  parameter   Parameter to define
 * \param [in]  value       Parameter value
 *
 */

void
PDM_ol_parameter_set
(
       PDM_ol_t           *ol,
 const PDM_ol_parameter_t  parameter,
 const double              value
)
{
  PDM_timer_resume(ol->timer);

  switch (parameter) {

  case PDM_OL_CAR_LENGTH_TOL :
    ol->vtxCarLengthTol = value;
    break;

  case PDM_OL_EXTENTS_TOL :
    ol->extentsTol = value;
    break;

  case PDM_OL_SAME_PLANE_TOL :
    ol->samePlaneTol = value;
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_parameter_set :"
            " Unavailable parameter to set\n");
    abort();

  }
  PDM_timer_hang_on(ol->timer);
}

/**
 * \brief Define input meshes properties
 *
 * This function defines the input meshes properties
 *
 * \param [in]  id          ol identifier
 * \param [in]  mesh        Input mesh to define
 *                          (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  i_part       Partition to define
 * \param [in]  n_face       Number of faces
 * \param [in]  face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]  face_vtx_idx  face -> vertex connectivity
 * \param [in]  face_ln_to_gn  Local face numbering to global face numbering
 * \param [in]  n_vtx        Number of vertices
 * \param [in]  coords      Coordinates
 * \param [in]  vtx_ln_to_gn   Local vertex numbering to global vertex numbering
 *
 */

void
PDM_ol_input_mesh_set
(
       PDM_ol_t      *ol,
 const PDM_ol_mesh_t  mesh,
 const int            i_part,
 const int            n_face,
 const int           *face_vtx_idx,
 const int           *face_vtx,
 const PDM_g_num_t    *face_ln_to_gn,
 const int            n_vtx,
 const double        *coords,
 const PDM_g_num_t    *vtx_ln_to_gn
)
{
  PDM_timer_resume(ol->timer);

  PDM_surf_mesh_t *_mesh = NULL;

  // printf("111 : %d %d / %d %d / %d %d %d / %ld // %d %12.5e %12.5e %ld\n",
  //        i_part, n_face, face_vtx_idx[0], face_vtx_idx[1], face_vtx[0], face_vtx[1], face_vtx[2],
  //        face_ln_to_gn[0], n_vtx, coords[0], coords[1], vtx_ln_to_gn[0]);

  switch(mesh) {

  case PDM_OL_MESH_A:
    _mesh = ol->meshA;
    break;

  case PDM_OL_MESH_B:
    _mesh = ol->meshB;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_input_mesh_set :"
            " unknown PDM_ol_mesh_t mesh\n");
    abort();
  }

  if (1 == 0) {
    printf("mesh %d g_num :\n", mesh);
    for (int i = 0; i < n_face; i++) {
      printf(" "PDM_FMT_G_NUM" :", face_ln_to_gn[i]);
      for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
        printf(" "PDM_FMT_G_NUM, vtx_ln_to_gn[face_vtx[j]-1]);
      }
      printf("\n");
    }
  }

  PDM_surf_mesh_part_input (_mesh,
                            i_part,
                            n_face,
                            face_vtx_idx,
                            face_vtx,
                            face_ln_to_gn,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn);
  PDM_timer_hang_on(ol->timer);
}


/**
 * \brief Overlaying the input surface meshes
 *
 * This function overlays the input surface meshes
 *
 * \param [in]  id       ol identifier
 *
 */

void
PDM_ol_compute
(
 PDM_ol_t *ol
)
{

  PDM_timer_resume(ol->timer);

  PDM_MPI_Barrier (ol->comm);
  PDM_timer_hang_on(ol->timer);
  ol->times_elapsed[INIT_DEF_DATA] = PDM_timer_elapsed(ol->timer);
  ol->times_cpu[INIT_DEF_DATA]     = PDM_timer_cpu(ol->timer);
  ol->times_cpu_u[INIT_DEF_DATA]   = PDM_timer_cpu_user(ol->timer);
  ol->times_cpu_s[INIT_DEF_DATA]   = PDM_timer_cpu_sys(ol->timer);
  PDM_timer_resume(ol->timer);

  /*
   * First computation
   */

  if (ol->olMeshA == NULL || ol->olMeshB== NULL) {

    /*
     * Computation of edges and inter partition exchange graph
     */
    _build_edges(ol);

    PDM_MPI_Barrier (ol->comm);
    PDM_timer_hang_on(ol->timer);
    ol->times_elapsed[INIT_BUILD_EDGES] = PDM_timer_elapsed(ol->timer);
    ol->times_cpu[INIT_BUILD_EDGES]     = PDM_timer_cpu(ol->timer);
    ol->times_cpu_u[INIT_BUILD_EDGES]   = PDM_timer_cpu_user(ol->timer);
    ol->times_cpu_s[INIT_BUILD_EDGES]   = PDM_timer_cpu_sys(ol->timer);
    PDM_timer_resume(ol->timer);

    /*
     * Inter partition exchange graph (exchange to vertices)
     */

    _build_exchange_graph(ol);

    PDM_MPI_Barrier (ol->comm);
    PDM_timer_hang_on(ol->timer);
    ol->times_elapsed[INIT_BUILD_EXCH_GRAPH] = PDM_timer_elapsed(ol->timer);
    ol->times_cpu[INIT_BUILD_EXCH_GRAPH]     = PDM_timer_cpu(ol->timer);
    ol->times_cpu_u[INIT_BUILD_EXCH_GRAPH]   = PDM_timer_cpu_user(ol->timer);
    ol->times_cpu_s[INIT_BUILD_EXCH_GRAPH]   = PDM_timer_cpu_sys(ol->timer);
    PDM_timer_resume(ol->timer);

    /*
     * Build ghost faces and edges (To compute caracteristic length
     * and normal to vertices)
     */

    _compute_carLgthVtx(ol);

    PDM_MPI_Barrier (ol->comm);
    PDM_timer_hang_on(ol->timer);
    ol->times_elapsed[INIT_COMPUTE_CAR_LGTH] = PDM_timer_elapsed(ol->timer);
    ol->times_cpu[INIT_COMPUTE_CAR_LGTH]     = PDM_timer_cpu(ol->timer);
    ol->times_cpu_u[INIT_COMPUTE_CAR_LGTH]   = PDM_timer_cpu_user(ol->timer);
    ol->times_cpu_s[INIT_COMPUTE_CAR_LGTH]   = PDM_timer_cpu_sys(ol->timer);
    PDM_timer_resume(ol->timer);

  }

  else {

    ol->times_elapsed[INIT_BUILD_EDGES] = ol->times_elapsed[INIT_DEF_DATA];
    ol->times_cpu[INIT_BUILD_EDGES]     = ol->times_cpu[INIT_DEF_DATA];
    ol->times_cpu_u[INIT_BUILD_EDGES]   = ol->times_cpu_u[INIT_DEF_DATA];
    ol->times_cpu_s[INIT_BUILD_EDGES]   = ol->times_cpu_s[INIT_DEF_DATA];

    ol->times_elapsed[INIT_BUILD_EXCH_GRAPH] = ol->times_elapsed[INIT_DEF_DATA];
    ol->times_cpu[INIT_BUILD_EXCH_GRAPH]     = ol->times_cpu[INIT_DEF_DATA];
    ol->times_cpu_u[INIT_BUILD_EXCH_GRAPH]   = ol->times_cpu_u[INIT_DEF_DATA];
    ol->times_cpu_s[INIT_BUILD_EXCH_GRAPH]   = ol->times_cpu_s[INIT_DEF_DATA];

    ol->times_elapsed[INIT_COMPUTE_CAR_LGTH] = ol->times_elapsed[INIT_DEF_DATA];
    ol->times_cpu[INIT_COMPUTE_CAR_LGTH]     = ol->times_cpu[INIT_DEF_DATA];
    ol->times_cpu_u[INIT_COMPUTE_CAR_LGTH]   = ol->times_cpu_u[INIT_DEF_DATA];
    ol->times_cpu_s[INIT_COMPUTE_CAR_LGTH]   = ol->times_cpu_s[INIT_DEF_DATA];

  }

  /*
   * Compute overlay
   */

  _compute_overlay(ol);

  PDM_timer_hang_on (ol->timer);
}



/**
 * \brief Define the type of a mesh moving
 *
 * This function defines the type of a mesh moving.
 * Only a mesh can move
 *
 * \param [in]  id       PDM_ol identifier
 * \param [in]  mesh     Moving mesh
 * \param [in]  mv       Type of moving
 *
 */

void
PDM_ol_moving_type_set
(
       PDM_ol_t      *ol,
 const PDM_ol_mesh_t  mesh,
 const PDM_ol_mv_t    mv
)
{

  PDM_timer_resume(ol->timer);
  PDM_UNUSED(mesh);
  PDM_UNUSED(mv);
  PDM_error(__FILE__, __LINE__, 0, "PDM ERROR : PDM_ol_moving_type_set not implemented yet");
  exit(1);
  PDM_timer_hang_on(ol->timer);
}


/**
 * \brief Define a translation
 *
 * This function defines a translation for the moving mesh
 *
 * \param [in]  id       PDM_ol identifier
 * \param [in]  vect     Translation vector
 * \param [in]  center   Translation center
 *
 */

void
PDM_ol_translation_set
(
       PDM_ol_t *ol,
 const double   *vect,
 const double   *center
 )
{

  PDM_timer_resume(ol->timer);

  PDM_UNUSED(vect);
  PDM_UNUSED(center);
  PDM_error(__FILE__, __LINE__, 0, "PDM ERROR : PDM_ol_translation_set not implemented yet");
  exit(1);

  PDM_timer_hang_on(ol->timer);
}



/**
 * \brief Define a rotation
 *
 * This function defines a rotation for the moving mesh
 *
 * \param [in]  id        PDM_ol identifier
 * \param [in]  direction Rotation direction
 * \param [in]  center    Rotation center
 * \param [in]  angle      Rotation center (degrees)
 *
 */

void
PDM_ol_rotation_set
(
       PDM_ol_t *ol,
 const double   *direction,
 const double   *center,
 const double    angle
)
{

  PDM_timer_resume(ol->timer);

  PDM_UNUSED(direction);
  PDM_UNUSED(center);
  PDM_UNUSED(angle);
  PDM_error(__FILE__, __LINE__, 0, "PDM ERROR : PDM_ol_rotation_set not implemented yet");
  exit(EXIT_FAILURE);
  PDM_timer_hang_on(ol->timer);
}


/**
 * \brief Return the entitie sizes of the overlay mesh
 *
 * This function returns the entities sizes of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id        PDM_ol identifier
 * \param [in]  mesh      Input mesh
 *                        (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [out] nGOlFace  Global number of faces of the overlay mesh
 * \param [out] nGOlVtx   Global number of vertices of the overlay mesh
 *
 */

void
PDM_ol_mesh_dim_get
(
 const PDM_ol_t      *ol,
 const PDM_ol_mesh_t  mesh,
       PDM_g_num_t   *nGOlFace,
       PDM_g_num_t   *nGOlVtx
)
{
  PDM_timer_resume(ol->timer);

  switch(mesh) {

  case PDM_OL_MESH_A:
    *nGOlFace = ol->olMeshA->nGFace;
    *nGOlVtx  = ol->olMeshA->nGVtx;
    break;

  case PDM_OL_MESH_B:
    *nGOlFace = ol->olMeshB->nGFace;
    *nGOlVtx  = ol->olMeshB->nGVtx;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_mesh_dim_get : unknown PDM_ol_mesh_t mesh\n");
    abort();
  }

  PDM_timer_hang_on(ol->timer);
}


/**
 * \brief Return the entitie sizes of the overlay mesh
 *
 * This function returns the entities sizes of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id            PDM_ol identifier
 * \param [in]  mesh          Input mesh
 *                            (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  i_part         Partition to define
 * \param [out] nOlFace       Number of faces of the overlay mesh
 * \param [out] nOlLinkedFace Number of linked faces
 * \param [out] nOlVtx        Number of vertices of the overlay mesh
 * \param [out] sOlface_vtx    Size of olface_vtx for each partition
 * \param [out] sInitToOlFace Size of initToOlFace for each partition
 *
 */

void
PDM_ol_part_mesh_dim_get
(
 const PDM_ol_t      *ol,
 const PDM_ol_mesh_t  mesh,
 const int            i_part,
       int           *nOlFace,
       int           *nOlLinkedFace,
       int           *nOlVtx,
       int           *sOlFaceIniVtx,
       int           *sOlface_vtx,
       int           *sInitToOlFace
)
{

  PDM_timer_resume(ol->timer);

  _ol_mesh_t *ol_mesh;

  switch(mesh) {

  case PDM_OL_MESH_A:
    ol_mesh = ol->olMeshA;
    break;

  case PDM_OL_MESH_B:
    ol_mesh = ol->olMeshB;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_mesh_dim_get : unknown PDM_ol_mesh_t mesh\n");
    abort();
  }

  PDM_surf_part_t *sp = ol_mesh->part[i_part]->part;

  *nOlFace       = sp->n_face;
  *nOlLinkedFace = ol_mesh->part[i_part]->nLinkedFace;
  *nOlVtx        = sp->n_vtx;
  *sOlFaceIniVtx = ol_mesh->part[i_part]->sFaceIniVtx;
  *sOlface_vtx   = sp->sface_vtx;
  *sInitToOlFace = ol_mesh->part[i_part]->sInitToOlFace;

  PDM_timer_hang_on(ol->timer);

}


/**
 * \brief Return the entitie of the overlay mesh
 *
 * This function returns the entities of the overlay mesh
 * for each partition of input meshA or input meshB
 *
 * \param [in]  id                   PDM_overlay identifier
 * \param [in]  mesh                 Input mesh
 *                                   (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
 * \param [in]  i_part                Mesh partition identifier
 * \param [out] olface_vtx_idx    Array adress of \ref olface_vtx index
 *                                   (size : \ref nOlFace + 1)
 * \param [out] olface_vtx       Array adress of face vertex connectivity
 *                                   (size : \ref sOlface_vtx[\ref i_part])
 * \param [out] olLinkedface_procIdx olLinkedFace Index (size = n_proc + 1)
 * \param [out] olLinkedFace  Array adress of linked face in other mesh
 *                              For each face, 4 link properties :
 *                                    - local face number
 *                                    - connected process,
 *                                    - connected part number,
 *                                    - connected local face number
 *                              (size : \ref 4 * nOlFace)
 * \param [out] olface_ln_to_gn    Array adress of local to global face numbering
 *                                   (size : \ref nOlFace)
 * \param [out] olCoords        Array adress of vertex coodinates
 *                                   (size : 3 * \ref nOlVtx)
 * \param [out] olvtx_ln_to_gn     Array adress of local to global vertex numbering array
 *                                   (size : \ref nOlVtx)
 * \param [out] initToOlFaceIdx Array adress of \ref initToOlFace index
 *                                   (size : \ref nOlVtx + 1)
 * \param [out] initToOlFace    Array adress of initial to ol faces
 * \param [out] initToOlVtx     Array adress of initial to ol vertices
 *
 */

void
PDM_ol_mesh_entities_get
(
const PDM_ol_t        *ol,
const PDM_ol_mesh_t    mesh,
const int              i_part,
      int            **olFaceIniVtxIdx,
      int            **olFaceIniVtx,
      int            **olface_vtx_idx,
      int            **olface_vtx,
      int            **olLinkedface_procIdx,
      int            **olLinkedFace,
      PDM_g_num_t    **olface_ln_to_gn,
      double         **olCoords,
      PDM_g_num_t    **olvtx_ln_to_gn,
      int            **initToOlFaceIdx,
      int            **initToOlFace
)
{

  PDM_timer_resume(ol->timer);

  _ol_mesh_t *ol_mesh;

  switch(mesh) {

  case PDM_OL_MESH_A:
    ol_mesh = ol->olMeshA;
    break;

  case PDM_OL_MESH_B:
    ol_mesh = ol->olMeshB;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_ol_mesh_dim_get : unknown PDM_ol_mesh_t mesh\n");
    abort();
  }

  PDM_surf_part_t *sp = ol_mesh->part[i_part]->part;

  *olFaceIniVtxIdx      = (int *) ol_mesh->part[i_part]->faceIniVtxIdx;
  *olFaceIniVtx         = (int *) ol_mesh->part[i_part]->faceIniVtx;
  *olface_vtx_idx       = (int *) sp->face_vtx_idx;
  *olface_vtx           = (int *) sp->face_vtx;
  *olLinkedface_procIdx = (int *) ol_mesh->part[i_part]->linkedFacesProcIdx;
  *olLinkedFace         = (int *) ol_mesh->part[i_part]->linkedFaces;
  *olface_ln_to_gn      = (PDM_g_num_t *) sp->face_ln_to_gn;
  *olCoords             = (double *) sp->coords;
  *olvtx_ln_to_gn       = (PDM_g_num_t *) sp->vtx_ln_to_gn;
  *initToOlFaceIdx      = (int *) ol_mesh->part[i_part]->initToOlFaceIdx;
  *initToOlFace         = (int *) ol_mesh->part[i_part]->initToOlFace;

  PDM_timer_hang_on(ol->timer);
}



/**
 * \brief Delete an overlay object
 *
 * This function deletes an overlay object
 *
 * \param [in]  id                PDM_ol identifier.
 *
 */

void
PDM_ol_del
(
 PDM_ol_t *ol
)
{
  /*
   * Delete
   */

  if (ol != NULL) {

    ol->meshA = PDM_surf_mesh_free(ol->meshA);
    ol->meshB = PDM_surf_mesh_free(ol->meshB);

    ol->olMeshA = _ol_mesh_free(ol->olMeshA);
    ol->olMeshB = _ol_mesh_free(ol->olMeshB);

    //PDM_dbbtree_free (ol->dbbtreeA);

    PDM_timer_free (ol->timer);

    /*
     * Update storaged array
     */

    free(ol);
  }
}


/**
 * \brief Dump elapsed an CPU time
 *
 *
 * \param [in]  id                PDM_ol identifier.
 *
 */

void
PDM_ol_dump_times
(
 PDM_ol_t *ol
)
{
  _ol_dump_times (ol);
}


#undef _DOT_PRODUCT
#undef _MODULE
#undef _MIN
#undef _MAX
#undef _CROSS_PRODUCT_3D

#ifdef __cplusplus
}
#endif /* __cplusplus */
