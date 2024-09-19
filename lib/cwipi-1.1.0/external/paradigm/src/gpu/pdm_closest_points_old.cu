// /*
//   This file is part of the ParaDiGM library.

//   Copyright (C) 2019       ONERA

//   This library is free software; you can redistribute it and/or
//   modify it under the terms of the GNU Lesser General Public
//   License as published by the Free Software Foundation; either
//   version 2.1 of the License, or (at your option) any later version.

//   This library is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.

//   You should have received a copy of the GNU Lesser General Public
//   License along with this library; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
// */
// /*----------------------------------------------------------------------------
//  * System headers
//  *----------------------------------------------------------------------------*/

// #include <assert.h>
// #include <stdio.h>
// #include <math.h>
// #include <string.h>
// #include <stdlib.h>

// /*----------------------------------------------------------------------------
//  * Local headers
//  *----------------------------------------------------------------------------*/

// #include "pdm_cuda_error.cuh"
// #include "pdm_cuda.cuh"
// #include "pdm_printf.h"
// #include "pdm_error.h"
// #include "pdm.h"
// #include "pdm_priv.h"
// #include "pdm_handles.h"
// #include "pdm_handles.cuh"
// #include "pdm_mpi.h"
// #include "pdm_timer.h"
// #include "pdm_timer.cuh"
// #include "pdm_closest_points.h"
// #include "pdm_closest_points.cuh"
// #include "pdm_para_octree.h"
// #include "pdm_para_octree.cuh"

// /*----------------------------------------------------------------------------*/

// #ifdef	__cplusplus
// extern "C" {
// #if 0
// } /* Fake brace to force Emacs auto-indentation back to column 0 */
// #endif
// #endif


// /*============================================================================
//  * Macro definitions
//  *============================================================================*/


// #define NTIMER 2

// /*============================================================================
//  * Type definitions
//  *============================================================================*/


// /**
//  * \enum _timer_step_t
//  *
//  */

// typedef enum {

//   BEGIN    = 0,
//   END      = 1,

// } _timer_step_t;


// /**
//  * \struct _tgt_point_cloud_t
//  * \brief  Target point cloud structure
//  *
//  */

// typedef struct {

//   int           n_part;            /*!< Number of partition */
//   int          *n_points;          /*!< Number of points of each partition */
//   double      **coords;            /*!< Point coordinates points of each partition */
//   PDM_g_num_t **gnum;              /*!< Point global numbering of each partition */
//   PDM_g_num_t **closest_src_gnum;  /*!< Global numbering of the n_closest source points
//                                      for each point of each partition  */
//   double      **closest_src_dist; /*!< Distance to the n_closest source points
//                                     for each point of each partition  */

// } _tgt_point_cloud_t;


// /**
//  * \struct _src_point_cloud_t
//  * \brief  Src point cloud structure
//  *
//  */

// typedef struct {

//   int           n_part;            /*!< Number of partition */
//   int          *n_points;          /*!< Number of points of each partition */
//   double      **coords;            /*!< Point coordinates points of each partition */
//   PDM_g_num_t **gnum;              /*!< Point global numbering of each partition */

// } _src_point_cloud_t;


// /**
//  * \struct _PDM_closest_t
//  * \brief  Closest points structure
//  *
//  */

// typedef struct {

//   PDM_MPI_Comm comm;  /*!< MPI communicator */

//   int n_closest;  /*!< Number of closest source points to find for each
//                     target point  */

//   _src_point_cloud_t *src_cloud; /*!< Source point cloud */

//   _tgt_point_cloud_t *tgt_cloud; /*!< Target point cloud */

//   PDM_timer_t *timer; /*!< Timer */

//   double times_elapsed[NTIMER]; /*!< Elapsed time */

//   double times_cpu[NTIMER];     /*!< CPU time */

//   double times_cpu_u[NTIMER];  /*!< User CPU time */

//   double times_cpu_s[NTIMER];  /*!< System CPU time */


// } _PDM_closest_t;


// /*============================================================================
//  * Global variable
//  *============================================================================*/

// static PDM_Handles_t *_closest_pts   = NULL;

// static int idebug = 0;

// /*=============================================================================
//  * Private function definitions
//  *============================================================================*/

// /**
//  *
//  * \brief Return ppart object from it identifier
//  *
//  * \param [in]   ppart_id        ppart identifier
//  *
//  */

// static _PDM_closest_t *
// _get_from_id
// (
//  int  id
//  )
// {
//   printf("before handle get\n");
//   printf("id = %d\n", id);
//   printf("s array = %d\n", _closest_pts->s_array);
//   _PDM_closest_t *closest = (_PDM_closest_t *) PDM_Handles_get (_closest_pts, id);

//   if (closest == NULL) {
//     PDM_error(__FILE__, __LINE__, 0, "PDM_closest_points error : Bad identifier\n");
//   }

//   return closest;
// }


// /*============================================================================
//  * Public function definitions
//  *============================================================================*/



// /**
//  *
//  * \brief Set the number of partitions of a point cloud
//  *
//  * \param [in]   id                Identifier
//  * \param [in]   n_part_cloud_src  Number of partitions of the source cloud
//  * \param [in]   n_part_cloud_tgt  Number of partitions of the target cloud
//  *
//  */

// /* void
// PDM_closest_points_n_part_cloud_set
// (
//  const int  id,
//  const int  n_part_cloud_src,
//  const int  n_part_cloud_tgt
//  )
// {
//   _PDM_closest_t *cls = _get_from_id (id);
//   assert(cls->src_cloud == NULL);
//   assert(cls->tgt_cloud == NULL);

//   cls->src_cloud = malloc (sizeof(_src_point_cloud_t));
//   cls->tgt_cloud = malloc (sizeof(_tgt_point_cloud_t));

//   cls->src_cloud->n_part = n_part_cloud_src;
//   cls->src_cloud->coords = malloc (sizeof(double *) * n_part_cloud_src);
//   cls->src_cloud->gnum = malloc (sizeof(int *) * n_part_cloud_src);
//   cls->src_cloud->n_points = malloc (sizeof(int) * n_part_cloud_src);

//   cls->tgt_cloud->n_part = n_part_cloud_tgt;
//   cls->tgt_cloud->coords = malloc (sizeof(double *) * n_part_cloud_tgt);
//   cls->tgt_cloud->gnum = malloc (sizeof(int *) * n_part_cloud_tgt);
//   cls->tgt_cloud->n_points = malloc (sizeof(int) * n_part_cloud_tgt);
//   cls->tgt_cloud->closest_src_gnum = NULL;
//   cls->tgt_cloud->closest_src_dist = NULL;
// } */


// /**
//  *
//  * \brief Look for closest points
//  *
//  * \param [in]   id  Identifier
//  *
//  */

// void
// PDM_closest_points_compute_GPU
// (
//  const int id,
//  PDM_Handles_t *var
//  )
// {
//   _closest_pts = var;
//   printf("s array = %d\n", _closest_pts->s_array);
//   _PDM_closest_t *cls = NULL;
//   int *octree_id = NULL;

//   //Allocate data on unified memory so it is accessible from CPU or GPU
//   gpuErrchk(cudaMallocManaged(&cls, sizeof(_PDM_closest_t)));
//   gpuErrchk(cudaMallocManaged(&octree_id, sizeof(int)));
//   printf("After malloc\n");

//   cls = _get_from_id (id);

//   printf("after cls\n");

//   double b_t_elapsed;
//   double b_t_cpu;
//   double b_t_cpu_u;
//   double b_t_cpu_s;

//   cls->times_elapsed[BEGIN] = PDM_timer_elapsed(cls->timer);
//   cls->times_cpu[BEGIN]     = PDM_timer_cpu(cls->timer);
//   cls->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(cls->timer);
//   cls->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(cls->timer);

//   b_t_elapsed = cls->times_elapsed[BEGIN];
//   b_t_cpu     = cls->times_cpu[BEGIN];
//   b_t_cpu_u   = cls->times_cpu_u[BEGIN];
//   b_t_cpu_s   = cls->times_cpu_s[BEGIN];
//   PDM_timer_resume(cls->timer);


//   int i_rank;
//   PDM_MPI_Comm_rank (cls->comm, &i_rank);

//   //-->GPU
//   const int depth_max = 31;//?
//   const int points_in_leaf_max = 1;//2*cls->n_closest;//?
//   const int build_leaf_neighbours = 1;

//   printf("Before create\n");
//   /* Create empty parallel octree structure */
//   *octree_id = PDM_para_octree_create_GPU (cls->src_cloud->n_part,
//                                           depth_max,
//                                           points_in_leaf_max,
//                                           build_leaf_neighbours,
//                                           cls->comm);

//   printf("gpu print:\n");
//   print_from_gpu<<<1,1>>>(*octree_id);
//   cudaDeviceSynchronize();
//   printf("after gpu print\n");

//   /* Set source point clouds */
//   for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
//     PDM_para_octree_point_cloud_set (*octree_id,
//                                      i_part,
//                                      cls->src_cloud->n_points[i_part],
//                                      cls->src_cloud->coords[i_part],
//                                      cls->src_cloud->gnum[i_part]);
//   }


//   /* Build parallel octree */
//   PDM_para_octree_build (*octree_id);
//   //PDM_para_octree_dump (octree_id);
//   PDM_para_octree_dump_times (*octree_id);
//   //<--


//   // /* Concatenate partitions */
//   // int n_tgt = 0;
//   // for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++)
//   //   n_tgt += cls->tgt_cloud->n_points[i_part];

//   // double      *tgt_coord = malloc (sizeof(double)      * n_tgt * 3);
//   // PDM_g_num_t *tgt_g_num = malloc (sizeof(PDM_g_num_t) * n_tgt);
//   // PDM_g_num_t *closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_tgt * cls->n_closest);
//   // double      *closest_src_dist = malloc (sizeof(double)      * n_tgt * cls->n_closest);

//   // n_tgt = 0;
//   // for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
//   //   for (int i = 0; i < cls->tgt_cloud->n_points[i_part]; i++) {
//   //     for (int j = 0; j < 3; j++)
//   //       tgt_coord[n_tgt + 3*i + j] = cls->tgt_cloud->coords[i_part][3*i + j];
//   //     tgt_g_num[n_tgt + i] = cls->tgt_cloud->gnum[i_part][i];
//   //   }
//   //   n_tgt += cls->tgt_cloud->n_points[i_part];
//   // }
  

//   // /* Search closest source points from target points */
//   // PDM_para_octree_closest_point (octree_id,
//   //                                cls->n_closest,
//   //                                n_tgt,
//   //                                tgt_coord,
//   //                                tgt_g_num,
//   //                                closest_src_gnum,
//   //                                closest_src_dist);


//   // /* Restore partitions */
//   // free (tgt_coord);
//   // free (tgt_g_num);
//   // n_tgt = 0;

//   // cls->tgt_cloud->closest_src_gnum = malloc (sizeof(PDM_g_num_t *) * cls->tgt_cloud->n_part);
//   // cls->tgt_cloud->closest_src_dist = malloc (sizeof(double *)      * cls->tgt_cloud->n_part);

//   // for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
//   //   int s_closest_src = cls->n_closest * cls->tgt_cloud->n_points[i_part];

//   //   cls->tgt_cloud->closest_src_gnum[i_part] = malloc (sizeof(PDM_g_num_t) * s_closest_src);
//   //   cls->tgt_cloud->closest_src_dist[i_part] = malloc (sizeof(double)      * s_closest_src);

//   //   for (int i = 0; i < cls->tgt_cloud->n_points[i_part]; i++) {
//   //     for (int j = 0; j < cls->n_closest; j++) {
//   //       cls->tgt_cloud->closest_src_gnum[i_part][cls->n_closest*i+j] =
//   //         closest_src_gnum[n_tgt + cls->n_closest*i + j];

//   //       cls->tgt_cloud->closest_src_dist[i_part][cls->n_closest*i+j] =
//   //         closest_src_dist[n_tgt + cls->n_closest*i + j];
//   //     }
//   //   }
//   //   n_tgt += cls->n_closest * cls->tgt_cloud->n_points[i_part];
//   // }
//   // free (closest_src_gnum);
//   // free (closest_src_dist);



//   //-->GPU
//   /* Free parallel octree */
//   PDM_para_octree_free_GPU (*octree_id);
//   gpuErrchk(cudaFree(octree_id));
//   //<--


//   PDM_timer_hang_on(cls->timer);

//   cls->times_elapsed[END] = PDM_timer_elapsed(cls->timer);
//   cls->times_cpu[END]     = PDM_timer_cpu(cls->timer);
//   cls->times_cpu_u[END]   = PDM_timer_cpu_user(cls->timer);
//   cls->times_cpu_s[END]   = PDM_timer_cpu_sys(cls->timer);

//   b_t_elapsed = cls->times_elapsed[END];
//   b_t_cpu     = cls->times_cpu[END];
//   b_t_cpu_u   = cls->times_cpu_u[END];
//   b_t_cpu_s   = cls->times_cpu_s[END];
//   PDM_timer_resume(cls->timer);
// }


// #ifdef	__cplusplus
// }
// #endif
// #undef NTIMER
