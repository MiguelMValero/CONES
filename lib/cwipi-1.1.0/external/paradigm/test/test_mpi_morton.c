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
#include "pdm_morton.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"








int
main
(
 int   argc,
 char *argv[]
 )
{

  //printf("%zu / %zu\n",sizeof(PDM_morton_code_t), 4*sizeof(PDM_morton_int_t));
  //printf("%zu / %zu\n",sizeof(_node_t), sizeof(PDM_morton_code_t) + sizeof(int) + 6*sizeof(double));
  //printf("%zu / %zu\n",sizeof(_explicit_node_t), (6+1+1+4+2)*8);


  PDM_MPI_Init (&argc, &argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int n_rank, i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  assert(n_rank == 3);

  int n;
  if (i_rank == 0) {
    n = 2;
  } else if (i_rank == 1) {
    n = 3;
  } else {
    n = 4;
  }


  int *all_n = malloc (sizeof(int) * n_rank);
  PDM_MPI_Allgather(&n, 1, PDM_MPI_INT, all_n, 1, PDM_MPI_INT, comm);

  PDM_log_trace_array_int(all_n, n_rank, "all_n : ");


  int *index = PDM_array_new_idx_from_sizes_int (all_n, n_rank);


  int *all_data = malloc (sizeof(int) * index[n_rank]);

  int *send_data = all_data + index[i_rank];
  for (int i = 0; i < n; i++) {
    send_data[i] = 100*i_rank + i;
  }

  PDM_log_trace_array_int (send_data, n, "send_data : ");

  /*int *recv_data = malloc (sizeof(int) * index[n_rank]);
  PDM_MPI_Allgatherv (send_data, n,            PDM_MPI_INT,
                      recv_data, all_n, index, PDM_MPI_INT,
                      comm);*/
  int *recv_data = all_data;
  PDM_MPI_Allgatherv (PDM_MPI_IN_PLACE, n,            PDM_MPI_INT,
                      recv_data,        all_n, index, PDM_MPI_INT,
                      comm);

  PDM_log_trace_array_int (recv_data, index[n_rank], "recv_data : ");

  free (all_n);
  free (index);
  free (all_data);

  /* Finalize */
  log_trace("-- End\n");
  PDM_MPI_Finalize ();

  return 0;
}













#if 0

typedef struct {
  int               n_points;
  double            pts_extents[6];
  PDM_morton_code_t code;

} _node_t;

typedef struct {
  int                n_nodes;
  PDM_morton_code_t *codes;
  int               *n_points;
} _l_octant_t;

typedef struct {
  PDM_morton_code_t code;
  int               n_points;
  int               range;
  int               ancestor_id;
  int               children_id[8];
  int               leaf_id; //-1 if internal, >=0 if leaf
  double            pts_extents[6];
} _explicit_node_t;

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

  //printf("%zu / %zu\n",sizeof(PDM_morton_code_t), 4*sizeof(PDM_morton_int_t));
  //printf("%zu / %zu\n",sizeof(_node_t), sizeof(PDM_morton_code_t) + sizeof(int) + 6*sizeof(double));
  //printf("%zu / %zu\n",sizeof(_explicit_node_t), (6+1+1+4+2)*8);


  PDM_MPI_Init (&argc, &argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int n_rank, i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  assert (n_rank == 2);


  PDM_MPI_Comm comm_node;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_SHARED, &comm_node);

  int n_rank_in_node, i_rank_in_node;
  PDM_MPI_Comm_rank (comm_node, &i_rank_in_node);
  PDM_MPI_Comm_size (comm_node, &n_rank_in_node);

  assert (n_rank_in_node == 2);

  int n_codes;
  PDM_morton_code_t *codes = NULL;

  log_trace("--- Morton codes ---\n");

  if (i_rank == 0) {
    n_codes = 3;
    PDM_MPI_Send (&n_codes, 1, PDM_MPI_INT, 1, 0, comm);

    codes = malloc (sizeof(PDM_morton_code_t) * n_codes);

    codes[0].L = 3;
    codes[0].X[0] = 0;
    codes[0].X[1] = 1;
    codes[0].X[2] = 2;

    codes[1].L = 6;
    codes[1].X[0] = 3;
    codes[1].X[1] = 4;
    codes[1].X[2] = 5;

    codes[2].L = 9;
    codes[2].X[0] = 6;
    codes[2].X[1] = 7;
    codes[2].X[2] = 8;

    PDM_MPI_Send (codes, n_codes*sizeof(PDM_morton_code_t), PDM_MPI_BYTE, 1, 1, comm);
  }

  else {
    PDM_MPI_Recv (&n_codes, 1, PDM_MPI_INT, 0, 0, comm);

    codes = malloc (sizeof(PDM_morton_code_t) * n_codes);

    PDM_MPI_Recv (codes, n_codes*sizeof(PDM_morton_code_t), PDM_MPI_BYTE, 0, 1, comm);
  }

  for (int i = 0; i < n_codes; i++) {
    log_trace("[%d] : L = %u, X = (%u, %u, %u)\n",
              i, codes[i].L, codes[i].X[0], codes[i].X[1], codes[i].X[2]);
  }

  free (codes);



  log_trace("--- Struct ---\n");



  int n_nodes;
  _node_t *nodes = NULL;

  if (i_rank == 0) {

    n_nodes = 3;
    PDM_MPI_Send (&n_nodes, 1, PDM_MPI_INT, 1, 0, comm);

    nodes = malloc (sizeof(_node_t) * n_nodes);

    nodes[0].code.L = 3;
    nodes[0].code.X[0] = 0;
    nodes[0].code.X[1] = 1;
    nodes[0].code.X[2] = 2;
    nodes[0].n_points  = 1;
    nodes[0].pts_extents[0] = 0.;
    nodes[0].pts_extents[1] = 1.;
    nodes[0].pts_extents[2] = 2.;
    nodes[0].pts_extents[3] = 4.;
    nodes[0].pts_extents[4] = 5.;
    nodes[0].pts_extents[5] = 6.;

    nodes[1].code.L = 6;
    nodes[1].code.X[0] = 3;
    nodes[1].code.X[1] = 4;
    nodes[1].code.X[2] = 5;
    nodes[1].n_points  = 2;
    nodes[1].pts_extents[0] =  7.;
    nodes[1].pts_extents[1] =  8.;
    nodes[1].pts_extents[2] =  9.;
    nodes[1].pts_extents[3] = 10.;
    nodes[1].pts_extents[4] = 11.;
    nodes[1].pts_extents[5] = 12.;

    nodes[2].code.L = 9;
    nodes[2].code.X[0] = 6;
    nodes[2].code.X[1] = 7;
    nodes[2].code.X[2] = 8;
    nodes[2].n_points  = 3;
    nodes[2].pts_extents[0] = 13.;
    nodes[2].pts_extents[1] = 14.;
    nodes[2].pts_extents[2] = 15.;
    nodes[2].pts_extents[3] = 16.;
    nodes[2].pts_extents[4] = 17.;
    nodes[2].pts_extents[5] = 18.;

    PDM_MPI_Send (nodes, n_nodes*sizeof(_node_t), PDM_MPI_BYTE, 1, 1, comm);
  }

  else {
    PDM_MPI_Recv (&n_nodes, 1, PDM_MPI_INT, 0, 0, comm);

    nodes = malloc (sizeof(_node_t) * n_nodes);

    PDM_MPI_Recv (nodes, n_nodes*sizeof(_node_t), PDM_MPI_BYTE, 0, 1, comm);
  }

  for (int i = 0; i < n_nodes; i++) {
    log_trace("[%d] : L = %u, X = (%u, %u, %u), n_pts = %d, ext = %.3f %.3f %.3f %.3f %.3f %.3f\n",
              i, nodes[i].code.L, nodes[i].code.X[0], nodes[i].code.X[1], nodes[i].code.X[2],
              nodes[i].n_points,
              nodes[i].pts_extents[0], nodes[i].pts_extents[1], nodes[i].pts_extents[2],
              nodes[i].pts_extents[3], nodes[i].pts_extents[4], nodes[i].pts_extents[5]);
  }

  free (nodes);








  log_trace("--- Shared Window (nodes) ---\n");

  int n_nodes_shared = 2;

  _node_t *nodes_to_copy = NULL;

  if (i_rank_in_node == 1) {
    nodes_to_copy = malloc (sizeof(_node_t) * n_nodes_shared);
    nodes_to_copy[0].code.L = 3;
    nodes_to_copy[0].code.X[0] = 0;
    nodes_to_copy[0].code.X[1] = 1;
    nodes_to_copy[0].code.X[2] = 2;
    nodes_to_copy[0].n_points  = 1;
    nodes_to_copy[0].pts_extents[0] = 0.;
    nodes_to_copy[0].pts_extents[1] = 1.;
    nodes_to_copy[0].pts_extents[2] = 2.;
    nodes_to_copy[0].pts_extents[3] = 4.;
    nodes_to_copy[0].pts_extents[4] = 5.;
    nodes_to_copy[0].pts_extents[5] = 6.;

    nodes_to_copy[1].code.L = 6;
    nodes_to_copy[1].code.X[0] = 3;
    nodes_to_copy[1].code.X[1] = 4;
    nodes_to_copy[1].code.X[2] = 5;
    nodes_to_copy[1].n_points  = 2;
    nodes_to_copy[1].pts_extents[0] =  7.;
    nodes_to_copy[1].pts_extents[1] =  8.;
    nodes_to_copy[1].pts_extents[2] =  9.;
    nodes_to_copy[1].pts_extents[3] = 10.;
    nodes_to_copy[1].pts_extents[4] = 11.;
    nodes_to_copy[1].pts_extents[5] = 12.;
  }


  PDM_mpi_win_shared_t *win = PDM_mpi_win_shared_create(n_nodes_shared, sizeof(_node_t), comm_node);



  _node_t *buff = PDM_mpi_win_shared_get(win);

  PDM_mpi_win_shared_lock_all(0, win);

  if (i_rank_in_node == 1) {
    memcpy(buff, nodes_to_copy, sizeof(_node_t) * n_nodes_shared);
  }

  PDM_mpi_win_shared_unlock_all(win);

  for (int i = 0; i < n_nodes_shared; i++) {
    log_trace("[%d] : L = %u, X = (%u, %u, %u), n_pts = %d, ext = %.3f %.3f %.3f %.3f %.3f %.3f\n",
              i, buff[i].code.L, buff[i].code.X[0], buff[i].code.X[1], buff[i].code.X[2],
              buff[i].n_points,
              buff[i].pts_extents[0], buff[i].pts_extents[1], buff[i].pts_extents[2],
              buff[i].pts_extents[3], buff[i].pts_extents[4], buff[i].pts_extents[5]);
  }

  PDM_mpi_win_shared_free(win);

  if (nodes_to_copy != NULL) free (nodes_to_copy);





  log_trace("--- Shared Window (octants) ---\n");


  _l_octant_t *octants_to_copy = NULL;

  //malloc (sizeof(l_octant_t));











  /* Finalize */

  PDM_MPI_Finalize ();

  return 0;
}
#endif
