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
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"

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

  int n_rank, i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_MPI_Comm comm_numa;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_numa);

  PDM_MPI_Comm comm_node;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_SHARED, &comm_node);

  int n_rank_in_numa, i_rank_in_numa;
  PDM_MPI_Comm_rank (comm_numa, &i_rank_in_numa);
  PDM_MPI_Comm_size (comm_numa, &n_rank_in_numa);

  int n_rank_in_node, i_rank_in_node;
  PDM_MPI_Comm_rank (comm_node, &i_rank_in_node);
  PDM_MPI_Comm_size (comm_node, &n_rank_in_node);

  log_trace(" i_rank_in_numa = %d/%d - i_rank_in_node = %d/%d  - i_rank = %d/%d \n", i_rank_in_numa, n_rank_in_numa,
                                                                                     i_rank_in_node, n_rank_in_node,
                                                                                     i_rank, n_rank);


  PDM_MPI_Comm comm_master_of_numa = PDM_MPI_get_group_of_master(comm, comm_numa);

  if(comm_master_of_numa != PDM_MPI_COMM_NULL) {
    int n_rank_node, i_rank_node;
    PDM_MPI_Comm_rank (comm_master_of_numa, &i_rank_node);
    PDM_MPI_Comm_size (comm_master_of_numa, &n_rank_node);

    log_trace(" i_rank_node = %d/%d - i_rank = %d/%d \n", i_rank_node, n_rank_node, i_rank, n_rank);
  }

  PDM_g_num_t n_g_data = 100; // Attention taille global !!!
  PDM_mpi_win_shared_t* win = PDM_mpi_win_shared_create(n_g_data, sizeof(int), comm_numa);
  int* buff = PDM_mpi_win_shared_get(win);

  PDM_g_num_t* distrib_numa = PDM_compute_uniform_entity_distribution(comm_numa, n_g_data);

  PDM_mpi_win_shared_lock_all(0, win);
  for(int i = distrib_numa[i_rank_in_numa]; i < distrib_numa[i_rank_in_numa+1]; ++i) {
    // buff[i] = i;
    buff[i] = i_rank;
  }
  PDM_mpi_win_shared_unlock_all(win);

  PDM_log_trace_array_int(buff, n_g_data, "buff::");

  free(distrib_numa);

  PDM_MPI_Comm_free(&comm_shared);
  PDM_MPI_Comm_free(&comm_node);
  PDM_mpi_win_shared_free(win);
  PDM_MPI_Finalize ();

  return 0;
}
