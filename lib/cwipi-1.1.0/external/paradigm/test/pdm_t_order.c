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
#include "pdm_distant_neighbor.h"
#include "pdm_order.h"
#include "pdm_printf.h"
#include "pdm_error.h"


static inline
int
is_same_triplet
(
int iproc1, int ipart1, int ielt1,
int iproc2, int ipart2, int ielt2
)
{
  if(iproc1 == iproc2){
    if(ipart1 == ipart2){
      if(ielt1 == ielt2){
        return 1;
      }
    }
  }
  return 0;
}

static
void
compute_unique_idx
(
 int order[],
 int order_unique[],
 int connect_triplet[],
 const size_t nb_ent
)
{
  if(nb_ent == 0){
    return;
  }

  int idx        = 0;
  int idx_unique = 0;
  int last_idx   = order[idx++];
  int last_proc  = connect_triplet[3*last_idx  ];
  int last_part  = connect_triplet[3*last_idx+1];
  int last_elmt  = connect_triplet[3*last_idx+2];
  order_unique[idx_unique] = 0;

  for(int i = 1; i < (int) nb_ent; i++){

    int curr_idx  = order[idx++];
    int curr_proc = connect_triplet[3*curr_idx  ];
    int curr_part = connect_triplet[3*curr_idx+1];
    int curr_elmt = connect_triplet[3*curr_idx+2];
    int is_same = is_same_triplet(last_proc, last_part, last_elmt,
                                  curr_proc, curr_part, curr_elmt);
    // printf(" curr:: ( %d / %d / %d ) | last:: ( %d / %d / %d ) \n",
    //        curr_proc, curr_part, curr_elmt,
    //        last_proc, last_part, last_elmt);

    if(is_same == 0){ // N'est pas le meme
      idx_unique++;
      last_proc = curr_proc;
      last_part = curr_part;
      last_elmt = curr_elmt;
    }
    order_unique[i] = idx_unique;
    // printf("[%d] = %d --> %d \n", i, is_same, idx_unique);
  }


  if(0 == 1){
    printf("order_unique:: \n");
    for(int i = 0; i < (int) nb_ent; i++){
      printf(" -------------------------- \n");
      // int pos_unique = order_unique_j1[i];
      // int curr_idx   = order_j1[pos_unique];
      int pos_unique = order_unique[i];
      int curr_idx   = order[i];

      int curr_proc  = connect_triplet[3*curr_idx  ];
      int curr_part  = connect_triplet[3*curr_idx+1];
      int curr_elmt  = connect_triplet[3*curr_idx+2];

      printf("\t pos_unique :: %d \n", pos_unique);
      printf("\t curr_idx   :: %d \n", curr_idx  );
      printf("\t triplet    :: ( %d / %d / %d ) \n", curr_proc, curr_part, curr_elmt);

    }
    printf("\n");
  }
}


static
void
compute_unique_idx2
(
 int order[],
 int order_unique[],
 int connect_triplet[],
 const size_t nb_ent
)
{
  int idx_unique = -1;
  int last_proc  = -1;
  int last_part  = -1;
  int last_elmt  = -1;

  for(int i = 0; i < (int) nb_ent; i++){

    int old_order = order[i];
    int curr_proc = connect_triplet[3*old_order  ];
    int curr_part = connect_triplet[3*old_order+1];
    int curr_elmt = connect_triplet[3*old_order+2];
    int is_same = is_same_triplet(last_proc, last_part, last_elmt,
                                  curr_proc, curr_part, curr_elmt);
    printf(" curr:: ( %d / %d / %d ) | last:: ( %d / %d / %d ) \n",
           curr_proc, curr_part, curr_elmt,
           last_proc, last_part, last_elmt);

    if(is_same == 0){ // N'est pas le meme
      idx_unique++;
      last_proc = curr_proc;
      last_part = curr_part;
      last_elmt = curr_elmt;
    }
    order_unique[i] = idx_unique;
    // printf("[%d] = %d --> %d \n", i, is_same, idx_unique);
  }


  if(0 == 1){
    printf("order_unique:: \n");
    for(int i = 0; i < (int) nb_ent; i++){
      printf(" -------------------------- \n");
      // int pos_unique = order_unique_j1[i];
      // int curr_idx   = order_j1[pos_unique];
      int pos_unique = order_unique[i];
      int curr_idx   = order[i];

      int curr_proc  = connect_triplet[3*curr_idx  ];
      int curr_part  = connect_triplet[3*curr_idx+1];
      int curr_elmt  = connect_triplet[3*curr_idx+2];

      printf("\t pos_unique :: %d \n", pos_unique);
      printf("\t curr_idx   :: %d \n", curr_idx  );
      printf("\t triplet    :: ( %d / %d / %d ) \n", curr_proc, curr_part, curr_elmt);

    }
    printf("\n");
  }
}

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
  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  /*
   * Triplet : (iproc, i_part, index)
   *
   * ----------------------------|           |----------------------------
   *     |                       |           |               |
   *     |       (0,0,2)         |           |               |
   *     |                       |           |  (1, 0, 1)    |
   * ----------------------------|           |               |
   *     |                       |           |               |
   *     |      (0,0,1)          |           |----------------------------
   *     |                       |           |               |
   * ----------------------------|           |               |
   *     |                       |           |  (1, 0, 0)    |
   *     |      (0,0,0)          |           |               |
   *     |                       |           |               |
   * ----------------------------|           |----------------------------
   *                         join_1        join_2
   */

  const int stride  = 3;

  /* Connection de join1 avec join2 */
  // const int n_faces_j1 = 3;
  int connect_triplet_j1[12] = {// Fisrt
                                1, 0, 0,
                                // Second
                                1, 0, 0,
                                1, 0, 1,
                                // Third
                                1, 0, 1};

  /* Connection de join2 avec join1 */
  // const int n_faces_j2 = 3;
  int connect_triplet_j2[12] = {// First
                                0, 0, 0,
                                0, 0, 1,
                                // Second
                                0, 0, 1,
                                0, 0, 2};


  // Ordering
  int nb_ent = 4;
  int order_j1[nb_ent];
  PDM_order_lnum_s(connect_triplet_j1, stride, order_j1, nb_ent);

  // printf("order_j1:: ");
  // for(int i = 0; i < nb_ent; i++){
  //   printf("%d ", order_j1[i]);
  // }
  // printf("\n");

  int order_j2[nb_ent];
  PDM_order_lnum_s(connect_triplet_j2, stride, order_j2, nb_ent);

  // printf("order_j2:: ");
  // for(int i = 0; i < nb_ent; i++){
  //   printf("%d ", order_j2[i]);
  // }
  // printf("\n");

  /*
   * Unique
   */
  int order_unique_j1[nb_ent]; // Très mauvais nom
  int order_unique_j2[nb_ent]; // Très mauvais nom
  compute_unique_idx(order_j1, order_unique_j1, connect_triplet_j1, nb_ent);
  compute_unique_idx(order_j2, order_unique_j2, connect_triplet_j2, nb_ent);

  compute_unique_idx2(order_j1, order_unique_j1, connect_triplet_j1, nb_ent);
  compute_unique_idx2(order_j2, order_unique_j2, connect_triplet_j2, nb_ent);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;

}
