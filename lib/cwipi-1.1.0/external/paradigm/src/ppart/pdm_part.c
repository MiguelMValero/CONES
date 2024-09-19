/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_multipart.h"
#include "pdm_part_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"
#include "pdm_para_graph_dual.h"
#include "pdm_array.h"

#include "pdm_part_geom.h"
#include "pdm_part_renum.h"
#include "pdm_fortran_to_c_string.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_quick_sort.h"
#include "pdm_binary_search.h"
#include "pdm_logging.h"


/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Call a couple MPI_Alltoall MPI_Alltoallv
 *
 * \param [in]   send_buff,            Sending buffer
 * \param [in]   send_buff_n,           Number of data to send to each process
 *                                    (size : communicator size)
 * \param [in]   send_buff_idx,         Index in send_buff for each process
 *                                    (size : communicator size)
 * \param [out]  recv_buff,            Receiving buffer
 * \param [out]  recv_buff_n           Receiving buffer size
 * \param [in]   exch_mpi_data_type   Data type to exchange
 * \param [in]   type_exch_size       Size of data type
 * \param [in]   comm                 Communicator
 *
 */

static void
_alltoall
(
 void              *send_buff,
 int               *send_buff_n,
 int               *send_buff_idx,
 void             **recv_buff,
 int               *recv_buff_n,
 int               *recv_buff_idx,
 PDM_MPI_Datatype   MPIDataType,
 size_t             MPIDataTypeSize,
 PDM_MPI_Comm       comm
)
{
  int n_rank = 0;
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Get number data to receive from each process */

  PDM_MPI_Alltoall(send_buff_n,
                   1,
                   PDM_MPI_INT,
                   recv_buff_n,
                   1,
                   PDM_MPI_INT,
                   comm);

  PDM_array_idx_from_sizes_int(recv_buff_n, n_rank, recv_buff_idx);

  *recv_buff = malloc(recv_buff_idx[n_rank] * MPIDataTypeSize);

  /* Receive data from each process */

  PDM_MPI_Alltoallv(send_buff,
                    send_buff_n,
                    send_buff_idx,
                    MPIDataType,
                    *recv_buff,
                    recv_buff_n,
                    recv_buff_idx,
                    MPIDataType,
                    comm);

}

/**
 *
 * \brief Builds dual graph face cell connectivity
 *
 * \param [inout] ppart       Ppart object
 *
 */

static void
_dual_graph_from_face_cell
(
 _PDM_part_t *ppart
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(ppart->comm, &i_rank);
  PDM_MPI_Comm_size(ppart->comm, &n_rank);

  /*
   * cell_to_send_n allocation
   */

  int *cell_to_send_n = PDM_array_zeros_int(n_rank);

  const int n_data = 3; /* Number data to send */

  /*
   * Set cell list to send to each process
   */

  for (int i = 0; i < ppart->dn_face; i++) {
    PDM_g_num_t i_cell1 = PDM_ABS (ppart->_dface_cell[2*i    ]);
    PDM_g_num_t i_cell2 = PDM_ABS (ppart->_dface_cell[2*i + 1]);

    int i_rank1 = PDM_search_rank(i_cell1, ppart->dcell_proc, 0, n_rank);
    cell_to_send_n[i_rank1] += n_data;

    if (i_cell2 > 0) {
      int i_rank2 = PDM_search_rank(i_cell2, ppart->dcell_proc, 0, n_rank);
      cell_to_send_n[i_rank2] += n_data;
    }
  }

  /*
   * Create index aray
   */

  int *cell_to_send_idx = PDM_array_new_idx_from_sizes_int(cell_to_send_n, n_rank);
  PDM_array_reset_int(cell_to_send_n, n_rank, 0);

  PDM_g_num_t *cell_to_send = (PDM_g_num_t *) malloc(cell_to_send_idx[n_rank] * sizeof(PDM_g_num_t));

  /*
   * Stores pair of cells to send to the others processes
   */

  for (int i = 0; i < ppart->dn_face ; i++) {
    PDM_g_num_t i_cell1 = PDM_ABS (ppart->_dface_cell[2*i    ]);
    PDM_g_num_t i_cell2 = PDM_ABS (ppart->_dface_cell[2*i + 1]);
    int i_rank1 = PDM_search_rank(i_cell1, ppart->dcell_proc, 0, n_rank);

    int idx1             = cell_to_send_idx[i_rank1] + cell_to_send_n[i_rank1];
    cell_to_send[idx1  ]   = i_cell1;
    cell_to_send[idx1+1]   = i_cell2;
    cell_to_send[idx1+2]   = ppart->dface_proc[i_rank] + i;
    cell_to_send_n[i_rank1] += n_data;

    if (i_cell2 > 0) {
      int i_rank2 = PDM_search_rank(i_cell2, ppart->dcell_proc, 0, n_rank);
      int idx2             = cell_to_send_idx[i_rank2] + cell_to_send_n[i_rank2];
      cell_to_send[idx2  ]   = i_cell2;
      cell_to_send[idx2+1]   = i_cell1;
      cell_to_send[idx2+2]   = ppart->dface_proc[i_rank] + i;
      cell_to_send_n[i_rank2] += n_data;
    }
  }

  /*
   * Receive pair of Cells from the others processes
   */

  int *cell_to_recv_n = (int *) malloc(n_rank * sizeof(int));

  PDM_MPI_Alltoall(cell_to_send_n,
                   1,
                   PDM_MPI_INT,
                   cell_to_recv_n,
                   1,
                   PDM_MPI_INT,
                   ppart->comm);

  int *cell_to_recv_idx =  PDM_array_new_idx_from_sizes_int(cell_to_recv_n, n_rank);

  PDM_g_num_t *cell_to_recv = (PDM_g_num_t *) malloc(cell_to_recv_idx[n_rank]*sizeof(PDM_g_num_t));

  PDM_MPI_Alltoallv(cell_to_send,
                    cell_to_send_n,
                    cell_to_send_idx,
                    PDM__PDM_MPI_G_NUM,
                    cell_to_recv,
                    cell_to_recv_n,
                    cell_to_recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    ppart->comm);

  int n_recv_pair = cell_to_recv_idx[n_rank]/n_data;

  /*
   * Free
   */

  free(cell_to_send_idx);
  free(cell_to_send_n);
  free(cell_to_send);
  free(cell_to_recv_idx);
  free(cell_to_recv_n);

  cell_to_send_idx  = NULL;
  cell_to_send_n    = NULL;
  cell_to_send      = NULL;
  cell_to_recv_idx  = NULL;
  cell_to_recv_n    = NULL;

  /*
   * Count neighbour cells for each cell
   */

  int *n_neighbour = (int *) malloc(ppart->dn_cell * sizeof(int));

  int have_dcell_face = 0;
  if (ppart->_dcell_face_idx != NULL)
    have_dcell_face = 1;

  int *dcell_face_n = NULL;
  if (!have_dcell_face) {
    dcell_face_n = PDM_array_zeros_int(ppart->dn_cell);
  }

  for (int i = 0; i < ppart->dn_cell; i++) {
    n_neighbour[i]= 0;
  }

  for (int i = 0; i < n_recv_pair; i++) {
    PDM_g_num_t  gelt1 = cell_to_recv[n_data*i  ];                         // Get global numbering
    PDM_g_num_t  gelt2 = cell_to_recv[n_data*i+1];                         // Get global numbering
    PDM_g_num_t  _lelt1 = gelt1 - ppart->dcell_proc[i_rank];
    int          lelt1 = (int) _lelt1;      // Switch to local numbering

    if (gelt2 > 0) {
      n_neighbour[lelt1] += 1;
    }
    if (!have_dcell_face) {
      dcell_face_n[lelt1] += 1;
    }
  }

  /*
   * Allocate dual graph from neighbour cells for each cell
   */

  if (!have_dcell_face) {
    ppart->dcell_face_idx = PDM_array_new_idx_from_sizes_int(dcell_face_n, ppart->dn_cell);
    ppart->dcell_face = PDM_array_const_gnum(ppart->dcell_face_idx[ppart->dn_cell], -1);
    PDM_array_reset_int(dcell_face_n, ppart->dn_cell, 0);

    ppart->_dcell_face_idx = ppart->dcell_face_idx;
    ppart->_dcell_face = ppart->dcell_face;
  }

  ppart->ddual_graph_idx = PDM_array_new_idx_from_sizes_gnum(n_neighbour, ppart->dn_cell);
  ppart->ddual_graph = PDM_array_const_gnum(ppart->ddual_graph_idx[ppart->dn_cell], -1);
  PDM_array_reset_int(n_neighbour, ppart->dn_cell, 0);

  /*
   * Complete dual graph
   */

  for (int i = 0; i < n_recv_pair; i++) {
    PDM_g_num_t  gcel1  = cell_to_recv[n_data*i];                      // global numbering
    PDM_g_num_t  _lcel1 = gcel1 - ppart->dcell_proc[i_rank];
    int           lcel1  = (int) _lcel1; // local numbering
    PDM_g_num_t  gcel2  = cell_to_recv[n_data*i+1];                    // global numbering
    PDM_g_num_t  gface2 = cell_to_recv[n_data*i+2];                    // global numbering

    if (!have_dcell_face) {
      ppart->dcell_face[ppart->dcell_face_idx[lcel1] + dcell_face_n[lcel1]] = gface2;
      dcell_face_n[lcel1] += 1;
    }

    /*
     * Search if cel2 is already stored (To optimize for polyhedra (lot of neighbours) ?)
     */

    if (gcel2 > 0) {

      PDM_g_num_t k;
      for (k = ppart->ddual_graph_idx[lcel1]; k < ppart->ddual_graph_idx[lcel1] + n_neighbour[lcel1]; k++) {
        if (ppart->ddual_graph[k] == gcel2 - 1)
          break;
      }

      if (k == ppart->ddual_graph_idx[lcel1] + n_neighbour[lcel1]) {
        ppart->ddual_graph[ppart->ddual_graph_idx[lcel1] + n_neighbour[lcel1]] = gcel2 - 1;
        n_neighbour[lcel1] += 1;
      }
    }
  }

  /*
   * Compress dual graph
   */

  int k = 0;
  int k1 = 0;
  while (k < ppart->ddual_graph_idx[ppart->dn_cell]) {
    if (ppart->ddual_graph[k] >= 0) {
      ppart->ddual_graph[k1] = ppart->ddual_graph[k];
      k++;
      k1++;
    }
    else
      k++;
  }

  /*
   * Reallocate to free unused memory
   */

  ppart->ddual_graph = realloc(ppart->ddual_graph, k1 * sizeof(PDM_g_num_t));
  PDM_array_idx_from_sizes_gnum(n_neighbour, ppart->dn_cell, ppart->ddual_graph_idx);

  /*
   * ppart->dcell_face_idx is ppart->ddual_graph_idx
   */

  if (1 == 0) {
    if (!have_dcell_face) {
      PDM_printf("ppart->_dcell_face : \n");
      for (int i = 0; i < ppart->dn_cell; i++) {
        for (int j = ppart->_dcell_face_idx[i]; j < ppart->_dcell_face_idx[i+1]; j++)
          PDM_printf(" "PDM_FMT_G_NUM, ppart->_dcell_face[j]);
        PDM_printf("\n");
      }
    }
  }

  free(cell_to_recv);
  free(n_neighbour);
  if (!have_dcell_face) {
    free(dcell_face_n);
  }
}

/**
 *
 * \brief Builds dual graph from face cell connectivity
 *
 * \param [inout] ppart       Ppart object
 *
 */

static void
_dual_graph_from_cell_face
(
 _PDM_part_t *ppart
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(ppart->comm, &i_rank);
  PDM_MPI_Comm_size(ppart->comm, &n_rank);

  /*
   * cell_to_send_n allocation
   */

  int *face_to_send_n = PDM_array_zeros_int(n_rank);

  const int n_data = 2; /* Number data to send */

  /*
   * Set cell list to send to each process
   */

  for (int i = 0; i < ppart->dn_cell; i++) {
    for (int j = ppart->_dcell_face_idx[i]; j < ppart->_dcell_face_idx[i+1]; j++) {
      PDM_g_num_t iface = PDM_ABS(ppart->_dcell_face[j]);

      int found_rank = PDM_search_rank(iface, ppart->dface_proc, 0, n_rank);
      face_to_send_n[found_rank] += n_data;
    }
  }

  /*
   * Create index aray
   */

  int *face_to_send_idx = PDM_array_new_idx_from_sizes_int(face_to_send_n, n_rank);
  PDM_array_reset_int(face_to_send_n, n_rank, 0);

  PDM_g_num_t *face_to_send =
    (PDM_g_num_t *) malloc(face_to_send_idx[n_rank] * sizeof(PDM_g_num_t));

  /*
   * Stores faces to send to the others processes
   */

  for (int i = 0; i < ppart->dn_cell; i++) {
    for (int j = ppart->_dcell_face_idx[i]; j < ppart->_dcell_face_idx[i+1]; j++) {
      PDM_g_num_t iface = PDM_ABS (ppart->_dcell_face[j]);

      int found_rank = PDM_search_rank(iface, ppart->dface_proc, 0, n_rank);
      int idx        = face_to_send_idx[found_rank] + face_to_send_n[found_rank];
      face_to_send[idx  ]   = iface;
      face_to_send[idx+1]   = ppart->dcell_proc[i_rank] + i;//ppart->dcell_proc[found_rank] + i;
      face_to_send_n[found_rank] += n_data;
    }
  }

  /*
   * Receive faces from the others processes
   */

  int *face_to_recv_n = (int *) malloc(n_rank * sizeof(int));

  PDM_MPI_Alltoall(face_to_send_n,
                   1,
                   PDM_MPI_INT,
                   face_to_recv_n,
                   1,
                   PDM_MPI_INT,
                   ppart->comm);

  int *face_to_recv_idx = PDM_array_new_idx_from_sizes_int(face_to_recv_n, n_rank);

  PDM_g_num_t *face_to_recv =
    (PDM_g_num_t *) malloc(face_to_recv_idx[n_rank]*sizeof(PDM_g_num_t));

  PDM_MPI_Alltoallv(face_to_send,
                    face_to_send_n,
                    face_to_send_idx,
                    PDM__PDM_MPI_G_NUM,
                    face_to_recv,
                    face_to_recv_n,
                    face_to_recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    ppart->comm);

  int n_recv_face = face_to_recv_idx[n_rank]/n_data;

  /*
   * Rename
   */

  int *cell_to_send_idx = face_to_recv_idx;
  int *cell_to_send_n   = face_to_recv_n;

  PDM_g_num_t *cell_to_send =
    (PDM_g_num_t *) malloc(face_to_recv_idx[n_rank]*sizeof(PDM_g_num_t));

  int         *cell_to_recv_idx = face_to_send_idx;
  PDM_g_num_t *cell_to_recv    = face_to_send;
  int         *cell_to_recv_n   = face_to_send_n;

  /*
   * Buid ppart->dface_cell
   */

  int have_dface_cell = 0;

  if (ppart->dface_cell != NULL) {
    have_dface_cell = 1;
  }

  if (!have_dface_cell) {
    ppart->dface_cell = PDM_array_const_gnum(2*ppart->dn_face, 0);

    for (int i = 0; i < n_recv_face; i++) {
      PDM_g_num_t  gface = face_to_recv[n_data*i  ];                    // Get global numbering
      PDM_g_num_t  gcell = face_to_recv[n_data*i+1];                    // Get global numbering
      PDM_g_num_t  _lface = gface - ppart->dface_proc[i_rank];
      int          lface = (int) _lface; // Switch to local numbering

      if (ppart->dface_cell[2*lface] == 0)
        ppart->dface_cell[2*lface] = gcell;
      else if (ppart->dface_cell[2*lface + 1] == 0)
        ppart->dface_cell[2*lface + 1] = gcell;
      else {
        printf("[%i] Failed on %i with gcell = "PDM_FMT_G_NUM" and gface = "PDM_FMT_G_NUM" --> "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" \n", i, lface, gcell, gface, ppart->dface_cell[2*lface], ppart->dface_cell[2*lface+1]);

        double xf = 0.;
        double yf = 0.;
        double zf = 0.;
        double count = 0.;
        for(int idx_vtx = ppart->_dface_vtx_idx[lface]; idx_vtx < ppart->_dface_vtx_idx[lface+1]; ++idx_vtx ){
          int ivtx = ppart->_dface_vtx[idx_vtx]-1;
          xf += ppart->_dvtx_coord[3*ivtx  ];
          yf += ppart->_dvtx_coord[3*ivtx+1];
          zf += ppart->_dvtx_coord[3*ivtx+2];
          count += 1;
        }
        xf = xf / count;
        yf = yf / count;
        zf = zf / count;
        printf("Face barycenter = %12.5e  %12.5e  %12.5e %12.5e \n", xf, yf, zf, count);

        PDM_printf("PPART internal error : Face already defined in ppart->dface_cell connectivity\n");
        exit(1);
      }
    }
    ppart->_dface_cell = ppart->dface_cell;
  }

  /*
   * Exchange cell neighbour
   */

  for (int i = 0; i < n_recv_face; i++) {
    PDM_g_num_t  gface  = face_to_recv[n_data*i  ];                    // Get global numbering
    PDM_g_num_t  gcell1 = face_to_recv[n_data*i+1];                    // Get global numbering
    PDM_g_num_t _lface = gface - ppart->dface_proc[i_rank]; // Switch to local numbering
    int          lface = (int) _lface;
    PDM_g_num_t gcell2;

    if (ppart->dface_cell[2*lface] == gcell1)
      gcell2 = PDM_ABS (ppart->dface_cell[2*lface + 1]);
    else if (ppart->dface_cell[2*lface + 1] == gcell1)
      gcell2 = PDM_ABS (ppart->dface_cell[2*lface]);
    else {
      PDM_printf("PPART internal error : Problem in dual grah building "
              PDM_FMT_G_NUM" "
              PDM_FMT_G_NUM" "
              PDM_FMT_G_NUM" \n",
             ppart->dface_cell[2*lface ], ppart->dface_cell[2*lface + 1], gcell1);
      exit(1);
    }

    cell_to_send[n_data*i    ] = gcell1;
    cell_to_send[n_data*i + 1] = gcell2;
  }

  free(face_to_recv);

  PDM_MPI_Alltoallv(cell_to_send,
                    cell_to_send_n,
                    cell_to_send_idx,
                    PDM__PDM_MPI_G_NUM,
                    cell_to_recv,
                    cell_to_recv_n,
                    cell_to_recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    ppart->comm);

  int n_recv_pair = cell_to_recv_idx[n_rank]/n_data;

  /*
   * Allocate dual graph
   */

  ppart->ddual_graph_idx = (PDM_g_num_t *) malloc((1+ppart->dn_cell) * sizeof(PDM_g_num_t));
  int *n_neighbour      = PDM_array_zeros_int(ppart->dn_cell);

  ppart->ddual_graph_idx[0] = 0;

  ppart->ddual_graph = PDM_array_const_gnum(ppart->_dcell_face_idx[ppart->dn_cell], -1);

  /*
   * Build dual graph
   */

  for (int i = 0; i < n_recv_pair; i++) {
    PDM_g_num_t   gcel1  = cell_to_recv[n_data*i];              // global numbering
    PDM_g_num_t  _lcel1  = gcel1 - ppart->dcell_proc[i_rank]; // local numbering
    int           lcel1  = (int) (_lcel1);
    PDM_g_num_t   gcel2  = cell_to_recv[n_data*i+1];            // global numbering

    /*
     * Search if cel2 is already stored (To optimize for polyhedra (lot of neighbours) ?)
     */

    if (gcel2 > 0) {

      int k;
      for (k = ppart->_dcell_face_idx[lcel1];
           k < ppart->_dcell_face_idx[lcel1] + n_neighbour[lcel1]; k++) {
        if (ppart->ddual_graph[k] == gcel2 - 1)
          break;
      }

      if (k == ppart->_dcell_face_idx[lcel1] + n_neighbour[lcel1]) {
        ppart->ddual_graph[ppart->_dcell_face_idx[lcel1] + n_neighbour[lcel1]] = gcel2 - 1;
        n_neighbour[lcel1] += 1;
      }
    }
  }

  /*
   * Compress dual graph
   */

  int k = 0;
  int k1 = 0;
  while (k < ppart->_dcell_face_idx[ppart->dn_cell]) {
    if (ppart->ddual_graph[k] >= 0) {
      ppart->ddual_graph[k1] = ppart->ddual_graph[k];
      k++;
      k1++;
    }
    else
      k++;
  }

  /*
   * Reallocate to free unused memory
   */

  ppart->ddual_graph = realloc(ppart->ddual_graph, k1 * sizeof(PDM_g_num_t));
  PDM_array_idx_from_sizes_gnum(n_neighbour,  ppart->dn_cell, ppart->ddual_graph_idx);

  /* Verifier tous les tableaux ..... */

  free(cell_to_recv);
  free(cell_to_recv_idx);
  free(cell_to_recv_n);
  free(cell_to_send);
  free(cell_to_send_idx);
  free(cell_to_send_n);
  free(n_neighbour);
}


/**
 *
 * \brief Splits the graph
 *
 * \param [in]  ppart     ppart object
 * \param [out] cell_part  Cell partitioning (size : dn_cell)
 *
 */

static void
_split
(
 _PDM_part_t  *ppart,
 int          *cell_part
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(ppart->comm, &i_rank);
  PDM_MPI_Comm_size(ppart->comm, &n_rank);

  PDM_array_reset_int(cell_part, ppart->dn_cell, 0);


  switch (ppart->split_method) {
  case PDM_PART_SPLIT_PARMETIS:
    {
#ifdef PDM_HAVE_PARMETIS

      /*
       * Define metis properties
       */

      int wgtflag    = 0;
      int numflag    = 0;        /* C or Fortran numbering (C = 0) */
      int edgecut;
      int ncon       = 1;

      double *ubvec = (double *) malloc(ncon * sizeof(double));
      for (int i = 0; i < ncon; i++) {
        ubvec[i] = 1.05;
      }

      double *tpwgts = (double *) malloc(ncon * ppart->tn_part * sizeof(double));

      for (int i = 0; i < ncon * ppart->tn_part; i++) {
        tpwgts[i] = (double) (1./ppart->tn_part);
      }

      /*
       * Call metis
       */

      PDM_g_num_t *_dcell_proc = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));

      for (int i = 0; i < n_rank + 1; i++) {
        _dcell_proc[i] = ppart->dcell_proc[i] - 1;
      }

      PDM_ParMETIS_V3_PartKway (_dcell_proc,
                                ppart->ddual_graph_idx,
                                ppart->ddual_graph,
                                (int *) ppart->_dcell_weight,
                                NULL,
                                &wgtflag,
                                &numflag,
                                &ncon,
                                &ppart->tn_part,
                                tpwgts,
                                ubvec,
                                &edgecut,
                                cell_part,
                                ppart->comm);
      free(ubvec);
      free(tpwgts);
      free(_dcell_proc);

#else
      if(i_rank == 0) {
        PDM_printf("PPART error : ParMETIS unavailable\n");
        exit(1);
      }
#endif
      break;
    }
  case PDM_PART_SPLIT_PTSCOTCH:
    {
#ifdef PDM_HAVE_PTSCOTCH
      int check = 0;
      // printf("chech : %i \n", check);
      int *edgeWeight = NULL;

      // int* dual_graph_n = (int*) malloc( sizeof(int) * ppart->dn_cell);
      // for(int i_entity = 0; i_entity < ppart->dn_cell; ++i_entity) {
      //   dual_graph_n[i_entity] = ppart->ddual_graph_idx[i_entity+1] - ppart->ddual_graph_idx[i_entity];
      // }
      // log_trace("ppart->ddual_graph_idx[ppart->dn_cell]:: %i \n ", ppart->ddual_graph_idx[ppart->dn_cell]);

      // PDM_para_graph_compress_connectivity(ppart->dn_cell, ppart->ddual_graph_idx, dual_graph_n, ppart->ddual_graph);
      // PDM_log_trace_array_long(ppart->ddual_graph_idx, ppart->dn_cell+1, "ddual_graph_idx:: ");
      // log_trace("------------------------- ");
      // PDM_log_trace_array_long(ppart->ddual_graph, ppart->ddual_graph_idx[ppart->dn_cell], "ddual_graph:: ");

      // log_trace("After ppart->ddual_graph_idx[ppart->dn_cell]:: %i \n ", ppart->ddual_graph_idx[ppart->dn_cell]);

      PDM_SCOTCH_dpart (ppart->dn_cell,
                        ppart->ddual_graph_idx,
                        ppart->ddual_graph,
                        ppart->_dcell_weight,
                        edgeWeight,
                        check,
                        ppart->comm,
                        ppart->tn_part,
                        cell_part);
      // free(dual_graph_n);

#else
      if(i_rank == 0) {
        PDM_printf("PPART error : PT-Scotch unavailable\n");
        exit(1);
      }
#endif
      break;
    }
  case PDM_PART_SPLIT_HILBERT:
    {

      for (int i = 0; i < n_rank + 1; ++i) {
        ppart->dface_proc[i] -= 1;
        ppart->dvtx_proc[i]  -= 1;
      }
      PDM_part_geom (PDM_PART_GEOM_HILBERT,
                     ppart->n_part,
                     ppart->comm,
                     ppart->dn_cell,
                     ppart->_dcell_face_idx,
                     ppart->_dcell_face,
                     ppart->_dcell_weight,
                     ppart->_dface_vtx_idx,
                     ppart->_dface_vtx,
                     ppart->dface_proc,
                     ppart->_dvtx_coord,
                     ppart->dvtx_proc,
                     cell_part);
      for (int i = 0; i < n_rank + 1; ++i) {
        ppart->dface_proc[i] += 1;
        ppart->dvtx_proc[i]  += 1;
      }
      break;
    }
  default:
    if(i_rank == 0) {
      PDM_printf("PPART error : '%i' unknown partioning choice\n", ppart->split_method);
      exit(1);
    }
  }
}


/**
 *
 * \brief Distributes cell arrays
 *
 * \param [in]  ppart      ppart object
 * \param [in] cell_part  Cell partitioning (size : 3*dn_cell)
 *
 */

static void
_distrib_cell
(
 _PDM_part_t  *ppart,
 int          *cell_part
)
{

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(ppart->comm, &i_rank);
  PDM_MPI_Comm_size(ppart->comm, &n_rank);

  /*
   *  For each cell face_to_send contains :
   *     - le numero de partition local
   *     - le numero global de l'element
   *     - le nombre de faces
   *     - la liste des faces en numerotation globale
   * on envoie aussi en parallele le tableau send_nbfac qui donne
   * le nombre de faces de chaque element envoye
   * ainsi que son numero de partition local (allant de 0 a nbMeshparProc)
   *
   */

  /* 1ere boucle pour compter le nombre d'elements qu'on envoie a chaque proc */

  int *face_to_send_idx = PDM_array_zeros_int(n_rank + 1);

  int n_data = 3; /* Num cell, Partition locale, nbFac */
  if (ppart->_dcell_tag != NULL)
    n_data += 1;

  for (int i = 0; i < ppart->dn_cell; i++) {
    int _part = (int) cell_part[i];
    int rank_to_send = ppart->gpart_to_lproc_part[2*_part];
    int nbfac      = ppart->_dcell_face_idx[i+1] - ppart->_dcell_face_idx[i];
    face_to_send_idx[rank_to_send+1]  += n_data + nbfac; /* Num cell,
                                                            Partition locale
                                                            nbFac,
                                                            liste des faces */
  }

  face_to_send_idx[0] = 0;
  PDM_array_accumulate_int(face_to_send_idx, n_rank+1);

  int         *face_to_send_n = PDM_array_zeros_int(n_rank);
  PDM_g_num_t *face_to_send  =
    (PDM_g_num_t *) malloc(face_to_send_idx[n_rank] * sizeof(PDM_g_num_t));


  /* 2nde boucle pour remplir le tableau a envoyer via alltoallv */

  for (int i = 0; i < ppart->dn_cell; i++) {
    int _part = (int) cell_part[i];
    int rank_to_send = ppart->gpart_to_lproc_part[2*_part    ];
    int lpart        = ppart->gpart_to_lproc_part[2*_part + 1];
    int nbfac        = ppart->_dcell_face_idx[i+1] - ppart->_dcell_face_idx[i];

    int place = face_to_send_idx[rank_to_send] + face_to_send_n[rank_to_send];

    face_to_send[place++] = lpart;  /* Partition locale */
    face_to_send_n[rank_to_send] += 1;

    face_to_send[place++] = ppart->dcell_proc[i_rank] + i;  /* Numero global de l'elt*/
    face_to_send_n[rank_to_send] += 1;

    face_to_send[place++] = nbfac;  /* Nombre de faces*/
    face_to_send_n[rank_to_send] += 1;

    for (int j = ppart->_dcell_face_idx[i]; j < ppart->_dcell_face_idx[i+1]; j++) {
      face_to_send[place++] = PDM_ABS(ppart->_dcell_face[j]);   /* Numero global de ses faces */
      face_to_send_n[rank_to_send] += 1;
    }

    if (ppart->_dcell_tag != NULL) {
      face_to_send[place++] = ppart->_dcell_tag[i];  /* Tag des cellules si elles existent */
      face_to_send_n[rank_to_send] += 1;
    }
  }

  PDM_g_num_t *face_to_recv     = NULL;
  int         *face_to_recv_n   = (int *) malloc(n_rank * sizeof(int));
  int         *face_to_recv_idx = (int *) malloc((n_rank + 1) * sizeof(int));

  _alltoall(face_to_send,
            face_to_send_n,
            face_to_send_idx,
            (void **) &face_to_recv,
            face_to_recv_n,
            face_to_recv_idx,
            PDM__PDM_MPI_G_NUM,
            sizeof(PDM_g_num_t),
            ppart->comm);

  int lface_to_recv = face_to_recv_idx[n_rank];

  free(face_to_send);
  free(face_to_send_n);
  free(face_to_send_idx);
  free(face_to_recv_n);
  free(face_to_recv_idx);

  /* Complete partitions */
  for (int i = 0; i < ppart->n_part; i++) {
    if (ppart->mesh_parts[i] == NULL){
      ppart->mesh_parts[i] = _part_create();
    }
    _part_t *mesh_part  = ppart->mesh_parts[i];
    mesh_part->n_vtx             = 0;
    mesh_part->n_face            = 0;
    mesh_part->n_cell            = 0;
    mesh_part->n_face_part_bound = 0;
  }

  /* First loop for counting */

  int k = 0;
  while (k < lface_to_recv) {

    _part_t *mesh_part      = ppart->mesh_parts[face_to_recv[k++]];
    k += 1;
    int          n_cell_face = (int) face_to_recv[k++];

    k += n_cell_face;
    if (ppart->_dcell_tag != NULL)
      k += 1;

    mesh_part->n_cell += 1;
    mesh_part->n_face += n_cell_face;  /* Utilisation temporaire de n_face */
  }

  /* Allocates arrays */

  for (int i = 0; i < ppart->n_part; i++) {

    _part_t *mesh_part  = ppart->mesh_parts[i];

    mesh_part->cell_face_idx    = (int *)          malloc((mesh_part->n_cell + 1) * sizeof(int));
    mesh_part->cell_face_idx[0] = 0;
    mesh_part->gcell_face       = (PDM_g_num_t *) malloc(mesh_part->n_face * sizeof(PDM_g_num_t));
    mesh_part->cell_ln_to_gn    = (PDM_g_num_t *) malloc(mesh_part->n_cell * sizeof(PDM_g_num_t));
    if (ppart->_dcell_tag != NULL)
      mesh_part->cell_tag      = (int *)          malloc(mesh_part->n_cell * sizeof(int));

    mesh_part->n_cell          = 0; /* reset temporary */

  }

  /* Second loop to complete arrays */

  k = 0;
  while (k < lface_to_recv) {

    _part_t *mesh_part  = ppart->mesh_parts[face_to_recv[k++]];

    PDM_g_num_t gn_cell    =       face_to_recv[k++];
    mesh_part->cell_ln_to_gn[mesh_part->n_cell] = gn_cell;

    int          n_cell_face = (int) face_to_recv[k++];
    int idx = mesh_part->cell_face_idx[mesh_part->n_cell];
    mesh_part->cell_face_idx[mesh_part->n_cell + 1] = idx + n_cell_face;

    for (int i = 0; i < n_cell_face; i++)
      mesh_part->gcell_face[idx + i] = face_to_recv[k++];

    if (ppart->_dcell_tag != NULL) {
      int tag = (int) face_to_recv[k++];
      mesh_part->cell_tag[mesh_part->n_cell] = tag;
    }

    mesh_part->n_cell += 1;
  }

  free(face_to_recv);

  /* Face local numbering */

  for (int i = 0; i < ppart->n_part; i++) {

    _part_t *mesh_part  = ppart->mesh_parts[i];

    int *initial_idx     = (int *) malloc(mesh_part->n_face * sizeof(int));
    mesh_part->cell_face  = (int *) malloc(mesh_part->n_face * sizeof(int));

    /* Map on gcell_face */

    mesh_part->face_ln_to_gn = mesh_part->gcell_face;

    for (int k1 = 0; k1 < mesh_part->n_face; k1++) {
      initial_idx[k1] = k1;
    }

    /* Sort face_ln_to_gn */

    if (1 == 0) {
      PDM_printf("mesh_part->n_face 1 : %i\n", mesh_part->n_face);
      PDM_printf("mesh_part->face_ln_to_gn 1 : ");
      for (int i1 = 0; i1 < mesh_part->n_face; i1++)
        PDM_printf(" "PDM_FMT_G_NUM, mesh_part->face_ln_to_gn[i1]);
      PDM_printf("\n");
    }

    PDM_quick_sort_long2(mesh_part->face_ln_to_gn, /* tableau a trier */
                         0,                        /* premier elt     */
                         mesh_part->n_face - 1,    /* dernier elt     */
                         initial_idx);

    /* Remove duplicate faces and build local cell face connectivity*/

    int n_dupl = 0;
    int k2 = 0;
    int k_compress = 0;

    while (k2 < mesh_part->n_face) {
      PDM_g_num_t iface = mesh_part->face_ln_to_gn[k2];
      mesh_part->face_ln_to_gn[k_compress]   = iface;
      mesh_part->cell_face[initial_idx[k2]] = k_compress + 1;
      k2 += 1;
      while (k2 < mesh_part->n_face) {
        if (mesh_part->face_ln_to_gn[k2] == iface) {
          mesh_part->cell_face[initial_idx[k2]] = k_compress + 1;
          k2 += 1;
          n_dupl += 1;
        }
        else
          break;
      }
      k_compress += 1;
    }
    PDM_UNUSED(n_dupl);

    mesh_part->n_face = k_compress;

    mesh_part->face_ln_to_gn = (PDM_g_num_t *) realloc(mesh_part->face_ln_to_gn,
                                                    mesh_part->n_face * sizeof(PDM_g_num_t));

    if (1 == 0) {
      PDM_printf("mesh_part->n_cell : %i\n", mesh_part->n_cell);

      PDM_printf("mesh_part->cell_ln_to_gn : ");
      for (int i1 = 0; i1 < mesh_part->n_cell; i1++)
        PDM_printf(" "PDM_FMT_G_NUM, mesh_part->cell_ln_to_gn[i1]);
      PDM_printf("\n");

      PDM_printf("mesh_part->cell_face : \n");
      for (int i1 = 0; i1 < mesh_part->n_cell; i1++) {
        for (int j = mesh_part->cell_face_idx[i1]; j < mesh_part->cell_face_idx[i1+1]; j++)
          PDM_printf(" %i", mesh_part->cell_face[j]);
        PDM_printf("\n");
      }

      PDM_printf("mesh_part->n_face : %i\n", mesh_part->n_face);

      PDM_printf("mesh_part->face_ln_to_gn : ");
      for (int i1 = 0; i1 < mesh_part->n_face; i1++)
        PDM_printf(" "PDM_FMT_G_NUM, mesh_part->face_ln_to_gn[i1]);
      PDM_printf("\n");
    }

    /* reordering cells and faces */
    mesh_part->new_to_old_order_cell = (int *) malloc (sizeof(int) * mesh_part->n_cell);
    for (int i1 = 0; i1 < mesh_part->n_cell; i1++){
      mesh_part->new_to_old_order_cell[i1] = i1;
    }
    mesh_part->new_to_old_order_face = (int *) malloc (sizeof(int) * mesh_part->n_face);
    for (int i1 = 0; i1 < mesh_part->n_face; i1++){
      mesh_part->new_to_old_order_face[i1] = i1;
    }

    /* Free */

    free(initial_idx);
    mesh_part->gcell_face = NULL;
  }

}


/**
 *
 * \brief Distributes face arrays
 *
 * \param [in]  ppart      ppart object
 *
 */

static void
_distrib_face
(
 _PDM_part_t      *ppart
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(ppart->comm, &i_rank);
  PDM_MPI_Comm_size(ppart->comm, &n_rank);

  const int n_data      = 1;
  int       n_data_face = 2;
  if (ppart->_dface_tag != NULL)
    n_data_face += 1;

  int          *face_to_send_idx = (int *) malloc((n_rank + 1) * sizeof(int));
  int          *face_to_send_n   = (int *) malloc(n_rank * sizeof(int));
  PDM_g_num_t *face_to_send    = NULL;

  int          *requested_face_n   = (int *) malloc(n_rank * sizeof(int));
  int          *requested_face_idx = (int *) malloc((n_rank + 1) * sizeof(int));

  for (int i_part = 0; i_part < ppart->mn_part; i_part++) {

    _part_t *mesh_part  = NULL;
    int *all_to_all_n_to_ln = NULL;

    PDM_array_reset_int(face_to_send_idx, n_rank+1, 0);
    PDM_array_reset_int(face_to_send_n, n_rank, 0);

    face_to_send = NULL;

    if (i_part < ppart->n_part) {

      mesh_part  = ppart->mesh_parts[i_part];
      all_to_all_n_to_ln = (int *) malloc(mesh_part->n_face * sizeof(int));

      /*
       *  Processes exchange list of faces which they want receive information
       */

      for (int i = 0; i < mesh_part->n_face; i++) {
        PDM_g_num_t iface = mesh_part->face_ln_to_gn[i];
        int          found_rank = PDM_search_rank(iface, ppart->dface_proc, 0, n_rank);
        face_to_send_idx[found_rank+1] += n_data;

      }

      for (int i = 0; i < n_rank; i++) {
        face_to_send_idx[i+1] += face_to_send_idx[i] ;
      }

      face_to_send = (PDM_g_num_t *) malloc(face_to_send_idx[n_rank] * sizeof(PDM_g_num_t));

      for (int i = 0; i < mesh_part->n_face; i++) {

        PDM_g_num_t iface = mesh_part->face_ln_to_gn[i];
        int found_rank = PDM_search_rank(iface, ppart->dface_proc, 0, n_rank);

        int idx = face_to_send_idx[found_rank] + face_to_send_n[found_rank];

        all_to_all_n_to_ln[idx/n_data] = i;

        face_to_send[idx++]   = iface;        /* Face global numbering */
        face_to_send_n[found_rank] += n_data;
      }
    }

    PDM_g_num_t *requested_face    = NULL;

    _alltoall(face_to_send,
              face_to_send_n,
              face_to_send_idx,
              (void **) &requested_face,
              requested_face_n,
              requested_face_idx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    if (face_to_send != NULL)
      free(face_to_send);

    /*
     *  Processes exchange information about requested faces
     *  For each face, information contains :
     *     - tag (if ppart->_dface_tag != NULL)
     *     - Number of vertices
     *     - Vertices
     *
     */

    int *sface_info_idx = face_to_send_idx;
    int *sface_info_n   = face_to_send_n;

    PDM_array_reset_int(sface_info_idx, n_rank+1, 0);
    PDM_array_reset_int(sface_info_n, n_rank, 0);

    for (int i = 0; i < n_rank; i++) {
      for (int k = requested_face_idx[i]; k < requested_face_idx[i+1]; k+=n_data) {
        PDM_g_num_t gface     = requested_face[k];
        PDM_g_num_t _lface    = gface - ppart->dface_proc[i_rank];
        int          lface     = (int) _lface;
        int          nb_vtx_face = (int) (ppart->_dface_vtx_idx[lface+1]
                                      - ppart->_dface_vtx_idx[lface]);
        sface_info_idx[i+1] += n_data_face + nb_vtx_face;
      }
    }

    for (int i = 0; i < n_rank; i++) {
      sface_info_idx[i+1] += sface_info_idx[i] ;
    }

    PDM_g_num_t *sface_info = (PDM_g_num_t *)
      malloc(sface_info_idx[n_rank] * sizeof(PDM_g_num_t));

    for (int i = 0; i < n_rank; i++) {
      for (int k = requested_face_idx[i]; k < requested_face_idx[i+1]; k+=n_data) {
        PDM_g_num_t gface     = requested_face[k];
        PDM_g_num_t _lface    = gface - ppart->dface_proc[i_rank];
        int          lface    = (int) _lface;
        int          nb_vtx_face = (int) (ppart->_dface_vtx_idx[lface+1]
                                      - ppart->_dface_vtx_idx[lface]);

        int idx = sface_info_idx[i] + sface_info_n[i];

        if (ppart->_dface_tag != NULL) {
          sface_info[idx++] = ppart->_dface_tag[lface];   /* Tag de la face */
          sface_info_n[i] += 1;
        }

        sface_info[idx++] = nb_vtx_face;                   /* Number of vertices */
        sface_info_n[i] += 1;

        for(int j = ppart->_dface_vtx_idx[lface]; j < ppart->_dface_vtx_idx[lface+1]; j++) {
          sface_info[idx++] =  ppart->_dface_vtx[j];  /*numero global du sommet qui compose la face*/
          sface_info_n[i] += 1;
        }
      }
    }

    free(requested_face);

    PDM_g_num_t  *rface_info    = NULL;
    int          *rface_info_n   = requested_face_n;
    int          *rface_info_idx = requested_face_idx;

    _alltoall(sface_info,
              sface_info_n,
              sface_info_idx,
              (void **) &rface_info,
              rface_info_n,
              rface_info_idx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    free(sface_info);
    sface_info = NULL;

    if (i_part < ppart->n_part) {

      /* Complete face_tag, face_vtx_idx gface_vtx */

      if (ppart->_dface_tag != NULL)
        mesh_part->face_tag  = (int *) malloc(mesh_part->n_face * sizeof(int));
      mesh_part->face_vtx_idx = (int *) malloc((mesh_part->n_face + 1) * sizeof(int));

      int k = 0;
      for (int i = 0; i < mesh_part->n_face; i++) {
        if (ppart->_dface_tag != NULL)
          mesh_part->face_tag[all_to_all_n_to_ln[i]] = (int) rface_info[k++];

        int n_vtx = (int) rface_info[k++];
        mesh_part->face_vtx_idx[all_to_all_n_to_ln[i]+1] = n_vtx;

        k += n_vtx;
      }

      mesh_part->face_vtx_idx[0] = 0;
      PDM_array_accumulate_int(mesh_part->face_vtx_idx, mesh_part->n_face+1);

      mesh_part->gface_vtx =
        (PDM_g_num_t *) malloc(mesh_part->face_vtx_idx[mesh_part->n_face] * sizeof(PDM_g_num_t));

      k = 0;
      for (int i = 0; i < mesh_part->n_face; i++) {
        if (ppart->_dface_tag != NULL)
          k += 1;

        int n_vtx = (int) rface_info[k++];
        int idx = mesh_part->face_vtx_idx[all_to_all_n_to_ln[i]];
        for (int j = 0; j < n_vtx; j++)
          mesh_part->gface_vtx[idx + j] = rface_info[k++];
      }

      if (rface_info != NULL)
        free(rface_info);
      rface_info = NULL;

      /* Vertex local numbering vtx_ln_to_gn */

      int *initial_idx   = (int *) malloc(mesh_part->face_vtx_idx[mesh_part->n_face] * sizeof(int));
      mesh_part->face_vtx = (int *) malloc(mesh_part->face_vtx_idx[mesh_part->n_face] * sizeof(int));

      /* Map on gface_vtx */

      mesh_part->vtx_ln_to_gn = mesh_part->gface_vtx;

      for (int k1 = 0; k1 <  mesh_part->face_vtx_idx[mesh_part->n_face]; k1++) {
        initial_idx[k1] = k1;
      }

      /* Sort face_ln_to_gn */

      PDM_quick_sort_long2(mesh_part->vtx_ln_to_gn,                        /* Array to sort */
                           0,                                              /* First face    */
                           mesh_part->face_vtx_idx[mesh_part->n_face] - 1, /* Latest face   */
                           initial_idx);

      /* Remove duplicate Vertex and build local face vertex connectivity*/

      int n_dupl = 0;
      int k2 = 0;
      int k_compress = 0;

      while (k2 < mesh_part->face_vtx_idx[mesh_part->n_face]) {
        PDM_g_num_t iVtx = mesh_part->vtx_ln_to_gn[k2];
        mesh_part->vtx_ln_to_gn[k_compress]   = iVtx;
        mesh_part->face_vtx[initial_idx[k2]] = k_compress + 1;
        k2 += 1;
        while (k2 < mesh_part->face_vtx_idx[mesh_part->n_face]) {
          if (mesh_part->vtx_ln_to_gn[k2] == iVtx) {
            mesh_part->face_vtx[initial_idx[k2]] = k_compress + 1;
            k2 += 1;
            n_dupl += 1;
          }
          else
            break;
        }
        k_compress += 1;
      }
      PDM_UNUSED(n_dupl);


      mesh_part->n_vtx = k_compress;

      mesh_part->vtx_ln_to_gn =
        (PDM_g_num_t *) realloc(mesh_part->vtx_ln_to_gn, mesh_part->n_vtx * sizeof(PDM_g_num_t));

      /* Free */

      free(initial_idx);
      mesh_part->gface_vtx = NULL;

      if (1 == 0) {
        PDM_printf("mesh_part->n_vtx 1 : %i\n", mesh_part->n_vtx);
        PDM_printf("mesh_part->vtx_ln_to_gn 1 : ");
        for (int i1 = 0; i1 < mesh_part->n_vtx; i1++)
          PDM_printf(" "PDM_FMT_G_NUM, mesh_part->vtx_ln_to_gn[i1]);
        PDM_printf("\n");

        PDM_printf("mesh_part->face_vtx : \n");
        for (int i1 = 0; i1 < mesh_part->n_face; i1++) {
          for (int j = mesh_part->face_vtx_idx[i1]; j < mesh_part->face_vtx_idx[i1+1]; j++)
            PDM_printf(" %i", mesh_part->face_vtx[j]);
          PDM_printf("\n");
        }
      }

    }

    if (rface_info != NULL)
      free(rface_info);
    rface_info = NULL;

    if (all_to_all_n_to_ln != NULL)
      free(all_to_all_n_to_ln);

  } /* For i_part */


  free(face_to_send_n);
  free(face_to_send_idx);
  free(requested_face_n);
  free(requested_face_idx);

}


/**
 *
 * \brief Distributes vertex arrays
 *
 * \param [in]  ppart      ppart object
 *
 */

static void
_distrib_vtx
(
 _PDM_part_t      *ppart
)
{

  const int n_data    = 1;
  int n_data_vtx = 0;
  if (ppart->_dvtx_tag != NULL)
    n_data_vtx += 1;

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(ppart->comm, &i_rank);
  PDM_MPI_Comm_size(ppart->comm, &n_rank);

  int          *vtx_to_send_idx = (int *) malloc((n_rank + 1) * sizeof(int));
  int          *vtx_to_send_n   = (int *) malloc(n_rank * sizeof(int));
  PDM_g_num_t  *vtx_to_send    = NULL;

  int          *requested_vtx_n   = (int *) malloc(n_rank * sizeof(int));
  int          *requested_vtx_idx = (int *) malloc((n_rank + 1) * sizeof(int));

  for (int i_part = 0; i_part < ppart->mn_part; i_part++) {

    _part_t *mesh_part  = NULL;
    int *all_to_all_n_to_ln = NULL;

    PDM_array_reset_int(vtx_to_send_idx, n_rank+1, 0);
    PDM_array_reset_int(vtx_to_send_n, n_rank, 0);

    vtx_to_send = NULL;

    if (i_part < ppart->n_part) {

      mesh_part  = ppart->mesh_parts[i_part];
      all_to_all_n_to_ln = (int *) malloc(mesh_part->n_vtx * sizeof(int));

      /*
       *  Processes exchange list of vtxs which they want receive information
       */

      for (int i = 0; i < mesh_part->n_vtx; i++) {
        PDM_g_num_t iVtx = mesh_part->vtx_ln_to_gn[i];
        int found_rank = PDM_search_rank(iVtx, ppart->dvtx_proc, 0, n_rank);
        vtx_to_send_idx[found_rank+1] += n_data;
      }

      for (int i = 0; i < n_rank; i++) {
        vtx_to_send_idx[i+1] += vtx_to_send_idx[i] ;
      }

      vtx_to_send = (PDM_g_num_t *) malloc(vtx_to_send_idx[n_rank] * sizeof(PDM_g_num_t));

      for (int i = 0; i < mesh_part->n_vtx; i++) {

        PDM_g_num_t iVtx = mesh_part->vtx_ln_to_gn[i];
        int found_rank = PDM_search_rank(iVtx, ppart->dvtx_proc, 0, n_rank);

        int idx = vtx_to_send_idx[found_rank] + vtx_to_send_n[found_rank];

        all_to_all_n_to_ln[idx/n_data] = i;

        vtx_to_send[idx++]   = iVtx;        /* Vtx global numbering */
        vtx_to_send_n[found_rank] += n_data;
      }
    }

    PDM_g_num_t *requested_vtx    = NULL;

    _alltoall(vtx_to_send,
              vtx_to_send_n,
              vtx_to_send_idx,
              (void **) &requested_vtx,
              requested_vtx_n,
              requested_vtx_idx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    if (vtx_to_send != NULL)
      free(vtx_to_send);

    /*
     *  Processes exchange information about requested vtxs
     *  For each vtx, information contains :
     *     - Tag (if ppart->_dvtx_tag != NULL)
     *     - Coordinates
     *
     */

    int *svtx_info_idx = vtx_to_send_idx;
    int *svtx_info_n   = vtx_to_send_n;

    PDM_array_reset_int(svtx_info_idx, n_rank+1, 0);
    PDM_array_reset_int(svtx_info_n, n_rank, 0);

    for (int i = 0; i < n_rank; i++) {
      for (int k = requested_vtx_idx[i]; k < requested_vtx_idx[i+1]; k += n_data) {
        svtx_info_idx[i+1] += n_data_vtx * (int) sizeof(int) + 3 * (int) sizeof(double);
      }
    }

    for (int i = 0; i < n_rank; i++) {
      svtx_info_idx[i+1] += svtx_info_idx[i];
    }

    unsigned char *svtx_info = (unsigned char *)
      malloc(svtx_info_idx[n_rank] * sizeof(unsigned char));

    for (int i = 0; i < n_rank; i++) {
      for (int k = requested_vtx_idx[i]; k < requested_vtx_idx[i+1]; k+=n_data) {
        PDM_g_num_t  gvtx     = requested_vtx[k];
        PDM_g_num_t _lvtx     = gvtx - ppart->dvtx_proc[i_rank];
        int          lvtx     = (int) _lvtx;

        int idx = svtx_info_idx[i] + svtx_info_n[i];

        if (ppart->_dvtx_tag != NULL) {
          int *_i_svtx_info = (int *) (svtx_info + idx);
          *_i_svtx_info     = ppart->_dvtx_tag[lvtx];   /* Tag de la vtx */
          svtx_info_n[i] += sizeof(int);
          idx += sizeof(int);
        }

        double *_d_svtx_info = (double *) (svtx_info + idx);
        _d_svtx_info[0] = ppart->_dvtx_coord[3*lvtx    ];
        _d_svtx_info[1] = ppart->_dvtx_coord[3*lvtx + 1];
        _d_svtx_info[2] = ppart->_dvtx_coord[3*lvtx + 2];
        svtx_info_n[i] += 3 * sizeof(double);

      }
    }

    free(requested_vtx);

    unsigned char *rvtx_info     = NULL;
    int           *rvtx_info_n   = requested_vtx_n;
    int           *rvtx_info_idx = requested_vtx_idx;

    _alltoall(svtx_info,
              svtx_info_n,
              svtx_info_idx,
              (void **) &rvtx_info,
              rvtx_info_n,
              rvtx_info_idx,
              PDM_MPI_UNSIGNED_CHAR,
              sizeof(PDM_MPI_UNSIGNED_CHAR),
              ppart->comm);

    if (svtx_info != NULL)
      free(svtx_info);

    if (i_part < ppart->n_part) {

      /* Complete vtx_tag, vtx */

      if (ppart->_dvtx_tag != NULL)
        mesh_part->vtx_tag  = (int *) malloc(mesh_part->n_vtx * sizeof(int));
      mesh_part->vtx = (double *) malloc(3 * mesh_part->n_vtx * sizeof(double));

      // mesh_part->new_to_old_order_vtx = (int *) malloc (sizeof(int) * mesh_part->n_vtx);
      // for (int i1 = 0; i1 < mesh_part->n_vtx; i1++){
      //   mesh_part->new_to_old_order_vtx[i1] = i1;
      // }

      int k = 0;
      for (int i = 0; i < mesh_part->n_vtx; i++) {
        if (ppart->_dvtx_tag != NULL) {
          int *_i_rVtxInfo = (int *) (rvtx_info + k);
          mesh_part->vtx_tag[all_to_all_n_to_ln[i]] = *_i_rVtxInfo;
          k += sizeof(int);
        }

        double *_d_rVtxInfo = (double *) (rvtx_info + k);
        mesh_part->vtx[3*all_to_all_n_to_ln[i]    ] = _d_rVtxInfo[0];
        mesh_part->vtx[3*all_to_all_n_to_ln[i] + 1] = _d_rVtxInfo[1];
        mesh_part->vtx[3*all_to_all_n_to_ln[i] + 2] = _d_rVtxInfo[2];
        k += 3*sizeof(double);
      }

      if (1 == 0) {
        PDM_printf("mesh_part->vtx : \n");
        for (int i1 = 0; i1 < mesh_part->n_vtx; i1++) {
          PDM_printf(" %12.5e %12.5e %12.5e", mesh_part->vtx[3*i1 ], mesh_part->vtx[3*i1+1], mesh_part->vtx[3*i1+2]);
          PDM_printf("\n");
        }
      }

    } /* if i_part */

    if (all_to_all_n_to_ln != NULL)
      free(all_to_all_n_to_ln);

    if (rvtx_info != NULL)
      free(rvtx_info);
    rvtx_info = NULL;

  } /* For i_part */

  free(vtx_to_send_n);
  free(vtx_to_send_idx);
  free(requested_vtx_n);
  free(requested_vtx_idx);

}


/**
 *
 * \brief Builds face-cell connectivity
 *
 * \param [in]  ppart      ppart object
 *
 */

static void
_build_faceCell
(
 _PDM_part_t  *ppart
)
{
  for (int i_part = 0; i_part < ppart->n_part; i_part++) {
    _part_t *mesh_part  = ppart->mesh_parts[i_part];

    mesh_part->face_cell = PDM_array_zeros_int(2*mesh_part->n_face);

    for (int i = 0; i < mesh_part->n_cell; i++) {
      for (int j = mesh_part->cell_face_idx[i]; j < mesh_part->cell_face_idx[i+1]; j++) {
        int idx = 2 * (PDM_ABS(mesh_part->cell_face[j])-1);
        if (mesh_part->face_cell[idx] == 0)
          mesh_part->face_cell[idx] = i + 1;
        else
          mesh_part->face_cell[idx + 1] = i + 1;
      }
    }
    if (1 == 0) {
      PDM_printf("mesh_part->face_cell : \n");
      for (int i1 = 0; i1 < mesh_part->n_face; i1++) {
        PDM_printf(" %i %i", mesh_part->face_cell[2*i1],  mesh_part->face_cell[2*i1+1]);
        PDM_printf("\n");
      }
    }
  }
}


/**
 *
 * \brief Search partitioning boundary faces
 *
 * \param [in]  ppart      ppart object
 *
 */

static void
_search_part_bound_face
(
 _PDM_part_t *ppart
)
{
  const int n_data  = 4;
  const int n_data2 = 5;

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(ppart->comm, &i_rank);
  PDM_MPI_Comm_size(ppart->comm, &n_rank);

  int          *face_to_send_idx = (int *) malloc((n_rank + 1) * sizeof(int));
  int          *face_to_send_n   = (int *) malloc(n_rank * sizeof(int));
  PDM_g_num_t  *face_to_send     = NULL;

  int          *requested_face_n   = (int *) malloc(n_rank * sizeof(int));
  int          *requested_face_idx = (int *) malloc((n_rank + 1) * sizeof(int));

  int n_data_pb = 6;

  ppart->dpart_bound = PDM_array_const_int(n_data_pb*ppart->dn_face, -1);

  /*
   * First loop on partitions to look for boundary faces
   */

  for (int i_part = 0; i_part < ppart->mn_part; i_part++) {

    _part_t *mesh_part  = NULL;

    PDM_array_reset_int(face_to_send_idx, n_rank+1, 0);
    PDM_array_reset_int(face_to_send_n, n_rank, 0);

    face_to_send = NULL;

    if (i_part < ppart->n_part) {

      mesh_part  = ppart->mesh_parts[i_part];

      /*
       *  Processes exchange list of faces which they want receive information
       */

      int nBoundFace = 0;

      for (int i = 0; i < mesh_part->n_face; i++) {
        int i_cell2 = PDM_ABS (mesh_part->face_cell[2*i + 1]);
        if (i_cell2 == 0) {
          PDM_g_num_t iface = mesh_part->face_ln_to_gn[i];
          int found_rank = PDM_search_rank(iface, ppart->dface_proc, 0, n_rank);
          face_to_send_idx[found_rank+1] += n_data;
          nBoundFace += 1;
        }
      }
      PDM_UNUSED(nBoundFace);

      for (int i = 0; i < n_rank; i++) {
        face_to_send_idx[i+1] += face_to_send_idx[i] ;
      }

      face_to_send = (PDM_g_num_t *) malloc(face_to_send_idx[n_rank] * sizeof(PDM_g_num_t));

      for (int i = 0; i < mesh_part->n_face; i++) {

        int i_cell2 = PDM_ABS (mesh_part->face_cell[2*i + 1]);
        if (i_cell2 == 0) {
          PDM_g_num_t gface = mesh_part->face_ln_to_gn[i];
          int found_rank = PDM_search_rank(gface, ppart->dface_proc, 0, n_rank);

          int idx = face_to_send_idx[found_rank] + face_to_send_n[found_rank];

          face_to_send[idx++]   = gface;        /* Face global numbering */
          face_to_send[idx++]   = i+1;          /* Face local numbering  */
          face_to_send[idx++]   = i_rank;       /* Rank                  */
          face_to_send[idx++]   = i_part;        /* Partition             */
          face_to_send_n[found_rank] += n_data;
        }
      }
    }

    PDM_g_num_t *requested_face    = NULL;

    _alltoall(face_to_send,
              face_to_send_n,
              face_to_send_idx,
              (void **) &requested_face,
              requested_face_n,
              requested_face_idx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    free(face_to_send);
    int n_face = requested_face_idx[n_rank]/n_data;


    int idx = 0;
    for(int i = 0; i < n_face; i++) {
      PDM_g_num_t  gface     = requested_face[idx++];
      PDM_g_num_t  _lface    = gface - ppart->dface_proc[i_rank];
      int          lface     = (int) _lface;
      int          lfaceRank = (int) requested_face[idx++];
      int          faceRank  = (int) requested_face[idx++];
      int          partition = (int) requested_face[idx++];

      int idx2 = 0;
      if (ppart->dpart_bound[n_data_pb * lface] != -1)
        idx2 += n_data_pb/2;

      ppart->dpart_bound[n_data_pb * lface + idx2++] = faceRank;
      ppart->dpart_bound[n_data_pb * lface + idx2++] = lfaceRank;
      ppart->dpart_bound[n_data_pb * lface + idx2  ] = partition;

    }

    free(requested_face);

  } /* mn_part */

  /* Exchange dpart_bound */

  PDM_array_reset_int(face_to_send_idx, n_rank+1, 0);
  PDM_array_reset_int(face_to_send_n, n_rank, 0);

  int idx = 0;
  for(int i = 0; i < ppart->dn_face; i++) {
    int faceRank1  = ppart->dpart_bound[idx++];
    idx += 2;

    int faceRank2  = ppart->dpart_bound[idx++];
    idx += 2;

    if ((faceRank1 > -1) && (faceRank2 > -1)) {
      face_to_send_idx[faceRank1 + 1] += n_data2;
      face_to_send_idx[faceRank2 + 1] += n_data2;
    }
  }

  for (int i = 0; i < n_rank; i++)
    face_to_send_idx[i + 1] += face_to_send_idx[i];

  int *face_to_sendInt = (int *) malloc(face_to_send_idx[n_rank] * sizeof(int));

  idx = 0;
  for(int i = 0; i < ppart->dn_face; i++) {
    int faceRank1  = ppart->dpart_bound[idx++];
    int lfaceRank1 = ppart->dpart_bound[idx++];
    int partition1 = ppart->dpart_bound[idx++];

    int faceRank2  = ppart->dpart_bound[idx++];
    int lfaceRank2 = ppart->dpart_bound[idx++];
    int partition2 = ppart->dpart_bound[idx++];

    if ((faceRank1 > -1) && (faceRank2 > -1)) {
      int idx2 = face_to_send_idx[faceRank1] + face_to_send_n[faceRank1];
      face_to_sendInt[idx2++]   = lfaceRank1;
      face_to_sendInt[idx2++]   = partition1;
      face_to_sendInt[idx2++]   = faceRank2;
      face_to_sendInt[idx2++]   = lfaceRank2;
      face_to_sendInt[idx2++]   = partition2;
      face_to_send_n[faceRank1] += n_data2;

      int idx3 = face_to_send_idx[faceRank2] + face_to_send_n[faceRank2];
      face_to_sendInt[idx3++]   = lfaceRank2;
      face_to_sendInt[idx3++]   = partition2;
      face_to_sendInt[idx3++]   = faceRank1;
      face_to_sendInt[idx3++]   = lfaceRank1;
      face_to_sendInt[idx3++]   = partition1;
      face_to_send_n[faceRank2] += n_data2;

    }
  }

  int *requested_face_int = NULL;

  _alltoall(face_to_sendInt,
            face_to_send_n,
            face_to_send_idx,
            (void **) &requested_face_int,
            requested_face_n,
            requested_face_idx,
            PDM_MPI_INT,
            sizeof(PDM_MPI_INT),
            ppart->comm);

  free(face_to_sendInt);

  /* Complete face_part_bound */

  int n_face_part_bound_rank = requested_face_idx[n_rank]/n_data2;

  idx = 0;
  for (int i = 0; i < n_face_part_bound_rank; i++) {
    idx += 1;
    int partition1 = requested_face_int[idx++];

    idx += 3;

    _part_t *mesh_part  = ppart->mesh_parts[partition1];
    mesh_part->n_face_part_bound += 1;
  }

  int n_data_face_part_bound = 4;

  for (int i = 0; i < ppart->n_part; i++) {
    _part_t *mesh_part  = ppart->mesh_parts[i];
    mesh_part->face_part_bound =
      (int *) malloc(n_data_face_part_bound * mesh_part->n_face_part_bound * sizeof(int));
    mesh_part->n_face_part_bound = 0;

    mesh_part->face_part_bound_proc_idx = PDM_array_zeros_int(n_rank+1);
    mesh_part->face_part_bound_part_idx = PDM_array_zeros_int(ppart->tn_part+1);
  }

  idx = 0;

  for (int i = 0; i < n_face_part_bound_rank; i++) {
    int lfaceRank1 = requested_face_int[idx++];
    int partition1 = requested_face_int[idx++];

    int faceRank2  = requested_face_int[idx++];
    int lfaceRank2 = requested_face_int[idx++];
    int partition2 = requested_face_int[idx++];

    _part_t *mesh_part  = ppart->mesh_parts[partition1];
    mesh_part->face_part_bound_proc_idx[faceRank2+1] += 1;
    mesh_part->face_part_bound_part_idx[ppart->dpart_proc[faceRank2] + partition2 + 1] += 1;
    mesh_part->face_part_bound[n_data_face_part_bound * mesh_part->n_face_part_bound    ] = lfaceRank1;
    mesh_part->face_part_bound[n_data_face_part_bound * mesh_part->n_face_part_bound + 1] = faceRank2;
    mesh_part->face_part_bound[n_data_face_part_bound * mesh_part->n_face_part_bound + 2] = partition2 + 1;
    mesh_part->face_part_bound[n_data_face_part_bound * mesh_part->n_face_part_bound + 3] = lfaceRank2;
    mesh_part->n_face_part_bound += 1;
  }

  for (int i = 0; i < ppart->n_part; i++) {
    _part_t *mesh_part  = ppart->mesh_parts[i];
    for (int j = 1; j < n_rank + 1; j++) {
      mesh_part->face_part_bound_proc_idx[j] = mesh_part->face_part_bound_proc_idx[j] + mesh_part->face_part_bound_proc_idx[j-1];
    }
    for (int j = 1; j <  ppart->tn_part + 1; j++) {
      mesh_part->face_part_bound_part_idx[j] = mesh_part->face_part_bound_part_idx[j] + mesh_part->face_part_bound_part_idx[j-1];
    }
  }

  for (int i = 0; i < ppart->n_part; i++) {
    _part_t *mesh_part  = ppart->mesh_parts[i];

    int *work_array  = (int *) malloc(mesh_part->n_face_part_bound * sizeof(int));
    PDM_g_num_t *work_array2;
    if (sizeof(PDM_g_num_t) == sizeof(int)) {
      work_array2 = (PDM_g_num_t *) work_array;
    }
    else {
      work_array2  = (PDM_g_num_t *) malloc(mesh_part->n_face_part_bound * sizeof(PDM_g_num_t));
    }

    int *ind         = (int *) malloc(mesh_part->n_face_part_bound * sizeof(int));
    int *copy_face_part_bound = (int *) malloc(n_data_face_part_bound * mesh_part->n_face_part_bound * sizeof(int));

    /* Sort by procs */

    int k = 0;
    for (int j = 0; j < mesh_part->n_face_part_bound; j++) {
      ind[j] = k;
      k += 1;
      work_array[j] =  mesh_part->face_part_bound[n_data_face_part_bound * j + 1];
    }

    PDM_quick_sort_int2(work_array,
                    0,
                    mesh_part->n_face_part_bound - 1,
                    ind);

    for (int j = 0; j < n_data_face_part_bound * mesh_part->n_face_part_bound; j++)
      copy_face_part_bound[j] =  mesh_part->face_part_bound[j];

    for (int j = 0; j <  mesh_part->n_face_part_bound; j++) {
      mesh_part->face_part_bound[n_data_face_part_bound * j    ] = copy_face_part_bound[n_data_face_part_bound * ind[j]    ];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 1] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 1];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 2] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 2];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 3] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 3];
    }

    /* Sort by part in procs */

    for (int j = 0; j < n_data_face_part_bound * mesh_part->n_face_part_bound; j++)
      copy_face_part_bound[j] =  mesh_part->face_part_bound[j];

    k = 0;
    for (int j = 0; j < mesh_part->n_face_part_bound; j++) {
      ind[j] = k;
      k += 1;
      work_array[j] =  mesh_part->face_part_bound[n_data_face_part_bound * j + 2];
    }

    for (int j = 0; j < n_rank; j++) {
      PDM_quick_sort_int2(work_array,
                      mesh_part->face_part_bound_proc_idx[j],
                      mesh_part->face_part_bound_proc_idx[j+1]-1,
                      ind);
    }

    for (int j = 0; j <  mesh_part->n_face_part_bound; j++) {
      mesh_part->face_part_bound[n_data_face_part_bound * j    ] = copy_face_part_bound[n_data_face_part_bound * ind[j]    ];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 1] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 1];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 2] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 2];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 3] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 3];
    }

    /* Sort by face absolute number in parts */

    for (int j = 0; j < n_data_face_part_bound * mesh_part->n_face_part_bound; j++)
      copy_face_part_bound[j] =  mesh_part->face_part_bound[j];

    k = 0;
    for (int j = 0; j < mesh_part->n_face_part_bound; j++) {
      ind[j] = k;
      k += 1;
      work_array2[j] =   mesh_part->face_ln_to_gn[mesh_part->face_part_bound[n_data_face_part_bound * j] - 1];
    }

    for (int j = 0; j < ppart->tn_part; j++) {
      PDM_quick_sort_long2(work_array2,
                           mesh_part->face_part_bound_part_idx[j],
                           mesh_part->face_part_bound_part_idx[j+1]-1,
                           ind);
    }

    for (int j = 0; j <  mesh_part->n_face_part_bound; j++) {
      mesh_part->face_part_bound[n_data_face_part_bound * j    ] = copy_face_part_bound[n_data_face_part_bound * ind[j]    ];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 1] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 1];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 2] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 2];
      mesh_part->face_part_bound[n_data_face_part_bound * j + 3] = copy_face_part_bound[n_data_face_part_bound * ind[j] + 3];
    }

    if (sizeof(PDM_g_num_t) != sizeof(int)) {
      free (work_array2);
    }
    free(work_array);
    free(ind);
    free(copy_face_part_bound);

  }

  /* Sort by absolute face in parts */

  if (0 == 1) {
    for (int i = 0; i < ppart->n_part; i++) {
      _part_t *mesh_part  = ppart->mesh_parts[i];
      PDM_printf("[%i] mesh_part->n_face_part_bound : %i\n",i_rank, mesh_part->n_face_part_bound);
      PDM_printf("[%i] mesh_part->face_part_bound : \n", i_rank);
      for (int i1 = 0; i1 < mesh_part->n_face_part_bound; i1++) {
        PDM_printf("[%i] %i %i %i %i", i_rank, mesh_part->face_part_bound[4*i1    ],
               mesh_part->face_part_bound[4*i1 + 1],
               mesh_part->face_part_bound[4*i1 + 2],
               mesh_part->face_part_bound[4*i1 + 3]);
        PDM_printf("\n");
      }
    }
  }

  free(requested_face_int);
  free(ppart->dpart_bound);
  ppart->dpart_bound = NULL;
  free(face_to_send_idx);
  free(face_to_send_n);

  free(requested_face_n);
  free(requested_face_idx);

}


/**
 *
 * \brief Distributes boundaries
 *
 * \param [in]  ppart      ppart object
 *
 */

static void
_distrib_face_groups
(
 _PDM_part_t      *ppart
)
{

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(ppart->comm, &i_rank);
  PDM_MPI_Comm_size(ppart->comm, &n_rank);

  int          *face_to_send_idx = (int *) malloc((n_rank + 1) * sizeof(int));
  int          *face_to_send_n   = (int *) malloc(n_rank * sizeof(int));
  PDM_g_num_t  *face_to_send    = NULL;

  int          *requested_face_n   = (int *) malloc(n_rank * sizeof(int));
  int          *requested_face_idx = (int *) malloc((n_rank + 1) * sizeof(int));

  PDM_g_num_t *dface_group_proc = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *dface_group     = (PDM_g_num_t *) malloc(ppart->dn_face * sizeof(PDM_g_num_t));

  for (int igroup = 0; igroup < ppart->n_face_group; igroup++) {

    /*
     *  Build dface_group_proc
     */

    PDM_g_num_t n_face_group = ppart->_dface_group_idx[igroup+1]
                             - ppart->_dface_group_idx[igroup];

    PDM_MPI_Allgather(&n_face_group,
                      1,
                      PDM__PDM_MPI_G_NUM,
                      (void *) (&dface_group_proc[1]),
                      1,
                      PDM__PDM_MPI_G_NUM,
                      ppart->comm);

    dface_group_proc[0] = 0;
    PDM_array_accumulate_gnum(dface_group_proc, n_rank+1);

    /*
     *  Build dface_group
     */

    const int n_data_g = 2;

    PDM_array_reset_int(face_to_send_idx, n_rank+1, 0);
    PDM_array_reset_int(face_to_send_n, n_rank, 0);
    PDM_array_reset_gnum(dface_group, ppart->dn_face, -1);

    for (int i = ppart->_dface_group_idx[igroup];
             i < ppart->_dface_group_idx[igroup+1];
             i++) {
      PDM_g_num_t iface =  ppart->_dface_group[i];
      int found_rank = PDM_search_rank(iface,  ppart->dface_proc, 0, n_rank);
      face_to_send_idx[found_rank+1] += n_data_g;
    }

    for (int i = 0; i < n_rank; i++) {
      face_to_send_idx[i+1] += face_to_send_idx[i] ;
    }

    face_to_send = (PDM_g_num_t *) malloc(face_to_send_idx[n_rank] * sizeof(PDM_g_num_t));
    for (int i = ppart->_dface_group_idx[igroup];
         i < ppart->_dface_group_idx[igroup+1];
         i++) {
      PDM_g_num_t iface = ppart->_dface_group[i];
      int found_rank = PDM_search_rank(iface, ppart->dface_proc, 0, n_rank);
      int idx = face_to_send_idx[found_rank] + face_to_send_n[found_rank];
      face_to_send[idx++] = dface_group_proc[i_rank] + (PDM_g_num_t) (i - ppart->_dface_group_idx[igroup] + 1);
      face_to_send[idx++] = iface;
      face_to_send_n[found_rank] += n_data_g;
    }

    PDM_g_num_t *requested_face = NULL;

    _alltoall(face_to_send,
              face_to_send_n,
              face_to_send_idx,
              (void **) &requested_face,
              requested_face_n,
              requested_face_idx,
              PDM__PDM_MPI_G_NUM,
              sizeof(PDM_g_num_t),
              ppart->comm);

    free(face_to_send);
    face_to_send = NULL;

    int idx = 0;
    for (int i = 0; i < requested_face_idx[n_rank]/n_data_g; i++) {
      PDM_g_num_t iface      = requested_face[idx++];
      PDM_g_num_t gfaceGroup = requested_face[idx++];
      PDM_g_num_t _lface = gfaceGroup - ppart->dface_proc[i_rank];
      int          lface = (int) _lface;
      dface_group[lface] = (PDM_g_num_t) iface;
    }

    /*
     *  As distributes_faces
     */

    const int n_data     = 1;
    const int n_data_face = 1;

    for (int i_part = 0; i_part < ppart->mn_part; i_part++) {

      _part_t *mesh_part  = NULL;
      int *all_to_all_n_to_ln = NULL;

      PDM_array_reset_int(face_to_send_idx, n_rank+1, 0);
      PDM_array_reset_int(face_to_send_n, n_rank, 0);

      face_to_send = NULL;

      if (i_part < ppart->n_part) {

        mesh_part  = ppart->mesh_parts[i_part];
        all_to_all_n_to_ln = (int *) malloc(mesh_part->n_face * sizeof(int));

        /*
         *  Processes exchange list of faces which they want receive information
         */

        for (int i = 0; i < mesh_part->n_face; i++) {
          PDM_g_num_t iface = mesh_part->face_ln_to_gn[i];
          int found_rank = PDM_search_rank(iface, ppart->dface_proc, 0, n_rank);
          face_to_send_idx[found_rank+1] += n_data;
        }

        for (int i = 0; i < n_rank; i++) {
          face_to_send_idx[i+1] += face_to_send_idx[i] ;
        }

        face_to_send = (PDM_g_num_t *) malloc(face_to_send_idx[n_rank] * sizeof(PDM_g_num_t));

        for (int i = 0; i < mesh_part->n_face; i++) {

          PDM_g_num_t iface = mesh_part->face_ln_to_gn[i];
          int found_rank = PDM_search_rank(iface, ppart->dface_proc, 0, n_rank);

          idx = face_to_send_idx[found_rank] + face_to_send_n[found_rank];

          all_to_all_n_to_ln[idx/n_data] = i;

          face_to_send[idx++]   = iface;        /* Face global numbering */
          face_to_send_n[found_rank] += n_data;
        }
      }

      if (requested_face != NULL)
        free(requested_face);
      requested_face    = NULL;

      PDM_g_num_t *requested_face2 = NULL;

      _alltoall(face_to_send,
                face_to_send_n,
                face_to_send_idx,
                (void **) &requested_face2,
                requested_face_n,
                requested_face_idx,
                PDM__PDM_MPI_G_NUM,
                sizeof(PDM_g_num_t),
                ppart->comm);

      if (face_to_send != NULL)
        free(face_to_send);

      /*
       *  Processes exchange information about requested faces
       *  For each face, information contains :
       *     - Face global number in the current group
       *
       */

      int *sface_info_idx = face_to_send_idx;
      int *sface_info_n   = face_to_send_n;

      PDM_array_reset_int(sface_info_idx, n_rank+1, 0);
      PDM_array_reset_int(sface_info_n, n_rank, 0);

      for (int i = 0; i < n_rank; i++) {
        for (int k = requested_face_idx[i]; k < requested_face_idx[i+1]; k+=n_data) {
          sface_info_idx[i+1] += n_data_face;
        }
      }

      for (int i = 0; i < n_rank; i++) {
        sface_info_idx[i+1] += sface_info_idx[i] ;
      }

      PDM_g_num_t *sface_info = (PDM_g_num_t *)
        malloc(sface_info_idx[n_rank] * sizeof(PDM_g_num_t));

      for (int i = 0; i < n_rank; i++) {
        for (int k = requested_face_idx[i]; k < requested_face_idx[i+1]; k+=n_data) {
          PDM_g_num_t gface     = requested_face2[k];
          PDM_g_num_t _lface    = gface - ppart->dface_proc[i_rank];
          int          lface     = (int) _lface;

          idx = sface_info_idx[i] + sface_info_n[i];

          sface_info[idx++] = dface_group[lface];                   /* Number of vertices */
          sface_info_n[i] += 1;

        }
      }

      free(requested_face2);

      PDM_g_num_t *rface_info      = NULL;
      int          *rface_info_n   = requested_face_n;
      int          *rface_info_idx = requested_face_idx;

      _alltoall(sface_info,
                sface_info_n,
                sface_info_idx,
                (void **) &rface_info,
                rface_info_n,
                rface_info_idx,
                PDM__PDM_MPI_G_NUM,
                sizeof(PDM_g_num_t),
                ppart->comm);

      free(sface_info);
      sface_info = NULL;

      if (i_part < ppart->n_part) {

        /* Complete face_group_idx faceGroupeFace */

        if (igroup == 0) {
          mesh_part->face_group_idx = PDM_array_zeros_int(ppart->n_face_group+1);
        }

        mesh_part->face_group_idx[igroup+1] = mesh_part->face_group_idx[igroup];
        for (int i = 0; i < rface_info_idx[n_rank]; i++)
          if (rface_info[i] != -1)
            mesh_part->face_group_idx[igroup+1] += 1;

        if (igroup == 0) {
          mesh_part->face_group = (int *) malloc(mesh_part->face_group_idx[igroup+1] * sizeof(int));
          mesh_part->face_group_ln_to_gn =
            (PDM_g_num_t *) malloc(mesh_part->face_group_idx[igroup+1]
                                    * sizeof(PDM_g_num_t));
        }
        else {
          mesh_part->face_group =
            (int *) realloc(mesh_part->face_group, mesh_part->face_group_idx[igroup+1] * sizeof(int));
          mesh_part->face_group_ln_to_gn =
            (PDM_g_num_t *) realloc(mesh_part->face_group_ln_to_gn, mesh_part->face_group_idx[igroup+1]
                                     * sizeof(PDM_g_num_t));
        }

        idx = mesh_part->face_group_idx[igroup];
        for (int i = 0; i < rface_info_idx[n_rank]; i++) {
          if (rface_info[i] != -1) {
            mesh_part->face_group[idx] = all_to_all_n_to_ln[i]+1;
            mesh_part->face_group_ln_to_gn[idx] = rface_info[i];
            idx += 1;
          }

        }
      }

      if (rface_info != NULL)
        free(rface_info);

      if (all_to_all_n_to_ln != NULL)
        free(all_to_all_n_to_ln);

    }  /* For i_part */

  } /* For n_face_group */

  if (1 == 0) {
    for (int i_part = 0; i_part < ppart->n_part; i_part++) {

      _part_t *mesh_part  = ppart->mesh_parts[i_part];

      PDM_printf("mesh_part->n_face_group : %i\n",  ppart->n_face_group);
      PDM_printf("mesh_part->face_group : \n");
      for (int i1 = 0; i1 < ppart->n_face_group; i1++) {
        for (int i2 = mesh_part->face_group_idx[i1]; i2 < mesh_part->face_group_idx[i1+1]; i2++)
          PDM_printf(" %i", mesh_part->face_group[i2]);
        PDM_printf(" --\n");
      }
      PDM_printf("mesh_part->face_group_ln_to_gn : \n");
      for (int i1 = 0; i1 < ppart->n_face_group; i1++) {
        for (int i2 = mesh_part->face_group_idx[i1]; i2 < mesh_part->face_group_idx[i1+1]; i2++)
          PDM_printf(" "PDM_FMT_G_NUM, mesh_part->face_group_ln_to_gn[i2]);
        PDM_printf(" --\n");
      }
    }
  }

  free(face_to_send_n);
  free(face_to_send_idx);
  free(requested_face_n);
  free(requested_face_idx);
  free(dface_group_proc);
  free(dface_group);
}

/**
 *
 * \brief Free partition
 *
 * \param [in]   part      partition
 *
 */

static void
_part_free
(
 _part_t *part
)
{
  if (part->cell_face_idx != NULL)
    free(part->cell_face_idx);
  part->cell_face_idx = NULL;

  if (part->gcell_face != NULL)
    free(part->gcell_face);
  part->gcell_face = NULL;

  if (part->cell_face != NULL)
    free(part->cell_face);
  part->cell_face = NULL;

  if (part->cell_ln_to_gn != NULL)
    free(part->cell_ln_to_gn);
  part->cell_ln_to_gn = NULL;

  if (part->cell_tag != NULL)
    free(part->cell_tag);
  part->cell_tag = NULL;

  if (part->face_cell != NULL)
    free(part->face_cell);
  part->face_cell = NULL;

  if (part->face_vtx_idx != NULL)
    free(part->face_vtx_idx);
  part->face_vtx_idx = NULL;

  if (part->gface_vtx != NULL)
    free(part->gface_vtx);
  part->gface_vtx = NULL;

  if (part->face_vtx != NULL)
    free(part->face_vtx);
  part->face_vtx = NULL;

  if (part->face_ln_to_gn != NULL)
    free(part->face_ln_to_gn);
  part->face_ln_to_gn = NULL;

  if (part->face_tag != NULL)
    free(part->face_tag);
  part->face_tag = NULL;

  if (part->edge_ln_to_gn != NULL)
    free(part->edge_ln_to_gn);
  part->edge_ln_to_gn = NULL;

  if (part->edge_tag != NULL)
    free(part->edge_tag);
  part->edge_tag = NULL;


  if (part->edge_face_idx != NULL)
    free(part->edge_face_idx);
  part->edge_face_idx = NULL;

  if (part->edge_face != NULL)
    free(part->edge_face);
  part->edge_face = NULL;

  if (part->face_edge_idx != NULL)
    free(part->face_edge_idx);
  part->face_edge_idx = NULL;

  if (part->face_edge != NULL)
    free(part->face_edge);
  part->face_edge = NULL;

  if (part->edge_vtx != NULL)
    free(part->edge_vtx);
  part->edge_vtx = NULL;

  if (part->face_part_bound_proc_idx != NULL)
    free(part->face_part_bound_proc_idx);
  part->face_part_bound_proc_idx = NULL;

  if (part->face_part_bound_part_idx != NULL)
    free(part->face_part_bound_part_idx);
  part->face_part_bound_part_idx = NULL;

  if (part->face_part_bound != NULL)
    free(part->face_part_bound);
  part->face_part_bound = NULL;

  if (part->face_group_idx != NULL)
    free(part->face_group_idx);
  part->face_group_idx = NULL;

  if (part->face_group != NULL)
    free(part->face_group);
  part->face_group = NULL;

  if (part->face_group_ln_to_gn != NULL)
    free(part->face_group_ln_to_gn);
  part->face_group_ln_to_gn = NULL;

  if (part->vtx != NULL)
    free(part->vtx);
  part->vtx = NULL;

  if (part->vtx_ln_to_gn != NULL)
    free(part->vtx_ln_to_gn);
  part->vtx_ln_to_gn = NULL;

  if (part->vtx_tag != NULL)
    free(part->vtx_tag);
  part->vtx_tag = NULL;

  if (part->cell_color != NULL)
    free(part->cell_color);
  part->cell_color = NULL;

  if (part->face_color != NULL)
    free(part->face_color);
  part->face_color = NULL;

  if (part->edge_color != NULL)
    free(part->edge_color);
  part->edge_color = NULL;

  if (part->vtx_color != NULL)
    free(part->vtx_color);
  part->vtx_color = NULL;

  if (part->thread_color != NULL)
    free(part->thread_color);
  part->thread_color = NULL;

  if (part->hyperplane_color != NULL)
    free(part->hyperplane_color);
  part->hyperplane_color = NULL;

  if (part->new_to_old_order_cell != NULL)
    free(part->new_to_old_order_cell);
  part->new_to_old_order_cell = NULL;

  if (part->new_to_old_order_face != NULL)
    free(part->new_to_old_order_face);
  part->new_to_old_order_face = NULL;

  if (part->new_to_old_order_edge != NULL)
    free(part->new_to_old_order_edge);
  part->new_to_old_order_edge = NULL;

  if (part->new_to_old_order_vtx != NULL)
    free(part->new_to_old_order_vtx);
  part->new_to_old_order_vtx = NULL;

  if(part->subpartlayout != NULL){
    if(part->subpartlayout->cell_tile_idx!= NULL)
      free(part->subpartlayout->cell_tile_idx);
    if(part->subpartlayout->face_tile_idx!= NULL)
      free(part->subpartlayout->face_tile_idx);
    if(part->subpartlayout->face_bnd_tile_idx!= NULL)
      free(part->subpartlayout->face_bnd_tile_idx);
    if(part->subpartlayout->mask_tile_idx!= NULL)
      free(part->subpartlayout->mask_tile_idx);
    if(part->subpartlayout->cell_vect_tile_idx!= NULL)
      free(part->subpartlayout->cell_vect_tile_idx);
    if(part->subpartlayout->mask_tile_n!= NULL)
      free(part->subpartlayout->mask_tile_n);
    if(part->subpartlayout->cell_vect_tile_n!= NULL)
      free(part->subpartlayout->cell_vect_tile_n);
    if(part->subpartlayout->mask_tile!= NULL)
      free(part->subpartlayout->mask_tile);
    free(part->subpartlayout);
  }

  free(part->n_elt);
  part->n_elt = NULL;
  for (int i=0; i<part->n_section; ++i) {
    free(part->elt_vtx_idx[i]);
    free(part->elt_vtx[i]);
    free(part->elt_section_ln_to_gn[i]);
  }
  free(part->elt_vtx_idx);
  free(part->elt_vtx);
  free(part->elt_section_ln_to_gn);
  part->elt_vtx_idx          = NULL;
  part->elt_vtx              = NULL;
  part->elt_section_ln_to_gn = NULL;
  part->n_section            = 0;

  free(part);
}

/**
 *
 * \brief Free partition
 *
 * \param [in]   part      partition
 *
 */

static void
_part_partial_free
(
 _part_t *part
)
{
  if (part->cell_face_idx != NULL)
    free(part->cell_face_idx);
  part->cell_face_idx = NULL;

  if (part->gcell_face != NULL)
    free(part->gcell_face);
  part->gcell_face = NULL;

  if (part->cell_face != NULL)
    free(part->cell_face);
  part->cell_face = NULL;

  if (part->face_cell != NULL)
    free(part->face_cell);
  part->face_cell = NULL;

  if (part->face_vtx_idx != NULL)
    free(part->face_vtx_idx);
  part->face_vtx_idx = NULL;

  if (part->gface_vtx != NULL)
    free(part->gface_vtx);
  part->gface_vtx = NULL;

  if (part->face_vtx != NULL)
    free(part->face_vtx);
  part->face_vtx = NULL;

  if (part->vtx != NULL)
    free(part->vtx);
  part->vtx = NULL;

  if (part->new_to_old_order_cell != NULL)
    free(part->new_to_old_order_cell);
  part->new_to_old_order_cell = NULL;

  if (part->new_to_old_order_face != NULL)
    free(part->new_to_old_order_face);
  part->new_to_old_order_face = NULL;

  if (part->new_to_old_order_vtx != NULL)
    free(part->new_to_old_order_vtx);
  part->new_to_old_order_vtx = NULL;
}

/* Wrappers of pdm_multipart for CEDRE */

static PDM_multipart_t *
_multipart_create
(
 const PDM_MPI_Comm           comm,
 const PDM_part_split_t       method,
 const char                  *renum_cell_method,
 const char                  *renum_face_method,
 const int                    n_property_cell,
 const int                   *renum_properties_cell,
 const int                    n_property_face,
 const int                   *renum_properties_face,
 const int                    n_part,
 const int                    dn_cell,
 const int                    dn_face,
 const int                    dn_vtx,
 const int                    n_face_group,
 const int                   *dcell_face_idx,
 const PDM_g_num_t           *dcell_face,
 const int                   *dcell_tag,
 const int                   *dcell_weight,
 const int                    have_dcell_part,
       int                   *dcell_part,
 const PDM_g_num_t           *dface_cell,
 const int                   *dface_vtx_idx,
 const PDM_g_num_t           *dface_vtx,
 const int                   *dface_tag,
 const double                *dvtx_coord,
 const int                   *dvtx_tag,
 const int                   *dface_group_idx,
 const PDM_g_num_t           *dface_group
 )
 {

  // unused for multipart
  PDM_UNUSED(n_property_cell      );
  PDM_UNUSED(n_property_face      );
  PDM_UNUSED(renum_properties_face);
  PDM_UNUSED(dcell_face_idx       );
  PDM_UNUSED(dcell_face           );
  PDM_UNUSED(dcell_tag            );
  PDM_UNUSED(dcell_weight         );
  PDM_UNUSED(have_dcell_part      );
  PDM_UNUSED(dcell_part           );
  PDM_UNUSED(dface_tag            );
  PDM_UNUSED(dvtx_tag             );

  PDM_split_dual_t split_method = (PDM_split_dual_t) method; // TO DO: uniformise method type

  // Partitioning
  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *multipart = PDM_multipart_create(n_domain,
                                                    &n_part_domains,
                                                    PDM_FALSE,
                                                    split_method,
                                                    PDM_PART_SIZE_HOMOGENEOUS,
                                                    NULL,
                                                    comm,
                                                    PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(multipart, -1, renum_cell_method, renum_properties_cell, renum_face_method);

  PDM_multipart_block_set(multipart,
                          0,
                          dn_cell,
                          dn_face,
                          dn_vtx,
                          n_face_group,
          (int         *) dcell_face_idx,
          (PDM_g_num_t *) dcell_face,
          (PDM_g_num_t *) dface_cell,
          (int         *) dface_vtx_idx,
          (PDM_g_num_t *) dface_vtx,
          (double      *) dvtx_coord,
          (int         *) dface_group_idx,
          (PDM_g_num_t *) dface_group);

  // Run
  PDM_multipart_compute(multipart);

  return multipart;
 }

static void
_dim_get
(
 PDM_multipart_t *multipart,
 const  int       i_part,
        int      *n_cell,
        int      *n_face,
        int      *n_face_part_bound,
        int      *n_vtx,
        int      *n_proc,
        int      *n_total_part,
        int      *scell_face,
        int      *sface_vtx,
        int      *sface_group,
        int      *n_face_group
)
{
  int  i_domain = 0;

  PDM_multipart_part_dim_get(multipart,
                             i_domain,
                             i_part,
                             n_cell,
                             n_face,
                             n_face_part_bound,
                             n_vtx,
                             n_proc,
                             n_total_part,
                             scell_face,
                             sface_vtx,
                             sface_group,
                             n_face_group);
}

static void
_val_get
(
PDM_multipart_t   *multipart,
const  int         i_part,
 int             **cell_tag,
 int             **cell_face_idx,
 int             **cell_face,
 PDM_g_num_t     **cell_ln_to_gn,
 int             **face_tag,
 int             **face_cell,
 int             **face_vtx_idx,
 int             **face_vtx,
 PDM_g_num_t     **face_ln_to_gn,
 int             **face_part_bound_proc_idx,
 int             **face_part_bound_part_idx,
 int             **face_part_bound,
 int             **vtx_tag,
 double          **vtx,
 PDM_g_num_t     **vtx_ln_to_gn,
 int             **face_group_idx,
 int             **face_group,
 PDM_g_num_t     **face_group_ln_to_gn
)
{
  int i_domain = 0;

  *cell_tag = NULL;
  *face_tag = NULL;
  *vtx_tag  = NULL;

  PDM_multipart_part_ln_to_gn_get(multipart,
                                  i_domain,
                                  i_part,
                                  PDM_MESH_ENTITY_CELL,
                                  cell_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_connectivity_get(multipart,
                                      i_domain,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                      cell_face_idx,
                                      cell_face,
                                      PDM_OWNERSHIP_KEEP);

  int *face_cell_idx = NULL;
  PDM_multipart_part_connectivity_get(multipart,
                                      i_domain,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                      &face_cell_idx,
                                      face_cell,
                                      PDM_OWNERSHIP_KEEP);
  assert(face_cell_idx == NULL);

  PDM_multipart_part_connectivity_get(multipart,
                                      i_domain,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                      face_vtx_idx,
                                      face_vtx,
                                      PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_ln_to_gn_get(multipart,
                                  i_domain,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  face_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_ln_to_gn_get(multipart,
                                  i_domain,
                                  i_part,
                                  PDM_MESH_ENTITY_VTX,
                                  vtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
  int pn_face_group = 0;
  PDM_multipart_group_get(multipart,
                          i_domain,
                          i_part,
                          PDM_MESH_ENTITY_FACE,
                          &pn_face_group,
                          face_group_idx,
                          face_group,
                          face_group_ln_to_gn,
                          PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_graph_comm_get(multipart,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    face_part_bound_proc_idx,
                                    face_part_bound_part_idx,
                                    face_part_bound,
                                    PDM_OWNERSHIP_KEEP);
  PDM_multipart_part_vtx_coord_get(multipart,
                                   i_domain,
                                   i_part,
                                   vtx,
                                   PDM_OWNERSHIP_KEEP);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a initial partitioning
 *  Build a initial partitioning from :
 *      - Cell block distribution with implicit global numbering
 *         (the first cell is the first cell of the first process and
 *          the latest cell is the latest cell of the latest process)
 *      - Face block distribution with implicit global numbering
 *      - Vertex block distribution with implicit global numbering
 *  To repart an existing partition use \ref PDM_part_repart function
 *
 * \param [in]   comm                   MPI Comminicator
 * \param [in]   split_method           Split method
 * \param [in]   renum_cell_method      Cell renumbering method
 * \param [in]   renum_face_method      Face renumbering method
 * \param [in]   n_property_cell        Number of cell properties
 * \param [in]   renum_properties_cell  For cache blocking [ n_cell_per_cache_wanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
 * \param [in]   n_property_face        Number of face properties
 * \param [in]   renum_properties_face  NOT USED
 * \param [in]   n_part                 Number of partition to build on this process
 * \param [in]   dn_cell                Number of distributed cells
 * \param [in]   dn_face                Number of distributed faces
 * \param [in]   dn_vtx                 Number of distributed vertices
 * \param [in]   n_face_group           Number of face groups
 * \param [in]   dcell_face_idx         Distributed cell face connectivity index or NULL
 *                                      (size : dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face             Distributed cell face connectivity or NULL
 *                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 * \param [in]   dcell_tag              Cell tag (size : n_cell) or NULL
 * \param [in]   dcell_weight           Cell weight (size : n_cell) or NULL
 * \param [in]   have_dcell_part        Presence of an array of cell part id
 * \param [in]   dcell_part             Distributed cell partitioning
 *                                      (size = dn_cell) or NULL (No partitioning if != NULL)
 * \param [in]   dface_cell             Distributed face cell connectivity or NULL
 *                                      (size : 2 * dn_face, numbering : 1 to n)
 * \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
 *                                      (size : dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx              Distributed face to vertex connectivity
 *                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 * \param [in]   dface_tag              Distributed face tag (size : dn_face)
 *                                      or NULL
 * \param [in]   dvtx_coord             Distributed vertex coordinates
 *                                      (size : 3*dn_vtx)
 * \param [in]   dvtx_tag               Distributed vertex tag (size : dn_vtx) or NULL
 * \param [in]   dface_group_idx        Index of distributed faces list of each group
 *                                      (size = n_face_group + 1) or NULL
 * \param [in]   dface_group            Distributed faces list of each group
 *                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
 *                                      or NULL
 *
 * \return    Pointer to \ref PDM_part object
 *
 */

PDM_part_t *
PDM_part_create
(
 const PDM_MPI_Comm           comm,
 const PDM_part_split_t       split_method,
 const char                  *renum_cell_method,
 const char                  *renum_face_method,
 const int                    n_property_cell,
 const int                   *renum_properties_cell,
 const int                    n_property_face,
 const int                   *renum_properties_face,
 const int                    n_part,
 const int                    dn_cell,
 const int                    dn_face,
 const int                    dn_vtx,
 const int                    n_face_group,
 const int                   *dcell_face_idx,
 const PDM_g_num_t           *dcell_face,
 const int                   *dcell_tag,
 const int                   *dcell_weight,
 const int                    have_dcell_part,
       int                   *dcell_part,
 const PDM_g_num_t           *dface_cell,
 const int                   *dface_vtx_idx,
 const PDM_g_num_t           *dface_vtx,
 const int                   *dface_tag,
 const double                *dvtx_coord,
 const int                   *dvtx_tag,
 const int                   *dface_group_idx,
 const PDM_g_num_t           *dface_group
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_part_renum_method_load_local();

  _PDM_part_t *_ppart = malloc(sizeof(_PDM_part_t));

  char *use_multipart_var;
  _ppart->use_multipart = 0;
  if (( use_multipart_var = getenv( "PDM_USE_MULTIPART" )) != NULL ) {
    _ppart->use_multipart = atoi(use_multipart_var);
  }

  if (_ppart->use_multipart) {
    _ppart->multipart = _multipart_create(comm,
                                          split_method,
                                          renum_cell_method,
                                          renum_face_method,
                                          n_property_cell,
                                          renum_properties_cell,
                                          n_property_face,
                                          renum_properties_face,
                                          n_part,
                                          dn_cell,
                                          dn_face,
                                          dn_vtx,
                                          n_face_group,
                                          dcell_face_idx,
                                          dcell_face,
                                          dcell_tag,
                                          dcell_weight,
                                          have_dcell_part,
                                          dcell_part,
                                          dface_cell,
                                          dface_vtx_idx,
                                          dface_vtx,
                                          dface_tag,
                                          dvtx_coord,
                                          dvtx_tag,
                                          dface_group_idx,
                                          dface_group);
  } else {

    /*
     * Build ppart structure
     */

    _ppart->timer = PDM_timer_create();
    for (int i = 0; i < 4; i++) {
      _ppart->times_elapsed[i] = 0.;
      _ppart->times_cpu[i] = 0.;
      _ppart->times_cpu_u[i] = 0.;
      _ppart->times_cpu_s[i] = 0.;
    }
    PDM_timer_resume(_ppart->timer);

    /* Local dimensions */

    _ppart->dn_vtx       = dn_vtx;
    _ppart->dn_cell      = dn_cell;
    _ppart->dn_face      = dn_face;
    _ppart->n_face_group = n_face_group;

    /* Cell definitions */

    _ppart->_dcell_face_idx = dcell_face_idx;
    _ppart->_dcell_face     = dcell_face;
    _ppart->_dcell_tag      = dcell_tag;
    _ppart->_dcell_weight   = dcell_weight;
    _ppart->_dcell_part     = dcell_part;
    _ppart->dcell_face_idx  = NULL;
    _ppart->dcell_face      = NULL;
    _ppart->dface_cell      = NULL;

    /* Set up for renumbering */
    _ppart->n_property_cell       = n_property_cell;
    _ppart->renum_properties_cell = renum_properties_cell;
    _ppart->n_property_face       = n_property_face;
    _ppart->renum_properties_face = renum_properties_face;

    _ppart->dcell_proc = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
    PDM_g_num_t _dn_cell = (PDM_g_num_t) dn_cell;
    PDM_MPI_Allgather((void *) &_dn_cell,
                      1,
                      PDM__PDM_MPI_G_NUM,
                      (void *) (&_ppart->dcell_proc[1]),
                      1,
                      PDM__PDM_MPI_G_NUM,
                      comm);

    _ppart->dcell_proc[0] = 1;

    for (int i = 1; i < n_rank+1; i++) {
      _ppart->dcell_proc[i] +=  _ppart->dcell_proc[i-1];
    }

    if (1 == 0) {
      PDM_printf("_ppart->dcell_proc : "PDM_FMT_G_NUM,  _ppart->dcell_proc[0]);
      for (int i = 1; i < n_rank+1; i++) {
        PDM_printf(" "PDM_FMT_G_NUM, _ppart->dcell_proc[i]);
      }
      PDM_printf("\n");
    }

    /* Face definitions */

    _ppart->_dface_tag     = dface_tag;
    _ppart->_dface_cell    = dface_cell;
    _ppart->_dface_vtx_idx = dface_vtx_idx;
    _ppart->_dface_vtx     = dface_vtx;

    _ppart->dface_proc = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
    int *dn_face_proc = (int *) malloc((n_rank) * sizeof(int));

    PDM_MPI_Allgather((void *) &dn_face,
                      1,
                      PDM_MPI_INT,
                      (void *) dn_face_proc,
                      1,
                      PDM_MPI_INT,
                      comm);
    _ppart->dface_proc[0] = 1;
    for (int i = 1; i < n_rank+1; i++) {
      _ppart->dface_proc[i] = (PDM_g_num_t) dn_face_proc[i-1] + _ppart->dface_proc[i-1];
    }

    free(dn_face_proc);

    /* Vertex definitions */

    _ppart->_dvtx_coord    = dvtx_coord;
    _ppart->_dvtx_tag      = dvtx_tag;

    _ppart->dvtx_proc = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
    int *dn_vtx_proc = (int *) malloc((n_rank) * sizeof(int));

    PDM_MPI_Allgather((void *) &dn_vtx,
                      1,
                      PDM_MPI_INT,
                      (void *) dn_vtx_proc,
                      1,
                      PDM_MPI_INT,
                      comm);
    _ppart->dvtx_proc[0] = 1;
    for (int i = 1; i < n_rank+1; i++) {
      _ppart->dvtx_proc[i] = dn_vtx_proc[i-1] + _ppart->dvtx_proc[i-1];
    }

    free(dn_vtx_proc);

    /* Boundaries definitions */

    _ppart->_dface_group_idx = dface_group_idx;
    _ppart->_dface_group     = dface_group;

    /* Dual graph */

    _ppart->ddual_graph_idx = NULL;
    _ppart->ddual_graph     = NULL;

    /* Partitions */

    _ppart->n_part = n_part;

    _ppart->mn_part = -1;

    _ppart->dpart_proc = (int *) malloc((n_rank + 1) * sizeof(int));
    PDM_MPI_Allgather((void *) &n_part,
                      1,
                      PDM_MPI_INT,
                      (void *) (&_ppart->dpart_proc[1]),
                      1,
                      PDM_MPI_INT,
                      comm);

    _ppart->dpart_proc[0] = 0;
    for (int i = 1; i < n_rank+1; i++) {
      _ppart->mn_part = PDM_MAX(_ppart->mn_part, _ppart->dpart_proc[i]);
      _ppart->dpart_proc[i] = _ppart->dpart_proc[i] + _ppart->dpart_proc[i-1];
    }

    _ppart->tn_part =  _ppart->dpart_proc[n_rank];

    _ppart->gpart_to_lproc_part = (int *) malloc(2*_ppart->tn_part * sizeof(int));

    for (int i = 0; i < n_rank; i++) {
      for (int j = _ppart->dpart_proc[i]; j < _ppart->dpart_proc[i+1]; j++) {
        _ppart->gpart_to_lproc_part[2*j    ] = i;
        _ppart->gpart_to_lproc_part[2*j + 1] = j - _ppart->dpart_proc[i];
      }
    }

    _ppart->mesh_parts = (_part_t **) malloc(_ppart->n_part * sizeof(_part_t *));

    for (int i = 0; i < _ppart->n_part; i++)
      _ppart->mesh_parts[i] = NULL;

    /* Communicator */

    _ppart->comm = comm;

    /* Method */

    _ppart->split_method = split_method;

    int _method = PDM_part_renum_method_cell_idx_get(renum_cell_method);

    if (_method == -1) {
      PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering cell method\n", renum_cell_method);
    }

    _ppart->renum_cell_method = _method;

    _method = PDM_part_renum_method_face_idx_get(renum_face_method);

    if (_method == -1) {
      PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering face method\n", renum_face_method);
    }
    _ppart->renum_face_method = _method;

    _ppart->dpart_bound = NULL;

    /*
     * Build dual graph
     */

    if (dcell_face != NULL)
      _dual_graph_from_cell_face(_ppart);
    else if (dface_cell != NULL)
      _dual_graph_from_face_cell(_ppart);
    else {
      PDM_printf("PDM_part_part_create error : dcell_face and dface_cell are undefined, define one of two\n");
      exit(1);
    }

    int itime = 1;
    PDM_timer_hang_on(_ppart->timer);
    _ppart->times_elapsed[itime] = PDM_timer_elapsed(_ppart->timer);
    _ppart->times_cpu[itime]     = PDM_timer_cpu(_ppart->timer);
    _ppart->times_cpu_u[itime]   = PDM_timer_cpu_user(_ppart->timer);
    _ppart->times_cpu_s[itime]   = PDM_timer_cpu_sys(_ppart->timer);
    itime += 1;

    /*
     * Graph partitioning
     */

    PDM_timer_resume(_ppart->timer);

    int *cell_part;

    if (have_dcell_part == 0) {
      cell_part = (int *) malloc(dn_cell * sizeof(int));
      _split(_ppart,
             cell_part);
      for (int i = 0; i < dn_cell; i++) {
        dcell_part[i] = cell_part[i];
      }
    }

    else {
      cell_part = (int *) _ppart->_dcell_part;
    }

    if (1 == 0) {
      PDM_printf("cell_part : ");
      for (int i = 0; i <dn_cell; i++)
        PDM_printf(" %d", cell_part[i]);
      PDM_printf("\n");
    }

    PDM_timer_hang_on(_ppart->timer);
    _ppart->times_elapsed[itime] = PDM_timer_elapsed(_ppart->timer);
    _ppart->times_cpu[itime]     = PDM_timer_cpu(_ppart->timer);
    _ppart->times_cpu_u[itime]   = PDM_timer_cpu_user(_ppart->timer);
    _ppart->times_cpu_s[itime]   = PDM_timer_cpu_sys(_ppart->timer);
    itime += 1;

    /*
     * Cell distribution to build local connectivities
     *     - cell_face_idx : ok
     *     - gcell_face   : ok
     *     - cell_face    : ok
     *     - cell_ln_to_gn  : ok
     *     - cell_tag     : ok
     *     - face_ln_to_gn  : ok
     */

    PDM_timer_resume(_ppart->timer);

    _distrib_cell(_ppart,
                  cell_part);

    if (have_dcell_part == 0) {
      free(cell_part);
      cell_part = NULL;
    }

    /*
     * Face distribution to build local connectivities
     *     - face_vtx_idx  : ok
     *     - gface_vtx    : ok
     *     - face_vtx     : ok
     *     - face_tag     : ok
     */

    _distrib_face(_ppart);
    _build_faceCell(_ppart);

    /*
     * Vertex distribution to build local connectivities
     */

    _distrib_vtx(_ppart);

    /*
     * Cell renumbering
     */

    PDM_part_renum_cell (        _ppart->mesh_parts,
                                 _ppart->n_part,
                                 _ppart->renum_cell_method,
                         (void*) _ppart->renum_properties_cell);

    /*
     * Face renumbering
     */
    PDM_part_renum_face (         _ppart->mesh_parts,
                                  _ppart->n_part,
                                  _ppart->renum_face_method,
                          (void*) _ppart->renum_properties_face);

    /*
     * Look for partitioning boundary faces
     */

    _search_part_bound_face(_ppart);

    /*
     * Face group distribution to build local connectivities
     */

    _distrib_face_groups(_ppart);

    PDM_timer_hang_on(_ppart->timer);
    _ppart->times_elapsed[itime] = PDM_timer_elapsed(_ppart->timer);
    _ppart->times_cpu[itime]     = PDM_timer_cpu(_ppart->timer);
    _ppart->times_cpu_u[itime]   = PDM_timer_cpu_user(_ppart->timer);
    _ppart->times_cpu_s[itime]   = PDM_timer_cpu_sys(_ppart->timer);

    _ppart->times_elapsed[0]     = _ppart->times_elapsed[itime];
    _ppart->times_cpu[0]         = _ppart->times_cpu[itime];
    _ppart->times_cpu_u[0]       = _ppart->times_cpu_u[itime];
    _ppart->times_cpu_s[0]       = _ppart->times_cpu_s[itime];

    for (int i = itime; i > 1; i--) {
      _ppart->times_elapsed[i] -= _ppart->times_elapsed[i-1];
      _ppart->times_cpu[i]     -= _ppart->times_cpu[i-1];
      _ppart->times_cpu_u[i]   -= _ppart->times_cpu_u[i-1];
      _ppart->times_cpu_s[i]   -= _ppart->times_cpu_s[i-1];
    }

    if (_ppart->dcell_face_idx != NULL)
      free(_ppart->dcell_face_idx);
    _ppart->dcell_face_idx = NULL;

    if (_ppart->dcell_face != NULL)
      free(_ppart->dcell_face);
    _ppart->dcell_face = NULL;

    if (_ppart->dcell_proc != NULL)
      free(_ppart->dcell_proc);
    _ppart->dcell_proc = NULL;

    if (_ppart->dface_proc != NULL)
      free(_ppart->dface_proc);
    _ppart->dface_proc = NULL;

    if (_ppart->dface_cell != NULL)
      free(_ppart->dface_cell);
    _ppart->dface_cell = NULL;

    if (_ppart->dvtx_proc != NULL)
      free(_ppart->dvtx_proc);
    _ppart->dvtx_proc = NULL;

    // if (_ppart->dpart_proc != NULL)
    //   free(_ppart->dpart_proc);
    // _ppart->dpart_proc = NULL;

    if (_ppart->gpart_to_lproc_part != NULL)
      free(_ppart->gpart_to_lproc_part);
    _ppart->gpart_to_lproc_part = NULL;

    if (_ppart->dpart_bound != NULL)
      free(_ppart->dpart_bound);
    _ppart->dpart_bound = NULL;

    if (_ppart->ddual_graph_idx != NULL)
      free(_ppart->ddual_graph_idx);
    _ppart->ddual_graph_idx = NULL;

    if (_ppart->ddual_graph != NULL)
      free(_ppart->ddual_graph);
    _ppart->ddual_graph = NULL;

  }

  return (PDM_part_t *) _ppart;
}


/**
 *
 * \brief Return a mesh partition dimensions
 *
 * \param [in]   ppart               Pointer to \ref PDM_part object
 * \param [in]   i_part              Current partition
 * \param [out]  n_cell              Number of cells
 * \param [out]  n_face              Number of faces
 * \param [out]  n_face_part_bound   Number of partitioning boundary faces
 * \param [out]  n_vtx               Number of vertices
 * \param [out]  n_proc              Number of processus
 * \param [out]  n_total_part        Number of partitions
 * \param [out]  scell_face          Size of cell-face connectivity
 * \param [out]  sface_vtx           Size of face-vertex connectivity
 * \param [out]  sFacePartBound      Size of face_part_bound array
 * \param [out]  sface_group         Size of face_group array
 * \param [out]  n_face_group        Number of face groups
 *
 */

void
PDM_part_part_dim_get
(
 PDM_part_t *ppart,
const   int  i_part,
        int *n_cell,
        int *n_face,
        int *n_face_part_bound,
        int *n_vtx,
        int *n_proc,
        int *n_total_part,
        int *scell_face,
        int *sface_vtx,
        int *sface_group,
        int *n_face_group
)
{

  if (ppart->use_multipart) {
    _dim_get(ppart->multipart,
             i_part,
             n_cell,
             n_face,
             n_face_part_bound,
             n_vtx,
             n_proc,
             n_total_part,
             scell_face,
             sface_vtx,
             sface_group,
             n_face_group);
  } else {
    _PDM_part_t *_ppart = (_PDM_part_t *) ppart;
    int n_rank;
    PDM_MPI_Comm_size(_ppart->comm, &n_rank);

    _part_t *mesh_part = NULL;
    if (i_part < _ppart->n_part)
      mesh_part  = _ppart->mesh_parts[i_part];

    if (mesh_part == NULL) {
      PDM_printf("PDM_part_part_get error : unknown partition\n");
      exit(1);
    }

    *n_cell            = mesh_part->n_cell;
    *n_face            = mesh_part->n_face;
    *n_face_part_bound = mesh_part->n_face_part_bound;
    *n_proc            = n_rank;
    *n_total_part      = _ppart->tn_part;
    *n_vtx             = mesh_part->n_vtx;
    *scell_face        = mesh_part->cell_face_idx[*n_cell];
    *sface_vtx         = mesh_part->face_vtx_idx[*n_face];
    *sface_group       = 0;
    if (_ppart->n_face_group > 0)
      *sface_group    = mesh_part->face_group_idx[_ppart->n_face_group];
    *n_face_group    = _ppart->n_face_group;
  }
}

/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   ppart                     Pointer to \ref PDM_part object
 * \param [in]   i_part                    Current partition
 * \param [out]  cell_tag                  Cell tag (size = n_cell)
 * \param [out]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 * \param [out]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 *                                                                   numbering : 1 to n)
 * \param [out]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 * \param [out]  face_tag                  Face tag (size = n_face)
 * \param [out]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 * \param [out]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 * \param [out]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 * \param [out]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 * \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
 * \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
 * \param [out]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
 *                                          sorted by processus, sorted by partition in each processus, and
 *                                          sorted by absolute face number in each partition
 *                                         For each face :
 *                                           - Face local number (numbering : 1 to n)
 *                                           - Connected process (numbering : 0 to n-1)
 *                                           - Connected Partition
 *                                             on the connected process (numbering :1 to n)
 *                                           - Connected face local number
 *                                             in the connected partition (numbering :1 to n)
 * \param [out]  vtx_tag                   Vertex tag (size = nVertex)
 * \param [out]  vtx                       Vertex coordinates (size = 3 * nVertex)
 * \param [out]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 * \param [out]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 * \param [out]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 * \param [out]  face_group_ln_to_gn       Faces global numbering for each group
 *                                         (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 *
 */

void PDM_part_part_val_get
(
PDM_part_t     *ppart,
const  int      i_part,
 int          **cell_tag,
 int          **cell_face_idx,
 int          **cell_face,
 PDM_g_num_t  **cell_ln_to_gn,
 int          **face_tag,
 int          **face_cell,
 int          **face_vtx_idx,
 int          **face_vtx,
 PDM_g_num_t  **face_ln_to_gn,
 int          **face_part_bound_proc_idx,
 int          **face_part_bound_part_idx,
 int          **face_part_bound,
 int          **vtx_tag,
 double       **vtx,
 PDM_g_num_t  **vtx_ln_to_gn,
 int          **face_group_idx,
 int          **face_group,
 PDM_g_num_t  **face_group_ln_to_gn
)
{
  if (ppart->use_multipart) {
    _val_get(ppart->multipart,
             i_part,
             cell_tag,
             cell_face_idx,
             cell_face,
             cell_ln_to_gn,
             face_tag,
             face_cell,
             face_vtx_idx,
             face_vtx,
             face_ln_to_gn,
             face_part_bound_proc_idx,
             face_part_bound_part_idx,
             face_part_bound,
             vtx_tag,
             vtx,
             vtx_ln_to_gn,
             face_group_idx,
             face_group,
             face_group_ln_to_gn);
  } else {
    _PDM_part_t *_ppart = (_PDM_part_t *) ppart;

    _part_t *mesh_part = NULL;
    if (i_part < _ppart->n_part)
      mesh_part  = _ppart->mesh_parts[i_part];

    if (mesh_part == NULL) {
      PDM_printf("PDM_part_part_val_get error : unknown partition\n");
      exit(1);
    }

    *cell_tag                 = mesh_part->cell_tag;
    *cell_face_idx            = mesh_part->cell_face_idx;
    *cell_face                = mesh_part->cell_face;
    *cell_ln_to_gn            = mesh_part->cell_ln_to_gn;
    *face_tag                 = mesh_part->face_tag;
    *face_cell                = mesh_part->face_cell;
    *face_vtx_idx             = mesh_part->face_vtx_idx;
    *face_vtx                 = mesh_part->face_vtx;
    *face_ln_to_gn            = mesh_part->face_ln_to_gn;
    *face_part_bound_proc_idx = mesh_part->face_part_bound_proc_idx;
    *face_part_bound_part_idx = mesh_part->face_part_bound_part_idx;
    *face_part_bound          = mesh_part->face_part_bound;
    *vtx_tag                  = mesh_part->vtx_tag;
    *vtx                      = mesh_part->vtx;
    *vtx_ln_to_gn             = mesh_part->vtx_ln_to_gn;
    *face_group_idx           = mesh_part->face_group_idx;
    *face_group               = mesh_part->face_group;
    *face_group_ln_to_gn      = mesh_part->face_group_ln_to_gn;
  }
}


/**
 *
 * \brief Return the coloring of a mesh partition
 *
 * \param [in]   ppart               Pointer to \ref PDM_part object
 * \param [in]   i_part              Current partition
 * \param [out]  cell_color          Cell color (size = n_cell)
 * \param [out]  face_color          Face color (size = n_face)
 * \param [out]  thread_color        Thread color (size = n_cell)
 * \param [out]  hyperplane_color    Hyperplane color (size = n_cell)
 *
 */

void PDM_part_part_color_get
(
PDM_part_t     *ppart,
const  int      i_part,
 int          **cell_color,
 int          **face_color,
 int          **thread_color,
 int          **hyperplane_color
)
{
  _PDM_part_t *_ppart = (_PDM_part_t *) ppart;

  if (_ppart->use_multipart) {
    int i_domain = 0;
    PDM_multipart_partition_color_get(ppart->multipart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      cell_color,
                                      PDM_OWNERSHIP_KEEP);
    PDM_multipart_partition_color_get(ppart->multipart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      face_color,
                                      PDM_OWNERSHIP_KEEP);
    PDM_multipart_part_hyperplane_color_get(ppart->multipart,
                                            i_domain,
                                            i_part,
                                            hyperplane_color,
                                            PDM_OWNERSHIP_KEEP);
    PDM_multipart_part_thread_color_get(ppart->multipart,
                                        i_domain,
                                        i_part,
                                        thread_color,
                                        PDM_OWNERSHIP_KEEP);
  } else {
    _part_t *mesh_part = NULL;
    if (i_part < _ppart->n_part)
      mesh_part  = _ppart->mesh_parts[i_part];

    if (mesh_part == NULL) {
      PDM_printf("PDM_part_part_color_get error : unknown partition\n");
      exit(1);
    }

    *cell_color       = mesh_part->cell_color;
    *face_color       = mesh_part->face_color;
    *thread_color     = mesh_part->thread_color;
    *hyperplane_color = mesh_part->hyperplane_color;
  }
}



/**
 *
 * \brief Free ppart
 *
 * \param [in]   ppart               Pointer to \ref PDM_part object
 *
 */

void
PDM_part_free
(
 PDM_part_t *ppart
)
{
  _PDM_part_t *_ppart = (_PDM_part_t *) ppart;

  if (_ppart->use_multipart) {

    if (ppart->multipart != NULL)
      PDM_multipart_free(ppart->multipart);
    _ppart->multipart = NULL;

  } else {

    if (_ppart->dcell_face_idx != NULL)
      free(_ppart->dcell_face_idx);
    _ppart->dcell_face_idx = NULL;

    if (_ppart->dcell_face != NULL)
      free(_ppart->dcell_face);
    _ppart->dcell_face = NULL;

    if (_ppart->dcell_proc != NULL)
      free(_ppart->dcell_proc);
    _ppart->dcell_proc = NULL;

    if (_ppart->dface_proc != NULL)
      free(_ppart->dface_proc);
    _ppart->dface_proc = NULL;

    if (_ppart->dface_cell != NULL)
      free(_ppart->dface_cell);
    _ppart->dface_cell = NULL;

    if (_ppart->dvtx_proc != NULL)
      free(_ppart->dvtx_proc);
    _ppart->dvtx_proc = NULL;

    if (_ppart->dpart_proc != NULL)
      free(_ppart->dpart_proc);
    _ppart->dpart_proc = NULL;

    if (_ppart->gpart_to_lproc_part != NULL)
      free(_ppart->gpart_to_lproc_part);
    _ppart->gpart_to_lproc_part = NULL;

    if (_ppart->dpart_bound != NULL)
      free(_ppart->dpart_bound);
    _ppart->dpart_bound = NULL;

    if (_ppart->ddual_graph_idx != NULL)
      free(_ppart->ddual_graph_idx);
    _ppart->ddual_graph_idx = NULL;

    if (_ppart->ddual_graph != NULL)
      free(_ppart->ddual_graph);
    _ppart->ddual_graph = NULL;

    for (int i = 0; i < _ppart->n_part; i++) {
      _part_free(_ppart->mesh_parts[i]);
      _ppart->mesh_parts[i] = NULL;
    }

    PDM_timer_free(_ppart->timer);
    _ppart->timer = NULL;

    if (_ppart->mesh_parts != NULL)
      free(_ppart->mesh_parts);
    _ppart->mesh_parts = NULL;

  }

  free (_ppart);

}



/**
 *
 * \brief Return times
 *
 * \param [in]   ppart       Pointer to \ref PDM_part object
 * \param [out]  elapsed     Elapsed times (size = 4)
 * \param [out]  cpu         Cpu times (size = 4)
 * \param [out]  cpu_user    User cpu times (size = 4)
 * \param [out]  cpu_sys     System cpu times (size = 4)
 *
 */

void PDM_part_time_get
(
 PDM_part_t  *ppart,
 double     **elapsed,
 double     **cpu,
 double     **cpu_user,
 double     **cpu_sys
)
{
  _PDM_part_t *_ppart = (_PDM_part_t *) ppart;

  *elapsed  = _ppart->times_elapsed;
  *cpu      = _ppart->times_cpu;
  *cpu_user = _ppart->times_cpu_u;
  *cpu_sys  = _ppart->times_cpu_s;
}


/**
 *
 * \brief Return statistics
 *
 * \param [in]   ppart                          Pointer to \ref PDM_part object
 * \param [out]  cells_average                  average of cells number
 * \param [out]  cells_median                   median of cells number
 * \param [out]  cells_std_deviation            standard deviation of cells number
 * \param [out]  cells_min                      minimum of cells nummber
 * \param [out]  cells_max                      maximum of cells nummber
 * \param [out]  bound_part_faces_average       average of partitioning boundary faces
 * \param [out]  bound_part_faces_median        median of partitioning boundary faces
 * \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 * \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 * \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 *
 */

void
PDM_part_stat_get
(
PDM_part_t  *ppart,
int         *cells_average,
int         *cells_median,
double      *cells_std_deviation,
int         *cells_min,
int         *cells_max,
int         *bound_part_faces_average,
int         *bound_part_faces_median,
double      *bound_part_faces_std_deviation,
int         *bound_part_faces_min,
int         *bound_part_faces_max,
int         *bound_part_faces_sum
)
{
  _PDM_part_t *_ppart = (_PDM_part_t *) ppart;
  if (_ppart->use_multipart) {
    PDM_multipart_stat_get(ppart->multipart,
                           0, // i_domain
                           cells_average,
                           cells_median,
                           cells_std_deviation,
                           cells_min,
                           cells_max,
                           bound_part_faces_average,
                           bound_part_faces_median,
                           bound_part_faces_std_deviation,
                           bound_part_faces_min,
                           bound_part_faces_max,
                           bound_part_faces_sum);
    return;
  }

  int n_rank;
  PDM_MPI_Comm_size(_ppart->comm, &n_rank);

  int *n_loc = (int *) malloc(_ppart->n_part * sizeof(int));
  int *n_tot = (int *) malloc(_ppart->dpart_proc[n_rank] * sizeof(int));

  int *s_loc = (int *) malloc(_ppart->n_part * sizeof(int));
  int *s_tot = (int *) malloc(_ppart->dpart_proc[n_rank] * sizeof(int));

  for (int i = 0; i < _ppart->n_part; i++) {
    n_loc[i] = 0;
    s_loc[i] = 0;
  }

  for (int i = 0; i < _ppart->dpart_proc[n_rank]; i++) {
    n_tot[i] = 0;
    s_tot[i] = 0;
  }

  for (int i = 0; i < _ppart->n_part; i++) {
    n_loc[i] = _ppart->mesh_parts[i]->n_cell;
    s_loc[i] = _ppart->mesh_parts[i]->n_face_part_bound;
  }

  int *n_partProc = (int *) malloc((n_rank) * sizeof(int));

  for (int i = 0; i < n_rank; i++) {
    n_partProc[i] = _ppart->dpart_proc[i+1] - _ppart->dpart_proc[i];
  }

  PDM_MPI_Allgatherv((void *) n_loc,
                     _ppart->n_part,
                     PDM_MPI_INT,
                     (void *) n_tot,
                     n_partProc,
                     _ppart->dpart_proc,
                     PDM_MPI_INT,
                     _ppart->comm);

  PDM_MPI_Allgatherv((void *) s_loc,
                     _ppart->n_part,
                     PDM_MPI_INT,
                     (void *) s_tot,
                     n_partProc,
                     _ppart->dpart_proc,
                     PDM_MPI_INT,
                     _ppart->comm);

  PDM_quick_sort_int(s_tot, 0, _ppart->dpart_proc[n_rank]-1);
  PDM_quick_sort_int(n_tot, 0, _ppart->dpart_proc[n_rank]-1);

  double   _cells_average;

  double   _bound_part_faces_average;

  *bound_part_faces_min = -1;
  *bound_part_faces_max = -1;
  *cells_min = -1;
  *cells_max = -1;
  _cells_average = 0;
  _bound_part_faces_average = 0;

  for (int i = 0; i < _ppart->dpart_proc[n_rank]; i++) {
    if (*bound_part_faces_min < 0)
      *bound_part_faces_min = s_tot[i];
    else
      *bound_part_faces_min = PDM_MIN(*bound_part_faces_min, s_tot[i]);
    if (*bound_part_faces_max < 0)
      *bound_part_faces_max = s_tot[i];
    else
      *bound_part_faces_max = PDM_MAX(*bound_part_faces_max, s_tot[i]);
    if (*cells_min < 0)
      *cells_min = n_tot[i];
    else
      *cells_min = PDM_MIN(*cells_min, n_tot[i]);
    if (*cells_max < 0)
      *cells_max = n_tot[i];
    else
      *cells_max = PDM_MAX(*cells_max, n_tot[i]);

    _cells_average += n_tot[i];
    _bound_part_faces_average += s_tot[i];
  }

  _cells_average = (_cells_average/((double) _ppart->dpart_proc[n_rank]));
  *bound_part_faces_sum = (int) _bound_part_faces_average;
  _bound_part_faces_average =
    _bound_part_faces_average/((double) _ppart->dpart_proc[n_rank]);

  *cells_average = (int) round(_cells_average);
  *bound_part_faces_average = (int) round(_bound_part_faces_average);

  *cells_std_deviation = 0.;
  *bound_part_faces_std_deviation = 0.;
  for (int i = 0; i < _ppart->dpart_proc[n_rank]; i++) {
    *cells_std_deviation += (n_tot[i] - _cells_average) * (n_tot[i] - _cells_average);
    *bound_part_faces_std_deviation += (s_tot[i] - _bound_part_faces_average) *
                                      (s_tot[i] - _bound_part_faces_average);
  }

  *cells_std_deviation = sqrt(*cells_std_deviation/_ppart->dpart_proc[n_rank]);
  *bound_part_faces_std_deviation =
    sqrt(*bound_part_faces_std_deviation/_ppart->dpart_proc[n_rank]);

  int mid = _ppart->dpart_proc[n_rank]/2;
  if (_ppart->dpart_proc[n_rank] % 2 == 1) {
    *cells_median = n_tot[mid];
    *bound_part_faces_median = s_tot[mid];
  }

  else {
    *cells_median =(int) round((n_tot[mid-1] + n_tot[mid])/2.);
    *bound_part_faces_median = (int) ((s_tot[mid-1] + s_tot[mid])/2.);
  }

  free(n_tot);
  free(s_tot);
  free(n_loc);
  free(s_loc);
  free(n_partProc);
}



/**
 *
 * \brief Free partially a \ref PDM_part_object
 *
 * \param [in]   ppart       Pointer to \ref PDM_part object
 *
 */

void
PDM_part_partial_free
(
 PDM_part_t  *ppart
)
{
  _PDM_part_t *_ppart = (_PDM_part_t *) ppart;
  PDM_MPI_Barrier(_ppart->comm);

  if (_ppart->dcell_proc != NULL)
    free(_ppart->dcell_proc);
  _ppart->dcell_proc = NULL;

  if (_ppart->dface_proc != NULL)
    free(_ppart->dface_proc);
  _ppart->dface_proc = NULL;

  if (_ppart->dvtx_proc != NULL)
    free(_ppart->dvtx_proc);
  _ppart->dvtx_proc = NULL;

  if (_ppart->dpart_proc != NULL)
    free(_ppart->dpart_proc);
  _ppart->dpart_proc = NULL;

  if (_ppart->ddual_graph_idx != NULL)
    free(_ppart->ddual_graph_idx);
  _ppart->ddual_graph_idx = NULL;

  if (_ppart->ddual_graph != NULL)
    free(_ppart->ddual_graph);
  _ppart->ddual_graph = NULL;

  // printf("PDM_part_partial_free f2 \n");
  for (int i = 0; i < _ppart->n_part; i++) {
    _part_partial_free(_ppart->mesh_parts[i]);
  }
  // printf("PDM_part_partial_free f3 \n");
  // printf("PDM_part_partial_free f8 \n");

}



#ifdef __cplusplus
}
#endif /* __cplusplus */
