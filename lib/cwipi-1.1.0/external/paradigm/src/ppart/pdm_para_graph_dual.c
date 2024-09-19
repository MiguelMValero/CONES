
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
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_para_graph_dual.h"
#include "pdm_array.h"
#include "pdm_dconnectivity_transform.h"


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

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Compress the connectivity of a graph, ie remove the multiple arcs connecting
 *        the same two nodes (if any) and remove occurence of the current line
 *
 * \param [in]    n_node            (local) number of nodes in the graph
 * \param [inout] dual_graph_idx    Node to node connectivity indexes (size=n_node+1)
 * \param [in] dual_graph_n         Original number of connected nodes (size=n_node)
 * \param [inout] dual_graph        Node to node connectivity (size=dual_graph_idx[n_node])
 *
 */
void
PDM_para_graph_compress_connectivity_dual
(
      int          n_node,
      PDM_g_num_t  shift_rank,
      PDM_g_num_t *dual_graph_idx,
const int         *dual_graph_n,
      PDM_g_num_t *dual_graph
)
{
  int idx_comp  = 0; /* Compressed index use to fill the buffer */
  int idx_block = 0; /* Index in the block to post-treat        */
  dual_graph_idx[0] = 0;
  int need_shift = 0;
  for(int i = 0; i < n_node; ++i){
    int n_cell_connect = dual_graph_n[i];
    /* Reshift next value in compressed block to avoid create a new shift */
    if(need_shift) {
      int idx_new = idx_comp;
      for(int j = idx_block; j < idx_block+n_cell_connect; ++j){
        dual_graph[idx_new++] = dual_graph[j];
      }
    }

    int end_connect         = idx_comp + n_cell_connect - 1;
    // printf(" idx_comp:: %d | end_connect:: %d \n", idx_comp, end_connect);

    int n_cell_connect_comp = PDM_inplace_unique_long(dual_graph, NULL, idx_comp, end_connect);
    // printf(" n_cell_connect:: %d | n_cell_connect_comp:: %d \n", n_cell_connect, n_cell_connect_comp);

    PDM_g_num_t g_num = i + shift_rank + 1;
    int ipos = PDM_binary_search_long(g_num, &dual_graph[idx_comp], n_cell_connect_comp);
    // printf("g_num : "PDM_FMT_G_NUM" | --> ipos::%i\n", g_num, ipos);

    if( ipos != -1 ) { /* So the global number of current line appears */
      // printf(" Suppress : "PDM_FMT_G_NUM"\n", dual_graph[idx_comp+ipos]);
      dual_graph[idx_comp+ipos] = dual_graph[idx_comp+n_cell_connect_comp-1];
      n_cell_connect_comp--;
    }

    //Once a shift is needed, need_shift must stay at one
    if(n_cell_connect_comp < n_cell_connect) {
      need_shift = 1;
    }

    dual_graph_idx[i+1] = dual_graph_idx[i] + n_cell_connect_comp;
    idx_comp  += n_cell_connect_comp;
    idx_block += n_cell_connect;

    // printf("idx_comp : %d | idx_block : %d \n", idx_comp, idx_block);
  }
}

/**
 *
 * \brief Compress the connectivity of a graph, ie remove the multiple arcs connecting
 *        the same two nodes (if any).
 *
 * \param [in]    n_node            (local) number of nodes in the graph
 * \param [inout] dual_graph_idx    Node to node connectivity indexes (size=n_node+1)
 * \param [in] dual_graph_n         Original number of connected nodes (size=n_node)
 * \param [inout] dual_graph        Node to node connectivity (size=dual_graph_idx[n_node])
 *
 */
void
PDM_para_graph_compress_connectivity
(
      int          n_node,
      int         *dual_graph_idx,
const int         *dual_graph_n,
      PDM_g_num_t *dual_graph
)
{
  int idx_comp  = 0; /* Compressed index use to fill the buffer */
  int idx_block = 0; /* Index in the block to post-treat        */
  dual_graph_idx[0] = 0;
  int need_shift = 0;
  for(int i = 0; i < n_node; ++i){
    int n_cell_connect = dual_graph_n[i];
    /* Reshift next value in compressed block to avoid create a new shift */
    if(need_shift) {
      int idx_new = idx_comp;
      for(int j = idx_block; j < idx_block+n_cell_connect; ++j){
        dual_graph[idx_new++] = dual_graph[j];
      }
    }

    int end_connect         = idx_comp + n_cell_connect - 1;
    // printf(" idx_comp:: %d | end_connect:: %d \n", idx_comp, end_connect);

    int n_cell_connect_comp = PDM_inplace_unique_long(dual_graph, NULL, idx_comp, end_connect);
    // printf(" n_cell_connect:: %d | n_cell_connect_comp:: %d \n", n_cell_connect, n_cell_connect_comp);

    //Once a shift is needed, need_shift must stay at one
    if(n_cell_connect_comp < n_cell_connect) {
      need_shift = 1;
    }

    dual_graph_idx[i+1] = dual_graph_idx[i] + n_cell_connect_comp;
    idx_comp  += n_cell_connect_comp;
    idx_block += n_cell_connect;

    // printf("idx_comp : %d | idx_block : %d \n", idx_comp, idx_block);
  }
}

/**
 *
 * \brief Compute in parallel the dual graph of an unstructured graph represented
 *        by its arc (edges of the graph) to node (vertices of the graph) connectivity.
 *        Arc and edge terminology is employed to avoid confusion with geometric entities
 *        such as vertices, edges, etc.
 *        Usually for a CFD mesh, the nodes of the graph are the cells of the mesh
 *        and the arcs of the graph are thus the faces of the mesh.
 *
 *        The dual graph computed by this function is the node to node connectivity
 *        as know as adjacency list, requested by graph partitioners.
 *
 *        Additively, this function computes the node to arc connectivity if
 *        compute_node_to_arc is true.
 *
 * \param [in]   comm               PDM_MPI communicator
 * \param [in]   graph_node_distrib distribution of nodes over the procs (size=n_rank+1)
 * \param [in]   graph_arc_distrib  distribution of arcs  over the procs (size=n_rank+1)
 * \param [in]   darc_to_node       Arc to node connectivity (size=2*dn_arc)
 * \param [out]  dual_graph_idx     Node to node connectivity indexes (size=dn_node+1)
 * \param [out]  dual_graph         Node to node connectivity (size=dual_graph_idx[dn_node])
 * \param [in]   compute_dnode_to_arc Compute or not node to arc connectivity
 * \param [out]  dnode_to_arc_idx   Node to arc connectivity indexes (size=dn_node+1)
 * \param [out]  dnode_to_arc       Node to arc connectivity (size=dnode_to_arc_idx[dn_node])
 *
 */
void
PDM_para_graph_dual_from_arc2node
(
const PDM_MPI_Comm     comm,
      PDM_g_num_t     *graph_node_distrib,
const PDM_g_num_t     *graph_arc_distrib,
const PDM_g_num_t     *darc_to_node,
      PDM_g_num_t    **dual_graph_idx,
      PDM_g_num_t    **dual_graph,
const int              compute_dnode_to_arc,
      int            **dnode_to_arc_idx,
      PDM_g_num_t    **dnode_to_arc
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_arc  = graph_arc_distrib[i_rank+1]  -  graph_arc_distrib[i_rank];
  int dn_node = graph_node_distrib[i_rank+1] -  graph_node_distrib[i_rank];
  // printf("dn_node : %d\n", dn_node);
  // printf("dn_arc  : %d\n", dn_arc);

  /*
   * We need for each node the connectivity with other nodes ( in global numbering )
   *    -> We can use part to block on arc2node to setup the correct connectivity
   *       because for each nodes we receive the contribution of each arc connectivity
   */
  PDM_g_num_t* dnode_ln_to_gn = (PDM_g_num_t *) malloc( 2 * dn_arc * sizeof(PDM_g_num_t));
  PDM_g_num_t* dopposite_node = (PDM_g_num_t *) malloc( 2 * dn_arc * sizeof(PDM_g_num_t));
  PDM_g_num_t* darc_g;

  int* arc_strid  = (int *) malloc(sizeof(int) * 2 * dn_arc);
  int* node_strid = (int *) malloc(sizeof(int) * 2 * dn_arc);

  if(compute_dnode_to_arc){
    darc_g = (PDM_g_num_t *) malloc( 2 * dn_arc * sizeof(PDM_g_num_t));
  }

  PDM_g_num_t shift_arc_g   = graph_arc_distrib[i_rank]+1; // Entre 1 et N
  // printf(" shift_arc_g = %i \n", shift_arc_g);
  int dn_arc_int    = 0;
  int idx_data_arc  = 0;
  int idx_data_node = 0;
  for(int i_arc = 0; i_arc < dn_arc; ++i_arc){

    PDM_g_num_t g_node1 = darc_to_node[2*i_arc  ];
    PDM_g_num_t g_node2 = darc_to_node[2*i_arc+1];

    dnode_ln_to_gn[dn_arc_int] = g_node1;

    if(compute_dnode_to_arc){
      arc_strid[dn_arc_int]  = 1;
      darc_g[idx_data_arc++] = shift_arc_g + i_arc ;
    }

    if(g_node2 > 0){
      node_strid[dn_arc_int++]        = 1;
      dopposite_node[idx_data_node++] = g_node2;

      if(compute_dnode_to_arc){
        arc_strid[dn_arc_int]  = 1;
        darc_g[idx_data_arc++] = -1*(shift_arc_g + i_arc);
      }

      node_strid[dn_arc_int]          = 1;
      dnode_ln_to_gn[dn_arc_int++]    = g_node2;
      dopposite_node[idx_data_node++] = g_node1;

    } else {
      node_strid[dn_arc_int++] = 0;
    }

  }

  if( 0 == 1){
    printf("dnode_ln_to_gn::");
    for(int i = 0; i < dn_arc_int; ++i){
      printf(PDM_FMT_G_NUM" ", dnode_ln_to_gn[i]);
    }
    printf("\n");

    printf("node_strid::");
    for(int i = 0; i < dn_arc_int; ++i){
      printf("%d ", node_strid[i]);
    }
    printf("\n");

    printf("dopposite_node::");
    for(int i = 0; i < idx_data_node ; ++i){
      printf(PDM_FMT_G_NUM" ", dopposite_node[i]);
    }
    printf("\n");

    printf("idx_data_node ::%d\n", idx_data_node );
    printf("dn_arc_int    ::%d\n", dn_arc_int);
    printf("dn_arc        ::%d\n", dn_arc);
  }


  node_strid     = realloc(node_strid,     dn_arc_int    * sizeof(int)         );
  dnode_ln_to_gn = realloc(dnode_ln_to_gn, dn_arc_int    * sizeof(PDM_g_num_t) );
  dopposite_node = realloc(dopposite_node, idx_data_node * sizeof(PDM_g_num_t) );

  /*
   * Initialize part_to_block for the computation of node_node
   *    -> Si on force une distribution utilisateur on devra passer le cell_distribution
   *           --> Semble necessaire pour parMetis mais pas scotch
   */
  PDM_part_to_block_t *ptb_dual =
   PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM,
                             1.,
                             &dnode_ln_to_gn,
                             graph_node_distrib,
                             &dn_arc_int,
                             1,
                             comm);

  const int n_node_block = PDM_part_to_block_n_elt_block_get (ptb_dual);
  const PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get(ptb_dual);

  /*
   * We exchange the dnode_ln_to_gn (almost asynchrone - send strid is synchrone )
   */
  int req_id_dual = PDM_part_to_block_async_exch(ptb_dual,
                                                 sizeof(PDM_g_num_t),
                                                 PDM_STRIDE_VAR_INTERLACED,
                                                 1,
                                                 &node_strid,
                                       (void **) &dopposite_node);



  int req_id_node_arc = -1;
  if(compute_dnode_to_arc){

    req_id_node_arc = PDM_part_to_block_async_exch(ptb_dual,
                                                   sizeof(PDM_g_num_t),
                                                   PDM_STRIDE_VAR_INTERLACED,
                                                   1,
                                                   &arc_strid,
                                         (void **) &darc_g);
  }

  /*
   * Reception asynchrone
   */
  // int*            node_node_n = NULL;
  int*             recv_strid = NULL;
  PDM_g_num_t* recv_node_node = NULL;
  PDM_part_to_block_async_wait(ptb_dual, req_id_dual);
  int n_data_recv = PDM_part_to_block_asyn_get_raw(ptb_dual,
                                                   req_id_dual,
                                                  &recv_strid,
                                        (void **) &recv_node_node);

  /*
   * The data is recv in raw format - We need to post-treat them as int
   */
  *dual_graph = (PDM_g_num_t *) malloc( n_data_recv * sizeof(PDM_g_num_t));
  PDM_g_num_t* _dual_graph     = (PDM_g_num_t  *) *dual_graph;

  /*
   * Allocate and setup convenient pointeur
   */

  *dual_graph_idx      = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_node_block+1));
  PDM_g_num_t* _dual_graph_idx = *dual_graph_idx;

  int* node_node_n   = (int *) malloc( (n_node_block+1) * sizeof(int)); /* Suralloc */

  /*
   * Count - In our case we know that recv_strid == 1 or 0
   */
  for(int i_recv = 0; i_recv < n_node_block+1; ++i_recv) {
    _dual_graph_idx[i_recv] = 0;
    node_node_n    [i_recv] = 0;
  }

  /* Stride can be 0 or 1 */
  for(int i_recv = 0; i_recv < n_data_recv; ++i_recv) {
    // if(recv_strid[i_recv] != 0){
      int ielmt = blk_gnum[i_recv] - graph_node_distrib[i_rank] - 1;
      // printf("ielmt::[%i] --> [%i]\n",(int)blk_gnum[i_recv], (int)ielmt );
      _dual_graph_idx[ielmt+1] += recv_strid[i_recv];
    // }
  }

  /* Index computation */
  for(int i_recv = 1; i_recv < n_node_block; ++i_recv) {
    _dual_graph_idx[i_recv+1] += _dual_graph_idx[i_recv];
  }

  /* Panic verbose */
  if( 0 == 1 ){
    printf("node_node_idx::");
    for(int i_recv = 0; i_recv < n_node_block+1; ++i_recv) {
      printf(" "PDM_FMT_G_NUM, _dual_graph_idx[i_recv]);
    }
    printf("\n");
  }

  /*
   * Fill buffer
   */
  int idx_recv = 0;
  for(int i_recv = 0; i_recv < n_data_recv; ++i_recv) {
    if(recv_strid[i_recv] != 0){
      int ielmt = blk_gnum[i_recv] - graph_node_distrib[i_rank] - 1;
      _dual_graph[_dual_graph_idx[ielmt] + node_node_n[ielmt]++] = recv_node_node[idx_recv++];
    }
  }

  free(recv_strid);
  free(recv_node_node);

  /*
   * Each block can have multiple same cell, we need to compress them
   *   We do it inplace cause unique will always have inferior size of complete array
   *
   */
  PDM_para_graph_compress_connectivity_dual(n_node_block,
                                            graph_node_distrib[i_rank],
                                            _dual_graph_idx,
                                            node_node_n,
                                            _dual_graph);

  free(node_node_n);

  /*
   * Realloc
   */
  *dual_graph = (PDM_g_num_t *) realloc(*dual_graph, sizeof(PDM_g_num_t) * _dual_graph_idx[n_node_block] );
  _dual_graph = (PDM_g_num_t *) *dual_graph;

  // For now we can change it later
  // printf("n_node_block::%d \n", n_node_block);
  // printf("dn_node      ::%d \n", dn_node );
  assert(n_node_block == dn_node );

  /*
   * Panic verbose
   */
  if( 0 == 1 ){
    printf("n_node_block:: %d \n", n_node_block);
    for(int i = 0; i < n_node_block; ++i){
      printf(" _dual_graph_idx = "PDM_FMT_G_NUM" ---> \n", _dual_graph_idx[i]);
      for(int i_data = _dual_graph_idx[i]; i_data < _dual_graph_idx[i+1]; ++i_data){
        // printf("%d ", _dual_graph[i_data]);
        printf("\t _dual_graph[%d] = "PDM_FMT_G_NUM" \n", i_data, _dual_graph[i_data]);
      }
      printf("\n");
    }
  }

  /* Scoth begin at 0 even if we put base value to 1 */
  for(int i = 0; i < n_node_block; ++i){
    for(int i_data = _dual_graph_idx[i]; i_data < _dual_graph_idx[i+1]; ++i_data){
      _dual_graph[i_data] = _dual_graph[i_data] - 1;
    }
  }

  /*
   * Async recv for cell_face
   */
  if(compute_dnode_to_arc){

    int*   recv_node2arc_strid = NULL;
    PDM_g_num_t* recv_node2arc = NULL;
    PDM_part_to_block_async_wait(ptb_dual, req_id_node_arc);
    int n_data_cf_recv = PDM_part_to_block_asyn_get_raw(ptb_dual,
                                                        req_id_node_arc,
                                                       &recv_node2arc_strid,
                                             (void **) &recv_node2arc);
    free(recv_node2arc_strid); // Always 1

    /*
     * Post treatment
     */
    *dnode_to_arc      = (PDM_g_num_t *) malloc(  n_data_cf_recv  * sizeof(PDM_g_num_t));
    *dnode_to_arc_idx  = PDM_array_zeros_int(n_node_block+1);
    int* node_to_arc_n = PDM_array_zeros_int(n_node_block+1);

    /* Short-cut */
    int*         _dnode_to_arc_idx = *dnode_to_arc_idx;
    PDM_g_num_t* _dnode_to_arc     = *dnode_to_arc;

    for(int i_recv = 0; i_recv < n_data_cf_recv; ++i_recv) {
      int ielmt = blk_gnum[i_recv] - graph_node_distrib[i_rank] - 1;
      // printf("ielmt = %i \n", ielmt);
      // printf("ielmt::[%d] --> [%d] - [%d] \n",blk_gnum[i_recv], ielmt, recv_node2arc_strid[i_recv]);
      _dnode_to_arc_idx[ielmt+1] += 1;
    }

    /* Index computation */
    for(int i = 1; i < n_node_block; ++i) {
      _dnode_to_arc_idx[i+1] += _dnode_to_arc_idx[i];
    }


    /*
     * Fill buffer -  Cas particulier ou recv_stride == 1
     */
    for(int i_recv = 0; i_recv < n_data_cf_recv; ++i_recv) {
      int ielmt = blk_gnum[i_recv] - graph_node_distrib[i_rank] - 1;
      // printf(" ielmt :: %d | i_recv :: %d | _dnode_to_arc_idx :: %d | node_to_arc_n :: %d \n", ielmt, i_recv, _dnode_to_arc_idx[ielmt], node_to_arc_n[ielmt]);
      _dnode_to_arc[_dnode_to_arc_idx[ielmt] + node_to_arc_n[ielmt]++] = recv_node2arc[i_recv];
    }

    free(recv_node2arc);

    if( 0 == 1 ){
      printf("n_node_block:: %d \n", n_node_block);
      int idx_block = 0;
      for(int i = 0; i < n_node_block; ++i){
        printf(" node_to_arc_n = %d ---> ", node_to_arc_n[i]);
        for(int i_data = 0; i_data < node_to_arc_n[i]; ++i_data){
          printf(PDM_FMT_G_NUM" ", _dnode_to_arc[idx_block]);
          idx_block++;
        }
        printf("\n");
      }
    }

    // _dnode_to_arc_idx[0] = 0;
    // for(int i_node = 0; i_node < n_node_block; ++i_node){
    //   _dnode_to_arc_idx[i_node+1] = _dnode_to_arc_idx[i_node] + node_to_arc_n[i_node];
    // }
    free(node_to_arc_n);

  }

  // abort();

  /*
   * Exchange is done we can free direclty memory
   */
  free(dnode_ln_to_gn);
  free(arc_strid);
  free(node_strid);
  free(dopposite_node);
  if(compute_dnode_to_arc){
    free(darc_g);
  }

  PDM_part_to_block_free (ptb_dual);
}

/**
 *
 * \brief Compute in parallel the dual graph of an unstructured graph represented
 *        by its node (vertices of the graph) to arc (edges of the graph) connectivity.
 *        Arc and edge terminology is employed to avoid confusion with geometric entities
 *        such as vertices, edges, etc.
 *        Usually for a CFD mesh, the nodes of the graph are the cells of the mesh
 *        and the arcs of the graph are thus the faces of the mesh.
 *
 * \param [in]   comm               PDM_MPI communicator
 * \param [in]   graph_node_distrib distribution of nodes over the procs (size=n_rank+1)
 * \param [in]   graph_arc_distrib  distribution of arcs  over the procs (size=n_rank+1)
 * \param [in]   dnode_arc_idx      Node to arc connectivity indexes (size=dn_node+1)
 * \param [in]   dnode_arc          Node to arc connectivity (size=dnode_to_arc_idx[dn_node])
 * \param [out]  dual_graph_idx     Node to node connectivity indexes (size=dn_node+1)
 * \param [out]  dual_graph         Node to node connectivity (size=dual_graph_idx[dn_node])
 */
void
PDM_para_graph_dual_from_node2arc
(
const PDM_MPI_Comm     comm,
      PDM_g_num_t     *graph_node_distrib,
const PDM_g_num_t     *graph_arc_distrib,
const int             *dnode_arc_idx,
const PDM_g_num_t     *dnode_arc,
      PDM_g_num_t    **dual_graph_idx,
      PDM_g_num_t    **dual_graph
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_arc  = graph_arc_distrib[i_rank+1]  -  graph_arc_distrib[i_rank];
  int dn_node = graph_node_distrib[i_rank+1] -  graph_node_distrib[i_rank];

  /* We reconstruct a arc_to_node connectivity from the node_to_arc input.
  This can be done with a part_to_block where the lntogn is the node_to_arc
  connectivity and the send value are the global number of the nodes.

  Each proc will received one (boundary) or two nodes in global num for his arcs
  and will thus be able to construct its arc_to_node
  */

  PDM_g_num_t* node_g       = (PDM_g_num_t *) malloc(dnode_arc_idx[dn_node] * sizeof(PDM_g_num_t));
  PDM_g_num_t* arc_ln_to_gn = (PDM_g_num_t *) malloc(dnode_arc_idx[dn_node] * sizeof(PDM_g_num_t));

  PDM_g_num_t shift_node_g = graph_node_distrib[i_rank]+1; // Entre 1 et N
  for (int i_node = 0; i_node < dn_node; i_node++) {
    for (int i_arc = dnode_arc_idx[i_node]; i_arc < dnode_arc_idx[i_node+1]; i_arc++) {
      int g_sign = PDM_SIGN(dnode_arc[i_arc]);
      node_g[i_arc] = g_sign * (i_node + shift_node_g);
      arc_ln_to_gn[i_arc] = PDM_ABS(dnode_arc[i_arc]);
    }
  }

  if(0 == 1) {
    PDM_printf("reverted cell face :");
    for (int i = 0; i < dnode_arc_idx[dn_node]; i++)
      printf(" "PDM_FMT_G_NUM, node_g[i]);
    printf("\n");
  }

  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
                            &arc_ln_to_gn,
                             graph_arc_distrib,
                 (int *)    &dnode_arc_idx[dn_node],
                             1,
                             comm);

  int* send_stride = PDM_array_const_int(dnode_arc_idx[dn_node], 1);
  free(arc_ln_to_gn);

  int        *recv_stride = NULL;
  PDM_g_num_t  *recv_data = NULL;

  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &send_stride,
                (void **) &node_g,
                          &recv_stride,
                (void **) &recv_data);

  int n_recv_block = PDM_part_to_block_n_elt_block_get(ptb);
  assert(n_recv_block == dn_arc);

  if (0 == 1) {
    int idx = 0;
    PDM_printf("[%d]Recv stride (%d): ", i_rank, n_recv_block);
    for (int i = 0; i < n_recv_block; i++)
      PDM_printf(" %d", recv_stride[i]);
    PDM_printf("\n[%d]Recv data:", i_rank);
    for (int i = 0; i < n_recv_block; i++) {
      for (int j = idx; j < idx + recv_stride[i]; j++)
        PDM_printf(" %d", recv_data[j]);
      idx += recv_stride[i];
    }
    PDM_printf("\n");
  }

  PDM_g_num_t *darc_to_node = (PDM_g_num_t *) malloc( 2*dn_arc * sizeof(PDM_g_num_t));

  int idx_recv_data = 0;
  for (int i_arc = 0; i_arc < dn_arc; i_arc++) {
    darc_to_node[2*i_arc] = recv_data[idx_recv_data++];
    if (recv_stride[i_arc] == 2) {
      darc_to_node[2*i_arc+1] = recv_data[idx_recv_data++];
      int sign_left = PDM_SIGN(darc_to_node[2*i_arc]);
      // Swap left and right face if left face is negative
      if (sign_left > 0) {
        darc_to_node[2*i_arc+1] = PDM_ABS(darc_to_node[2*i_arc+1]);
      } else {
        PDM_g_num_t tmp_face    = PDM_ABS(darc_to_node[2*i_arc  ]);
        darc_to_node[2*i_arc  ] = PDM_ABS(darc_to_node[2*i_arc+1]); //normally useless
        darc_to_node[2*i_arc+1] = tmp_face;
      }
    } else {
      darc_to_node[2*i_arc+1] = 0;
    }
  }

  if (0 == 1) {
    PDM_printf("[%d] Generated arc_to_node ::", i_rank);
    for (int i = 0; i < dn_arc; i++){
      PDM_printf(" %d %d", darc_to_node[2*i], darc_to_node[2*i+1]);
    }
    PDM_printf("\n");
  }

  PDM_part_to_block_free(ptb);
  free(node_g);
  free(send_stride);
  free(recv_stride);
  free(recv_data);

  /*
   * Now we have a arc_to_node connectivity, we can call graph_dual_from_arc2node
   */
  PDM_para_graph_dual_from_arc2node(comm,
                                    graph_node_distrib,
                                    graph_arc_distrib,
                                    darc_to_node,
                                    dual_graph_idx,
                                    dual_graph,
                                    0,
                                    NULL,
                                    NULL);
  free(darc_to_node);
}

/**
 *
 * \brief Compute in parallel the dual graph of an unstructured graph represented
 *        by its node (vertices of the graph) to arc (edges of the graph) connectivity.
 *        Arc and edge terminology is employed to avoid confusion with geometric entities
 *        such as vertices, edges, etc.
 *        Usually for a CFD mesh, the nodes of the graph are the cells of the mesh
 *        and the arcs of the graph are thus the faces of the mesh.
 *
 * \param [in]   comm               PDM_MPI communicator
 * \param [in]   graph_node_distrib distribution of nodes over the procs (size=n_rank+1)
 * \param [in]   graph_arc_distrib  distribution of arcs  over the procs (size=n_rank+1)
 * \param [in]   dnode_arc_idx      Node to arc connectivity indexes (size=dn_node+1)
 * \param [in]   dnode_arc          Node to arc connectivity (size=dnode_to_arc_idx[dn_node])
 * \param [out]  dual_graph_idx     Node to node connectivity indexes (size=dn_node+1)
 * \param [out]  dual_graph         Node to node connectivity (size=dual_graph_idx[dn_node])
 */
void
PDM_para_graph_dual_from_combine_connectivity
(
const PDM_MPI_Comm   comm,
const PDM_g_num_t   *cell_distrib,
const PDM_g_num_t   *face_distrib,
      PDM_g_num_t   *vtx_distrib,
const int           *dcell_face_idx,
const PDM_g_num_t   *dcell_face,
const int           *dface_vtx_idx,
const PDM_g_num_t   *dface_vtx,
      PDM_g_num_t  **dual_graph_idx,
      PDM_g_num_t  **dual_graph
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_cell = cell_distrib[i_rank+1] - cell_distrib[i_rank];
  // int dn_face = face_distrib[i_rank+1] - face_distrib[i_rank];

  int*         dcell_vtx_idx;
  PDM_g_num_t* dcell_vtx;

  /*
   *  Call generic function to deduce the induce connectivity
   */
  PDM_deduce_combine_connectivity(comm,
                                  cell_distrib,
                                  face_distrib,
                                  dcell_face_idx,
                                  dcell_face,
                                  dface_vtx_idx,
                                  dface_vtx,
                                  1,
                ( int         **) &dcell_vtx_idx,
                ( PDM_g_num_t **) &dcell_vtx);

  // int* dcell_vtx_n = (int *) malloc( dn_cell * sizeof(int));
  // for(int i_entity = 0; i_entity < dn_cell; ++i_entity) {
  //   dcell_vtx_n[i_entity] = dcell_vtx_idx[i_entity+1] - dcell_vtx_idx[i_entity];
  // }
  // PDM_log_trace_array_int (dcell_vtx_n  , dn_cell               , "dcell_vtx_n::");
  // PDM_log_trace_array_int (dcell_vtx_idx, dn_cell+1             , "dcell_vtx_idx::");
  // PDM_log_trace_array_long(dcell_vtx    , dcell_vtx_idx[dn_cell], "dcell_vtx::");
  // free(dcell_vtx_n);

  int*         dvtx_cell_idx;
  PDM_g_num_t* dvtx_cell;
  PDM_dconnectivity_transpose(comm,
                               cell_distrib,
                               vtx_distrib,
                               dcell_vtx_idx,
                               dcell_vtx,
                               1,
                               &dvtx_cell_idx,
                               &dvtx_cell);

  // int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];
  // int* dvtx_cell_n = (int *) malloc( dn_vtx * sizeof(int));
  // for(int i_entity = 0; i_entity < dn_vtx; ++i_entity) {
  //   dvtx_cell_n[i_entity] = dvtx_cell_idx[i_entity+1] - dvtx_cell_idx[i_entity];
  // }
  // PDM_log_trace_array_int (dvtx_cell_n  , dn_vtx               , "dvtx_cell_n::");
  // PDM_log_trace_array_int (dvtx_cell_idx, dn_vtx+1             , "dvtx_cell_idx::");
  // PDM_log_trace_array_long(dvtx_cell    , dvtx_cell_idx[dn_vtx], "dvtx_cell::");
  // free(dvtx_cell_n);

  /*
   * Call the standard fonction : arc = vtx , node = cell
   */
  PDM_deduce_combine_connectivity_dual(comm,
                                       cell_distrib,
                                       vtx_distrib,
                                       dcell_vtx_idx,
                                       dcell_vtx,
                                       dvtx_cell_idx,
                                       dvtx_cell,
                                       1,
                                       dual_graph_idx,
                                       dual_graph);

  free(dcell_vtx);
  free(dcell_vtx_idx);
  free(dvtx_cell_idx);
  free(dvtx_cell);


  // Realloc
  *dual_graph = (PDM_g_num_t *) realloc(*dual_graph, (*dual_graph_idx)[dn_cell] * sizeof(PDM_g_num_t));

  PDM_g_num_t* _dual_graph_idx = *dual_graph_idx;
  PDM_g_num_t* _dual_graph     = *dual_graph;

  // -> Shift by one
  for(int i_entity = 0; i_entity < _dual_graph_idx[dn_cell]; ++i_entity) {
    _dual_graph[i_entity] = _dual_graph[i_entity] - 1;
    // assert(_dual_graph[i_entity] >= 0);
    // assert(_dual_graph[i_entity] <  cell_distrib[n_rank]-1);
  }
  // int* _dual_graph_n = (int *) malloc( dn_cell * sizeof(int));
  // for(int i_entity = 0; i_entity < dn_cell; ++i_entity) {
  //   _dual_graph_n[i_entity] = _dual_graph_idx[i_entity+1] - _dual_graph_idx[i_entity];
  // }

  // PDM_log_trace_array_int (_dual_graph_n  , dn_cell                 , "PDM_para_graph_dual_from_combine_connectivity::_dual_graph_n::");
  // PDM_log_trace_array_int (_dual_graph_idx, dn_cell+1               , "PDM_para_graph_dual_from_combine_connectivity::dual_graph_idx::");
  // PDM_log_trace_array_long(_dual_graph    , _dual_graph_idx[dn_cell], "PDM_para_graph_dual_from_combine_connectivity::dual_graph::");

  /*
   * Patch
   */
  // int* _dual_comp_graph_idx = (int *) malloc( dn_cell * sizeof(int));
  // PDM_g_num_t* _dual_comp_graph = (PDM_g_num_t *) malloc( _dual_graph_idx[dn_cell] * sizeof(PDM_g_num_t));

  // _dual_comp_graph_idx[0] = 0;
  // for(int i_entity = 0; i_entity < dn_cell; ++i_entity) {
  //   PDM_g_num_t g_num = i_entity + cell_distrib[i_rank];
  //   _dual_comp_graph_idx[i_entity+1] = _dual_comp_graph_idx[i_entity];
  //   for(int j = _dual_graph_idx[i_entity]; j < _dual_graph_idx[i_entity+1]; ++j) {
  //     if(_dual_graph[j] != g_num-1) {
  //       _dual_comp_graph[_dual_comp_graph_idx[i_entity+1]++] = _dual_graph[j];
  //     }
  //   }
  // }

  // PDM_log_trace_array_int (_dual_graph_n  , dn_cell                 , "PDM_para_graph_dual_from_combine_connectivity::_dual_graph_n::");
  // PDM_log_trace_array_int (_dual_comp_graph_idx, dn_cell+1               , "PDM_para_graph_dual_from_combine_connectivity::_dual_comp_graph_idx::");
  // PDM_log_trace_array_long(_dual_comp_graph    , _dual_comp_graph_idx[dn_cell], "PDM_para_graph_dual_from_combine_connectivity::_dual_comp_graph::");

  /*
   * Realloc
   */
  // free(_dual_graph);
  // free(_dual_graph_idx);

  // _dual_comp_graph = (PDM_g_num_t *) realloc(_dual_comp_graph, _dual_comp_graph_idx[dn_cell] * sizeof(PDM_g_num_t));

  // *dual_graph_idx = _dual_comp_graph_idx;
  // *dual_graph     = _dual_comp_graph;


  // free(_dual_graph_n);
}

/**
 *
 * \brief Call the chosen graph partitioner to split the dual graph
 *
 * \param [in]   split_method       Choice of the graph partitioner
 * \param [in]   graph_node_distrib distribution of nodes over the procs (size=n_rank+1)
 * \param [in]   dual_graph_idx     Node to node connectivity indexes (size=dn_node+1)
 * \param [in]   dual_graph         Node to node connectivity (size=dual_graph_idx[dn_node])
 * \param [in]   node_weight        Weight associated to each node of the graph or NULL
 * \param [in]   arc_weight         Weight associated to each arc of the graph or NULL
 * \param [in]   n_part             Total number of partitions to produce
 * \param [in]   part_fraction      Fraction of (weighted) vertex wanted on each part (Metis only)
                                    or NULL for homogeneous sizes (size = n_part)
 * \param [out]  node_part_id       Attributed partition number for each node (size=dn_node)
 * \param [in]   comm               PDM_MPI communicator
 */
void
PDM_para_graph_split
(
const PDM_split_dual_t  split_method,
const PDM_g_num_t      *graph_node_distrib,
const PDM_g_num_t      *dual_graph_idx,
const PDM_g_num_t      *dual_graph,
const int              *node_weight,
const int              *arc_weight,
const int               n_part,
const double           *part_fraction,
      int              *node_part_id,
const PDM_MPI_Comm      comm
)
{
  int i_rank;
  int n_rank;

  PDM_UNUSED(part_fraction);
  PDM_UNUSED(dual_graph_idx);
  PDM_UNUSED(dual_graph);
  PDM_UNUSED(node_weight);
  PDM_UNUSED(arc_weight);
  PDM_UNUSED(n_part);

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_elmt = graph_node_distrib[i_rank+1] - graph_node_distrib[i_rank];

  // printf("dn_elmt = %i \n", dn_elmt);
  PDM_array_reset_int(node_part_id, dn_elmt, 0);

  switch (split_method) {

    case PDM_SPLIT_DUAL_WITH_PARMETIS:
    {
      #ifndef PDM_HAVE_PARMETIS
        PDM_printf("PPART error : ParMETIS unavailable\n");
        exit(1);
      #else

        // Define metis properties
        int numflag    = 0;   /* C or Fortran numbering (C = 0)                    */
        int wgtflag    = 0;   /* Indicate if graph is weighted                     */
        int ncon       = 1;   /* Number of weight to consider per vertex           */
        double *ubvec;        /* Imbalance tolerance for vertex weights            */
        double *tpwgts;       /* Fraction of (weighted) vertex wanted on each part */
        int edgecut;          /* Number of edges cutted by metis                   */

        ubvec = (double *) malloc(ncon * sizeof(double));
        for (int i = 0; i < ncon; i++) {
          ubvec[i] = 1.05;
        }

        tpwgts = (double *) malloc(ncon * n_part * sizeof(double));
        if (part_fraction != NULL) {
          for (int i = 0; i < n_part; i++) {
            for (int j = 0; j < ncon; j++) {
              tpwgts[ncon*i + j] = part_fraction[i];
            }
          }
        }
        else {
          for (int i = 0; i < ncon * n_part; i++) {
            tpwgts[i] = (double) (1./n_part);
          }
        }

        if (arc_weight != NULL) {
          if (node_weight != NULL) wgtflag = 3;     //Weights on both the vertices and edges
          else wgtflag = 1;                         //Weights on the edges only
        }
        else if (node_weight != NULL) wgtflag = 2;  //Weights on the vertices only

        // printf("PDM_ParMETIS_dpart %d | %d \n", n_part, dn_elmt);
        PDM_ParMETIS_V3_PartKway (graph_node_distrib,
                                  dual_graph_idx,
                                  dual_graph,
                                  node_weight,
                                  arc_weight,
                                 &wgtflag,
                                 &numflag,
                                 &ncon,
                                 &n_part,
                                  tpwgts,
                                  ubvec,
                                 &edgecut,
                                  node_part_id,
                                  comm);
        // printf("PDM_ParMETIS_dpart %d | %d  END\n", n_part, dn_elmt);

        free(ubvec);
        free(tpwgts);
      #endif
        break;
    }

    case PDM_SPLIT_DUAL_WITH_PTSCOTCH:
    {
      #ifndef PDM_HAVE_PTSCOTCH
        PDM_printf("PPART error : PT-Scotch unavailable\n");
        exit(1);
      #else
        int check = 1;
        // printf("PDM_SCOTCH_dpart %d | %d \n", n_part, dn_elmt);
        PDM_SCOTCH_dpart (dn_elmt,
                          dual_graph_idx,
                          dual_graph,
                          node_weight,
                          arc_weight,
                          check,
                          comm,
                          n_part,
                          node_part_id);
        // printf("PDM_SCOTCH_dpart end \n");
      #endif
      break;
    }

    default:
      PDM_printf("PPART error : unknown partioning choice '%i'\n", split_method);
      exit(1);
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
