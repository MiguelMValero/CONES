/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"

#include "pdm_domain_interface.h"
#include "pdm_part_domain_interface.h"
#include "pdm_domain_interface_priv.h"
#include "pdm_part_domain_interface_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/
static int _unique_pairs(int n_pairs, PDM_g_num_t *ids, int *dom_ids) {
  //Rebuild n_occurences
  int n_read = 0;
  PDM_g_num_t last = -1;
  for (int i = 0; i < n_pairs; i++) {
    if (ids[2*i] != last) {
      n_read++;
      last = ids[2*i];
    }
  }
  int *n_occurences = PDM_array_const_int(n_read, 1);
  n_read = -1;
  last   = -1;
  for (int i = 0; i < n_pairs; i++) {
    if (ids[2*i] != last) {
      n_read++;
      last = ids[2*i];
    }
    else {
      n_occurences[n_read]++;
    }
  }
  n_read++; //Add Last
  int max_occur = 0;
  for (int j = 0; j < n_read; j++)
    max_occur = PDM_MAX(max_occur, n_occurences[j]);

  PDM_g_num_t *working_array     = (PDM_g_num_t *) malloc(max_occur*sizeof(PDM_g_num_t));
  PDM_g_num_t *backup_array_gnum = (PDM_g_num_t *) malloc(2*max_occur*sizeof(PDM_g_num_t));
  int         *backup_array_int  = (int *)         malloc(2*max_occur*sizeof(int));
  int         *ordering_array    = (int *)         malloc(max_occur*sizeof(int));

  int start = 0;
  int rewrite_start = 0;
  for (int j = 0; j < n_read; j++) {
    //Fill working array
    for (int k = 0; k < n_occurences[j]; k++)
      working_array[k] = ids[2*(start+k)+1];
    int n_unique = PDM_inplace_unique_long2(working_array, ordering_array, 0, n_occurences[j]-1);
    //Copy into array
    memcpy(backup_array_gnum, &ids[2*start],     2*n_occurences[j]*sizeof(PDM_g_num_t));
    memcpy(backup_array_int,  &dom_ids[2*start], 2*n_occurences[j]*sizeof(int));
    for (int k = 0; k < n_occurences[j]; k++) {
      int pos = ordering_array[k];
      ids[2*(rewrite_start+pos)]   = backup_array_gnum[2*k];
      ids[2*(rewrite_start+pos)+1] = backup_array_gnum[2*k+1];
      dom_ids[2*(rewrite_start+pos)]   = backup_array_int[2*k];
      dom_ids[2*(rewrite_start+pos)+1] = backup_array_int[2*k+1];
    }
    start += n_occurences[j];
    rewrite_start += n_unique;
  }

  free(n_occurences);
  free(working_array);
  free(ordering_array);
  free(backup_array_gnum);
  free(backup_array_int);
  return rewrite_start;
}

static
PDM_g_num_t*
_per_block_offset
(
 int           n_block,
 int          *sizes,
 PDM_MPI_Comm  comm
)
{
  PDM_g_num_t *sizes_as_gn = (PDM_g_num_t *) malloc(n_block*sizeof(PDM_g_num_t));
  for (int i = 0; i < n_block; i ++) {
    sizes_as_gn[i] = sizes[i];
  }
  PDM_g_num_t *per_block_offset = (PDM_g_num_t *) malloc((n_block+1) * sizeof(PDM_g_num_t));
  per_block_offset[0] = 0;
  PDM_MPI_Allreduce(sizes_as_gn, &per_block_offset[1], n_block, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_array_accumulate_gnum(per_block_offset, n_block+1);
  free(sizes_as_gn);
  return per_block_offset;
}

// Translate interfaces list to a graph representation
static int _interface_to_graph
(
  const int                    n_interface,
  PDM_domain_interface_mult_t  multidomain_intrf,
  int                         *interface_dn,
  PDM_g_num_t                **interface_ids,
  int                        **interface_dom,
  int                        **graph_idx,
  PDM_g_num_t                **graph_ids,
  int                        **graph_dom,
  PDM_MPI_Comm                 comm
)
{
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  // Step 0 : retrieve some data. We need (to offset gnums)
  //   - the number of involved blocks
  //   - the max id occuring in each block
  int n_domain = -1;
  if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
    int max_domain_loc = 0;
    for (int itrf = 0; itrf < n_interface; itrf++) {
      for (int k = 0; k < 2*interface_dn[itrf]; k++) {
        max_domain_loc = PDM_MAX(max_domain_loc, interface_dom[itrf][k]);
      }
    }
    PDM_MPI_Allreduce(&max_domain_loc, &n_domain, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  }
  else {
    for (int itrf = 0; itrf < n_interface; itrf++) {
      n_domain = PDM_MAX(n_domain, interface_dom[itrf][0]);
      n_domain = PDM_MAX(n_domain, interface_dom[itrf][1]);
    }
  }
  n_domain++; //Because domain numbering start at 0

  PDM_g_num_t *max_per_domain_loc = PDM_array_const_gnum(n_domain, 0);
  PDM_g_num_t *max_per_domain     = (PDM_g_num_t *) malloc((n_domain+1) * sizeof(PDM_g_num_t));
  for (int itrf = 0; itrf < n_interface; itrf++) {
    int dom    = -1;
    int domopp = -1;
    if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      dom    = interface_dom[itrf][0];
      domopp = interface_dom[itrf][1];
    }
    for (int k = 0; k < interface_dn[itrf]; k++) {
      if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        dom    = interface_dom[itrf][2*k];
        domopp = interface_dom[itrf][2*k+1];
      }
      max_per_domain_loc[dom   ] = PDM_MAX(max_per_domain_loc[dom   ], interface_ids[itrf][2*k  ]);
      max_per_domain_loc[domopp] = PDM_MAX(max_per_domain_loc[domopp], interface_ids[itrf][2*k+1]);
    }
  }
  max_per_domain[0] = 0;
  PDM_MPI_Allreduce(max_per_domain_loc, &max_per_domain[1], n_domain, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_array_accumulate_gnum(max_per_domain, n_domain+1);
  if (0 == 1) {
    PDM_log_trace_array_long(max_per_domain, n_domain+1, "max per domain");
  }
  free(max_per_domain_loc);
  
  // Prepare first PtB with multiple partitions.
  // Use (shifted) ids as gnum and send tuple (shited) id, opp_id
  PDM_g_num_t **interface_ids_shifted = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t*));
  PDM_g_num_t **send_data             = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t*));
  double      **weight                = (double      **) malloc(n_interface * sizeof(double*     ));
  int         **stride_one            = (int         **) malloc(n_interface * sizeof(int*        ));
  int          *interface_dn_twice    = (int          *) malloc(n_interface * sizeof(int         ));
  for (int itrf = 0; itrf < n_interface; itrf++) {
    stride_one[itrf]            = (int         *) malloc(2*interface_dn[itrf]*sizeof(int        ));
    interface_ids_shifted[itrf] = (PDM_g_num_t *) malloc(2*interface_dn[itrf]*sizeof(PDM_g_num_t));
    send_data[itrf]             = (PDM_g_num_t *) malloc(2*interface_dn[itrf]*sizeof(PDM_g_num_t));
    weight[itrf]                = (double      *) malloc(2*interface_dn[itrf]*sizeof(double     ));
    interface_dn_twice[itrf]    = 2*interface_dn[itrf];
    int dom    = -1;
    int domopp = -1;
    if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      dom    = interface_dom[itrf][0];
      domopp = interface_dom[itrf][1];
    }
    for (int k = 0; k < interface_dn[itrf]; k++) {
      if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        dom    = interface_dom[itrf][2*k];
        domopp = interface_dom[itrf][2*k+1];
      }
      interface_ids_shifted[itrf][2*k]   = interface_ids[itrf][2*k] + max_per_domain[dom];
      interface_ids_shifted[itrf][2*k+1] = interface_ids[itrf][2*k+1] + max_per_domain[domopp];
      send_data[itrf][2*k]   = interface_ids[itrf][2*k+1] + max_per_domain[domopp];
      send_data[itrf][2*k+1] = interface_ids[itrf][2*k] + max_per_domain[dom];
      weight[itrf][2*k]   = 1.;
      weight[itrf][2*k+1] = 1.;
      stride_one[itrf][2*k]   = 1;
      stride_one[itrf][2*k+1] = 1;
    }
    if (0 == 1) {
      log_trace("Interface %d\n", itrf);
      PDM_log_trace_array_long(interface_ids_shifted[itrf], 2*interface_dn[itrf], "  shifted gnum");
      PDM_log_trace_array_long(send_data[itrf], 2*interface_dn[itrf], "  send");
    }
  }
  
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      interface_ids_shifted,
                                                      weight,
                                                      interface_dn_twice,
                                                      n_interface,
                                                      comm);
  // Save distribution & gnum from first PtB. We will use it for following PtBs
  int n_gnum = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *distri = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *gnum   = (PDM_g_num_t *) malloc(n_gnum    * sizeof(PDM_g_num_t));
  memcpy(gnum,   PDM_part_to_block_block_gnum_get(ptb),        n_gnum*sizeof(PDM_g_num_t));
  memcpy(distri, PDM_part_to_block_distrib_index_get(ptb), (n_rank+1)*sizeof(PDM_g_num_t));

  int         *recv_stride = NULL;
  PDM_g_num_t *recv_data   = NULL;
  int n_connected_l = PDM_part_to_block_exch(ptb,
                                             sizeof(PDM_g_num_t),
                                             PDM_STRIDE_VAR_INTERLACED,
                                             -1,
                                             stride_one,
                                   (void **) send_data,
                                             &recv_stride,
                                   (void **) &recv_data);
  if (0 == 1) {
    PDM_log_trace_array_long(PDM_part_to_block_block_gnum_get(ptb), n_gnum, "gnum");
    PDM_log_trace_array_int(recv_stride, n_gnum, "recv stride");
    PDM_log_trace_array_long(recv_data, n_connected_l, "recv data");
  }

  PDM_part_to_block_free(ptb);
  for (int itrf = 0; itrf < n_interface; itrf++) {
    free(stride_one           [itrf]);
    free(send_data            [itrf]);
    free(weight               [itrf]);
    free(interface_ids_shifted[itrf]);
  }
  free(stride_one);
  free(weight);
  free(send_data);
  free(interface_dn_twice);
  free(interface_ids_shifted);

  int n_connected;
  PDM_MPI_Allreduce(&n_connected_l, &n_connected, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  // log_trace("Initial size of graph : %d \n", n_connected);


  /* After first exchange, we received for each gnum a list (usually of size one) of connected
   * gnum. The idea is to consider this received list as a partition, and to send to each
   * entity of this partition : the neighbors in the received list (if stride was > 1) + the original gnum.
   *
   * After some iterations this will group the related ids together
  */

  int n_connected_prev = 0;
  while(n_connected_prev != n_connected) {
    // log_trace("\nSize of graph : %d (prev %d), start new it\n", n_connected, n_connected_prev);
    n_connected_prev = n_connected;

    ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                    1.,
                                   &recv_data,
                                    distri,
                                   &n_connected_l,
                                    1,
                                    comm);

    
    int *send_stride = (int *) malloc(n_connected_l*sizeof(int));
    int w_idx = 0;
    int n_data = 0;
    for (int k = 0; k < n_gnum; k++) {
      for (int j = 0; j < recv_stride[k]; j++) {
        send_stride[w_idx++] = recv_stride[k];
        n_data += recv_stride[k];
      }
    }
    assert (w_idx == n_connected_l);
    PDM_g_num_t *send_data2 = (PDM_g_num_t *) malloc(n_data*sizeof(PDM_g_num_t));
    w_idx = 0;
    int r_idx = 0;
    for (int k = 0; k < n_gnum; k++) {
      for (int j = 0; j < recv_stride[k]; j++) {
        //Add gnum
        send_data2[w_idx++] = gnum[k];
        //Add others
        for (int i = 0; i < j; i++)
          send_data2[w_idx++] = recv_data[r_idx + i];
        for (int i = j+1; i < recv_stride[k]; i++)
          send_data2[w_idx++] = recv_data[r_idx + i];
      }
      r_idx += recv_stride[k];
    }
    assert (r_idx == n_connected_l);
    assert (w_idx == n_data);
    if (0 == 1) {
      PDM_log_trace_array_int(send_stride, n_connected_l, "  send stride");
      PDM_log_trace_array_long(send_data2, n_data, "  send_data2");
    }

    int         *recv_stride_next = NULL;
    PDM_g_num_t *recv_data_next   = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           &send_stride,
                 (void **) &send_data2,
                           &recv_stride_next,
                 (void **) &recv_data_next);

    free(send_stride);
    free(send_data2);
    PDM_part_to_block_free(ptb);
    
    // Post treat recv data to remove duplicated per gnum and count size of graph
    int start = 0;
    n_connected_l = 0;
    for (int i = 0; i < n_gnum; i ++) {
      int n_unique = PDM_inplace_unique_long(recv_data_next, NULL, start, start+recv_stride_next[i]-1);
      for(int k = 0; k < n_unique; ++k) {
        recv_data_next[n_connected_l+k] = recv_data_next[start+k];
      }
      start += recv_stride_next[i];
      recv_stride_next[i] = n_unique;
      n_connected_l += n_unique;
    }

    PDM_MPI_Allreduce(&n_connected_l, &n_connected, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    if (0 == 1) {
      PDM_log_trace_array_long(gnum, n_gnum, "  gnum");
      PDM_log_trace_array_int(recv_stride_next, n_gnum, "  recv stride");
      PDM_log_trace_array_long(recv_data_next, n_connected_l, "  recv data");
      log_trace("  Total size of graph : %d \n", n_connected);
    }

    free(recv_stride);
    free(recv_data); // To free after PTB because it was used as lngn
    recv_data = recv_data_next;
    recv_stride = recv_stride_next;
  }
  
  // When iteration are completed, all the connections are known by every id.
  // Last step is to compress the graph and to redistribute it
  // To do that we take for each group of related id the min of it as lngn
  int n_keys = 0;
  n_connected_l = 0;
  int r_idx = 0;
  /* is_key_gr : -1 if gnum is not a key (only min of group is the key),
   *              1 if gnum is a key and gnum is not included in received ids
   *              0 otherwise */
  int *is_key_gr = PDM_array_const_int(n_gnum, 1);
  for (int k = 0; k < n_gnum; k++) {
    for (int j = 0; j < recv_stride[k]; j++) {
      if (recv_data[r_idx+j] == gnum[k]) {
        is_key_gr[k] = 0;
      }
      else if (recv_data[r_idx+j] < gnum[k]) {
        is_key_gr[k] = -1;
        break;
      }
    }
    if (is_key_gr[k] != -1) {
      n_keys++;
      n_connected_l += recv_stride[k] + is_key_gr[k];
    }
    r_idx  += recv_stride[k];
  }

  PDM_g_num_t *lngn_gr        = (PDM_g_num_t *) malloc(n_keys*sizeof(PDM_g_num_t));
  int         *send_stride_gr = (int         *) malloc(n_keys*sizeof(int        ));
  double      *weight_gr      = (double      *) malloc(n_keys*sizeof(double     ));
  PDM_g_num_t *send_data_gr   = (PDM_g_num_t *) malloc(n_connected_l *sizeof(PDM_g_num_t));
  int w_idx = 0;
  int w_idx2 = 0;
  r_idx = 0;
  for (int k = 0; k < n_gnum; k++) {
    if (is_key_gr[k] != -1) {
      lngn_gr[w_idx]        = gnum[k];
      send_stride_gr[w_idx] = recv_stride[k] + is_key_gr[k]; //Include gnum in send data (if needed) so we have directly graph
      weight_gr[w_idx]      = (double) (recv_stride[k] + is_key_gr[k]);
      w_idx++;
      if (is_key_gr[k] == 1) {
        send_data_gr[w_idx2++] = gnum[k];
      }
      memcpy(&send_data_gr[w_idx2], &recv_data[r_idx], recv_stride[k]*sizeof(PDM_g_num_t));
      w_idx2 += recv_stride[k];
    }
    r_idx += recv_stride[k];
  }
  if (0 == 1) {
    log_trace("Build graph\n");
    PDM_log_trace_array_long(lngn_gr, n_keys, "  keys graph");
    PDM_log_trace_array_int(send_stride_gr, n_keys, "  send stride");
    PDM_log_trace_array_long(send_data_gr, n_connected_l, "  send data");
  }
  
  // Data of previous iteration is not usefull anymore
  free(gnum);
  free(distri);
  free(recv_stride);
  free(recv_data);
  
  //TODO In fact we just want to do a block to block, but PTB + weights compute distribution for us
  //and we are lazy
  ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                 PDM_PART_TO_BLOCK_POST_MERGE,
                                 1.,
                                &lngn_gr,
                                &weight_gr,
                                &n_keys,
                                 1,
                                 comm);
  int graph_dn = PDM_part_to_block_n_elt_block_get(ptb);

  int         *graph_size = NULL;
  PDM_g_num_t *graph_gnum = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                         &send_stride_gr,
               (void **) &send_data_gr,
                         &graph_size,
               (void **) &graph_gnum);

  PDM_part_to_block_free(ptb);
  free(send_stride_gr);
  free(send_data_gr);
  free(weight_gr);
  free(lngn_gr);
  free(is_key_gr);

  int* _graph_idx = PDM_array_new_idx_from_sizes_int(graph_size, graph_dn);
  int *_graph_dom = (int *) malloc(_graph_idx[graph_dn]*sizeof(int));

  if (0 == 1) {
    PDM_log_trace_array_int(graph_size, graph_dn, "  recv stride");
    PDM_log_trace_array_long(graph_gnum, _graph_idx[graph_dn], "  recv data");
  }

  // Retrieve domain and local gnum in domain
  for (int i = 0; i < _graph_idx[graph_dn]; i++) {
    _graph_dom[i] = PDM_binary_search_gap_long(graph_gnum[i]-1, max_per_domain, n_domain+1);
    graph_gnum[i] -= max_per_domain[_graph_dom[i]];
  }
  free(graph_size);
  free(max_per_domain);

  *graph_idx = _graph_idx;
  *graph_ids =  graph_gnum;
  *graph_dom = _graph_dom;
  return graph_dn;
}

static int _extract_and_shift_jn_faces
(
 int           n_domain,
 int          *dn_face,
 int          *dn_vtx,
 int           n_interface,
 int          *interfaces_size,
 PDM_g_num_t **interface_face_ids,
 int         **interface_domains_ids,
 int         **dface_vtx_idx,
 PDM_g_num_t **dface_vtx,
 int         **face_vtx_both_idx,
 PDM_g_num_t **face_vtx_both,
 int         **dextract_face_group_id,
 int         **dextract_face_dom_id,
 PDM_g_num_t **dextract_face_id,
 PDM_MPI_Comm  comm
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *face_per_block_offset = _per_block_offset(n_domain, dn_face, comm);
  PDM_g_num_t *vtx_per_block_offset  = _per_block_offset(n_domain, dn_vtx,  comm);
  
  int n_face_join = 0; // Each interface comes with a pair of faces
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    n_face_join += 2*interfaces_size[i_interface];
  }

  PDM_g_num_t *_dextract_face_id_tmp       = (PDM_g_num_t *) malloc(n_face_join * sizeof(PDM_g_num_t));
  // Also transport some data to the extracted faces
  int         *_dextract_face_group_id_tmp = (int *) malloc(n_face_join*sizeof(int));
  int         *_dextract_face_dom_id_tmp   = (int *) malloc(n_face_join*sizeof(int));

  int idx = 0;
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {

    PDM_g_num_t *_interface_ids = interface_face_ids[i_interface]; //Shortcuts
    int         *_interface_dom = interface_domains_ids[i_interface];
    for (int i_pair = 0; i_pair < interfaces_size[i_interface]; i_pair++) {
      int i_domain_cur = _interface_dom[2*i_pair];
      int i_domain_opp = _interface_dom[2*i_pair+1];

      // _dextract_face_id_tmp[idx]         = _interface_ids[2*i_pair] + face_per_block_offset[i_domain_cur];
      _dextract_face_id_tmp[idx]         = PDM_ABS(_interface_ids[2*i_pair]) + face_per_block_offset[i_domain_cur];
      _dextract_face_dom_id_tmp  [idx]   = i_domain_cur;
      _dextract_face_group_id_tmp[idx++] = i_interface;

      // _dextract_face_id_tmp[idx]         = _interface_ids[2*i_pair+1] + face_per_block_offset[i_domain_opp];
      _dextract_face_id_tmp[idx]         = PDM_ABS(_interface_ids[2*i_pair+1]) + face_per_block_offset[i_domain_opp];
      _dextract_face_dom_id_tmp  [idx]   = i_domain_opp;
      _dextract_face_group_id_tmp[idx++] = i_interface;
    }
  }
  assert (idx == n_face_join);
  
  // Multi gnum is not equilibrated, we have to redistribute it but we want to keep the face/face_opp groups
  PDM_g_num_t *cur_distri   = PDM_compute_entity_distribution(comm, n_face_join/2);
  PDM_g_num_t *ideal_distri = PDM_compute_uniform_entity_distribution(comm, cur_distri[n_rank]);
  PDM_block_to_block_t *btb = PDM_block_to_block_create(cur_distri, ideal_distri, comm);

  //dextract_face_id will be used to get the face->vtx of extracted faces, and also
  //to build the hash table later so return it
  PDM_block_to_block_exch(btb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          2,
                          NULL,
                          _dextract_face_id_tmp,
                          NULL,
              (void **)   dextract_face_id);
  PDM_block_to_block_exch(btb,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          2,
                          NULL,
                          _dextract_face_group_id_tmp,
                          NULL,
              (void **)  dextract_face_group_id);
  PDM_block_to_block_exch(btb,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          2,
                          NULL,
                          _dextract_face_dom_id_tmp,
                          NULL,
              (void **)  dextract_face_dom_id);

  
  PDM_block_to_block_free(btb);
  // Update n_face_join before freeing distribution
  n_face_join = 2*(ideal_distri[i_rank+1]-ideal_distri[i_rank]);
  free(ideal_distri);
  free(cur_distri);
  free(_dextract_face_id_tmp);
  free(_dextract_face_group_id_tmp);
  free(_dextract_face_dom_id_tmp);

  if (0 == 1) {
    PDM_log_trace_array_long(face_per_block_offset, n_domain+1, "face_per_block_offset :: ");
    PDM_log_trace_array_long(vtx_per_block_offset,  n_domain+1, "vtx_per_block_offset :: ");
    PDM_log_trace_array_long(*dextract_face_id, n_face_join, "dextract_face_id :: ");
  }

  PDM_g_num_t **all_face_distribution = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t*));
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    all_face_distribution[i_domain] = PDM_compute_entity_distribution(comm, dn_face[i_domain]);
  }

  PDM_multi_block_to_part_t *mptb = PDM_multi_block_to_part_create(face_per_block_offset,
                                                                   n_domain,
                                            (const PDM_g_num_t **) all_face_distribution,
                                            (const PDM_g_num_t **) dextract_face_id,
                                                                  &n_face_join,
                                                                   1,
                                                                   comm);
  //Prepare data to send : face -> vtx connectivity 
  int         **face_vtx_n       = (int         **) malloc(n_domain * sizeof(int*));
  PDM_g_num_t **face_vtx_shifted = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t*));
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    face_vtx_n[i_domain]       = (int         *) malloc(dn_face[i_domain] * sizeof(int));
    face_vtx_shifted[i_domain] = (PDM_g_num_t *) malloc(dface_vtx_idx[i_domain][dn_face[i_domain]] * sizeof(PDM_g_num_t));

    for (int i_face = 0; i_face < dn_face[i_domain]; i_face++) {
      face_vtx_n[i_domain][i_face] = dface_vtx_idx[i_domain][i_face+1] - dface_vtx_idx[i_domain][i_face];
      for (int j = dface_vtx_idx[i_domain][i_face]; j < dface_vtx_idx[i_domain][i_face+1]; j++) {
        face_vtx_shifted[i_domain][j] = dface_vtx[i_domain][j] + vtx_per_block_offset[i_domain];
      }
    }
  }

  int         **part_stride = NULL;
  PDM_g_num_t **part_data   = NULL;
  PDM_multi_block_to_part_exch2(mptb,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR_INTERLACED,
                                face_vtx_n,
                 (void **)      face_vtx_shifted,
                               &part_stride,
                 (void ***)    &part_data);

  *face_vtx_both_idx = PDM_array_new_idx_from_sizes_int(part_stride[0], n_face_join);
  *face_vtx_both     = part_data[0];

  free(part_data);
  free(part_stride[0]);
  free(part_stride);

  PDM_multi_block_to_part_free(mptb);

  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    free(all_face_distribution[i_domain]);
    free(face_vtx_n[i_domain]);
    free(face_vtx_shifted[i_domain]);
  }
  free(all_face_distribution);
  free(face_vtx_n);
  free(face_vtx_shifted);

  free(face_per_block_offset);
  free(vtx_per_block_offset);

  return n_face_join;
}
static int _generate_edge_face
(
 int            n_face,
 int           *face_vtx_idx,
 PDM_g_num_t   *face_vtx,
 PDM_g_num_t  **dedge_distrib,
 int          **dedge_vtx_idx,
 PDM_g_num_t  **dedge_vtx,
 int          **dedge_face_idx,
 PDM_g_num_t  **dedge_face,
 PDM_MPI_Comm  comm
)
{

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  
  // 1. Get the number of unique vertex
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &face_vtx,
                                                      NULL,
                                                     &face_vtx_idx[n_face],
                                                      1,
                                                      comm);
  int n_vtx;
  int n_vtx_loc = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_MPI_Allreduce(&n_vtx_loc, &n_vtx, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  PDM_part_to_block_free(ptb);

  // 2. Get (unmerged) edges
  int n_edge_unmerged = face_vtx_idx[n_face];
  int n_elmt_current  = 0;
  int n_edge_current  = 0;

  PDM_g_num_t *face_distri         = PDM_compute_entity_distribution(comm, n_face);
  PDM_g_num_t *tmp_edge_face       = (PDM_g_num_t *) malloc(n_edge_unmerged       * sizeof(PDM_g_num_t));
  int         *tmp_parent_elmt_pos = (int         *) malloc(n_edge_unmerged       * sizeof(int        ));
  int         *tmp_edge_vtx_idx    = (int         *) malloc((n_edge_unmerged + 1) * sizeof(int        ));
  PDM_g_num_t *tmp_edge_vtx        = (PDM_g_num_t *) malloc(2*n_edge_unmerged     * sizeof(PDM_g_num_t));

  tmp_edge_vtx_idx[0] = 0;
  PDM_poly2d_decomposes_edges(n_face,
                             &n_elmt_current,
                             &n_edge_current,
                              face_distri[i_rank],
                              -1,
                              face_vtx,
                              face_vtx_idx,
                              tmp_edge_vtx_idx, //Numéro de sommet des edges
                              tmp_edge_vtx,
                              tmp_edge_face,         // Numéro des edges pour une face
                              NULL,
                              NULL,
                              tmp_parent_elmt_pos);
  assert(n_edge_current == n_edge_unmerged);

  /*PDM_log_trace_array_long(tmp_dface_edge, n_edge_unmerged, "dface_edge :: ");*/
  /*PDM_log_trace_connectivity_long(tmp_dface_edge_vtx_idx, tmp_dface_edge_vtx, n_edge_unmerged, "dface_edge :: ");*/

  // 3. Merge shared edges
  int dn_edge;
  PDM_generate_entitiy_connectivity_raw(comm,
                                        n_vtx,
                                        n_edge_unmerged,
                                        tmp_edge_face,
                                        tmp_edge_vtx_idx,
                                        tmp_edge_vtx,
                                       &dn_edge,
                                        dedge_distrib,
                                        dedge_vtx_idx,
                                        dedge_vtx,
                                        dedge_face_idx,
                                        dedge_face);
  free(tmp_parent_elmt_pos);
  free(face_distri);
  return dn_edge;
}

static int _match_internal_edges
(
 int            dn_edge,
 PDM_g_num_t   *dedge_distrib,
 int           *dedge_face_idx,
 PDM_g_num_t   *dedge_face,
 PDM_g_num_t   *dedge_face_join,
 PDM_g_num_t   *dedge_face_join_opp,
 PDM_g_num_t  **dedge_gnum,
 PDM_g_num_t  **dedge_gnum_opp,
 PDM_MPI_Comm  comm
)
{
  PDM_UNUSED(dedge_face);
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  
  // 0. Count the number of internal edges
  int dn_internal_edge = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {
    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    dn_internal_edge += (int) (n_face_this_edge > 1);
  }
  if (0 == 1) {
    log_trace("dn internal edges is %i \n", dn_internal_edge);
    log_trace("dn external edges is %i \n", dn_edge - dn_internal_edge);
  }

  // 1. Build hash keys
  PDM_g_num_t *key_ln_to_gn = (PDM_g_num_t *) malloc(dn_internal_edge * sizeof(PDM_g_num_t)); 
  int         *stride_one   = (int         *) malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_two   = (int         *) malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_four  = (int         *) malloc(dn_internal_edge * sizeof(int        ));
  double      *weight       = (double      *) malloc(dn_internal_edge * sizeof(double     ));

  PDM_g_num_t *data_send_connect    = (PDM_g_num_t *) malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_edge_g_num = (PDM_g_num_t *) malloc(  dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_group      = (PDM_g_num_t *) malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_sens       = (PDM_g_num_t *) malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_face_g_num = (PDM_g_num_t *) malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));

  int i_int_edge = 0;
  int idx_write2 = 0;
  int idx_write4 = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {

    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    if (n_face_this_edge == 1) {
      continue;
    }
    assert (n_face_this_edge == 2);

    stride_one [i_int_edge] = 1;
    stride_two [i_int_edge] = 2;
    stride_four[i_int_edge] = 4;
    weight[i_int_edge] = 1.;
    //Retrive domain id using group id of any of two faces
    data_send_edge_g_num[i_int_edge] = dedge_distrib[i_rank] + i_edge + 1;

    int key = 0;
    for(int j = dedge_face_idx[i_edge]; j < dedge_face_idx[i_edge+1]; ++j) { //Do it for the two faces data
      key += (dedge_face_join[j] + dedge_face_join_opp[j]);

      //data_send_face_g_num[idx_write2]   = dedge_face[j];
      //data_send_sens      [idx_write2++] = dedge_face_group_sens[j];
      idx_write2++;

      data_send_connect[idx_write4++] = dedge_face_join    [j];

      data_send_connect[idx_write4++] = dedge_face_join_opp[j];
    }
    key_ln_to_gn[i_int_edge] = key;

    i_int_edge++;
  }
  assert(idx_write2 == 2*dn_internal_edge);
  assert(idx_write4 == 4*dn_internal_edge);

  // 2. Exchange data over hash key
  //Attention, pb d'équilibrage car les clés sont réparties vers la fin ... Un proc risque
  // de se retrouver avec tt les clés
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                    (PDM_g_num_t **) &key_ln_to_gn,
                                                     &weight,
                                                     &dn_internal_edge,
                                                      1,
                                                      comm);
  // Get protocol data
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  // This one will contain 2 for non conflicting hash key, more otherwise
  int *gnum_n_occurences = PDM_part_to_block_block_gnum_count_get(ptb);

  int gnum_n_occurences_tot = 0;
  for (int k = 0; k < blk_size; k++) {
    gnum_n_occurences_tot += gnum_n_occurences[k];
  }
  
  int         *unused_recv_stride = NULL;
  PDM_g_num_t *blk_edge_g_num     = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int **)  &stride_one,
                           (void **) &data_send_edge_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_edge_g_num);
  free(unused_recv_stride); // Same as gnum_n_occurences 
  assert (exch_size == gnum_n_occurences_tot);

  /*
  PDM_g_num_t* blk_data_sens   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int **)  &stride_two,
                           (void **) &data_send_sens,
                                     &unused_recv_stride,
                           (void **) &blk_data_sens);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences 
  assert (exch_size == 2*gnum_n_occurences_tot);
  */

  /*
  PDM_g_num_t* blk_data_face_g_num   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int **)  &stride_two,
                           (void **) &data_send_face_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_data_face_g_num);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences
  assert (exch_size == 2*gnum_n_occurences_tot);
  */

  PDM_g_num_t *blk_data_connect = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                          (int **)  &stride_four,
                          (void **) &data_send_connect,
                                    &unused_recv_stride,
                          (void **) &blk_data_connect);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);

  if (0 == 1) {
    PDM_log_trace_array_long(key_ln_to_gn, dn_internal_edge, "key_ln_to_gn :: ");
    PDM_log_trace_array_int(gnum_n_occurences   , blk_size               , "gnum_n_occurences   :: ");
    PDM_log_trace_array_long(blk_edge_g_num     , gnum_n_occurences_tot  , "blk_edge_g_num      :: ");
    //PDM_log_trace_array_long(blk_data_face_g_num, 2*gnum_n_occurences_tot, "blk_data_face_g_num :: ");
    //PDM_log_trace_array_long(blk_data_sens      , 2*gnum_n_occurences_tot, "blk_data_sens       :: ");
    PDM_log_trace_array_long(blk_data_connect   , 4*gnum_n_occurences_tot, "blk_data_connect    :: ");
  }


  free(key_ln_to_gn        );
  free(data_send_connect   );
  free(data_send_group     );
  free(data_send_sens      );
  free(data_send_face_g_num);
  free(data_send_edge_g_num);
  free(stride_one          );
  free(stride_two          );
  free(stride_four         );
  free(weight              );



  // 3. Post treatemement : resolve conflicting keys
  PDM_g_num_t *results_edge     = (PDM_g_num_t *) malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *results_edge_opp = (PDM_g_num_t *) malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));

  int n_max_entity_per_key = 0;
  for(int i = 0; i < blk_size; ++i) {
    n_max_entity_per_key = PDM_MAX(gnum_n_occurences   [i], n_max_entity_per_key);
  }

  int *already_treat   = (int *) malloc(n_max_entity_per_key * sizeof(int));
  int *same_entity_idx = (int *) malloc(n_max_entity_per_key * sizeof(int));
  // int *sens_entity     = (int *) malloc(n_max_entity_per_key * sizeof(int));


  int idx  = 0;
  int idx_w = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {

    int n_matching_edge = gnum_n_occurences[i_key];

    /* Reset */
    PDM_array_reset_int(already_treat, n_matching_edge, -1);

    /* Loop over all entitys in conflict and sort all */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {

      //Each internal edge comes with 2 faces -> 4 data
      int beg    = 4*(idx+i_entity);

      // Caution inplace sort !!!!
      PDM_sort_long(&blk_data_connect[beg], NULL, 4);

      if(0 == 1) {
        PDM_log_trace_array_long(&blk_data_connect[beg], 4, "blk_data_connect (sort) :: ");
      }
    }

    /*
     *  Identify pair or invalid other //todo : shortcut only 2 occurences or let for debug
     */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {
      int beg1    = 4*(idx+i_entity);

      int n_same         = 0;
      int i_entity2_same = -1;

      if(already_treat[i_entity] == 1) {
        continue;
      }

      for(int i_entity2 = i_entity+1; i_entity2 < n_matching_edge; ++i_entity2) {

        if(already_treat[i_entity2] == -1) {
          int beg2    = 4*(idx+i_entity2);

          if(!PDM_array_are_equal_gnum(&blk_data_connect[beg1], &blk_data_connect[beg2], 4)) {
            continue;
          }

          already_treat[i_entity2] = 1;
          same_entity_idx[n_same++] = i_entity2;
          i_entity2_same = i_entity2;

        }
      } /* End for i_entity2 */
      assert(n_same <= 1);


      if (n_same == 1) {
        int edge_idx     = (idx+i_entity);
        int edge_idx_opp = (idx+i_entity2_same);

        // Set data for edge
        results_edge    [idx_w] = blk_edge_g_num[edge_idx];
        results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx_opp];

        // Set data for opposite edge
        results_edge    [idx_w] = blk_edge_g_num[edge_idx_opp];
        results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx];
      }

      already_treat[i_entity] = 1;
    }
    idx  += n_matching_edge;
  }
  free(already_treat);
  free(same_entity_idx);
  // Some pairs can be still unresolved, eg if a edge is internal from one interface point of view but
  // external for the other
  int rsvd_gnum_n_occurences_tot = idx_w;
  results_edge     = (PDM_g_num_t *) realloc(results_edge,     rsvd_gnum_n_occurences_tot*sizeof(PDM_g_num_t));
  results_edge_opp = (PDM_g_num_t *) realloc(results_edge_opp, rsvd_gnum_n_occurences_tot*sizeof(PDM_g_num_t));


  free(blk_edge_g_num);
  free(blk_data_connect);
  PDM_part_to_block_free(ptb); // Needed until here for gnum_n_occurences

  if (0 == 1) {
    log_trace("Conflict resolved, gnum are\n");
    PDM_log_trace_array_long(results_edge, rsvd_gnum_n_occurences_tot, "edge gnum ::");
    PDM_log_trace_array_long(results_edge_opp, rsvd_gnum_n_occurences_tot, "edge gnum opp::");
  }



  // 4. Send back result on edge distribution
                       ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_NOTHING,
                                                       1.,
                                                      &results_edge,
                                                       dedge_distrib,
                                                      &rsvd_gnum_n_occurences_tot,
                                                       1,
                                                       comm);

  int resolved_dn_internal_edge = PDM_part_to_block_n_elt_block_get(ptb);
  *dedge_gnum = (PDM_g_num_t *) malloc(resolved_dn_internal_edge * sizeof(PDM_g_num_t));

  PDM_g_num_t *dedge_gnum_tmp = PDM_part_to_block_block_gnum_get(ptb);
  memcpy(*dedge_gnum, dedge_gnum_tmp, resolved_dn_internal_edge*sizeof(PDM_g_num_t));

  PDM_part_to_block_exch(ptb,
                        sizeof(PDM_g_num_t),
                        PDM_STRIDE_CST_INTERLACED,
                        1,
                        NULL,
             (void **) &(results_edge_opp),
                        NULL,
              (void **) dedge_gnum_opp);
  PDM_part_to_block_free(ptb);

  free(results_edge);
  free(results_edge_opp);

  return resolved_dn_internal_edge;
}

static void _match_all_edges_from_faces
(
  int          dn_face,
  int         *face_edge_idx,
  PDM_g_num_t *face_edge,
  PDM_g_num_t *face_edge_wopp,
  PDM_g_num_t *pedge_vtx,
  PDM_g_num_t *p_all_vtx,
  PDM_g_num_t *p_all_vtx_opp
)
{
  // Avoid multiple malloc using max size
  int max_face_len = 0;
  for (int i_face = 0; i_face < dn_face/2; i_face++) {
    max_face_len = PDM_MAX(max_face_len, face_edge_idx[2*i_face+1] - face_edge_idx[2*i_face]);
  }
  PDM_g_num_t *ordered_edge     = (PDM_g_num_t *) malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_edge_opp = (PDM_g_num_t *) malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_vtx      = (PDM_g_num_t *) malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_vtx_opp  = (PDM_g_num_t *) malloc(max_face_len * sizeof(PDM_g_num_t));

  // Avec la construction des faces de bord, on a des paires faces / face opp
  assert (dn_face%2 == 0);
  assert (face_edge_idx[dn_face]%2 == 0);
  int glob_idx     = 0;
  for (int i_face = 0; i_face < dn_face/2; i_face++) {
    int start_idx     = face_edge_idx[2*i_face];
    int start_idx_opp = face_edge_idx[2*i_face+1];
    int face_len = face_edge_idx[2*i_face+1] - face_edge_idx[2*i_face];
    if (0 == 1) {
      log_trace("iface %d :\n", i_face);
      log_trace("face data\n");
      // PDM_log_trace_array_long(&face_vtx[face_edge_idx[i_face]], face_len, "face_vtx");
      PDM_log_trace_array_long(&face_edge[start_idx], face_len, "face_edge");
      PDM_log_trace_array_long(&pedge_vtx[2*start_idx], 2*face_len, "edge vertices");
      PDM_log_trace_array_long(&face_edge_wopp[start_idx], face_len, "face_edge_wopp");
      log_trace("face opp data\n");
      // PDM_log_trace_array_long(&face_vtx_opp[face_edge_idx[i_face]], face_len, "face_vtx");
      PDM_log_trace_array_long(&face_edge[start_idx_opp], face_len, "face_opp_edge_opp");
      PDM_log_trace_array_long(&pedge_vtx[2*start_idx_opp], 2*face_len, "edge vertices opp");
    }
    
    // Search any received edge (we should have at least one, but not necessary the first
    // since we can have aliases edges
    PDM_g_num_t opp_edge_key = 0;
    int idx = -1;
    int opp_idx = -1;
    int found = 0;
    while (!found && idx < face_len) {
      opp_edge_key = face_edge_wopp[start_idx + 1 + idx++];
      if (opp_edge_key != 0) { //Skip undetermined edges
        opp_idx = -1;
        PDM_g_num_t candidate = 0;
        while (candidate != opp_edge_key && opp_idx < face_len) {
          candidate = face_edge[start_idx_opp + 1 + opp_idx++];
          candidate = PDM_ABS(candidate);
        }
        found = (int) (candidate == opp_edge_key); //Break first while if OK
      }
    }
    // log_trace("Match pos %d and %d using shared edge %d\n", idx, opp_idx, opp_edge_key);
    assert (found);

    // Find starting point using the two edge indices and signs
    int sign     = PDM_SIGN(face_edge[start_idx + idx]);
    int opp_sign = PDM_SIGN(face_edge[start_idx_opp + opp_idx]);

    int next_vtx, next_vtx_opp, cur_vtx, cur_vtx_opp;
    if (sign == 1) {
      cur_vtx  = pedge_vtx[2*(start_idx + idx)];
      next_vtx = pedge_vtx[2*(start_idx + idx)+1];
    }
    else {
      cur_vtx  = pedge_vtx[2*(start_idx + idx)+1];
      next_vtx = pedge_vtx[2*(start_idx + idx)];
    }
    if (opp_sign == 1) {//Invert looping order for opposite face
      cur_vtx_opp  = pedge_vtx[2*(start_idx_opp + opp_idx)+1];
      next_vtx_opp = pedge_vtx[2*(start_idx_opp + opp_idx)];
    }
    else {
      cur_vtx_opp  = pedge_vtx[2*(start_idx_opp + opp_idx)];
      next_vtx_opp = pedge_vtx[2*(start_idx_opp + opp_idx)+1];
    }

    for (int i = 0; i < face_len; i++) {
      //Fill
      // log_trace("Cur vtx %d and opp %d \n", cur_vtx, cur_vtx_opp);
      ordered_edge[i]     = PDM_ABS(face_edge[start_idx+idx]);
      ordered_edge_opp[i] = PDM_ABS(face_edge[start_idx_opp+opp_idx]);
      ordered_vtx[i]      = cur_vtx;
      ordered_vtx_opp[i]  = cur_vtx_opp;
      //This is for face
      for (int j = 0; j < face_len; j++) {
        if (j != idx) {
          int vtx1 = pedge_vtx[2*(start_idx + j)];
          int vtx2 = pedge_vtx[2*(start_idx + j)+1];
          if (vtx1 == next_vtx) {
            idx = j;
            cur_vtx  = vtx1;
            next_vtx = vtx2;
            break;
          }
          else if (vtx2 == next_vtx) {
            idx = j;
            cur_vtx  = vtx2;
            next_vtx = vtx1;
            break;
          }
        }
      }
      //This is for opposite face
      for (int j = 0; j < face_len; j++) {
        if (j != opp_idx) {
          int vtx1 = pedge_vtx[2*(start_idx_opp+j)];
          int vtx2 = pedge_vtx[2*(start_idx_opp+j)+1];
          if (vtx1 == next_vtx_opp) {
            opp_idx = j;
            cur_vtx_opp  = vtx1;
            next_vtx_opp = vtx2;
            break;
          }
          else if (vtx2 == next_vtx_opp) {
            opp_idx = j;
            cur_vtx_opp  = vtx2;
            next_vtx_opp = vtx1;
            break;
          }
        }
      }
    }
    if (0 == 1) {
      PDM_log_trace_array_long(ordered_edge,     face_len, "ordered edges");
      PDM_log_trace_array_long(ordered_edge_opp, face_len, "ordered edges_opp");
      PDM_log_trace_array_long(ordered_vtx,     face_len, "ordered vtx");
      PDM_log_trace_array_long(ordered_vtx_opp, face_len, "ordered vtx_opp");
    }
    //Copy results for this face
    memcpy(&p_all_vtx[glob_idx],      ordered_vtx,      face_len*sizeof(PDM_g_num_t));
    memcpy(&p_all_vtx_opp[glob_idx],  ordered_vtx_opp,  face_len*sizeof(PDM_g_num_t));
    //memcpy(&p_all_edge[glob_idx],     ordered_edge      face_len*sizeof(PDM_g_num_t));
    //memcpy(&p_all_edge_opp[glob_idx], ordered_edge_opp, face_len*sizeof(PDM_g_num_t));
    glob_idx     += face_len;
  }
  free(ordered_edge);
  free(ordered_edge_opp);
  free(ordered_vtx);
  free(ordered_vtx_opp);
}

static void
_create_vtx_join
(
int            n_interface,
int            p_all_vtx_n,
PDM_g_num_t   *p_all_vtx,
PDM_g_num_t   *p_all_vtx_opp,
int           *p_all_vtx_group,
int           *p_all_vtx_dom_id,
int           *p_all_vtx_domopp_id,
int           *vtx_interface_size,
PDM_g_num_t  **interface_vtx_ids,
int          **interface_vtx_dom_ids,
PDM_MPI_Comm   comm
)
{
  //Todo : we could shift back to position of vtx in extraction to have a better
  //balance of edge distribution
  int *stride_one = PDM_array_const_int(p_all_vtx_n, 1);
  //Merge vertices using gnum
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &p_all_vtx,
                                                      NULL,
                                                     &p_all_vtx_n,
                                                      1,
                                                      comm);

  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *dall_vtx   = PDM_part_to_block_block_gnum_get(ptb);

  int         *recv_stride      = NULL;
  int         *dall_vtx_dom     = NULL;
  int         *dall_vtx_dom_opp = NULL;
  int         *dall_vtx_group   = NULL;
  PDM_g_num_t *dall_vtx_opp     = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_opp,
                                     &recv_stride,
                           (void **) &dall_vtx_opp);

  int *unused_recv_stride = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_group,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_group);
  free(unused_recv_stride);

  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_dom_id,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_dom);
  free(unused_recv_stride);

  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_domopp_id,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_dom_opp);
  free(unused_recv_stride);
  free(stride_one);

  if (0 == 1) {
    PDM_log_trace_array_long(dall_vtx,       blk_size,  "dall_vtx            :");
    PDM_log_trace_array_int (recv_stride,    blk_size,  "recv stride         :");
    PDM_log_trace_array_long(dall_vtx_opp,   exch_size, "recv dall_vtx_opp   :");
    PDM_log_trace_array_int (dall_vtx_group, exch_size, "recv dall_vtx_group :");
  }

  //First, dispatch vertices depending of the original interface
  PDM_array_reset_int(vtx_interface_size, n_interface, 0);
  int start_vtx = 0;
  for (int i_vtx =  0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    for (int i = 0; i < n_recv; i++) {
      vtx_interface_size[dall_vtx_group[start_vtx + i]]++;
    }
    start_vtx += n_recv;
  }
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    interface_vtx_ids[i_interface]     = (PDM_g_num_t *) malloc(2*vtx_interface_size[i_interface]*sizeof(PDM_g_num_t));
    interface_vtx_dom_ids[i_interface] = (int         *) malloc(2*vtx_interface_size[i_interface]*sizeof(int));
  }
  PDM_array_reset_int(vtx_interface_size, n_interface, 0);
  start_vtx = 0;
  for (int i_vtx =  0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    for (int i = 0; i < n_recv; i++) {
      int i_interface = dall_vtx_group[start_vtx + i];
      interface_vtx_ids[i_interface][2*vtx_interface_size[i_interface]]   = dall_vtx[i_vtx];
      interface_vtx_ids[i_interface][2*vtx_interface_size[i_interface]+1] = dall_vtx_opp[start_vtx+i];
      interface_vtx_dom_ids[i_interface][2*vtx_interface_size[i_interface]]   = dall_vtx_dom    [start_vtx+i];
      interface_vtx_dom_ids[i_interface][2*vtx_interface_size[i_interface]+1] = dall_vtx_dom_opp[start_vtx+i];
      vtx_interface_size[i_interface]++;
    }
    start_vtx += n_recv;
  }

  //Then, for each interface, eliminate pairs of vertices occuring more than once
  //If dom_id & dom_opp_id differs, vtx_id & opp should also differ because of the shift
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    
    int n_pairs_u = _unique_pairs(vtx_interface_size[i_interface],
                                  interface_vtx_ids[i_interface],
                                  interface_vtx_dom_ids[i_interface]);

    //Update
    vtx_interface_size[i_interface] = n_pairs_u;
    interface_vtx_ids[i_interface] = (PDM_g_num_t *) realloc(interface_vtx_ids[i_interface], 2*n_pairs_u*sizeof(PDM_g_num_t));
    interface_vtx_dom_ids[i_interface] = (int *) realloc(interface_vtx_dom_ids[i_interface], 2*n_pairs_u*sizeof(int));
  }

  PDM_part_to_block_free(ptb);
  free(recv_stride);
  free(dall_vtx_dom);
  free(dall_vtx_dom_opp);
  free(dall_vtx_opp);
  free(dall_vtx_group);
}


static void _connect_additional_edges
(
 int           n_extr_face,
 int          *face_vtx_both_idx,
 PDM_g_num_t  *face_vtx_both,
 PDM_g_num_t  *dface_edge,
 PDM_g_num_t  *dextract_face_join,
 PDM_g_num_t  *pedge_vtx,
 int          *face_status,
 PDM_g_num_t  *face_edge_wopp,
 PDM_MPI_Comm  comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  //TODO : equilibrate with weights
  //1. Use face_vtx connectivity to send to the vertex the id of face to which they belong
  PDM_part_to_block_t *ptb_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                         &face_vtx_both,
                                                          NULL,
                                                         &face_vtx_both_idx[n_extr_face],
                                                          1,
                                                          comm);

  PDM_g_num_t *send_data = (PDM_g_num_t *) malloc(2*face_vtx_both_idx[n_extr_face]*sizeof(PDM_g_num_t));
  for (int i_face = 0; i_face < n_extr_face / 2; i_face++) { //Split loop because of face,face_opp order
    for (int i_vtx = face_vtx_both_idx[2*i_face]; i_vtx < face_vtx_both_idx[2*i_face+1]; i_vtx++) {
      send_data[2*i_vtx] = dextract_face_join[2*i_face];
      send_data[2*i_vtx+1] = dextract_face_join[2*i_face+1];
    }
    for (int i_vtx = face_vtx_both_idx[2*i_face+1]; i_vtx < face_vtx_both_idx[2*(i_face+1)]; i_vtx++) {
      send_data[2*i_vtx] = dextract_face_join[2*i_face+1];
      send_data[2*i_vtx+1] = dextract_face_join[2*i_face];
    }
  }
  //PDM_log_trace_array_long(send_data, 2*face_vtx_both_idx[n_extr_face], "send data");
  int *stride_two = PDM_array_const_int(face_vtx_both_idx[n_extr_face], 2);

  int *vtx_stride_two = NULL;
  PDM_g_num_t *vtx_face_ids = NULL;
  int n_recv = PDM_part_to_block_exch(ptb_vtx,
                                      sizeof(PDM_g_num_t),
                                      PDM_STRIDE_VAR_INTERLACED,
                                      -1,
                            (int **)  &stride_two,
                           (void **)  &send_data,
                                      &vtx_stride_two,
                            (void **) &vtx_face_ids);

  int n_vtx_blk = PDM_part_to_block_n_elt_block_get(ptb_vtx);
  PDM_g_num_t *vtx_distri = (PDM_g_num_t *) malloc((n_rank+1)*sizeof(PDM_g_num_t));
  memcpy(vtx_distri, PDM_part_to_block_distrib_index_get(ptb_vtx), (n_rank+1)*sizeof(PDM_g_num_t));

  PDM_g_num_t *vtx_gnum = PDM_part_to_block_block_gnum_get(ptb_vtx);
  
  if (0 == 1) {
    PDM_log_trace_array_long(vtx_gnum, n_vtx_blk, "block gnum");
    PDM_log_trace_array_int(vtx_stride_two, n_vtx_blk, "recv stride 2");
    PDM_log_trace_array_long(vtx_face_ids, n_recv, "recv data");
  }

  int *vtx_face_n = (int *) malloc(n_vtx_blk * sizeof(int)); //We received 2 data per face connected to each vtx
  for (int i=0; i < n_vtx_blk; i++) {
    vtx_face_n[i] = vtx_stride_two[i] / 2;
  }

  free(stride_two);
  free(send_data);

  // 2. Build key and send face ids & gnum using key numbering
  PDM_g_num_t *vtx_key = PDM_array_const_gnum(n_vtx_blk, 0);
  int read_idx = 0;
  for (int j=0; j < n_vtx_blk; j++) {
    for (int k = 0; k < vtx_stride_two[j]; k++)
      vtx_key[j] += vtx_face_ids[read_idx++];
  }

  // PDM_log_trace_array_long(vtx_key, n_vtx_blk, "vtx_key");

  PDM_part_to_block_t *ptb_key = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                        (PDM_g_num_t **) &vtx_key,
                                                          NULL,
                                                         &n_vtx_blk,
                                                          1,
                                                          comm);

  int *stride_one = PDM_array_const_int(n_vtx_blk, 1);

  int         *unused_recv_stride = NULL;
  PDM_g_num_t *key_vtx_gnum       = NULL;
  // n_key_vtx is the number of vertices involved in keys (counted multiple times) managed by this proc
  int n_key_vtx =  PDM_part_to_block_exch(ptb_key,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int **)  &stride_one,
                           (void **) &vtx_gnum,
                                     &unused_recv_stride, //Same as key count
                           (void **) &key_vtx_gnum);
  free(unused_recv_stride);
  int *key_recv_face_n = NULL; // For each key (unmerged), number of face related to key
  PDM_part_to_block_exch(ptb_key,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
               (int **)  &stride_one,
               (void **) &vtx_face_n,
                         &unused_recv_stride,  //Same as key count
               (void **) &key_recv_face_n);
  free(unused_recv_stride);

  int         *key_recv_stride = NULL;
  PDM_g_num_t *key_recv_data = NULL; //For each key (unmerged), tuples (face/face_opp)  * nb of face related to key
  PDM_part_to_block_exch(ptb_key,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
               (int **)  &vtx_stride_two,
               (void **) &vtx_face_ids,
                         &key_recv_stride,
               (void **) &key_recv_data);


  free(stride_one);
  PDM_part_to_block_free(ptb_vtx); // Needed until here for vtx gnum

  int n_keys = PDM_part_to_block_n_elt_block_get(ptb_key);
  PDM_g_num_t *keys_ids = PDM_part_to_block_block_gnum_get(ptb_key);
  int         *keys_cnt = PDM_part_to_block_block_gnum_count_get(ptb_key);

  if (0 == 1) {
    PDM_log_trace_array_long(keys_ids, n_keys, "key to treat");
    PDM_log_trace_array_int(keys_cnt, n_keys, "n recept");
    PDM_log_trace_array_int(key_recv_stride, n_keys, "key recv stride");
    PDM_log_trace_array_int(key_recv_face_n, n_key_vtx, "key recv facen"); 
    //PDM_log_trace_array_long(key_recv_data, n_recv, "key recv data"); 
    PDM_log_trace_array_long(key_vtx_gnum, n_key_vtx, "key recv gnum");
  }

  free(vtx_face_ids);
  free(vtx_face_n);
  free(vtx_stride_two);
  free(vtx_key);
  
  //3. Match data on key distribution
  PDM_g_num_t *key_vtx_gnum_opp = PDM_array_const_gnum(n_key_vtx, 0);
  int count_idx = 0; //Start of data in key_recv_face_n
  int data_idx = 0; //Start of data in key_recv_data
  for (int i_key = 0; i_key < n_keys; i_key++) {
    int n_vtx_this_key = keys_cnt[i_key];
    // Each key occurs n_vtx_this_key times, each time with 2 face ids * nb of face connected to vertex
    // First we sort these sections
    for (int k=0; k < n_vtx_this_key; k++) {
      PDM_sort_long(&key_recv_data[data_idx], NULL, 2*key_recv_face_n[count_idx]);
      data_idx += 2*key_recv_face_n[count_idx];
      count_idx++;
    }
    // Now search matches
    int idx1 = data_idx - key_recv_stride[i_key]; //Reset idx1 for comparaison
    count_idx -= n_vtx_this_key;
    for (int k = 0; k < n_vtx_this_key; k++) {
      int n_match = 0;
      int i_match;
      int idx2 = data_idx - key_recv_stride[i_key]; //Reset idx2 for comparaison

      for (int k2 = 0; k2 < n_vtx_this_key; k2++) {
        if (k2 != k) { //Skip myself
          if (key_recv_face_n[count_idx+k] == key_recv_face_n[count_idx+k2]) {
            if (PDM_array_are_equal_gnum(&key_recv_data[idx1], &key_recv_data[idx2], 2*key_recv_face_n[count_idx+k])) {
              n_match++;
              i_match = k2;
            }
          }
        }
        idx2 += 2*key_recv_face_n[count_idx+k2];
      }
      if (n_match == 1) { //Register match
        key_vtx_gnum_opp[count_idx + k] = key_vtx_gnum[count_idx+i_match];
      }
      idx1 += 2*key_recv_face_n[count_idx+k];
    }
    count_idx += n_vtx_this_key;
  }

  PDM_part_to_block_free(ptb_key);
  free(key_recv_stride);
  free(key_recv_data);
  free(key_recv_face_n);

  // 4. We send back the matches to vertex distribution to have block property
  PDM_part_to_block_t *ptb_vtx2 = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                            PDM_PART_TO_BLOCK_POST_NOTHING,
                                                            1.,
                                          (PDM_g_num_t **) &key_vtx_gnum,
                                                            vtx_distri,
                                                           &n_key_vtx,
                                                            1,
                                                            comm);
  assert (PDM_part_to_block_n_elt_block_get(ptb_vtx2) == n_vtx_blk);
  PDM_g_num_t *matched_gnum = (PDM_g_num_t *) malloc(n_vtx_blk*sizeof(PDM_g_num_t));
  memcpy(matched_gnum, PDM_part_to_block_block_gnum_get(ptb_vtx2), n_vtx_blk*sizeof(PDM_g_num_t));

  PDM_g_num_t *matched_gnum_opp = NULL;
  PDM_part_to_block_exch(ptb_vtx2,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &key_vtx_gnum_opp,
                         NULL,
               (void **) &matched_gnum_opp);

  PDM_part_to_block_free(ptb_vtx2);
  free(key_vtx_gnum);
  free(key_vtx_gnum_opp);


  // 5a. Prepare edge matching : get the vtx_gnum_opp only for the untreated faces
  int unsolvable_edge = 0;
  for (int i_face=0; i_face < n_extr_face; i_face++) {
    if (face_status[i_face] == 0)
      unsolvable_edge += face_vtx_both_idx[i_face+1] - face_vtx_both_idx[i_face];
  }

  PDM_g_num_t *requested_gnum = (PDM_g_num_t *) malloc(unsolvable_edge*sizeof(PDM_g_num_t));
  int idx = 0;
  for (int i_face = 0; i_face < n_extr_face; i_face++) {
    if (face_status[i_face] == 0) {
      int n_vtx_face = face_vtx_both_idx[i_face+1] - face_vtx_both_idx[i_face];
      memcpy(&requested_gnum[idx], &face_vtx_both[face_vtx_both_idx[i_face]], n_vtx_face*sizeof(PDM_g_num_t));
      idx += n_vtx_face;
    }
  }
  PDM_block_to_part_t *btp = PDM_block_to_part_create(vtx_distri,
                               (const PDM_g_num_t **) &requested_gnum,
                                                      &unsolvable_edge,
                                                      1,
                                                      comm);
  int *blk_stride = PDM_array_zeros_int(vtx_distri[i_rank+1] - vtx_distri[i_rank]);
  for (int i = 0; i < n_vtx_blk; i++) {
    blk_stride[matched_gnum[i] - vtx_distri[i_rank] - 1] = 1;
  }
  int         **recv_stride;
  PDM_g_num_t **recv_data;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          blk_stride,
                          matched_gnum_opp, 
                         &recv_stride,
              (void ***) &recv_data);
  free(recv_stride[0]);
  free(recv_stride);
  PDM_g_num_t *requested_gnum_opp = recv_data[0];
  free(recv_data);
  free(blk_stride);
  PDM_block_to_part_free(btp);

  free(matched_gnum_opp);
  free(matched_gnum);
  free(vtx_distri);

  //5b. Now we can match edges and update face_edge_wopp !
  int extracted_face_id = 0;
  for (int i_face = 0; i_face < n_extr_face/2; i_face++) {
    if (face_status[2*i_face] == 0) {
      int edge_start     = face_vtx_both_idx[2*i_face];
      int edge_opp_start = face_vtx_both_idx[2*i_face+1];
      int n_vtx_face = face_vtx_both_idx[2*i_face+1] - face_vtx_both_idx[2*i_face];
      
      //Find any opp vtx gnum to have a gnum / gnum opp couple
      PDM_g_num_t opp_vtx_gnum = 0;
      int         opp_vtx_pos  = 0;
      for (int j=0; j < n_vtx_face; j++) {
        opp_vtx_gnum = requested_gnum_opp[extracted_face_id + j];
        if (opp_vtx_gnum != 0) {
          opp_vtx_pos = j;
          break;
        }
      }
      assert (opp_vtx_gnum != 0);
      PDM_g_num_t my_vtx_gnum = requested_gnum[extracted_face_id+opp_vtx_pos];

      // Search my vtx gnum in edges
      int edge_pos = -1;
      int edge_sens;
      for (int j = 0; j < n_vtx_face; j++) {
        if (pedge_vtx[2*edge_start + 2*j] == my_vtx_gnum) {
          edge_pos = j;
          edge_sens = 0;
          break;
        }
        else if (pedge_vtx[2*edge_start + 2*j + 1] == my_vtx_gnum) {
          edge_pos = j;
          edge_sens = 1;
          break;
        }
      }
      assert (edge_pos >= 0);
      
      //Search opp vtx in edge opp, reverse order
      int opp_pos = -1;
      int opp_sens = 1 - edge_sens;
      for (int j = 0; j < n_vtx_face; j++) {
        if (pedge_vtx[2*edge_opp_start + 2*j + opp_sens] == opp_vtx_gnum) {
          opp_pos = j;
          break;
        }
      }
      assert (opp_pos >= 0);

      //Complete face_edge_wopp for face & opp at the same time
      face_edge_wopp[edge_start + edge_pos] = dface_edge[edge_opp_start + opp_pos];
      face_edge_wopp[edge_opp_start + opp_pos] = dface_edge[edge_start + edge_pos];

      face_status[2*i_face]   = 1;
      face_status[2*i_face+1] = 1;
      extracted_face_id += 2*n_vtx_face;
    }
  }

  free(requested_gnum);
  free(requested_gnum_opp);
}



static void _domain_interface_face_to_vertex
(
 int            n_interface,             /* Total number of interfaces */
 int           *interfaces_size,         /* Number of face pairs in each interface */
 PDM_g_num_t  **interface_face_ids,      /* For each interface, list of pairs face,face_opp */
 int          **interface_domains_ids,   /* For each interface, list of domains dom,dom_opp */
 int            n_domain,                  /* Number of domains */
 int           *dn_vtx,                  /* Number of vertex in each domain (distributed) */
 int           *dn_face,                 /* Number of face in each domain (distributed) */
 int          **dface_vtx_idx,           /* Face->vertex connectivity for each domain */
 PDM_g_num_t  **dface_vtx,
 int           *vtx_interface_size,      /* [OUT] Number of vtx pairs in each interface */
 PDM_g_num_t  **interface_vtx_ids,       /* [OUT] For each interface, list of pairs vtx,vtx_opp */
 int          **interface_vtx_dom_ids,   /* [OUT] For each interface, list of domains dom,dom_opp */
 PDM_MPI_Comm   comm                     /* Mpi comunicator */
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *face_per_block_offset = _per_block_offset(n_domain, dn_face, comm);
  PDM_g_num_t *vtx_per_block_offset  = _per_block_offset(n_domain, dn_vtx,  comm);

  int         *face_vtx_both_idx = NULL;
  PDM_g_num_t *face_vtx_both     = NULL;
  int         *dextract_face_dom_id   = NULL;
  int         *dextract_face_group_id = NULL;
  PDM_g_num_t *dextract_face_join     = NULL;
  int n_extr_face = _extract_and_shift_jn_faces(n_domain,
                                                dn_face,
                                                dn_vtx,
                                                n_interface,
                                                interfaces_size,
                                                interface_face_ids,
                                                interface_domains_ids,
                                                dface_vtx_idx,
                                                dface_vtx,
                                               &face_vtx_both_idx,
                                               &face_vtx_both,
                                               &dextract_face_group_id,
                                               &dextract_face_dom_id,
                                               &dextract_face_join,
                                                comm);

  PDM_g_num_t *extracted_face_distri = PDM_compute_entity_distribution(comm, n_extr_face);

  //Duplicate this data for easier send to edges
  PDM_g_num_t *dextract_face_join_opp   = (PDM_g_num_t *) malloc(n_extr_face*sizeof(PDM_g_num_t));
  for (int i = 0; i < n_extr_face/2; i++) {
    dextract_face_join_opp[2*i]     = dextract_face_join[2*i+1];
    dextract_face_join_opp[2*i+1]   = dextract_face_join[2*i];
  }

  if (0 == 1) {
    log_trace("Face vtx received after MBTP\n");
    PDM_log_trace_array_int (face_vtx_both_idx, n_extr_face+1, "face_vtx_idx :: ");
    PDM_log_trace_array_long(face_vtx_both, face_vtx_both_idx[n_extr_face], "face_vtx :: ");
  }


  //Now we have some pairs of faces (each pair appears only one) + face_vtx for this pairs

  //Generate edge numbering

  PDM_g_num_t *dedge_distrib  = NULL;
  int         *dedge_vtx_idx  = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;
  int         *dedge_face_idx = NULL;
  PDM_g_num_t *dedge_face     = NULL;
  int dn_edge = _generate_edge_face(n_extr_face,
                                    face_vtx_both_idx,
                                    face_vtx_both,
                                   &dedge_distrib,
                                   &dedge_vtx_idx,
                                   &dedge_vtx,
                                   &dedge_face_idx,
                                   &dedge_face,
                                    comm);
  
  if (0 == 1) {
    log_trace("Edges rebuild\n");
    PDM_log_trace_array_long(dedge_distrib, n_rank+1, "dedge_distri ::");
    //PDM_log_trace_array_long(dedge_vtx, 2*dn_edge, "dedge_vtx ::");
    //PDM_log_trace_connectivity_long(dedge_face_idx, dedge_face, dn_edge, "dedge_face ::");
  }

  // Transport face data to edges

  //Prepare numbering
  PDM_g_num_t *dedge_face_abs = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  int         *dedge_face_sgn = (int         *) malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  for(int i = 0; i < dedge_face_idx[dn_edge]; ++i) {
    dedge_face_abs[i] = PDM_ABS (dedge_face[i]);
    dedge_face_sgn[i] = PDM_SIGN(dedge_face[i]);
  }


  PDM_g_num_t *dedge_face_join         = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_g_num_t *dedge_face_join_opp     = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_block_to_part_t *btp = PDM_block_to_part_create(extracted_face_distri,
                               (const PDM_g_num_t **) &dedge_face_abs,
                                                      &dedge_face_idx[dn_edge],
                                                      1,
                                                      comm);
  int cst_stride = 1;
  PDM_block_to_part_exch_in_place(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &cst_stride,
                (void *) dextract_face_join,
                         NULL,
             (void ** ) &dedge_face_join);

  PDM_block_to_part_exch_in_place(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &cst_stride,
                (void *) dextract_face_join_opp,
                         NULL,
             (void ** ) &dedge_face_join_opp);

  PDM_block_to_part_free(btp);
  free(dedge_face_abs);
  free(dedge_face_sgn);
  free(dextract_face_join_opp);


  if (0 == 1) {
    log_trace("Transport data on edges\n");
    PDM_log_trace_array_int (dedge_face_idx,        dn_edge+1,               "dedge_face_idx        ::");
    PDM_log_trace_array_long(dedge_face_join    ,   dedge_face_idx[dn_edge], "dedge_face_join       ::");
    PDM_log_trace_array_long(dedge_face_join_opp,   dedge_face_idx[dn_edge], "dedge_face_join_opp   ::");
  }

  // Match internal edges
  PDM_g_num_t *dedge_gnum     = NULL;
  PDM_g_num_t *dedge_gnum_opp = NULL;
  int dn_internal_edge = _match_internal_edges(dn_edge,
                                               dedge_distrib,
                                               dedge_face_idx,
                                               dedge_face,
                                               dedge_face_join,
                                               dedge_face_join_opp,
                                              &dedge_gnum,
                                              &dedge_gnum_opp,
                                               comm);

  if (0 == 1) {
    log_trace("Internal edge matches after conflict resolution \n");
    PDM_log_trace_array_long(dedge_gnum, dn_internal_edge, "dedge gnum :: ");
    PDM_log_trace_array_long(dedge_gnum_opp, dn_internal_edge, "dedge gnum_opp :: ");
  }

  // To match external edges, we will work on face distribution. We need face->edge connectivity

  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  PDM_dconnectivity_transpose(comm,
                              dedge_distrib,
                              extracted_face_distri,
                              dedge_face_idx,
                              dedge_face,
                              1,
                             &dface_edge_idx,
                             &dface_edge);

  if (0 == 1) {
    log_trace("Generate dface->edge from dedge->face\n");
    PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, n_extr_face, "dface_edge :: ");
  }

  //Transfert some data from the edges to the edges know by the faces

  //Prepare gnum and stride
  PDM_g_num_t *dface_edge_abs = (PDM_g_num_t *) malloc(dface_edge_idx[n_extr_face] * sizeof(PDM_g_num_t));
  for (int i = 0; i < dface_edge_idx[n_extr_face]; i++) {
    dface_edge_abs[i] = PDM_ABS(dface_edge[i]);
  }
  int idx = 0;
  int *dedge_gnum_n = PDM_array_zeros_int(dn_edge);
  for (int i = 0; i < dn_edge; i++) {
    if (i + dedge_distrib[i_rank] + 1 == dedge_gnum[idx]) {
      dedge_gnum_n[i] = 1;
      idx++;
    }
    if (idx == dn_internal_edge) { //End of internal edge reached, no more comparaison is needed
      break;
    }
  }

                       btp = PDM_block_to_part_create(dedge_distrib,
                               (const PDM_g_num_t **) &dface_edge_abs,
                                                      &dface_edge_idx[n_extr_face],
                                                      1,
                                                      comm);


  int         **recv_stride_tmp = NULL;
  PDM_g_num_t **recv_data_tmp   = NULL;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          dedge_gnum_n,
                          dedge_gnum_opp,
                         &recv_stride_tmp,
              (void ***) &recv_data_tmp);
  int         *pedge_gnum_n   = recv_stride_tmp[0];
  PDM_g_num_t *pedge_gnum_opp = recv_data_tmp[0];
  free(recv_stride_tmp);
  free(recv_data_tmp);

  int stride2 = 2;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                         &stride2,
                          dedge_vtx,
                          NULL,
              (void ***) &recv_data_tmp);
  PDM_g_num_t *pedge_vtx = recv_data_tmp[0];
  free(recv_data_tmp);
  //PDM_log_trace_array_long(pedge_vtx, 2*dface_edge_idx[n_extr_face], "pedge_vtx");



  //Attention, on devrait pouvoir travailler sur face externes uniquement (filtre le dface_edge_abs)
  PDM_g_num_t *face_edge_wopp = (PDM_g_num_t *) malloc(dface_edge_idx[n_extr_face]*sizeof(PDM_g_num_t));
  idx = 0;
  for (int i = 0; i < dface_edge_idx[n_extr_face]; i++) {
    if (pedge_gnum_n[i] == 1)
      face_edge_wopp[i] = pedge_gnum_opp[idx++];
    else
      face_edge_wopp[i] = 0;
  }

  //PDM_log_trace_connectivity_long(dface_edge_idx, face_edge_wopp, n_extr_face, "dface_edge :: ");
  int *face_status = (int *) malloc(n_extr_face*sizeof(int));
  for (int i_face=0; i_face < n_extr_face; i_face++) {
    int n_treated_edge = 0;
    for (int i_edge = dface_edge_idx[i_face]; i_edge < dface_edge_idx[i_face+1]; i_edge++)
      n_treated_edge += pedge_gnum_n[i_edge];
    face_status[i_face] = (int) (n_treated_edge > 0);
  }


  free(dface_edge_abs);
  free(dedge_gnum_n);
  free(pedge_gnum_n);
  free(pedge_gnum_opp);
  PDM_block_to_part_free(btp);

  int need_more_edge_l = 0;
  int need_more_edge;
  for (int i_face=0; i_face < n_extr_face; i_face++) {
    if (face_status[i_face] == 0) {
      need_more_edge_l = 1; 
      break;
    }
  }
  PDM_MPI_Allreduce(&need_more_edge_l, &need_more_edge, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  if (need_more_edge > 0) {
    log_trace("Warning -- Some face have not shared edges. Try to retrieve it using vertex connectivity\n");
    _connect_additional_edges(n_extr_face,
                              face_vtx_both_idx,
                              face_vtx_both,
                              dface_edge,
                              dextract_face_join,
                              pedge_vtx,
                              face_status,
                              face_edge_wopp,
                              comm);
  }
  free(face_status);

  //Match external edges
  assert (dface_edge_idx[n_extr_face] % 2 == 0);
  int n_vtx_interface_tot = dface_edge_idx[n_extr_face] / 2;
  PDM_g_num_t *p_all_vtx      = (PDM_g_num_t *) malloc(n_vtx_interface_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *p_all_vtx_opp  = (PDM_g_num_t *) malloc(n_vtx_interface_tot * sizeof(PDM_g_num_t));
  /*PDM_g_num_t *p_all_edge_gnum     = malloc(dface_edge_idx[n_extr_face] * sizeof(PDM_g_num_t));*/
  /*PDM_g_num_t *p_all_edge_gnum_opp = malloc(dface_edge_idx[n_extr_face] * sizeof(PDM_g_num_t));*/
  _match_all_edges_from_faces(n_extr_face,
                              dface_edge_idx,
                              dface_edge,
                              face_edge_wopp,
                              pedge_vtx,
                              p_all_vtx,
                              p_all_vtx_opp);

  //Copy group from face to vertices
  int *p_all_vtx_group     = (int *) malloc(n_vtx_interface_tot * sizeof(int));
  int *p_all_vtx_dom_id    = (int *) malloc(n_vtx_interface_tot * sizeof(int));
  int *p_all_vtx_domopp_id = (int *) malloc(n_vtx_interface_tot * sizeof(int));
  int glob_idx = 0;
  for (int i_face = 0; i_face < n_extr_face/2; i_face++) {
    int face_len = dface_edge_idx[2*i_face+1] - dface_edge_idx[2*i_face];
    for (int i = 0; i < face_len; i++) {
      p_all_vtx_group[glob_idx]     = dextract_face_group_id[2*i_face];
      p_all_vtx_dom_id[glob_idx]    = dextract_face_dom_id  [2*i_face];
      p_all_vtx_domopp_id[glob_idx] = dextract_face_dom_id  [2*i_face+1];
      glob_idx ++;
    }
  }


  free(face_edge_wopp);
  if (0 == 1) {
    log_trace("Vtx matching on face distribution\n");
    PDM_log_trace_array_long(p_all_vtx,     n_vtx_interface_tot, "p_all_vtx     ::");
    PDM_log_trace_array_long(p_all_vtx_opp, n_vtx_interface_tot, "p_all_vtx_opp ::");
    PDM_log_trace_array_int (p_all_vtx_group, n_vtx_interface_tot, "p_all_vtx_group ::");
    PDM_log_trace_array_int (p_all_vtx_dom_id, n_vtx_interface_tot, "p_all_vtx_dom_id ::");
    PDM_log_trace_array_int (p_all_vtx_domopp_id, n_vtx_interface_tot, "p_all_vtx_domopp_id ::");
  }

  _create_vtx_join(n_interface,
                   n_vtx_interface_tot,
                   p_all_vtx,
                   p_all_vtx_opp,
                   p_all_vtx_group,
                   p_all_vtx_dom_id,
                   p_all_vtx_domopp_id,
                   vtx_interface_size,
                   interface_vtx_ids,
                   interface_vtx_dom_ids,
                   comm);


  //Ultimate step : go back to original vtx numbering. All we have to do is retrieve domain
  // and substract domain offset
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    for (int i_vtx = 0; i_vtx < 2*vtx_interface_size[i_interface]; i_vtx++) {
      interface_vtx_ids[i_interface][i_vtx] -= vtx_per_block_offset[interface_vtx_dom_ids[i_interface][i_vtx]];
    }
  }
  if (0 == 1) {
    PDM_log_trace_array_int(vtx_interface_size, n_interface, "Vtx interfaces sizes");
    for (int i = 0; i < n_interface; i++) {
      log_trace("Vtx & vtx opp in vtx global numbering for interface %d \n", i);
      PDM_log_trace_array_long(interface_vtx_ids[i], 2*vtx_interface_size[i], "vertex ids     ::");
      PDM_log_trace_array_int(interface_vtx_dom_ids[i], 2*vtx_interface_size[i], "vertex doms     ::");
    }
  }



  free(dface_edge_idx);
  free(dface_edge);

  free(dextract_face_join);
  free(dextract_face_group_id);
  free(dextract_face_dom_id);

  free(dedge_face_join);
  free(dedge_face_join_opp);


  free(extracted_face_distri);

  free(face_vtx_both_idx);
  free(face_vtx_both);

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);

  free(dedge_gnum);
  free(dedge_gnum_opp);

  free(p_all_vtx);
  free(p_all_vtx_opp);
  free(p_all_vtx_group);
  free(p_all_vtx_dom_id);
  free(p_all_vtx_domopp_id);

  free(pedge_vtx);
  free(face_per_block_offset);
  free(vtx_per_block_offset);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

PDM_domain_interface_t*
PDM_domain_interface_create
(
const int                         n_interface,
const int                         n_domain,
      PDM_domain_interface_mult_t multidomain_interface,
      PDM_ownership_t             ownership,
      PDM_MPI_Comm                comm
)
{
  PDM_domain_interface_t *dom_intrf = (PDM_domain_interface_t *) malloc (sizeof(PDM_domain_interface_t));
  dom_intrf->n_interface       = n_interface;
  dom_intrf->n_domain          = n_domain;
  dom_intrf->multidomain_intrf = multidomain_interface;
  dom_intrf->ownership         = ownership;
  dom_intrf->comm              = comm;

  dom_intrf->interface_dn_face   = NULL;
  dom_intrf->interface_ids_face  = NULL;
  dom_intrf->interface_dom_face  = NULL;
  dom_intrf->interface_dn_edge   = NULL;
  dom_intrf->interface_ids_edge  = NULL;
  dom_intrf->interface_dom_edge  = NULL;
  dom_intrf->interface_dn_vtx    = NULL;
  dom_intrf->interface_ids_vtx   = NULL;
  dom_intrf->interface_dom_vtx   = NULL;

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    dom_intrf->is_result[i] = 0;
  }

  dom_intrf->translation_vect   = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_direction = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_center    = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_angle     = (double  *) malloc(n_interface * sizeof(double  ));

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    dom_intrf->translation_vect  [i_interface] = NULL;
    dom_intrf->rotation_direction[i_interface] = NULL;
    dom_intrf->rotation_center   [i_interface] = NULL;
    dom_intrf->rotation_angle    [i_interface] = 0.;
  }

  return dom_intrf;
}

void PDM_domain_interface_set
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                    *interface_dn,
 PDM_g_num_t           **interface_ids,
 int                   **interface_dom
)
{
  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    dom_intrf->interface_dn_vtx  = interface_dn;
    dom_intrf->interface_ids_vtx = interface_ids;
    dom_intrf->interface_dom_vtx = interface_dom;
  } else if (interface_kind == PDM_BOUND_TYPE_EDGE) {
    dom_intrf->interface_dn_edge  = interface_dn;
    dom_intrf->interface_ids_edge = interface_ids;
    dom_intrf->interface_dom_edge = interface_dom;
  } else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    dom_intrf->interface_dn_face  = interface_dn;
    dom_intrf->interface_ids_face = interface_ids;
    dom_intrf->interface_dom_face = interface_dom;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
  }

}

void
PDM_domain_interface_translate_face2vtx
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *dn_vtx,
 int                     *dn_face,
 int                    **dface_vtx_idx,
 PDM_g_num_t            **dface_vtx
)
{
  assert (dom_intrf != NULL);
  assert (dom_intrf->interface_dn_face != NULL);
  assert (dom_intrf->interface_dn_vtx  == NULL);
  dom_intrf->interface_dn_vtx  = (int *)          malloc(dom_intrf->n_interface * sizeof(int));
  dom_intrf->interface_ids_vtx = (PDM_g_num_t **) malloc(dom_intrf->n_interface * sizeof(PDM_g_num_t*));
  dom_intrf->interface_dom_vtx = (int         **) malloc(dom_intrf->n_interface * sizeof(int*));

  // Simple case is not yet managed, copy to go back to full case
  int **_interface_dom_face = NULL;
  if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
    _interface_dom_face = (int **) malloc(dom_intrf->n_interface*sizeof(int*));
    for (int i_intrf = 0; i_intrf < dom_intrf->n_interface; i_intrf++) {
      _interface_dom_face[i_intrf] = (int *) malloc(2*dom_intrf->interface_dn_face[i_intrf]*sizeof(int));
      for (int j = 0; j < dom_intrf->interface_dn_face[i_intrf]; j++) {
        _interface_dom_face[i_intrf][2*j]   = dom_intrf->interface_dom_face[i_intrf][0];
        _interface_dom_face[i_intrf][2*j+1] = dom_intrf->interface_dom_face[i_intrf][1];
      }
    }
  }
  else {
    _interface_dom_face = dom_intrf->interface_dom_face;
  }

  _domain_interface_face_to_vertex(dom_intrf->n_interface,
                                   dom_intrf->interface_dn_face,
                                   dom_intrf->interface_ids_face,
                                   _interface_dom_face,
                                   dom_intrf->n_domain,
                                   dn_vtx,
                                   dn_face,
                                   dface_vtx_idx,
                                   dface_vtx,
                                   dom_intrf->interface_dn_vtx,
                                   dom_intrf->interface_ids_vtx,
                                   dom_intrf->interface_dom_vtx,
                                   dom_intrf->comm);

  dom_intrf->is_result[PDM_BOUND_TYPE_VTX] = 1;

  // Simple case is not yet managed, free working arrays
  if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
    for (int i_intrf = 0; i_intrf < dom_intrf->n_interface; i_intrf++) {
      for (int j = 0; j < dom_intrf->interface_dn_vtx[i_intrf]; j++) {
        assert(dom_intrf->interface_dom_vtx[i_intrf][2*j]   == dom_intrf->interface_dom_face[i_intrf][0]);
        assert(dom_intrf->interface_dom_vtx[i_intrf][2*j+1] == dom_intrf->interface_dom_face[i_intrf][1]);
      }
      free(dom_intrf->interface_dom_vtx[i_intrf]);
      free(_interface_dom_face[i_intrf]);
    }
    free(_interface_dom_face);
    free(dom_intrf->interface_dom_vtx);
    dom_intrf->interface_dom_vtx = dom_intrf->interface_dom_face;
  }
}


void
PDM_domain_interface_translate_entity1_entity2
(
 int                      n_domain,
 int                      n_interface,
 int                     *dn_entity1,
 int                     *dn_entity2,
 int                     *dn_interface,
 int                    **interface_dom,
 PDM_g_num_t            **interface_ids,
 int                    **dentity2_entity1_idx,
 PDM_g_num_t            **dentity2_entity1,
 int                      connectivity_is_signed,
 PDM_MPI_Comm             comm,
 int                    **interface_dn_entity2,
 PDM_g_num_t           ***interface_ids_entity2,
 int                   ***interface_dom_entity2
)
{
  // log_trace("PDM_domain_interface_translate_entity1_entity2 beg \n");
  // TODO :
  //  -> reduce time by extracting dentity2_entity1_idx for only concerns interfaces
  //  -> Pour l'insant la reduction est faite en dehors (via PDM_dmesh_extract )

  if(0 == 1) {
    for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
      PDM_log_trace_connectivity_long(dentity2_entity1_idx[i_domain],
                                      dentity2_entity1[i_domain],dn_entity2[i_domain], "dentity2_entity1 ::" );
    }
  }

  int i_rank = -1;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *entity1_per_block_offset = _per_block_offset(n_domain, dn_entity1, comm);
  PDM_g_num_t *entity2_per_block_offset = _per_block_offset(n_domain, dn_entity2, comm);

  // Prepare first PtB with multiple partitions.
  // Use (shifted) ids as gnum and send tuple (shited) id, opp_id
  PDM_g_num_t **interface_ids_shifted = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
  PDM_g_num_t **send_data_gnum        = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
  int         **send_data_dom         = (int         **) malloc(n_interface * sizeof(int         *));
  int         **send_data_sens        = (int         **) malloc(n_interface * sizeof(int         *));
  int         **send_data_intno       = (int         **) malloc(n_interface * sizeof(int         *));
  double      **weight                = (double      **) malloc(n_interface * sizeof(double      *));
  int         **stride_one            = (int         **) malloc(n_interface * sizeof(int         *));
  int          *dn_interface_twice    = (int          *) malloc(n_interface * sizeof(int          ));

  if(n_domain > 1) {
    for(int i_domain = 1; i_domain < n_domain; ++i_domain) {
      for(int i = 0; i < dentity2_entity1_idx[i_domain][dn_entity2[i_domain]]; ++i) {
        int sgn = PDM_SIGN(dentity2_entity1[i_domain][i]);
        dentity2_entity1[i_domain][i] = sgn * ( PDM_ABS(dentity2_entity1[i_domain][i]) + entity1_per_block_offset[i_domain]);
      }
    }
  }

  for (int itrf = 0; itrf < n_interface; itrf++) {
    stride_one           [itrf] = (int         *) malloc( 2 * dn_interface[itrf] * sizeof(int        ));
    interface_ids_shifted[itrf] = (PDM_g_num_t *) malloc( 2 * dn_interface[itrf] * sizeof(PDM_g_num_t));
    send_data_gnum       [itrf] = (PDM_g_num_t *) malloc( 2 * dn_interface[itrf] * sizeof(PDM_g_num_t));
    send_data_dom        [itrf] = (int         *) malloc( 2 * dn_interface[itrf] * sizeof(int        ));
    send_data_sens       [itrf] = (int         *) malloc( 2 * dn_interface[itrf] * sizeof(int        ));
    send_data_intno      [itrf] = (int         *) malloc( 2 * dn_interface[itrf] * sizeof(int        ));
    weight               [itrf] = (double      *) malloc( 2 * dn_interface[itrf] * sizeof(double     ));
    dn_interface_twice   [itrf] = 2*dn_interface[itrf];

    if (0 == 1) {
      PDM_log_trace_array_long(interface_ids[itrf], 2 * dn_interface[itrf], "interface_ids:: ");
    }

    for (int k = 0; k < dn_interface[itrf]; k++) {
      int dom    = interface_dom[itrf][2*k  ];
      int domopp = interface_dom[itrf][2*k+1];
      // interface_ids_shifted[itrf][2*k  ] = interface_ids[itrf][2*k  ] + entity1_per_block_offset[dom   ];
      // interface_ids_shifted[itrf][2*k+1] = interface_ids[itrf][2*k+1] + entity1_per_block_offset[domopp];
      // send_data_gnum       [itrf][2*k  ] = interface_ids[itrf][2*k+1] + entity1_per_block_offset[domopp];
      // send_data_gnum       [itrf][2*k+1] = interface_ids[itrf][2*k  ] + entity1_per_block_offset[dom   ];
      interface_ids_shifted[itrf][2*k  ] = PDM_ABS(interface_ids[itrf][2*k  ]) + entity1_per_block_offset[dom   ];
      interface_ids_shifted[itrf][2*k+1] = PDM_ABS(interface_ids[itrf][2*k+1]) + entity1_per_block_offset[domopp];

      int sgn_cur = PDM_SIGN(interface_ids[itrf][2*k  ]);
      int sgn_opp = PDM_SIGN(interface_ids[itrf][2*k+1]);
      send_data_gnum       [itrf][2*k  ] = PDM_ABS(interface_ids[itrf][2*k+1]) + entity1_per_block_offset[domopp];
      send_data_gnum       [itrf][2*k+1] = PDM_ABS(interface_ids[itrf][2*k  ]) + entity1_per_block_offset[dom   ];
      send_data_dom        [itrf][2*k  ] = domopp;
      send_data_dom        [itrf][2*k+1] = dom   ;

      /* Send it separatly */
      // send_data_sens       [itrf][2*k  ] = 1;
      // send_data_sens       [itrf][2*k+1] = sgn_opp;

      // send_data_sens       [itrf][2*k  ] = sgn_opp;
      // send_data_sens       [itrf][2*k+1] = 1;

      // send_data_sens       [itrf][2*k  ] = sgn_cur;
      // send_data_sens       [itrf][2*k+1] = sgn_opp;
      // send_data_sens       [itrf][2*k  ] = sgn_opp;
      // send_data_sens       [itrf][2*k+1] = sgn_cur;

      send_data_sens       [itrf][2*k  ] = sgn_opp;
      send_data_sens       [itrf][2*k+1] = sgn_cur;


      send_data_intno      [itrf][2*k  ] =  (itrf+1);
      send_data_intno      [itrf][2*k+1] = -(itrf+1);
      weight               [itrf][2*k  ] = 1.;
      weight               [itrf][2*k+1] = 1.;
      stride_one           [itrf][2*k  ] = 1;
      stride_one           [itrf][2*k+1] = 1;
    }

    if (0 == 1) {
      log_trace("Interface %d\n", itrf);
      PDM_log_trace_array_long(interface_ids_shifted[itrf], 2*dn_interface[itrf], "shifted gnum    :: ");
      PDM_log_trace_array_int (send_data_dom        [itrf], 2*dn_interface[itrf], "send_data_dom   :: ");
      PDM_log_trace_array_int (send_data_sens       [itrf], 2*dn_interface[itrf], "send_data_sens  :: ");
      PDM_log_trace_array_int (send_data_intno      [itrf], 2*dn_interface[itrf], "send_data_intno :: ");
      PDM_log_trace_array_long(send_data_gnum       [itrf], 2*dn_interface[itrf], "send_data_gnum  :: ");
    }
  }

  /*
   * part_to_block to equilibrate / fit by block
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      interface_ids_shifted,
                                                      weight,
                                                      dn_interface_twice,
                                                      n_interface,
                                                      comm);

  // Save distribution & gnum from first PtB. We will use it for following PtBs
  int n_gnum = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *gnum   = PDM_part_to_block_block_gnum_get(ptb);
  // PDM_g_num_t *distri = PDM_part_to_block_distrib_index_get(ptb);

  if(0 == 1) {
    PDM_log_trace_array_long(gnum, n_gnum, "gnum ::");
  }


  int *recv_stride   = NULL;
  int *recv_data_dom = NULL;
  int n_connected_l = PDM_part_to_block_exch(ptb,
                                             sizeof(int),
                                             PDM_STRIDE_VAR_INTERLACED,
                                             -1,
                                             stride_one,
                                   (void **) send_data_dom,
                                             &recv_stride,
                                   (void **) &recv_data_dom);

  free(recv_stride);
  recv_stride = NULL;
  int *recv_data_sens = NULL;
  n_connected_l = PDM_part_to_block_exch(ptb,
                                         sizeof(int),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         stride_one,
                               (void **) send_data_sens,
                                         &recv_stride,
                               (void **) &recv_data_sens);

  free(recv_stride);
  recv_stride = NULL;
  int *recv_data_intno = NULL;
  n_connected_l = PDM_part_to_block_exch(ptb,
                                         sizeof(int),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         stride_one,
                               (void **) send_data_intno,
                                         &recv_stride,
                               (void **) &recv_data_intno);

  free(recv_stride);
  recv_stride = NULL;
  PDM_g_num_t *recv_data_gnum = NULL;
  n_connected_l = PDM_part_to_block_exch(ptb,
                                         sizeof(PDM_g_num_t),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         stride_one,
                               (void **) send_data_gnum,
                                         &recv_stride,
                               (void **) &recv_data_gnum);

  if (0 == 1) {
    PDM_log_trace_array_long(gnum           , n_gnum       , "gnum           ::");
    PDM_log_trace_array_int (recv_stride    , n_gnum       , "recv stride    ::");
    PDM_log_trace_array_int (recv_data_dom  , n_connected_l, "recv_data_dom  ::");
    PDM_log_trace_array_int (recv_data_sens , n_connected_l, "recv_data_sens ::");
    PDM_log_trace_array_int (recv_data_intno, n_connected_l, "recv_data_intno::");
    PDM_log_trace_array_long(recv_data_gnum , n_connected_l, "recv_data_gnum ::");
  }

  /*
   * At this stage we have in block frame all entity1 and all occurence for all interface
   */
  int* n_dentity2_entity1 = (int *) malloc( n_domain * sizeof(int));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    n_dentity2_entity1[i_domain] = dentity2_entity1_idx[i_domain][dn_entity2[i_domain]];
    // PDM_log_trace_array_long(dentity2_entity1[i_domain], n_dentity2_entity1[i_domain], "dentity2_entity1");
  }

  for (int itrf = 0; itrf < n_interface; itrf++) {
    free(interface_ids_shifted[itrf]);
    free(send_data_dom        [itrf]);
    free(send_data_sens       [itrf]);
    free(send_data_intno      [itrf]);
    free(send_data_gnum       [itrf]);
    free(weight               [itrf]);
  }
  free(interface_ids_shifted);
  free(send_data_dom        );
  free(send_data_sens       );
  free(send_data_intno      );
  free(send_data_gnum       );
  free(weight               );
  free(dn_interface_twice   );

  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(gnum,
                                                                        n_gnum,
                                             (const PDM_g_num_t    **)  dentity2_entity1,
                                                                        n_dentity2_entity1,
                                                                        n_domain,
                                                                        comm);

  int **part_stride   = NULL;
  int **part_data_dom = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_dom,
                         &part_stride,
             (void ***)  &part_data_dom);
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    free(part_stride[i_domain]);
  }
  free(part_stride);

  int **part_data_sens   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_sens,
                         &part_stride,
             (void ***)  &part_data_sens);

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    free(part_stride[i_domain]);
  }
  free(part_stride);

  int **part_data_intno   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_intno,
                         &part_stride,
             (void ***)  &part_data_intno);

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    free(part_stride[i_domain]);
  }
  free(part_stride);

  PDM_g_num_t **part_data_gnum   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_gnum,
                         &part_stride,
             (void ***)  &part_data_gnum);

  free(recv_data_dom);

  PDM_block_to_part_free(btp);

  /*
   * Post-Treated
   */
  int **part_stride_idx             = (int ** ) malloc( n_domain * sizeof(int *));
  // int  *n_max_possible_entity2_intf = (int *  ) malloc( n_domain * sizeof(int  ));
  int max_size                    = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    part_stride_idx[i_domain] = PDM_array_new_idx_from_sizes_int(part_stride[i_domain], n_dentity2_entity1[i_domain]);
    max_size = PDM_MAX(max_size, part_stride_idx[i_domain][n_dentity2_entity1[i_domain]]);
  }

  PDM_g_num_t** distrib_entity2 = (PDM_g_num_t**) malloc( n_domain * sizeof(PDM_g_num_t *));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    distrib_entity2[i_domain] = PDM_compute_entity_distribution(comm, dn_entity2[i_domain]);
  }


  /*
   * For each face we receive for each vertex all possible interface
   *   - Condition 1 : All vtx of the should have the same interface domain
   */
  int  *l_interface_n        = (int *  ) malloc( n_interface * sizeof(int  ));
  int  *l_interface_sgn      = (int *  ) malloc( n_interface * sizeof(int  ));
  int **entity2_lids         = (int ** ) malloc( n_domain    * sizeof(int *));
  int **entity2_intf_no_idx  = (int ** ) malloc( n_domain    * sizeof(int *));
  int **entity2_intf_no      = (int ** ) malloc( n_domain    * sizeof(int *));
  int  *n_entity2_intf       = (int *  ) malloc( n_domain    * sizeof(int  ));
  int  *key_data_size_approx = (int *  ) malloc( n_domain    * sizeof(int  ));

  if(0 == 1) {
    PDM_log_trace_array_int(dn_entity1, n_domain, "dn_entity1"    );
    PDM_log_trace_array_int(dn_entity2, n_domain, "dn_entity2"    );
  }

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {

    if(0 == 1) {
      PDM_log_trace_array_int(part_stride     [i_domain], n_dentity2_entity1[i_domain], "part_stride"    );
      PDM_log_trace_array_int(part_stride_idx [i_domain], n_dentity2_entity1[i_domain], "part_stride_idx");
      PDM_log_trace_array_int(part_data_sens  [i_domain], part_stride_idx[i_domain][n_dentity2_entity1[i_domain]], "part_data_sens  ::");
      PDM_log_trace_array_int(part_data_intno [i_domain], part_stride_idx[i_domain][n_dentity2_entity1[i_domain]], "part_data_intno ::");
      PDM_log_trace_array_long(part_data_gnum [i_domain], part_stride_idx[i_domain][n_dentity2_entity1[i_domain]], "part_data_gnum  ::");
    }

    key_data_size_approx[i_domain] = 0;

    int *_dentity2_entity1_idx = dentity2_entity1_idx[i_domain];

    entity2_lids       [i_domain] = (int * ) malloc(  dn_entity2[i_domain]    * sizeof(int));
    entity2_intf_no_idx[i_domain] = (int * ) malloc( (dn_entity2[i_domain]+1) * sizeof(int));
    entity2_intf_no    [i_domain] = (int * ) malloc( (max_size              ) * sizeof(int));

    n_entity2_intf     [i_domain] = 0;
    entity2_intf_no_idx[i_domain][0] = 0;

    for(int i_entity2 = 0; i_entity2 < dn_entity2[i_domain]; ++i_entity2) {

      for(int i = 0; i < n_interface; ++i) {
        l_interface_n  [i] = 0;
        l_interface_sgn[i] = 0;
      }

      if(_dentity2_entity1_idx[i_entity2] == _dentity2_entity1_idx[i_entity2+1]) {
        continue;
      }

      for(int j = _dentity2_entity1_idx[i_entity2]; j < _dentity2_entity1_idx[i_entity2+1]; ++j) {
        int idx = part_stride_idx[i_domain][j];
        for(int k = 0; k < part_stride[i_domain][j]; ++k) {
          int itrf = PDM_ABS(part_data_intno[i_domain][idx+k])-1;
          l_interface_n  [itrf]++;
          if(l_interface_sgn[itrf] == 0) {
            l_interface_sgn[itrf] = PDM_SIGN(part_data_intno[i_domain][idx+k]);
          }
          // Ce asset peut planter pour de mauvaise raison sur le cas d'un maillage monocellule par exemple
          // Pour corriger il suffit de dire que si tout les éléments recu sont de sign différents
          // assert(l_interface_sgn[itrf] == PDM_SIGN(part_data_intno[i_domain][idx+k]));
        }
      }

      if( 0 == 1) {
        // PDM_log_trace_array_int(work         , n_unique   , "work          :: ");
        PDM_log_trace_array_int(l_interface_n, n_interface, "l_interface_n :: ");
      }

      // Count for all interface
      int n_connect = _dentity2_entity1_idx[i_entity2+1] - _dentity2_entity1_idx[i_entity2];

      int n_new_interface = 0;
      for(int i = 0; i < n_interface; ++i) {
        if(l_interface_n[i] == n_connect) {
          int idx_write = entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]] + n_new_interface++;
          // entity2_intf_no[i_domain][idx_write] = i;
          assert(l_interface_sgn[i] != 0);
          entity2_intf_no[i_domain][idx_write] = l_interface_sgn[i] * (i+1);
          key_data_size_approx[i_domain] += n_connect * 2;
        }
      }

      if(n_new_interface > 0) {
        entity2_intf_no_idx [i_domain][n_entity2_intf[i_domain]+1] = entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]] + n_new_interface;
        entity2_lids        [i_domain][n_entity2_intf[i_domain]++] = i_entity2;
      }
    }

    entity2_lids       [i_domain] = (int * ) realloc(entity2_lids       [i_domain], ( n_entity2_intf[i_domain]   )                          * sizeof(int) );
    entity2_intf_no_idx[i_domain] = (int * ) realloc(entity2_intf_no_idx[i_domain], ( n_entity2_intf[i_domain]+1 )                          * sizeof(int) );
    entity2_intf_no    [i_domain] = (int * ) realloc(entity2_intf_no    [i_domain], entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]] * sizeof(int) );

    if(0 == 1) {
      PDM_log_trace_array_int(entity2_lids       [i_domain], n_entity2_intf[i_domain]  , " n_entity2_intf ::");
      PDM_log_trace_array_int(entity2_intf_no_idx[i_domain], n_entity2_intf[i_domain]+1, " entity2_intf_no_idx ::");
      PDM_log_trace_array_int(entity2_intf_no    [i_domain], entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]], " n_entity2_intf ::");
    }

  }

  free(l_interface_n);
  free(l_interface_sgn);

  // PDM_log_trace_array_int(key_data_size_approx, n_domain, "key_data_size_approx ::");
  for (int itrf = 0; itrf < n_interface; itrf++) {
    free(stride_one[itrf]);
  }
  free(stride_one);

  /*
   * At this stage we identify all gnum of faces that concern by a domain interface
   * We need now to setup link between this faces
   * For this we hash by connectivity
   */
  PDM_g_num_t mod_g_num_entity1 = entity1_per_block_offset[n_domain] / 4;
  int          *n_lkey       = (int          *) malloc( n_domain * sizeof(int          ));
  PDM_g_num_t **key_ln_to_gn = (PDM_g_num_t **) malloc( n_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t **key_data     = (PDM_g_num_t **) malloc( n_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t **gnum_entity2 = (PDM_g_num_t **) malloc( n_domain * sizeof(PDM_g_num_t *));
  int         **key_data_n   = (int         **) malloc( n_domain * sizeof(int         *));
  double      **key_weight   = (double      **) malloc( n_domain * sizeof(double      *));

  stride_one   = (int         **) malloc( n_domain * sizeof(int         *));

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {

    key_ln_to_gn[i_domain] = (PDM_g_num_t *) malloc( entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]] * sizeof(PDM_g_num_t));
    key_data_n  [i_domain] = (int         *) malloc( entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]] * sizeof(int        ));
    key_weight  [i_domain] = (double      *) malloc( entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]] * sizeof(double     ));
    stride_one  [i_domain] = (int         *) malloc( entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]] * sizeof(int        ));
    gnum_entity2[i_domain] = (PDM_g_num_t *) malloc( entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]] * sizeof(PDM_g_num_t));

    key_data    [i_domain] = (PDM_g_num_t *) malloc( key_data_size_approx[i_domain] * sizeof(PDM_g_num_t));

    PDM_g_num_t *_key_ln_to_gn         = key_ln_to_gn        [i_domain];
    int         *_key_data_n           = key_data_n          [i_domain];
    int         *_stride_one           = stride_one          [i_domain];
    PDM_g_num_t *_key_data             = key_data            [i_domain];
    PDM_g_num_t *_gnum_entity2         = gnum_entity2        [i_domain];
    double      *_key_weight           = key_weight          [i_domain];

    int         *_dentity2_entity1_idx = dentity2_entity1_idx[i_domain];
    PDM_g_num_t *_dentity2_entity1     = dentity2_entity1    [i_domain];

    // PDM_log_trace_connectivity_long(_dentity2_entity1_idx, _dentity2_entity1, dn_entity2[i_domain], "_dentity2_entity1 :: ");

    int idx_write      = 0;
    int idx_write_data = 0;
    for(int idx_entity2 = 0; idx_entity2 < n_entity2_intf[i_domain]; idx_entity2++) {
      int i_entity2 = entity2_lids[i_domain][idx_entity2];

      // log_trace(" -----------------------  %i \n", i_entity2);

      for(int idx_interf = entity2_intf_no_idx[i_domain][idx_entity2]; idx_interf < entity2_intf_no_idx[i_domain][idx_entity2+1]; ++idx_interf) {

        int i_interf = PDM_ABS(entity2_intf_no[i_domain][idx_interf])-1;

        _stride_one  [idx_write] = 1;
        _key_data_n  [idx_write] = 0;
        _key_ln_to_gn[idx_write] = 0;
        _key_weight  [idx_write] = 1.;
        _gnum_entity2[idx_write] = i_entity2 + distrib_entity2[i_domain][i_rank] + entity2_per_block_offset[i_domain] + 1;

        // int first_idx_write_data = idx_write_data;

        for(int idx_entity1 = _dentity2_entity1_idx[i_entity2]; idx_entity1 < _dentity2_entity1_idx[i_entity2+1]; ++idx_entity1) {
          PDM_g_num_t entity1_gnum = _dentity2_entity1[idx_entity1];
          _key_ln_to_gn[idx_write] += PDM_ABS(entity1_gnum);

          _key_data_n[idx_write       ] += 1;
          _key_data  [idx_write_data++]  = entity1_gnum;

          // log_trace(" \t Current part   %i \n", entity1_gnum);
        }

        /* Add opposite contribution */
        for(int j = _dentity2_entity1_idx[i_entity2]; j < _dentity2_entity1_idx[i_entity2+1]; ++j) {
          int idx = part_stride_idx[i_domain][j];
          for(int k = 0; k < part_stride[i_domain][j]; ++k) {
            int t_itrf = PDM_ABS(part_data_intno[i_domain][idx+k])-1;
            if(t_itrf == i_interf) {

              int sens = part_data_sens[i_domain][idx+k];
              // log_trace("sens = %i \n", sens);
              _key_ln_to_gn[idx_write] += PDM_ABS(part_data_gnum[i_domain][idx+k]);

              _key_data_n[idx_write       ] += 1;
              _key_data  [idx_write_data++]  = sens * part_data_gnum[i_domain][idx+k];


              // log_trace(" \t Opposit part part sens = %i / part_data_gnum = %i / key_data = %i \n", sens, part_data_gnum[i_domain][idx+k], sens * part_data_gnum[i_domain][idx+k]);

            } else {

              // log_trace(" \t Skip with sens = %i \n", part_data_sens[i_domain][idx+k]);

            }
          }
        }

        _key_ln_to_gn[idx_write] = _key_ln_to_gn[idx_write] % mod_g_num_entity1 + 1;

        // Sort it
        // PDM_sort_long(&_key_data[first_idx_write_data], NULL, _key_data_n[idx_write]);

        idx_write++;
      }
    }
    assert(idx_write == entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]]);

    n_lkey[i_domain] = idx_write;
    // PDM_log_trace_array_long(_key_ln_to_gn, entity2_intf_no_idx[i_domain][n_entity2_intf[i_domain]], " _key_ln_to_gn ::");

  }

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {

    free(part_data_dom  [i_domain]);
    free(part_data_intno[i_domain]);
    free(part_data_sens [i_domain]);
    free(part_data_gnum [i_domain]);
    free(part_stride    [i_domain]);
    free(part_stride_idx[i_domain]);
  }
  free(part_data_dom);
  free(part_data_intno);
  free(part_data_sens);
  free(part_data_gnum);
  free(part_stride);
  free(part_stride_idx);

  /*
   * Create hash table to find for all entity2 the connection between them
   */
  PDM_part_to_block_t *ptb_hash = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           key_ln_to_gn,
                                                           key_weight,
                                                           n_lkey,
                                                           n_domain,
                                                           comm);

  /*
   * Exchange key data
   */
  int         *dkey_data_n = NULL;
  PDM_g_num_t *dkey_data   = NULL;
  int data_size = PDM_part_to_block_exch(ptb_hash,
                                         sizeof(PDM_g_num_t),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         key_data_n,
                               (void **) key_data,
                                         &dkey_data_n,
                               (void **) &dkey_data);
  free(dkey_data_n); // Useless cause is concatenate

  int *dkey_strid   = NULL;
  int *dkey_intf_no = NULL;
  PDM_part_to_block_exch(ptb_hash,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                         stride_one,
               (void **) entity2_intf_no,
                         &dkey_strid,
               (void **) &dkey_intf_no);
  free(dkey_strid);

  PDM_g_num_t *dkey_gnum_entity2 = NULL;
  PDM_part_to_block_exch(ptb_hash,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                         stride_one,
               (void **) gnum_entity2,
                         &dkey_strid,
               (void **) &dkey_gnum_entity2);
  free(dkey_strid);

  int data_size_n = PDM_part_to_block_exch(ptb_hash,
                                         sizeof(int),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         stride_one,
                               (void **) key_data_n,
                                         &dkey_strid,
                               (void **) &dkey_data_n);

  PDM_UNUSED(data_size);
  PDM_UNUSED(data_size_n);

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    free(key_ln_to_gn[i_domain]);
    free(key_data    [i_domain]);
    free(key_data_n  [i_domain]);
    free(gnum_entity2[i_domain]);
    free(key_weight  [i_domain]);
  }
  free(key_ln_to_gn);
  free(key_data    );
  free(key_data_n  );
  free(gnum_entity2);
  free(key_weight);
  free(n_lkey);


  /*
   * Post-Treatment :
   *   - Dans chaque bucket de clé on a dkey_strid conflict
   *   - Yon need to unified them
   */
  int dn_key = PDM_part_to_block_n_elt_block_get(ptb_hash);

  int* dkey_data_idx = (int * ) malloc( (data_size_n+1) * sizeof(int));
  dkey_data_idx[0] = 0;
  int tmp_idx   = 0;
  dkey_data_idx[0] = 0;
  int max_conflict = 0;
  int max_n_data_conflict = 0;
  int *key_by_interface = PDM_array_zeros_int(n_interface);
  for(int i = 0; i < dn_key; ++i) {
    max_conflict = PDM_MAX(max_conflict, dkey_strid[i]);
    int n_data_conflict = 0;
    for(int j = 0; j < dkey_strid[i]; ++j) {
      dkey_data_idx[tmp_idx+1] = dkey_data_idx[tmp_idx] + dkey_data_n[tmp_idx];
      n_data_conflict += dkey_data_n[tmp_idx];
      int itrf = PDM_ABS(dkey_intf_no[tmp_idx])-1;
      key_by_interface[itrf]++;
      tmp_idx++;
    }
    max_n_data_conflict = PDM_MAX(max_n_data_conflict, n_data_conflict);
  }

  if(0 == 1) {
    log_trace("data_size_n = %i \n", data_size_n);
    PDM_log_trace_array_int(dkey_strid, dn_key, "dkey_strid : ");
    PDM_log_trace_array_int(dkey_data_idx, data_size_n+1, "dkey_data_idx : ");
    PDM_log_trace_array_int(key_by_interface, n_interface, "key_by_interface : ");
  }


  /*
   * Solve conflict
   */
  int         idx_read            = 0;
  int         idx_read_data       = 0;
  int         *is_treated         = (int         *) malloc( max_conflict         * sizeof(int        ));
  int         *conflict_data_idx  = (int         *) malloc((max_conflict+1)      * sizeof(int        ));
  PDM_g_num_t *conflict_sort_data = (PDM_g_num_t *) malloc((max_n_data_conflict) * sizeof(PDM_g_num_t));

  int          *_interface_dn_entity2  = (int          *) malloc( n_interface * sizeof(int           ));
  PDM_g_num_t **_interface_ids_entity2 = (PDM_g_num_t **) malloc( n_interface * sizeof(PDM_g_num_t * ));
  int         **_interface_dom_entity2 = (int         **) malloc( n_interface * sizeof(int         * ));

  /*
   * Allocation (sur alloc )
   */
  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    _interface_dn_entity2 [i_interface] = 0;
    _interface_ids_entity2[i_interface] = (PDM_g_num_t *) malloc( key_by_interface[i_interface] * sizeof(PDM_g_num_t));
    _interface_dom_entity2[i_interface] = (int         *) malloc( key_by_interface[i_interface] * sizeof(int        ));
  }


  for(int i = 0; i < dn_key; ++i) {

    int n_conflict_keys = dkey_strid[i];

    conflict_data_idx[0] = 0;
    for(int j = 0; j < n_conflict_keys; ++j) {
      is_treated[j] = 0;
      conflict_data_idx[j+1] = conflict_data_idx[j] + dkey_data_n[idx_read+j];
    }

    // log_trace("----------------------------------------------- \n");
    // log_trace("%i --> n_conflict = %i \n", i, n_conflict_keys);

    /*
     * On trie tout dans conflict_sort_data
     */
    for(int j = 0; j < n_conflict_keys; ++j) {
      for(int k = conflict_data_idx[j]; k < conflict_data_idx[j+1]; ++k) {
        conflict_sort_data[k] = PDM_ABS(dkey_data[idx_read_data+k]);
      }

      PDM_sort_long(&conflict_sort_data[conflict_data_idx[j]], 0, dkey_data_n[idx_read+j]);

      // PDM_log_trace_array_long(&conflict_sort_data[conflict_data_idx[j]], dkey_data_n[idx_read+j], "conflict_sort_data ::");
    }


    for(int i_conflict = 0; i_conflict < n_conflict_keys; ++i_conflict) {

      if(is_treated[i_conflict] == 1) {
        continue;
      }

      // PDM_g_num_t* data1  = &dkey_data[idx_read_data+conflict_data_idx[i_conflict]];
      PDM_g_num_t* data1  = &conflict_sort_data[conflict_data_idx[i_conflict]];
      int          n_val1 = dkey_data_n[idx_read+i_conflict];

      // PDM_log_trace_array_long(data1, n_val1, "data1 :: ");
      // PDM_g_num_t key1 = 0;
      // for(int k = 0; k < n_val1; ++k) {
      //   key1 += data1[k];
      // }
      // log_trace("key1 = %i \n", key1);

      for(int i_conflict2 = i_conflict+1; i_conflict2 < n_conflict_keys; ++i_conflict2) {

        if(is_treated[i_conflict2] == 1) {
          continue;
        }

        // PDM_g_num_t* data2  = &dkey_data[idx_read_data+conflict_data_idx[i_conflict2]];
        PDM_g_num_t* data2  = &conflict_sort_data[conflict_data_idx[i_conflict2]];
        int          n_val2 = dkey_data_n[idx_read+i_conflict2];

        // PDM_log_trace_array_long(data2, n_val2, "data2 :: ");
        // PDM_g_num_t key2 = 0;
        // for(int k = 0; k < n_val2; ++k) {
        //   key2 += data2[k];
        // }
        // log_trace("key2 = %i \n", key2);

        if(n_val1 != n_val2) {
          continue;
        }

        int is_same = 1;
        for(int k = 0; k < n_val1; k++) {
          if(data1[k] != data2[k]) {
            is_same = 0;
            break;
          }
        }

        if(is_same == 0) {
          continue;
        }

        int int_no1    = dkey_intf_no[idx_read+i_conflict ];
        int int_no2    = dkey_intf_no[idx_read+i_conflict2];

        if(PDM_ABS(int_no1) != PDM_ABS(int_no2)) {
          continue;
        }

        /*
         * Si on arrive ici on a un match !!
         */
        is_treated[i_conflict ] = 1;
        is_treated[i_conflict2] = 1;

        /*
         * Build for each interface the correct array
         */
        // log_trace("int_no1 = %i | int_no2 = %i \n", int_no1, int_no2);
        assert(int_no1 == -int_no2);
        int int_no    = PDM_ABS(dkey_intf_no[idx_read+i_conflict])-1;
        int idx_write = _interface_dn_entity2 [int_no]++;

        PDM_g_num_t gnum1 = dkey_gnum_entity2[idx_read+i_conflict ];
        PDM_g_num_t gnum2 = dkey_gnum_entity2[idx_read+i_conflict2];

        int dom1 = PDM_binary_search_gap_long(gnum1-1, entity2_per_block_offset, n_domain+1);
        int dom2 = PDM_binary_search_gap_long(gnum2-1, entity2_per_block_offset, n_domain+1);

        assert(dom1 != -1);
        assert(dom2 != -1);

        //
        // log_trace("Match !!!! "PDM_FMT_G_NUM" |  "PDM_FMT_G_NUM"\n", gnum1, gnum2);

        // PDM_g_num_t* raw_data1  = &dkey_data[idx_read_data+conflict_data_idx[i_conflict ]];
        // PDM_g_num_t* raw_data2  = &dkey_data[idx_read_data+conflict_data_idx[i_conflict2]];

        /*
         *  Mandatory reminder
         *    data is organized as follow : n_data * 2 : [original data form connecvtivty (signed), opposite data using connectivity between domain (Unnsigned)]
         *    So in case of we identify the two entity to be the same they can have a different sign
         */
        // int n_connect = n_val1/2;
        // PDM_g_num_t* raw_data1_opp = raw_data1 + n_val1/2;
        // PDM_g_num_t* raw_data2_opp = raw_data2 + n_val2/2;


        // PDM_log_trace_array_long(raw_data1    , n_val1, "raw_data1_cur (all) :: ");
        // PDM_log_trace_array_long(raw_data2    , n_val1, "raw_data2_cur (all) :: ");

        /*
         * entity2_entity1 can be implicitly ordered or signed, we need to manage all cases :
         *   - ex1 : face_vtx  -> connectivity_is_signed = 0 and implicit order
         *   - ex2 : edge_vtx  -> connectivity_is_signed = 0 and implicit order
         *   - ex3 : edge_vtx  -> connectivity_is_signed = 1  (an edge can be -6 4 )
         *   - ex3 : face_edge -> connectivity_is_signed = 1
         *
         */
        // int sens = 0;
        // if(connectivity_is_signed == 0) {

        //   // PDM_log_trace_array_long(raw_data1_opp, n_connect, "raw_data1_opp :: ");
        //   // PDM_log_trace_array_long(raw_data2_opp, n_connect, "raw_data2_opp :: ");

        //   if(n_connect == 2) { // Cas particulier pour les edges -_-'

        //     // Comparaison opp1 with cur2
        //     PDM_g_num_t i1 = raw_data1_opp[0];
        //     PDM_g_num_t i2 = raw_data1_opp[1];

        //     PDM_g_num_t j1 = raw_data2[0];
        //     PDM_g_num_t j2 = raw_data2[1];

        //     int lsens1 = 0;
        //     if(i1 == j1) {
        //       lsens1 = 1;
        //     } else {
        //       assert(i1 == j2);
        //       assert(i2 == j1);
        //       lsens1 = -1;
        //     }

        //     // Comparaison opp2 with cur1
        //     i1 = raw_data2_opp[0];
        //     i2 = raw_data2_opp[1];

        //     j1 = raw_data1[0];
        //     j2 = raw_data1[1];

        //     int lsens2 = 0;
        //     if(i1 == j1) {
        //       lsens2 = 1;
        //     } else {
        //       assert(i1 == j2);
        //       assert(i2 == j1);
        //       lsens2 = -1;
        //     }

        //     assert(lsens1 == lsens2);
        //     sens = lsens1;


        //   } else {

        //     /*
        //      * In this case we search permutation inside raw_data
        //      */
        //     int idx_min_cur1 = 0;
        //     int idx_min_cur2 = 0;
        //     int idx_min_opp1 = 0;
        //     int idx_min_opp2 = 0;
        //     PDM_g_num_t first_gnum_cur1 = raw_data1    [0];
        //     PDM_g_num_t first_gnum_cur2 = raw_data2    [0];
        //     PDM_g_num_t first_gnum_opp1 = raw_data1_opp[0];
        //     PDM_g_num_t first_gnum_opp2 = raw_data2_opp[0];


        //     for(int p = 0; p < n_connect; ++p) {
        //       if(raw_data1[p] < first_gnum_cur1) {
        //         idx_min_cur1    = p;
        //         first_gnum_cur1 = raw_data1[p];
        //       }
        //       if(raw_data2[p] < first_gnum_cur2) {
        //         idx_min_cur2    = p;
        //         first_gnum_cur2 = raw_data2[p];
        //       }
        //       if(raw_data1_opp[p] < first_gnum_opp1) {
        //         idx_min_opp1    = p;
        //         first_gnum_opp1 = raw_data1_opp[p];
        //       }
        //       if(raw_data2_opp[p] < first_gnum_opp2) {
        //         idx_min_opp2    = p;
        //         first_gnum_opp2 = raw_data2_opp[p];
        //       }
        //     }

        //     // log_trace("idx_min_cur1 = %i | idx_min_cur2 = %i \n", idx_min_cur1, idx_min_cur2);
        //     // log_trace("idx_min_opp1 = %i | idx_min_opp2 = %i \n", idx_min_opp1, idx_min_opp2);

        //     int next_idx_cur1 = (idx_min_cur1 + 1            ) % n_connect;
        //     // int prev_idx_cur1 = (idx_min_cur1 - 1 + n_connect) % n_connect;
        //     int next_idx_cur2 = (idx_min_cur2 + 1            ) % n_connect;
        //     int prev_idx_cur2 = (idx_min_cur2 - 1 + n_connect) % n_connect;
        //     int next_idx_opp1 = (idx_min_opp1 + 1            ) % n_connect;
        //     // int prev_idx_opp1 = (idx_min_opp1 - 1 + n_connect) % n_connect;
        //     int next_idx_opp2 = (idx_min_opp2 + 1            ) % n_connect;
        //     int prev_idx_opp2 = (idx_min_opp2 - 1 + n_connect) % n_connect;

        //     // Comparaison opp1 with cur2
        //     int lsens1 = 0;
        //     if(raw_data2[next_idx_cur2] == raw_data1_opp[next_idx_opp1]){
        //       lsens1 = 1;
        //     } else {
        //       assert(raw_data2[prev_idx_cur2] == raw_data1_opp[next_idx_opp1]);
        //       lsens1 = -1;
        //     }

        //     // Comparaison opp2 with cur1
        //     int lsens2 = 0;
        //     if(raw_data2_opp[next_idx_opp2] == raw_data1[next_idx_cur1]){
        //       lsens2 = 1;
        //     } else {
        //       assert(raw_data2_opp[prev_idx_opp2] == raw_data1[next_idx_cur1]);
        //       lsens2 = -1;
        //     }

        //     assert(lsens1 == lsens2);
        //     sens = lsens1;

        //     // log_trace("sens = %i \n", sens);
        //   }

        // } else {
        //   // Mise a jour des donnés opposés avec l'orientation courante
        //   // for(int p = 0; p < n_connect; ++p) {
        //   //   int sgn1 = PDM_SIGN(raw_data1[p]);
        //   //   raw_data1_opp[p] = sgn1 * raw_data1_opp[p];
        //   //   int sgn2 = PDM_SIGN(raw_data2[p]);
        //   //   raw_data2_opp[p] = sgn2 * raw_data2_opp[p];
        //   // }

        //   // PDM_log_trace_array_long(raw_data1    , n_connect, "raw_data1_cur :: ");
        //   // PDM_log_trace_array_long(raw_data1_opp, n_connect, "raw_data1_opp :: ");
        //   // PDM_log_trace_array_long(raw_data2    , n_connect, "raw_data2_cur :: ");
        //   // PDM_log_trace_array_long(raw_data2_opp, n_connect, "raw_data2_opp :: ");

        //   int sens1 = 0;
        //   int sens2 = 0;
        //   for(int p = 0; p < n_connect; ++p) {
        //     PDM_g_num_t gnum_opp1 = PDM_ABS (raw_data1_opp[p]);
        //     int         sgn_opp1  = PDM_SIGN(raw_data1_opp[p]);
        //     PDM_g_num_t gnum_opp2 = PDM_ABS (raw_data2_opp[p]);
        //     int         sgn_opp2  = PDM_SIGN(raw_data2_opp[p]);
        //     int lsens1 = 0;
        //     int lsens2 = 0;

        //     // Search brute force
        //     for(int q = 0; q < n_connect; ++q) {
        //       PDM_g_num_t gnum_cur1 = PDM_ABS (raw_data2[q]);
        //       PDM_g_num_t gnum_cur2 = PDM_ABS (raw_data1[q]);
        //       if(gnum_cur1 == gnum_opp1){
        //         int sgn_cur2 = PDM_SIGN (raw_data2[q]);
        //         if(sgn_opp1 != sgn_cur2) { lsens1 = -1;}
        //         else {                     lsens1 =  1;}
        //       }
        //       if(gnum_cur2 == gnum_opp2){
        //         int sgn_cur1 = PDM_SIGN (raw_data1[q]);
        //         if(sgn_opp2 != sgn_cur1) { lsens2 = -1;}
        //         else {                     lsens2 =  1;}
        //       }

        //       // In realease  this test is mandatory (sinon on plante dans l'assert plus bas)
        //       if(lsens1 != 0 && lsens2 != 0) {
        //         break; // Car on a tout trouvé, si on veut tout check il faut commenté cette ligne
        //       }

        //       // log_trace("\t\t gnum_cur1 = %i | gnum_cur2 = %i \n", gnum_cur1, gnum_cur2);
        //       // log_trace("\t\t lsens1 = %i | lsens2 = %i \n", lsens1, lsens2);
        //     }

        //     // log_trace("--> lsens1 = %i | lsens2 = %i \n", lsens1, lsens2);
        //     // if(sens1 != 0) {
        //     //   assert(sens1 == lsens1);
        //     // }
        //     // if(sens2 != 0) {
        //     //   assert(sens2 == lsens2);
        //     // }
        //     sens1 = lsens1;
        //     sens2 = lsens2;
        //   }

        //   log_trace("gnum1 = %i | gnum2 = %i | sens1 = %i | sens2 = %i\n", gnum1, gnum2, sens1, sens2);
        //   // assert(sens1 != 0);
        //   // assert(sens2 != 0);
        //   // assert(sens1 == sens2);
        //   sens = sens1;
        // }
        // assert(sens != 0);

        // Temporary patch - Sens are find after all
        // int sens = 1;
        if(int_no1 > 0) {
          // _interface_ids_entity2[int_no][2*idx_write  ] =         gnum1 - entity2_per_block_offset[dom1];
          // _interface_ids_entity2[int_no][2*idx_write+1] = sens * (gnum2 - entity2_per_block_offset[dom2]);

          _interface_ids_entity2[int_no][2*idx_write  ] = gnum1;
          _interface_ids_entity2[int_no][2*idx_write+1] = gnum2;

          _interface_dom_entity2[int_no][2*idx_write  ] = dom1;
          _interface_dom_entity2[int_no][2*idx_write+1] = dom2;
        } else {
          // _interface_ids_entity2[int_no][2*idx_write  ] =         gnum2 - entity2_per_block_offset[dom2];
          // _interface_ids_entity2[int_no][2*idx_write+1] = sens * (gnum1 - entity2_per_block_offset[dom1]);

          _interface_ids_entity2[int_no][2*idx_write  ] = gnum2;
          _interface_ids_entity2[int_no][2*idx_write+1] = gnum1;

          _interface_dom_entity2[int_no][2*idx_write  ] = dom2;
          _interface_dom_entity2[int_no][2*idx_write+1] = dom1;
        }

      }
    }

    idx_read      += n_conflict_keys;
    idx_read_data += conflict_data_idx[n_conflict_keys];
  }


  free(dkey_data_idx);
  free(is_treated);
  free(conflict_data_idx);
  free(conflict_sort_data);
  free(key_by_interface);


  free(dkey_data_n);
  free(dkey_data);
  free(dkey_strid);
  free(dkey_intf_no);
  free(dkey_gnum_entity2);

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    _interface_ids_entity2[i_interface] = (PDM_g_num_t *) realloc( _interface_ids_entity2[i_interface], 2 * _interface_dn_entity2 [i_interface] * sizeof(PDM_g_num_t));
    _interface_dom_entity2[i_interface] = (int         *) realloc( _interface_dom_entity2[i_interface], 2 * _interface_dn_entity2 [i_interface] * sizeof(int        ));

    if(0 == 1) {
      PDM_log_trace_array_long(_interface_ids_entity2[i_interface], 2 * _interface_dn_entity2 [i_interface], "_interface_ids_entity2 : " );
      PDM_log_trace_array_int (_interface_dom_entity2[i_interface], 2 * _interface_dn_entity2 [i_interface], "_interface_dom_entity2 : " );
    }

  }

  if(0 == 1) {
    PDM_log_trace_connectivity_long(dentity2_entity1_idx[0], dentity2_entity1[0], dn_entity2[0], "dentity2_entity1 :: ");
  }

  PDM_part_to_block_free(ptb_hash);


  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    free(entity2_lids       [i_domain]);
    free(entity2_intf_no_idx[i_domain]);
    free(entity2_intf_no    [i_domain]);
  }
  free(entity2_lids        );
  free(entity2_intf_no_idx );
  free(entity2_intf_no     );
  free(n_entity2_intf      );
  free(key_data_size_approx);


  free(n_dentity2_entity1);

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    free(distrib_entity2[i_domain]);
    free(stride_one[i_domain]);
  }
  free(stride_one);
  free(distrib_entity2);

  /*
   * Update interface_id with the good sign
   *   - Hook dentity2_entity1 for each ids in interface -> block_to_part
   *   - We have a partition of entity1 (stride by entity2 )
   *   - Use the initial ptb of entity1 to get "sens" of interface of entity1
   */

  PDM_g_num_t **all_entity2_distribution = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t *));
  int         **dentity2_entity1_n       = (int         **) malloc(n_domain * sizeof(int         *));
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    all_entity2_distribution[i_domain] = PDM_compute_entity_distribution(comm, dn_entity2[i_domain]);

    dentity2_entity1_n[i_domain] = malloc(dn_entity2[i_domain] * sizeof(int));

    for(int i = 0; i < dn_entity2[i_domain]; ++i) {
      dentity2_entity1_n[i_domain][i] = dentity2_entity1_idx[i_domain][i+1] - dentity2_entity1_idx[i_domain][i];
    }
  }

  int *_interface_dn_entity2_twice = (int  *) malloc(n_interface * sizeof(int ));
  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    _interface_dn_entity2_twice[i_interface] = 2 * _interface_dn_entity2[i_interface];
  }



  PDM_multi_block_to_part_t* mbtp = PDM_multi_block_to_part_create(entity2_per_block_offset,
                                                                   n_domain,
                                            (const PDM_g_num_t **) all_entity2_distribution,
                                            (const PDM_g_num_t **) _interface_ids_entity2,
                                            (const int          *) _interface_dn_entity2_twice,
                                                                   n_interface,
                                                                   comm);

  int         **pentity2_entity1_n = NULL;
  PDM_g_num_t **pentity2_entity1   = NULL;
  PDM_multi_block_to_part_exch2(mbtp,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR_INTERLACED,
                                dentity2_entity1_n,
                 (void **)      dentity2_entity1,
                               &pentity2_entity1_n,
                 (void ***)    &pentity2_entity1);

  int *pn_entity2_entity1 = malloc(n_interface * sizeof(int));

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

    // PDM_log_trace_array_int(pentity2_entity1_n[i_interface], _interface_dn_entity2_twice[i_interface], "pentity2_entity1_n ::");
    int n_data = 0;
    for(int i = 0; i < _interface_dn_entity2_twice[i_interface]; ++i) {
      n_data += pentity2_entity1_n[i_interface][i];
    }

    pn_entity2_entity1[i_interface] = n_data;

    // PDM_log_trace_array_long(pentity2_entity1[i_interface], n_data, "pentity2_entity1 ::");

  }


  PDM_block_to_part_t* btp_sens = PDM_block_to_part_create_from_sparse_block(gnum,
                                                                             n_gnum,
                                                  (const PDM_g_num_t    **)  pentity2_entity1,
                                                                             pn_entity2_entity1,
                                                                             n_interface,
                                                                             comm);

  int **pentity2_entity1_data_n = NULL;
  int **pentity2_entity1_sens   = NULL;
  PDM_block_to_part_exch(btp_sens,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_sens,
                         &pentity2_entity1_data_n,
             (void ***)  &pentity2_entity1_sens);

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    free(pentity2_entity1_data_n[i_interface]);
  }
  free(pentity2_entity1_data_n);

  int **pentity2_entity1_intno = NULL;
  PDM_block_to_part_exch(btp_sens,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_intno,
                         &pentity2_entity1_data_n,
             (void ***)  &pentity2_entity1_intno);

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    free(pentity2_entity1_data_n[i_interface]);
  }
  free(pentity2_entity1_data_n);

  PDM_g_num_t **pentity2_entity1_gnum_opp   = NULL;
  PDM_block_to_part_exch(btp_sens,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_gnum,
                         &pentity2_entity1_data_n,
             (void ***)  &pentity2_entity1_gnum_opp);

  /*
   * Free exchange protocol
   */
  free(recv_stride);
  free(recv_data_sens);
  free(recv_data_intno);
  free(recv_data_gnum);
  PDM_block_to_part_free(btp_sens);
  PDM_part_to_block_free(ptb);
  PDM_multi_block_to_part_free(mbtp);
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    free(all_entity2_distribution[i_domain]);
    free(dentity2_entity1_n      [i_domain]);
  }
  free(all_entity2_distribution);
  free(dentity2_entity1_n);

  /*
   * Post-processing
   */
  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

    int         *linterface_dom_entity2 = _interface_dom_entity2[i_interface];
    PDM_g_num_t *linterface_ids_entity2 = _interface_ids_entity2[i_interface];

    // Compute idx for entity2_entity1
    int *pentity2_entity1_idx      = PDM_array_new_idx_from_sizes_int(pentity2_entity1_n     [i_interface], _interface_dn_entity2_twice [i_interface]);
    int *pentity2_entity1_data_idx = PDM_array_new_idx_from_sizes_int(pentity2_entity1_data_n[i_interface], pn_entity2_entity1[i_interface]);

    PDM_g_num_t *_pentity2_entity1          = pentity2_entity1         [i_interface];
    int         *_pentity2_entity1_n        = pentity2_entity1_n       [i_interface];
    int         *_pentity2_entity1_data_n   = pentity2_entity1_data_n  [i_interface];
    int         *_pentity2_entity1_intno    = pentity2_entity1_intno   [i_interface];
    int         *_pentity2_entity1_sens     = pentity2_entity1_sens    [i_interface];
    PDM_g_num_t *_pentity2_entity1_gnum_opp = pentity2_entity1_gnum_opp[i_interface];

    if(0 == 1) {
      PDM_log_trace_array_int(pentity2_entity1_idx     , _interface_dn_entity2_twice[i_interface]+1, "pentity2_entity1_idx      ::");
      PDM_log_trace_array_int(pentity2_entity1_data_idx, pn_entity2_entity1         [i_interface]+1, "pentity2_entity1_data_idx ::");
      PDM_log_trace_array_int(_pentity2_entity1_data_n , pn_entity2_entity1         [i_interface]  , "pentity2_entity1_data_n ::");
      PDM_log_trace_array_int(_pentity2_entity1_n      , _interface_dn_entity2_twice[i_interface]  , "_pentity2_entity1_n ::");

      int n_data = pentity2_entity1_data_idx[pn_entity2_entity1[i_interface]];
      PDM_log_trace_array_int (_pentity2_entity1_sens    , n_data, "pentity2_entity1_sens      ::");
      PDM_log_trace_array_int (_pentity2_entity1_intno   , n_data, "pentity2_entity1_intno     ::");
      PDM_log_trace_array_long(_pentity2_entity1_gnum_opp, n_data, "_pentity2_entity1_gnum_opp ::");
    }

    int max_connect = 0;
    for(int i = 0; i < _interface_dn_entity2_twice[i_interface]; ++i) {
      max_connect = PDM_MAX(max_connect, pentity2_entity1_idx[i+1] - pentity2_entity1_idx[i]);
    }

    PDM_g_num_t *lentity2_entity1_cur = malloc(max_connect * sizeof(PDM_g_num_t));
    PDM_g_num_t *lentity2_entity1_opp = malloc(max_connect * sizeof(PDM_g_num_t));


    // Traitement des pairs
    for(int i = 0; i < _interface_dn_entity2 [i_interface]; ++i) {

      // log_trace(" i = %i - %i \n", i, 2*(i+1));

      // PDM_g_num_t gnum1 = linterface_ids_entity2[2*i  ];
      // PDM_g_num_t gnum2 = linterface_ids_entity2[2*i+1];

      int idx_read1  = pentity2_entity1_idx[2*i  ];
      int idx_read2  = pentity2_entity1_idx[2*i+1];

      int n_data_cur = _pentity2_entity1_n[2*i  ];
      int n_data_opp = _pentity2_entity1_n[2*i+1];

      // Copy
      for(int j = 0; j < n_data_cur; ++j) {
        lentity2_entity1_cur[j] = _pentity2_entity1[idx_read1+j];
      }
      for(int j = 0; j < n_data_opp; ++j) {
        lentity2_entity1_opp[j] = _pentity2_entity1[idx_read2+j];
      }

      // Apply sens of opposite on current
      assert(n_data_cur == n_data_opp);

      int i_data_cur = 0;
      for(int j = 0; j < n_data_opp; ++j) {
        for(int p = pentity2_entity1_data_idx[idx_read1+j]; p < pentity2_entity1_data_idx[idx_read1+j+1]; ++p) {
          if(PDM_ABS(_pentity2_entity1_intno[p]) == i_interface+1) {
            lentity2_entity1_cur[i_data_cur] = lentity2_entity1_cur[i_data_cur] * _pentity2_entity1_sens[p];
            i_data_cur++;
          }
        }
      }
      assert(i_data_cur == n_data_cur);

      int i_data_opp = 0;
      for(int j = 0; j < n_data_opp; ++j) {
        for(int p = pentity2_entity1_data_idx[idx_read2+j]; p < pentity2_entity1_data_idx[idx_read2+j+1]; ++p) {
          if(PDM_ABS(_pentity2_entity1_intno[p]) == i_interface+1) {
            lentity2_entity1_opp[i_data_opp] = PDM_SIGN(lentity2_entity1_opp[i_data_opp]) * _pentity2_entity1_gnum_opp[p];
            i_data_opp++;
          }
        }
      }
      // log_trace("i_data_opp = %i  | n_data_cur = %i \n", i_data_opp, n_data_cur);
      assert(i_data_opp == n_data_cur);

      if(0 == 1) {
        PDM_log_trace_array_long(lentity2_entity1_cur, n_data_cur, "lentity2_entity1_cur ::");
        PDM_log_trace_array_long(lentity2_entity1_opp, n_data_opp, "lentity2_entity1_opp ::");
      }

      int sens = 0;
      if(connectivity_is_signed == 0) {
        if(n_data_cur == 2) { // Cas particulier pour les edges -_-'

          // Comparaison opp1 with cur2
          PDM_g_num_t i1 = lentity2_entity1_cur[0];
          PDM_g_num_t i2 = lentity2_entity1_cur[1];

          PDM_g_num_t j1 = lentity2_entity1_opp[0];
          PDM_g_num_t j2 = lentity2_entity1_opp[1];

          if(i1 == j1) {
            sens = 1;
          } else {
            assert(i1 == j2);
            assert(i2 == j1);
            sens = -1;
          }
          assert(sens != 0);

        } else {

          /*
           * In this case we search permutation inside raw_data
           */
          int idx_min_cur1 = 0;
          int idx_min_opp1 = 0;
          PDM_g_num_t first_gnum_cur1 = lentity2_entity1_cur[0];
          PDM_g_num_t first_gnum_opp1 = lentity2_entity1_opp[0];


          for(int p = 0; p < n_data_cur; ++p) {
            if(lentity2_entity1_cur[p] < first_gnum_cur1) {
              idx_min_cur1    = p;
              first_gnum_cur1 = lentity2_entity1_cur[p];
            }
            if(lentity2_entity1_opp[p] < first_gnum_opp1) {
              idx_min_opp1    = p;
              first_gnum_opp1 = lentity2_entity1_opp[p];
            }
          }

          // log_trace("idx_min_cur1 = %i | idx_min_opp1 = %i \n", idx_min_cur1, idx_min_opp1);
          // log_trace("idx_min_opp1 = %i | idx_min_opp2 = %i \n", idx_min_opp1, idx_min_opp2);

          int next_idx_cur1 = (idx_min_cur1 + 1             ) % n_data_cur;
          int prev_idx_cur1 = (idx_min_cur1 - 1 + n_data_cur) % n_data_cur;
          int next_idx_opp1 = (idx_min_opp1 + 1             ) % n_data_cur;
          // int prev_idx_opp1 = (idx_min_opp1 - 1 + n_data_cur) % n_data_cur;

          // Comparaison opp1 with cur2
          if(lentity2_entity1_cur[next_idx_cur1] == lentity2_entity1_opp[next_idx_opp1]){
            sens = 1;
          } else {
            assert(lentity2_entity1_cur[prev_idx_cur1] == lentity2_entity1_opp[next_idx_opp1]);
            sens = -1;
          }

          // log_trace("sens = %i \n", sens);
        }
      } else {

        // Brute force to identify sens
        for(int p = 0; p < n_data_cur; ++p) {
          PDM_g_num_t gnum_cur = PDM_ABS (lentity2_entity1_cur[p]);
          // int         sgn_cur  = PDM_SIGN(lentity2_entity1_cur[p]);
          for(int q = 0; q < n_data_cur; ++q) {
            PDM_g_num_t gnum_opp = PDM_ABS (lentity2_entity1_opp[q]);
            int         sgn_opp  = PDM_SIGN(lentity2_entity1_opp[q]);
            if(gnum_cur == gnum_opp && sens == 0) {
              sens = sgn_opp;
            }
          }
        }

        // log_trace("\t \t sens = %i \n", sens);

      }

      assert(sens != 0);

      // Update
      PDM_g_num_t gnum1 = PDM_ABS(linterface_ids_entity2[2*i  ]);
      PDM_g_num_t gnum2 = PDM_ABS(linterface_ids_entity2[2*i+1]);

      // Shift
      int dom1 = linterface_dom_entity2[2*i  ];
      int dom2 = linterface_dom_entity2[2*i+1];

      linterface_ids_entity2[2*i  ] =         gnum1 - entity2_per_block_offset[dom1];
      linterface_ids_entity2[2*i+1] = sens * (gnum2 - entity2_per_block_offset[dom2]);

    }

    if(0 == 1) {
      PDM_log_trace_array_long(linterface_ids_entity2, 2 * _interface_dn_entity2 [i_interface], "linterface_ids_entity2 ::");
      PDM_log_trace_array_int (linterface_dom_entity2, 2 * _interface_dn_entity2 [i_interface], "linterface_dom_entity2 ::");
    }

    free(lentity2_entity1_cur);
    free(lentity2_entity1_opp);

    free(pentity2_entity1_idx);
    free(pentity2_entity1_data_idx);
  }

  free(pn_entity2_entity1);

  /*
   * Free all post-processing array
   */
  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    free(pentity2_entity1_n       [i_interface]);
    free(pentity2_entity1         [i_interface]);
    free(pentity2_entity1_data_n  [i_interface]);
    free(pentity2_entity1_sens    [i_interface]);
    free(pentity2_entity1_intno   [i_interface]);
    free(pentity2_entity1_gnum_opp[i_interface]);
  }
  free(pentity2_entity1_n       );
  free(pentity2_entity1         );
  free(pentity2_entity1_data_n  );
  free(pentity2_entity1_intno   );
  free(pentity2_entity1_gnum_opp);
  free(pentity2_entity1_sens    );


  /*
   * Unshift
   */
  if(n_domain > 1) {
    for(int i_domain = 1; i_domain < n_domain; ++i_domain) {
      for(int i = 0; i < dentity2_entity1_idx[i_domain][dn_entity2[i_domain]]; ++i) {
        int sgn = PDM_SIGN(dentity2_entity1[i_domain][i]);
        dentity2_entity1[i_domain][i] = sgn * ( PDM_ABS(dentity2_entity1[i_domain][i]) - entity1_per_block_offset[i_domain]);
      }
    }
  }


  free(entity1_per_block_offset);
  free(entity2_per_block_offset);
  free(_interface_dn_entity2_twice);

  /*
   * Result assignement
   */
  *interface_dn_entity2  = _interface_dn_entity2;
  *interface_ids_entity2 = _interface_ids_entity2;
  *interface_dom_entity2 = _interface_dom_entity2;

  // exit(1);
  // log_trace("PDM_domain_interface_translate_entity1_entity2 end \n");
}



void
PDM_domain_interface_translate_vtx2face
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *dn_vtx,
 int                     *dn_face,
 int                    **dface_vtx_idx,
 PDM_g_num_t            **dface_vtx
)
{
  PDM_UNUSED(dom_intrf);
  PDM_UNUSED(dn_vtx);
  PDM_UNUSED(dn_face);
  PDM_UNUSED(dface_vtx_idx);
  PDM_UNUSED(dface_vtx);

  assert (dom_intrf != NULL);
  assert (dom_intrf->interface_dn_face == NULL);
  assert (dom_intrf->interface_dn_vtx  != NULL);

  // Simple case is not yet managed, copy to go back to full case
  int **_interface_dom_vtx = NULL;
  if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
    _interface_dom_vtx = (int **) malloc(dom_intrf->n_interface * sizeof(int*));
    for (int i_intrf = 0; i_intrf < dom_intrf->n_interface; i_intrf++) {
      _interface_dom_vtx[i_intrf] = (int *) malloc(2*dom_intrf->interface_dn_vtx[i_intrf]*sizeof(int));
      for (int j = 0; j < dom_intrf->interface_dn_vtx[i_intrf]; j++) {
        _interface_dom_vtx[i_intrf][2*j]   = dom_intrf->interface_dom_vtx[i_intrf][0];
        _interface_dom_vtx[i_intrf][2*j+1] = dom_intrf->interface_dom_vtx[i_intrf][1];
      }

      if(0 == 1) {
        PDM_log_trace_array_int (_interface_dom_vtx[i_intrf], 2 * dom_intrf->interface_dn_vtx[i_intrf], "interf_by_vtx");
        PDM_log_trace_array_long(dom_intrf->interface_ids_vtx[i_intrf], 2 * dom_intrf->interface_dn_vtx[i_intrf], "interf_by_vtx");
      }
    }
  } else {
    _interface_dom_vtx = dom_intrf->interface_dom_vtx;
  }

  /*
   * - part_to_block on interfaces ids (with global shift on vtx)
   */
  // int n_domain = dom_intrf->n_domain;

  PDM_domain_interface_translate_entity1_entity2(dom_intrf->n_domain,
                                                 dom_intrf->n_interface,
                                                 dn_vtx,
                                                 dn_face,
                                                 dom_intrf->interface_dn_vtx,
                                                 _interface_dom_vtx,
                                                 dom_intrf->interface_ids_vtx,
                                                 dface_vtx_idx,
                                                 dface_vtx,
                                                 0, // Connectivity is signed
                                                 dom_intrf->comm,
                                                 &dom_intrf->interface_dn_face,
                                                 &dom_intrf->interface_ids_face,
                                                 &dom_intrf->interface_dom_face);
  dom_intrf->is_result[PDM_BOUND_TYPE_FACE] = 1;

  // Simple case is not yet managed, free working arrays
  if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
    for (int i_intrf = 0; i_intrf < dom_intrf->n_interface; i_intrf++) {
      free(_interface_dom_vtx[i_intrf]);
    }
    free(_interface_dom_vtx);
  }
}



void
PDM_domain_interface_translate_vtx2edge
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *dn_vtx,
 int                     *dn_edge,
 int                    **dedge_vtx_idx,
 PDM_g_num_t            **dedge_vtx
)
{
  PDM_UNUSED(dom_intrf);
  PDM_UNUSED(dn_vtx);
  PDM_UNUSED(dn_edge);
  PDM_UNUSED(dedge_vtx_idx);
  PDM_UNUSED(dedge_vtx);

  assert (dom_intrf != NULL);
  assert (dom_intrf->interface_dn_edge == NULL);
  assert (dom_intrf->interface_dn_vtx  != NULL);

  // Simple case is not yet managed, copy to go back to full case
  int **_interface_dom_vtx = NULL;
  if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
    _interface_dom_vtx = (int **) malloc(dom_intrf->n_interface * sizeof(int*));
    for (int i_intrf = 0; i_intrf < dom_intrf->n_interface; i_intrf++) {
      _interface_dom_vtx[i_intrf] = (int *) malloc(2*dom_intrf->interface_dn_vtx[i_intrf]*sizeof(int));
      for (int j = 0; j < dom_intrf->interface_dn_vtx[i_intrf]; j++) {
        _interface_dom_vtx[i_intrf][2*j]   = dom_intrf->interface_dom_vtx[i_intrf][0];
        _interface_dom_vtx[i_intrf][2*j+1] = dom_intrf->interface_dom_vtx[i_intrf][1];
      }

      if(0 == 1) {
        PDM_log_trace_array_int (_interface_dom_vtx[i_intrf], 2 * dom_intrf->interface_dn_vtx[i_intrf], "interf_by_vtx");
        PDM_log_trace_array_long(dom_intrf->interface_ids_vtx[i_intrf], 2 * dom_intrf->interface_dn_vtx[i_intrf], "interf_by_vtx");
      }
    }
  } else {
    _interface_dom_vtx = dom_intrf->interface_dom_vtx;
  }

  /*
   * - part_to_block on interfaces ids (with global shift on vtx)
   */
  // int n_domain = dom_intrf->n_domain;

  PDM_domain_interface_translate_entity1_entity2(dom_intrf->n_domain,
                                                 dom_intrf->n_interface,
                                                 dn_vtx,
                                                 dn_edge,
                                                 dom_intrf->interface_dn_vtx,
                                                 _interface_dom_vtx,
                                                 dom_intrf->interface_ids_vtx,
                                                 dedge_vtx_idx,
                                                 dedge_vtx,
                                                 0,  // connectivity_is_signed
                                                 dom_intrf->comm,
                                                 &dom_intrf->interface_dn_edge,
                                                 &dom_intrf->interface_ids_edge,
                                                 &dom_intrf->interface_dom_edge);
  dom_intrf->is_result[PDM_BOUND_TYPE_EDGE] = 1;

  // Simple case is not yet managed, free working arrays
  if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
    for (int i_intrf = 0; i_intrf < dom_intrf->n_interface; i_intrf++) {
      free(_interface_dom_vtx[i_intrf]);
    }
    free(_interface_dom_vtx);
  }
}



void
PDM_domain_interface_get
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                   **interface_dn,
 PDM_g_num_t          ***interface_ids,
 int                  ***interface_dom
)
{
  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_FACE) {
    assert (dom_intrf->interface_dn_face != NULL);
    *interface_dn  = dom_intrf->interface_dn_face;
    *interface_ids = dom_intrf->interface_ids_face;
    *interface_dom = dom_intrf->interface_dom_face;
  } else if (interface_kind == PDM_BOUND_TYPE_EDGE) {
    assert (dom_intrf->interface_dn_edge != NULL);
    *interface_dn  = dom_intrf->interface_dn_edge;
    *interface_ids = dom_intrf->interface_ids_edge;
    *interface_dom = dom_intrf->interface_dom_edge;
  } else if (interface_kind == PDM_BOUND_TYPE_VTX) {
    assert (dom_intrf->interface_dn_vtx != NULL);
    *interface_dn  = dom_intrf->interface_dn_vtx;
    *interface_ids = dom_intrf->interface_ids_vtx;
    *interface_dom = dom_intrf->interface_dom_vtx;
  } else  {
    PDM_error(__FILE__, __LINE__, 0, "This kind of entity is not yet supported\n");
  }
}

int
PDM_domain_interface_get_as_graph
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                   **interface_graph_idx,
 PDM_g_num_t           **interface_graph_ids,
 int                   **interface_graph_dom
)
{
  assert (dom_intrf != NULL);
  int          *interface_dn = NULL;
  PDM_g_num_t **interface_ids = NULL;
  int         **interface_dom = NULL;
  if (interface_kind == PDM_BOUND_TYPE_FACE) {
    assert (dom_intrf->interface_dn_face != NULL);
    interface_dn  = dom_intrf->interface_dn_face;
    interface_ids = dom_intrf->interface_ids_face;
    interface_dom = dom_intrf->interface_dom_face;
  } else if (interface_kind == PDM_BOUND_TYPE_EDGE) {
    assert (dom_intrf->interface_dn_edge != NULL);
    interface_dn  = dom_intrf->interface_dn_edge;
    interface_ids = dom_intrf->interface_ids_edge;
    interface_dom = dom_intrf->interface_dom_edge;
  } else if (interface_kind == PDM_BOUND_TYPE_VTX) {
    assert (dom_intrf->interface_dn_vtx != NULL);
    interface_dn  = dom_intrf->interface_dn_vtx;
    interface_ids = dom_intrf->interface_ids_vtx;
    interface_dom = dom_intrf->interface_dom_vtx;
  } else  {
    PDM_error(__FILE__, __LINE__, 0, "This kind of entity is not yet supported\n");
  }

  int graph_dn = _interface_to_graph(dom_intrf->n_interface,
                                     dom_intrf->multidomain_intrf,
                                     interface_dn,
                                     interface_ids,
                                     interface_dom,
                                     interface_graph_idx,
                                     interface_graph_ids,
                                     interface_graph_dom,
                                     dom_intrf->comm);
  return graph_dn;
}

void PDM_domain_interface_free
(
 PDM_domain_interface_t *dom_intrf
)
{
  assert (dom_intrf != NULL);
  if (dom_intrf->ownership == PDM_OWNERSHIP_KEEP) {
    if (dom_intrf->is_result[PDM_BOUND_TYPE_VTX]) {
      for (int i_interface = 0; i_interface < dom_intrf->n_interface; i_interface++) {
        free(dom_intrf->interface_ids_vtx[i_interface]);
        // if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
          free(dom_intrf->interface_dom_vtx[i_interface]);
        // }
      }
      free(dom_intrf->interface_dn_vtx);
      free(dom_intrf->interface_ids_vtx);
      // if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        free(dom_intrf->interface_dom_vtx);
      // }
    }

    if (dom_intrf->is_result[PDM_BOUND_TYPE_EDGE]) {
      for (int i_interface = 0; i_interface < dom_intrf->n_interface; i_interface++) {
        free(dom_intrf->interface_ids_edge[i_interface]);
        // if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
          free(dom_intrf->interface_dom_edge[i_interface]);
        // }
      }
      free(dom_intrf->interface_dn_edge);
      free(dom_intrf->interface_ids_edge);
      // if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        free(dom_intrf->interface_dom_edge);
      // }
    }

    if (dom_intrf->is_result[PDM_BOUND_TYPE_FACE]) {
      for (int i_interface = 0; i_interface < dom_intrf->n_interface; i_interface++) {
        free(dom_intrf->interface_ids_face[i_interface]);
        // if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
          free(dom_intrf->interface_dom_face[i_interface]);
        // }
      }
      free(dom_intrf->interface_dn_face);
      free(dom_intrf->interface_ids_face);
      // if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        free(dom_intrf->interface_dom_face);
      // }
    }
  }

  for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface) {
    if(dom_intrf->translation_vect[i_interface]   != NULL) {
      free(dom_intrf->translation_vect[i_interface]);
      dom_intrf->translation_vect[i_interface] = NULL;
    }
    if(dom_intrf->rotation_direction[i_interface]   != NULL) {
      free(dom_intrf->rotation_direction[i_interface]);
      dom_intrf->rotation_direction[i_interface] = NULL;
    }
    if(dom_intrf->rotation_center[i_interface]   != NULL) {
      free(dom_intrf->rotation_center[i_interface]);
      dom_intrf->rotation_center[i_interface] = NULL;
    }
  }

  free(dom_intrf->translation_vect  );
  free(dom_intrf->rotation_direction);
  free(dom_intrf->rotation_center   );
  free(dom_intrf->rotation_angle    );


  free(dom_intrf);
}

void
PDM_ddomain_interface_to_pdomain_interface
(
 PDM_MPI_Comm                   comm,
 int                            n_interface,
 int                            n_domain,
 PDM_domain_interface_mult_t    multidomain_intrf,
 PDM_bound_type_t               interface_kind,
 int                           *dn_interface,
 PDM_g_num_t                  **interface_ids,
 int                          **interface_dom,
 int                           *n_part,
 int                          **pn_entity,
 PDM_g_num_t                 ***entity_ln_to_gn,
 PDM_part_domain_interface_t   *pditrf
)
{
  PDM_UNUSED(comm);
  PDM_UNUSED(n_interface);
  PDM_UNUSED(n_domain);
  PDM_UNUSED(dn_interface);
  PDM_UNUSED(interface_ids);
  PDM_UNUSED(interface_dom);
  PDM_UNUSED(n_part);
  PDM_UNUSED(pn_entity);
  PDM_UNUSED(entity_ln_to_gn);
  PDM_UNUSED(pditrf);

  // log_trace("PDM_ddomain_interface_to_pdomain_interface \n");

  int i_rank = -1;
  PDM_MPI_Comm_rank(comm, &i_rank);

  /*
   * Find max on interface or partition to rebuild a total distribution
   */
  PDM_g_num_t *max_per_domain_loc = PDM_array_const_gnum(n_domain, 0);
  PDM_g_num_t *max_per_domain     = (PDM_g_num_t *) malloc((n_domain+1) * sizeof(PDM_g_num_t));
  for(int itrf = 0; itrf < n_interface; ++itrf) {
    int dom    = -1;
    int domopp = -1;
    if (dn_interface[itrf] > 0 && multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      dom    = interface_dom[itrf][0];
      domopp = interface_dom[itrf][1];
    }
    for (int k = 0; k < dn_interface[itrf]; k++) {
      if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        dom    = interface_dom[itrf][2*k];
        domopp = interface_dom[itrf][2*k+1];
      }
      max_per_domain_loc[dom   ] = PDM_MAX(max_per_domain_loc[dom   ], PDM_ABS(interface_ids[itrf][2*k  ]));
      max_per_domain_loc[domopp] = PDM_MAX(max_per_domain_loc[domopp], PDM_ABS(interface_ids[itrf][2*k+1]));
    }
  }

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      for(int i_entity = 0; i_entity < pn_entity[i_domain][i_part]; ++i_entity) {
        max_per_domain_loc[i_domain] = PDM_MAX(max_per_domain_loc[i_domain], entity_ln_to_gn[i_domain][i_part][i_entity]);
      }
    }
  }
  max_per_domain[0] = 0;
  PDM_MPI_Allreduce(max_per_domain_loc, &max_per_domain[1], n_domain, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_array_accumulate_gnum(max_per_domain, n_domain+1);
  if (0 == 1) {
    PDM_log_trace_array_long(max_per_domain, n_domain+1, "max per domain");
  }

  /*
   * Shift all ln_to_gn
   */
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      for(int i_entity = 0; i_entity < pn_entity[i_domain][i_part]; ++i_entity) {
        entity_ln_to_gn[i_domain][i_part][i_entity] += max_per_domain[i_domain];
      }
      if(0 == 1) {
        PDM_log_trace_array_long(entity_ln_to_gn[i_domain][i_part], pn_entity[i_domain][i_part], "entity_ln_to_gn  ::" );
      }
    }
  }


  // Prepare first PtB with multiple partitions.
  // Use (shifted) ids as gnum and send tuple (shited) id, opp_id
  PDM_g_num_t **interface_ids_shifted = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
  PDM_g_num_t **send_data_itrf_gnum   = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
  PDM_g_num_t **send_data_gnum        = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
  int         **send_data_itrf_sgn    = (int         **) malloc(n_interface * sizeof(int         *));
  // int         **send_data_dom         = (int         **) malloc(n_interface * sizeof(int         *));
  int         **send_data_intno       = (int         **) malloc(n_interface * sizeof(int         *));
  double      **weight                = (double      **) malloc(n_interface * sizeof(double      *));
  int         **stride_one            = (int         **) malloc(n_interface * sizeof(int         *));
  int          *dn_interface_twice    = (int          *) malloc(n_interface * sizeof(int          ));
  PDM_g_num_t **distrib_itrf          = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
  for (int itrf = 0; itrf < n_interface; itrf++) {
    stride_one           [itrf] = (int         *) malloc( 2 * dn_interface[itrf] * sizeof(int        ));
    interface_ids_shifted[itrf] = (PDM_g_num_t *) malloc( 2 * dn_interface[itrf] * sizeof(PDM_g_num_t));
    send_data_gnum       [itrf] = (PDM_g_num_t *) malloc( 2 * dn_interface[itrf] * sizeof(PDM_g_num_t));
    send_data_itrf_gnum  [itrf] = (PDM_g_num_t *) malloc( 2 * dn_interface[itrf] * sizeof(PDM_g_num_t));
    // send_data_dom        [itrf] = (int         *) malloc( 2 * dn_interface[itrf] * sizeof(int        ));
    send_data_itrf_sgn   [itrf] = (int         *) malloc( 2 * dn_interface[itrf] * sizeof(int        ));
    send_data_intno      [itrf] = (int         *) malloc( 2 * dn_interface[itrf] * sizeof(int        ));
    weight               [itrf] = (double      *) malloc( 2 * dn_interface[itrf] * sizeof(double     ));

    distrib_itrf[itrf] = PDM_compute_entity_distribution(comm, dn_interface[itrf]);

    dn_interface_twice   [itrf] = 2*dn_interface[itrf];
    int dom    = -1;
    int domopp = -1;
    if (dn_interface[itrf] > 0 && multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      dom    = interface_dom[itrf][0];
      domopp = interface_dom[itrf][1];
    }
    for (int k = 0; k < dn_interface[itrf]; k++) {
      if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        dom    = interface_dom[itrf][2*k  ];
        domopp = interface_dom[itrf][2*k+1];
      }
      interface_ids_shifted[itrf][2*k  ] =    PDM_ABS(interface_ids[itrf][2*k  ]) + max_per_domain[dom   ];
      interface_ids_shifted[itrf][2*k+1] = - (PDM_ABS(interface_ids[itrf][2*k+1]) + max_per_domain[domopp]);

      int sgn1 = PDM_SIGN(interface_ids[itrf][2*k  ]);
      int sgn2 = PDM_SIGN(interface_ids[itrf][2*k+1]);


      // send_data_gnum       [itrf][2*k  ] =    PDM_ABS(interface_ids[itrf][2*k  ]) + max_per_domain[dom   ];
      // send_data_gnum       [itrf][2*k+1] = - (PDM_ABS(interface_ids[itrf][2*k+1]) + max_per_domain[domopp]);

      send_data_gnum       [itrf][2*k  ] =  sgn1 * (PDM_ABS(interface_ids[itrf][2*k  ]) + max_per_domain[dom   ]);
      send_data_gnum       [itrf][2*k+1] =  sgn2 * (PDM_ABS(interface_ids[itrf][2*k+1]) + max_per_domain[domopp]);

      // send_data_itrf_sgn       [itrf][2*k  ] =  PDM_SIGN(interface_ids[itrf][2*k+1]); // On stcoke le sens sur le deuxieme uniquement dans l'autre
      // send_data_itrf_sgn       [itrf][2*k+1] = -PDM_SIGN(interface_ids[itrf][2*k+1]);

      send_data_itrf_sgn       [itrf][2*k  ] =  1;
      send_data_itrf_sgn       [itrf][2*k+1] = -1;

      /*
       * Reminder :
       *   - sgn  :  de gauche à droite ou l'inverse --> C'est implicite en block !!!!
       *   - sens : sens de l'entité connecté par rapport à l'autre (Donc edge ou face retourné )
       */

      // send_data_itrf_sgn       [itrf][2*k  ] = PDM_SIGN(interface_ids[itrf][2*k  ]); // On stcoke le sens sur le deuxieme uniquement dans l'autre
      // send_data_itrf_sgn       [itrf][2*k+1] = PDM_SIGN(interface_ids[itrf][2*k+1]);

      // send_data_dom        [itrf][2*k  ] = domopp;
      // send_data_dom        [itrf][2*k+1] = dom   ;
      send_data_itrf_gnum  [itrf][2*k  ] = distrib_itrf[itrf][i_rank] + k + 1;
      send_data_itrf_gnum  [itrf][2*k+1] = distrib_itrf[itrf][i_rank] + k + 1;
      send_data_intno      [itrf][2*k  ] = itrf;
      send_data_intno      [itrf][2*k+1] = itrf;
      weight               [itrf][2*k  ] = 1.;
      weight               [itrf][2*k+1] = 1.;
      stride_one           [itrf][2*k  ] = 1;
      stride_one           [itrf][2*k+1] = 1;
    }

    if (0 == 1) {
      log_trace("Interface %d\n", itrf);
      PDM_log_trace_array_long(interface_ids[itrf], 2 * dn_interface[itrf], "shifted interface_ids    :: ");
      PDM_log_trace_array_long(interface_ids_shifted[itrf], 2*dn_interface[itrf], "shifted gnum    :: ");
      // PDM_log_trace_array_int (send_data_dom        [itrf], 2*dn_interface[itrf], "send_data_dom   :: ");
      PDM_log_trace_array_int (send_data_itrf_sgn   [itrf], 2*dn_interface[itrf], "send_data_itrf_sgn  :: ");
      PDM_log_trace_array_int (send_data_intno      [itrf], 2*dn_interface[itrf], "send_data_intno :: ");
      PDM_log_trace_array_long(send_data_gnum       [itrf], 2*dn_interface[itrf], "send_data_gnum  :: ");
    }

  }

  /*
   *   part_to_block des n_interfaces to have all interface in block frame
   *   send relative position to have the entity_ln_to_gn and i_interface
   *   Caution : We can have multiple interface for the same gnum (exemple for vtx )
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      interface_ids_shifted,
                                                      weight,
                                                      dn_interface_twice,
                                                      n_interface,
                                                      comm);

  int          n_gnum     = PDM_part_to_block_n_elt_block_get  (ptb);
  PDM_g_num_t* block_gnum = PDM_part_to_block_block_gnum_get   (ptb);
  // PDM_g_num_t* distrib    = PDM_part_to_block_distrib_index_get(ptb);

  if (0 == 1) {
    PDM_log_trace_array_long(block_gnum, n_gnum, "gnum");
  }
  int *recv_stride   = NULL;
  int n_connected_l  = 0;
  // int *recv_data_dom = NULL;
  // n_connected_l = PDM_part_to_block_exch(ptb,
  //                                        sizeof(int),
  //                                        PDM_STRIDE_VAR_INTERLACED,
  //                                        -1,
  //                                        stride_one,
  //                              (void **) send_data_dom,
  //                                        &recv_stride,
  //                              (void **) &recv_data_dom);

  // free(recv_stride);
  // recv_stride = NULL;
  int *recv_data_itrf_sgn = NULL;
  n_connected_l = PDM_part_to_block_exch(ptb,
                                         sizeof(int),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         stride_one,
                               (void **) send_data_itrf_sgn,
                                         &recv_stride,
                               (void **) &recv_data_itrf_sgn);

  free(recv_stride);
  recv_stride = NULL;
  int *recv_data_intno = NULL;
  n_connected_l = PDM_part_to_block_exch(ptb,
                                         sizeof(int),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         stride_one,
                               (void **) send_data_intno,
                                         &recv_stride,
                               (void **) &recv_data_intno);

  free(recv_stride);
  recv_stride = NULL;
  PDM_g_num_t *recv_data_gnum = NULL;
  n_connected_l = PDM_part_to_block_exch(ptb,
                                         sizeof(PDM_g_num_t),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         stride_one,
                               (void **) send_data_gnum,
                                         &recv_stride,
                               (void **) &recv_data_gnum);

  free(recv_stride);
  recv_stride = NULL;
  PDM_g_num_t *recv_data_itrf_gnum = NULL;
  n_connected_l = PDM_part_to_block_exch(ptb,
                                         sizeof(PDM_g_num_t),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         -1,
                                         stride_one,
                               (void **) send_data_itrf_gnum,
                                         &recv_stride,
                               (void **) &recv_data_itrf_gnum);

  if (0 == 1) {
    PDM_log_trace_array_long(block_gnum         , n_gnum       , "block_gnum"         );
    PDM_log_trace_array_int (recv_stride        , n_gnum       , "recv stride"        );
    // PDM_log_trace_array_int (recv_data_dom      , n_connected_l, "recv_data_dom"      );
    PDM_log_trace_array_int (recv_data_itrf_sgn , n_connected_l, "recv_data_itrf_sgn  ::");
    PDM_log_trace_array_int (recv_data_intno    , n_connected_l, "recv_data_intno     ::");
    PDM_log_trace_array_long(recv_data_gnum     , n_connected_l, "recv_data_gnum      ::");
    PDM_log_trace_array_long(recv_data_itrf_gnum, n_connected_l, "recv_data_itrf_gnum ::");
  }

  for (int itrf = 0; itrf < n_interface; itrf++) {
    free(interface_ids_shifted[itrf]);
    // free(send_data_dom        [itrf]);
    free(send_data_itrf_sgn   [itrf]);
    free(send_data_intno      [itrf]);
    free(send_data_gnum       [itrf]);
    free(send_data_itrf_gnum  [itrf]);
    free(weight               [itrf]);
    free(stride_one           [itrf]);
  }
  free(interface_ids_shifted);
  // free(send_data_dom        );
  free(send_data_itrf_sgn   );
  free(send_data_intno      );
  free(send_data_gnum       );
  free(send_data_itrf_gnum  );
  free(weight               );
  free(stride_one           );
  free(dn_interface_twice   );

  /*
   *
   */
  int n_part_tot = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    n_part_tot += n_part[i_domain];
  }

  int          *pn_entity_all        = (int          *) malloc( n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity_ln_to_gn_all = (PDM_g_num_t **) malloc( n_part_tot * sizeof(PDM_g_num_t *));

  int shift_domain = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      pn_entity_all       [shift_domain+i_part] = pn_entity      [i_domain][i_part];
      pentity_ln_to_gn_all[shift_domain+i_part] = entity_ln_to_gn[i_domain][i_part];
    }
    shift_domain += n_part[i_domain];
  }


  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(block_gnum,
                                                                        n_gnum,
                                             (const PDM_g_num_t    **)  pentity_ln_to_gn_all,
                                                                        pn_entity_all,
                                                                        n_part_tot,
                                                                        comm);

  /*
   * Exchange to have for each partitions all information
   */
  int **part_stride   = NULL;
  // int **part_data_dom = NULL;
  // PDM_block_to_part_exch(btp,
  //                        sizeof(int),
  //                        PDM_STRIDE_VAR_INTERLACED,
  //                        recv_stride,
  //                        recv_data_dom,
  //                        &part_stride,
  //            (void ***)  &part_data_dom);
  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   free(part_stride[i_part]);
  // }
  // free(part_stride);

  int **part_data_itrf_sgn   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_itrf_sgn,
                         &part_stride,
             (void ***)  &part_data_itrf_sgn);

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(part_stride[i_part]);
  }
  free(part_stride);
  free(recv_data_itrf_sgn);

  int **part_data_intno   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_intno,
                         &part_stride,
             (void ***)  &part_data_intno);

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(part_stride[i_part]);
  }
  free(part_stride);


  PDM_g_num_t **part_data_gnum   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_gnum,
                         &part_stride,
             (void ***)  &part_data_gnum);

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(part_stride[i_part]);
  }
  free(part_stride);

  PDM_g_num_t **part_data_itrf_gnum   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         recv_stride,
                         recv_data_itrf_gnum,
                         &part_stride,
             (void ***)  &part_data_itrf_gnum);


  if(0 == 1) {
    for(int i_part = 0; i_part < n_part_tot; ++i_part) {
      PDM_log_trace_array_int (part_stride    [i_part], pn_entity_all[i_part], "part_stride"    );
      int n_data = 0;
      for(int i = 0; i < pn_entity_all[i_part]; ++i) {
        n_data += part_stride    [i_part][i];
      }
      // PDM_log_trace_array_int (part_data_dom      [i_part], n_data, "part_data_dom"      );
      PDM_log_trace_array_int (part_data_itrf_sgn [i_part], n_data, "part_data_itrf_sgn"     );
      PDM_log_trace_array_int (part_data_intno    [i_part], n_data, "part_data_intno"    );
      PDM_log_trace_array_long(part_data_gnum     [i_part], n_data, "part_data_gnum"     );
      PDM_log_trace_array_long(part_data_itrf_gnum[i_part], n_data, "part_data_itrf_gnum");
    }
  }

  /*
   *  Post-treatment - The first one
   */
  int         **pn_interface        = (int         **) malloc( n_part_tot * sizeof(int         **));
  int         **pn_interface_idx    = (int         **) malloc( n_part_tot * sizeof(int         **));
  int         **pinterface_triplet  = (int         **) malloc( n_part_tot * sizeof(int          *));
  int         **pinterface_itrf_sgn = (int         **) malloc( n_part_tot * sizeof(int          *));
  // int         **pinterface_dom     = (int         **) malloc( n_part_tot * sizeof(int          *));
  int         **pinterface_sens     = (int         **) malloc( n_part_tot * sizeof(int          *));
  PDM_g_num_t **pinterface_gnum     = (PDM_g_num_t **) malloc( n_part_tot * sizeof(PDM_g_num_t  *));

  shift_domain = 0;
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    // int *current_desc = (int *) malloc(3 * n_data * sizeof(int));

    pn_interface      [i_part] = (int          * ) malloc(  n_interface    * sizeof(int          ));
    pn_interface_idx  [i_part] = (int          * ) malloc( (n_interface+1) * sizeof(int          ));

    for(int i = 0; i < n_interface; ++i) {
      pn_interface[i_part][i] = 0;
    }

    // Sort buffer by interface first
    int idx_read = 0;
    for(int i = 0; i < pn_entity_all[i_part]; ++i) {
      for(int k = 0; k < part_stride    [i_part][i]; ++k) {
        int i_interface = part_data_intno[i_part][idx_read++];
        pn_interface[i_part][i_interface]++;
      }
    }

    pn_interface_idx[i_part][0] = 0;
    for(int i = 0; i < n_interface; ++i) {
      pn_interface_idx[i_part][i+1] = pn_interface_idx[i_part][i] + pn_interface[i_part][i];
      pn_interface[i_part][i] = 0;
    }

    // PDM_log_trace_array_int(pn_interface_idx[i_part], n_interface+1, "pn_interface_idx[i_part]");

    pinterface_triplet[i_part] = (int         * ) malloc( pn_interface_idx[i_part][n_interface] * sizeof(int         ));
    // pinterface_dom    [i_part] = (int         * ) malloc( pn_interface_idx[i_part][n_interface] * sizeof(int         ));
    pinterface_sens   [i_part] = (int         * ) malloc( pn_interface_idx[i_part][n_interface] * sizeof(int         ));

    pinterface_itrf_sgn[i_part] = (int         * ) malloc( pn_interface_idx[i_part][n_interface] * sizeof(int         ));


    pinterface_gnum   [i_part] = (PDM_g_num_t * ) malloc( pn_interface_idx[i_part][n_interface] * sizeof(PDM_g_num_t ));

    idx_read = 0;
    for(int i = 0; i < pn_entity_all[i_part]; ++i) {
      for(int k = 0; k < part_stride    [i_part][i]; ++k) {
        int i_interface = part_data_intno[i_part][idx_read];
        int idx_write   = pn_interface_idx[i_part][i_interface] + pn_interface[i_part][i_interface]++;

        // pinterface_dom    [i_part][idx_write] = part_data_dom [i_part][idx_read];
        // pinterface_sens   [i_part][idx_write] = part_data_itrf_sgn[i_part][idx_read];
        pinterface_sens   [i_part][idx_write] = PDM_SIGN(part_data_gnum[i_part][idx_read]);
        pinterface_triplet[i_part][idx_write] = i;

        pinterface_itrf_sgn[i_part][idx_write] = part_data_itrf_sgn[i_part][idx_read];

        PDM_g_num_t check_g_num = PDM_ABS(part_data_gnum[i_part][idx_read]);
        assert(check_g_num == pentity_ln_to_gn_all[i_part][i]);

        pinterface_gnum[i_part][idx_write] = part_data_itrf_gnum[i_part][idx_read];

        idx_read++;
      }
    }

    if(0 == 1) {
      // PDM_log_trace_array_int (part_data_dom     [i_part], pn_interface_idx[i_part][n_interface], "part_data_dom      :: ");
      PDM_log_trace_array_int (part_data_itrf_sgn [i_part], pn_interface_idx[i_part][n_interface], "part_data_itrf_sgn :: ");
      PDM_log_trace_array_int (pinterface_sens    [i_part], pn_interface_idx[i_part][n_interface], "pinterface_sens    :: ");
      PDM_log_trace_array_int (pinterface_itrf_sgn[i_part], pn_interface_idx[i_part][n_interface], "pinterface_itrf_sgn:: ");
      PDM_log_trace_array_int (pinterface_triplet [i_part], pn_interface_idx[i_part][n_interface], "pinterface_triplet :: ");
      PDM_log_trace_array_long(pinterface_gnum    [i_part], pn_interface_idx[i_part][n_interface], "pinterface_gnum    :: ");
    }

  }

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    // free(part_data_dom      [i_part]);
    free(part_data_itrf_sgn [i_part]);
    free(part_data_intno    [i_part]);
    free(part_data_gnum     [i_part]);
    free(part_data_itrf_gnum[i_part]);
    free(part_stride        [i_part]);
  }
  // free(part_data_dom      );
  free(part_data_itrf_sgn     );
  free(part_data_intno    );
  free(part_data_gnum     );
  free(part_data_itrf_gnum);
  free(part_stride        );


  free(pn_entity_all);
  free(pentity_ln_to_gn_all);

  PDM_block_to_part_free(btp);
  free(recv_stride);
  // free(recv_data_dom);
  // free(recv_data_itrf_sgn);
  free(recv_data_intno);
  free(recv_data_gnum);
  free(recv_data_itrf_gnum);

  PDM_part_to_block_free(ptb);

  /*
   *  Get opposite information by rexchange all data of partition interface
   */
  int          ***pres_interface_pn       = ( int          ***) malloc( n_domain * sizeof(int          **));
  PDM_g_num_t ****pres_interface_ln_to_gn = ( PDM_g_num_t ****) malloc( n_domain * sizeof(PDM_g_num_t ***));
  int         ****pres_interface_ids      = ( int         ****) malloc( n_domain * sizeof(int         ***));
  int         ****pres_interface_sgn      = ( int         ****) malloc( n_domain * sizeof(int         ***));
  int         ****pres_interface_sens     = ( int         ****) malloc( n_domain * sizeof(int         ***));
  int         ****pres_interface_ids_idx  = ( int         ****) malloc( n_domain * sizeof(int         ***));
  int         ****pres_interface_dom      = ( int         ****) malloc( n_domain * sizeof(int         ***));

  for( int i_domain = 0; i_domain < n_domain; ++i_domain) {
    pres_interface_pn      [i_domain] = ( int          **) malloc( n_part[i_domain] * sizeof(int          *));
    pres_interface_ln_to_gn[i_domain] = ( PDM_g_num_t ***) malloc( n_part[i_domain] * sizeof(PDM_g_num_t **));
    pres_interface_ids     [i_domain] = ( int         ***) malloc( n_part[i_domain] * sizeof(int         **));
    pres_interface_sgn     [i_domain] = ( int         ***) malloc( n_part[i_domain] * sizeof(int         **));
    pres_interface_sens    [i_domain] = ( int         ***) malloc( n_part[i_domain] * sizeof(int         **));
    pres_interface_ids_idx [i_domain] = ( int         ***) malloc( n_part[i_domain] * sizeof(int         **));
    pres_interface_dom     [i_domain] = ( int         ***) malloc( n_part[i_domain] * sizeof(int         **));

    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {

      pres_interface_pn      [i_domain][i_part] = ( int          *) malloc( n_interface * sizeof(int          ));
      pres_interface_ln_to_gn[i_domain][i_part] = ( PDM_g_num_t **) malloc( n_interface * sizeof(PDM_g_num_t *));
      pres_interface_sgn     [i_domain][i_part] = ( int         **) malloc( n_interface * sizeof(int         *));
      pres_interface_sens    [i_domain][i_part] = ( int         **) malloc( n_interface * sizeof(int         *));
      pres_interface_ids     [i_domain][i_part] = ( int         **) malloc( n_interface * sizeof(int         *));
      pres_interface_ids_idx [i_domain][i_part] = ( int         **) malloc( n_interface * sizeof(int         *));
      pres_interface_dom     [i_domain][i_part] = ( int         **) malloc( n_interface * sizeof(int         *));
    }
  }



  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

    int          *_ln_interface    = (int          *) malloc( n_part_tot * sizeof(int          ));
    PDM_g_num_t **_linterface_gnum = (PDM_g_num_t **) malloc( n_part_tot * sizeof(PDM_g_num_t *));
    for(int i_part = 0; i_part < n_part_tot; ++i_part) {
      int beg = pn_interface_idx[i_part][i_interface];
      _ln_interface   [i_part] = pn_interface_idx[i_part][i_interface+1] - beg;
      _linterface_gnum[i_part] = &pinterface_gnum[i_part][beg];
    }

    PDM_part_to_block_t* ptb_sync_part = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                               PDM_PART_TO_BLOCK_POST_MERGE,
                                                                               1.,
                                                                               _linterface_gnum,
                                                                               distrib_itrf[i_interface],
                                                                               _ln_interface,
                                                                               n_part_tot,
                                                                               comm);

    /*
     * Le bloc est plein normalement
     */
    int          n_gnum_interf     = PDM_part_to_block_n_elt_block_get  (ptb_sync_part);
    // PDM_g_num_t* block_gnum = PDM_part_to_block_block_gnum_get   (ptb);
    PDM_g_num_t* distrib_interf    = PDM_part_to_block_distrib_index_get(ptb_sync_part);
    assert(n_gnum_interf == distrib_interf[i_rank+1] - distrib_interf[i_rank]);

    /*
     * Exch i_part / i_proc
     */
    int **entity_desc      = (int ** ) malloc( n_part_tot * sizeof(int *));
    int **entity_sens      = (int ** ) malloc( n_part_tot * sizeof(int *));
    int **entity_itrf_sgn = (int ** ) malloc( n_part_tot * sizeof(int *));
    int **pstride_one      = (int ** ) malloc( n_part_tot * sizeof(int *));

    int shift = 0;
    for( int i_domain = 0; i_domain < n_domain; ++i_domain) {
      for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
        entity_desc    [shift+i_part] = (int *) malloc( 3 * _ln_interface[shift+i_part] * sizeof(int));
        entity_sens    [shift+i_part] = (int *) malloc(     _ln_interface[shift+i_part] * sizeof(int));
        entity_itrf_sgn[shift+i_part] = (int *) malloc(     _ln_interface[shift+i_part] * sizeof(int));
        pstride_one    [shift+i_part] = (int *) malloc(     _ln_interface[shift+i_part] * sizeof(int));
        int* _entity_desc     = entity_desc     [shift+i_part];
        int* _entity_sens     = entity_sens     [shift+i_part];
        int* _entity_itrf_sgn = entity_itrf_sgn[shift+i_part];
        int beg = pn_interface_idx[shift+i_part][i_interface];
        for(int i = 0; i < _ln_interface[shift+i_part]; ++i) {
          _entity_desc[3*i  ] = i_rank;
          _entity_desc[3*i+1] = i_part + shift;
          _entity_desc[3*i+2] = pinterface_triplet[shift+i_part][beg+i];

          _entity_sens    [i] = pinterface_sens    [shift+i_part][beg+i];
          _entity_itrf_sgn[i] = pinterface_itrf_sgn[shift+i_part][beg+i];

          pstride_one[shift+i_part][i] = 1;
        }
      }
      shift += n_part[i_domain];
    }

    int *blk_strid       = NULL;
    int *blk_entity_desc = NULL;
    int exch_size = PDM_part_to_block_exch(ptb_sync_part,
                                           3 * sizeof(int),
                                           PDM_STRIDE_VAR_INTERLACED,
                                           -1,
                                           pstride_one,
                                 (void **) entity_desc,
                                           &blk_strid,
                                 (void **) &blk_entity_desc);

    free(blk_strid);
    int *blk_entity_sens = NULL;
    int exch_size2 = PDM_part_to_block_exch(ptb_sync_part,
                                            sizeof(int),
                                            PDM_STRIDE_VAR_INTERLACED,
                                            -1,
                                            pstride_one,
                                  (void **) entity_sens,
                                            &blk_strid,
                                  (void **) &blk_entity_sens);
    PDM_UNUSED(exch_size2);

    free(blk_strid);
    int *blk_entity_itrf_sgn = NULL;
    int exch_size3 = PDM_part_to_block_exch(ptb_sync_part,
                                            sizeof(int),
                                            PDM_STRIDE_VAR_INTERLACED,
                                            -1,
                                            pstride_one,
                                  (void **) entity_itrf_sgn,
                                            &blk_strid,
                                  (void **) &blk_entity_itrf_sgn);
    PDM_UNUSED(exch_size3);



    for(int i_part = 0; i_part < n_part_tot; ++i_part) {
      free(entity_desc    [i_part]);
      free(entity_sens    [i_part]);
      free(entity_itrf_sgn[i_part]);
      free(pstride_one    [i_part]);
    }
    free(entity_desc    );
    free(entity_sens    );
    free(entity_itrf_sgn);
    free(pstride_one    );

    if(0 == 1) {
      PDM_log_trace_array_int(blk_strid          ,     n_gnum_interf, "blk_strid           ::");
      PDM_log_trace_array_int(blk_entity_desc    , 3 * exch_size    , "blk_entity_desc     ::");
      PDM_log_trace_array_int(blk_entity_sens    ,     exch_size    , "blk_entity_sens     ::");
      PDM_log_trace_array_int(blk_entity_itrf_sgn,     exch_size    , "blk_entity_itrf_sgn ::");
    }

    /*
     * Renvoi vers les partitions
     */
    PDM_block_to_part_t *btp_sync_part = PDM_block_to_part_create(distrib_itrf[i_interface],
                                           (const PDM_g_num_t **) _linterface_gnum,
                                                                  _ln_interface,
                                                                  n_part_tot,
                                                                  comm);

    int **precv_stride      = NULL;
    int **precv_entity_desc = NULL;
    PDM_block_to_part_exch(btp_sync_part,
                           3 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           blk_strid,
                           blk_entity_desc,
                          &precv_stride,
               (void ***) &precv_entity_desc);


    int **precv_stride_sens  = NULL;
    int **precv_sens         = NULL;
    PDM_block_to_part_exch(btp_sync_part,
                           1 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           blk_strid,
                           blk_entity_sens,
                          &precv_stride_sens,
               (void ***) &precv_sens);

    int **precv_stride_itrf_sgn  = NULL;
    int **precv_itrf_sgn         = NULL;
    PDM_block_to_part_exch(btp_sync_part,
                           1 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           blk_strid,
                           blk_entity_itrf_sgn,
                          &precv_stride_itrf_sgn,
               (void ***) &precv_itrf_sgn);

    /* Change stride */
    for(int i = 0; i < n_gnum_interf; ++i) {
      blk_strid[i] = 1;
    }

    int* _linterface_dom = NULL;
    if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      _linterface_dom = (int * ) malloc(2 * n_gnum_interf * sizeof(int));
      for(int i = 0; i < n_gnum_interf; ++i) {
        _linterface_dom[2*i  ] = interface_dom[i_interface][0];
        _linterface_dom[2*i+1] = interface_dom[i_interface][1];
      }
    } else {
      _linterface_dom = interface_dom[i_interface];
    }

    int **precv_stride_dom  = NULL;
    int **precv_dom         = NULL;
    PDM_block_to_part_exch(btp_sync_part,
                           2 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           blk_strid,
                           _linterface_dom,
                          &precv_stride_dom,
               (void ***) &precv_dom);

    PDM_g_num_t* _lpart_gnum = (PDM_g_num_t * ) malloc(2 * n_gnum_interf * sizeof(PDM_g_num_t));
    for(int j = 0; j < n_gnum_interf; ++j) {
      int dom    = _linterface_dom[2*j  ];
      int domopp = _linterface_dom[2*j+1];

      int sgn1 = PDM_SIGN(interface_ids[i_interface][2*j  ]);
      int sgn2 = PDM_SIGN(interface_ids[i_interface][2*j+1]);

      _lpart_gnum[2*j  ] = sgn1 * (PDM_ABS(interface_ids[i_interface][2*j  ]) + max_per_domain[dom   ]);
      _lpart_gnum[2*j+1] = sgn2 * (PDM_ABS(interface_ids[i_interface][2*j+1]) + max_per_domain[domopp]);

      // ODL
      // _lpart_gnum[2*j  ] =    PDM_ABS(interface_ids[i_interface][2*j  ]) + max_per_domain[dom   ];
      // _lpart_gnum[2*j+1] = - (PDM_ABS(interface_ids[i_interface][2*j+1]) + max_per_domain[domopp]);
    }

    int         **precv_stride_gnum  = NULL;
    PDM_g_num_t **precv_gnum         = NULL;
    PDM_block_to_part_exch(btp_sync_part,
                           2 * sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           blk_strid,
                           _lpart_gnum,
                          &precv_stride_gnum,
               (void ***) &precv_gnum);

    if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      free(_linterface_dom);
    }
    // free(_linterface_sens);
    free(_lpart_gnum);

    free(blk_strid);
    free(blk_entity_desc);
    free(blk_entity_sens);
    free(blk_entity_itrf_sgn);

    PDM_part_to_block_free(ptb_sync_part);
    PDM_block_to_part_free(btp_sync_part);

    /*
     * Post-Treatment
     */

    shift_domain = 0;
    for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
      for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {

        int s_i_part = shift_domain + i_part;
        int n_data = 0;

        int *interface_ids_idx = (int * ) malloc( (_ln_interface[s_i_part] + 1) * sizeof(int ));

        interface_ids_idx[0] = 0;
        for(int i = 0; i < _ln_interface[s_i_part]; ++i) {
          n_data += precv_stride[s_i_part][i];
          interface_ids_idx[i+1] = interface_ids_idx[i] + precv_stride[s_i_part][i];
        }
        int n_data_dom = 0;
        for(int i = 0; i < _ln_interface[s_i_part]; ++i) {
          n_data_dom += precv_stride_dom[s_i_part][i];
        }

        PDM_g_num_t* _pentity_ln_to_gn = (PDM_g_num_t *) malloc(_ln_interface[s_i_part] * sizeof(PDM_g_num_t));
        for(int i = 0; i < _ln_interface[s_i_part]; ++i) {
          int l_num = pinterface_triplet[s_i_part][i];
          _pentity_ln_to_gn[i] = entity_ln_to_gn[i_domain][i_part][l_num];
        }


        if(0 == 1) {
          PDM_log_trace_array_int (precv_stride     [s_i_part] ,     _ln_interface[s_i_part], "precv_stride      ::");
          PDM_log_trace_array_long(_pentity_ln_to_gn           ,     _ln_interface[s_i_part], "_pentity_ln_to_gn ::");
          PDM_log_trace_array_int (precv_entity_desc[s_i_part] , 3 * n_data                 , "precv_entity_desc ::");
          PDM_log_trace_array_int (precv_dom        [s_i_part] , 2 * n_data_dom             , "precv_dom         ::");
          // PDM_log_trace_array_int (precv_sens       [s_i_part] , 2 * n_data_dom             , "precv_sens        ::");
          PDM_log_trace_array_int (precv_itrf_sgn   [s_i_part] , 1 * n_data                 , "precv_itrf_sgn    ::");
          PDM_log_trace_array_int (precv_sens       [s_i_part] , 1 * n_data                 , "precv_sens        ::");
          PDM_log_trace_array_long(precv_gnum       [s_i_part] , 2 * n_data_dom             , "precv_gnum        ::");
          PDM_log_trace_array_int (precv_stride_gnum[s_i_part] , _ln_interface[s_i_part]    , "precv_stride_gnum ::");
          PDM_log_trace_array_int (pinterface_triplet[s_i_part], _ln_interface[s_i_part]    , "pinterface_triplet ::");
        }

        /*
         *  Post-traitement :
         *    - We refind the current entity concers by the domain_interace
         *    - Keep sens information
         */
        int* precv_entity_desc_post = (int * ) malloc( 2 * 3 * n_data          * sizeof(int));
        int* precv_sens_post        = (int * ) malloc( _ln_interface[s_i_part] * sizeof(int));
        int* precv_sgn              = (int * ) malloc( _ln_interface[s_i_part] * sizeof(int));

        int idx_read      = 0;
        int idx_read_desc = 0;
        int idx_write     = 0;
        int shift_beg = pn_interface_idx[shift_domain+i_part][i_interface];

        interface_ids_idx[0] = 0;
        for(int i = 0; i < _ln_interface[s_i_part]; ++i) {

          int         lnum         = pinterface_triplet[s_i_part][shift_beg+i];
          PDM_g_num_t gnum_to_find = entity_ln_to_gn[i_domain][i_part][lnum];

          int pos = -1;
          int sgn =  0;

          // Brute force because useless to sort (list is small =2 )
          for(int k = 0; k < 2 * precv_stride_gnum[s_i_part][i]; ++k) {
            PDM_g_num_t gnum = PDM_ABS(precv_gnum[s_i_part][idx_read+k]);
            if(gnum_to_find == gnum) {
              pos = k;
              sgn = PDM_SIGN(precv_gnum[s_i_part][idx_read+k]);
            }
          }
          // precv_sgn[i] = sgn;
          precv_sens_post[i] = sgn;
          assert(sgn != 0);

          // log_trace(" ------------------------------------------------------------- \n");
          // log_trace(" lnum = %i | gnum_to_find = %i | idx_read = %i | idx_write = %i \n", lnum, (int) gnum_to_find, idx_read, idx_write);
          // log_trace(" i = %i | pos = %i | sgn = %i  \n", i, pos, sgn);
          // PDM_log_trace_array_int(order, 2 * precv_stride_gnum[s_i_part][i], "order :: ");
          assert(pos != -1);

          int idx_first = idx_write;

          /* Determine the sens */
          int sens_cur = 0;
          for(int k = 0; k < precv_stride[s_i_part][i]; ++k) {
            int idx_read2 = idx_read_desc + k;
            int i_cur_proc   = precv_entity_desc[s_i_part][3*(idx_read2)  ];
            int i_cur_part   = precv_entity_desc[s_i_part][3*(idx_read2)+1];
            int i_cur_entity = precv_entity_desc[s_i_part][3*(idx_read2)+2];
            if(i_cur_proc == i_rank && i_cur_part == i_part + shift_domain && i_cur_entity == lnum ) {
              // sens_cur = precv_sens[s_i_part][idx_read2];
              sens_cur = precv_itrf_sgn[s_i_part][idx_read2];
            }
          }
          assert(sens_cur != 0);
          // precv_sens_post[i] = sens_cur;
          precv_sgn[i] = sens_cur;

          interface_ids_idx[i+1] = interface_ids_idx[i];

          int found = 0;
          for(int k = 0; k < precv_stride[s_i_part][i]; ++k) {
            int idx_read2 = idx_read_desc + k;
            int i_cur_proc   = precv_entity_desc[s_i_part][3*(idx_read2)  ];
            int i_cur_part   = precv_entity_desc[s_i_part][3*(idx_read2)+1];
            int i_cur_entity = precv_entity_desc[s_i_part][3*(idx_read2)+2];
            // int lsens        = precv_sens[s_i_part][idx_read2];
            int lsens        = precv_itrf_sgn[s_i_part][idx_read2];

            // log_trace(" i_cur_proc = %i | i_cur_part = %i | i_cur_entity = %i | lnum = %i \n", i_cur_proc, i_cur_part, i_cur_entity, lnum);
            if(i_cur_proc == i_rank && i_cur_part == i_part + shift_domain && i_cur_entity == lnum ) {
              precv_entity_desc_post[3*idx_first  ] = precv_entity_desc[s_i_part][3*(idx_read2)  ];
              precv_entity_desc_post[3*idx_first+1] = precv_entity_desc[s_i_part][3*(idx_read2)+1];
              precv_entity_desc_post[3*idx_first+2] = precv_entity_desc[s_i_part][3*(idx_read2)+2];
              interface_ids_idx[i+1]++;
              found = 1;
            } else if(lsens == -sens_cur){
              precv_entity_desc_post[3*(idx_write+1)  ] = precv_entity_desc[s_i_part][3*(idx_read2)  ];
              precv_entity_desc_post[3*(idx_write+1)+1] = precv_entity_desc[s_i_part][3*(idx_read2)+1];
              precv_entity_desc_post[3*(idx_write+1)+2] = precv_entity_desc[s_i_part][3*(idx_read2)+2];
              idx_write++;
              interface_ids_idx[i+1]++;
            }
          }
          // log_trace("found = %i \n", found);
          // log_trace(" ------------------------------------------------------------- \n");

          assert(found == 1);
          idx_write++;
          idx_read_desc +=     precv_stride     [s_i_part][i];
          idx_read      += 2 * precv_stride_gnum[s_i_part][i];
        }

        // PDM_log_trace_array_int(precv_entity_desc_post, 3 * n_data, "precv_entity_desc ::");
        free(_pentity_ln_to_gn);
        // free(order);

        // interface_ids_idx is transfer to PDM_part_domain_interface
        // Set ptr properly :
        pres_interface_pn      [i_domain][i_part][i_interface] = _ln_interface[s_i_part];
        pres_interface_ln_to_gn[i_domain][i_part][i_interface] = (PDM_g_num_t * ) malloc( _ln_interface[s_i_part] * sizeof(PDM_g_num_t));

        for(int i = 0; i < _ln_interface[s_i_part]; ++i) {
          pres_interface_ln_to_gn[i_domain][i_part][i_interface][i] = _linterface_gnum[s_i_part][i];
        }

        precv_entity_desc_post = (int*) realloc(precv_entity_desc_post, 3 * interface_ids_idx[_ln_interface[s_i_part]] * sizeof(int) );

        // PDM_log_trace_array_int(precv_sgn, _ln_interface[s_i_part], "precv_sgn ::");

        // pres_interface_ids     [i_domain][i_part][i_interface] = precv_entity_desc[s_i_part];
        pres_interface_sgn     [i_domain][i_part][i_interface] = precv_sgn;
        pres_interface_ids     [i_domain][i_part][i_interface] = precv_entity_desc_post;
        pres_interface_ids_idx [i_domain][i_part][i_interface] = interface_ids_idx;
        pres_interface_dom     [i_domain][i_part][i_interface] = precv_dom [s_i_part];
        pres_interface_sens    [i_domain][i_part][i_interface] = precv_sens_post; //precv_sens[s_i_part];

        free(precv_stride         [s_i_part]);
        free(precv_stride_dom     [s_i_part]);
        free(precv_stride_sens    [s_i_part]);
        free(precv_stride_itrf_sgn[s_i_part]);
        // free(precv_dom        [s_i_part]);
        free(precv_sens       [s_i_part]);
        free(precv_itrf_sgn   [s_i_part]);
        free(precv_entity_desc[s_i_part]);
        free(precv_gnum       [s_i_part]);
        free(precv_stride_gnum[s_i_part]);
      }
      shift_domain += n_part[i_domain];
    }

    free(precv_stride     );
    free(precv_stride_dom );
    free(precv_stride_sens);
    free(precv_dom        );
    free(precv_stride_gnum);
    free(precv_gnum       );
    free(precv_sens       );
    free(precv_stride_itrf_sgn);
    free(precv_itrf_sgn   );
    free(precv_entity_desc);

    free(_ln_interface);
    free(_linterface_gnum);

  }

  shift_domain = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {

      for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
        PDM_part_domain_interface_set(pditrf,
                                      interface_kind,
                                      i_domain,
                                      i_part,
                                      i_interface,
                                      pres_interface_pn      [i_domain][i_part][i_interface],
                                      pres_interface_ln_to_gn[i_domain][i_part][i_interface],
                                      pres_interface_sgn     [i_domain][i_part][i_interface],
                                      pres_interface_sens    [i_domain][i_part][i_interface],
                                      pres_interface_ids     [i_domain][i_part][i_interface],
                                      pres_interface_ids_idx [i_domain][i_part][i_interface],
                                      pres_interface_dom     [i_domain][i_part][i_interface]);

        if(0 == 1) {
          PDM_log_trace_array_long(pres_interface_ln_to_gn[i_domain][i_part][i_interface], pres_interface_pn[i_domain][i_part][i_interface], "pres_interface_ln_to_gn ::");
          PDM_log_trace_array_int (pres_interface_sgn     [i_domain][i_part][i_interface], pres_interface_pn[i_domain][i_part][i_interface], "pres_interface_sgn ::");
          PDM_log_trace_array_int (pres_interface_sens    [i_domain][i_part][i_interface], pres_interface_pn[i_domain][i_part][i_interface], "pres_interface_sens ::");
          PDM_log_trace_array_int (pres_interface_dom     [i_domain][i_part][i_interface], pres_interface_pn[i_domain][i_part][i_interface], "pres_interface_dom ::");
          PDM_log_trace_graph_nuplet_int (pres_interface_ids_idx [i_domain][i_part][i_interface],
                                          pres_interface_ids[i_domain][i_part][i_interface], 3, pres_interface_pn[i_domain][i_part][i_interface], "pres_interface_ids ::");
        }
      }
    }
    shift_domain += n_part[i_domain];
  }

  shift_domain = 0;
  for( int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      free(pres_interface_pn      [i_domain][i_part]);
      free(pres_interface_ln_to_gn[i_domain][i_part]);
      free(pres_interface_sgn     [i_domain][i_part]);
      free(pres_interface_sens    [i_domain][i_part]);
      free(pres_interface_ids     [i_domain][i_part]);
      free(pres_interface_ids_idx [i_domain][i_part]);
      free(pres_interface_dom     [i_domain][i_part]);
    }
    free(pres_interface_pn      [i_domain]);
    free(pres_interface_ln_to_gn[i_domain]);
    free(pres_interface_sgn     [i_domain]);
    free(pres_interface_sens    [i_domain]);
    free(pres_interface_ids     [i_domain]);
    free(pres_interface_ids_idx [i_domain]);
    free(pres_interface_dom     [i_domain]);
    shift_domain += n_part[i_domain];
  }
  free(pres_interface_pn      );
  free(pres_interface_ln_to_gn);
  free(pres_interface_sgn     );
  free(pres_interface_sens    );
  free(pres_interface_ids     );
  free(pres_interface_ids_idx );
  free(pres_interface_dom     );

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(pinterface_triplet [i_part]);
    free(pinterface_sens    [i_part]);
    free(pinterface_itrf_sgn[i_part]);
    // free(pinterface_dom    [i_part]);
    free(pinterface_gnum   [i_part]);
    free(pn_interface      [i_part]);
    free(pn_interface_idx  [i_part]);
  }
  free(pn_interface      );
  free(pn_interface_idx  );
  free(pinterface_triplet);
  free(pinterface_sens);
  free(pinterface_itrf_sgn);
  // free(pinterface_dom    );
  free(pinterface_gnum   );


  for (int itrf = 0; itrf < n_interface; itrf++) {
    free(distrib_itrf[itrf]);
  }
  free(distrib_itrf);


  /*
   * Unshift all ln_to_gn
   */
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      for(int i_entity = 0; i_entity < pn_entity[i_domain][i_part]; ++i_entity) {
        entity_ln_to_gn[i_domain][i_part][i_entity] -= max_per_domain[i_domain];
      }
    }
  }

  free(max_per_domain);
  free(max_per_domain_loc);
  // log_trace("PDM_ddomain_interface_to_pdomain_interface end \n");

}

PDM_part_domain_interface_t*
PDM_domain_interface_to_part_domain_interface
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *n_part,
 int                    **pn_face,
 int                    **pn_edge,
 int                    **pn_vtx,
 PDM_g_num_t           ***face_ln_to_gn,
 PDM_g_num_t           ***edge_ln_to_gn,
 PDM_g_num_t           ***vtx_ln_to_gn
)
{
  PDM_UNUSED(dom_intrf);
  PDM_UNUSED(n_part);
  PDM_UNUSED(pn_face);
  PDM_UNUSED(pn_edge);
  PDM_UNUSED(pn_vtx);
  PDM_UNUSED(face_ln_to_gn);
  PDM_UNUSED(edge_ln_to_gn);
  PDM_UNUSED(vtx_ln_to_gn);

  PDM_part_domain_interface_t* pditrf = PDM_part_domain_interface_create(dom_intrf->n_interface,
                                                                         dom_intrf->n_domain,
                                                                         n_part,
                                                                         dom_intrf->multidomain_intrf,
                                                                         PDM_OWNERSHIP_KEEP,
                                                                         dom_intrf->comm);

  if(dom_intrf->interface_dn_vtx != NULL) {
    PDM_ddomain_interface_to_pdomain_interface(dom_intrf->comm,
                                               dom_intrf->n_interface,
                                               dom_intrf->n_domain,
                                               dom_intrf->multidomain_intrf,
                                               PDM_BOUND_TYPE_VTX,
                                               dom_intrf->interface_dn_vtx,
                                               dom_intrf->interface_ids_vtx,
                                               dom_intrf->interface_dom_vtx,
                                               n_part,
                                               pn_vtx,
                                               vtx_ln_to_gn,
                                               pditrf);

  }

  if(dom_intrf->interface_dn_edge != NULL) {
    PDM_ddomain_interface_to_pdomain_interface(dom_intrf->comm,
                                               dom_intrf->n_interface,
                                               dom_intrf->n_domain,
                                               dom_intrf->multidomain_intrf,
                                               PDM_BOUND_TYPE_EDGE,
                                               dom_intrf->interface_dn_edge,
                                               dom_intrf->interface_ids_edge,
                                               dom_intrf->interface_dom_edge,
                                               n_part,
                                               pn_edge,
                                               edge_ln_to_gn,
                                               pditrf);

  }

  if(dom_intrf->interface_dn_face != NULL) {
    PDM_ddomain_interface_to_pdomain_interface(dom_intrf->comm,
                                               dom_intrf->n_interface,
                                               dom_intrf->n_domain,
                                               dom_intrf->multidomain_intrf,
                                               PDM_BOUND_TYPE_FACE,
                                               dom_intrf->interface_dn_face,
                                               dom_intrf->interface_ids_face,
                                               dom_intrf->interface_dom_face,
                                               n_part,
                                               pn_face,
                                               face_ln_to_gn,
                                               pditrf);

  }

  /* Copy information of translation/rotation */
  for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface) {
    double* translation_vect   = NULL;
    double* rotation_direction = NULL;
    double* rotation_center    = NULL;
    double  rotation_angle     = 0;
    PDM_domain_interface_translation_get(dom_intrf, i_interface, &translation_vect);
    PDM_domain_interface_rotation_get   (dom_intrf, i_interface, &rotation_direction, &rotation_center, &rotation_angle);

    if(translation_vect != NULL) {
      PDM_part_domain_interface_translation_set(pditrf, i_interface, translation_vect);
      free(translation_vect);
    }
    if(rotation_direction != NULL){
      PDM_part_domain_interface_rotation_set(pditrf, i_interface, rotation_direction, rotation_center, rotation_angle);
      free(rotation_direction);
      free(rotation_center);
    }
  }

  return pditrf;
}



void
PDM_domain_interface_translation_set
(
        PDM_domain_interface_t  *dom_intrf,
        int                      i_interface,
  const double                  *vect
)
{
  assert(i_interface < dom_intrf->n_interface);
  assert(dom_intrf->translation_vect[i_interface] == NULL);

  dom_intrf->translation_vect[i_interface] = (double *) malloc( 3 * sizeof(double));

  for(int i = 0; i < 3; ++i) {
    dom_intrf->translation_vect[i_interface][i] = vect[i];
  }

}

void
PDM_domain_interface_rotation_set
(
        PDM_domain_interface_t  *dom_intrf,
  const int                      i_interface,
  const double                  *direction,
  const double                  *center,
  const double                   angle
)
{
  assert(i_interface < dom_intrf->n_interface);
  assert(dom_intrf->rotation_direction[i_interface] == NULL);
  assert(dom_intrf->rotation_center   [i_interface] == NULL);

  for(int i = 0; i < 3; ++i) {
    dom_intrf->rotation_direction[i_interface][i] = direction[i];
    dom_intrf->rotation_center   [i_interface][i] = center   [i];
  }
  dom_intrf->rotation_angle[i_interface] = angle;
}


void
PDM_domain_interface_translation_get
(
        PDM_domain_interface_t       *dom_intrf,
        int                           i_interface,
        double                      **vect
)
{
  assert(i_interface < dom_intrf->n_interface);
  if(dom_intrf->translation_vect[i_interface] != NULL){

    *vect = (double *) malloc( 3 * sizeof(double));
    double* _vect = *vect;

    for(int i = 0; i < 3; ++i) {
      _vect[i] = dom_intrf->translation_vect[i_interface][i];
    }
  } else {
    *vect = NULL;
  }
}

void
PDM_domain_interface_rotation_get
(
        PDM_domain_interface_t       *dom_intrf,
  const int                           i_interface,
        double                      **direction,
        double                      **center,
        double                       *angle
)
{
  assert(i_interface < dom_intrf->n_interface);
  if(dom_intrf->rotation_direction[i_interface] != NULL) {
    assert(dom_intrf->rotation_center   [i_interface] != NULL);

    *direction = (double *) malloc( 3 * sizeof(double));
    *center    = (double *) malloc( 3 * sizeof(double));
    double *_direction = *direction;
    double *_center    = *center   ;

    for(int i = 0; i < 3; ++i) {
      _direction[i] = dom_intrf->rotation_direction[i_interface][i];
      _center   [i] = dom_intrf->rotation_center   [i_interface][i];
    }
    *angle = dom_intrf->rotation_angle[i_interface];
  } else {
    *direction = NULL;
    *center    = NULL;
    *angle     = 0;
  }
}


void
PDM_domain_interface_make_flat_view
(
  PDM_domain_interface_t  *dom_intrf,
  PDM_bound_type_t         interface_kind,
  PDM_g_num_t             *shift_by_domain,
  PDM_part_to_block_t   ***ptb_interface_out,
  PDM_g_num_t           ***entity_opp_gnum_out
)
{

  int          *interface_dn  = NULL;
  PDM_g_num_t **interface_ids = NULL;
  int         **interface_dom = NULL;
  PDM_domain_interface_get(dom_intrf,
                           interface_kind,
                           &interface_dn,
                           &interface_ids,
                           &interface_dom);

  int n_interface = dom_intrf->n_interface;

  PDM_g_num_t **interface_ids_shifted = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
  PDM_g_num_t **send_data             = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
  int         **stride_one            = (int         **) malloc(n_interface * sizeof(int         *));
  PDM_g_num_t **entity_opp_gnum       = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));

  PDM_part_to_block_t  **ptb_interface = (PDM_part_to_block_t ** ) malloc(n_interface * sizeof(PDM_part_to_block_t *));

  for (int itrf = 0; itrf < n_interface; itrf++) {
    // stride_one           [itrf] = (int         *) malloc(2*interface_dn[itrf]*sizeof(int        ));
    interface_ids_shifted[itrf] = (PDM_g_num_t *) malloc( 2 * interface_dn[itrf] * sizeof(PDM_g_num_t));
    send_data            [itrf] = (PDM_g_num_t *) malloc( 2 * interface_dn[itrf] * sizeof(PDM_g_num_t));
    stride_one           [itrf] = (int         *) malloc( 2 * interface_dn[itrf] * sizeof(int        ));

    int dom    = -1;
    int domopp = -1;
    if(dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      dom    = interface_dom[itrf][0];
      domopp = interface_dom[itrf][1];
    }

    for (int k = 0; k < interface_dn[itrf]; k++) {
      if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        dom    = interface_dom[itrf][2*k  ];
        domopp = interface_dom[itrf][2*k+1];
      }

      PDM_g_num_t gnum1 = PDM_ABS(interface_ids[itrf][2*k  ]) + shift_by_domain[dom   ];
      PDM_g_num_t gnum2 = PDM_ABS(interface_ids[itrf][2*k+1]) + shift_by_domain[domopp];
      int         sgn1  = PDM_SIGN(interface_ids[itrf][2*k  ]);
      int         sgn2  = PDM_SIGN(interface_ids[itrf][2*k+1]);

      interface_ids_shifted[itrf][2*k  ] = gnum1;
      interface_ids_shifted[itrf][2*k+1] = gnum2;
      send_data            [itrf][2*k  ] = sgn2 * gnum2;
      send_data            [itrf][2*k+1] = sgn1 * gnum1;
      stride_one           [itrf][2*k  ] = 1;
      stride_one           [itrf][2*k+1] = 1;
    }

    int dn_interface_twice = 2 * interface_dn[itrf];

    PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        &interface_ids_shifted[itrf],
                                                        NULL,
                                                        &dn_interface_twice,
                                                        1,
                                                        dom_intrf->comm);

    int         *recv_stride = NULL;
    PDM_g_num_t *recv_data   = NULL;
    int n_connected_l = PDM_part_to_block_exch(ptb,
                                               sizeof(PDM_g_num_t),
                                               PDM_STRIDE_VAR_INTERLACED,
                                               -1,
                                               &stride_one[itrf],
                                     (void **) &send_data[itrf],
                                               &recv_stride,
                                     (void **) &recv_data);

    int n_gnum = PDM_part_to_block_n_elt_block_get(ptb);

    assert(n_gnum == n_connected_l); // ie all recv_stride == 1

    if (0 == 1) {
      PDM_log_trace_array_long(PDM_part_to_block_block_gnum_get(ptb), n_gnum, "gnum");
      PDM_log_trace_array_int (recv_stride, n_gnum       , "recv_stride ::");
      PDM_log_trace_array_long(recv_data  , n_connected_l, "recv_data   ::");
    }

    free(recv_stride);
    free(send_data[itrf]);

    ptb_interface  [itrf] = ptb;
    entity_opp_gnum[itrf] = recv_data;

    // PDM_part_to_block_free(ptb);

  }

  *ptb_interface_out   = ptb_interface;
  *entity_opp_gnum_out = entity_opp_gnum;


  for (int itrf = 0; itrf < n_interface; itrf++) {
    free(interface_ids_shifted[itrf]);
    free(stride_one           [itrf]);
  }
  free(interface_ids_shifted);
  free(stride_one           );
  free(send_data           );



}







#ifdef __cplusplus
}
#endif /* __cplusplus */
