/*
 * \file
 */

#ifndef __PDM_PARA_GRAPH_DUAL_H__
#define __PDM_PARA_GRAPH_DUAL_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"

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
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

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
PDM_para_graph_compress_connectivity_dual
(
      int          n_node,
      PDM_g_num_t  shift_rank,
      PDM_g_num_t *dual_graph_idx,
const int         *dual_graph_n,
      PDM_g_num_t *dual_graph
);

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
);

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
);


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
);

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
);

void
PDM_para_graph_dual_from_combine_connectivity
(
const PDM_MPI_Comm     comm,
const PDM_g_num_t     *cell_distrib,
const PDM_g_num_t     *face_distrib,
      PDM_g_num_t     *vtx_distrib,
const int             *dcell_face_idx,
const PDM_g_num_t     *dcell_face,
const int             *dface_vtx_idx,
const PDM_g_num_t     *dface_vtx,
      PDM_g_num_t    **dual_graph_idx,
      PDM_g_num_t    **dual_graph
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARA_GRAPH_DUAL_H__ */
