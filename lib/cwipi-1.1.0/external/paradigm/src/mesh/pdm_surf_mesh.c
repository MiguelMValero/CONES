/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_surf_mesh.h"
#include "pdm_surf_mesh_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_surf_part.h"
#include "pdm_surf_part_priv.h"
#include "pdm_binary_search.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 * \brief Compute global number of entities
 *
 * \param [in]  mesh         Mesh 
 *
 */

static 
void
_n_g_enttities_compute
(
PDM_surf_mesh_t *mesh
)
{

  if (mesh->nGVtx < 0) {

    PDM_g_num_t n_g_num_vtx = 0;

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_surf_part_t *part = mesh->part[i];
      for (int j = 0; j < part->n_vtx; j++) {
        n_g_num_vtx = PDM_MAX (n_g_num_vtx, part->vtx_ln_to_gn[j]);
      }
    }

    PDM_MPI_Allreduce(&n_g_num_vtx, &(mesh->nGVtx), 1,
                      PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, mesh->comm);
  }

  if (mesh->nGFace < 0) {

    PDM_g_num_t n_g_num_face = 0;

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_surf_part_t *part = mesh->part[i];
      for (int j = 0; j < part->n_face; j++) {
        n_g_num_face = PDM_MAX (n_g_num_face, part->face_ln_to_gn[j]);
      }
    }

    PDM_MPI_Allreduce(&n_g_num_face, &(mesh->nGFace), 1,
                      PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, mesh->comm);
  }
 
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Return an intialized \ref PDM_surf_mesh_t structure
 *
 * This function returns an initialized \ref PDM_surf_mesh_t structure
 *
 * \param [in]  n_part        Number of partition
 * \param [in]  comm         MSG communicator of mesh
 *
 * \return      A new initialized \ref PDM_surf_mesh_t structure
 *
 */

PDM_surf_mesh_t *
PDM_surf_mesh_create
(
const int    n_part,
PDM_MPI_Comm comm
)
{

  PDM_surf_mesh_t *mesh = malloc (sizeof(PDM_surf_mesh_t));

  mesh->nGFace  = -1;
  mesh->nGVtx   = -1;

  mesh->comm    = comm;
  mesh->n_part   = n_part;
  mesh->part    = (PDM_surf_part_t **) malloc(n_part * sizeof(PDM_surf_part_t *));

  PDM_MPI_Allreduce ((void *)&n_part, (void *)&(mesh->nGPart), 1,
                     PDM_MPI_INT, PDM_MPI_SUM, mesh->comm);

  mesh->gMinCarLgthVtx = DBL_MAX;
  mesh->gMaxCarLgthVtx = -DBL_MAX;
  mesh->interPartEdgeGraph = NULL;
  mesh->interPartVtxGraph = NULL;
  mesh->vtxPartBound = NULL;
  mesh->edgePartBound = NULL;

  return (PDM_surf_mesh_t *) mesh;
}


/**
 * \brief Delete a \ref _mesh_t structure
 *
 * This function returns an initialized \ref _mesh_t structure
 *
 * \param [in]  mesh         Mesh to delete
 *
 * \return     Null pointer
 *
 */

PDM_surf_mesh_t *
PDM_surf_mesh_free
(
PDM_surf_mesh_t *mesh
)
{

  if (mesh != NULL) {

    if (mesh->part != NULL) {
      for (int i = 0; i < mesh->n_part; i++)
        mesh->part[i] = PDM_surf_part_free(mesh->part[i]);
      free(mesh->part);
      mesh->part = NULL;

    }

    mesh->interPartEdgeGraph = PDM_graph_bound_free (mesh->interPartEdgeGraph);
    mesh->interPartVtxGraph = PDM_graph_bound_free (mesh->interPartVtxGraph);

    free (mesh->vtxPartBound);
    free (mesh->edgePartBound);

    free (mesh);
  }

  return NULL;
}


/**
 * \brief Compute edge global numbering
 *
 * This function computes edge global numbering
 *
 * \param [in]  mesh    mesh to compute global numbering
 *
 */

void
PDM_surf_mesh_build_edges_gn_and_edge_part_bound
(
 PDM_surf_mesh_t *mesh
)
{
  assert (mesh != NULL);

  int n_part = mesh->n_part;

  int myRank;
  int lComm;

  PDM_MPI_Comm_rank(mesh->comm, &myRank);
  PDM_MPI_Comm_size(mesh->comm, &lComm);

  /*
   * Compute global numbering for internal edge of each partition and
   * define hash table to look for inter partition boundary edges
   * (intra processus)
   */

  int nEdgeProc = 0;
  int n_vtxProc = 0;
  int *nIntEdgePart = (int *) malloc(sizeof(int) * (n_part + 1));
  nIntEdgePart[0] = 0;
  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part =  mesh->part[i];
    nIntEdgePart[i+1] = 0;
    part->edgeLnToGn = PDM_array_const_gnum(part->nEdge, -1);
    n_vtxProc +=  part->n_vtx;
    nEdgeProc += part->nEdge;
  }

  int  keyMax        = 2 * n_vtxProc;
  int  lHashTableIdx = keyMax + 1;
  int *hashTableIdx  = PDM_array_zeros_int(lHashTableIdx);
  int *hashTable     = (int *) malloc(sizeof(int) * 2 * nEdgeProc);
  int *nHashTable    = PDM_array_zeros_int(keyMax);

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    for (int j = 0; j < part->nEdge; j++) {
      if ((part->edgeFace[2*j] == 0) ||
          (part->edgeFace[2*j + 1] == 0)) {
        int vtx1 = part->edgeVtx[2*j];
        int vtx2 = part->edgeVtx[2*j + 1];
        PDM_g_num_t _keyEdge = (part->vtx_ln_to_gn[vtx1-1] + part->vtx_ln_to_gn[vtx2-1]) % keyMax;
        int keyEdge = (int) _keyEdge;
        hashTableIdx[keyEdge+1] += 2;
      }
      else
        nIntEdgePart[i+1] += 1;
    }
  }

  for (int i = 1; i < lHashTableIdx; i++) {
    hashTableIdx[i] = hashTableIdx[i] + hashTableIdx[i-1];
  }

  for (int i = 1; i < n_part + 1; i++) {
    nIntEdgePart[i] = nIntEdgePart[i] + nIntEdgePart[i-1];
  }

  int nIntEdgeProc = nIntEdgePart[n_part];

  int *nEdgeBoundPart = (int *) malloc(sizeof(int) * n_part);

  for (int i = 0; i < 2*nEdgeProc; i++) {
    hashTable[i]=-1;
  }

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    nEdgeBoundPart[i] = 0;
    for (int j = 0; j < part->nEdge; j++) {
      if ((part->edgeFace[2*j] == 0) ||
          (part->edgeFace[2*j + 1] == 0)) {
        int vtx1 = part->edgeVtx[2*j];
        int vtx2 = part->edgeVtx[2*j + 1];
        PDM_g_num_t _keyEdge =
           (part->vtx_ln_to_gn[vtx1-1] + part->vtx_ln_to_gn[vtx2-1]) % keyMax;
        int keyEdge = (int) _keyEdge;
        hashTable[hashTableIdx[keyEdge] + (nHashTable[keyEdge]++)] = i;
        hashTable[hashTableIdx[keyEdge] + (nHashTable[keyEdge]++)] = j;
        nEdgeBoundPart[i] += 1;
      }
      else {
        part->edgeLnToGn[j] = ++(nIntEdgePart[i]);
      }
    }
  }

  free(nHashTable);
  free(nIntEdgePart);

  /*
   * Compute global numbering for inter partition boundary edges (intra processus)
   * Build exchange communication graph between boundary edges
   */

  int *edgePartCur = (int *) malloc(sizeof(int) * mesh->n_part);

  for (int i = 0; i < n_part; i++) {
    edgePartCur[i] = 0;
    PDM_surf_part_t *part =  mesh->part[i];
    const int nConnectedElt = 1;
    const int nOfferElt = 1;
    part->edgePartBound = PDM_part_bound_create (lComm,
                                                 part->nEdge,
                                                 nEdgeBoundPart[i],
                                                 PDM_PART_BOUND_SIMPLE,
                                                 &nConnectedElt,
                                                 &nOfferElt,
                                                 mesh->nGFace,
                                                 part->n_face,
                                                 part->face_ln_to_gn);
  }

  /* for (int i = 0; i < keyMax; i++) { */
  /*   int idx          = hashTableIdx[i]; */
  /*   int nEdgeSameSum = (hashTableIdx[i+1] - idx); */
  /*   for (int j = idx; j < idx + nEdgeSameSum; j+=2) { */
  /*     int i_part = hashTable[j]; */
  /*     int iEdge = hashTable[j+1]; */
  /*     _part_t *part = mesh->part[i_part]; */
  /*   } */
  /* } */

  for (int i = 0; i < keyMax; i++) {
    int idx          = hashTableIdx[i];
    int nEdgeSameSum = (hashTableIdx[i+1] - idx);
    for (int j = idx; j < idx + nEdgeSameSum; j+=2) {
      int i_part = hashTable[j];
      int iEdge = hashTable[j+1];
      if (iEdge != -1) {
        PDM_surf_part_t *part = mesh->part[i_part];
        PDM_part_bound_t *edgePartBound = part->edgePartBound;
        PDM_g_num_t gVtx1 = part->vtx_ln_to_gn[part->edgeVtx[2*iEdge    ] - 1];
        PDM_g_num_t gVtx2 = part->vtx_ln_to_gn[part->edgeVtx[2*iEdge + 1] - 1];
        for (int k = j + 2; k < idx + nEdgeSameSum; k+=2) {
          int iPart1 = hashTable[k];
          int iEdge1 = hashTable[k+1];
          if (iEdge1 != -1) {
            PDM_surf_part_t *part1 = mesh->part[iPart1];
            PDM_part_bound_t *edgePartBound1 = part1->edgePartBound;
            PDM_g_num_t gVtx11 = part1->vtx_ln_to_gn[part1->edgeVtx[2*iEdge1  ] - 1];
            PDM_g_num_t gVtx12 = part1->vtx_ln_to_gn[part1->edgeVtx[2*iEdge1+1] - 1];

            if (((gVtx1 == gVtx11) && (gVtx2 == gVtx12)) ||
                ((gVtx2 == gVtx11) && (gVtx1 == gVtx12))) {

              hashTable[k] = -1;
              hashTable[k+1] = -1;

              PDM_g_num_t gEdge = ++nIntEdgeProc;
              part->edgeLnToGn[iEdge] = gEdge;
              part1->edgeLnToGn[iEdge1] = gEdge;

              edgePartCur[i_part] += 1;
              PDM_part_bound_local_elt_set (edgePartBound,
                                            edgePartCur[i_part],
                                            iEdge+1);

              PDM_part_bound_offer_elt_set (edgePartBound,
                                            edgePartCur[i_part],
                                            0,
                                            PDM_ABS (part->edgeFace[2*iEdge]),
                                            part->face_ln_to_gn[PDM_ABS(part->edgeFace[2*iEdge])-1]);

              PDM_part_bound_distant_elt_set (edgePartBound,
                                              edgePartCur[i_part],
                                              0,
                                              myRank,
                                              iPart1,
                                              iEdge1+1);

              edgePartCur[iPart1] += 1;
              PDM_part_bound_local_elt_set (edgePartBound1,
                                            edgePartCur[iPart1],
                                            iEdge1+1);

              PDM_part_bound_offer_elt_set (edgePartBound1,
                                            edgePartCur[iPart1],
                                            0,
                                            PDM_ABS (part1->edgeFace[2*iEdge1]),
                                            part1->face_ln_to_gn[PDM_ABS (part1->edgeFace[2*iEdge1])-1]);

              PDM_part_bound_distant_elt_set (edgePartBound1,
                                              edgePartCur[iPart1],
                                              0,
                                              myRank,
                                              i_part,
                                              iEdge+1);
              break;
            }
          }
        }
      }
      hashTable[j+1] = -1;
    }
  }
  free(hashTable);
  free(hashTableIdx);
  free(nEdgeBoundPart);

  /*
   * Update global numbering
   */

  PDM_g_num_t *nIntEdgeProcs =
    (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (lComm + 1));
  nIntEdgeProcs[0] = 0;
  PDM_g_num_t *_nIntEdgeProcs = nIntEdgeProcs + 1;

  PDM_g_num_t _nIntEdgeProc = (PDM_g_num_t) nIntEdgeProc;
  PDM_MPI_Allgather((void *) &_nIntEdgeProc, 1, PDM__PDM_MPI_G_NUM,
                (void *) _nIntEdgeProcs, 1, PDM__PDM_MPI_G_NUM, mesh->comm);

  int nEdgeWithoutNG = 0;
  int *edgeWithoutNG = (int *) malloc(sizeof(int) * (2 * nEdgeProc));

  for (int i = 1; i < lComm + 1; i++)
    nIntEdgeProcs[i] = nIntEdgeProcs[i] +  nIntEdgeProcs[i-1];

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    for (int j = 0; j < part->nEdge; j++) {
      if (part->edgeLnToGn[j] != -1)
        part->edgeLnToGn[j] += nIntEdgeProcs[myRank];
      else {
        edgeWithoutNG[nEdgeWithoutNG++] = i;
        edgeWithoutNG[nEdgeWithoutNG++] = j;
      }
    }
  }
  nEdgeWithoutNG = nEdgeWithoutNG/2;

  /*
   * Define distributed hash table to look for
   * inter partition boundary edges (inter processus)
   */

  PDM_g_num_t *nKeyProcs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (lComm + 1));

  _n_g_enttities_compute(mesh);
  
  nKeyProcs[0] = 1;
  PDM_g_num_t gKeyMax = 2 * mesh->nGVtx;

  PDM_g_num_t _nKeyProc = gKeyMax / lComm;
  int nKeyProc = (int) _nKeyProc;

  for (int i = 0; i < lComm; i++) {
    nKeyProcs[i+1] = nKeyProc;
  }

  PDM_g_num_t _rest = gKeyMax % lComm;
  int rest = (int ) _rest;
  for (int i = 0; i < rest; i++)
    nKeyProcs[i+1] += 1;

  for (int i = 1; i < lComm + 1; i++)
    nKeyProcs[i] = nKeyProcs[i] + nKeyProcs[i-1];
  _nKeyProc = nKeyProcs[myRank+1] - nKeyProcs[myRank];
  nKeyProc = (int) _nKeyProc;

  int lDHashTableIdx  = nKeyProc + 1;
  int *dHashTableIdx  = PDM_array_zeros_int(lDHashTableIdx);
  int *dNHashTable    = PDM_array_zeros_int(nKeyProc);


  const int nDataToSend = 6;

  //FIXMe : AJOUTER UNE VALEUR LE NUMERO DE PROC

  int *edgeToSendN = PDM_array_zeros_int(lComm);


  for (int i = 0; i < nEdgeWithoutNG; i++) {
    int i_part = edgeWithoutNG[2*i];
    PDM_surf_part_t *part = mesh->part[i_part];
    int iEdge = edgeWithoutNG[2*i + 1];
    int vtx1 = part->edgeVtx[2*iEdge];
    int vtx2 = part->edgeVtx[2*iEdge + 1];
    PDM_g_num_t gVtx1 = part->vtx_ln_to_gn[vtx1 - 1];
    PDM_g_num_t gVtx2 = part->vtx_ln_to_gn[vtx2 - 1];

    PDM_g_num_t key = gVtx1 + gVtx2;
    edgeToSendN[PDM_binary_search_gap_long (key, nKeyProcs, lComm + 1)] += nDataToSend;
  }

  int *edgeToSendIdx = PDM_array_new_idx_from_sizes_int(edgeToSendN, lComm);
  PDM_array_reset_int(edgeToSendN, lComm, 0);

  /*
   * Store keys to send to the others processes
   */

  PDM_g_num_t *edgeToSend =
    (PDM_g_num_t *) malloc(edgeToSendIdx[lComm] * sizeof(PDM_g_num_t));

  for (int i = 0; i < nEdgeWithoutNG; i++) {
    int i_part = edgeWithoutNG[2*i];
    PDM_surf_part_t *part = mesh->part[i_part];
    int iEdge = edgeWithoutNG[2*i + 1];
    int vtx1 = part->edgeVtx[2*iEdge];
    int vtx2 = part->edgeVtx[2*iEdge + 1];
    PDM_g_num_t gVtx1 = part->vtx_ln_to_gn[vtx1 - 1];
    PDM_g_num_t gVtx2 = part->vtx_ln_to_gn[vtx2 - 1];

    PDM_g_num_t key = gVtx1 + gVtx2;

    int i_rank1 = PDM_binary_search_gap_long (key, nKeyProcs, lComm + 1);

    int idx1             = edgeToSendIdx[i_rank1] + edgeToSendN[i_rank1];
    edgeToSend[idx1  ]   = key;

    edgeToSend[idx1+1]   = gVtx1;
    edgeToSend[idx1+2]   = gVtx2;

    edgeToSend[idx1+3]   = i_part;
    edgeToSend[idx1+4]   = iEdge;
    edgeToSend[idx1+5]   = -1;
    edgeToSendN[i_rank1] += nDataToSend;
  }

  free(edgeWithoutNG);

  /*
   * Receive keys from the others processes
   */
  int *edgeToRecvN = (int *) malloc(lComm * sizeof(int));

  PDM_MPI_Alltoall(edgeToSendN,
               1,
               PDM_MPI_INT,
               edgeToRecvN,
               1,
               PDM_MPI_INT,
               mesh->comm);

  int *edgeToRecvIdx =  PDM_array_new_idx_from_sizes_int(edgeToRecvN, lComm);

  PDM_g_num_t *edgeToRecv =
    (PDM_g_num_t *) malloc(edgeToRecvIdx[lComm]*sizeof(PDM_g_num_t));

  PDM_MPI_Alltoallv(edgeToSend,
                edgeToSendN,
                edgeToSendIdx,
                PDM__PDM_MPI_G_NUM,
                edgeToRecv,
                edgeToRecvN,
                edgeToRecvIdx,
                PDM__PDM_MPI_G_NUM,
                mesh->comm);

  int nRecvKey = edgeToRecvIdx[lComm]/nDataToSend;

  /*
   * Fill hash table
   */

  const int nDatadHashTable = 2;

  PDM_g_num_t *dHashTable =
    (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * nDatadHashTable * nRecvKey);


  for (int i = 0; i < nRecvKey; i++) {
    PDM_g_num_t _keyLoc = edgeToRecv[i*nDataToSend] - nKeyProcs[myRank];
    int keyLoc       = (int) (_keyLoc);
    dHashTableIdx[keyLoc+1] += nDatadHashTable;
  }

  dHashTableIdx[0] = 0;
  PDM_array_accumulate_int(dHashTableIdx, nKeyProc + 1);

  int cptEdgeToRecv = 0;
  for (int i = 0; i < lComm; i++) {
    for (int j = edgeToRecvIdx[i]; j <  edgeToRecvIdx[i+1]; j += nDataToSend) {
      PDM_g_num_t _keyLoc   = edgeToRecv[j] - nKeyProcs[myRank];
      cptEdgeToRecv++;

      int keyLoc       = (int) _keyLoc;
      dHashTable[dHashTableIdx[keyLoc]+(dNHashTable[keyLoc]++)] =
        cptEdgeToRecv - 1;
	  dHashTable[dHashTableIdx[keyLoc]+(dNHashTable[keyLoc]++)] = i;
   }
  }

  free(dNHashTable);

  PDM_g_num_t *gNBoundPartEdge =
    (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 1 * nRecvKey);
  PDM_g_num_t gNCurrent = 0;

  /*
   * Give a global number to partition boundary edges
   */

  for (int i = 0; i < nKeyProc; i++) {
    int idx          = dHashTableIdx[i];
    int nEdgeSameSum = (dHashTableIdx[i+1] - idx);
    for (int j = idx; j < idx + nEdgeSameSum; j+=2) {
      int iEdge = (int) dHashTable[j];
      int iProc = (int) dHashTable[j+1];
      if (iEdge != -1) {
        int        fusion = 0;
        PDM_g_num_t gVtx1 = edgeToRecv[iEdge*nDataToSend + 1];
        PDM_g_num_t gVtx2 = edgeToRecv[iEdge*nDataToSend + 2];
        int        i_part = (int) edgeToRecv[iEdge*nDataToSend + 3];
        int        iEdgePart = (int) edgeToRecv[iEdge*nDataToSend + 4];

        for (int k = j + 2; k < idx + nEdgeSameSum; k+=2) {
          int        iEdge1 = (int) dHashTable[k];
          int        iProc1 = (int) dHashTable[k+1];

          if (iEdge1 != -1) {
            PDM_g_num_t gVtx11 = edgeToRecv[iEdge1*nDataToSend + 1];
            PDM_g_num_t gVtx12 = edgeToRecv[iEdge1*nDataToSend + 2];
            int        iPart1 = (int) edgeToRecv[iEdge1*nDataToSend + 3];
            int        iEdgePart1 = (int) edgeToRecv[iEdge1*nDataToSend + 4];

            if (((gVtx1 == gVtx11) && (gVtx2 == gVtx12)) ||
                ((gVtx2 == gVtx11) && (gVtx1 == gVtx12))) {

              gNBoundPartEdge[iEdge]  = ++gNCurrent;
              gNBoundPartEdge[iEdge1] =   gNCurrent;

              edgeToRecv[iEdge1*nDataToSend    ] = -1;
              edgeToRecv[iEdge1*nDataToSend + 1] = iPart1;
              edgeToRecv[iEdge1*nDataToSend + 2] = iEdgePart1;
              edgeToRecv[iEdge1*nDataToSend + 3] = i_part;
              edgeToRecv[iEdge1*nDataToSend + 4] = iEdgePart;
              edgeToRecv[iEdge1*nDataToSend + 5] = iProc;

              edgeToRecv[iEdge*nDataToSend    ] = -1;
              edgeToRecv[iEdge*nDataToSend + 1] = i_part;
              edgeToRecv[iEdge*nDataToSend + 2] = iEdgePart;
              edgeToRecv[iEdge*nDataToSend + 3] = iPart1;
              edgeToRecv[iEdge*nDataToSend + 4] = iEdgePart1;
              edgeToRecv[iEdge*nDataToSend + 5] = iProc1;

              dHashTable[k]   = -1;
              dHashTable[k+1] = -1;
              dHashTable[j]   = -1;
              dHashTable[j+1] = -1;
              fusion = 1;
              break;
            }
          }
        }
        if (!fusion) {
          edgeToRecv[iEdge*nDataToSend    ] = -1;
          edgeToRecv[iEdge*nDataToSend + 1] = i_part;
          edgeToRecv[iEdge*nDataToSend + 2] = iEdgePart;
          edgeToRecv[iEdge*nDataToSend + 3] = -1;
          edgeToRecv[iEdge*nDataToSend + 4] = -1;
          edgeToRecv[iEdge*nDataToSend + 5] = -1;
        }
      }
    }
  }

  /*
   * Give a global number to real boundary edges
   */

  for (int i = 0; i < nKeyProc; i++) {
  int idx          = dHashTableIdx[i];
    int nEdgeSameSum = (dHashTableIdx[i+1] - idx);
    for (int j = idx; j < idx + nEdgeSameSum; j+=2) {
      int iEdge = (int) dHashTable[j];
      if (iEdge != -1) {
        gNBoundPartEdge[iEdge] = ++gNCurrent;
      }
    }
  }

  /*
   * Update global numbering with taking into account o
   * of global numbering of internal edges
   */

  PDM_g_num_t *gNCurrentProcs =
    (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (lComm + 1));

  PDM_MPI_Allgather((void *) &gNCurrent, 1, PDM__PDM_MPI_G_NUM,
                (void *) (&gNCurrentProcs[1]), 1, PDM__PDM_MPI_G_NUM, mesh->comm);

  gNCurrentProcs[0] = 0;
  PDM_array_accumulate_gnum(gNCurrentProcs, lComm + 1);

  for (int i = 0; i < lComm + 1; i++)
    gNCurrentProcs[i] +=  nIntEdgeProcs[lComm];

  /*
   * Update global numbering with taking into account
   * the global numbering of internal edges
   */

  for (int i = 0; i < nRecvKey; i++) {
    edgeToRecv[i*nDataToSend    ] = gNBoundPartEdge[i] + gNCurrentProcs[myRank];
  }

  free(nIntEdgeProcs);
  free(dHashTable);
  free(dHashTableIdx);

  /*
   * Return to sender of results
   */

  PDM_MPI_Alltoallv(edgeToRecv,
                edgeToRecvN,
                edgeToRecvIdx,
                PDM__PDM_MPI_G_NUM,
                edgeToSend,
                edgeToSendN,
                edgeToSendIdx,
                PDM__PDM_MPI_G_NUM,
                mesh->comm);

  free(edgeToRecvN);
  free(edgeToRecvIdx);
  free(edgeToRecv);

  /*
   * Copy data already computed (intra processus boundary partition)
   */

  for (int i = 0; i < lComm; i++) {

    int id1 = edgeToSendIdx[i];
    int id2 = (edgeToSendIdx[i+1] - edgeToSendIdx[i]) / nDataToSend;

    int idx = id1;
    for (int j = 0; j < id2; j++) {

      PDM_g_num_t edgeGn     =       edgeToSend[idx++];
      int        i_part      = (int) edgeToSend[idx++];
      int        iEdgePart  = (int) edgeToSend[idx++];
      int        iPart1     = (int) edgeToSend[idx++];
      int        iEdgePart1 = (int) edgeToSend[idx++];
      int        iProc1     = (int) edgeToSend[idx++];
      PDM_surf_part_t *part =  mesh->part[i_part];

      part->edgeLnToGn[iEdgePart] = edgeGn;

      if (iPart1 != -1) {
        //        _part_t *part = mesh->part[i_part];
        PDM_part_bound_t *edgePartBound = part->edgePartBound;

        edgePartCur[i_part] += 1;
        PDM_part_bound_local_elt_set(edgePartBound,
                                     edgePartCur[i_part],
                                     iEdgePart+1);

        PDM_part_bound_offer_elt_set (edgePartBound,
                                      edgePartCur[i_part],
                                      0,
                                      PDM_ABS (part->edgeFace[2*iEdgePart]),
                                      part->face_ln_to_gn[PDM_ABS (part->edgeFace[2*iEdgePart])-1]);

        PDM_part_bound_distant_elt_set(edgePartBound,
                                       edgePartCur[i_part],
                                       0,
                                       iProc1,
                                       iPart1,
                                       iEdgePart1+1);
      }
    }
  }

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part =  mesh->part[i];
    PDM_part_bound_adjust_size(part->edgePartBound,
                               edgePartCur[i]);
  }

  PDM_g_num_t nGEdgeProc = 0;
  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    for (int j = 0; j < part->nEdge; j++) {
      nGEdgeProc = PDM_MAX (nGEdgeProc, part->edgeLnToGn[j]);
    }
  }

  PDM_MPI_Allreduce(&nGEdgeProc, &(mesh->nGEdge), 1,
                PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, mesh->comm);

  free(edgePartCur);
  free(nKeyProcs);
  free(edgeToSendN);
  free(edgeToSendIdx);
  free(edgeToSend);
  free(gNBoundPartEdge);
  free(gNCurrentProcs);

}


/**
 * \brief Compute inter partition communication graph between vertices
 *
 * This function computes edge global numbering
 *
 * \param [in]  mesh    mesh to compute global numbering
 *
 */

void
PDM_surf_mesh_build_vtx_part_bound
(
PDM_surf_mesh_t *mesh
)
{
  assert (mesh != NULL);

  int n_part = mesh->n_part;

  int myRank;
  int lComm;

  PDM_MPI_Comm_rank(mesh->comm, &myRank);
  PDM_MPI_Comm_size(mesh->comm, &lComm);

  /*
   * Define distributed hash table to look for
   * inter partition boundary vertices (inter processus)
   *
   */

  int *nKeyProcs = (int *) malloc(sizeof(int) * (lComm + 1));

  _n_g_enttities_compute(mesh);

  nKeyProcs[0] = 1;
  PDM_g_num_t gKeyMax = mesh->nGVtx + 1;

  PDM_g_num_t _nKeyProc = gKeyMax / lComm;
  int nKeyProc = (int) _nKeyProc;

  for (int i = 0; i < lComm; i++)
    nKeyProcs[i+1] = nKeyProc;

  PDM_g_num_t _rest = gKeyMax % lComm;
  int rest = (int) _rest;
  for (int i = 0; i < rest; i++)
    nKeyProcs[i+1] += 1;

  for (int i = 1; i < lComm + 1; i++)
    nKeyProcs[i] = nKeyProcs[i] + nKeyProcs[i-1];

  nKeyProc = nKeyProcs[myRank+1] - nKeyProcs[myRank];

  int lDHashTableIdx = nKeyProc + 1;
  int *dHashTableIdx = PDM_array_zeros_int(lDHashTableIdx);
  int *dNHashTable   = PDM_array_zeros_int(nKeyProc);

  const int nDataToSend = 3;
  int *vtx_to_sendN = PDM_array_zeros_int(lComm);

  int maxn_vtx = 0;
  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    maxn_vtx = PDM_MAX(maxn_vtx, part->n_vtx);
  }

  int *tagVtx = PDM_array_zeros_int(maxn_vtx);


  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    for (int j = 0; j < part->nEdge; j++) {
      if ((part->edgeFace[2*j] == 0) ||
          (part->edgeFace[2*j + 1] == 0)) {
        int vtx1 = part->edgeVtx[2*j] - 1;
        int vtx2 = part->edgeVtx[2*j + 1] -1;
        if (tagVtx[vtx1] == 0) {
          PDM_g_num_t _keyVtx = part->vtx_ln_to_gn[vtx1] % gKeyMax;
          int keyVtx = (int) _keyVtx;
          vtx_to_sendN[PDM_binary_search_gap_int (keyVtx, nKeyProcs, lComm + 1)] +=
            nDataToSend;
          tagVtx[vtx1] = 1;
        }

        if (tagVtx[vtx2] == 0) {
          PDM_g_num_t _keyVtx = part->vtx_ln_to_gn[vtx2] % gKeyMax;
          int keyVtx = (int) _keyVtx;
          vtx_to_sendN[PDM_binary_search_gap_int (keyVtx, nKeyProcs, lComm + 1)] +=
            nDataToSend;
          tagVtx[vtx2] = 1;
        }
      }
    }
    PDM_array_reset_int(tagVtx, maxn_vtx, 0);
  }

  int *vtx_to_sendIdx = PDM_array_new_idx_from_sizes_int(vtx_to_sendN, lComm);
  PDM_array_reset_int(vtx_to_sendN, lComm, 0);

  /*
   * Store keys to send to the others processes
   */

  PDM_g_num_t *vtx_to_send =
    (PDM_g_num_t *) malloc(vtx_to_sendIdx[lComm] * sizeof(PDM_g_num_t));

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    for (int j = 0; j < part->nEdge; j++) {
      if ((part->edgeFace[2*j    ] == 0) ||
          (part->edgeFace[2*j + 1] == 0)) {
        int vtx1 = part->edgeVtx[2*j] - 1;
        int vtx2 = part->edgeVtx[2*j + 1] -1;

        if (tagVtx[vtx1] == 0) {
          PDM_g_num_t _keyVtx = part->vtx_ln_to_gn[vtx1] % gKeyMax;
          int keyVtx = (int) _keyVtx;
          int i_rank  = PDM_binary_search_gap_int (keyVtx, nKeyProcs, lComm + 1);
          int idx    = vtx_to_sendIdx[i_rank] + vtx_to_sendN[i_rank];

          vtx_to_send[idx  ]   = keyVtx;
          vtx_to_send[idx+1]   = i;
          vtx_to_send[idx+2]   = vtx1;

          vtx_to_sendN[i_rank] += nDataToSend;
          tagVtx[vtx1]       = 1;
        }

        if (tagVtx[vtx2] == 0) {
          PDM_g_num_t _keyVtx = part->vtx_ln_to_gn[vtx2] % gKeyMax;
          int keyVtx = (int) _keyVtx;
          int i_rank  = PDM_binary_search_gap_int (keyVtx, nKeyProcs, lComm + 1);
          int idx    = vtx_to_sendIdx[i_rank] + vtx_to_sendN[i_rank];

          vtx_to_send[idx  ]   = keyVtx;
          vtx_to_send[idx+1]   = i;
          vtx_to_send[idx+2]   = vtx2;

          vtx_to_sendN[i_rank] += nDataToSend;
          tagVtx[vtx2]       = 1;
        }
      }
    }
    PDM_array_reset_int(tagVtx, maxn_vtx, 0);
  }
  free(tagVtx);

  /*
   * Receive keys from the others processes
   */

  int *vtxToRecvN = PDM_array_const_int(lComm, -100);


  PDM_MPI_Alltoall(vtx_to_sendN,
               1,
               PDM_MPI_INT,
               vtxToRecvN,
               1,
               PDM_MPI_INT,
               mesh->comm);

  int *vtxToRecvIdx =  PDM_array_new_idx_from_sizes_int(vtxToRecvN, lComm);

  PDM_g_num_t *vtxToRecv =
    (PDM_g_num_t *) malloc(vtxToRecvIdx[lComm]*sizeof(PDM_g_num_t));


  PDM_MPI_Alltoallv(vtx_to_send,
                vtx_to_sendN,
                vtx_to_sendIdx,
                PDM__PDM_MPI_G_NUM,
                vtxToRecv,
                vtxToRecvN,
                vtxToRecvIdx,
                PDM__PDM_MPI_G_NUM,
                mesh->comm);

  int n_vtxProc = vtxToRecvIdx[lComm]/nDataToSend;

  /*
   * Fill hash table
   */

  const int nDatadHashTable = 2;
  PDM_g_num_t *dHashTable =
    (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * nDatadHashTable * n_vtxProc);

  for (int i = 0; i < n_vtxProc; i++) {
    PDM_g_num_t key = vtxToRecv[i*nDataToSend] - nKeyProcs[myRank];

    int keyLoc       = (int) key;
    dHashTableIdx[keyLoc+1] += nDatadHashTable;
  }

  dHashTableIdx[0] = 0;
  PDM_array_accumulate_int(dHashTableIdx, nKeyProc+1);

  int cptVtxToRecv = 0;
  for (int i = 0; i < lComm; i++) {
    for (int j = vtxToRecvIdx[i]; j < vtxToRecvIdx[i+1]; j += nDataToSend) {
      PDM_g_num_t key = vtxToRecv[nDataToSend*cptVtxToRecv] - nKeyProcs[myRank];
      cptVtxToRecv++;

      int keyLoc       = (int) key;
      dHashTable[dHashTableIdx[keyLoc]+(dNHashTable[keyLoc]++)] =
        cptVtxToRecv - 1;
      dHashTable[dHashTableIdx[keyLoc]+(dNHashTable[keyLoc]++)] = i;
    }
  }

  free(dNHashTable);
  free(nKeyProcs);

  /*
   * Look for
   */

  const int nDataToSend1 = 3;
  const int nDataToSend2 = 3;

  PDM_array_reset_int(vtxToRecvN, lComm, 0);

  for (int i = 0; i < nKeyProc; i++) {
    int idx        = dHashTableIdx[i];
    int n_vtxSameGN = (dHashTableIdx[i+1] - idx);
    int nn = n_vtxSameGN/2 - 1;

    for (int j = 0; j < n_vtxSameGN; j+=2) {

      int iProc = (int) dHashTable[idx + j + 1];

      vtxToRecvN[iProc] += nDataToSend2 * nn + nDataToSend1;

    }
  }

  PDM_MPI_Alltoall(vtxToRecvN,
               1,
               PDM_MPI_INT,
               vtx_to_sendN,
               1,
               PDM_MPI_INT,
               mesh->comm);

  /*
   *  Backup vtxToRecv
   */

  PDM_g_num_t *_cpVtxToRecv =
    (PDM_g_num_t *) malloc(vtxToRecvIdx[lComm]*sizeof(PDM_g_num_t));
  memcpy(_cpVtxToRecv, vtxToRecv, vtxToRecvIdx[lComm]*sizeof(PDM_g_num_t));

  /*
   * Adjust Size of vtxToRecvN, vtxToRecvIdx and vtxToRecv,
   * vtx_to_sendN, vtx_to_sendIdx and vtx_to_send
   */

  vtxToRecvIdx[0] = 0;
  vtx_to_sendIdx[0] = 0;

  for (int i = 0; i < lComm; i++) {
    vtxToRecvIdx[i+1] = vtxToRecvIdx[i] + vtxToRecvN[i];
    vtx_to_sendIdx[i+1] = vtx_to_sendIdx[i] + vtx_to_sendN[i];
    vtxToRecvN[i] = 0;
  }

  vtxToRecv = (PDM_g_num_t *) realloc((void *) vtxToRecv,
                                     vtxToRecvIdx[lComm]*sizeof(PDM_g_num_t));
  vtx_to_send = (PDM_g_num_t *) realloc((void *) vtx_to_send,
                                     vtx_to_sendIdx[lComm]*sizeof(PDM_g_num_t));

  for (int i = 0; i < nKeyProc; i++) {
    int idx        = dHashTableIdx[i];
    int n_vtxSameGN = (dHashTableIdx[i+1] - idx);
    int nn        = n_vtxSameGN/2 - 1;

    for (int j = 0; j < n_vtxSameGN; j+=2) {

      int iVtx  = (int) dHashTable[idx + j    ];
      int iProc = (int) dHashTable[idx + j + 1];

      int        i_part     = (int) _cpVtxToRecv[iVtx*nDataToSend + 1];
      int        iVtxPart  = (int) _cpVtxToRecv[iVtx*nDataToSend + 2];

      vtxToRecv[vtxToRecvIdx[iProc]+vtxToRecvN[iProc]++] = i_part;
      vtxToRecv[vtxToRecvIdx[iProc]+vtxToRecvN[iProc]++] = iVtxPart;
      vtxToRecv[vtxToRecvIdx[iProc]+vtxToRecvN[iProc]++] = nn;

      for (int k = 0; k < n_vtxSameGN; k+=2) {
        if (k != j) {
          int iVtx1  = (int) dHashTable[idx + k    ];
          int iProc1 = (int) dHashTable[idx + k + 1];

          int        iPart1    = (int) _cpVtxToRecv[iVtx1*nDataToSend + 1];
          int        iVtxPart1 = (int) _cpVtxToRecv[iVtx1*nDataToSend + 2];

          vtxToRecv[vtxToRecvIdx[iProc]+vtxToRecvN[iProc]++] = iProc1;
          vtxToRecv[vtxToRecvIdx[iProc]+vtxToRecvN[iProc]++] = iPart1;
          vtxToRecv[vtxToRecvIdx[iProc]+vtxToRecvN[iProc]++] = iVtxPart1;
        }
      }
    }
  }

  free(dHashTable);
  free(dHashTableIdx);
  free(_cpVtxToRecv);

  PDM_MPI_Alltoallv(vtxToRecv,
                vtxToRecvN,
                vtxToRecvIdx,
                PDM__PDM_MPI_G_NUM,
                vtx_to_send,
                vtx_to_sendN,
                vtx_to_sendIdx,
                PDM__PDM_MPI_G_NUM,
                mesh->comm);

  free(vtxToRecv);
  free(vtxToRecvN);
  free(vtxToRecvIdx);

  int **nConnectedElt = (int **) malloc(sizeof(int *) * mesh->n_part);
  int **nOfferElt = (int **) malloc(sizeof(int *) * mesh->n_part);
  int *nEltPartBound = (int *) malloc(sizeof(int) * mesh->n_part);
  for (int i = 0; i < n_part; i++) {
    nConnectedElt[i] = NULL;
    nEltPartBound[i] = 0;
  }


  //FIXME: 3 prochaines Boucles a verifier lors de la validation !!!!!

  for (int i = 0; i < lComm; i++) {

    int id1 = vtx_to_sendIdx[i];
    int id2 = vtx_to_sendIdx[i+1];

    int j = id1;
    while (j < id2) {
      int i_part         = (int) vtx_to_send[j++];
      j += 1;
      int n_vtxConnected = (int) vtx_to_send[j++];


      if (n_vtxConnected > 0) {

        nEltPartBound[i_part] += 1;
        j += 3 * n_vtxConnected;
      }
    }
  }

  for (int i = 0; i < n_part; i++) {
    nConnectedElt[i] = (int *) malloc(sizeof(int) * nEltPartBound[i]);
    nOfferElt[i] = (int *) malloc(sizeof(int) * nEltPartBound[i]);
    nEltPartBound[i] = 0;
  }

  for (int i = 0; i < lComm; i++) {

    int id1 = vtx_to_sendIdx[i];
    int id2 = vtx_to_sendIdx[i+1];

    int j = id1;
    while (j < id2) {

      int i_part         = (int) vtx_to_send[j++];
      int iVtxPart      = (int) vtx_to_send[j++];
      int n_vtxConnected = (int) vtx_to_send[j++];

      if (n_vtxConnected > 0) {

        PDM_surf_part_t *part = mesh->part[i_part];
        nOfferElt[i_part][nEltPartBound[i_part]] =
          part->vtxEdgeIdx[iVtxPart+1] -
          part->vtxEdgeIdx[iVtxPart];

        nConnectedElt[i_part][nEltPartBound[i_part]++] = n_vtxConnected;

        j += 3 * n_vtxConnected;
      }
    }
  }

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    part->vtxPartBound = PDM_part_bound_create (lComm,
                                                part->n_vtx,
                                                nEltPartBound[i],
                                                PDM_PART_BOUND_CPLX,
                                                nConnectedElt[i],
                                                nOfferElt[i],
                                                mesh->nGEdge,
                                                part->nEdge,
                                                part->edgeLnToGn);

  }

  for (int i = 0; i < n_part; i++) {
    nEltPartBound[i] = 0;
    free (nConnectedElt[i]);
    free (nOfferElt[i]);
  }

  free (nOfferElt);
  free (nConnectedElt);

  for (int i = 0; i < lComm; i++) {

    int id1 = vtx_to_sendIdx[i];
    int id2 = vtx_to_sendIdx[i+1];

    int j = id1;
    while (j < id2) {

      int i_part         = (int) vtx_to_send[j++];
      int iVtxPart      = (int) vtx_to_send[j++];
      int n_vtxConnected = (int) vtx_to_send[j++];

      if (n_vtxConnected > 0) {
        PDM_surf_part_t *part = mesh->part[i_part];
        PDM_part_bound_t *part_bound = part->vtxPartBound;
        int idxEdge1 = part->vtxEdgeIdx[iVtxPart];
        int idxEdge2 = part->vtxEdgeIdx[iVtxPart+1];

        PDM_part_bound_local_elt_set (part_bound,
                                      ++nEltPartBound[i_part],
                                      iVtxPart+1);

        int nOffert = PDM_part_bound_n_offer_elt_get (part_bound,
                                                      nEltPartBound[i_part]);

        assert(nOffert == idxEdge2 - idxEdge1);
        for (int iOfferElt = 0; iOfferElt < nOffert; ++iOfferElt) {

          int lNum = part->vtxEdge[idxEdge1 + iOfferElt];
          PDM_g_num_t gNum = part->edgeLnToGn[lNum-1];

          PDM_part_bound_offer_elt_set (part_bound,
                                        nEltPartBound[i_part],
                                        iOfferElt,
                                        lNum,
                                        gNum);
        }

        for (int k = 0; k < n_vtxConnected; k++) {

          int iProc1        = (int) vtx_to_send[j++];
          int iPart1        = (int) vtx_to_send[j++];
          int iVtxPart1     = (int) vtx_to_send[j++];

          PDM_part_bound_distant_elt_set (part_bound,
                                          nEltPartBound[i_part],
                                          k,
                                          iProc1,
                                          iPart1,
                                          iVtxPart1+1);
        }
      }
    }
  }

  free(nEltPartBound);
  free(vtx_to_sendN);
  free(vtx_to_sendIdx);
  free(vtx_to_send);

}


/**
 * \brief Build ghost faces and edges
 *
 * This function computes ghost edges and ghost faces
 *
 * \param [in]  mesh    mesh to compute global numbering
 *
 */

void
PDM_surf_mesh_build_ghost_element
(
 PDM_surf_mesh_t *mesh
)
{
  assert (mesh != NULL);

  int n_part = mesh->n_part;

  int myRank;
  int lComm;

  PDM_MPI_Comm_rank(mesh->comm, &myRank);
  PDM_MPI_Comm_size(mesh->comm, &lComm);

  /*
   * Ghost faces
   */

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    PDM_graph_bound_t *graph_bound = mesh->interPartEdgeGraph;


    int nGhostFace = PDM_graph_bound_n_ghost_elt_get (graph_bound, i);
    part->nGhostFace = nGhostFace;

    part->nTotalFace = nGhostFace + part->n_face;

    /*
     * Update allocation
     */

    part->faceEdgeIdx = (int *) realloc (part->faceEdgeIdx,
                                         (part->nTotalFace + 1) * sizeof(int));

    int *ghostFaceEdgeIdx = part->faceEdgeIdx + part->n_face + 1;

    for (int ghostElt = 0; ghostElt < nGhostFace; ghostElt++) {
      int nTouchElt = PDM_graph_bound_ghost_elt_n_touch_elt_get (graph_bound,
                                                                 i,
                                                                 ghostElt);
      ghostFaceEdgeIdx[ghostElt] = nTouchElt;
    }

    for (int j = part->n_face; j < part->nTotalFace; j++) {
      part->faceEdgeIdx[j+1] += part->faceEdgeIdx[j];
    }

    part->faceEdge =
      (int *) realloc (part->faceEdge,
                       part->faceEdgeIdx[part->nTotalFace] * sizeof(int));

    /*
     * Update faceEdge
     */

    for (int ghostElt = 0; ghostElt < nGhostFace; ghostElt++) {
      int nTouchElt = PDM_graph_bound_ghost_elt_n_touch_elt_get (graph_bound,
                                                                 i,
                                                                 ghostElt);

      int *touchElt = PDM_graph_bound_ghost_elt_touch_elt_get (graph_bound,
                                                               i,
                                                               ghostElt);
      int deb = part->faceEdgeIdx[part->n_face];
      for (int iTouchElt = 0; iTouchElt < nTouchElt; iTouchElt++) {
        part->faceEdge[deb++] = touchElt[iTouchElt];
        part->edgeFace[2 * (touchElt[iTouchElt] - 1) + 1] =
          part->n_face + ghostElt + 1;
      }
    }
  }

  /*
   * Ghost edges
   */

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    PDM_graph_bound_t *graph_bound = mesh->interPartVtxGraph;

    int nGhostEdge = PDM_graph_bound_n_ghost_elt_get (graph_bound, i);
    part->nGhostEdge = nGhostEdge;

    part->nTotalEdge = nGhostEdge + part->nEdge;

    /*
     * Update allocation
     */

    part->edgeVtx =
      (int *) realloc (part->edgeVtx, 2 * part->nTotalEdge * sizeof(int));

    /*
     * Update edgeVtx
     */

    int *newVtxEdgeIdx = PDM_array_zeros_int(part->n_vtx + 1);

    for (int ghostElt = 0; ghostElt < nGhostEdge; ghostElt++) {
      int nTouchElt = PDM_graph_bound_ghost_elt_n_touch_elt_get (graph_bound,
                                                                 i,
                                                                 ghostElt);

      int *touchElt = PDM_graph_bound_ghost_elt_touch_elt_get (graph_bound,
                                                               i,
                                                               ghostElt);

      int deb = 2 * part->nEdge;
      for (int iTouchElt = 0; iTouchElt < nTouchElt; iTouchElt++) {
        int iVtx = touchElt[iTouchElt];
        part->edgeVtx[deb++] = iVtx;
        part->edgeVtx[deb++] = -1;
        newVtxEdgeIdx[iVtx]++;
      }
    }

    for (int i1 = 0; i1 < part->n_vtx; i1++) {
      newVtxEdgeIdx[i1+1] += part->vtxEdgeIdx[i1+1] - part->vtxEdgeIdx[i1];
    }

    for (int i1 = 0; i1 < part->n_vtx; i1++) {
      newVtxEdgeIdx[i1+1] += newVtxEdgeIdx[i1];
    }

    int *newVtxEdge = (int *) malloc (newVtxEdgeIdx[part->n_vtx] * sizeof(int));

    int *cptVtxEdge = PDM_array_zeros_int(part->n_vtx);

    for (int i1 = 0; i1 < part->n_vtx; i1++) {
      int idx = newVtxEdgeIdx[i1];
      for (int j = part->vtxEdgeIdx[i1]; j < part->vtxEdgeIdx[i1+1]; j++) {
        newVtxEdge[idx + cptVtxEdge[i1]++] = part->vtxEdge[j];
      }
    }

    free (part->vtxEdgeIdx);
    free (part->vtxEdge);

    part->vtxEdgeIdx = newVtxEdgeIdx;
    part->vtxEdge = newVtxEdge;

    for (int ghostElt = 0; ghostElt < nGhostEdge; ghostElt++) {
      int nTouchElt = PDM_graph_bound_ghost_elt_n_touch_elt_get (graph_bound,
                                                                 i,
                                                                 ghostElt);

      int *touchElt = PDM_graph_bound_ghost_elt_touch_elt_get (graph_bound,
                                                               i,
                                                               ghostElt);

      for (int iTouchElt = 0; iTouchElt < nTouchElt; iTouchElt++) {
        int iVtx = touchElt[iTouchElt] - 1;
        int idx = newVtxEdgeIdx[iVtx]+cptVtxEdge[iVtx]++;
        part->vtxEdge[idx] = ghostElt + part->nEdge + 1;
      }
    }

    free (cptVtxEdge);
  }
}


/**
 * \brief Build communication graph between internal partitions of any
 * initial mesh
 *
 * This function builds the communication graph between internal partitions
 * of each initial mesh
 *
 * \param [in]  mesh    mesh to compute global numbering
 *
 */

void
PDM_surf_mesh_build_exchange_graph
(
 PDM_surf_mesh_t *mesh
)
{

  assert (mesh != NULL);

  PDM_surf_mesh_build_vtx_part_bound (mesh);

  mesh->vtxPartBound =
    (PDM_part_bound_t **) malloc (mesh->n_part * sizeof(PDM_part_bound_t *));
  mesh->edgePartBound =
    (PDM_part_bound_t **) malloc (mesh->n_part * sizeof(PDM_part_bound_t *));

  for (int i = 0; i < mesh->n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    mesh->vtxPartBound[i]  = part->vtxPartBound;
    mesh->edgePartBound[i] = part->edgePartBound;
  }

  mesh->interPartVtxGraph = PDM_graph_bound_create (mesh->comm,
                                                    mesh->n_part,
                                                    mesh->vtxPartBound);

  mesh->interPartEdgeGraph = PDM_graph_bound_create (mesh->comm,
                                                     mesh->n_part,
                                                     mesh->edgePartBound);

  PDM_surf_mesh_build_ghost_element (mesh);

}


/**
 * \brief Build ghost faces and ghost edges
 *
 * This function builds ghost faces and ghost edges
 * of each initial mesh
 *
 * \param [in]  mesh        Mesh
 *
 */

void
PDM_surf_mesh_compute_carLgthVtx
(
PDM_surf_mesh_t *mesh
)
{
  assert (mesh != NULL);

  const int n_part = mesh->n_part;
  double **lEdge = (double **) malloc (sizeof(double *) * n_part);
  double **lEdgeGhost = (double **) malloc (sizeof(double *) * n_part);

  /*
   * Compute edge length
   */

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    const int nEdge = part->nEdge;
    const int nTotalEdge = part->nTotalEdge;

    if (part->vtxEdgeIdx == NULL) {
      PDM_error(__FILE__, __LINE__, 0,
              "Error _compute_carLgthVtx_mesh :"
              " Edges are not computed\n");
      abort();
    }

    if (part->carLgthVtx != NULL){
      PDM_error(__FILE__, __LINE__, 0,
              "Error _compute_carLgthVtx_mesh :"
              " Caracteristic length already computed\n");
      abort();
    }

    lEdge[i] = (double *) malloc (sizeof(double) * nTotalEdge);
    lEdgeGhost[i] = lEdge[i] + nEdge;

    for (int j = 0; j < nEdge; j++) {
      int vtx1 = part->edgeVtx[2*j    ] - 1;
      int vtx2 = part->edgeVtx[2*j + 1] - 1;

      const double *coordVtx1 = part->coords + 3 * vtx1;
      const double *coordVtx2 = part->coords + 3 * vtx2;

      double vEdge[3];
      for (int j1 = 0; j1 < 3; j1++) {
        vEdge[j1] = coordVtx2[j1] - coordVtx1[j1];
      }
      lEdge[i][j] = _MODULE (vEdge);
    }
  }

  /*
   * Initialize exchange to ghost egdes
   */

  PDM_graph_bound_exch_data_init (mesh->interPartVtxGraph,
                                  1,
                                  PDM_DOUBLE,
                                  (void **) lEdge,
                                  (void **) lEdgeGhost);

  /*
   * TODO : Vertex renumbering to store boundary vertices behind internal vertices
   *        To cover exchange by computation
   */

  /*
   * Finalize exchange to ghost egdes
   */

  PDM_graph_bound_exch_data_wait (mesh->interPartVtxGraph);

  free (lEdgeGhost);

  /*
   * Compute caracteristic length with local contribution
   */

  for (int i = 0; i < n_part; i++) {

    PDM_surf_part_t *part = mesh->part[i];
    const int n_vtx = part->n_vtx;
    double *_lEdge = lEdge[i];

    part->carLgthVtx = (double *) malloc (sizeof(double) * n_vtx);
    for (int j = 0; j < n_vtx; j++) {
      part->carLgthVtx[j] = DBL_MAX;
    }

    /*
     * TODO : Optimization : Split this boucle by blocks (vectorization)
     */

    for (int j = 0; j < n_vtx; j++) {
      for (int k = part->vtxEdgeIdx[j]; k < part->vtxEdgeIdx[j+1]; k++) {
        int edge = part->vtxEdge[k] - 1;
        part->carLgthVtx[j] = _MIN (part->carLgthVtx[j], _lEdge[edge]);
      }
      mesh->gMinCarLgthVtx = _MIN (mesh->gMinCarLgthVtx, part->carLgthVtx[j]);
      mesh->gMaxCarLgthVtx = _MAX (mesh->gMaxCarLgthVtx, part->carLgthVtx[j]);
    }

  }

  double gMin;
  double gMax;

  PDM_MPI_Allreduce(&mesh->gMinCarLgthVtx, &gMin, 1,
                PDM__PDM_MPI_REAL, PDM_MPI_MIN, mesh->comm);

  PDM_MPI_Allreduce(&mesh->gMaxCarLgthVtx, &gMax, 1,
                PDM__PDM_MPI_REAL, PDM_MPI_MIN, mesh->comm);

  mesh->gMinCarLgthVtx = gMin;
  mesh->gMaxCarLgthVtx = gMax;

  /*
   * Clean up
   */

  for (int i = 0; i < n_part; i++) {
    free (lEdge[i]);
  }
  free (lEdge);
}


/**
 * \brief Return face normal
 *
 * This function returns face normal after computation
 *
 * \param [in]  mesh       mesh
 *
 */

const double *
PDM_surf_mesh_face_normal_get
(
 PDM_surf_mesh_t  *mesh,
 int              i_part
)
{
  assert (mesh != NULL);

  int n_part = mesh->n_part;

  if ((i_part < 0) || (i_part > n_part)) {
    PDM_error(__FILE__, __LINE__, 0, "Error _part_face_normal_get : undefined part\n");
    abort();
  }
  PDM_surf_part_t *part = mesh->part[i_part];

  int n_face = part->n_face;
  const double *coords = part->coords;
  const int *face_vtx_idx = part->face_vtx_idx;
  const int *face_vtx = part->face_vtx;

  if (part->faceNormal  == NULL) {
    part->faceNormal = (double *) malloc (sizeof(double) * 3 * n_face);

    //TODO : vectorisation par paquet

    double *_faceNormal = part->faceNormal;
    for (int j = 0; j < n_face; j++) {
      double v1[3];
      double v2[3];
      double subNorm[3];
      for (int k = 0; k < 3; k++) {
        _faceNormal[k] = 0;
      }
      int n_vtx = face_vtx_idx[j+1] - face_vtx_idx[j];
      int ideb = face_vtx_idx[j];

      int idx0 = face_vtx[ideb] - 1;
      for (int k = 1; k < n_vtx - 1; k++) {
        int idx1 = face_vtx[ideb + k] - 1;
        int idx2 = face_vtx[ideb + (k+1) % n_vtx] - 1;
        for (int k1 = 0; k1 < 3; k1++) {
          v1[k1] = coords[3*idx1+k1] - coords[3*idx0+k1];
          v2[k1] = coords[3*idx2+k1] - coords[3*idx0+k1];
        }
        _CROSS_PRODUCT_3D (subNorm, v1, v2);
        for (int k1 = 0; k1 < 3; k1++) {
          _faceNormal[k1] += subNorm[k1];
        }
      }
        _faceNormal += 3;
    }
  }
  return part->faceNormal;
}


/**
 * \brief Chek if mesh is a plane surfece
 *
 * This function cheks if the mesh is a plane surface
 * and returns plane equation ant vertices barycenter
 *
 * \param [in]  ol             Overlay object
 * \param [in]  tolerance      Tolerance to accept surface as plane
 * \param [out] planeEquation  Plane equation
 * \param [out] barycenter     Vertices barycenter
 *
 */

int
PDM_surf_mesh_is_plane_surface
(
 PDM_surf_mesh_t  *mesh,
 double            tolerance,
 double            planeEquation[4],
 double            barycenter[3]
)
{
  assert (mesh != NULL);

  int myRank;
  PDM_MPI_Comm_rank (mesh->comm, &myRank);

  for (int k = 0; k < 4; k++) {
    planeEquation[k] = 0;
  }
  int isPlane = 1;
  int n_part = mesh->n_part;
  const double epsilonAbs = 1e-15;  // Absolute epsilon to protect / 0

  /*
   * Compute normal surface (not oriented)
   */

  double normalSum[3] = {0, 0, 0};

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    const int n_face = part->n_face;
    const double* faceNormal = PDM_surf_mesh_face_normal_get (mesh, i);

    for (int j = 0; j < n_face; j++) {

      if (_DOT_PRODUCT(faceNormal, normalSum) > 0) {
        for (int k = 0; k < 3; k++) {
          normalSum[k] += faceNormal[3*j+k];
        }
      }
      else {
        for (int k = 0; k < 3; k++) {
          normalSum[k] += -faceNormal[3*j+k];
        }
      }
    }
  }

  /*
   * Normalize
   */

  double _mod = _MAX(_MODULE (normalSum), epsilonAbs);
  for (int k = 0; k < 3; k++) {
    normalSum[k] = normalSum[k]/_mod;
  }

  /*
   * Check if surface is plane
   */

  for (int i = 0; i < n_part; i++) {
    PDM_surf_part_t *part = mesh->part[i];
    const int n_face = part->n_face;
    const double* faceNormal = PDM_surf_mesh_face_normal_get (mesh, i);

    const double* _faceNormal = faceNormal;
    for (int j = 0; j < n_face; j++) {
      double faceNormalNorm[3];
      double _mod1 = _MAX (_MODULE (_faceNormal), epsilonAbs);
      for (int k = 0; k < 3; k++) {
        faceNormalNorm[k] = _faceNormal[k]/_mod1;
      }
      if ((1 - fabs (_DOT_PRODUCT(faceNormalNorm, normalSum))) > tolerance) {
        isPlane = 0;
        break;
      }
      _faceNormal += 3;
    }
  }

  int isPlaneS;

  PDM_MPI_Allreduce (&isPlane, &isPlaneS, 1,
                 PDM_MPI_INT, PDM_MPI_SUM, mesh->comm);

  isPlane = isPlaneS;

  int lComm;
  PDM_MPI_Comm_size (mesh->comm, &lComm);

  /*
   * Compare plane equations between all ranks
   */

  if (isPlane == lComm) {

    isPlane = 1;
    double *normalSums = (double *) malloc (sizeof(double) * 3 * lComm);

    PDM_MPI_Allgather(normalSum, 3, PDM__PDM_MPI_REAL,
                  normalSums, 3, PDM__PDM_MPI_REAL,
                  mesh->comm);

    double *_normalSumRef = normalSums;
    int idx = 0;
    while (   (fabs(_normalSumRef[0]) < epsilonAbs)
           && (fabs(_normalSumRef[1]) < epsilonAbs)
           && (fabs(_normalSumRef[2]) < epsilonAbs)) {
      _normalSumRef += 3;
      idx++;
    }

    for (int k = 0; k < 3; k++) {
      planeEquation[k] = _normalSumRef[k];
    }

    for (int j = idx+1; j < lComm; j++) {
      double *_normalSumRank =  normalSums + 3 * j;

      /*
       * Pairwise comparaison
       */

      if (   (fabs(_normalSumRank[0]) < epsilonAbs)
          && (fabs(_normalSumRank[1]) < epsilonAbs)
          && (fabs(_normalSumRank[2]) < epsilonAbs)) {
        continue;
      }

      if (  (1 - fabs (_DOT_PRODUCT (_normalSumRef, _normalSumRank)))
          > tolerance) {
        isPlane = 0;
        break;
      }

    }
    free (normalSums);
  }
  else {
    isPlane = 0;
  }

  double center[3] = {0, 0, 0};


  if (isPlane) {

    PDM_g_num_t **_gnums = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);
    int *_n_elts = (int *) malloc (sizeof(int) * n_part);
    double **_coords = (double **) malloc (sizeof(double *) * n_part);

    for (int i = 0; i < n_part; i++) {
      PDM_surf_part_t *part = mesh->part[i];
      const int n_vtx = part->n_vtx;
      const double* coords = part->coords;
      const PDM_g_num_t* gnum = part->vtx_ln_to_gn;

      _gnums[i] = (PDM_g_num_t *) gnum;
      _coords[i] = (double *) coords;
      _n_elts[i] = n_vtx;

    }

    PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                         1,
                                                         _gnums,
                                                         NULL,
                                                         _n_elts,
                                                         n_part,
                                                         mesh->comm);

    double *_dcoords;
    int * block_stride;

    PDM_part_to_block_exch (ptb,
                            sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            3,
                            NULL,
                            (void **) _coords,
                            &block_stride,
                            (void **) &_dcoords);


    int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb);

    PDM_part_to_block_free (ptb);

    free (_gnums);
    free (_n_elts);
    free (_coords);

    for (int j = 0; j < n_elt_block; j++) {
      for (int k = 0; k < 3; k++) {
        center[k] += _dcoords[3*j+k];
      }
    }

    free (_dcoords);

    PDM_MPI_Allreduce (center, barycenter, 3,
                       PDM__PDM_MPI_REAL, PDM_MPI_SUM, mesh->comm);

    _n_g_enttities_compute(mesh);

    for (int k = 0; k < 3; k++) {
      barycenter[k] = barycenter[k] / mesh->nGVtx;
    }

    planeEquation[3] = - (  barycenter[0] * planeEquation[0]
                          + barycenter[1] * planeEquation[1]
                          + barycenter[2] * planeEquation[2]);

  }

  else {
    for (int k = 0; k < 4; k++) {
      planeEquation[k] = 0;
    }
  }

  return isPlane;
}


/**
 * \brief Compute face extents
 *
 * This function computes face extents
 *
 * \param [in]  ol       overlay object
 *
 */

void
PDM_surf_mesh_compute_faceExtentsMesh
(
 PDM_surf_mesh_t *mesh,
 double           tolerance
)
{
  assert (mesh != NULL);

  for (int i = 0; i < mesh->n_part; i++) {

    PDM_surf_part_t *part = mesh->part[i];
    const int n_face = part->n_face;
    const int *face_vtx = part->face_vtx;
    const int *face_vtx_idx = part->face_vtx_idx;
    const double *coords = part->coords;

    part->extents = (double *) malloc (sizeof(double) * 6 * n_face);

    /*
     * TODO : Optimization : Split this boucle by blocks (vectorization)
     */

    double *_extents = part->extents;
    for (int j = 0; j < n_face; j++) {

      for (int k1 = 0; k1 < 3; k1++) {
        _extents[k1]   = DBL_MAX;
        _extents[3+k1] = -DBL_MAX;
      }

      for (int k = face_vtx_idx[j]; k < face_vtx_idx[j+1]; k++) {
        int iVtx = face_vtx[k] - 1;
        double *_coords = (double *) coords + 3 * iVtx;

        for (int k1 = 0; k1 < 3; k1++) {
          _extents[k1]   = _MIN (_coords[k1], _extents[k1]);
          _extents[3+k1] = _MAX (_coords[k1], _extents[3+k1]);
        }

      }

      double delta = -DBL_MAX;

      for (int k1 = 0; k1 < 3; k1++) {
        delta = _MAX (delta, fabs (_extents[k1+3] - _extents[k1]));
      }

      delta *= tolerance;

      for (int k1 = 0; k1 < 3; k1++) {
        _extents[k1]   +=  - delta;
        _extents[3+k1] +=    delta;
      }

      _extents += 6;
    }
  }
}


/**
 * \brief Compute face extents
 *
 * This function computes face extents
 *
 * \param [in]  mesh       Mesh object
 *
 */

void
PDM_surf_mesh_build_edges
(
 PDM_surf_mesh_t *mesh
)
{

  assert (mesh != NULL);

  for (int i = 0; i < mesh->n_part; i++)
    PDM_surf_part_build_edges (mesh->part[i]);

}

/**
 * \brief Return global minimum of caracteristic length vertex
 *
 * This function returns global minimum of caracteristic length vertex
 *
 * \param [in]  mesh       Mesh object
 *
 */

double
PDM_surf_mesh_gMinCarLgthVtx_get
(
 PDM_surf_mesh_t *mesh
)
{
  assert (mesh != NULL);

  return mesh->gMinCarLgthVtx;

}


/**
 * \brief Input a partition
 *
 * This function inputs a partition
 *
 * \param [in]  mesh       Mesh object
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
PDM_surf_mesh_part_input
(
 PDM_surf_mesh_t      *mesh,
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
  assert (mesh != NULL);

  mesh->part[i_part] = PDM_surf_part_create (n_face,
                                             face_vtx_idx,
                                             face_vtx,
                                             face_ln_to_gn,
                                             n_vtx,
                                             coords,
                                             vtx_ln_to_gn);

}


/**
 * \brief Return number of partitions
 *
 * \param [in]  mesh       Mesh object
 *
 * \return    Number of partitions
 */

int
PDM_surf_mesh_n_part_get
(
 PDM_surf_mesh_t      *mesh
)
{
  assert (mesh != NULL);

  return mesh->n_part;
}


/**
 * \brief Return number of faces
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Number of faces
 */

int
PDM_surf_mesh_part_n_face_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->n_face;
}


/**
 * \brief Return number of vertices
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Number of faces
 */

int
PDM_surf_mesh_part_n_vtx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->n_vtx;
}


/**
 * \brief Return extents for any face
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Extents
 */

const double *
PDM_surf_mesh_part_extents_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);


  PDM_surf_part_t *part  =  mesh->part[i_part];

  return part->extents;
}


/**
 * \brief Return face global number
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Extents
 */

const PDM_g_num_t *
PDM_surf_mesh_part_face_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->face_ln_to_gn;
}

/**
 * \brief Return Vertex global number
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Vertex global number
 */

const PDM_g_num_t *
PDM_surf_mesh_part_vtx_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->vtx_ln_to_gn;
}


/**
 * \brief Return Edge global number
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return  Edge global number
 */

const PDM_g_num_t *
PDM_surf_mesh_part_edge_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->edgeLnToGn;
}


/**
 * \brief Return Face to edge connectivity
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Face to edge connectivity
 */

const int *
PDM_surf_mesh_part_face_edge_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->faceEdge;
}


/**
 * \brief Return Face to vertex connectivity
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Face to vertex connectivity
 */

const int *
PDM_surf_mesh_part_face_vtx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->face_vtx;
}


/**
 * \brief Return Face to vertex connectivity index
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Face to vertex connectivity index
 */

const int *
PDM_surf_mesh_part_face_vtx_idx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  =  mesh->part[i_part];

  return part->face_vtx_idx;
}

/**
 * \brief Return Face to edge connectivity index
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Face to edge connectivity index
 */

const int *
PDM_surf_mesh_part_face_edge_idx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->faceEdgeIdx;
}



/**
 * \brief Return vertex coordinates
 *
 * \param [in]  mesh       Mesh object
 * \param [in]  i_part      Part number
 *
 * \return    Vertex coordinates
 */

const double *
PDM_surf_mesh_part_vtx_get
(
 PDM_surf_mesh_t *mesh,
 int              i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->coords;
}


/**
 * \brief Return Vertex caracteristic length
 *
 * This function returns global minimum of caracteristic length of vertex
 *
 * \param [in]  mesh       Mesh object
 *
 * \return  Vertex caracteristic length
 */

double *
PDM_surf_mesh_part_carLgthVtx_get
(
 PDM_surf_mesh_t *mesh,
 int              i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->carLgthVtx;
}

/**
 * \brief Return global number of edges
 *
 * \param [in]  mesh       Mesh object
 *
 * \return   Global number of edges
 */

PDM_g_num_t
PDM_surf_mesh_n_g_edge_get
(
 PDM_surf_mesh_t      *mesh
)
{
  assert (mesh != NULL);

  return mesh->nGEdge;
}


/**
 * \brief Return global number of vertices
 *
 * \param [in]  mesh       Mesh object
 *
 * \return   Global number of vertices
 */

PDM_g_num_t
PDM_surf_mesh_n_g_vtx_get
(
 PDM_surf_mesh_t      *mesh
)
{
  assert (mesh != NULL);
  
  _n_g_enttities_compute(mesh);

  return mesh->nGVtx;
}


/**
 * \brief Return global number of faces
 *
 * \param [in]  mesh       Mesh object
 *
 * \return   Global number of faces
 */

PDM_g_num_t
PDM_surf_mesh_n_g_face_get
(
 PDM_surf_mesh_t      *mesh
)
{
  assert (mesh != NULL);

  _n_g_enttities_compute(mesh);

  return mesh->nGFace;
}


#undef _DOT_PRODUCT
#undef _MODULE
#undef _MIN
#undef _MAX
#undef _CROSS_PRODUCT_3D

#ifdef __cplusplus
}
#endif /* __cplusplus */
