/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_surf_part.h"
#include "pdm_surf_part_priv.h"
#include "pdm_part_bound.h"
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

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_surf_part_t structure
 *
 * This function returns an initialized \ref PDM_surf_part_t structure
 *
 * \param [in]  n_face       Number of faces
 * \param [in]  face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]  face_vtx_idx  face -> vertex connectivity
 * \param [in]  face_ln_to_gn  Local face numbering to global face numbering
 * \param [in]  n_vtx        Number of vertices
 * \param [in]  coords      Coordinates
 * \param [in]  vtx_ln_to_gn   Local vertex numbering to global vertex numbering
 *
 * \return      A new initialized \ref PDM_surf_part_t structure
 *
 */

PDM_surf_part_t *
PDM_surf_part_create
(
const int         n_face,
const int        *face_vtx_idx,
const int        *face_vtx,
const PDM_g_num_t *face_ln_to_gn,
const int         n_vtx,
const double     *coords,
const PDM_g_num_t *vtx_ln_to_gn
)
{
  PDM_surf_part_t *_part = (PDM_surf_part_t *) malloc(sizeof(PDM_surf_part_t));

  _part->n_face                = n_face;
  _part->nGhostFace           = 0;
  _part->nTotalFace           = n_face;
  _part->sface_vtx             = face_vtx_idx[n_face];
  _part->face_vtx_idx           = face_vtx_idx;
  _part->face_vtx              = face_vtx;
  _part->faceEdgeIdx          = NULL;
  _part->faceEdge             = NULL;
  _part->face_ln_to_gn           = face_ln_to_gn;
  _part->n_vtx                 = n_vtx;
  _part->coords               = coords;
  _part->vtxEdgeIdx           = NULL;
  _part->vtxEdge              = NULL;
  _part->vtx_ln_to_gn            = vtx_ln_to_gn;
  _part->nEdge                = -1;
  _part->nGhostEdge           = -1;
  _part->nTotalEdge           = -1;
  _part->edgeFace             = NULL;
  _part->edgeVtx              = NULL;
  _part->edgeLnToGn           = NULL;
  _part->edgePartBound        = NULL;
  _part->vtxPartBound         = NULL;
  _part->carLgthVtx           = NULL;
  _part->faceNormal           = NULL;
  _part->extents              = NULL;

  return _part;
}


/**
 * \brief Delete a \ref PDM_surf_part_t structure
 *
 * This function deletes a  PDM_surf_part_t structure
 *
 * \param [in]  part      part to delete
 *
 * \return     Null pointer
 */

PDM_surf_part_t *
PDM_surf_part_free
(
 PDM_surf_part_t * part
)
{
  assert (part != NULL);

  if (part != NULL) {
    part->face_vtx_idx = NULL;
    part->face_vtx = NULL;
    if (part->faceEdgeIdx != NULL)
      free(part->faceEdgeIdx);
    if (part->faceEdge != NULL)
      free(part->faceEdge);

    part->face_ln_to_gn = NULL;
    part->coords = NULL;
    if (part->vtxEdgeIdx != NULL)
      free(part->vtxEdgeIdx);
    if (part->vtxEdge != NULL)
      free(part->vtxEdge);
    part->vtx_ln_to_gn = NULL;

    if (part->edgeFace != NULL)
      free(part->edgeFace);

    if (part->edgeVtx != NULL)
      free(part->edgeVtx);

    if (part->edgeLnToGn != NULL)
      free(part->edgeLnToGn);

    if (part->edgePartBound != NULL)
      part->edgePartBound = PDM_part_bound_free(part->edgePartBound);
    if (part->vtxPartBound != NULL)
      part->vtxPartBound = PDM_part_bound_free(part->vtxPartBound);

    if (part->carLgthVtx != NULL)
      free (part->carLgthVtx);

    if (part->faceNormal != NULL)
      free (part->faceNormal);

    if (part->extents != NULL)
      free (part->extents);

    free(part);
  }

  return NULL;
}


/**
 * \brief Compute partition edge entities
 *
 * This function defines edges of an initial partitiob and
 * computes edge connectivities
 *
 * \param [in]  _part      Partition to compute
 *
 */

void
PDM_surf_part_build_edges
(
PDM_surf_part_t *part
)
{

  assert (part != NULL);

  /*
   * build hash table with key = vertices sum and list of edges
   */

  int  lHashTableIdx = 2 * part->n_vtx + 1;
  int *hashTableIdx  = PDM_array_zeros_int(lHashTableIdx);
  int *hashTable     = PDM_array_const_int(part->face_vtx_idx[part->n_face], -1);
  int *nHashTable    = PDM_array_zeros_int(2 * part->n_vtx);

  int l_edges = 0;
  for (int i = 0; i < part->n_face; i++) {
    int idxFace = part->face_vtx_idx[i];
    for (int j = idxFace; j < part->face_vtx_idx[i+1]; j++) {
      int k = (j == part->face_vtx_idx[i+1] - 1) ? idxFace : (j + 1);
      int s_vtx = part->face_vtx[j] + part->face_vtx[k];
      hashTableIdx[s_vtx+1] += 1;
      l_edges += 2;
    }
  }

  hashTableIdx[0] = 0;
  PDM_array_accumulate_int(hashTableIdx, lHashTableIdx);

  int *listEdges = (int *) malloc(sizeof(int) * l_edges);
  int *edgeFaceUncompress = (int *) malloc(sizeof(int) * l_edges/2);
  l_edges = 0;

  int max_n_vtx_face = 0;
  for (int i = 0; i < part->n_face; i++) {
    int idxFace = part->face_vtx_idx[i];
    int n_vtx_face = part->face_vtx_idx[i+1] - idxFace;
    max_n_vtx_face = PDM_MAX(max_n_vtx_face, n_vtx_face);
    for (int j = idxFace; j < idxFace + n_vtx_face; j++) {
      int k = (j == part->face_vtx_idx[i+1] - 1) ? idxFace : (j + 1);
      int vtx1 = part->face_vtx[j];
      int vtx2 = part->face_vtx[k];
      int s_vtx = vtx1 + vtx2;

      hashTable[hashTableIdx[s_vtx] + (nHashTable[s_vtx]++)] = l_edges/2;
      edgeFaceUncompress[l_edges/2] = i;

      listEdges[l_edges++] = vtx1;
      listEdges[l_edges++] = vtx2;
    }
  }

  free(nHashTable);

  /*
   * Compress edges
   */

  int nEdgesFaces = l_edges/2;
  int *edgesToCompressEdges = PDM_array_const_int(nEdgesFaces, -1);

  /*
   * First allocation to l_edges
   */

  part->edgeVtx  = (int *) malloc(sizeof(int) * l_edges);
  part->nEdge    = 0;

  /*
   * Compute real edges
   */

  int nEdge = 0;

  part->vtxEdgeIdx  = PDM_array_zeros_int(part->n_vtx + 1);

  for (int i = 0; i < 2 * part->n_vtx; i++) {
    int idx          = hashTableIdx[i];
    int nEdgeSameSum = (hashTableIdx[i+1] - idx);
    for (int j = idx; j < idx + nEdgeSameSum; j++) {
      int iEdge = hashTable[j];
      if (iEdge != -1) {
        edgesToCompressEdges[iEdge] = nEdge;
        int vtx1 = listEdges[2*iEdge];
        int vtx2 = listEdges[2*iEdge+1];
        part->edgeVtx[2*nEdge]      = vtx1;
        part->edgeVtx[2*nEdge + 1]  = vtx2;
        part->vtxEdgeIdx[vtx1] += 1;
        part->vtxEdgeIdx[vtx2] += 1;
        for (int k = j + 1; k < idx + nEdgeSameSum; k++) {
          int iEdge1 = hashTable[k];
          if (iEdge1 != -1) {
            int vtx11 = listEdges[2*iEdge1];
            int vtx12 = listEdges[2*iEdge1+1];

            if (((vtx1 == vtx11) && (vtx2 == vtx12)) ||
                ((vtx2 == vtx11) && (vtx1 == vtx12))) {
              hashTable[k] = -1;

              edgesToCompressEdges[iEdge1] = nEdge;
            }
          }
        }
        nEdge += 1;
      }
    }
  }

  part->nEdge = nEdge;

  for (int i = 0; i < part->n_vtx; i++) {
    part->vtxEdgeIdx[i+1] = part->vtxEdgeIdx[i+1] + part->vtxEdgeIdx[i];
  }
  part->vtxEdge  = (int *) malloc(sizeof(int) * part->vtxEdgeIdx[part->n_vtx]);
  int *n_vtxEdge  = PDM_array_zeros_int(part->n_vtx);

  for (int i = 0; i < part->nEdge; i++) {
    int vtx1 = part->edgeVtx[2*i    ] - 1;
    int vtx2 = part->edgeVtx[2*i + 1] - 1;
    part->vtxEdge[part->vtxEdgeIdx[vtx1] + n_vtxEdge[vtx1]++] = i+1;
    part->vtxEdge[part->vtxEdgeIdx[vtx2] + n_vtxEdge[vtx2]++] = i+1;
  }
  free(n_vtxEdge);

  free(hashTable);
  free(hashTableIdx);
  free(listEdges);

  /*
   * Re-allocation to real size
   */

  part->edgeVtx = (int *) realloc(part->edgeVtx, sizeof(int) * 2 * nEdge);

  /*
   * edge -> face connectivity
   */

  part->edgeFace = PDM_array_zeros_int(2 * nEdge);

  for (int i = 0; i < nEdgesFaces; i++) {
    int iEdge = edgesToCompressEdges[i];
    int ifac1 = 2*iEdge;
    int ifac2 = 2*iEdge + 1;

    if (part->edgeFace[ifac1] == 0) {
      part->edgeFace[ifac1] = edgeFaceUncompress[i] + 1;
    }
    else if (part->edgeFace[ifac2] == 0) {
      part->edgeFace[ifac2] = edgeFaceUncompress[i] + 1;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Error _build_edges_part : Error in edgeFace computing\n");
      abort();
    }
  }

  if (0 == 1) {
    PDM_printf("edgeface : PDM_surf_part_build_edges\n");
    for (int i = 0; i < part->nEdge; i++) {
      PDM_printf("%d: %d %d\n", i, part->edgeFace[2*i], part->edgeFace[2*i+1]);
    }
  }

  free(edgesToCompressEdges);

  /*
   * face -> edge connectivity
   */

  part->faceEdgeIdx = (int *) malloc(sizeof(int) * (part->n_face + 1));
  memcpy(part->faceEdgeIdx, part->face_vtx_idx, sizeof(int) * (part->n_face + 1));
  const int *faceEdgeIdx = part->faceEdgeIdx;
  part->faceEdge = (int *) malloc(sizeof(int) * faceEdgeIdx[part->n_face]);
  int *n_faceEdge = PDM_array_zeros_int(part->n_face);

  for (int i = 0; i < nEdge; i++) {
    int face1 = PDM_ABS (part->edgeFace[2*i]);
    int face2 = PDM_ABS (part->edgeFace[2*i+1]);
    if (face1 > 0) {
      part->faceEdge[faceEdgeIdx[face1-1] + n_faceEdge[face1-1]++] = i + 1;
    }
    if (face2 > 0) {
      part->faceEdge[faceEdgeIdx[face2-1] + n_faceEdge[face2-1]++] = i + 1;
    }
  }

  /*
   * face -> edge orientation (same than face -> vtx)
   */

  int *vtxEdge = PDM_array_const_int(2 * part->n_vtx, -1);


  for (int i = 0; i < part->n_face; i++) {
    int idx = 0;
    for (int j = part->faceEdgeIdx[i]; j < part->faceEdgeIdx[i+1]; j++) {
      int _edge = part->faceEdge[j];
      int _iedge = _edge - 1;

      int ivtx1 = part->edgeVtx[2*_iedge] - 1;
      int ivtx2 = part->edgeVtx[2*_iedge+1] -1;

      if (vtxEdge[2*ivtx1] == -1) {
        vtxEdge[2*ivtx1] = _iedge;
      }
      else {
        vtxEdge[2*ivtx1+1] = _iedge;
      }

      if (vtxEdge[2*ivtx2] == -1) {
        vtxEdge[2*ivtx2] = _iedge;
      }
      else {
        vtxEdge[2*ivtx2+1] = _iedge;
      }
    }

    int n_vtx_face = part->face_vtx_idx[i+1] - part->face_vtx_idx[i];
    idx = part->face_vtx_idx[i];
    for (int j = part->face_vtx_idx[i]; j < part->face_vtx_idx[i+1]; j++) {
      int ivtx = part->face_vtx[j] - 1;
      int ivtx_next = part->face_vtx[idx + (j+1-idx) % n_vtx_face] -1;

      int iedge1;
      int iedge2;
      int k;
      int k1;

      for (k = 0; k < 2; k++) {
        iedge1 = vtxEdge[2*ivtx + k];
        assert (iedge1 != -1);
        for (k1 = 0; k1 < 2; k1++) {
          iedge2 = vtxEdge[2*ivtx_next + k1];
          assert (iedge2 != -1);
          if (iedge1 == iedge2) {
            break;
          }
        }
        if (iedge1 == iedge2) {
          break;
        }
      }
      assert (iedge1 == iedge2);

      if (part->edgeVtx[2*iedge1] == (ivtx + 1)) {
        part->faceEdge[j] = iedge1 + 1;
      }
      else {
        part->faceEdge[j] = -(iedge1 + 1);
      }
    }

    for (int j = part->faceEdgeIdx[i]; j < part->faceEdgeIdx[i+1]; j++) {
      int _edge = PDM_ABS(part->faceEdge[j]);
      int _iedge = _edge - 1;

      assert(_edge <= part->nEdge);

      int ivtx1 = part->edgeVtx[2*_iedge] - 1;
      int ivtx2 = part->edgeVtx[2*_iedge+1] -1;

      vtxEdge[2*ivtx1]   = -1;
      vtxEdge[2*ivtx1+1] = -1;

      vtxEdge[2*ivtx2]   = -1;
      vtxEdge[2*ivtx2+1] = -1;

    }
  }

  free (vtxEdge);

  free(edgeFaceUncompress);
  free(n_faceEdge);

}


/**
 * \brief Return face_ln_to_gn
 *
 *
 * \param [in]  part      Partition to compute
 *
 */

const PDM_g_num_t *
PDM_surf_part_faceLnToGn_get
(
PDM_surf_part_t *part
)
{

  return part->face_ln_to_gn;
}

/**
 * \brief Dump a PDM_surf_part_t object
 *
 * This function dumps a surf part structure
 *
 * \param [in]  part     surf part to dump
 *
 */

void
PDM_surf_part_dump
(
 PDM_surf_part_t *part
 )
{
  if (part != NULL) {
	  PDM_printf ("PDM_surf_part_t : \n");
    PDM_printf ("  - edgeLnToGn :");
    if (part->edgeLnToGn != NULL) {
      for (int j=0; j<part->nEdge; j++)
        PDM_printf (" %ld", j+1, part->edgeLnToGn[j]);
    }
    else {
      printf(" NULL");
    }
    printf("\n");

    PDM_printf ("  - n_face : %d \n", part->n_face);
    PDM_printf ("  - faceNormal :");
    if (part->faceNormal == NULL)
      PDM_printf (" NULL\n");
    else{
      for (int j=0; j<3*part->n_face; j++)
        PDM_printf ("%d ", part->faceNormal[j]);
      PDM_printf ("\n");
    }
    PDM_printf ("  - nGhostFace : %d \n", part->nGhostFace);
    PDM_printf ("  - nTotalFace : %d \n", part->nTotalFace);
    PDM_printf ("  - sface_vtx : %d \n", part->sface_vtx);
		PDM_printf ("  - face_vtx : \n");
		for (int j = 0; j < part->n_face; j++) {
		  for (int k = (part->face_vtx_idx)[j]; k < (part->face_vtx_idx)[j+1]; k++)
        PDM_printf (" %d", part->face_vtx[k]);
		  PDM_printf ("\n");
		}

		PDM_printf ("  - face_ln_to_gn : \n");
		for (int j = 0; j < part->n_face; j++) {
		  PDM_printf (" "PDM_FMT_G_NUM, part->face_ln_to_gn[j]);
    }
    PDM_printf ("\n");

    PDM_printf ("  - faceEdge : \n");
    if (part->faceEdge != NULL) {
      for (int j = 0; j <part->n_face; j++) {
        PDM_printf ("FG2EG %ld-> ", part->face_ln_to_gn[j]);
        for (int k = (part->faceEdgeIdx)[j]; k < (part->faceEdgeIdx)[j+1]; k++)
          PDM_printf (" "PDM_FMT_G_NUM,
                      (PDM_ABS(part->faceEdge[k])/ part->faceEdge[k])*  part->edgeLnToGn[PDM_ABS(part->faceEdge[k])-1]);
        PDM_printf ("\n");
      }
    }
    else {
      printf(" NULL\n");
    }

    PDM_printf ("  - n_vtx : %d \n", part->n_vtx);
		PDM_printf ("  - vtx_ln_to_gn : \n");
		for (int j = 0; j < part->n_vtx; j++) {
		  PDM_printf (" "PDM_FMT_G_NUM, part->vtx_ln_to_gn[j]);
    }
    PDM_printf ("\n");
    PDM_printf ("  - coords : \n");
    for (int j=0; j<part->n_vtx; j++){
      PDM_printf ("%d :", j);
      PDM_printf (" %12.5e", part->coords[3*j]);
      PDM_printf (" %12.5e", part->coords[3*j+1]);
      PDM_printf (" %12.5e", part->coords[3*j+2]);
      PDM_printf ("\n");
    }
    PDM_printf ("  - vtxEdge loc with  ghost edges: \n");
    if (part->vtxEdge != NULL) {
      for (int j = 0; j <part->n_vtx; j++) {
        PDM_printf ("VG2EG %ld-> ", part->vtx_ln_to_gn[j]);
        for (int k = (part->vtxEdgeIdx[j]); k < (part->vtxEdgeIdx)[j+1]; k++)
          PDM_printf (" %d",part->vtxEdge[k]);
        PDM_printf ("\n");
      }
    }
    else {
      printf(" NULL\n");
    }

    PDM_printf ("  - vtxEdge without ghost edges: \n");
    if (part->vtxEdge != NULL) {
      for (int j = 0; j <part->n_vtx; j++) {
        PDM_printf ("VG2VE %ld-> ", part->vtx_ln_to_gn[j]);
        for (int k = (part->vtxEdgeIdx[j]); k < (part->vtxEdgeIdx)[j+1]; k++)
          if (part->vtxEdge[k] <= part->nEdge) {
            PDM_printf (" "PDM_FMT_G_NUM, part->edgeLnToGn[part->vtxEdge[k]-1]);
          }
        PDM_printf ("\n");
      }
    }
    else {
      printf(" NULL\n");
    }

    PDM_printf ("  - nEdge : %d \n", part->nEdge);
    PDM_printf ("  - nGhostEdge : %d \n", part->nGhostEdge);
    PDM_printf ("  - nTotalEdge : %d \n", part->nTotalEdge);
    PDM_printf ("  - edgeFace loc : ");
    for (int j=0; j<part->nEdge; j++){
      PDM_printf ("\nEG2FG %ld : ", part->edgeLnToGn[j]);
      PDM_printf (" %d", part->edgeFace[2*j]);
      if (part->edgeFace[2*j+1] > 0)
        PDM_printf (" %d", part->edgeFace[2*j+1]);
    }
    PDM_printf ("\n");
    PDM_printf ("  - edgeFace with ghost face : ");
    if (part->edgeFace != NULL) {
      for (int j=0; j<part->nEdge; j++){
        PDM_printf ("\n%ld : ", part->edgeLnToGn[j]);
        if (part->edgeFace[2*j] <= part->n_face) {
          PDM_printf (" "PDM_FMT_G_NUM, part->face_ln_to_gn[part->edgeFace[2*j]-1]);
        }
        if (part->edgeFace[2*j+1] <= part->n_face) {
          if (part->edgeFace[2*j+1] > 0)
            PDM_printf (" "PDM_FMT_G_NUM, part->face_ln_to_gn[part->edgeFace[2*j+1]-1]);
        }
      }
    }
    else {
      printf(" NULL\n");
    }

    PDM_printf ("\n");
    PDM_printf ("  - edgeVtx loc without ghost faces : ");
    if (part->edgeVtx != NULL) {
      for (int j=0; j<part->nEdge; j++){
        PDM_printf ("\n%ld : ",  part->edgeLnToGn[j]);
        PDM_printf (" %d", part->edgeVtx[2*j]);
        PDM_printf (" %d", part->edgeVtx[2*j+1]);
      }
      PDM_printf ("\n");
    }
    else {
      printf(" NULL\n");
    }
    PDM_printf ("  - edgeVtx : ");
    if (part->edgeVtx != NULL) {
      for (int j=0; j<part->nEdge; j++){
        PDM_printf ("\nEG2VG %ld : ",  part->edgeLnToGn[j]);
        PDM_printf (" "PDM_FMT_G_NUM, part->vtx_ln_to_gn[part->edgeVtx[2*j]-1]);
        PDM_printf (" "PDM_FMT_G_NUM, part->vtx_ln_to_gn[part->edgeVtx[2*j+1]-1]);
      }
      PDM_printf ("\n");
    }
    else {
      printf(" NULL\n");
    }
    PDM_printf ("  - edgePartBound : %d \n", part->edgePartBound);
    PDM_printf ("  - vtxPartBound : %d \n", part->vtxPartBound);
    PDM_printf ("  - carLgthVtx : %d : ", part->carLgthVtx);
    if (part->carLgthVtx == NULL) {
      PDM_printf ("\n");
    }
    else {
      for (int j=0; j<part->n_vtx; j++) {
        PDM_printf ("%d ", part->carLgthVtx[j]);
      }
      PDM_printf ("\n");
    }
    PDM_printf ("  - extents : %d : ", part->extents);
    if (part->extents == NULL) {
      PDM_printf ("\n");
    }
    else {
      for (int j=0; j<part->n_vtx; j++) {
        PDM_printf ("%d ", part->extents[j]);
      }
      PDM_printf ("\n");
    }
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
