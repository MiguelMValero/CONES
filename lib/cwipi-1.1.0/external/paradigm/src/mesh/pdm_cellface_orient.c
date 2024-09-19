
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
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_array.h"
#include "pdm_geom_elem.h"
#include "pdm_cellface_orient.h"
#include "pdm_hash_tab.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_logging.h"

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


enum {false, true};

typedef enum {
  FACE_Un_procESSED,
  FACE_IN_STACK,
  FACE_UNCHANGED_CYCLE,
  FACE_CHANGED_CYCLE
} _face_state_t;


typedef enum {
  CELL_Un_procESSED,
  CELL_IN_STACK,
  CELL_COMPLETED
} _cell_state_t;

/*============================================================================
 * Static global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/



/**
 * \brief Orient cell->face connectivity
 *
 * At the output of th function, a face number in \ref cell_face is positive
 * if surface normal is inside the cell, negative otherwise. \ref face_cell is
 * oriented in the same way
 *
 * \param [in]      n_cell         Number of cells
 * \param [in]      n_face         Number of faces
 * \param [in]      n_vtx          Number of vertices
 * \param [in]      coords         Vertices coordinates
 * \param [in]      cell_face_idx  Cell to face connectivity index (size = \ref n_cell + 1)
 * \param [in, out] cell_face      Cell to face connectivity (size = cell_face_idx[n_cell])
 * \param [in, out] face_cell      Face to cell connectivity (size = 2 * \ref n_face) or NULL
 * \param [in]      face_vtx_idx   Face to vertex connectivity index (size = \ref n_face + 1)
 * \param [in]      face_vtx       Face to vertex connectivity (size = face_vtx_idx[n_face])
 *
 */

void
PDM_cellface_orient
(
const int      n_cell,
const int      n_face,
const int      n_vtx,
const double  *coords,
const int     *cell_face_idx,
int           *cell_face,
int           *face_cell,
const int     *face_vtx_idx,
const int     *face_vtx
)
{




  if (n_cell == 0) {
    return;
  }

  PDM_timer_t *t1 = PDM_timer_create();
  PDM_timer_resume(t1);

  int *_faceCell = NULL;

  if (face_cell != NULL) {
    _faceCell = face_cell;
  }
  else {
    _faceCell = PDM_array_zeros_int(2*n_face);

    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = PDM_ABS (cell_face[j]);
        int iFace = 2 * (face - 1);
        if (_faceCell[iFace] == 0) {
          _faceCell[iFace] = i + 1;
        }
        else {
          _faceCell[iFace+1] = i + 1;
        }
      }
    }
  }

  if (1 == 0) {
    printf("_faceCell : ");
    for (int i = 0; i < n_face; i++) {
      printf("%d : %d %d\n", i+1, _faceCell[2*i], _faceCell[2*i+1]);
    }
    printf("\n");
  }

  int *orientedface_cell = PDM_array_zeros_int(2*n_face);


  int keyMax = 2 * n_vtx;
  PDM_hash_tab_t *hashOrient = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT, &keyMax);

  int maxNPolyFace = -1;
  int maxEdges = -1;

  for (int ipoly = 0; ipoly < n_cell; ipoly++) {
    const int polyIdx   = cell_face_idx[ipoly];
    const int nPolyFace = cell_face_idx[ipoly + 1] - polyIdx;
    maxNPolyFace = PDM_MAX (maxNPolyFace, nPolyFace);
    int nEdgeCell = 0;
    for (int i = polyIdx; i < polyIdx + nPolyFace; i++) {
      int face = PDM_ABS(cell_face[i]) - 1;
      nEdgeCell += face_vtx_idx[face+1] - face_vtx_idx[face];

    }
    maxEdges = PDM_MAX (maxEdges, nEdgeCell);
  }

  int *stackFace = (int *) malloc (sizeof(int) * maxNPolyFace);

  int nStackCell = -1;
  int *stackCell = (int *) malloc (sizeof(int) * n_cell);
  int n_processedFace = 0;
  int *processedFace = (int *) malloc (sizeof(int) * maxNPolyFace);
  int *tagCell = PDM_array_const_int(n_cell, CELL_Un_procESSED);
  int *tagFace = PDM_array_const_int(maxNPolyFace, FACE_Un_procESSED);


  tagCell[0] = CELL_COMPLETED;

  int nEdges = 0;
  const int nDataEdge = 3;
  int *edges = malloc (sizeof(int) * maxEdges * nDataEdge);

 /*
  * Orient the first cell of the first component
  * --------------------------------------------
  *
  * As the oriented volume is positive, the face normals
  * are outside of the element
  *
  */

  int fistCellComp = 0;

  while (fistCellComp != -1) {

    int     isOriented = 0;
    int     nPolyhedra = 1;
    double  volume[3];
    double  center[3];

    int *_cell_face_idx = (int *) cell_face_idx + fistCellComp;

    PDM_geom_elem_polyhedra_properties (isOriented,
                                        nPolyhedra,
                                        n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        _cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        coords,
                                        volume,
                                        center,
                                        NULL,
                                        NULL);

   /*
    * Initialize an oriented face cell
    * --------------------------------
    *
    * left cell : normal is inside the cell
    * right cell : normal is outside the cell
    *
    * The orientation of the first cell is taking into account
    *
    */

    tagCell[fistCellComp] = CELL_COMPLETED;

    for (int i = cell_face_idx[fistCellComp]; i < cell_face_idx[fistCellComp+1]; i++) {
      int Face = cell_face[i];

      if (Face > 0) {
        orientedface_cell[2 * (Face - 1)] = fistCellComp + 1;
      }
      else {
        orientedface_cell[2 * (PDM_ABS (Face) - 1) + 1] = fistCellComp + 1;
      }
    }


   /*
    * Other cells are oriented from the faces of the first cell
    * ---------------------------------------------------------
    *
    */

    /* Add neighbours of the first cell in the stack */

    for (int i = cell_face_idx[fistCellComp]; i < cell_face_idx[fistCellComp+1]; i++) {
      int iFace = 2 * (PDM_ABS (cell_face[i]) - 1);

      // if (_faceCell[iFace] == 1) {
      if (_faceCell[iFace] == fistCellComp+1) {
        if (_faceCell[iFace + 1] > 0) {
          int cell = PDM_ABS (_faceCell[iFace + 1]);
          if (tagCell[cell - 1] == CELL_Un_procESSED) {
            stackCell[++nStackCell] = cell;
            tagCell[cell - 1] = CELL_IN_STACK;
          }
        }
      }
      else {
        int cell = PDM_ABS (_faceCell[iFace]);
        if (tagCell[cell - 1] == CELL_Un_procESSED) {
          stackCell[++nStackCell] = cell;
        }
        tagCell[cell - 1] = CELL_IN_STACK;
      }
    }

    /* Orientation process */

    while (nStackCell >= 0) {
      if (1 == 0) {
        printf("orientedface_cell : ");
        for (int i = 0; i < n_face; i++) {
          printf("%d : %d %d\n", i+1,orientedface_cell[2*i], orientedface_cell[2*i+1]);
        }
        printf("\n");

        printf("\nstackCell : ");
        for (int i = 0; i < nStackCell; i++) {
          printf(" %d", stackCell[i]);
        }
        printf("\n");
      }

      nEdges = 0;

      int iCell = stackCell[nStackCell--] - 1;

      if (1 == 0) {
        printf("iCell : %d\n", iCell);
      }

      if (tagCell[iCell] == CELL_COMPLETED) {
        continue;
      }

      const int polyIdx   = cell_face_idx[iCell];
      const int nPolyFace = cell_face_idx[iCell + 1] - polyIdx;

      /* Build pseudo edges of the current cell and store them into a hash table */

      n_processedFace = 0;

      for (int iface = 0; iface < nPolyFace; iface++) {

        tagFace[iface] = FACE_Un_procESSED;

        const int face          = PDM_ABS (cell_face[polyIdx + iface]) - 1;

        if (orientedface_cell[2*face] != 0) {
          assert (orientedface_cell[2*face  ] != iCell + 1);
          assert (orientedface_cell[2*face+1] != iCell + 1);
          assert (orientedface_cell[2*face+1] == 0        );
          tagFace[iface] = FACE_CHANGED_CYCLE;
          orientedface_cell[2*face+1] = iCell + 1;
          processedFace[n_processedFace++] = iface;
        }
        else if (orientedface_cell[2*face + 1] != 0) {
          assert (orientedface_cell[2*face  ] != iCell + 1);
          assert (orientedface_cell[2*face+1] != iCell + 1);
          assert (orientedface_cell[2*face  ] == 0        );
          tagFace[iface] = FACE_UNCHANGED_CYCLE;
          orientedface_cell[2*face] = iCell + 1;
          processedFace[n_processedFace++] = iface;
        }
      }

      if (n_processedFace == 0) {
        PDM_error (__FILE__, __LINE__, 0, "Error reorient : no processed face found\n");
      }

      for (int iface = 0; iface < nPolyFace; iface++) {

        const int face          = PDM_ABS (cell_face[polyIdx + iface]) - 1;
        const int faceIdx       = face_vtx_idx[face];
        const int n_faceVertices = face_vtx_idx[face+1] - faceIdx;

        for (int ivert = 0; ivert < n_faceVertices; ivert++) {
          const int vertex = face_vtx[faceIdx + ivert] - 1;

          const int inext = (ivert + 1) % n_faceVertices;
          const int vertexNext = face_vtx[faceIdx + inext] - 1;
          const int key = vertex + vertexNext;

          int *edge = edges + nDataEdge * nEdges;
          edge[0] = vertex;
          edge[1] = vertexNext;
          edge[2] = iface;

          nEdges += 1;

          PDM_hash_tab_data_add (hashOrient, (void *) &key, edge);

        }
      }

      /* FIXME Check edges connexion (ignore edges connected with more 2 faces) */

      int nStackFace = -1;

      /* Look for a neighbour of this face */

      for (int i = 0; i < n_processedFace; i++) {
        int _currentProcessedFace = processedFace[i];
        const int face          = PDM_ABS (cell_face[polyIdx + _currentProcessedFace]) - 1;
        const int faceIdx       = face_vtx_idx[face];
        const int n_faceVertices = face_vtx_idx[face+1] - faceIdx;

        for (int ivert = 0; ivert < n_faceVertices; ivert++) {
          const int inext = (ivert + 1) % n_faceVertices;

          const int vertex = face_vtx[faceIdx + ivert] - 1;
          const int vertexNext = face_vtx[faceIdx + inext] - 1;
          int key = vertex + vertexNext;

          int nData = PDM_hash_tab_n_data_get (hashOrient, &key);

          int **data = (int **) PDM_hash_tab_data_get (hashOrient, &key);

          for (int j = 0; j < nData; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int isInverseEdge = (vertex == _edge[1]) && (vertexNext == _edge[0]);
              int isSameEdge    = (vertex == _edge[0]) && (vertexNext == _edge[1]);
              int isSameFace    = _currentProcessedFace == _edge[2];
              if (!isSameFace) {
                if (isSameEdge || isInverseEdge) {
                  if (tagFace[ _edge[2]] == FACE_Un_procESSED) {
                    stackFace[++nStackFace] = _edge[2];
                    tagFace[ _edge[2]] = FACE_IN_STACK;
                  }
                  break;
                }
              }
            }
          }
        }
      }
      if (1 == 0) {
        printf("\nstackFace : ");
        for (int i = 0; i < nStackFace+1; i++) {
          printf(" %d", stackFace[i]);
        }
        printf("\n");
      }

      while (nStackFace >= 0) {

        int iFace = stackFace[nStackFace--];

        if ((tagFace[iFace] == FACE_UNCHANGED_CYCLE) ||
            (tagFace[iFace] == FACE_CHANGED_CYCLE)) {
          continue;
        }

        const int face          = PDM_ABS (cell_face[polyIdx + iFace]) - 1;
        const int faceIdx       = face_vtx_idx[face];
        const int n_faceVertices = face_vtx_idx[face+1] - faceIdx;

        for (int ivert = 0; ivert < n_faceVertices; ivert++) {
          const int inext = (ivert + 1) % n_faceVertices;

          const int vertex = face_vtx[faceIdx + ivert] - 1;
          const int vertexNext = face_vtx[faceIdx + inext] - 1;
          int key = vertex + vertexNext;

          int nData = PDM_hash_tab_n_data_get (hashOrient, &key);
          int **data = (int **) PDM_hash_tab_data_get (hashOrient, &key);

          int jCurrentEdge = -1;
          for (int j = 0; j < nData; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int isSameEdge    = (vertex == _edge[0]) && (vertexNext == _edge[1]);
              int isSameFace    = iFace == _edge[2];
              if (isSameEdge && isSameFace) {
                jCurrentEdge = j;
                break;
              }
            }
          }

          assert (jCurrentEdge > -1);

          for (int j = 0; j < nData; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int isInverseEdge = (vertex == _edge[1]) && (vertexNext == _edge[0]);
              int isSameEdge    = (vertex == _edge[0]) && (vertexNext == _edge[1]);
              int isSameFace    = iFace == _edge[2];

              int neighbour = _edge[2];

              if (!isSameFace) {
                if (isSameEdge || isInverseEdge) {

                  if (tagFace[iFace] < FACE_UNCHANGED_CYCLE) {

                    if (tagFace[neighbour] >= FACE_UNCHANGED_CYCLE) {
                      if (tagFace[neighbour] == FACE_UNCHANGED_CYCLE) {
                        if (isSameEdge) {
                          tagFace[iFace] = FACE_CHANGED_CYCLE;
                          orientedface_cell[2*face+1] = iCell+1;
                        }
                        else  {
                          tagFace[iFace] = FACE_UNCHANGED_CYCLE;
                          orientedface_cell[2*face] = iCell+1;
                        }
                      }
                      else {
                        if (isSameEdge) {
                          tagFace[iFace] = FACE_UNCHANGED_CYCLE;
                          orientedface_cell[2*face] = iCell+1;
                        }
                        else  {
                          tagFace[iFace] = FACE_CHANGED_CYCLE;
                          orientedface_cell[2*face+1] = iCell+1;
                        }
                      }
                    }
                  }

                  if (tagFace[neighbour] == FACE_Un_procESSED) {
                    stackFace[++nStackFace] = neighbour;
                    tagFace[neighbour] = FACE_IN_STACK;
                  }

                  break;

                }
              }
            }
          }
        }

        if (tagFace[iFace] == FACE_IN_STACK) {
          // printf(" oooo %i \n", iFace);
          // printf(" oooo %i \n", face);
          PDM_error (__FILE__, __LINE__, 0, "Error reorient : no neighbour processed face found\n");
        }
      }

      /* Add cell neighbours in the stack */

      for (int iface = 0; iface < nPolyFace; iface++) {

        if (!((tagFace[iface] == FACE_UNCHANGED_CYCLE) ||
              (tagFace[iface] == FACE_CHANGED_CYCLE))) {
          printf(" oooo %i \n", iface);
          printf(" Link to : %i %i \n", face_cell[2*iface], face_cell[2*iface+1]);
          PDM_error (__FILE__, __LINE__, 0, "Error reorient : a face of polyhedron is not processed\n");
        }

        if (tagFace[iface] == FACE_CHANGED_CYCLE) {
          cell_face[polyIdx + iface] = -cell_face[polyIdx + iface];
        }

        tagFace[iface] = FACE_Un_procESSED;

        const int face          = PDM_ABS (cell_face[polyIdx + iface]) - 1;

        int nextCell = -1;
        if (_faceCell[2 * face] == (iCell + 1)) {
          if (_faceCell[2 * face + 1] != 0) {
            nextCell = _faceCell[2 * face + 1] - 1;
          }
        }
        else {
          nextCell = _faceCell[2 * face] - 1;
        }

        if (nextCell != -1) {
          if (tagCell[nextCell] == CELL_Un_procESSED) {
            stackCell[++nStackCell] = nextCell + 1;
            tagCell[nextCell] = CELL_IN_STACK;
          }
        }
      }

      PDM_hash_tab_purge(hashOrient, PDM_FALSE);

      tagCell[iCell] = CELL_COMPLETED;

    }

    int icheck = fistCellComp;
    fistCellComp = -1;
    for (int k = icheck; k < n_cell; k++) {
      if (tagCell[k] == CELL_Un_procESSED) {
        fistCellComp = k;
        break;
      }
    }
  }

  PDM_hash_tab_free(hashOrient);

  /* Orient face_cell */

  if (face_cell != NULL) {
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = cell_face[j];
        int iFace = 2 * (PDM_ABS(face) - 1);
        if (PDM_ABS (face_cell[iFace]) == i+1) {
          if (face < 0) {
            face_cell[iFace] = -(i+1);
          }
          else {
            face_cell[iFace] = i+1;
          }
        }
        else {
          if (face < 0) {
            face_cell[iFace+1] = -(i+1);
          }
          else {
            face_cell[iFace+1] = i+1;
          }
        }
      }
    }
  }

  free (edges);
  free (stackFace);
  free (tagFace);
  free (stackCell);
  free (tagCell);
  free (orientedface_cell);
  free (processedFace);
  if (face_cell == NULL) {
    free (_faceCell);
  }

  PDM_timer_hang_on (t1);
  double et1 = PDM_timer_elapsed (t1);
  PDM_timer_free (t1);

  if (0 == 1) {
    printf("elapsed time cell_face_orient : %12.5e\n", et1);
  }
}

#ifdef	__cplusplus
}
#endif
