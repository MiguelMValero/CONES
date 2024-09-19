/*============================================================================
 * Hilbert encoding for 2D or 3D coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_part_graph.h"
#include "pdm_hilbert.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_ext_wrapper.h"
#include "pdm_array.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */

static void
_quickSort_int
(
 int a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    int t;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    _quickSort_int(a, l  , j-1);
    _quickSort_int(a, j+1,   r);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Splits the graph
 *
 * \param [in]  method       Method to be used: choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]  n_part        Number of partitions
 * \param [in]  part_ini     Part object fine mesh
 *
 * \param [in]  cell_cell                  Dual graph (size : cell_cell_idx[n_cell])
 * \param [in]  cell_cell_idx               Array of indexes of the dual graph (size : n_cell + 1)
 * \param [in]  cell_weight         Cell weight (size = n_cell)
 * \param [in]  face_weight         Face weight (size = n_face)
 *
 * \param [inout] cell_part  Cell partitioning (size : n_cell)
 *
 */

void
PDM_part_graph_split
(
 int         method,
 int         n_part,
 _part_t    *part_ini,
 int        *cell_cell_idx,
 int        *cell_cell,
 int        *cell_weight,
 int        *face_weight,
 int       **cell_part
)
{
  *cell_part = PDM_array_zeros_int(part_ini->n_cell);

  if (n_part > 1) {
    switch(method) {
    case 1:
      {
#ifdef PDM_HAVE_PARMETIS
        //            Define Metis properties

        //          int flag_weights = 0; //0 = False -> weights are unused

        int flag_weights = 0; //0 = False -> weights are unused

        int ncon = 1; //The number of balancing constraints

        int *vwgt = cell_weight; //Weights of the vertices of the graph (NULL if unused)

        int *adjwgt = face_weight; //Weights of the edges of the graph (NULL if unused)

        double *tpwgts = NULL;
        if (flag_weights != 0) {
          tpwgts = (double *) malloc(ncon * n_part * sizeof(double));
          for (int i = 0; i < ncon * n_part; i++){
            tpwgts[i] = (double) (1./n_part);
          }
        }

        double *ubvec = NULL;
        if (flag_weights != 0) {
          ubvec = (double *) malloc(ncon * sizeof(double));
          for (int i = 0; i < ncon; i++) {
            ubvec[i] = 1.05;
          }
        }

        //TO ADD: USE OF ADJWGT IN AN IF STATEMENT

        //This value is solely a memory space to be filled by METIS

        int edgecut;
        printf("PDM_part_graph_split \n");
        if (n_part < 8) {

          PDM_METIS_PartGraphRecursive (&(part_ini->n_cell),
                                        &ncon,
                                        cell_cell_idx,
                                        cell_cell,
                                        vwgt,
                                        adjwgt,
                                        &n_part,
                                        tpwgts,
                                        ubvec,
                                        &edgecut,
                                        *cell_part);
        }

        else {

          PDM_METIS_PartGraphKway (&(part_ini->n_cell),
                                   &ncon,
                                   cell_cell_idx,
                                   cell_cell,
                                   vwgt,
                                   adjwgt,
                                   &n_part,
                                   tpwgts,
                                   ubvec,
                                   &edgecut,
                                   *cell_part);
        }
        // double inbalance = 0.03;
        // double balance   = 0;
        // edgecut   = 0;
        // // bool   suppress_output = False;
        // // bool   graph_partitioned = False;
        // int time_limit = 0;
        // int seed  = 0;
        // int mode = 2;
        // PDM_kaffpa(&(part_ini->n_cell),
        //            NULL,
        //            cell_cell_idx,
        //            NULL,
        //            cell_cell,
        //            &n_part,
        //             &inbalance,
        //             seed,
        //             mode,
        //             &edgecut,
        //             *cell_part );

        if (0 == 1) {
          PDM_printf("\n Contenu de cell_part : \n");
          for (int i = 0; i < part_ini->n_cell; i++) {
            PDM_printf(" %d ", (*cell_part)[i]);
          }
          PDM_printf("\n");
        }

        if (flag_weights != 0) {
          if(ubvec!= NULL)
            free(ubvec);
          if(tpwgts!= NULL)
            free(tpwgts);
          // if(adjwgt!= NULL)
          //   free(adjwgt);
        }

#else
        PDM_UNUSED(method);
        PDM_UNUSED(n_part);
        PDM_UNUSED(part_ini);
        PDM_UNUSED(cell_cell_idx);
        PDM_UNUSED(cell_cell);
        PDM_UNUSED(cell_weight);
        PDM_UNUSED(face_weight);
        PDM_UNUSED(cell_part);
        PDM_printf("PDM_part error : METIS unavailable\n");
        exit(1);

#endif
        break;
      }
    case 2:
      {
#ifdef PDM_HAVE_PTSCOTCH

        int check = 0;

        PDM_SCOTCH_part (part_ini->n_cell,
                         cell_cell_idx,
                         cell_cell,
                         cell_weight,
                         face_weight,
                         check,
                         n_part,
                         *cell_part);

#else
        PDM_UNUSED(method);
        PDM_UNUSED(n_part);
        PDM_UNUSED(part_ini);
        PDM_UNUSED(cell_cell_idx);
        PDM_UNUSED(cell_cell);
        PDM_UNUSED(cell_weight);
        PDM_UNUSED(face_weight);
        PDM_UNUSED(cell_part);
        PDM_printf("PDM_part error : Scotch unavailable\n");
        exit(1);
#endif

        break;
      }
      case 3:
      {

        // To see with eric ...
        // abort();
        /* Allocation */
        double *cellCenter = (double *) malloc (part_ini->n_cell * 3 * sizeof(double ));

        PDM_hilbert_code_t *hilbert_codes = (PDM_hilbert_code_t *) malloc (part_ini->n_cell * sizeof(PDM_hilbert_code_t));

        /** Barycentre computation **/

        /* Allocate */
        double *cellPond = (double *) malloc (part_ini->n_cell * sizeof(double));

        /* Nulliffy cellCenterArray */
        for(int iCell = 0; iCell < part_ini->n_cell; iCell++) {
          cellCenter[3*iCell  ] = 0.;
          cellCenter[3*iCell+1] = 0.;
          cellCenter[3*iCell+2] = 0.;
          cellPond[iCell]     = 0.;
        }

        /* Compute */
        for(int iCell = 0; iCell < part_ini->n_cell; iCell++) {

          /* Cellule composé de n_face */
          int aFac = part_ini->cell_face_idx[iCell];
          int nFac = part_ini->cell_face_idx[iCell+1] - aFac;

          for(int iFac = 0; iFac < nFac; iFac++) {

            /* Face composé de n_vtx */
            int lFac = PDM_ABS(part_ini->cell_face[aFac + iFac]) - 1;

            int aVtx = part_ini->face_vtx_idx[lFac];
            int n_vtx = part_ini->face_vtx_idx[lFac+1] - aVtx;

            for(int iVtx = 0; iVtx < n_vtx; iVtx++) {

              /* Face composé de n_vtx */
              int lVtx = part_ini->face_vtx[aVtx + iVtx] - 1;

              /* Add to current cell and stack weight */
              cellCenter[3*iCell  ] += part_ini->vtx[3*lVtx  ];
              cellCenter[3*iCell+1] += part_ini->vtx[3*lVtx+1];
              cellCenter[3*iCell+2] += part_ini->vtx[3*lVtx+2];

              cellPond[iCell] += 1.;
            }
          }
        }

        /* Nulliffy cellCenterArray */
        for(int iCell = 0; iCell < part_ini->n_cell; iCell++) {
          cellCenter[3*iCell  ] = cellCenter[3*iCell  ]/cellPond[iCell];
          cellCenter[3*iCell+1] = cellCenter[3*iCell+1]/cellPond[iCell];
          cellCenter[3*iCell+2] = cellCenter[3*iCell+2]/cellPond[iCell];
        }


        double extents[3 * 2];

        /** Get EXTENTS LOCAL **/

        PDM_hilbert_get_coord_extents_seq(3, part_ini->n_cell, cellCenter, extents);

        /** Hilbert Coordinates Computation **/

        PDM_hilbert_encode_coords(3, PDM_HILBERT_CS, extents, part_ini->n_cell, cellCenter, hilbert_codes);

        /** CHECK H_CODES **/

        free(cellCenter);
        free(cellPond);

        int *newToOldOrder = (int *) malloc (part_ini->n_cell * sizeof(int));
        for(int i = 0; i < part_ini->n_cell; ++i) {
          newToOldOrder [i] = i;
        }

        PDM_sort_double (hilbert_codes, *cell_part, part_ini->n_cell);

        /* Free */
        free (hilbert_codes);
        free (newToOldOrder);

        break;
      }
    default:
      PDM_printf("PART error : '%i' unknown partitioning method\n", method);
      exit(1);
    }
  }
  else if (n_part < 0) {
    PDM_printf("PART error : n_part must be > 0\n");
    exit(1);
  }
}


/**
 *
 * \brief Builds dual graph face cell connectivity
 *
 * \param [inout] part_ini                 Part object - fine mesh partition
 *
 * \param [inout] cell_cell_idxCompressed    Array of indexes of the dual graph
 * \param [inout] cell_cellCompressed       Dual graph
 *
 */
void
PDM_part_graph_compute_from_face_cell
(
  _part_t        *part_ini,
  int           **cell_cell_idxCompressed,
  int           **cell_cellCompressed
)
{
  int *cell_cell_n = PDM_array_zeros_int(part_ini->n_cell);

  int *cell_cell = PDM_array_const_int(part_ini->cell_face_idx[part_ini->n_cell], -1);

  int *cell_cell_idx = (int *) malloc((part_ini->n_cell + 1) * sizeof(int));
  for(int i = 0; i < part_ini->n_cell + 1; i++) {
    cell_cell_idx[i] = part_ini->cell_face_idx[i];
  }

  for (int i = 0; i < part_ini->n_face; i++) {
    int i_cell1 = PDM_ABS (part_ini->face_cell[2*i    ]);
    int i_cell2 = PDM_ABS (part_ini->face_cell[2*i + 1]);
    //Only the non-boundary faces are stored
    if (i_cell1 > 0 && i_cell2 > 0) {
      int idx1 = cell_cell_idx[i_cell1-1] + cell_cell_n[i_cell1-1];
      cell_cell[idx1] = i_cell2 ;
      cell_cell_n[i_cell1-1] += 1;

      int idx2 = cell_cell_idx[i_cell2-1] + cell_cell_n[i_cell2-1];
      cell_cell[idx2] = i_cell1 ;
      cell_cell_n[i_cell2-1] += 1;
    }

  }

  if (0 == 1) {
    PDM_printf("Content of cell_cell_n after looping over cell_face: ");
    for(int i = 0; i < part_ini->n_cell; i++) {
      PDM_printf(" %d ", cell_cell_n[i]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cell_cell after looping over cell_face: ");
    for(int i = 0; i < part_ini->cell_face_idx[part_ini->n_cell]; i++) {
      PDM_printf(" %d ", cell_cell[i]);
    }
    PDM_printf("\n");
  }

  //cell_cell_idx is rebuilt
  assert((*cell_cell_idxCompressed) == NULL);
  *cell_cell_idxCompressed = PDM_array_new_idx_from_sizes_int(cell_cell_n, part_ini->n_cell);

  //We compress the dual graph since cell_cell_idx was built from cell_face_idx
  //We have then n_face elements in cell_cell whereas it needs to be composed of n_cell elements

  //    PDM_printf("(*cell_cell_idxCompressed)[part_ini->n_cell] : %d \n", (*cell_cell_idxCompressed)[part_ini->n_cell]);
  //
  assert( (*cell_cellCompressed) == NULL);
  (*cell_cellCompressed) = (int *) malloc((*cell_cell_idxCompressed)[part_ini->n_cell] * sizeof(int));

  int cpt_cell_cellCompressed = 0;
  for(int i = 0; i < part_ini->cell_face_idx[part_ini->n_cell]; i++) {
    //        PDM_printf("I am testing a value for the %d time! \n", i);

    //We have an information to store when a neighboring cell exists
    if(cell_cell[i] > -1){
      //            PDM_printf("I am storing a value for the %d time! \n", i);
      //We add a -1 to have the graph vertices numbered from 0 to n (C numbering)
      (*cell_cellCompressed)[cpt_cell_cellCompressed++] = cell_cell[i] - 1;
      //            PDM_printf("Valeur stockee : %d \n ", (*cell_cellCompressed)[cpt_cell_cellCompressed - 1]);
    }
  }

  if( 0 == 1) {
    PDM_printf("Content of cell_cellCompressed after compression and renumbering: ");
    for(int i = 0; i < (*cell_cell_idxCompressed)[part_ini->n_cell]; i++) {
      PDM_printf(" %d ", (*cell_cellCompressed)[i]);
    }
    PDM_printf("\n");
  }

  /* Free temporary arrays*/

  free(cell_cell_n);
  free(cell_cell);
  free(cell_cell_idx);

  //Remove duplicate cells of the dual graph
  //We use the following scheme:
  //We loop over the indexes for the whole array to subdivide it into subarrays
  //We sort locally each subarray (determined thanks to cell_cell_idxCompressed)
  //We loop over each subarray
  //We store the first value of each subarray anyway
  //We store each non-duplicated value and increment the writing index
  //We update the index array at each iteration

  int idx_write = 0;
  int tabIdxTemp = 0;

  for (int i = 0; i < part_ini->n_cell; i++) {
    _quickSort_int((*cell_cellCompressed), tabIdxTemp, (*cell_cell_idxCompressed)[i + 1] - 1);

    int last_value = -1;

    for (int j = tabIdxTemp; j < (*cell_cell_idxCompressed)[i + 1]; j++) {
      //We need to have a local index (between 0 and n_face)
      //If the value is different from the previous one (higher than is the same as different since the array is sorted)

      if(last_value != (*cell_cellCompressed)[j]) {
        (*cell_cellCompressed)[idx_write++] = (*cell_cellCompressed)[j];
        last_value = (*cell_cellCompressed)[j];
      }
    }

    if (0 == 1) {
      PDM_printf("\n Contenu de cell_cellCompressed apres reecriture: \n");
      for(int i1 = 0; i1 < (*cell_cell_idxCompressed)[part_ini->n_cell]; i1++) {
        PDM_printf(" %d ", (*cell_cellCompressed)[i1]);
      }
      PDM_printf("\n");
    }

    tabIdxTemp = (*cell_cell_idxCompressed)[i + 1];
    (*cell_cell_idxCompressed)[i + 1] = idx_write;

    if (0 == 1) {
      PDM_printf("\n Contenu de cell_cell_idxCompressed apres reecriture: \n");
      for(int i1 = 0; i1 < part_ini->n_cell + 1; i1++) {
        PDM_printf(" %d ", (*cell_cell_idxCompressed)[i1]);
      }
      PDM_printf("\n");
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cell_cell_idxCompressed after compression: ");
    for(int i1 = 0; i1 < part_ini->n_cell + 1; i1++) {
      PDM_printf(" %d ", (*cell_cell_idxCompressed)[i1]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cell_cellCompressed after compression: ");
    for(int i1 = 0; i1 < (*cell_cell_idxCompressed)[part_ini->n_cell]; i1++) {
      PDM_printf(" %d ", (*cell_cellCompressed)[i1]);
    }
    PDM_printf("\n");
  }

  //We reallocate the memory in case of duplicated values removed
  //The new array size is idx_write (stored in (*cell_cell_idxCompressed)[part_ini->n_cell])
  *cell_cellCompressed = realloc(*cell_cellCompressed,
                                (*cell_cell_idxCompressed)[part_ini->n_cell] * sizeof(int));

}




void
PDM_part_graph_split_bis
(
 int         method,
 int         n_part,
 int         graph_size,
 int        *cell_cell_idx,
 int        *cell_cell,
 int        *cell_weight,
 int        *face_weight,
 int       **cell_part
)
{
  *cell_part = PDM_array_zeros_int(graph_size);

  PDM_UNUSED (n_part);
  PDM_UNUSED (cell_cell_idx);
  PDM_UNUSED (cell_cell);
  PDM_UNUSED (cell_weight);
  PDM_UNUSED (face_weight);

  switch(method) {
  case 1:
    {
#ifdef PDM_HAVE_PARMETIS
      //            Define Metis properties

      //          int flag_weights = 0; //0 = False -> weights are unused

      int flag_weights = 1; //0 = False -> weights are unused

      int ncon = 1; //The number of balancing constraints

      int *vwgt = cell_weight; //Weights of the vertices of the graph (NULL if unused)

      int *adjwgt = face_weight; //Weights of the edges of the graph (NULL if unused)

      double *tpwgts = NULL;
      if (flag_weights != 0) {
        tpwgts = (double *) malloc(ncon * n_part * sizeof(double));
        for (int i = 0; i < ncon * n_part; i++){
          tpwgts[i] = (double) (1./n_part);
        }
      }

      double *ubvec = NULL;
      if (flag_weights != 0) {
        ubvec = (double *) malloc(ncon * sizeof(double));
        for (int i = 0; i < ncon; i++) {
          ubvec[i] = 1.05;
        }
      }

      //TO ADD: USE OF ADJWGT IN AN IF STATEMENT

      //This value is solely a memory space to be filled by METIS

      int edgecut;

      if (n_part < 8) {

        PDM_METIS_PartGraphRecursive (&(graph_size),
                                      &ncon,
                                      cell_cell_idx,
                                      cell_cell,
                                      vwgt,
                                      adjwgt,
                                      &n_part,
                                      tpwgts,
                                      ubvec,
                                      &edgecut,
                                      *cell_part);
      }

      else {

        PDM_METIS_PartGraphKway (&(graph_size),
                                 &ncon,
                                 cell_cell_idx,
                                 cell_cell,
                                 vwgt,
                                 adjwgt,
                                 &n_part,
                                 tpwgts,
                                 ubvec,
                                 &edgecut,
                                 *cell_part);
      }

      // double inbalance = 0.03;
      // double balance   = 0;
      // edgecut   = 0;
      // // bool   suppress_output = False;
      // // bool   graph_partitioned = False;
      // int time_limit = 0;
      // int seed  = 0;
      // int mode = 2;
      // PDM_kaffpa(&(graph_size),
      //            NULL,
      //            cell_cell_idx,
      //            NULL,
      //            cell_cell,
      //            &n_part,
      //             &inbalance,
      //             seed,
      //             mode,
      //             &edgecut,
      //             *cell_part );

      if (0 == 1) {
        PDM_printf("\n Contenu de cell_part : \n");
        for (int i = 0; i < graph_size; i++) {
          PDM_printf(" %d ", (*cell_part)[i]);
        }
        PDM_printf("\n");
      }

      if (flag_weights != 0) {
        if(ubvec!= NULL)
          free(ubvec);
        if(tpwgts!= NULL)
          free(tpwgts);
        // if(adjwgt!= NULL)
        //   free(adjwgt);
      }

#else
      PDM_printf("PDM_part error : METIS unavailable\n");
      exit(1);

#endif
      break;
    }
  case 2:
    {
#ifdef PDM_HAVE_PTSCOTCH

      int check = 0;

      PDM_SCOTCH_part (graph_size,
                       cell_cell_idx,
                       cell_cell,
                       cell_weight,
                       face_weight,
                       check,
                       n_part,
                       *cell_part);

#else
      PDM_printf("PDM_part error : Scotch unavailable\n");
      exit(1);
#endif

      break;
    }
    case 3:
    {

      // To see with eric ...
      abort();
      break;
    }
  default:
    PDM_printf("PART error : '%i' unknown partitioning method\n", method);
    exit(1);
  }

}


#ifdef  __cplusplus
}
#endif
