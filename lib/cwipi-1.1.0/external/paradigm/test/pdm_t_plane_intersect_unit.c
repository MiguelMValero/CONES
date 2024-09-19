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
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part.h"

#include "pdm_writer.h"
#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_overlay.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_edges_intersect.h"

/*============================================================================
 * Cas test
 *============================================================================*/

/*  /\* CAS 1 : un triangle et un quadrangle */
/* 0,2;5,4;5,0 */
/* 2,5;6,3;6,1;2,0 *\/ */
/* static int t1( */
/*   const int myRank, */
/*   double **faceVtxCooA */
/* /\*  PDM_edges_get_t *get_t, */
/*   PDM_g_num_t *nGvtxA, */
/*   PDM_g_num_t *nGvtxB, */
/*   PDM_g_num_t *maxGNEdgeA, */
/*   PDM_g_num_t *maxGNEdgeB, */
/*   int *nVtxA, */
/*   PDM_g_num_t **faceToEdgeA, */
/*   PDM_g_num_t **faceToVtxA, */
/*   double **faceVtxCooA, */
/*   double **faceVtxEpsA, */
/*   int *nVtxB, */
/*   PDM_g_num_t **faceToEdgeB, */
/*   PDM_g_num_t **faceToVtxB, */
/*   double **faceVtxCooB, */
/*   double **faceVtxEpsB */
/* //const double *pts, */
/* //double        n[3]*\/ */
/* ) */
/* { */
/* PDM_printf ("myRank : %d\n", myRank); */
/*   (*faceVtxCooA) = malloc(sizeof(double) * 5); */
/*   for (int i=0; i<5; i++) */
/* 	  (*faceVtxCooA)[i] = i * 2.5; */
/*  PDM_printf ("faceVtxCooA[2] : %e\n", (*faceVtxCooA)[2]); */

/*   return 0; */
/* } */

static int init_cas1(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{

  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 4;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = 3;
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 5.;
    (*faceVtxCooA)[3*j+1] = 4.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 5.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = 4;
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 5.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 3;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 4;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = 4;
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 5.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 3;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = 3;
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 5.;
    (*faceVtxCooB)[3*j+1] = 4.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 5.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

  }
  return 0;
}


/* CAS 2 : deux triangles identiques
-0.5,0;4,-1;2,3
-0.5,0;4,-1;2,3*/
static int init_cas2(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = -0.5;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = -1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = -0.5;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = -1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = -0.5;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = -1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = -0.5;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = -1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 3 : deux triangles avec une arete commune
-0.5,0;4,-1;3,2
-0.5,0;3,0;3,2*/
static int init_cas3(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = -0.5;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = -1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = -0.5;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 3.5;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = -0.5;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 3.5;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = -0.5;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = -1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 4 : deux triangles avec une arete commune et le dernier sommet dans l arrete
-0.5,0;4,-1;2,3
-0.5,0;3.5,0;2,3
*/
static int init_cas4(
  const int myRank,
  PDM_edges_get_t *get_t, PDM_g_num_t *nGvtxA, PDM_g_num_t *nGvtxB, PDM_g_num_t *maxGNEdgeA, PDM_g_num_t *maxGNEdgeB,
  int *nVtxA, PDM_g_num_t **faceToEdgeA, PDM_g_num_t **faceToVtxA, double **faceVtxCooA, double **faceVtxEpsA,
  int *nVtxB, PDM_g_num_t **faceToEdgeB, PDM_g_num_t **faceToVtxB, double **faceVtxCooB, double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = -0.5;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = -1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = -0.5;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 3.5;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = -0.5;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 3.5;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = -0.5;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = -1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  return 0;
}

/* CAS 5 : deux triangles avec une arete commune et le dernier sommet en dehors : croisement d aretes
0,0;4,0;2,3
0,0;5,2;2,3
*/
static int init_cas5(
  const int myRank,
  PDM_edges_get_t *get_t, PDM_g_num_t *nGvtxA, PDM_g_num_t *nGvtxB, PDM_g_num_t *maxGNEdgeA, PDM_g_num_t *maxGNEdgeB,
  int *nVtxA, PDM_g_num_t **faceToEdgeA, PDM_g_num_t **faceToVtxA, double **faceVtxCooA, double **faceVtxEpsA,
  int *nVtxB, PDM_g_num_t **faceToEdgeB, PDM_g_num_t **faceToVtxB, double **faceVtxCooB, double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 5.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 5.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  return 0;
}

/* CAS 6 : deux triangles partageant un morceau d arete
0,0;4,1;3,3
2,0.5;4,1;3,-1
*/
static int init_cas6(
  const int myRank,
  PDM_edges_get_t *get_t, PDM_g_num_t *nGvtxA, PDM_g_num_t *nGvtxB, PDM_g_num_t *maxGNEdgeA, PDM_g_num_t *maxGNEdgeB,
  int *nVtxA, PDM_g_num_t **faceToEdgeA, PDM_g_num_t **faceToVtxA, double **faceVtxCooA, double **faceVtxEpsA,
  int *nVtxB, PDM_g_num_t **faceToEdgeB, PDM_g_num_t **faceToVtxB, double **faceVtxCooB, double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 0.5;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = -1.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 0.5;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = -1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  return 0;
}

/* CAS 7 : le sommet d un des triangles et un point de l arete de la utre
0,0;4,0;2,4
3.5,1;6,2;5,4
*/
static int init_cas7(
  const int myRank,
  PDM_edges_get_t *get_t, PDM_g_num_t *nGvtxA, PDM_g_num_t *nGvtxB, PDM_g_num_t *maxGNEdgeA, PDM_g_num_t *maxGNEdgeB,
  int *nVtxA, PDM_g_num_t **faceToEdgeA, PDM_g_num_t **faceToVtxA, double **faceVtxCooA, double **faceVtxEpsA,
  int *nVtxB, PDM_g_num_t **faceToEdgeB, PDM_g_num_t **faceToVtxB, double **faceVtxCooB, double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 4.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 3.5;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 5.;
    (*faceVtxCooB)[3*j+1] = 4.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 3.5;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 5.;
    (*faceVtxCooA)[3*j+1] = 4.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 4.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  return 0;
}

/* CAS 8 : deux triangles inverses
0,0;4,0.5;2.5,5
2,-1;4,3;0.5,3
*/
static int init_cas8(
  const int myRank,
  PDM_edges_get_t *get_t, PDM_g_num_t *nGvtxA, PDM_g_num_t *nGvtxB, PDM_g_num_t *maxGNEdgeA, PDM_g_num_t *maxGNEdgeB,
  int *nVtxA, PDM_g_num_t **faceToEdgeA, PDM_g_num_t **faceToVtxA, double **faceVtxCooA, double **faceVtxEpsA,
  int *nVtxB, PDM_g_num_t **faceToEdgeB, PDM_g_num_t **faceToVtxB, double **faceVtxCooB, double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.5;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.5;
    (*faceVtxCooA)[3*j+1] = 5.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = -1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 0.5;
    (*faceVtxCooB)[3*j+1] = 3.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = -1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 0.5;
    (*faceVtxCooA)[3*j+1] = 3.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.5;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.5;
    (*faceVtxCooB)[3*j+1] = 5.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  return 0;
}

/* CAS 9 : deux triangles avec un sommet commun
0,2;0,0;3,1
3,1;5,2;5,0
*/
static int init_cas9(
  const int myRank,
  PDM_edges_get_t *get_t, PDM_g_num_t *nGvtxA, PDM_g_num_t *nGvtxB, PDM_g_num_t *maxGNEdgeA, PDM_g_num_t *maxGNEdgeB,
  int *nVtxA, PDM_g_num_t **faceToEdgeA, PDM_g_num_t **faceToVtxA, double **faceVtxCooA, double **faceVtxEpsA,
  int *nVtxB, PDM_g_num_t **faceToEdgeB, PDM_g_num_t **faceToVtxB, double **faceVtxCooB, double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 5.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 5.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 5.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 5.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  return 0;
}

/* CAS 10 : deux carres identiques mais un contenant 4 mailles MAIS NON!!! notre test ne doit avoir qu une maille par maillage!!!!!
un triangle et un quadranle avec A1OnEdgeB
0,0.8;2,2;2,0
0,0;0,2;2,2;2,0
*/
static int init_cas10(
  const int myRank,
  PDM_edges_get_t *get_t, PDM_g_num_t *nGvtxA, PDM_g_num_t *nGvtxB, PDM_g_num_t *maxGNEdgeA, PDM_g_num_t *maxGNEdgeB,
  int *nVtxA, PDM_g_num_t **faceToEdgeA, PDM_g_num_t **faceToVtxA, double **faceVtxCooA, double **faceVtxEpsA,
  int *nVtxB, PDM_g_num_t **faceToEdgeB, PDM_g_num_t **faceToVtxB, double **faceVtxCooB, double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 4;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.8;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 3;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 4;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 3;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.8;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

  }
  return 0;
}

/* CAS 11 : deux triangles avec A1=B1 mais A2 != B2 et A1 B1 A2 B2 alignes
4,0;8,0;6,2
4,0;0,0;2,2*/
static int init_cas11(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 8.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 8.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 12 : deux triangles avec A1=B1 mais A2 != B2 et A1 B1 A2 B2 alignes
4,0;8,0;6,2
4,0;6,0;6,2*/
static int init_cas12(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 8.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 8.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

static int init_cas13(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 8.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 8.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 6.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 6.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 14 : deux triangles dont un a une rete incluse ds l arete de l autre
0,0;4,0;2,2
1,0;3,0;2,-2*/
static int init_cas14(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 1.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = -2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 1.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = -2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}


/* CAS 15 : deux triangles (isA2B1 && isA1B2)
0,0;4,0;2,2
0,0;2,-2.4,0*/
static int init_cas15(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = -2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = -2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 16 : deux triangles (!isA1B1 && isA2B2) mais alignes
0,0;4,0;0,2
2,0;4,0;2,2*/
static int init_cas16(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 17 : deux triangles (!isA1B1 && isA2B2) mais alignes ds lautre sens pour que nNewPointsB == 2
2,0;4,0;2,2
0,0;4,0;0,2*/
static int init_cas17(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 18 : deux triangles (!isA2B1 && isA1B2) mais alignes et pour que nNewPointsA == 2
0,0;4,0;2,0
2,0;0,0;2,-2*/
static int init_cas18(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = -2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = -2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 19 : deux triangles (!isA2B1 && isA1B2) mais alignes et pour que nNewPointsB == 2
0,0;2,0;0,2
4,0;0,0;2,-2*/
static int init_cas19(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 2.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = -2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 2.;
    (*faceVtxCooA)[3*j+1] = -2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 2.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 2.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 20 : deux triangles deux aretes alignees mais aucun point commun et mA1B1<vB
0,0;3,0;0,1
1,0;4,0;1,1*/
static int init_cas20(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 1.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 1.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 1.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 1.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 21 : deux triangles deux aretes alignees mais aucun point commun et mA1B1<vB et a1 entre B1B2
1,0;3,0;1,1
0,0;4,0;0,1*/
static int init_cas21(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 1.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 1.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 1.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 1.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/* CAS 22 : un triangle tres ters plat avec un sommet de l autre
0,0;5,0.;5,0.00001
3,0.000005;4,-1;4,1*/
static int init_cas22(
  const int myRank,
  PDM_edges_get_t *get_t,
  PDM_g_num_t *nGvtxA,
  PDM_g_num_t *nGvtxB,
  PDM_g_num_t *maxGNEdgeA,
  PDM_g_num_t *maxGNEdgeB,
  int *nVtxA,
  PDM_g_num_t **faceToEdgeA,
  PDM_g_num_t **faceToVtxA,
  double **faceVtxCooA,
  double **faceVtxEpsA,
  int *nVtxB,
  PDM_g_num_t **faceToEdgeB,
  PDM_g_num_t **faceToVtxB,
  double **faceVtxCooB,
  double **faceVtxEpsB
)
{
  if (myRank == 0) {
    (*get_t) = PDM_EDGES_GET_FROM_A;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 0.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 5.;
    (*faceVtxCooA)[3*j+1] = 0.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 5.;
    (*faceVtxCooA)[3*j+1] = 1.e-16;
    (*faceVtxCooA)[3*j+2] = 0.;


    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 3.;
    (*faceVtxCooB)[3*j+1] = 1.e-16/2.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = -1.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 4.;
    (*faceVtxCooB)[3*j+1] = 1.;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
  else {
    (*get_t) = PDM_EDGES_GET_FROM_B;
    (*nGvtxA) = 3;
    (*nGvtxB) = 3;
    (*maxGNEdgeA) = (*nGvtxA);
    (*maxGNEdgeB) = (*nGvtxB);

    (*nVtxA) = (*nGvtxA);
    (*faceToEdgeA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceToVtxA) = malloc (sizeof(PDM_g_num_t) * (*nVtxA));
    (*faceVtxCooA) = malloc (sizeof(double) * 3 * (*nVtxA));
    (*faceVtxEpsA)  = malloc (sizeof(double) * (*nVtxA));

    for (int i = 0; i < (*nVtxA); i++) {
      (*faceToEdgeA)[i] = i+1;
      (*faceToVtxA)[i] = i+1;
      (*faceVtxEpsA)[i] = 1e-3;
    }

    int j;

    j = 0;
    (*faceVtxCooA)[3*j+0] = 3.;
    (*faceVtxCooA)[3*j+1] = 1.e-16/2.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = -1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooA)[3*j+0] = 4.;
    (*faceVtxCooA)[3*j+1] = 1.;
    (*faceVtxCooA)[3*j+2] = 0.;

    (*nVtxB) = (*nGvtxB);
    (*faceToEdgeB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceToVtxB) = malloc (sizeof(PDM_g_num_t) * (*nVtxB));
    (*faceVtxCooB) = malloc (sizeof(double) * 3 * (*nVtxB));
    (*faceVtxEpsB)  = malloc (sizeof(double) * (*nVtxB));

    for (int i = 0; i < (*nVtxB); i++) {
      (*faceToEdgeB)[i] = i+1;
      (*faceToVtxB)[i] = i+1;
      (*faceVtxEpsB)[i] = 1e-3;
    }

    j = 0;
    (*faceVtxCooB)[3*j+0] = 0.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 1;
    (*faceVtxCooB)[3*j+0] = 5.;
    (*faceVtxCooB)[3*j+1] = 0.;
    (*faceVtxCooB)[3*j+2] = 0.;

    j = 2;
    (*faceVtxCooB)[3*j+0] = 5.;
    (*faceVtxCooB)[3*j+1] = 1.e-16;
    (*faceVtxCooB)[3*j+2] = 0.;
  }
    return 0;
}

/*============================================================================
 * Type definitions
 *============================================================================*/

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

  int myRank;
  int numProcs;

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  assert (numProcs == 2);

  const double     vtxCarLengthTol = 1e-3;

  PDM_edges_get_t get_t;
  PDM_g_num_t nGvtxA;
  PDM_g_num_t nGvtxB;
  PDM_g_num_t maxGNEdgeA;
  PDM_g_num_t maxGNEdgeB;

  int nVtxA;
  PDM_g_num_t *faceToEdgeA;
  PDM_g_num_t *faceToVtxA;
  double      *faceVtxCooA;
  double      *faceVtxEpsA;

  int nVtxB;
  PDM_g_num_t *faceToEdgeB;
  PDM_g_num_t *faceToVtxB;
  double      *faceVtxCooB;
  double      *faceVtxEpsB;

  int cas = 1;

  if (cas == 1) {
    init_cas1(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 2) {
    init_cas2(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 3) {
    init_cas3(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 4) {
    init_cas4(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 5) {
    init_cas5(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 6) {
    init_cas6(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 7) {
    init_cas7(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 8) {
    init_cas8(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 9) {
    init_cas9(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 10) {
    init_cas10(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 11) {
    init_cas11(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 12) {
    init_cas12(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 13) {
    init_cas13(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 14) {
    init_cas14(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 15) {
    init_cas15(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 16) {
    init_cas16(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 17) {
    init_cas17(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 18) {
    init_cas18(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 19) {
    init_cas19(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 20) {
    init_cas20(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  else if (cas == 21) {
    init_cas21(
              myRank,
              &get_t,
              &nGvtxA,
              &nGvtxB,
              &maxGNEdgeA,
              &maxGNEdgeB,
              &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
              &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
  // double *f;
  // ret = t1(myRank, &f);
   // PDM_printf ("f[2] : %e\n", f[2]);
  else {
  init_cas22(
    myRank,
    &get_t,
    &nGvtxA,
    &nGvtxB,
    &maxGNEdgeA,
    &maxGNEdgeB,
    &nVtxA, &faceToEdgeA, &faceToVtxA, &faceVtxCooA, &faceVtxEpsA,
    &nVtxB, &faceToEdgeB, &faceToVtxB, &faceVtxCooB, &faceVtxEpsB);
  }
	PDM_edges_intersect_t *ei = PDM_edges_intersect_create (maxGNEdgeA,
							  maxGNEdgeB,
							  vtxCarLengthTol,
							  PDM_MPI_COMM_WORLD);

//  PDM_edges_intersect_dump(ei);
PDM_printf ("****  PDM_edges_intersect_create  OK  ******************************************************************\n");

  PDM_edges_intersect_poly_add (ei,
				nVtxA,
				faceToEdgeA,
				faceToVtxA,
				faceVtxCooA,
				faceVtxEpsA,
				nVtxB,
				faceToEdgeB,
				faceToVtxB,
				faceVtxCooB,
				faceVtxEpsB);
  PDM_edges_intersect_dump(ei);
PDM_printf ("****  PDM_edges_intersect_poly_add  OK  ******************************************************************\n");

  int *vtxAOnEdgeB = malloc(sizeof(int) * nVtxA);
  PDM_edges_intersect_res_t **vtxAOnEdgeBEir = malloc(sizeof(PDM_edges_intersect_res_t *) * nVtxA);
  for (int i = 0; i < nVtxA; i++) {
    vtxAOnEdgeB[i] = -1;
    vtxAOnEdgeBEir[i] = NULL;
  }

  int *vtxBOnEdgeA = malloc(sizeof(int) * nVtxB);
  PDM_edges_intersect_res_t **vtxBOnEdgeAEir = malloc(sizeof(PDM_edges_intersect_res_t *) * nVtxB);
  for (int i = 0; i < nVtxB; i++) {
    vtxBOnEdgeA[i] = -1;
    vtxBOnEdgeAEir[i] = NULL;
  }
  PDM_printf ("oNewPoints : \n  0 : PDM_EDGES_INTERSECT_POINT_NEW\n  1 : PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB\n  2 : PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA\n  3 : PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB\n");
  free(vtxAOnEdgeBEir);
  free(vtxBOnEdgeAEir);

  for (int iedgeA = 0; iedgeA < nVtxA; iedgeA++) {

    int n_intersect;
    PDM_edges_intersect_res_t **eir;
    for (int iedgeB = 0; iedgeB < nVtxB; iedgeB++) {
  		PDM_printf ("iedgeA : %d     iedgeB : %d", iedgeA, iedgeB);
  		eir = PDM_edges_intersect_get (ei,
  			                             get_t,
  			                             (PDM_g_num_t) (iedgeA + 1),
  			                             (PDM_g_num_t) (iedgeB + 1),
  			                             &n_intersect);
  		PDM_printf ("  ===>>   n_intersect : %d \n",n_intersect);
  		for (int iinter = 0; iinter < n_intersect; iinter++) {

  			PDM_printf ("iinter : %d \n", iinter);
  			PDM_line_intersect_t tIntersect;
  			PDM_g_num_t nGEdgeA;
  			PDM_g_num_t originEdgeA;
  			PDM_g_num_t endEdgeA;
  			int nNewPointsA;
  			PDM_edges_intersect_point_t *oNewPointsA;
  			PDM_g_num_t *linkA;
  			PDM_g_num_t *gNumA;
  			double *coordsA;
  			double *uA;

  			PDM_edges_intersect_res_data_get (eir[iinter],
  											                  PDM_EDGES_INTERSECT_MESHA,
  											                  &nGEdgeA,
  											                  &originEdgeA,
  											                  &endEdgeA,
  											                  &tIntersect,
  											                  &nNewPointsA,
  											                  &oNewPointsA,
  											                  &linkA,
  											                  &gNumA,
  											                  &coordsA,
  											                  &uA);

  			//if(tIntersect)
  			{
  			  PDM_printf ("**********************************************************************\n");
  			  PDM_printf ("******* MESH A *******************************************************\n");
  			  PDM_printf ("nGEdge : %d   ", nGEdgeA);
  			  PDM_printf ("originEdge : %d   ", originEdgeA);
  			  PDM_printf ("tIntersect : %d  (UNDEF=-1, NO=0, YES=1, ON_LINE=2)   ", tIntersect);
  			  PDM_printf ("nNewPointsA : %d\n", nNewPointsA);
  			  for (int ip=0; ip<nNewPointsA; ip++) {
  			  	PDM_printf ("oNewPoints[%d] : %d   ", ip, oNewPointsA[ip]);
  			  	PDM_printf ("link[%d] : %d   ", ip, linkA[ip]);
  			  	PDM_printf ("gNum[%d] : %d   ", ip, gNumA[ip]);
  			  	PDM_printf ("u[%d] : %e   ", ip, uA[ip]);
  			  	PDM_printf ("coords : %ld   \n",coordsA);
  			  	PDM_printf ("coords[%d] : %e    ", 3*ip+0, coordsA[3*ip+0]);
  			  	PDM_printf ("coords[%d] : %e    ", 3*ip+1, coordsA[3*ip+1]);
  			  	PDM_printf ("coords[%d] : %e  \n", 3*ip+2, coordsA[3*ip+2]);
  			  }
  			  PDM_printf ("**********************************************************************\n");
  			}
  		}
      if(eir != NULL){
        free(eir);
      }
  	}
  }

  for (int iedgeA = 0; iedgeA < nVtxA; iedgeA++) {
	  int n_intersect;
    PDM_edges_intersect_res_t **eir;
    for (int iedgeB = 0; iedgeB < nVtxB; iedgeB++) {
	  	PDM_printf ("iedgeA : %d     iedgeB : %d", iedgeA, iedgeB);
	  	eir = PDM_edges_intersect_get (ei,
	  		                             get_t,
	  		                             (PDM_g_num_t) (iedgeA + 1),
	  		                             (PDM_g_num_t) (iedgeB + 1),
	  		                             &n_intersect);
	  	PDM_printf ("  ===>>   n_intersect : %d \n",n_intersect);

	  	for (int iinter = 0; iinter < n_intersect; iinter++) {

	  		PDM_printf ("iinter : %d \n", iinter);
	  		PDM_line_intersect_t tIntersect;
	  		PDM_g_num_t nGEdgeB;
	  		PDM_g_num_t originEdgeB;
	  		PDM_g_num_t endEdgeB;
	  		int nNewPointsB;
	  		PDM_edges_intersect_point_t *oNewPointsB;
	  		PDM_g_num_t *linkB;
	  		PDM_g_num_t *gNumB;
	  		double *coordsB;
	  		double *uB;

	  		PDM_edges_intersect_res_data_get (eir[iinter],
	  										                  PDM_EDGES_INTERSECT_MESHB,
	  										                  &nGEdgeB,
	  										                  &originEdgeB,
	  										                  &endEdgeB,
	  										                  &tIntersect,
	  										                  &nNewPointsB,
	  										                  &oNewPointsB,
	  										                  &linkB,
	  										                  &gNumB,
	  										                  &coordsB,
	  										                  &uB);

	  		//if(tIntersect)
	  		{
	  		  PDM_printf ("**********************************************************************\n");
	  		  PDM_printf ("******* MESH B *******************************************************\n");
	  		  PDM_printf ("nGEdge : %d   ", nGEdgeB);
	  		  PDM_printf ("originEdge : %d   ", originEdgeB);
	  		  PDM_printf ("tIntersect : %d  (UNDEF=-1, NO=0, YES=1, ON_LINE=2)   ", tIntersect);
	  		  PDM_printf ("nNewPointsB : %d\n", nNewPointsB);
	  		  for (int ip=0; ip<nNewPointsB; ip++) {
	  		  	PDM_printf ("oNewPoints[%d] : %d   ", ip, oNewPointsB[ip]);
	  		  	PDM_printf ("link[%d] : %d   ", ip, linkB[ip]);
	  		  	PDM_printf ("gNum[%d] : %d   ", ip, gNumB[ip]);
	  		  	PDM_printf ("u[%d] : %e   ", ip, uB[ip]);
	  		  	PDM_printf ("coords : %d\n", coordsB);
	  		  	PDM_printf ("coords[%d] : %e    ", 3*ip+0, coordsB[3*ip+0]);
	  		  	PDM_printf ("coords[%d] : %e    ", 3*ip+1, coordsB[3*ip+1]);
	  		  	PDM_printf ("coords[%d] : %e  \n", 3*ip+2, coordsB[3*ip+2]);
	  		  }
	  		  PDM_printf ("**********************************************************************\n");
	  		}
	  	}
      if(eir != NULL){
        free(eir);
      }
	  }
  }


  for (int iedgeA = 0; iedgeA < nVtxA; iedgeA++) {
   PDM_printf ("vtxAOnEdgeB[%d] : %d \n", iedgeA, vtxAOnEdgeB[iedgeA]);
  }
  for (int iedgeB = 0; iedgeB < nVtxB; iedgeB++) {
   PDM_printf ("vtxBOnEdgeA[%d] : %d \n", iedgeB, vtxBOnEdgeA[iedgeB]);
  }

  free (faceToEdgeA);
  free (faceToVtxA);
  free (faceVtxCooA);
  free (faceVtxEpsA);

  free (faceToEdgeB);
  free (faceToVtxB);
  free (faceVtxCooB);
  free (faceVtxEpsB);

  free(vtxBOnEdgeA);
  free(vtxAOnEdgeB);

  PDM_edges_intersect_free (ei);
  PDM_MPI_Finalize ();

  return 0;

}
