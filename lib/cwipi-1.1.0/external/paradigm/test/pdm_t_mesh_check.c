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

#include "pdm_mesh_check.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"



int main(int argc, char *argv[])
{
  PDM_MPI_Init (&argc, &argv);

  const int n_i = 3;
  const int n_j = 3;

  PDM_g_num_t n_vtx      = n_i * n_j;
  PDM_g_num_t l_face_vtx = 8;

  double *coords = (double *) malloc(sizeof(double) * 3 * n_vtx);
  int k = 0;
  for (int j = 0; j < n_j; j++) {
    for (int i = 0; i < n_i; i++) {
      if (k < 3*n_vtx) {
        coords[k++] = (double) i;
        coords[k++] = (double) j;
        coords[k++] = 0.;
      }
    }
  }


  PDM_g_num_t *face_vtx = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * l_face_vtx);
  face_vtx[0] = 1;
  face_vtx[1] = 2;
  face_vtx[2] = 5;
  face_vtx[3] = 4;
  // face_vtx[4] = 2;
  // face_vtx[5] = 3;
  // face_vtx[6] = 6;
  // face_vtx[7] = 5;
  face_vtx[4] = 4;
  face_vtx[5] = 5;
  face_vtx[6] = 8;
  face_vtx[7] = 7;

  int n_holes;
  PDM_mesh_check_unconnected_vertex(&n_vtx,
                                    &l_face_vtx,
                                    face_vtx,
                                    coords,
                                    &n_holes);

  printf(PDM_FMT_G_NUM" vertices, %d holes\n", n_vtx, n_holes);

  free(face_vtx);
  free(coords);

  PDM_MPI_Finalize();

  return 0;
}
