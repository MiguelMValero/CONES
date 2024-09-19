#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_ho_bezier.h"


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init(&argc, &argv);

  double point_coord[3] = {
    -1.7516283337977998e+03, -2.4366417578898277e+03, -1.6904612902513327e+03
  };

  double node_coord[3*10] = {
    -8.0296495277832003e+02, -1.3890805515241000e+03, -9.6136786496847003e+02,
    -7.7952673129348398e+02, -1.3930219640567957e+03, -9.9049326381078606e+02,
    -7.5524779978206868e+02, -1.3962348401410800e+03, -1.0191639157580662e+03,
    -7.3017588396461997e+02, -1.3984962345714000e+03, -1.0470787613708001e+03,
    -7.7634577363059907e+02, -1.4168863752472719e+03, -9.6070418259120015e+02,
    -7.5232512212535721e+02, -1.4203639364413416e+03, -9.8951284758796032e+02,
    -7.2749705788034714e+02, -1.4229139556528405e+03, -1.0177074408379543e+03,
    -7.4877630148282299e+02, -1.4439180852325642e+03, -9.5946534676684917e+02,
    -7.2426425149072077e+02, -1.4467704408092804e+03, -9.8784935035719172e+02,
    -7.2039530384498005e+02, -1.4698580869662001e+03, -9.5755153831604002e+02
  };

  double proj_coord[3];
  double uvw[3];
  PDM_ho_bezier_triangle_location(3,
                                  10,
                                  node_coord,
                                  point_coord,
                                  proj_coord,
                                  uvw);

  if (0) {
    PDM_vtk_write_point_cloud("point.vtk",
                              1,
                              point_coord,
                              NULL,
                              NULL);

    int connec[10] = {1, 4, 10, 2, 3, 7, 9, 8, 5, 6};
    PDM_vtk_write_std_elements_ho("bezier_triangle.vtk",
                                  3,
                                  10,
                                  node_coord,
                                  NULL,
                                  PDM_MESH_NODAL_TRIAHO_BEZIER,
                                  1,
                                  connec,
                                  NULL,
                                  0,
                                  NULL,
                                  NULL);
  }
  PDM_MPI_Finalize();

  return 0;
}
