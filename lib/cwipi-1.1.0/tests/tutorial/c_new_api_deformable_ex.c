/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_generate_mesh.h"

/*----------------------------------------------------------------------
 *
 * Main : advanced test : Deformable (code 2)
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  // Initialize MPI :
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // Initialize CWIPI :
  // Use CWP_Init for code2 written in C.
  // ------------------------------------------------------- To fill in
  int n_code = 1;
  const char  **code_name  = malloc(sizeof(char *) * n_code);
  MPI_Comm     *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

  code_name[0]      = "code2";

  // ---------------------------------------------------- End To fill in

  // Create the coupling :
  // In this tutorial the mesh will be deformed over the coupling
  // itreations. Those iterations mimic solver steps. Use CWP_Cpl_create
  // to couple code1 with code2. A rectangular mesh is used with triangle
  // elements. Operate the localization with the octree method.
  // ------------------------------------------------------- To fill in
  int n_part = 1;
  CWP_UNUSED(n_part); // For remove compilation warning

  const char  *coupling_name     = "coupling";
  CWP_UNUSED(coupling_name); // For remove compilation warning

  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  coupled_code_name[0] = "code1";

  // ---------------------------------------------------- End To fill in

  // Set coupling visualisation:
  // Use CWP_Visu_set to output ASCII Ensight format files
  // at each iteration step. Open the CHR.case file in
  // paraview to visualize the deformed mesh and exchanged fields.
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Create mesh :
  // A ParaDiGM library function is used to generate the rectangular mesh.
  // PDM_MPI_mpi_2_pdm_mpi_comm is used to convert the MPI communicator into
  // the expected ParaDiGM format. It is the users responsability to
  // free arrays from _simplified mesh generation functions. In a real life
  // coupling case here a user generated mesh would be read/loaded/given.
  int     n_vtx = 0;
  int     n_elt = 0;
  CWP_UNUSED(n_elt); // For remove compilation warning

  double *coords      = NULL;
  int    *elt_vtx_idx = NULL;
  int    *elt_vtx     = NULL;
  // TO UNCOMMENT -->>
  // PDM_generate_mesh_rectangle_simplified(PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[0]),
  //                                        10,
  //                                        &n_vtx,
  //                                        &n_elt,
  //                                        &coords,
  //                                        &elt_vtx_idx,
  //                                        &elt_vtx);
  // <<--

  // Interations :
  // At each iteration the mesh coordinates and the exchanged fields are modified.
  const char *send_field_name = "chinchilla";
  CWP_UNUSED(send_field_name); // For remove compilation warning

  const char *recv_field_name = "girafe";
  CWP_UNUSED(recv_field_name); // For remove compilation warning

  int         n_components    = 1;
  CWP_UNUSED(n_components); // For remove compilation warning
  
  double     *send_field_data = malloc(sizeof(double) * n_vtx);
  CWP_UNUSED(send_field_data); // For remove compilation warning

  double     *recv_field_data = malloc(sizeof(double) * n_vtx);
  CWP_UNUSED(recv_field_data); // For remove compilation warning

  const int    itdeb = 1;
  const int    itend = 10;
  const double freq  = 0.20;
  const double ampl  = 0.012;
  const double phi   = 0.1;
  double       ttime = 0.0;
  double       dt    = 0.1;

  double omega = 2.0*acos(-1.0)*freq;

  for (int it = itdeb; it <= itend; it ++) {

    ttime = (it-itdeb)*dt;

    // Start time step
    // ------------------------------------------------------- To fill in

    // ---------------------------------------------------- End To fill in

    // Deform mesh and set send field :
    // The field that will be sent, "chinchilla", is set to the value of
    // the z-coordinate of the mesh nodes.
    for (int i = 0; i < n_vtx; i++) {
      coords[3 * i + 2]  = ampl * (coords[3 * i]*coords[3 * i]+coords[1 + 3 * i]*coords[1 + 3 * i])*cos(omega*ttime+phi);
      send_field_data[i] = coords[3 * i + 2];
    }

    // First iteration :
    // Since array pointers are given to CWIPI, they are only set
    // at the first iteration.
    if (it == itdeb) {

      // Set the mesh :
      // Here the appropriate CWIPI functions should be called to
      // set the arrays difining the mesh and finalize it.
      // ------------------------------------------------------- To fill in

      // ---------------------------------------------------- End To fill in

      // Set fields :
      // code2 sends the "chinchilla" field and receives "girafe".
      // ------------------------------------------------------- To fill in

      // ---------------------------------------------------- End To fill in

      // Set tolerance to 10%.
      // ------------------------------------------------------- To fill in

      // ---------------------------------------------------- End To fill in

    } else {

      // Find the correct CWIPI function to inform CWIPI that the
      // the time step has been updated.
      // ------------------------------------------------------- To fill in

      // ---------------------------------------------------- End To fill in
    }

    // Compute interpolation weights :
    // ------------------------------------------------------- To fill in

    // ---------------------------------------------------- End To fill in

    // Exchange field values between codes :
    // Use the CWIPI exchange functions similar to the MPI ones
    // for code2 to send "chinchilla" and receive "girafe".
    // ------------------------------------------------------- To fill in

    // ---------------------------------------------------- End To fill in

    // Check interpolation :
    // For the receiving field "girafe", check if all  points have been located.
    // ------------------------------------------------------- To fill in
    int n_uncomputed_tgts = 0;
    const int *uncomputed_tgts = NULL;

    if (n_uncomputed_tgts != 0) {

    }

    // ---------------------------------------------------- End To fill in

    PDM_UNUSED(n_uncomputed_tgts);
    PDM_UNUSED(uncomputed_tgts);

    // End time step
    // ------------------------------------------------------- To fill in

    // ---------------------------------------------------- End To fill in

  } // end interations

  // Delete fields :
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Delete Mesh :
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Delete the coupling :
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // free
  if (intra_comm        != NULL) free(intra_comm);
  if (code_name         != NULL) free(code_name);
  if (coupled_code_name != NULL) free(coupled_code_name);
  if (coords            != NULL) free(coords);
  if (elt_vtx_idx       != NULL) free(elt_vtx_idx);
  if (elt_vtx           != NULL) free(elt_vtx);
  if (send_field_data   != NULL) free(send_field_data);
  if (recv_field_data   != NULL) free(recv_field_data);

  // Finalize CWIPI :
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;
}
