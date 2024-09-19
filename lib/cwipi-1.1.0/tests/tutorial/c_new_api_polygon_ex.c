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

#include "cwp.h"

/*----------------------------------------------------------------------
 *
 * Main : base test : Polygon
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  // Initialize MPI
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // This test mimics the coupling between 2 codes running each
  // on one processor.
  assert (n_rank == 2);

  // Initialize CWIPI :
  // Use CWP_Init for code1 running on MPI rank 0 and code2
  // running on MPU rank 1.
  // ------------------------------------------------------- To fill in
  int n_code = 1;
  const char **code_name = malloc(sizeof(char *) * n_code);

  int I_am_code1 = 0;
  int I_am_code2 = 0;

  if (i_rank == 0) {
    code_name[0] = "code1";
    I_am_code1   = 1;
  }

  if (i_rank == 1) {
    code_name[0] = "code2";
    I_am_code2   = 1;
  }

  // ---------------------------------------------------- End To fill in

  // Create the coupling :
  // Use CWP_Cpl_create to couple code1 with code2 on a surface
  // interface. Operate the localization with the octree method.
  // In this tutorial, one coupling iteration is done on a mesh
  // partitionned over all processors. This induces that the mesh
  // is static.
  // ------------------------------------------------------- To fill in
//  int n_part = 1;
//  const char  *coupling_name     = "code1_code2";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);

  if (I_am_code1) {
    coupled_code_name[0] = "code2";
  }

  if (I_am_code2) {
    coupled_code_name[0] = "code1";
  }

  // ---------------------------------------------------- End To fill in

  // Set coupling visualisation:
  // Use CWP_Visu_set to output ASCII Ensight format files
  // at each iteration step.
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Set the mesh vertices coordinates :
  // Use CWP_Mesh_interf_vtx_set to set the mesh vertex coordinates,
  // no global numbering of the vertices will be given. In this
  // simple setting, there is only one partition per processor.
  // ------------------------------------------------------- To fill in
  int    n_vtx  = 11;
  double coords[33] = {0,0,0,  1,0,0,  2,0,0,  3,0,0,  0,1,0,  2,1,0,
                       3,1,0,  1,2,0,  0,3,0,  2,3,0,  3,3,0};

  // ---------------------------------------------------- End To fill in

  // Set the mesh polygons connectivity :
  // Use CWP_Mesh_interf_block_add to create a block of
  // of polygons. Choose the correct CWIPI function
  // to set a polygonal mesh, no need to give the elements
  // global numbering.
  // ------------------------------------------------------- To fill in

//  int n_elts = 5;
//  int connec_idx[6] = {0,3,7,11,16,21};
//  int connec[21]    = {1,2,5,   3,4,7,6,   5,8,10,9   ,5,2,3,6,8,   6,7,11,10,8};

  // ---------------------------------------------------- End To fill in

  // Finalize mesh :
  // Use the correct CWIPI function to generate the
  // mesh global numbering and the underlying mesh data structure
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Create and set the field values :
  // Use CWP_Field_create and CWP_Field_data_set to create and set
  // a field onto the mesh. code1 will send its field which code2
  // will receive. The field is located at the mesh nodes.
  // There is only one mesh partition in this tutorial. Activate
  // visualization for this field if you wish it to be in the
  // Ensight file.
  // Do not forget to begin the time step AFTER creating the fields,
  // but BEFORE setting the fields data!
  // ------------------------------------------------------- To fill in
//  const char *field_name      = "a super fancy field";
  int         n_components    = 1;
  double     *send_field_data = malloc(sizeof(double) * n_vtx * n_components);
  double     *recv_field_data = malloc(sizeof(double) * n_vtx * n_components);

  // Create the field :
  //...

  // Begin the time step :
  //...

  // Set the field pointers :
  // for code1
  if (I_am_code1) {

    for (int i = 0; i < n_vtx; i++) {
      send_field_data[i] = coords[3*i];
    }

    //...
  }

  // for code2
  if (I_am_code2) {
    //...
  }
  // ---------------------------------------------------- End To fill in

  // Compute interpolation weights :
  // Choose the two correct CWIPI functions to set the geometric
  // tolerance to 10% of an element size for point localisation
  // and to compute the interpolation weights.
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Exchange field values between codes :
  // Use the CWIPI exchange functions similar to the MPI ones
  // for code1 to send a field and code2 to receive that field.
  // The reason for the double if statements is that the exchange
  // operation is non-blocking. That means that work can be done
  // before calling the wait function where the field need to be
  // used or changed.
  // ------------------------------------------------------- To fill in
  if (I_am_code1) {

  }

  // for code2
  if (I_am_code2) {

  }

  // for code1
  if (I_am_code1) {

  }

  // for code2
  if (I_am_code2) {

  }
  // ---------------------------------------------------- End To fill in

  // End the time step :
  // ------------------------------------------------------- To fill in
  // ---------------------------------------------------- End To fill in

  // Check interpolation :
  // The field that has been sent will be interpolated on the vertices
  // of the mesh of the receiving code. Depending on the geometric
  // tolerance and the interpolation algorithm used, some points might
  // not be located. Users often want all points on the receiving mesh
  // to have a field value. Thus, for the receiving code, check if there
  // are vertices for which the interpolation has been unsuccessful.
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Delete field :
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Delete Mesh :
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Delete the coupling :
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // free
  if (send_field_data != NULL) free(send_field_data);
  if (recv_field_data != NULL) free(recv_field_data);
  free(code_name);
  free(coupled_code_name);

  // Finalize CWIPI :
  // ------------------------------------------------------- To fill in

  // ---------------------------------------------------- End To fill in

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;
}
