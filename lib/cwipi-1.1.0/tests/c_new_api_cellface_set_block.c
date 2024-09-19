/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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
#include <string.h>
#include <math.h>

#include "cwp.h"
#include "cwipi_config.h"

/*----------------------------------------------------------------------
 *
 * Read mesh dimension
 *
 * parameters:
 *   f                   <-- Mesh file
 *   dimension           --> Dimension
 *   n_vtx             --> number of vertices
 *   n_elmt           --> number of elements
 *   nConnecVertex       --> size of connectivity
 *---------------------------------------------------------------------*/

static int read_mesh_dim(FILE *f,
                         int *dimension,
                         int *n_vtx,
                         int *n_face,
                         int *n_elt,
                         int *lface_connec,
                         int *lcell_connec) {
  int r;
  r = fscanf(f, "%d %d %d %d %d %d",
             dimension,
             n_vtx,
             n_face,
             n_elt,
             lface_connec,
             lcell_connec);
  if (r == EOF)
    return 0;
  else return 1;
}


/*----------------------------------------------------------------------
 *
 * Read mesh
 *
 * parameters:
 *   f                   <-- Mesh file
 *   dimension           --> Dimension
 *   n_vtx               <-- number of vertices
 *   n_elmt              <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/

static int read_mesh(FILE    *f,
                     int     dimension,
                     int     n_vtx,
                     int     n_face,
                     int     n_elt,
                     int     lface_connec,
                     int     lcell_connec,
                     double *coords,
                     int    *face_vtx_idx,
                     int    *face_vtx,
                     int    *cell_face_idx,
                     int    *cell_face) {
  int i, j, r;

  // Read coordinates
  for (i = 0 ; i < n_vtx ; i++) {
    for (j = 0 ; j < dimension ; j++) {
      r = fscanf(f, "%lf", coords + i * dimension + j);
      if (r == EOF)
        return EXIT_FAILURE;
    }
  }

  // Read face -> vertex connectivity index
  for (i = 0 ; i < n_face + 1 ; i++) {
    r = fscanf(f, "%d", face_vtx_idx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read face -> vertex connectivity
  for (i = 0 ; i < lface_connec ; i++) {
    r = fscanf(f, "%d", face_vtx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity index
  for (i = 0 ; i < n_elt + 1 ; i++) {
    r = fscanf(f, "%d", cell_face_idx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity
  for (i = 0 ; i < lcell_connec ; i++) {
    r = fscanf(f, "%d", cell_face + i);
    //if(cell_face[i]<0) printf("cell_face[%i] %i\n",i,cell_face[i]);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


/*----------------------------------------------------------------------
 *
 * Main : linear coupling test
 *
 *---------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  FILE *mesh_file;
  mesh_file = fopen(CWP_MESH_DIR"mesh_poly_d1", "r");

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);


  int n_partition = 0;
  const int two = 2;
  while (two * pow(n_partition, two) < comm_world_size) n_partition++;

  int n2 = (int) (two * pow(n_partition, two));

  if (n2 != comm_world_size) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  char *src_name = (char *) malloc(sizeof(char) * (strlen(__FILE__) + 1));
  strcpy(src_name, __FILE__);
  char *src_base_name = NULL;
  src_base_name = strrchr(src_name, '.');
  if (src_base_name != NULL)
    *src_base_name = '\0';
  src_base_name = NULL;
  src_base_name = strrchr(src_name, '/');
  if (src_base_name != NULL)
    src_base_name += 1;
  else
    src_base_name = src_name;

  if (rank == 0)
    printf("\nSTART: %s\n", src_base_name);


  /* Initialization
   * -------------- */

  int n_code = 0;
  const char **code_names = NULL;
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  if (rank == 0) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1_cell_faces";
  }

  if (rank == 1) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code2";
  }


  MPI_Comm *local_comm = malloc(sizeof(MPI_Comm) * n_code);

  printf("CWIPI Initialization rank %i\n", rank);
  CWP_Init(MPI_COMM_WORLD,
           n_code,
           code_names,
           is_active_rank,
           local_comm);


  /* Output redirection
   * ------------------ */

  int current_rank;
  int local_comm_size;

  for (int i = 0 ; i < n_code ; i++) {
    MPI_Comm_rank(local_comm[i], &current_rank);
    MPI_Comm_size(local_comm[i], &local_comm_size);
    printf("Size of local_comm[%i]=%i et rang du proc=%i.\n", i, local_comm_size, current_rank);
  }

  /* Finalize
   * -------- */

  char cpl_id1[] = "cpl_code1_code2";

  printf("Coupling creation\n");

  if (rank == 0) {
    CWP_Cpl_create("code1_cell_faces", cpl_id1, "code2", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE, 1,
                   CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (rank == 1) {
    CWP_Cpl_create("code2", cpl_id1, "code1_cell_faces", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE, 1,
                   CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }
  printf("Coupling created %i\n", current_rank);

  /* Building of the local mesh */

  int dimension    = 3;             // Dimension of the space
  int n_vtx        = 0;             // Number of points in the mesh
  int n_face       = 0;             // Number of face
  int n_elmt       = 0;             // Number of cells
  int lface_connec = 0;
  int lcell_connec = 0;

  double *coords        = NULL;         // Coordinates of the points
  int    *face_vtx_idx  = NULL;
  int    *face_vtx      = NULL;
  int    *cell_face_idx = NULL;
  int    *cell_face     = NULL;

  if (rank == 0)
    printf("        Read mesh\n");

  read_mesh_dim(mesh_file, &dimension, &n_vtx, &n_face, &n_elmt, &lface_connec, &lcell_connec);

  coords        = (double *) malloc(dimension * n_vtx * sizeof(double));
  face_vtx_idx  = (int    *) malloc((n_face + 1)      * sizeof(int   ));
  face_vtx      = (int    *) malloc(lface_connec      * sizeof(int   ));
  cell_face_idx = (int    *) malloc((n_elmt + 1)      * sizeof(int   ));
  cell_face     = (int    *) malloc(lcell_connec      * sizeof(int   ));

  read_mesh(mesh_file,
            dimension,
            n_vtx,
            n_face,
            n_elmt,
            lface_connec,
            lcell_connec,
            coords,
            face_vtx_idx,
            face_vtx,
            cell_face_idx,
            cell_face);

  fclose(mesh_file);

  if (rank == 0) {

    printf("vtx_set\n");
    CWP_Mesh_interf_vtx_set("code1_cell_faces", cpl_id1, 0, n_vtx, coords, NULL);

    printf("Cell_face Add and Setting\n");
    CWP_Mesh_interf_from_cellface_set("code1_cell_faces",
                                      cpl_id1,
                                      0,
                                      n_elmt,
                                      cell_face_idx,
                                      cell_face,
                                      n_face,
                                      face_vtx_idx,
                                      face_vtx,
                                      NULL);


    printf("Interface Mesh deletion\n");
    CWP_Mesh_interf_del("code1_cell_faces", cpl_id1);
    printf("Interface Mesh deleted\n");
  }


  fflush(stdout);

  CWP_Finalize();

  MPI_Finalize();

  free(src_name);
  free(local_comm);
  free(code_names);
  free(coords       );
  free(face_vtx_idx );
  free(face_vtx     );
  free(cell_face_idx);
  free(cell_face    );

  return 0;
}
