/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

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
#include <unistd.h>
#include <assert.h>

#include "cwp.h"
#include "pdm_io.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "client_server/client.h"

#include "cwp_priv.h"
#include "cwipi_config.h"

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

/*=============================================================================
 * Util functions
 *============================================================================*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c     Filename of the server configuration file.\n\n"
     "  -h     This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int            argc,
 char         **argv,
 char         **config  // filename for server ip adresses + ports
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *config = argv[i];
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}

/*=============================================================================
 * Main
 *============================================================================*/

int
main
(
 int   argc,
 char *argv[]
)
{
  // default
  char *config     = NULL;

  _read_args(argc,
             argv,
             &config);

  // mpi
  int rank;
  int comm_world_size;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &comm_world_size);

  if (config == NULL) {
    if (rank == 0) {
      config = (char *) "client_new_api_cellface_set_block_o/code1/cwp_config_srv.txt";
    }
    else {
      config = (char *) "client_new_api_cellface_set_block_o/code2/cwp_config_srv.txt";
    }
  }

  assert (comm_world_size == 2);

  // launch server

  if (rank == 0) {
    system("mkdir -p client_new_api_cellface_set_block_o/code1");
    system("mkdir -p client_new_api_cellface_set_block_o/code2");
    system("rm -f ./client_new_api_cellface_set_block_o/code1/cwp_config_srv.txt");
    system("rm -f ./client_new_api_cellface_set_block_o/code2/cwp_config_srv.txt");
    system("mpirun -n 1 cwp_server -cn code0 -p 49104 49104 -c \"client_new_api_cellface_set_block_o/code1/cwp_config_srv.txt\" : -n 1  cwp_server -cn code1 -p 49105 49105 -c \"client_new_api_cellface_set_block_o/code2/cwp_config_srv.txt\" &");
  }

  while (access(config, R_OK) != 0) {
    sleep(1);
  }
  sleep(5);

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


  // Init
  const char *code_name = NULL;
  int id_code = 0;

  if (rank == 0) {
    id_code = 0;
    code_name = "code1_cell_faces";
  }

  if (rank == 1) {
    id_code = 1;
    code_name = "code2";
  }


  CWP_Status_t is_coupled_rank = CWP_STATUS_ON;

  MPI_Comm intra_comm;
  MPI_Comm_split(comm, id_code, rank, &intra_comm);

  CWP_client_Init(intra_comm,
                  config,
                  code_name,
                  is_coupled_rank);

  char cpl_id1[] = "cpl_code1_code2";

  if (rank == 0) {
    CWP_client_Cpl_create("code1_cell_faces", cpl_id1, "code2", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                          CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE, 1,
                          CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (rank == 1) {
    CWP_client_Cpl_create("code2", cpl_id1, "code1_cell_faces", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                          CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE, 1,
                          CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }

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

  FILE *mesh_file;
  mesh_file = fopen(CWP_MESH_DIR"mesh_poly_d1", "r"); // WARNING: adapt depending on where client is launched

  fscanf(mesh_file, "%d %d %d %d %d %d",
         &dimension,
         &n_vtx,
         &n_face,
         &n_elmt,
         &lface_connec,
         &lcell_connec);

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

    CWP_client_Visu_set("code1_cell_faces", cpl_id1, 1.0, CWP_VISU_FORMAT_ENSIGHT, "binary");

    CWP_g_num_t *global_num_vtx = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_vtx);
    for (int i = 0; i < n_vtx; i++) {
      global_num_vtx[i] = i + 1;
    }

    CWP_client_Mesh_interf_vtx_set("code1_cell_faces", cpl_id1, 0, n_vtx, coords, global_num_vtx);

    CWP_g_num_t *global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elmt);
    for (int i = 0; i < n_elmt; i++) {
      global_num[i] = i + 1;
    }

    CWP_client_Mesh_interf_from_cellface_set("code1_cell_faces",
                                             cpl_id1,
                                             0,
                                             n_elmt,
                                             cell_face_idx,
                                             cell_face,
                                             n_face,
                                             face_vtx_idx,
                                             face_vtx,
                                             global_num);

    // free
    free(global_num_vtx);
    free(global_num);

    CWP_client_Mesh_interf_del("code1_cell_faces", cpl_id1);
  }

  if (rank == 0) {
    CWP_client_Cpl_del("code1_cell_faces", cpl_id1);
  }

  if (rank == 1) {
    CWP_client_Cpl_del("code2", cpl_id1);
  }

  fflush(stdout);

  CWP_client_Finalize();

  MPI_Comm_free(&intra_comm);

  MPI_Finalize();

  free(coords       );
  free(face_vtx_idx );
  free(face_vtx     );
  free(cell_face_idx);
  free(cell_face    );

  return 0;
}
