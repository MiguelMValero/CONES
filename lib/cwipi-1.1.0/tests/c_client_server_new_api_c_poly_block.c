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
 * Read mesh dimension
 *
 * parameters:
 *   f                   <-- Mesh file
 *   dimension           --> Dimension
 *   nvertex             --> number of vertices
 *   nElements           --> number of elements
 *   nConnecVertex       --> size of connectivity
 *---------------------------------------------------------------------*/

static int read_mesh_dim(FILE *f,
                         int *dimension,
                         int *nVertex,
                         int *nFace,
                         int *nElt,
                         int *lFaceConnec,
                         int *lCellConnec) {
  int r;
  r = fscanf(f, "%d %d %d %d %d %d",
             dimension,
             nVertex,
             nFace,
             nElt,
             lFaceConnec,
             lCellConnec);
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
 *   nvertex             <-- number of vertices
 *   nElements           <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/

static int read_mesh(FILE *f,
                     int dimension,
                     int nVertex,
                     int nFace,
                     int nElt,
                     int lFaceConnec,
                     int lCellConnec,
                     double *coords,
                     int *faceVertexIdx,
                     int *faceVertex,
                     int *cellFaceIdx,
                     int *cellFace) {
  int i, j, r;

  // Read coordinates
  for (i = 0 ; i < nVertex ; i++) {
    for (j = 0 ; j < dimension ; j++) {
      r = fscanf(f, "%lf", coords + i * dimension + j);
      if (r == EOF)
        return EXIT_FAILURE;
    }
  }

  // Read face -> vertex connectivity index
  for (i = 0 ; i < nFace + 1 ; i++) {
    r = fscanf(f, "%d", faceVertexIdx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read face -> vertex connectivity
  for (i = 0 ; i < lFaceConnec ; i++) {
    r = fscanf(f, "%d", faceVertex + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity index
  for (i = 0 ; i < nElt + 1 ; i++) {
    r = fscanf(f, "%d", cellFaceIdx + i);
    if (r == EOF)
      return EXIT_FAILURE;
  }

  // Read cell -> face connectivity
  for (i = 0 ; i < lCellConnec ; i++) {
    r = fscanf(f, "%d", cellFace + i);
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
      config = (char *) "client_new_api_c_poly_block_o/code1/cwp_config_srv.txt";
    }
    else {
      config = (char *) "client_new_api_c_poly_block_o/code2/cwp_config_srv.txt";
    }
  }

  assert (comm_world_size == 2);

  // launch server

  if (rank == 0) {
    system("mkdir -p client_new_api_c_poly_block_o/code1");
    system("mkdir -p client_new_api_c_poly_block_o/code2");
    system("rm -f ./client_new_api_c_poly_block_o/code1/cwp_config_srv.txt");
    system("rm -f ./client_new_api_c_poly_block_o/code2/cwp_config_srv.txt");
    system("mpirun -n 1 cwp_server -cn code1 -p 48102 48102 -c \"client_new_api_c_poly_block_o/code1/cwp_config_srv.txt\" : -n 1 cwp_server -cn code2 -p 48103 48103 -c \"client_new_api_c_poly_block_o/code2/cwp_config_srv.txt\" &");
  }

  while (access(config, R_OK) != 0) {
    sleep(1);
  }
  sleep(10);

  FILE *meshFile;

  meshFile = fopen(CWP_MESH_DIR"mesh_poly_d1", "r"); // WARNING: adapt depending on where client is launched

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

  char *srcName = (char *) malloc(sizeof(char) * (strlen(__FILE__) + 1));
  strcpy(srcName, __FILE__);
  char *srcBaseName = NULL;
  srcBaseName = strrchr(srcName, '.');
  if (srcBaseName != NULL)
    *srcBaseName = '\0';
  srcBaseName = NULL;
  srcBaseName = strrchr(srcName, '/');
  if (srcBaseName != NULL)
    srcBaseName += 1;
  else
    srcBaseName = srcName;


  /* Initialization
   * -------------- */

  const char *codeName = NULL;

  if (rank == 0) {
    codeName = "cpoly";
  }

  if (rank == 1) {
    codeName = "code2";
  }

  CWP_Status_t is_coupled_rank = CWP_STATUS_ON;

  // Outputfile
  if (rank == 0) {
    FILE *f = fopen("output_file_code1.txt", "w");
    CWP_client_Output_file_set(f);
  }

  if (rank == 1) {
    FILE *f = fopen("output_file_code2.txt", "w");
    CWP_client_Output_file_set(f);
  }

  MPI_Comm intra_comm;
  int color = rank;
  MPI_Comm_split(comm, color, rank, &intra_comm);

  CWP_client_Init(intra_comm,
                  config,
                  codeName,
                  is_coupled_rank);

  // EXIT_SUCCESS ?
  int exit_check = 0;

  char cpl_id1[] = "cpl_code1_code2";

  if (rank == 0) {
    CWP_client_Cpl_create("cpoly", cpl_id1, "code2", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                          CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, 1,
                          CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (rank == 1) {
    CWP_client_Cpl_create("code2", cpl_id1, "cpoly", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                          CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, 1,
                          CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }

  /* Building of the local mesh */

  int dimension = 0;             // Dimension of the space
  int nVertex = 0;               // Number of points in the mesh
  int nFace = 0;                 // Number of face
  int nElements = 0;             // Number of cells
  int lFaceConnec = 0;
  int lCellConnec = 0;

  double *coords        = NULL;         // Coordinates of the points
  int    *faceVertexIdx = NULL;
  int    *faceVertex    = NULL;
  int    *cellFaceIdx   = NULL;
  int    *cellFace      = NULL;

  read_mesh_dim(meshFile, &dimension, &nVertex, &nFace, &nElements, &lFaceConnec, &lCellConnec);

  coords        = (double *) malloc(dimension * nVertex * sizeof(double));
  faceVertexIdx = (int    *) malloc((nFace + 1)         * sizeof(int   ));
  faceVertex    = (int    *) malloc(lFaceConnec         * sizeof(int   ));
  cellFaceIdx   = (int    *) malloc((nElements + 1)     * sizeof(int   ));
  cellFace      = (int    *) malloc(lCellConnec         * sizeof(int   ));

  read_mesh(meshFile,
            dimension,
            nVertex,
            nFace,
            nElements,
            lFaceConnec,
            lCellConnec,
            coords,
            faceVertexIdx,
            faceVertex,
            cellFaceIdx,
            cellFace);

  fclose(meshFile);

  if (rank == 0) {
    CWP_client_Visu_set("cpoly", cpl_id1, 1, CWP_VISU_FORMAT_ENSIGHT, "binary");

    CWP_client_Mesh_interf_vtx_set("cpoly", cpl_id1, 0, nVertex, coords, NULL);

    int block_id = CWP_client_Mesh_interf_block_add("cpoly", cpl_id1, CWP_BLOCK_CELL_POLY);

    CWP_client_Mesh_interf_c_poly_block_set("cpoly", cpl_id1, 0, block_id,
                                            nElements,
                                            nFace,
                                            faceVertexIdx,
                                            faceVertex,
                                            cellFaceIdx,
                                            cellFace,
                                            NULL); // global_num or try with gnum NULL

    CWP_g_num_t *cellGnum = NULL;
    int getNElements = -1;
    int getNFace = -1;
    int *getFaceVertexIdx = NULL;
    int *getFaceVertex    = NULL;
    int *getCellFaceIdx   = NULL;
    int *getCellFace      = NULL;
    CWP_client_Mesh_interf_c_poly_block_get("cpoly", cpl_id1, 0, block_id,
                                            &getNElements,
                                            &getNFace,
                                            &getFaceVertexIdx,
                                            &getFaceVertex,
                                            &getCellFaceIdx,
                                            &getCellFace,
                                            &cellGnum);

    // --> check
    if (!(getNElements == nElements)) {
      exit_check = 1;
    }

    if (!(getNFace == nFace)) {
      exit_check = 1;
    }

    for (int i = 0; i < nFace + 1; i++) {
      if (!(getFaceVertexIdx[i] == faceVertexIdx[i])) {
        exit_check = 1;
        break;
      }
    }

    for (int i = 0; i < faceVertexIdx[nFace]; i++) {
      if (!(getFaceVertex[i] == faceVertex[i])) {
        exit_check = 1;
        break;
      }
    }

    for (int i = 0; i < nElements + 1; i++) {
      if (!(getCellFaceIdx[i] == cellFaceIdx[i])) {
        exit_check = 1;
        break;
      }
    }

    for (int i = 0; i < cellFaceIdx[nElements]; i++) {
      if (!(getCellFace[i] == cellFace[i])) {
        exit_check = 1;
        break;
      }
    }

    if (cellGnum != NULL) free(cellGnum    );

    CWP_client_Mesh_interf_del("cpoly", cpl_id1);  }

  if (rank == 0) {
    CWP_client_Cpl_del("cpoly", cpl_id1);
  }

  if (rank == 1) {
    CWP_client_Cpl_del("code2", cpl_id1);
  }

  free(coords       );
  free(faceVertexIdx);
  free(faceVertex   );
  free(cellFaceIdx  );
  free(cellFace     );

  CWP_client_Finalize();

  MPI_Comm_free(&intra_comm);

  MPI_Finalize();

  free(srcName);

  return exit_check;
}
