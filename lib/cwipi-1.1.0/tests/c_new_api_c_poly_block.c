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
#include <assert.h>

#include "cwp.h"
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
  assert(f!=NULL);
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


/*----------------------------------------------------------------------
 *
 * Main : linear coupling test
 *
 *---------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  FILE *meshFile;


  meshFile = fopen(CWP_MESH_DIR"mesh_poly_d1", "r");

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

  if (rank == 0)
    printf("\nSTART: %s\n", srcBaseName);


  /* Initialization
   * -------------- */

  int n_code_name = 0;
  const char **codeNames = NULL;
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  if (rank == 0) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "cpoly";
  }

  if (rank == 1) {
    n_code_name = 1;
    codeNames = malloc(sizeof(char *) * n_code_name);
    codeNames[0] = "code2";
  }

  //CWP_Output_file_set (outputFile);

  MPI_Comm *localComm = malloc(sizeof(MPI_Comm) * n_code_name);

  CWP_Init(MPI_COMM_WORLD,
           n_code_name,
           codeNames,
           is_active_rank,
           localComm);


  int currentRank;
  int localCommSize;

  for (int i = 0 ; i < n_code_name ; i++) {
    MPI_Comm_rank(localComm[i], &currentRank);
    MPI_Comm_size(localComm[i], &localCommSize);
    printf("Size of localComm[%i]=%i et rang du proc=%i.\n", i, localCommSize, currentRank);
  }

  char cpl_id1[] = "cpl_code1_code2";

  printf("Coupling creation\n");
  if (rank == 0) {
    CWP_Cpl_create("cpoly", cpl_id1, "code2", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, 1,
                   CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (rank == 1) {
    CWP_Cpl_create("code2", cpl_id1, "cpoly", CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITHOUT_PART,
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, 1,
                   CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
  }
  printf("Coupling created %i\n", currentRank);

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

  if (rank == 0)
    printf("        Read mesh\n");

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

    printf("vtx_set\n");
    CWP_Mesh_interf_vtx_set("cpoly", cpl_id1, 0, nVertex, coords, NULL);

    printf("3D Cell Polyhedra Block Add\n");
    int block_id = CWP_Mesh_interf_block_add("cpoly", cpl_id1, CWP_BLOCK_CELL_POLY);

    printf("3D Cell Polyhedra Block Set\n");
    CWP_Mesh_interf_c_poly_block_set("cpoly", cpl_id1, 0, block_id,
                                     nElements,
                                     nFace,
                                     faceVertexIdx,
                                     faceVertex,
                                     cellFaceIdx,
                                     cellFace,
                                     NULL);


    printf("Interface Mesh deletion\n");
    CWP_Mesh_interf_del("cpoly", cpl_id1);
    printf("Interface Mesh deleted\n");
    fflush(stdout);
  }

  if (rank == 0) {
    CWP_Cpl_del("cpoly", cpl_id1);
  }

  if (rank == 1) {
    CWP_Cpl_del("code2", cpl_id1);
  }

  free(coords       );
  free(faceVertexIdx);
  free(faceVertex   );
  free(cellFaceIdx  );
  free(cellFace     );

  CWP_Finalize();
  MPI_Finalize();

  free(srcName);
  free(localComm);
  free(codeNames);

  return 0;
}
