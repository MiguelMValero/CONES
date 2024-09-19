/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <mpi.h>

#include "cwipi.h"
#include "grid_mesh.h"


/*----------------------------------------------------------------------
 *                                                                     
 * Dump status exchange                                                
 *                                                                     
 * parameters:
 *   status              <-- Exchange status           
 *---------------------------------------------------------------------*/

static void _dumpStatus(FILE* outputFile, cwipi_exchange_status_t status)
{
  switch(status) {
  case CWIPI_EXCHANGE_OK :
    fprintf(outputFile, "Exchange Ok\n");
    break;
  case CWIPI_EXCHANGE_BAD_RECEIVING :
    fprintf(outputFile, "Bad receiving\n");
    break;
  default :
    printf("Error : bad exchange status\n");
    exit(1);
  }
}

/*----------------------------------------------------------------------
 *                                                                     
 * Dump not located points                                             
 *                                                                     
 * parameters:
 *   coupling_id         <-- Coupling id               
 *   nNotLocatedPoints   <-- Number of not located points
 *---------------------------------------------------------------------*/

static void _dumpNotLocatedPoints(FILE* outputFile,
                                  const char *coupling_id,
                                  const int nNotLocatedPoints)
{
  if ( nNotLocatedPoints > 0) {
    fprintf(outputFile, "Not located points :\n");
    const int* notLocatedPoints = cwipi_get_not_located_points(coupling_id);
    for(int i = 0; i < nNotLocatedPoints; i++)
     fprintf(outputFile, "%i ", notLocatedPoints[i]);
    fprintf(outputFile, "\n");
  }
}

/*----------------------------------------------------------------------
 *                                                                     
 * Display usage                                             
 *                                                                     
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code)
{
  printf
    ("\n"
     "  Usage: \n\n"
     "  -n     <level>  Number of vertices in band width.\n\n"
     "  -rand  <level>  Random level ( > 0 and < 0.4) \n\n"
     "  -h             this message.\n\n");

  exit(exit_code);
}

/*----------------------------------------------------------------------
 *                                                                     
 * Read args from the command line                           
 *                                                                     
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth                         
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

static void
_read_args(int            argc,
           char         **argv,
           int          *nVertex,
           double       *randLevel)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *nVertex = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *randLevel = atof(argv[i]);
    }
    i++;
  }
}

double *shapef = NULL;

static void _userInterpolation2(const int entities_dim,
                               const int n_local_vertex,
                               const int n_local_element,
                               const int n_local_polhyedra,
                               const int n_distant_point,
                               const double local_coordinates[],
                               const int local_connectivity_index[],
                               const int local_connectivity[],
                               const int local_polyhedra_face_index[],
                               const int local_polyhedra_cell_to_face_connectivity[],
                               const int local_polyhedra_face_connectivity_index[],
                               const int local_polyhedra_face_connectivity[],
                               const double distant_points_coordinates[],
                               const int distant_points_location[],
                               const double distant_points_distance[],
                               const int distant_points_barycentric_coordinates_index[],
                               const double distant_points_barycentric_coordinates[],
                               const int stride,
                               const cwipi_solver_type_t  solver_type,
                               const void *local_field,
                               void *distant_field)
{

  // Compute shapef

  int compute_shape_f = 0;

  if (solver_type == CWIPI_SOLVER_CELL_VERTEX) {
    for (int i = 0; i < n_distant_point; i++) {

      int ielt = distant_points_location[i] - 1;
      int ivertex[4];

      ivertex[0] = local_connectivity[local_connectivity_index[ielt]  ] - 1;
      ivertex[1] = local_connectivity[local_connectivity_index[ielt]+1] - 1;
      ivertex[2] = local_connectivity[local_connectivity_index[ielt]+2] - 1;
      ivertex[3] = local_connectivity[local_connectivity_index[ielt]+3] - 1;

      if (shapef == NULL) {
        shapef = (double *) malloc(4 * n_local_element * sizeof(double));
        compute_shape_f = 1;
      }

      double *shapef_elt = shapef + 4 * ielt;

      //
      // Compute shape function
      //

      if (compute_shape_f == 1) {
        double deriv[4][2];
        double uv[2];
        double a[2][2];
        double b[2];
        double det_a;
        double x[2];
        double inv_a[2][2];

        for (int k = 0; k < 2; k++)
          uv[k] = 0.5;


        const int it_max = 100;
        for (int it = 0; it < it_max; it++) {

          shapef_elt[0] = (1 - uv[0]) * (1 - uv[1]);
          shapef_elt[1] = uv[0] * (1 - uv[1]);
          shapef_elt[2] = uv[0] * uv[1];
          shapef_elt[3] = (1 - uv[0]) * uv[1];
        
          deriv[0][0] = - (1 - uv[1]);
          deriv[0][1] = - (1 - uv[0]);
          deriv[1][0] =   (1 - uv[1]);
          deriv[1][1] = - uv[0];
          deriv[2][0] =   uv[1];
          deriv[2][1] =   uv[0];
          deriv[3][0] = - uv[1];
          deriv[3][1] =   (1 - uv[0]);

          for (int k = 0; k < 2; k++) {
            for (int l = 0; l < 2; l++)
              a[k][l] = 0.0;
          }
        
          b[0] = - distant_points_coordinates[3 * i    ];
          b[1] = - distant_points_coordinates[3 * i + 1];
        
          for (int k = 0; k < 4; k++) {
            
            b[0] += (shapef_elt[k] * local_coordinates[3 * ivertex[k]   ]); 
            b[1] += (shapef_elt[k] * local_coordinates[3 * ivertex[k] +1]); 
            
            for (int l = 0; l < 2; l++) { 
              a[0][l]  -=  (deriv[k][l] * local_coordinates[3 * ivertex[k]    ]);
              a[1][l]  -=  (deriv[k][l] * local_coordinates[3 * ivertex[k] + 1]);
            }
          }

          det_a = a[0][0] * a[1][1] - a[0][1] * a[1][0];  
          if (fabs(det_a) < 1e-12) {
            printf("matrice non inversible\n");
            exit(1);
          }

          double det_inv = 1./det_a;

          inv_a[0][0] =   det_inv * a[1][1];
          inv_a[0][1] = - det_inv * a[0][1];
          inv_a[1][0] =   det_inv * a[1][0];
          inv_a[1][1] =   det_inv * a[0][0];

          x[0] = inv_a[0][0] * b[0] + inv_a[0][1] * b[1]; 
          x[1] = inv_a[1][0] * b[0] + inv_a[1][1] * b[1]; 

          double dist = 0.0;

          for (int k = 0; k < 2; k++) {
            dist += x[k] * x[k];
            uv[k] += x[k];
          }

          if (dist <= 1e-5)
            break;
        }

        shapef_elt[0] = (1 - uv[0]) * (1 - uv[1]);
        shapef_elt[1] = uv[0] * (1 - uv[1]);
        shapef_elt[2] = uv[0] * uv[1];
        shapef_elt[3] = (1 - uv[0]) * uv[1];

      }

      //
      // Insterpolation
      //

      ((double *) distant_field)[i] = 0;
      
      for (int k = 0; k < 4; k++) {
        ((double *) distant_field)[i] += shapef_elt[k] * ((double *) local_field)[ivertex[k]];
      }
    }
  }
  else {
    printf("Error in _userInterpolation : bad solver_type\n");
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------
 *                                                                     
 * Main : surface coupling test : P1P1 
 *
 *---------------------------------------------------------------------*/
 
int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{

  FILE *outputFile;

  MPI_Init(&argc, &argv);

  int rank;
  int commWorldSize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  char *srcName = (char *) malloc (sizeof(char) * (strlen(__FILE__) + 1));
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

  srand(rank + time(0));

  int n_partition = 0;
  const int two = 2;
  while(two * pow(n_partition, two) < commWorldSize) n_partition++;

  int n2 = (int) (two * pow(n_partition, two));

  if (n2 != commWorldSize) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  /* Read args from command line
   * --------------------------- */

/*   int nVertexSeg = 10; */
/*   double randLevel = 0.4; */

  int nVertexSeg = 10;
  double randLevel = 0.4;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  /* Initialization
   * -------------- */

  char *codeName;
  int codeId;
  char *codeCoupledName;

  if (rank < commWorldSize / 2) {
    codeName = "code1";
    codeId = 1;
    codeCoupledName = "code2";
  }
  else {
    codeName = "code2";
    codeId = 2;
    codeCoupledName = "code1";
  }

  char* fileName = (char *) malloc(sizeof(char) * 45);
  sprintf(fileName,"c_surf_coupling_user_interpolation2%4.4d.txt",rank);

  outputFile = fopen(fileName,"w");

  free(fileName);

  cwipi_set_output_listing(outputFile);

  MPI_Comm localComm;
  cwipi_init(MPI_COMM_WORLD,
             codeName ,
             &localComm);

  /* Output redirection
   * ------------------ */

  int currentRank;
  int localCommSize;

  MPI_Comm_rank(localComm, &currentRank);
  MPI_Comm_size(localComm, &localCommSize);

  fprintf(outputFile, "  Surface coupling test : P1P1 with polygon\n");
  fprintf(outputFile, "\n");

  fprintf(outputFile, "\nDump after initialization\n");
  fprintf(outputFile, "-------------------------\n");
  cwipi_dump_application_properties();

  if (rank == 0)
    printf("        Create coupling\n");
  
  cwipi_solver_type_t solver_type;
  
  solver_type = CWIPI_SOLVER_CELL_VERTEX;
  
  /* Coupling creation
   * ----------------- */

  cwipi_create_coupling("c_surf_cpl_usr_interpolation2",                                // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        2,                                         // Geometric entities dimension
                        0.1,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        1,                                         // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option
  
  /* Mesh definition
   * --------------- */

  if (rank == 0)
    printf("        Create mesh\n");

  int nVertex = 0;               // Number of vertex
  double *coords = NULL;         // Vertex coordinates
  int nElts = 0;                 // Number of elements
  int *eltsConnecPointer = NULL; // Connectivity index
  int *eltsConnec = NULL;        // Connectivity
  
  /* Domain bounds */

  const double xmin = -10;
  const double xmax =  10;
  const double ymin = -10;
  const double ymax =  10;

  nVertex = nVertexSeg * nVertexSeg;
  nElts = (nVertexSeg - 1) * (nVertexSeg - 1);

  coords = (double *) malloc(sizeof(double) * 3 * nVertex );
  eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec = (int *) malloc(sizeof(int) * 4 * nElts);
  
  grid_mesh(xmin, 
            xmax, 
            ymin, 
            ymax, 
            randLevel,
            nVertexSeg,
            n_partition, 
            coords, 
            eltsConnecPointer,
            eltsConnec,
            localComm); 


  fprintf(outputFile, "   Number of vertex   : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  cwipi_define_mesh("c_surf_cpl_usr_interpolation2",
                    nVertex,
                    nElts,
                    coords,
                    eltsConnecPointer,
                    eltsConnec);

  /* Fields exchange
   *     - Proc 0 : Send X coordinates
   *                Recv Y coordinates
   *     - Proc 1 : Send Y coordinates
   *                Recv X coordinates
   * --------------------------------- */

  if (rank == 0)
    printf("        Exchange Code1 <-> Code2\n");

  double *sendValues = NULL;
  double *recvValues = NULL;
  
  sendValues = (double *) malloc(sizeof(double) * nVertex);
  recvValues = (double *) malloc(sizeof(double) * nVertex);

  /* Define fields to send (X coordinate or Y coordinate) */

  for (int i = 0; i < nVertex; i++) {
    if (codeId == 1)
      sendValues[i] = coords[3 * i];
    else
      sendValues[i] = coords[3 * i + 1];
  }

  /* Exchange */

  int nNotLocatedPoints = 0;
  char *sendValuesName;
  char *recvValuesName;
  if (codeId == 1) {
    sendValuesName = "cooX";
    recvValuesName = "cooY";
  }
  else{
    sendValuesName = "cooY";
    recvValuesName = "cooX";
  }

  cwipi_set_interpolation_function("c_surf_cpl_usr_interpolation2", _userInterpolation2);

  cwipi_exchange_status_t status = cwipi_exchange("c_surf_cpl_usr_interpolation2",
                                                  "ech",
                                                  1,
                                                  1,     // n_step
                                                  0.1,   // physical_time
                                                  sendValuesName,
                                                  sendValues,
                                                  recvValuesName,
                                                  recvValues,
                                                  &nNotLocatedPoints); 

  _dumpStatus(outputFile, status);
  _dumpNotLocatedPoints(outputFile, "c_surf_cpl_usr_interpolation2", nNotLocatedPoints);

  /* Coupling deletion
   * ----------------- */

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("c_surf_cpl_usr_interpolation2");

  /* Check barycentric coordinates */

  if (rank == 0)
    printf("        Check results\n");    

  double err;
  if (codeId == 1)
    err = fabs(recvValues[0] - coords[3 * 0 + 1]);
  else
    err = fabs(recvValues[0] - coords[3 * 0    ]);
 
  for (int i = 1; i < nVertex; i++) {
    if (codeId == 1)
      err = ((fabs(recvValues[i] - coords[3 * i + 1])) < (err) ? (err) : 
             (fabs(recvValues[i] - coords[3 * i + 1])));
    else
      err = ((fabs(recvValues[i] - coords[3 * i    ])) < (err) ? (err) : 
             (fabs(recvValues[i] - coords[3 * i    ])));
  }

  // TODO: Optimiser la fonction utilisateur !! 

  if (err >= 1e-2) {
    if (rank == 0) {
      printf("        !!! Error = %12.5e\n", err);
      return EXIT_FAILURE;
    }
  }

  /* Free memory
   * ----------- */

  free(coords);
  free(eltsConnecPointer);
  free(eltsConnec);
  free(sendValues);
  free(recvValues);
  free(srcName);

  /* Finalize
   * -------- */

  cwipi_finalize();

  MPI_Finalize();

  fclose(outputFile);

  return EXIT_SUCCESS;
}
