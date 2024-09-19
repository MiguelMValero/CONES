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
#include "cwipi_priv.h"
#include "grid_mesh.h"


/*----------------------------------------------------------------------
 *                                                                     
 * Dump status exchange                                                
 *                                                                     
 * parameters:
 *   status              <-- Exchange status           
 *---------------------------------------------------------------------*/

/* static void _dumpStatus(FILE* outputFile, cwipi_exchange_status_t status) */
/* { */
/*   switch(status) { */
/*   case CWIPI_EXCHANGE_OK : */
/*     fprintf(outputFile, "Exchange Ok\n"); */
/*     break; */
/*   case CWIPI_EXCHANGE_BAD_RECEIVING : */
/*     fprintf(outputFile, "Bad receiving\n"); */
/*     break; */
/*   default : */
/*     printf("Error : bad exchange status\n"); */
/*     exit(1); */
/*   } */
/* } */

/*----------------------------------------------------------------------
 *                                                                     
 *                                             
 *                                                                     
 * parameters:
 *   status              <-- Exchange status           
 *---------------------------------------------------------------------*/

static double _f(double x, double y, double z)
{
  CWIPI_UNUSED(y);
  return x*x + z*z - x*z + z - x + 2. + 3*z; 

}


/*----------------------------------------------------------------------
 *                                                                     
 * Display usage                                             
 *                                                                     
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

/* static void */
/* _usage(int exit_code) */
/* { */
/*   printf */
/*     ("\n" */
/*      "  Usage: \n\n" */
/*      "  -n     <level>  Number of vertices in band width.\n\n" */
/*      "  -rand  <level>  Random level ( > 0 and < 0.4) \n\n" */
/*      "  -visu           CWP_VISU_FORMAT_ENSIGHT outputs \n\n" */
/*      "  -a              Unlocking communication \n\n" */
/*      "  -stdout         Standard output \n\n" */
/*      "  -h             this message.\n\n"); */

/*   exit(exit_code); */
/* } */

/*----------------------------------------------------------------------
 *                                                                     
 * Read args from the command line                           
 *                                                                     
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth                         
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

/* static void */
/* _read_args(int            argc, */
/*            char         **argv, */
/*            int          *nVertex, */
/*            double       *randLevel, */
/* 	   int          *randFromClock, */
/* 	   int          *postFreq, */
/* 	   int          *t_com, */
/* 	   int          *tostdout) */

/* { */
/*   int i = 1; */

/*   /\* Parse and check command line *\/ */

/*   while (i < argc) { */

/*     if (strcmp(argv[i], "-h") == 0) */
/*       _usage(EXIT_SUCCESS); */

/*     else if (strcmp(argv[i], "-n") == 0) { */
/*       i++; */
/*       if (i >= argc) */
/*         _usage(EXIT_FAILURE); */
/*       else */
/*         *nVertex = atoi(argv[i]); */
/*     } */
/*     else if (strcmp(argv[i], "-rand") == 0) { */
/*       i++; */
/*       if (i >= argc) */
/*         _usage(EXIT_FAILURE); */
/*       else */
/*         *randLevel = atof(argv[i]); */
/*     } */
/*     else if (strcmp(argv[i], "-randFromClock") == 0) { */
/*       *randFromClock = 1; */
/*     } */
/*     else if (strcmp(argv[i], "-a") == 0) { */
/*       *t_com = 1; */
/*     } */
/*     else if (strcmp(argv[i], "-visu") == 0) { */
/*       *postFreq = 1; */
/*     } */
/*     else if (strcmp(argv[i], "-stdout") == 0) { */
/*       *tostdout = 1; */
/*     } */
/*     else */
/*       _usage(EXIT_FAILURE); */
/*     i++; */
/*   } */
/* } */


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

  if (two != commWorldSize) {
    if (rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  /* Initialization
   * -------------- */

  const char *codeName;
  int codeId;
  const char *codeCoupledName;

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

  char* fileName = (char *) malloc(sizeof(char) * 32);
  sprintf(fileName,"c_surf_location_triaP2_%4.4d.txt",rank);

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

  fprintf(outputFile, "  Surface coupling test : location in tria P2\n");
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

  const int postFreq = -1;
  
  cwipi_create_coupling("c_surf_cpl_location_triaP2",                                // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        2,                                         // Geometric entities dimension
                        1e-1,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        postFreq,                                         // Postprocessing frequency
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

  const double xmin = -0.1;
  const double xmax =  0.1;
  const double zmin = -0.1;
  const double zmax =  0.1;

  nVertex = 6;
  nElts = 1;

  coords = (double *) malloc(sizeof(double) * 3 * nVertex );
  eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1));
  eltsConnec = (int *) malloc(sizeof(int) * nVertex);

  eltsConnecPointer[0] = 0;
  eltsConnecPointer[1] = 6;
  
  eltsConnec[0] = 1;
  eltsConnec[1] = 2;
  eltsConnec[2] = 3;
  eltsConnec[3] = 4;
  eltsConnec[4] = 5;
  eltsConnec[5] = 6;
  
  coords[0] = xmin;
  coords[1] = 0.;
  coords[2] = zmin;

  coords[3] = xmax;
  coords[4] = 0.;
  coords[5] = 0.;

  coords[6] = 0;
  coords[7] = 0.;
  coords[8] = zmax;
  
  coords[9] = (xmin + xmax) / 2.;
  coords[10] = 0.;
  coords[11] = zmin/2.;

  coords[12] = xmax/2.;
  coords[13] = 0.;
  coords[14] = zmax/2.;

  coords[15] = xmin/2.;
  coords[16] = 0.;
  coords[17] = (zmax+zmin)/2.;

  fprintf(outputFile, "   Number of vertex   : %i\n", nVertex);
  fprintf(outputFile, "   Number of elements : %i\n", nElts);

  const int order = 2;
  
  cwipi_ho_define_mesh("c_surf_cpl_location_triaP2",
                       nVertex,
                       nElts,
                       order,
                       coords,
                       eltsConnecPointer,
                       eltsConnec);

  const int n_node = 6;
  
  int *ijk = malloc(sizeof(int)*2*n_node);

  ijk[0] = 0;
  ijk[1] = 0;

  ijk[2] = 2;
  ijk[3] = 0;

  ijk[4] = 0;
  ijk[5] = 2;
  
  ijk[6] = 1;
  ijk[7] = 0;

  ijk[8] = 1;
  ijk[9] = 1;

  ijk[10] = 0;
  ijk[11] = 1;
  
  cwipi_ho_ordering_from_IJK_set ("c_surf_cpl_location_triaP2",
                                  CWIPI_FACE_TRIAHO,
                                  n_node,
                                  ijk);
  
  int n_pts_to_locate = 11;

  double *pts_to_locate = (double *) malloc(sizeof(double) * 3 * n_pts_to_locate);

  pts_to_locate[0] = xmin;
  pts_to_locate[1] = 0.;
  pts_to_locate[2] = zmin;

  pts_to_locate[0] = xmax;
  pts_to_locate[1] = 0.;
  pts_to_locate[2] = 0.;

  pts_to_locate[3] = xmax;
  pts_to_locate[4] = 0.;
  pts_to_locate[5] = 0.;

  pts_to_locate[6] = 0;
  pts_to_locate[7] = 0.;
  pts_to_locate[8] = zmax;
  
  pts_to_locate[9] = (xmin + xmax) / 2.;
  pts_to_locate[10] = 0.;
  pts_to_locate[11] = zmin/2.;

  pts_to_locate[12] = xmax/2.;
  pts_to_locate[13] = 0.;
  pts_to_locate[14] = zmax/2.;

  pts_to_locate[15] = xmin/2.;
  pts_to_locate[16] = 0.;
  pts_to_locate[17] = (zmax+zmin)/2.;
  
  pts_to_locate[18] = 0.;
  pts_to_locate[19] = 0.;
  pts_to_locate[20] = 0.;

  pts_to_locate[21] = 0.;
  pts_to_locate[22] = 0.005;
  pts_to_locate[23] = 0.;

  pts_to_locate[24] = xmin - 0.005;
  pts_to_locate[25] = 0.;
  pts_to_locate[26] = zmin - 0.005;

  pts_to_locate[27] = (xmax + xmax/2.) / 2 + 0.005;
  pts_to_locate[28] = 0.;
  pts_to_locate[29] = zmax/4. + 0.005;

  pts_to_locate[30] = xmax/4.;
  pts_to_locate[31] = 0.;
  pts_to_locate[32] = zmax/4.;

  
  //  n_pts_to_locate = 1;
  cwipi_set_points_to_locate ("c_surf_cpl_location_triaP2",
                              n_pts_to_locate,
                              pts_to_locate);
  
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
  recvValues = (double *) malloc(sizeof(double) * n_pts_to_locate);

  /* Define fields to send (X coordinate or Y coordinate) */

  for (int i = 0; i < nVertex; i++) {
    sendValues[i] = _f(coords[3 * i], coords[3 * i+1], coords[3 * i+2]);
  }

  /* Exchange */

  int nNotLocatedPoints = 0;
  const char *sendValuesName;
  const char *recvValuesName;

  sendValuesName = "_fs";
  recvValuesName = "_fr";

  cwipi_locate("c_surf_cpl_location_triaP2");

  nNotLocatedPoints = cwipi_get_n_not_located_points("c_surf_cpl_location_triaP2");
  if (nNotLocatedPoints > 0) {
    printf("--- Error --- : %d not located points found\n", nNotLocatedPoints);
    exit(1);
  }

  int sRequest, rRequest;
  int tag = 1;
  
  cwipi_irecv("c_surf_cpl_location_triaP2",
              "ech",
              tag,
              1,
              1,
              0.1,
              recvValuesName,
              recvValues,
              &rRequest);
  

  cwipi_issend("c_surf_cpl_location_triaP2",
               "ech",
               tag,
               1,
               1,
               0.1,
               sendValuesName,
               sendValues,
               &sRequest);
  
  cwipi_wait_irecv("c_surf_cpl_location_triaP2", rRequest);
  cwipi_wait_issend("c_surf_cpl_location_triaP2", sRequest);


  /* Coupling deletion
   * ----------------- */

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("c_surf_cpl_location_triaP2");


  /* Check barycentric coordinates */

  if (rank == 0)
    printf("        Check results\n");    

  double *res = (double *) malloc(sizeof(double) *  n_pts_to_locate);

  for (int i = 0; i < nVertex; i++) {
    res[i] = sendValues[i];
  }
  
  res[nVertex    ] = _f(0.           , 0.             , 0.             );
  res[nVertex + 1] = _f(0.           , 0.             , 0.             );
  res[nVertex + 2] = res[0];  
  res[nVertex + 3] = _f( (xmax + xmax/2.) / 2, 0., zmax/4.);
  res[nVertex + 4] = _f( xmax/4., 0., zmax/4.);

  double err;

  for (int i = 0; i < n_pts_to_locate; i++) {
    printf (" %12.5e ", res[i]);

  }
  for (int i = 0; i < n_pts_to_locate; i++) {
    printf (" %12.5e ", recvValues[i]);

  }

  for (int i = 0; i < n_pts_to_locate; i++) {
    err = fabs(recvValues[i] - res[i]);
    //    if (err > 1e-6) {
    printf ("[%d] err %d : %12.5e %12.5e %12.5e\n", codeId, i, err, recvValues[i], res[i]);
      // }
  }

  double err_max;
  MPI_Allreduce(&err, &err_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  if (err_max >= 1e-6) {
    if (rank == 0) {
      printf("        !!! Error = %12.5e\n", err_max);
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  /* Free
   * ---- */
  
  free(pts_to_locate);
  free(res);
  free(ijk);
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
