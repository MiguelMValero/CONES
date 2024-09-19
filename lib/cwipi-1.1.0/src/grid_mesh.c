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

#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "grid_mesh.h"

#define ABS(a)     ((a) <  0  ? -(a) : (a))


/*----------------------------------------------------------------------
 *                                                                     
 * Return a random number in [-1, 1]                                    
 *                                                                     
 * parameters:
 *   coupling_id         <-- Coupling id               
 *   nNotLocatedPoints   <-- Number of not located points
 *---------------------------------------------------------------------*/

static double _random01(void)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / ABS(rsigna - rsignb);
  double resultat =   sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}


/*----------------------------------------------------------------------
 *                                                                     
 * Dump status exchange                                                
 *                                                                     
 * parameters:
 *   xmin                <-- grid xmin (global) 
 *   xmax                <-- grid xmax (global)
 *   ymin                <-- grid ymin (global)
 *   ymax                <-- grid ymax (global)
 *   randLevel           <-- random level
 *   nVertexSeg          <-- number of vertices in X and Y
 *   nPartitionSeg       <-- number of partitions in X and Y
 *   nVertex             --> number of vertices of the local partition
 *   coords              --> coordinates of vertices of the local partition
 *   nElts               --> number of elements of the local partition
 *   eltsConnecPointer   --> connectivity index of the local partition 
 *   eltsConnec          --> connectivity of the local partition
 *   localComm           <-- MPI comm of the global grid
 *---------------------------------------------------------------------*/

void  grid_mesh(double xmin, 
                double xmax, 
                double ymin, 
                double ymax, 
                double randLevel,
                int nVertexSeg,
                int nPartitionSeg, 
                double *coords, 
                int *eltsConnecPointer,
                int *eltsConnec,
                MPI_Comm localComm)
{

  int currentRank;
  int localCommSize;

  MPI_Comm_rank(localComm, &currentRank);
  MPI_Comm_size(localComm, &localCommSize);

  /* Bandwidth */

  double lX = (xmax - xmin) / nPartitionSeg; 
  double lY = (ymax - ymin) / nPartitionSeg;

  /* Compute local partition bounds with random level */

  int nBoundVerticesSeg = nPartitionSeg + 1;
  int nBoundVertices = nBoundVerticesSeg * nBoundVerticesSeg;
  double *boundRanks = (double *) malloc(3 * sizeof(double) * nBoundVertices);

  if (currentRank == 0) {
    
    for (int j = 0; j < nBoundVerticesSeg; j++) {
      for (int i = 0; i < nBoundVerticesSeg; i++) {
        boundRanks[3 * (j * nBoundVerticesSeg + i)] = xmin + i * lX;
        boundRanks[3 * (j * nBoundVerticesSeg + i) + 1] = ymin + j * lY;
        boundRanks[3 * (j * nBoundVerticesSeg + i) + 2] = 0.;
        if (j != 0 && j != (nBoundVerticesSeg - 1))
          boundRanks[3 * (j * nBoundVerticesSeg + i) + 1] += _random01() * randLevel * lY;  
        if (i != 0 && i != (nBoundVerticesSeg - 1))
          boundRanks[3 * (j * nBoundVerticesSeg + i)] += _random01() * randLevel * lX;  
      }
    }
  }

  MPI_Bcast(boundRanks, 3 * nBoundVertices, MPI_DOUBLE, 0, localComm);

  double *boundRank = (double *) malloc(3 * sizeof(double) * 4);

  int p1, p2, p3, p4;

  int ii = currentRank % nPartitionSeg;
  int jj = currentRank / nPartitionSeg;

  p1 = (nBoundVerticesSeg * jj)       + ii;
  p2 = (nBoundVerticesSeg * jj)       + ii + 1;
  p3 = (nBoundVerticesSeg * (jj + 1)) + ii + 1;
  p4 = (nBoundVerticesSeg * (jj + 1)) + ii;

  boundRank[0 * 3 + 0] = boundRanks[3 * p1    ]; 
  boundRank[0 * 3 + 1] = boundRanks[3 * p1 + 1];
  boundRank[0 * 3 + 2] = boundRanks[3 * p1 + 2];

  boundRank[1 * 3 + 0] = boundRanks[3 * p2    ];
  boundRank[1 * 3 + 1] = boundRanks[3 * p2 + 1];
  boundRank[1 * 3 + 2] = boundRanks[3 * p2 + 2];

  boundRank[2 * 3 + 0] = boundRanks[3 * p3    ];
  boundRank[2 * 3 + 1] = boundRanks[3 * p3 + 1];
  boundRank[2 * 3 + 2] = boundRanks[3 * p3 + 2];

  boundRank[3 * 3 + 0] = boundRanks[3 * p4    ];
  boundRank[3 * 3 + 1] = boundRanks[3 * p4 + 1];
  boundRank[3 * 3 + 2] = boundRanks[3 * p4 + 2];

  free(boundRanks);

  /* Number of vertices and elements in the partition */

  const int nElts   = (nVertexSeg - 1) * (nVertexSeg - 1);

  /* Define coordinates */

  double deltaU = 2.0/(nVertexSeg - 1);
  double deltaV = 2.0/(nVertexSeg - 1);
  double u = -1;
  double v = -1;
  for (int j = 0; j < nVertexSeg; j++) {
    for (int i = 0; i < nVertexSeg; i++) {
      double randU = u;
      double randV = v;
      if ((i != 0) && (j != 0) && (j != nVertexSeg - 1) && (i != nVertexSeg - 1)) {
        randU +=  _random01() * randLevel * deltaU;
        randV +=  _random01() * randLevel * deltaV;
      }
      
      coords[3 * (j * nVertexSeg + i) + 0] = 
        0.25 * ((1 - randU - randV + randU * randV) * boundRank[0 * 3 + 0] +
                (1 + randU - randV - randU * randV) * boundRank[1 * 3 + 0] +
                (1 + randU + randV + randU * randV) * boundRank[2 * 3 + 0] +
                (1 - randU + randV - randU * randV) * boundRank[3 * 3 + 0] );
      
      coords[3 * (j * nVertexSeg + i) + 1] = 
        0.25 * ((1 - randU - randV + randU * randV) * boundRank[0 * 3 + 1] +
                (1 + randU - randV - randU * randV) * boundRank[1 * 3 + 1] +
                (1 + randU + randV + randU * randV) * boundRank[2 * 3 + 1] +
                (1 - randU + randV - randU * randV) * boundRank[3 * 3 + 1] );

      coords[3 * (j * nVertexSeg + i) + 2] = 0.;

      u += deltaU;
    }
    v += deltaV;
    u = -1;
  }

  free(boundRank);

  /* Define connectivity */

  eltsConnecPointer[0] = 0;
  for (int i = 1; i < nElts + 1; i++) 
    eltsConnecPointer[i] = eltsConnecPointer[i-1] + 4;

  int k = 0;
  for (int j = 0; j < (nVertexSeg - 1); j++) {
    for (int i = 0; i < (nVertexSeg - 1); i++) {
      eltsConnec[4 * k]     =       j * nVertexSeg + i     + 1;
      eltsConnec[4 * k + 1] =       j * nVertexSeg + i + 1 + 1;
      eltsConnec[4 * k + 2] = (j + 1) * nVertexSeg + i + 1 + 1;
      eltsConnec[4 * k + 3] = (j + 1) * nVertexSeg + i     + 1;
      k++;
    }
  }
}

void PROCF(grid_mesh_f, GRID_MESH_F)
     (double *xmin, 
      double *xmax, 
      double *ymin, 
      double *ymax, 
      double *randLevel,
      int *nVertexSeg,
      int *nPartitionSeg, 
      double *coords, 
      int *eltsConnecPointer,
      int *eltsConnec,
      MPI_Fint *localComm)
{

  MPI_Comm localComm_c =  MPI_Comm_f2c(*localComm);

  grid_mesh(*xmin, *xmax, *ymin, *ymax, *randLevel,
            *nVertexSeg, *nPartitionSeg, coords,
            eltsConnecPointer, eltsConnec, localComm_c);

}

/* rotation d'un angle alpha en degre */
void mesh_rotate(double *coords, int npt, double alpha) {
  double degrad = acos(-1.0)/180.;
  double sina, cosa, x, y;
  
  alpha = alpha * degrad;
  sina = sin(alpha);
  cosa = cos(alpha);

  for (int i = 0; i < npt; i++) {
    x = coords[i*3];
    y = coords[i*3 +1];
    coords[i*3]   = cosa*x - sina*y;
    coords[i*3+1] = sina*x + cosa*y;
  }
}

void carre2rond(double xmin, double xmax, double ymin, double ymax, double *coords, int nVertex)  {
  double xc, yc , x , y;
  double hyp, l;
  // tolerance pour eviter division par ZERO
  const double toler = 0.00000001;

  // remarque : ne fonctionne que dans un plan XY
  // coord centre
  xc = (xmax + xmin)/2.;
  yc = (ymax + ymin)/2.;

  for (int i = 0; i < nVertex ; i++) {
    x = coords[i*3] - xc;
    y = coords[i*3+1] - yc;

    hyp = sqrt(x*x + y*y);
    if (hyp > toler) {
      if (fabs(x) > fabs(y)) {
    l = fabs(cos(acos(x / hyp)));
      } else {
    l = fabs(sin(asin(y / hyp)));
      }
      coords[i*3] = xc + x*l;
      coords[i*3+1] = yc + y*l;
    }
  }
}

