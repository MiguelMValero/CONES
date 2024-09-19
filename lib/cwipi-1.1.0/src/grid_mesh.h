#ifndef __GRID_MESH_H__
#define __GRID_MESH_H__
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

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#include <mpi.h>

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

void grid_mesh(double xmin, 
               double xmax, 
               double ymin, 
               double ymax, 
               double randLevel,
               int nVertexSeg,
               int nPartitionSeg, 
               double *coords, 
               int *eltsConnecPointer,
               int *eltsConnec,
               MPI_Comm localComm);

void  mesh_rotate(double *coords, int npt, double alpha);
void carre2rond(double xmin, double xmax, double ymin, double ymax, double *coords, int nVertex);

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
      MPI_Fint *localComm);

#endif
