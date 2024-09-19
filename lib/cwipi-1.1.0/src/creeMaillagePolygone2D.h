#ifndef __CREE_MAIL_POLY_H__
#define __CREE_MAIL_POLY_H__
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

/* #include <fvm_defs.h> */
/* #include <fvm_nodal.h> */
#include <mpi.h>

#include "cwipi_config.h"
#include "grid_mesh.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

//#define PROCF(x, y) x
#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

void creeMaillagePolygone2D(int order,
                            MPI_Comm localComm,
                            double xmin,
                            double xmax,
                            double ymin,
                            double ymax,
                            int initRandom,
                            int nx,
                            int ny,
                            int *nVertex,
                            double **meshCoords,
                            int *nElts,
                            int **eltsConnecPointer,
                            int **eltsConnec
                            );

void PROCF(creemaillagepolygone2d_f, CREEMAILLAGEPOLYGONE2D_F)(int *order,
							       MPI_Fint* localFComm,
                                                               double *xmin,
                                                               double *xmax,
                                                               double *ymin,
                                                               double *ymax,
                                                               int *initRandom,
                                                               int *nx,
                                                               int *ny,
                                                               int *nVertex,
                                                               double *meshCoords_f,
                                                               int *nElts,
                                                               int *lEltsConnecPointer_f,
                                                               int *eltsConnecPointer_f,
                                                               int *eltsConnec_f
                                                               );

#ifdef __cplusplus
}
#endif
#endif
