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

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>
#include <cfloat>

#include <bftc_printf.h>
#include "fvmc_point_location.h"
#include "fvmc_ho_basis.h"

#include "locationToLocalMesh.hxx"
#include "locationToDistantMesh.hxx"
#include "vectorUtilities.hxx"
#include "geomUtilities.hxx"
#include "applicationProperties.hxx"

#define _DOT_PRODUCT(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

#define _MODULE(vect) \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])

#define _CROSS_PRODUCT(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define _PI 3.1415926535897931


namespace cwipi
{


/* Compute the polygon normal from an array of points. This version assumes */
/* that the polygon is convex, and looks for the first valid normal. */

void LocationToLocalMesh::_computeBary (int numPts, double *pts, double bary[3])
{
  bary[0] = 0.;
  bary[1] = 0.;
  bary[2] = 0.;

  for (int i = 0; i < 3; i++) {
    for (int ipt = 0; ipt < numPts; ipt++) {
      bary[i] += pts[3*ipt+i];
    }
    bary[i] /= numPts;
  }

}


void LocationToLocalMesh::_project_point2(double x[3], double pt_plan[3],
                            double normal[3], double xproj[3])
{

  double cst   = _DOT_PRODUCT(normal, pt_plan);
  double cst1  = _DOT_PRODUCT(normal, x);
  double norm2 = _DOT_PRODUCT(normal, normal);

  double t = - (cst1 - cst)/ norm2;

  xproj[0] = x[0] + t * normal[0];
  xproj[1] = x[1] + t * normal[1];
  xproj[2] = x[2] + t * normal[2];
}

/* ---------------------------------------------------------------------------- */
/* Compute the polygon normal from an array of points. This version assumes */
/* that the polygon is convex, and looks for the first valid normal. */

void LocationToLocalMesh::_computeNormal (int numPts, double *pts, double n[3])
{
  double length = 0.;
  double bary[3]= {0., 0., 0.};

  n[0] = 0.;
  n[1] = 0.;
  n[2] = 0.;

  _computeBary (numPts, pts, bary);

  for (int ipt = 0; ipt < numPts; ipt++) {

    double *pt1 = pts + 3 * ipt;
    double *pt2 = pts + 3 * ((ipt+1)%numPts);
    double vect1[3];
    double vect2[3];

    for (int i = 0; i < 3; i++) {
      vect1[i] = pt1[i] - bary[i];
      vect2[i] = pt2[i] - bary[i];
    }

    n[0] += vect1[1] * vect2[2] - vect1[2] * vect2[1];
    n[1] += vect1[2] * vect2[0] - vect1[0] * vect2[2];
    n[2] += vect1[0] * vect2[1] - vect1[1] * vect2[0];

  } //over all points

  length = sqrt (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (length != 0.0) {
    n[0] /= length;
    n[1] /= length;
    n[2] /= length;

  }
  return;
}


LocationToLocalMesh::LocationToLocalMesh(
                                         const cwipi_solver_type_t  &solverType,
                                         const int  optBboxStep,
                                         const double &tolerance,
                                         const MPI_Comm& couplingComm,
                                         const int &coupledApplicationNRankCouplingComm,
                                         const int &coupledApplicationBeginningRankCouplingComm,
                                         const bool isCoupledRank,
                                         const int entitiesDim,
                                         const ApplicationProperties& localApplicationProperties,
                                         LocationToDistantMesh &locationToDistantMesh)

: _supportMesh(NULL), _solverType(solverType), _optBboxStep(optBboxStep), _tolerance(tolerance),
  _couplingComm(couplingComm),
  _coupledApplicationNRankCouplingComm(coupledApplicationNRankCouplingComm),
  _coupledApplicationBeginningRankCouplingComm(coupledApplicationBeginningRankCouplingComm),
  _isCoupledRank(isCoupledRank), _entitiesDim(entitiesDim), _localApplicationProperties(localApplicationProperties),
  _locationToDistantMesh(locationToDistantMesh),
  _locatedPointsDistribution(NULL), _distantDistribution(NULL)
{
  _fvmLocator = NULL;
  _barycentricCoordinatesIndex = NULL;
  _barycentricCoordinates = NULL;
  _uvw = NULL;
  _nDistantPoint = 0;
  _location = NULL;
  _toLocate = true;
  _maxElementContainingNVertex = -1;
  _nVertex = NULL;
  _supportMesh = NULL;

}

LocationToLocalMesh::~LocationToLocalMesh()
{

  if (_barycentricCoordinatesIndex != NULL) {
    delete _barycentricCoordinatesIndex;
    delete _barycentricCoordinates;
    _barycentricCoordinatesIndex = NULL;
    _barycentricCoordinates = NULL;
  }

  if (_uvw != NULL) {
    delete _uvw;
    _uvw = NULL;
  }
  if (_fvmLocator != NULL)
    fvmc_locator_destroy(_fvmLocator);

  if (_nVertex != NULL)
    delete _nVertex;

}

size_t LocationToLocalMesh::locationSize()
{
  size_t il_size = 0;

  if ( _isCoupledRank) {
    il_size += sizeof(int);
    il_size += fvmc_locator_size(_fvmLocator);
    const int nDistantPoint = fvmc_locator_get_n_dist_points(_fvmLocator);

    il_size += sizeof(int);
    if (_barycentricCoordinatesIndex != NULL) {
      il_size += (nDistantPoint + 1)*sizeof(int);
    }

    il_size += sizeof(int);
    if (_barycentricCoordinates != NULL) {
      std::vector <double> &  _refBarycentricCoordinates = *_barycentricCoordinates;
      il_size += _refBarycentricCoordinates.size()*sizeof(double);
    }

    il_size += sizeof(int);
    if (_uvw != NULL) {
      std::vector <double> &  _refUvw = *_uvw;
      il_size += _refUvw.size()*sizeof(double);
    }

    il_size += sizeof(int);
    if (_nVertex != NULL) {
      il_size +=  _nDistantPoint*sizeof(int);
    }
  }
  return il_size;
}




void LocationToLocalMesh::packLocation(unsigned char *buff)
{
  int s;
  void *p;
  size_t s_pack;

  p = (void *)buff;

  if ( _isCoupledRank) {

    if (_fvmLocator != NULL) {
      s = 1;
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
      p = fvmc_locator_pack(p, _fvmLocator);
    }
    else {
      s = 0;
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
    }

    const int nDistantPoint      = fvmc_locator_get_n_dist_points(_fvmLocator);
    if (_barycentricCoordinatesIndex != NULL) {
      s = 1;
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
      std::vector <int> &  _refBarycentricCoordinatesIndex = *_barycentricCoordinatesIndex;
      s_pack = (nDistantPoint + 1)*sizeof(int);
      memcpy(p,(void *)&_refBarycentricCoordinatesIndex[0], s_pack);
      p = (void *) ((char *) p + s_pack);
    }
    else {
      s = 0;
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
    }

    if (_barycentricCoordinates != NULL) {
      std::vector <double> &  _refBarycentricCoordinates = *_barycentricCoordinates;
      // calcul de la taille de _barycentricCoordinates
      s = _refBarycentricCoordinates.size();
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
      s_pack = s * sizeof(double);
      memcpy(p,(void *)&_refBarycentricCoordinates[0], s_pack);
      p = (void *) ((char *) p + s_pack);
    }
    else {
      s = 0;
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
    }

    if (_uvw != NULL) {
      std::vector <double> &  _refUvw = *_uvw;
      // calcul de la taille de _uvw
      s = _refUvw.size();
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
      s_pack = s * sizeof(double);
      memcpy(p,(void *)&_refUvw[0], s_pack);
      p = (void *) ((char *) p + s_pack);
    }
    else {
      s = 0;
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
    }

    if (_nVertex != NULL) {
      s = 1;
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
      std::vector <int> & _nVertexRef = *_nVertex;
      s_pack = _nDistantPoint*sizeof(int);
      memcpy(p,(void *)&_nVertexRef[0], s_pack);
      p = (void *) ((char *) p + s_pack);
    }
    else {
      s = 0;
      s_pack = sizeof(int);
      memcpy(p,(void *)&s, s_pack);
      p = (void *) ((char *) p + s_pack);
    }
  }
}


void LocationToLocalMesh::unpackLocation(unsigned char *buff)
{
  int s;
  size_t cur_pos;
  cur_pos = 0;

  if ( _isCoupledRank) {
    // read the locator
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
    if (s == 1) {
      if (_fvmLocator != NULL) {
        fvmc_locator_destroy(_fvmLocator);
      }

      _fvmLocator = fvmc_locator_create(_optBboxStep,
                                        _tolerance,
                                        _couplingComm,
                                        _coupledApplicationNRankCouplingComm,
                                        _coupledApplicationBeginningRankCouplingComm);
      cur_pos += fvmc_locator_unpack(&buff[cur_pos],_fvmLocator);

      // mise à jour de l'objet locationToLocalMesh

      const int* locationList = fvmc_locator_get_dist_locations(_fvmLocator);
      const float* distanceList = fvmc_locator_get_dist_distances(_fvmLocator);
      _location = const_cast<int *> (locationList);
      _distance = const_cast<float *> (distanceList);

      _distantDistribution =  const_cast<int *> (fvmc_locator_get_dist_distrib(_fvmLocator));
      _locatedPointsDistribution =  const_cast<int *> (fvmc_locator_get_loc_distrib(_fvmLocator));

      _nDistantPoint = fvmc_locator_get_n_dist_points(_fvmLocator);

      // mise à jour de l'objet locationToDistantMesh

      const int* exteriorList = fvmc_locator_get_exterior_list(_fvmLocator);
      const int* interiorList = fvmc_locator_get_interior_list(_fvmLocator);
      _locationToDistantMesh._unlocatedPoint = const_cast<int *> (exteriorList);
      _locationToDistantMesh._locatedPoint = const_cast<int *> (interiorList);


      // printf("LocationToLocalMesh::load ATTENTION DUMP \n");
      // fvmc_locator_dump(_fvmLocator);

    }

    // load the barycentric coordinates
    const int nDistantPoint      = fvmc_locator_get_n_dist_points(_fvmLocator);

    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
    if (s == 1) {
      if (_barycentricCoordinatesIndex != NULL) delete _barycentricCoordinatesIndex;
      _barycentricCoordinatesIndex = new std::vector <int> (nDistantPoint + 1);
      std::vector <int> &  _refBarycentricCoordinatesIndex = *_barycentricCoordinatesIndex;
      cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_refBarycentricCoordinatesIndex[0], (nDistantPoint + 1)*sizeof(int));
    }

    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
    if (s != 0) {
      if (_barycentricCoordinates != NULL) delete _barycentricCoordinates;
      _barycentricCoordinates = new std::vector <double> (s);
      std::vector <double> &  _refBarycentricCoordinates = *_barycentricCoordinates;
      cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_refBarycentricCoordinates[0], s*sizeof(double));
    }

    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
    if (s != 0) {
      if (_uvw != NULL) delete _uvw;
      _uvw = new std::vector <double> (s);
      std::vector <double> &  _refUvw = *_uvw;
      cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_refUvw[0], s*sizeof(double));
    }

    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
    if (s == 1) {
      if(_nVertex != NULL) delete _nVertex;
      _nVertex = new std::vector <int> (_nDistantPoint, 0);
      std::vector <int> & _nVertexRef = *_nVertex;
      cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_nVertexRef[0], _nDistantPoint*sizeof(int));
    }
    _toLocate = false;
    _locationToDistantMesh._toLocate = false;
  }
}


void LocationToLocalMesh::locate()
{
  //
  // Create locator

  // TODO: Exchange MPI of _toLocate param between master rank !!!
  // TODO: Attention la fonction clear ne devrait-elle pas etre dans le if suivant !

  //  _locationToDistantMesh.clear();

  if (_toLocate || _locationToDistantMesh.getToLocateStatus()) {


    _locationToDistantMesh.clear();

    if (_isCoupledRank) {


      if (_supportMesh == NULL)
        bftc_error(__FILE__, __LINE__, 0, "undefined support mesh\n");

      if (_fvmLocator == NULL)
        _fvmLocator = fvmc_locator_create(_optBboxStep,
                                          _tolerance,
                                          _couplingComm,
                                          _coupledApplicationNRankCouplingComm,
                                          _coupledApplicationBeginningRankCouplingComm);

      // TODO: Revoir les coordonnees des points a localiser
      // (cas centres sommets + centres faces + autres points)
      // TODO: Ajouter un locator pour les sommets pour les centres faces,...

      double* coords = NULL;
      if (_locationToDistantMesh._coordsPointsToLocate != NULL)
        coords = _locationToDistantMesh._coordsPointsToLocate;

      else if(_solverType == CWIPI_SOLVER_CELL_CENTER) {
        _locationToDistantMesh._nPointsToLocate = _supportMesh->getNElts();
        coords = const_cast <double*> (&(_supportMesh->getCellCenterCoords()[0]));
      }

      else if(_solverType == CWIPI_SOLVER_CELL_VERTEX) {
        _locationToDistantMesh._nPointsToLocate = _supportMesh->getNVertex();
        coords = const_cast <double*> (_supportMesh->getVertexCoords());
      }

      fvmc_locator_set_nodal(_fvmLocator,
                             &_supportMesh->getFvmNodal(),
                             0,
                             3,
                             _locationToDistantMesh._nPointsToLocate,
                             NULL,
                             coords);

      _toLocate = false;
      const int nLocatedPoint = fvmc_locator_get_n_interior(_fvmLocator);
      const int nNotLocatedPoint = _locationToDistantMesh._nPointsToLocate - nLocatedPoint;
      const int* exteriorList = fvmc_locator_get_exterior_list(_fvmLocator);
      const int* interiorList = fvmc_locator_get_interior_list(_fvmLocator);
      const int* locationList = fvmc_locator_get_dist_locations(_fvmLocator);
      const float* distanceList = fvmc_locator_get_dist_distances(_fvmLocator);
      const int nExterior = fvmc_locator_get_n_exterior(_fvmLocator);
      assert(nNotLocatedPoint == nExterior);
      _locationToDistantMesh._unlocatedPoint = const_cast<int *> (exteriorList);
      _locationToDistantMesh._locatedPoint = const_cast<int *> (interiorList);
      _locationToDistantMesh._nUnlocatedPoint = nNotLocatedPoint;
      _locationToDistantMesh._nLocatedPoint = nLocatedPoint;
      _locationToDistantMesh._toLocate = false;

      _location = const_cast<int *> (locationList);
      _distance = const_cast<float *> (distanceList);

      _distantDistribution =  const_cast<int *> (fvmc_locator_get_dist_distrib(_fvmLocator));
      _locatedPointsDistribution =  const_cast<int *> (fvmc_locator_get_loc_distrib(_fvmLocator));

      _nDistantPoint = fvmc_locator_get_n_dist_points(_fvmLocator);

      if (_barycentricCoordinatesIndex != NULL) {
        delete _barycentricCoordinatesIndex;
        delete _barycentricCoordinates;
        _barycentricCoordinatesIndex = NULL;
        _barycentricCoordinates = NULL;
      }

      if (_uvw != NULL) {
        delete _uvw;
        _uvw = NULL;
      }

      //
      // TODO: Prevoir une fabrique pour supprimer les tests if sur _entitiesDim
      //       Le calcul des coordonnees barycentriques se fera dans cette fabrique

      if (_barycentricCoordinatesIndex == NULL) {
        if (_entitiesDim == 1) {
          const int nDistantPoint      = fvmc_locator_get_n_dist_points(_fvmLocator);
          const int *distantLocation   = fvmc_locator_get_dist_locations(_fvmLocator);
          const double *distantCoords   = fvmc_locator_get_dist_coords(_fvmLocator);

          const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
          const int *eltsConnec = _supportMesh->getEltConnectivity();
          const double *localCoords    = _supportMesh->getVertexCoords();

          const int order = _supportMesh->getOrder();
          //          const fvmc_coord_t *proj_coords = fvmc_locator_get_dist_projected_coords(_fvmLocator);

          if (order == -1) {

            if ( nDistantPoint > 0 ) {
              _barycentricCoordinatesIndex = new std::vector <int> (nDistantPoint + 1);
              _barycentricCoordinates = new std::vector <double> (2 * nDistantPoint);
              std::vector <int> &  _refBarycentricCoordinatesIndex = *_barycentricCoordinatesIndex;
              std::vector <double> &  _refBarycentricCoordinates = *_barycentricCoordinates;

              _refBarycentricCoordinatesIndex[0] = 0;

              for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
                int iel = distantLocation[ipoint] - 1;
                _refBarycentricCoordinatesIndex[ipoint+1] = _refBarycentricCoordinatesIndex[ipoint] + 2;
                int index = eltsConnecPointer[iel];
                int nVertex = eltsConnecPointer[iel+1] - eltsConnecPointer[iel];
                assert(nVertex == 2);
                int pt1 = eltsConnec[index] - 1;
                int pt2 = eltsConnec[index+1] - 1;
                double coef1 = sqrt((localCoords[3*pt1]-distantCoords[3*ipoint])
                                    *(localCoords[3*pt1]-distantCoords[3*ipoint])+
                                    (localCoords[3*pt1+1]-distantCoords[3*ipoint+1])
                                    *(localCoords[3*pt1+1]-distantCoords[3*ipoint+1])+
                                    (localCoords[3*pt1+2]-distantCoords[3*ipoint+2])
                                    *(localCoords[3*pt1+2]-distantCoords[3*ipoint+2]));
                double coef2 = sqrt((localCoords[3*pt2]-distantCoords[3*ipoint])
                                    *(localCoords[3*pt2]-distantCoords[3*ipoint])+
                                    (localCoords[3*pt2+1]-distantCoords[3*ipoint+1])
                                    *(localCoords[3*pt2+1]-distantCoords[3*ipoint+1])+
                                    (localCoords[3*pt2+2]-distantCoords[3*ipoint+2])
                                    *(localCoords[3*pt2+2]-distantCoords[3*ipoint+2]));
                _refBarycentricCoordinates[_refBarycentricCoordinatesIndex[ipoint]] =
                  coef2/(coef1+coef2);
                _refBarycentricCoordinates[_refBarycentricCoordinatesIndex[ipoint]+1] =
                  coef1/(coef1+coef2);
              }
            }
          }

          else {
            _barycentricCoordinatesIndex = new std::vector <int> (nDistantPoint + 1);
            _barycentricCoordinates = new std::vector <double> ((order+ 1) * nDistantPoint);
            std::vector <int> &  _refBarycentricCoordinatesIndex = *_barycentricCoordinatesIndex;
            std::vector <double> &  _refBarycentricCoordinates = *_barycentricCoordinates;

            int max_entity_dim = fvmc_nodal_get_max_entity_dim (&(_supportMesh->getFvmNodal()));

            _refBarycentricCoordinatesIndex[0] = 0;
            for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
              _refBarycentricCoordinatesIndex[ipoint + 1] =
                _refBarycentricCoordinatesIndex[ipoint] + (order + 1);
            }

            const int *meshConnectivityIndex = _supportMesh->getEltConnectivityIndex();

            const double *uvw = fvmc_locator_get_dist_uvw (_fvmLocator);

            for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {

              int ielt = locationList[ipoint] - 1;

              int n_node = meshConnectivityIndex[ielt+1] - meshConnectivityIndex[ielt];

              double *_weights = &(_refBarycentricCoordinates[0]) +
                                 _refBarycentricCoordinatesIndex[ipoint];

              const double *_uvw_point_fvm = uvw + ipoint * max_entity_dim;

              FVMC_ho_basis (FVMC_EDGE, order, n_node, 1, _uvw_point_fvm, _weights);

            }
          }
        }

        else if (_entitiesDim == 2) {
          compute2DMeanValues();
        }

        else if (_entitiesDim == 3) {
          compute3DMeanValues();
        }
      }

      //
      // Exchange info status between root ranks

      int currentRank;
      MPI_Comm_rank(_couplingComm, &currentRank);
      cwipi_located_point_info_t distantInfo;
      cwipi_located_point_info_t localInfo = _locationToDistantMesh._locationInfo;
      MPI_Status MPIStatus;

      const bool isRootRank = (currentRank == 0 ||
                               currentRank == _coupledApplicationBeginningRankCouplingComm +
                               _coupledApplicationNRankCouplingComm);

      int rootRank;
      if (isRootRank) {

        MPI_Sendrecv((int*) &localInfo,   1, MPI_INT,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     (int*) &distantInfo, 1, MPI_INT,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     _couplingComm, &MPIStatus);
      }

      //
      // application rootRank send the value of distantInfo to the other applications in the couplingComm

      int sizeComm;
      MPI_Comm_size(_couplingComm, &sizeComm);
      MPI_Status status;

      if ( _coupledApplicationBeginningRankCouplingComm == 0)
        rootRank = _coupledApplicationNRankCouplingComm;
      else
        rootRank = 0;

      if (! isRootRank)
        MPI_Recv(&distantInfo, 1, MPI_INT, rootRank, 0, _couplingComm, &status);
      else {
        int rank1;
        int rank2;

        if (rootRank == 0) {
          rank1 = 1;
          rank2 = _coupledApplicationBeginningRankCouplingComm;
        }
        else {
          rank1 = _coupledApplicationNRankCouplingComm+1;
          rank2 = sizeComm;
        }

        for (int i = rank1; i < rank2; i++)
          MPI_Send(&distantInfo, 1, MPI_INT, i, 0, _couplingComm);

      }

      //
      // Exchange info about distant mesh

      if (distantInfo == CWIPI_DISTANT_MESH_INFO ||
          localInfo == CWIPI_DISTANT_MESH_INFO) {

        //
        // Location

        if (localInfo == CWIPI_DISTANT_MESH_INFO)
          _locationToDistantMesh._elementContaining =  (int *) malloc (sizeof(int) * (_locationToDistantMesh._nLocatedPoint));

        int *pLocation = NULL;
        if (distantInfo == CWIPI_DISTANT_MESH_INFO)
          pLocation = _location;

        fvmc_locator_exchange_point_var(_fvmLocator,
                                        (void *) pLocation,
                                        (void *) _locationToDistantMesh._elementContaining,
                                        NULL,
                                        sizeof(int),
                                        1,
                                        0);

        //
        // elementContainingMPIrankContaining

        if (localInfo == CWIPI_DISTANT_MESH_INFO)
          _locationToDistantMesh._elementContainingMPIrankContaining =
             (int *) malloc (sizeof(int) * (_locationToDistantMesh._nLocatedPoint));

        int *pMPIrank = NULL;
        std::vector <int> *MPIrank = NULL;

        if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
          int rankInGlobalComm;
          MPI_Comm_rank(_localApplicationProperties.getGlobalComm(), &rankInGlobalComm);
          MPIrank = new std::vector <int> (_nDistantPoint, rankInGlobalComm);
          pMPIrank = &((*MPIrank)[0]);
        }

        fvmc_locator_exchange_point_var(_fvmLocator,
                                        (void *) pMPIrank,
                                        (void *) _locationToDistantMesh._elementContainingMPIrankContaining,
                                        NULL,
                                        sizeof(int),
                                        1,
                                        0);

        if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
          delete MPIrank;
        }

        //
        // ElementContainingNVertex

        if (localInfo == CWIPI_DISTANT_MESH_INFO)
          _locationToDistantMesh._elementContainingNVertex =  (int *) malloc (sizeof(int) * (_locationToDistantMesh._nLocatedPoint + 1));

        int *p_nVertex = NULL;
        _maxElementContainingNVertex = 0;

        if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
          _nVertex = new std::vector <int> (_nDistantPoint, 0);
          std::vector <int> & _nVertexRef = *_nVertex;

          for (int i = 0; i < _nDistantPoint; i++)
            _nVertexRef[i] = _supportMesh->getEltConnectivityIndex()[_location[i]] -
              _supportMesh->getEltConnectivityIndex()[_location[i]-1];

          p_nVertex = &(*_nVertex)[0];

          int localMaxElementContainingNVertex = *std::max_element(_nVertexRef.begin(), _nVertexRef.end());
          MPI_Allreduce (&localMaxElementContainingNVertex,
                         &_maxElementContainingNVertex,
                         1, MPI_INT, MPI_MAX,
                         _localApplicationProperties.getLocalComm());
        }

        int distantMaxElementContainingNVertex = 0;
        if (isRootRank) {
          MPI_Sendrecv((int*) &_maxElementContainingNVertex,   1, MPI_INT,
                       _coupledApplicationBeginningRankCouplingComm, 0,
                       (int*) &distantMaxElementContainingNVertex, 1, MPI_INT,
                       _coupledApplicationBeginningRankCouplingComm, 0,
                       _couplingComm, &MPIStatus);
        }

        MPI_Bcast(&distantMaxElementContainingNVertex, 1, MPI_INT, 0, _localApplicationProperties.getLocalComm());

        _maxElementContainingNVertex = std::max(_maxElementContainingNVertex, distantMaxElementContainingNVertex);

        fvmc_locator_exchange_point_var(_fvmLocator,
                                        (void *) p_nVertex,
                                        (void *) _locationToDistantMesh._elementContainingNVertex,
                                        NULL,
                                        sizeof(int),
                                        1,
                                        0);

        if (localInfo == CWIPI_DISTANT_MESH_INFO) {
          int previous1 = _locationToDistantMesh._elementContainingNVertex[0];
          _locationToDistantMesh._elementContainingNVertex[0] = 0;
          for (int i = 1; i < _locationToDistantMesh._nLocatedPoint + 1; i++) {
            int previous2 = _locationToDistantMesh._elementContainingNVertex[i];
            _locationToDistantMesh._elementContainingNVertex[i] = previous1 + _locationToDistantMesh._elementContainingNVertex[i-1];
            previous1 = previous2;
          }
        }

        //
        // ElementContainingVertex

        int *tmpLocal = NULL;
        int *tmpDistant = NULL;

        if (localInfo == CWIPI_DISTANT_MESH_INFO) {
          tmpDistant =  (int *) malloc (sizeof(int) * (_locationToDistantMesh._nLocatedPoint));
          _locationToDistantMesh._elementContainingVertex =  (int *) malloc (sizeof(int) * (_locationToDistantMesh._elementContainingNVertex[_locationToDistantMesh._nLocatedPoint]));
        }

        if (distantInfo == CWIPI_DISTANT_MESH_INFO)
          tmpLocal =  (int *) malloc (sizeof(int) * (_nDistantPoint));

        for (int i = 0; i < _maxElementContainingNVertex; i++) {
          if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
            std::vector <int> & _nVertexRef = *_nVertex;
            for (int j = 0; j < _nDistantPoint; j++) {
              if (i < _nVertexRef[j]) {
                int index = _supportMesh->getEltConnectivityIndex()[_location[j]-1] + i;
                tmpLocal[j] = _supportMesh->getEltConnectivity()[index];
              }
              else
                tmpLocal[j] = -1;
            }
          }

          fvmc_locator_exchange_point_var(_fvmLocator,
                                          (void *) tmpLocal,
                                          (void *) tmpDistant,
                                          NULL,
                                          sizeof(int),
                                          1,
                                          0);

          if (localInfo == CWIPI_DISTANT_MESH_INFO) {
            for (int j = 0; j < _locationToDistantMesh._nLocatedPoint; j++) {
              int nVertices = _locationToDistantMesh._elementContainingNVertex[j+1] - _locationToDistantMesh._elementContainingNVertex[j];
              if (i < nVertices) {
                assert(tmpDistant[j] != -1);
                _locationToDistantMesh._elementContainingVertex[_locationToDistantMesh._elementContainingNVertex[j] + i] = tmpDistant[j];
              }
            }
          }
        }

        if (localInfo == CWIPI_DISTANT_MESH_INFO)
          free ( tmpDistant);

        if (distantInfo == CWIPI_DISTANT_MESH_INFO)
          free ( tmpLocal);

        //
        // ElementContainingBarycentricCoordinates

        double *tmpLocal1 = NULL;
        double *tmpDistant1 = NULL;

        if (localInfo == CWIPI_DISTANT_MESH_INFO) {
          tmpDistant1 =  (double *) malloc (sizeof(double) * (_locationToDistantMesh._nLocatedPoint));
          _locationToDistantMesh._elementContainingBarycentricCoordinates =  (double *) malloc (sizeof(double) * (_locationToDistantMesh._elementContainingNVertex[_locationToDistantMesh._nLocatedPoint]));
        }

        if (distantInfo == CWIPI_DISTANT_MESH_INFO)
          tmpLocal1 =  (double *) malloc (sizeof(double) * (_nDistantPoint));

        for (int i = 0; i < _maxElementContainingNVertex; i++) {
          std::vector <double> &  _refBarycentricCoordinates = *_barycentricCoordinates;

          if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
            std::vector <int> & _nVertexRef = *_nVertex;
            for (int j = 0; j < _nDistantPoint; j++) {
              if (i < _nVertexRef[j]) {
                int index = _supportMesh->getEltConnectivityIndex()[_location[j]-1] + i;
                tmpLocal1[j] = _refBarycentricCoordinates[index];
              }
            }
          }

          fvmc_locator_exchange_point_var(_fvmLocator,
                                          (void *) tmpLocal1,
                                          (void *) tmpDistant1,
                                          NULL,
                                          sizeof(double),
                                          1,
                                          0);

          if (localInfo == CWIPI_DISTANT_MESH_INFO) {
            for (int j = 0; j < _locationToDistantMesh._nLocatedPoint; j++) {
              int nVertices = _locationToDistantMesh._elementContainingNVertex[j+1] - _locationToDistantMesh._elementContainingNVertex[j];
              if (i < nVertices) {
                _locationToDistantMesh._elementContainingBarycentricCoordinates[_locationToDistantMesh._elementContainingNVertex[j] + i] = tmpDistant1[j];
              }
            }
          }
        }

        if (localInfo == CWIPI_DISTANT_MESH_INFO)
          free ( tmpDistant1);

        if (distantInfo == CWIPI_DISTANT_MESH_INFO)
          free ( tmpLocal1);

        // TODO: Optimisation a réaliser dans fvm en mettant un stride[dim]
        //       pour l'instant on calcule le max des nombre de sommets
        //       ou bien creer son propre graphe de communication MPI

        //
        // ElementContainingVertexCoords
        tmpLocal1 = NULL;
        tmpDistant1 = NULL;

        int stride = 3;
        if (localInfo == CWIPI_DISTANT_MESH_INFO) {
          tmpDistant1 =  (double *) malloc (sizeof(double) * (stride * _locationToDistantMesh._nLocatedPoint));
          _locationToDistantMesh._elementContainingVertexCoords =  (double *) malloc (sizeof(double) * (stride * _locationToDistantMesh._elementContainingNVertex[_locationToDistantMesh._nLocatedPoint]));
        }

        if (distantInfo == CWIPI_DISTANT_MESH_INFO)
          tmpLocal1 =  (double *) malloc (sizeof(double) * (stride * _nDistantPoint));

        for (int i = 0; i < _maxElementContainingNVertex; i++) {
          std::vector <int> & _nVertexRef = *_nVertex;
          if (distantInfo == CWIPI_DISTANT_MESH_INFO) {
            for (int j = 0; j < _nDistantPoint; j++) {
              if (i < _nVertexRef[j]) {
                int index = _supportMesh->getEltConnectivityIndex()[_location[j]-1] + i;
                int icel = _supportMesh->getEltConnectivity()[index];
                for (int k = 0; k < stride; k++)
                  tmpLocal1[stride *j + k] = _supportMesh->getVertexCoords()[stride * (icel - 1) + k];
              }
            }
          }

          fvmc_locator_exchange_point_var(_fvmLocator,
                                          (void *) tmpLocal1,
                                          (void *) tmpDistant1,
                                          NULL,
                                          sizeof(double),
                                          stride,
                                          0);

          if (localInfo == CWIPI_DISTANT_MESH_INFO) {
            for (int j = 0; j < _locationToDistantMesh._nLocatedPoint; j++) {
              int nVertices = _locationToDistantMesh._elementContainingNVertex[j+1] - _locationToDistantMesh._elementContainingNVertex[j];
              if (i < nVertices) {
                for (int k = 0; k < stride; k++) {
                  _locationToDistantMesh._elementContainingVertexCoords[stride * (_locationToDistantMesh._elementContainingNVertex[j] + i) + k] =
                    tmpDistant1[stride * j + k];
                }
              }
            }
          }
        }

        if (localInfo == CWIPI_DISTANT_MESH_INFO)
          free ( tmpDistant1);

        if (distantInfo == CWIPI_DISTANT_MESH_INFO)
          free ( tmpLocal1);
      }

      // TODO: Attention la fonction synchronise ne devrait-elle pas etre dans le if precedent !

    }


    _locationToDistantMesh.synchronize();
    _toLocate = false;

  }
  //  _locationToDistantMesh.synchronize();

}

///
/// \brief Exchange field on vertices of cells that contain each located points
///

void LocationToLocalMesh::exchangeCellVertexFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride)
{

  // TODO: Pas optimal : une fois FVM amélioré ou remplacé : transmettre le champ  a une liste de sommets et plus sur la conectivité
  // TODO: exchangeCellVertexFieldOfElementContaining : doute pour un fonctionnement en parallele sans partitionnement
  double *tmpDistant1 = NULL;
  double *tmpLocal1 = NULL;

  if (receivingField != NULL) {
    tmpDistant1 =  (double *) malloc (sizeof(double) * (stride * _locationToDistantMesh._nLocatedPoint));
    _locationToDistantMesh._elementContainingVertexCoords =  (double *) malloc (sizeof(double) * (stride * _locationToDistantMesh._elementContainingNVertex[_locationToDistantMesh._nLocatedPoint]));
  }

  if (sendingField != NULL)
    tmpLocal1 =  (double *) malloc (sizeof(double) * (stride * _nDistantPoint));

  for (int i = 0; i < _maxElementContainingNVertex; i++) {
    std::vector <int> & _nVertexRef = *_nVertex;
    if (sendingField != NULL) {
      for (int j = 0; j < _nDistantPoint; j++) {
        if (i < _nVertexRef[j]) {
          int index = _supportMesh->getEltConnectivityIndex()[_location[j]-1] + i;
          int ivertex = _supportMesh->getEltConnectivity()[index];
          for (int k = 0; k < stride; k++)
            tmpLocal1[stride *j + k] = sendingField[stride*(ivertex - 1) + k];
        }
      }
    }

    fvmc_locator_exchange_point_var(_fvmLocator,
                                   (void *) tmpLocal1,
                                   (void *) tmpDistant1,
                                   NULL,
                                   sizeof(double),
                                   stride,
                                   0);

    if (_locationToDistantMesh._elementContainingNVertex != NULL) {
      for (int j = 0; j < _locationToDistantMesh._nLocatedPoint; j++) {
        int nVertices = _locationToDistantMesh._elementContainingNVertex[j+1] - _locationToDistantMesh._elementContainingNVertex[j];
        if (i < nVertices) {
          for (int k = 0; k < stride; k++)
            receivingField[stride * (_locationToDistantMesh._elementContainingNVertex[j] + i) + k] = tmpDistant1[stride * j + k];
        }
      }
    }
  }

  free ( tmpLocal1);
  free ( tmpDistant1);
}

///
/// \brief Exchange field on cells that contain each located points
///

void LocationToLocalMesh::exchangeCellCenterFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride)
{

  // TODO: exchangeCellCenterFieldOfElementContaining : doute pour un fonctionnement en parallele sans partitionnement

  fvmc_locator_exchange_point_var(_fvmLocator,
                                 (void *) sendingField,
                                 (void *) receivingField,
                                 NULL,
                                 sizeof(double),
                                 stride,
                                 0);
}

///
/// \brief Compute Mean Values
///
///

void  LocationToLocalMesh::midplaneProjection
(
 const int     nbr_som_fac,
 double *const coo_som_fac,
 double *const coo_point_dist
)

{

  int    icoo;
  int    isom;
  int    itri;

  double   cost;
  double   sint;

  double   coo_tmp;

  double   vect1[3];
  double   vect2[3];

  double  prod_vect[3];

  double  barycentre_fac[3];
  double  normale_fac[3];

  double *coo_som_fac_tmp = NULL;
  double  coo_point_dist_tmp[3];

  const double eps = 1e-15;


  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Calcul des coordonnées du barycentre B du polygone P */
  /*======================================================*/

  for (icoo = 0; icoo < 3; icoo++) {

    barycentre_fac[icoo] = 0.;

    for (isom = 0; isom < nbr_som_fac; isom++)
      barycentre_fac[icoo] += coo_som_fac[3*isom+icoo];

    barycentre_fac[icoo] /= nbr_som_fac;

  }

  for (icoo = 0; icoo < 3; icoo++)
    normale_fac[icoo] = 0.;

  /* Calcul de la normale */
  /*======================*/

  for (itri = 0; itri < nbr_som_fac; itri++) {

    for (icoo = 0; icoo < 3; icoo++) {

      vect1[icoo] = coo_som_fac[3*itri+icoo] - barycentre_fac[icoo];

      if (itri < nbr_som_fac - 1)
        vect2[icoo] = coo_som_fac[3*(itri+1)+icoo] - barycentre_fac[icoo];
      else
        vect2[icoo] = coo_som_fac[icoo] - barycentre_fac[icoo];

    }

    normale_fac[0] += vect1[1] * vect2[2] - vect2[1] * vect1[2];
    normale_fac[1] += vect2[0] * vect1[2] - vect1[0] * vect2[2];
    normale_fac[2] += vect1[0] * vect2[1] - vect2[0] * vect1[1];

  }

  /* Projection dans un plan parallèle à la face */
  /*=============================================*/

  /* On ramène l'origine au centre de gravité de la fac */

  for (isom = 0; isom < nbr_som_fac; isom++) {
    for (icoo = 0; icoo < 3; icoo++)
      coo_som_fac[3*isom+icoo] -= barycentre_fac[icoo];
  }

  for (icoo = 0; icoo < 3; icoo++)
    coo_point_dist[icoo] -= barycentre_fac[icoo];


  if (fabs(normale_fac[0]) > eps
      || fabs(normale_fac[1]) > eps) {

    /* Première rotation d'axe (Oz) et d'angle (Ox, proj normale sur Oxy) */

    coo_som_fac_tmp =  (double *) malloc (sizeof(double) * (3 * nbr_som_fac));

    vect1[0] = 1.;
    vect1[1] = 0.;
    vect1[2] = 0.;

    vect2[0] = normale_fac[0];
    vect2[1] = normale_fac[1];
    vect2[2] = 0.;

    crossProduct(vect1, vect2, prod_vect);

    cost = dotProduct(vect1, vect2) / norm(vect2);

    if (prod_vect[2] > 0.)
      sint =  norm(prod_vect) / norm(vect2);
    else
      sint = -norm(prod_vect) / norm(vect2);

    for (isom = 0; isom < nbr_som_fac; isom++) {

      coo_som_fac_tmp[3*isom] =
         cost*coo_som_fac[3*isom] + sint*coo_som_fac[3*isom+1];
      coo_som_fac_tmp[3*isom+1] =
        -sint*coo_som_fac[3*isom] + cost*coo_som_fac[3*isom+1];
      coo_som_fac_tmp[3*isom+2] = coo_som_fac[3*isom+2];

    }

    coo_point_dist_tmp[0] = cost*coo_point_dist[0] + sint*coo_point_dist[1];
    coo_point_dist_tmp[1] = -sint*coo_point_dist[0] + cost*coo_point_dist[1];
    coo_point_dist_tmp[2] = coo_point_dist[2];

    /* Deuxième rotation d'axe (Oy) et d'angle (Oz', proj normale sur Ox'z) */

    vect1[0] =  0.;
    vect1[1] =  0.;
    vect1[2] =  1.;

    vect2[0] =
      sqrt(normale_fac[0]*normale_fac[0] + normale_fac[1]*normale_fac[1]);
    vect2[1] = 0.;
    vect2[2] = normale_fac[2];

    crossProduct(vect1, vect2, prod_vect);

    cost = dotProduct(vect1, vect2) / norm(vect2);

    if (prod_vect[2] > 0.)
      sint =  norm(prod_vect) / norm(vect2);
    else
      sint = -norm(prod_vect) / norm(vect2);


    for (isom = 0; isom < nbr_som_fac; isom++) {

      coo_som_fac[3*isom]   = cost*coo_som_fac_tmp[3*isom] + sint*coo_som_fac_tmp[3*isom + 2];
      coo_som_fac[3*isom+1] = coo_som_fac_tmp[3*isom+1];
      coo_som_fac[3*isom+2] = 0.;

    }

    coo_point_dist[0] = cost*coo_point_dist_tmp[0] + sint*coo_point_dist_tmp[2];
    coo_point_dist[1] = coo_point_dist_tmp[1];
    coo_point_dist[2] = 0.;


    free ( (coo_som_fac_tmp));

  }
  else {

    /* On écrase seulement la coordonnée z du sommet, en intervertissant
       éventuellement les coordonnées dans le plan de projection(Oxy).  */

    if (normale_fac[2] > 0.) {
      for (isom = 0; isom < nbr_som_fac; isom++)
        coo_som_fac[3*isom+2] = 0.;

      coo_point_dist[2] = 0.;
    }
    else {
      for (isom = 0; isom < nbr_som_fac; isom++) {
        coo_tmp = coo_som_fac[3*isom];
        coo_som_fac[3*isom] = coo_som_fac[3*isom+1];
        coo_som_fac[3*isom+1] = coo_tmp;
        coo_som_fac[3*isom+2] = 0.;
      }
      coo_tmp = coo_point_dist[0];
      coo_point_dist[0] = coo_point_dist[1];
      coo_point_dist[1] = coo_tmp;
      coo_point_dist[2] = 0.;
    }
  }
}

///
/// \brief Compute Mean Values
///

void LocationToLocalMesh::compute2DMeanValues()
{
  const int n_dist_points = fvmc_locator_get_n_dist_points(_fvmLocator);
  const fvmc_lnum_t *dist_locations = fvmc_locator_get_dist_locations(_fvmLocator);
  const fvmc_coord_t *dist_coords = fvmc_locator_get_dist_coords(_fvmLocator);
  // const fvmc_coord_t *proj_coords = fvmc_locator_get_dist_projected_coords(_fvmLocator);

  const int *meshConnectivityIndex = _supportMesh->getEltConnectivityIndex();
  const int *meshConnectivity = _supportMesh->getEltConnectivity();
  const double *meshVertexCoords = _supportMesh->getVertexCoords();

  const int order = _supportMesh->getOrder();

  _barycentricCoordinatesIndex = new std::vector <int> (n_dist_points + 1);
  std::vector <int>& nDistBarCoords = *_barycentricCoordinatesIndex;

  //
  // Allocation
  //

  nDistBarCoords[0] = 0;
  for (int ipoint =  0; ipoint < n_dist_points; ipoint++ ) {
    int ielt = dist_locations[ipoint] - 1;

    int nbr_som_fac =  meshConnectivityIndex[ielt+1] -
                       meshConnectivityIndex[ielt];
    nDistBarCoords[ipoint+1] = nDistBarCoords[ipoint] + nbr_som_fac;
  }

  _barycentricCoordinates = new std::vector <double> (nDistBarCoords[n_dist_points]);
  std::vector <double>& distBarCoords = *_barycentricCoordinates;

  if (order == -1) {

    computePolygonMeanValues(n_dist_points,
                             dist_locations,
                             dist_coords,
                             meshConnectivityIndex,
                             meshConnectivity,
                             meshVertexCoords,
                             nDistBarCoords,
                             distBarCoords);
  }

  else {

    const double *uvw = fvmc_locator_get_dist_uvw (_fvmLocator);

    int max_entity_dim = fvmc_nodal_get_max_entity_dim (&(_supportMesh->getFvmNodal()));

    for (int ipoint = 0; ipoint < n_dist_points; ipoint++) {

      int ielt = dist_locations[ipoint] - 1;

      int n_node = meshConnectivityIndex[ielt+1] - meshConnectivityIndex[ielt];

      double *_weights = &(distBarCoords[0]) + nDistBarCoords[ipoint];
      const double *_uvw_point_fvm = uvw + ipoint * max_entity_dim;

      const int n_node_tria = (order+1)*(order+2)/2;

      fvmc_element_t type;

      if (n_node == n_node_tria) {
        type = FVMC_FACE_TRIA;
      }
      else {
        type = FVMC_FACE_QUAD;
      }

      FVMC_ho_basis (type, order, n_node, 1, _uvw_point_fvm, _weights);

    }
  }
}

///
/// \brief Compute Polygo Mean Values
///

void LocationToLocalMesh::computePolygonMeanValues(const int           n_dist_points,
                                                   const fvmc_lnum_t  *dist_locations,
                                                   const fvmc_coord_t *dist_coords,
                                                   const int          *meshConnectivityIndex,
                                                   const int          *meshConnectivity,
                                                   const double       *meshVertexCoords,
                                                   const std::vector <int>& nDistBarCoords,
                                                   std::vector <double>& distBarCoords)

{
  /* Boucle sur les points distants */

  fvmc_coord_t coo_point_dist[3];

  /* Tableaux locaux */

  const double eps_base = 1e-10;
  std::vector <double> coo_som_fac;
  std::vector <double> s;
  std::vector <double> dist;
  std::vector <double> aire;
  std::vector <double> proScal;

  for (int ipoint =  0; ipoint < n_dist_points; ipoint++ ) {

    double *_distBarCoords = &(distBarCoords[nDistBarCoords[ipoint]]);

    /* Initialisation - Copie locale */

    int isOnEdge = 0;
    int isVertex = 0;
    int ielt = dist_locations[ipoint] - 1;

    int nbr_som_fac =  meshConnectivityIndex[ielt+1] -
                       meshConnectivityIndex[ielt];
    coo_point_dist[0] = dist_coords[3*ipoint];
    coo_point_dist[1] = dist_coords[3*ipoint + 1];
    coo_point_dist[2] = dist_coords[3*ipoint + 2];

    if (ipoint == 0) {
      coo_som_fac.resize(3 * nbr_som_fac);
      s.resize(3 * nbr_som_fac);
      dist.resize(nbr_som_fac);
      aire.resize(nbr_som_fac);
      proScal.resize(nbr_som_fac);
    }
    else {
      if (proScal.size() < nbr_som_fac) {
        coo_som_fac.resize(3 * nbr_som_fac);
        s.resize(3 * nbr_som_fac);
        dist.resize(nbr_som_fac);
        aire.resize(nbr_som_fac);
        proScal.resize(nbr_som_fac);
      }
    }

    for (int isom = 0; isom < nbr_som_fac; isom++) {
      coo_som_fac[3*isom]   =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)];

      coo_som_fac[3*isom+1] =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)+1];

      coo_som_fac[3*isom+2] =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)+2];
    }

    /* Projection sur un plan moyen */

    double bary[3];
    _computeBary (nbr_som_fac, &(coo_som_fac[0]), bary);

    double n[3] = {0, 0, 1};
    _computeNormal (nbr_som_fac, &(coo_som_fac[0]), n);

    _project_point2(coo_point_dist, bary, n, coo_point_dist);

    for (int isom = 0; isom < nbr_som_fac; isom++) {

      double *pt1 = &(coo_som_fac[0]) + 3 *isom;
      _project_point2 (pt1, bary, n, pt1);

    }

    double bounds[6] = {DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX};

    for (int isom = 0; isom < nbr_som_fac; isom++) {
      bounds[0] = std::min(bounds[0], coo_som_fac[3*isom]);
      bounds[1] = std::max(bounds[1], coo_som_fac[3*isom]);

      bounds[2] = std::min(bounds[2], coo_som_fac[3*isom + 1]);
      bounds[3] = std::max(bounds[3], coo_som_fac[3*isom + 1]);

      bounds[4] = std::min(bounds[4], coo_som_fac[3*isom + 2]);
      bounds[5] = std::max(bounds[5], coo_som_fac[3*isom + 2]);
    }


    /* Verification que le point est dans l'element */

    if (fvmc_point_in_polygon (coo_point_dist,
                               nbr_som_fac,
                               &(coo_som_fac[0]),
                               bounds,
                               n) != 1) {

      double closestPoint[3];
      double dist_min = DBL_MAX;

      for (int k = 0; k < nbr_som_fac; k++) {
        double *p1 = &(coo_som_fac[3 * k]);
        double *p2 = &(coo_som_fac[3 * ((k+1) % nbr_som_fac)]);
        double closest[3];
        double t;

        double dist2 = fvmc_distance_to_line (coo_point_dist,
                                             p1,
                                             p2,
                                             &t,
                                             closest);

        if (dist2 < dist_min) {
          dist_min = dist2;
          closestPoint[0] = closest[0];
          closestPoint[1] = closest[1];
          closestPoint[2] = closest[2];
        }
      }

      coo_point_dist[0] = closestPoint[0];
      coo_point_dist[1] = closestPoint[1];
      coo_point_dist[2] = closestPoint[2];

    }

    /* Calcul des coordonnnees barycentriques */

    double min_dist = DBL_MAX;
    for (int isom = 0; isom < nbr_som_fac; isom++) {

      int inext = (isom + 1) % nbr_som_fac;
      double *vect = &s[0] + 3*isom;
      double l_edge;
      vect[0] = coo_som_fac[3*inext]   - coo_som_fac[3*isom];
      vect[1] = coo_som_fac[3*inext+1] - coo_som_fac[3*isom+1];
      vect[2] = coo_som_fac[3*inext+2] - coo_som_fac[3*isom+2];
      l_edge  = _MODULE (vect);

      min_dist = std::min(l_edge, min_dist);
    }

    double eps = std::max(min_dist * eps_base, 1.e-30);

    for (int isom = 0; isom < nbr_som_fac; isom++) {

      double *vect = &s[0] + 3*isom;
      vect[0] = coo_som_fac[3*isom]   - coo_point_dist[0];
      vect[1] = coo_som_fac[3*isom+1] - coo_point_dist[1];
      vect[2] = coo_som_fac[3*isom+2] - coo_point_dist[2];
      dist[isom] = _MODULE (vect);

    }

    int currentVertex;
    for (int isom = 0; isom < nbr_som_fac; isom++) {
      int inext = (isom + 1) % nbr_som_fac;
      double *vect1 = &s[0] + 3 * isom;
      double *vect2 = &s[0] + 3 * inext;
      double pvect[3];

      proScal[isom] = _DOT_PRODUCT (vect1, vect2);
      _CROSS_PRODUCT(pvect, vect1, vect2);

      double sign = _DOT_PRODUCT (pvect, n);
      aire[isom] = _MODULE(pvect);

      if (sign < 0) {
        aire[isom] = -aire[isom];
      }

      if (dist[isom] <= eps) {

        isVertex = 1;
        currentVertex = isom;
        break;
      }

      else if ((fabs(aire[isom]) <= eps)  && (proScal[isom] < 0)) {

        isOnEdge = 1;
        currentVertex = isom;
        break;

      }

    }

    /* Le point distant est un sommet */

    if (isVertex) {
      for (int isom = 0; isom < nbr_som_fac; isom++)
        _distBarCoords[isom] = 0.;
      _distBarCoords[currentVertex] = 1.;
    }

    /* Le point distant est sur arete */

    else if (isOnEdge) {

      for (int isom = 0; isom < nbr_som_fac; isom++)
        _distBarCoords[isom] = 0.;

      int nextPoint = (currentVertex + 1) % nbr_som_fac;

      _distBarCoords[currentVertex] =
        dist[nextPoint]     / (dist[nextPoint]+dist[currentVertex]);
      _distBarCoords[nextPoint]     =
        dist[currentVertex] / (dist[nextPoint]+dist[currentVertex]);

    }

    /* Cas general */

    else {

      double sigma = 0;
      for (int isom = 0; isom < nbr_som_fac; isom++) {
        double coef = 0.;
        int previousVertex = (isom - 1 + nbr_som_fac) % nbr_som_fac;
        int nextVertex = (isom + 1) % nbr_som_fac;

        if (fabs(aire[previousVertex]) > eps)
          coef += (dist[previousVertex] - proScal[previousVertex]/dist[isom]) / aire[previousVertex];
        if (fabs(aire[isom]) > eps)
          coef += (dist[nextVertex]     - proScal[isom]/dist[isom])           / aire[isom];
        sigma += coef;
        _distBarCoords[isom] = coef;

      }

      if (fabs(sigma) >= eps ) {
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          _distBarCoords[isom] /= sigma;
        }
      }

      else {

        double abs_sigma = fabs(sigma);
        printf("Warning : Mise à NAN %f %f\n", abs_sigma,  eps);
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          _distBarCoords[isom] = NAN;
        }
      }

      /* Check Result */

      for (int isom = 0; isom <  nbr_som_fac; isom++) {
        if ( _distBarCoords[isom] != _distBarCoords[isom] ||
             _distBarCoords[isom] < 0. ||
             _distBarCoords[isom] > 1. ) {

          double dist_min = DBL_MAX;
          int k_min = 0;
          double t_min;

          for (int k = 0; k < nbr_som_fac; k++) {
            _distBarCoords[k] = 0.0;
          }

          for (int k = 0; k < nbr_som_fac; k++) {
            double *p1 = &(coo_som_fac[3 * k]);
            double *p2 = &(coo_som_fac[3 * ((k+1) % nbr_som_fac)]);
            double closest[3];
            double t;

            double dist2 = fvmc_distance_to_line (coo_point_dist,
                                                 p1,
                                                 p2,
                                                 &t,
                                                 closest);
            if (dist2 < dist_min) {
              t_min = t;
              k_min = k;
            }
          }

          _distBarCoords[k_min] = 1 - t_min;
          _distBarCoords[(k_min + 1) % nbr_som_fac] = t_min;

          break;

        }

      }

    }

    if (0 == 1) {
      if ((n_dist_points == 1) && (dist_locations[0] == 1)) {

        bftc_printf("coord %i %i :", ipoint+1, ielt+1);
        bftc_printf(" %12.5e %12.5e %12.5e", dist_coords[3*ipoint],
                    dist_coords[3*ipoint+1],
                    dist_coords[3*ipoint+2] );
        bftc_printf("\n");

        bftc_printf("coo b %i :", ipoint+1);
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          bftc_printf(" %f", _distBarCoords[isom]);
        }
        bftc_printf("\n");
      }
    }
  }

  coo_som_fac.clear();
  s.clear();
  aire.clear();
  dist.clear();
  proScal.clear();
}


void LocationToLocalMesh::compute3DMeanValues()
{

  const int n_dist_points           = fvmc_locator_get_n_dist_points(_fvmLocator);
  const fvmc_lnum_t *dist_locations = fvmc_locator_get_dist_locations(_fvmLocator);
  const float *dist_distances       = fvmc_locator_get_dist_distances(_fvmLocator);
  const fvmc_coord_t *dist_coords   = fvmc_locator_get_dist_coords(_fvmLocator);
  //const fvmc_coord_t *proj_coords = fvmc_locator_get_dist_projected_coords(_fvmLocator);


  /**** Tableaux barycentriques ****/

  int tailleDistBarCoords =  4 * n_dist_points;

  if (_barycentricCoordinatesIndex == NULL) {
    _barycentricCoordinatesIndex = new std::vector <int> (n_dist_points + 1);
    _barycentricCoordinates = new std::vector <double> (tailleDistBarCoords);
  }
  else {
    _barycentricCoordinatesIndex->resize(n_dist_points + 1);
    _barycentricCoordinates->resize(tailleDistBarCoords);
  }

  std::vector <int>& nDistBarCoords = *_barycentricCoordinatesIndex;
  std::vector <double>& distBarCoords = *_barycentricCoordinates;

  /* Const */

  const int *meshConnectivityIndex = _supportMesh->getEltConnectivityIndex();
  const int *meshConnectivity = _supportMesh->getEltConnectivity();
  const int *polyMeshConnectivityIndex     = NULL;
  const int *polyMeshConnectivity          = NULL;
  if (_supportMesh->getNPolyhedra() > 0) {
    polyMeshConnectivityIndex = &(_supportMesh->getPolyhedraCellToVertexConnectivityIndex()[0]);
    polyMeshConnectivity      = &(_supportMesh->getPolyhedraCellToVertexConnectivity()[0]);
  }

  const double *meshVertexCoords = _supportMesh->getVertexCoords();

  const std::vector<int>& isDegenerated  = _supportMesh->getIsDegenerated();
  const std::vector<double>& characteristicLength  = _supportMesh->getCharacteristicLength();

  nDistBarCoords[0] = 0;

  /* Constantes */

  const int nEltStd = _supportMesh->getNElts() - _supportMesh->getNPolyhedra();

  /* Boucle sur les points distants */

  const int n_face_max = 6;
  const int n_vtx_max = 24;

  int faceConnectivityIndex[n_face_max + 1];
  int faceConnectivity[n_vtx_max];
  double faceVtxCoords[3*n_vtx_max];
  int localConnec[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int localFaceConnectivity[n_vtx_max];
  std::vector <int> nFaceDistBarCoords(2);
  std::vector <double> faceDistBarCoords(8);

  //  for (int ipoint =  0; ipoint < n_dist_points; ipoint++ ) {

  const int order = _supportMesh->getOrder();

  if (order == -1) {

    int ipoint = 0;
    int ipoint_old = -1;

    while (ipoint < n_dist_points) {
      int ielt = dist_locations[ipoint] - 1; // numero de l'element le plus proche du point

      float dist;

      if (ipoint != ipoint_old)

        dist = dist_distances[ipoint]; // numero de l'element le plus proche du point

      else

        dist = 1.1;  //

      ipoint_old = ipoint;

      /* Coordonnees du point distant */

      fvmc_coord_t coo_point_dist[3];

      coo_point_dist[0] = dist_coords[3 * ipoint    ];
      coo_point_dist[1] = dist_coords[3 * ipoint + 1];
      coo_point_dist[2] = dist_coords[3 * ipoint + 2];


      //
      // Adjust table length
      //

      int nbr_som;
      if (ielt < nEltStd)
        nbr_som = meshConnectivityIndex[ielt + 1] - meshConnectivityIndex[ielt];
      else {
        int i = ielt - nEltStd;
        nbr_som = polyMeshConnectivityIndex[i+1] - polyMeshConnectivityIndex[i];
      }

      //
      // If element is degenerated
      //

      nDistBarCoords[ipoint + 1] = nDistBarCoords[ipoint] + nbr_som;

      while (distBarCoords.size() <= nDistBarCoords[ipoint + 1]) {
        distBarCoords.resize(2 * distBarCoords.size());
      }

      if (isDegenerated[ielt]) {

        for (int ivertex = 0; ivertex < nbr_som; ivertex++)
          distBarCoords[nDistBarCoords[ipoint] + ivertex] = 1./nbr_som;

        ipoint++;
      }

      else {

        if (ielt < nEltStd) {

          //
          // Standard element
          //

          int ierr = 0;

          if (dist > 1.) {

            int *_cellConnectivity = const_cast<int *> (&(meshConnectivity[meshConnectivityIndex[ielt]]));

            //
            // If point is not in element : Closest face projection
            //

            int n_face = 0;

            switch(nbr_som){

            case 4 :

              //
              // Tetraedra :
              //

              n_face = 4;

              tetrahedraFaces (1,
                               0,
                               _cellConnectivity,
                               faceConnectivityIndex,
                               faceConnectivity);

              tetrahedraFaces (1,
                               0,
                               localConnec,
                               faceConnectivityIndex,
                               localFaceConnectivity);

              break;

            case 5 :

              //
              // Pyramid
              //

              n_face = 5;

              pyramidFaces (1,
                            0,
                            _cellConnectivity,
                            faceConnectivityIndex,
                            faceConnectivity);

              pyramidFaces (1,
                            0,
                            localConnec,
                            faceConnectivityIndex,
                            localFaceConnectivity);

              break;

            case 6 :

              //
              // Prism
              //

              n_face = 5;

              prismFaces (1,
                          0,
                          _cellConnectivity,
                          faceConnectivityIndex,
                          faceConnectivity);

              prismFaces (1,
                          0,
                          localConnec,
                          faceConnectivityIndex,
                          localFaceConnectivity);

              break;

            case 8 :

              //
              // Hexahedron
              //

              n_face = 6;

              hexahedraFaces (1,
                              0,
                              _cellConnectivity,
                              faceConnectivityIndex,
                              faceConnectivity);

              hexahedraFaces (1,
                              0,
                              localConnec,
                              faceConnectivityIndex,
                              localFaceConnectivity);

              break;

            }

            for (int k = 0; k < n_face; k++) {
              for (int k1 = faceConnectivityIndex[k]; k1 < faceConnectivityIndex[k+1]; k1++) {
                int iVtx = faceConnectivity[k1] - 1;
                faceVtxCoords[3*k1  ] = meshVertexCoords[3*iVtx  ];
                faceVtxCoords[3*k1+1] = meshVertexCoords[3*iVtx+1];
                faceVtxCoords[3*k1+2] = meshVertexCoords[3*iVtx+2];
              }
            }

            double *_faceVtxCoords = faceVtxCoords;
            double dist_min = DBL_MAX;
            int k_min = 0;

            for (int k = 0; k < n_face; k++) {
              const int nVtx = faceConnectivityIndex[k+1] - faceConnectivityIndex[k];
              double closest[3];
              double pcoords[3];
              double dist_face;

              fvmc_polygon_evaluate_Position(coo_point_dist, nVtx, _faceVtxCoords, closest,
                                             pcoords, &dist_face);

              _faceVtxCoords += 3*nVtx;

              if (dist_face < dist_min) {
                dist_min = dist_face;
                k_min = k;
              }

            }

            int dist_locations1 = k_min + 1;
            const int nVtx = faceConnectivityIndex[k_min+1] - faceConnectivityIndex[k_min];

            nFaceDistBarCoords[0] = 0;
            nFaceDistBarCoords[1] = nVtx;


            computePolygonMeanValues (1,
                                      &dist_locations1,
                                      coo_point_dist,
                                      faceConnectivityIndex,
                                      faceConnectivity,
                                      meshVertexCoords,
                                      nFaceDistBarCoords,
                                      faceDistBarCoords);

            int *_localFaceConnectivity = localFaceConnectivity + faceConnectivityIndex[k_min];

            for (int isom = 0; isom < nbr_som; isom++) {
              distBarCoords[nDistBarCoords[ipoint] + isom] = 0.;
            }

            for (int k = 0; k < nVtx; k++) {
              distBarCoords[nDistBarCoords[ipoint] + _localFaceConnectivity[k]] = faceDistBarCoords[k];
            }

            ipoint++;

          }

          else {

            double uvw[3];
            double vertex_coords[8][3];
            double deriv[8][3];

            for (int ivertex = 0; ivertex < nbr_som; ivertex++) {
              vertex_coords[ivertex][0] =
                meshVertexCoords[3 * (meshConnectivity[meshConnectivityIndex[ielt] + ivertex] - 1) ];
              vertex_coords[ivertex][1] =
                meshVertexCoords[3 * (meshConnectivity[meshConnectivityIndex[ielt] + ivertex] - 1) + 1];
              vertex_coords[ivertex][2] =
                meshVertexCoords[3 * (meshConnectivity[meshConnectivityIndex[ielt] + ivertex] - 1) + 2];
            }

            switch(nbr_som){

            case 4 :

              //
              // Tetraedra :
              //

              ierr = compute_uvw(CWIPI_CELL_TETRA4,
                                 coo_point_dist,
                                 vertex_coords,
                                 1e-6,
                                 uvw);

              compute_shapef_3d(CWIPI_CELL_TETRA4,
                                uvw,
                                &distBarCoords[0] + nDistBarCoords[ipoint],
                                deriv);



              break;

            case 5 :
              {
                //
                // Pyramid
                //

                // Check if not a vertex for stability

                int onVtx = 0;
                for (int k1 = 0; k1 < nbr_som; k1++) {
                  double v[3] = {vertex_coords[k1][0] - coo_point_dist[0],
                                 vertex_coords[k1][1] - coo_point_dist[1],
                                 vertex_coords[k1][2] - coo_point_dist[2]};

                  double _dist = _MODULE(v);

                  if (_dist < 1e-6 * _tolerance) {
                    double *_bar = &(distBarCoords[0]) + nDistBarCoords[ipoint];
                    for (int k2 = 0; k2 < nbr_som; k2++) {
                      _bar[k2] = 0.;
                    }
                    _bar[k1] = 1.;
                    onVtx = 1;
                    ierr = 1;
                    break;
                  }
                }

                if (!onVtx) {

                  ierr = compute_uvw(CWIPI_CELL_PYRAM5,
                                     coo_point_dist,
                                     vertex_coords,
                                     1e-6,
                                     uvw);

                  compute_shapef_3d(CWIPI_CELL_PYRAM5,
                                    uvw,
                                    &(distBarCoords[0]) + nDistBarCoords[ipoint],
                                    deriv);
                }
                break;
              }
            case 6 :

              //
              // Prism
              //

              ierr = compute_uvw(CWIPI_CELL_PRISM6,
                                 coo_point_dist,
                                 vertex_coords,
                                 1e-6,
                                 uvw);

              compute_shapef_3d(CWIPI_CELL_PRISM6,
                                uvw,
                                &(distBarCoords[0]) + nDistBarCoords[ipoint],
                                deriv);
              break;

            case 8 :

              //
              // Hexahedron
              //

              ierr = compute_uvw(CWIPI_CELL_HEXA8,
                                 coo_point_dist,
                                 vertex_coords,
                                 1e-6,
                                 uvw);

              compute_shapef_3d(CWIPI_CELL_HEXA8,
                                uvw,
                                &(distBarCoords[0]) + nDistBarCoords[ipoint],
                                deriv);
              break;

            default:
              bftc_error(__FILE__, __LINE__, 0,
                         "compute3DMeanValues: unhandled element type\n");


            }

            if (ierr != 0) {
              ipoint++;
            }
          }
        }

        else {

          //
          // Polyhedron
          //

          const int ipoly = ielt - nEltStd;

          const int *polyCellToFaceIdx   =  _supportMesh->getPolyhedraFaceIndex();
          const int *polyCellToFace      = _supportMesh->getPolyhedraCellToFaceConnectivity();
          const int *polyFaceToVertex    =  _supportMesh->getPolyhedraFaceConnectivity();
          const int *polyFaceToVertexIdx = _supportMesh->getPolyhedraFaceConnectivityIndex();
          const int  n_poly_vertex       =  polyMeshConnectivityIndex[ipoly + 1]
                                          - polyMeshConnectivityIndex[ipoly];
          const int  vertexIdx           =  polyMeshConnectivityIndex[ipoly];
          const int  n_poly_face         =  polyCellToFaceIdx[ipoly + 1] - polyCellToFaceIdx[ipoly];
          const int  faceIdx             =  polyCellToFaceIdx[ipoly];

          std::map <int, int> indirection;
          for (int i = 0; i < n_poly_vertex; i++)
            indirection[polyMeshConnectivity[vertexIdx + i]] = i+1;


          // TODO : A optimiser : faire une allocation unique sur le nombre de faces max
          //  d'un polyedre et le nombre de sommet max par face

          int *faceToVertexEltIdx =  (int *) malloc (sizeof(int) * (n_poly_face + 1));

          faceToVertexEltIdx[0] = 0;
          for (int i = 0; i < n_poly_face; i++) {
            int iface = abs(polyCellToFace[faceIdx + i]) - 1;
            faceToVertexEltIdx[i+1] = faceToVertexEltIdx[i]
              + polyFaceToVertexIdx[iface+1]
              - polyFaceToVertexIdx[iface];
          }

          int *faceToVertexElt =  (int *) malloc (sizeof(int) * (faceToVertexEltIdx[n_poly_face]));
          int *faceDirection =  (int *) malloc (sizeof(int) * (n_poly_face));
          int k = 0;
          for (int i = 0; i < n_poly_face; i++) {
            const int iface          = abs(polyCellToFace[faceIdx + i]) - 1;
            const int direction     =     (polyCellToFace[faceIdx + i] < 0) ? -1 : 1;
            faceDirection[i] = direction;
            for (int j = polyFaceToVertexIdx[iface]; j < polyFaceToVertexIdx[iface+1]; j++) {
              faceToVertexElt[k++] = indirection[polyFaceToVertex[j]];
            }
          }

          double* vertex_coords_Elts =  (double *) malloc (sizeof(double) * (3*n_poly_vertex));
          for (int i = 0; i < n_poly_vertex; i++) {
            int ivertex = polyMeshConnectivity[vertexIdx + i] - 1;
            for (int j = 0; j < 3; j++) {
              vertex_coords_Elts[3 * i + j] = meshVertexCoords[3 * ivertex + j];
            }
          }

          indirection.clear();

          compute3DMeanValuesPoly(coo_point_dist,
                                  n_poly_face,
                                  n_poly_vertex,
                                  faceDirection,
                                  faceToVertexEltIdx,
                                  faceToVertexElt,
                                  vertex_coords_Elts,
                                  characteristicLength[ielt],
                                  dist_distances[ipoint],
                                  &distBarCoords[0] + nDistBarCoords[ipoint]);

          free ( vertex_coords_Elts);
          free ( faceToVertexElt);
          free ( faceToVertexEltIdx);
          free ( faceDirection);

          ipoint++;
        }
      }
    }
  }

  else {

    int max_entity_dim = fvmc_nodal_get_max_entity_dim (&(_supportMesh->getFvmNodal()));

    nDistBarCoords.resize(n_dist_points + 1);

    nDistBarCoords[0] = 0;
    for (int ipoint =  0; ipoint < n_dist_points; ipoint++ ) {
      int ielt = dist_locations[ipoint] - 1;
      int nbr_som_fac =  meshConnectivityIndex[ielt+1] - meshConnectivityIndex[ielt];
      nDistBarCoords[ipoint+1] = nDistBarCoords[ipoint] + nbr_som_fac;
    }

    distBarCoords.resize(nDistBarCoords[n_dist_points]);

    const double *uvw = fvmc_locator_get_dist_uvw (_fvmLocator);

    for (int ipoint = 0; ipoint < n_dist_points; ipoint++) {

      int ielt = dist_locations[ipoint] - 1;

      int n_node = meshConnectivityIndex[ielt+1] - meshConnectivityIndex[ielt];

      double *_weights = &(distBarCoords[0]) + nDistBarCoords[ipoint];
      const double *_uvw_point_fvm = uvw + ipoint * max_entity_dim;

      fvmc_element_t type;

      const int n_node_tetra = (order+1)*(order+2)*(order+3)/6;
      const int n_node_hexa = (order+1)*(order+1)*(order+1);
      //const int n_node_prism = (order+1)*(order+1)*(order+2)/2;
      const int n_node_pyramid = (order+1)*(order+2)*(2*order+3)/6;

      if (n_node == n_node_tetra) {
        type = FVMC_CELL_TETRA;
      }
      else if (n_node == n_node_hexa) {
        type = FVMC_CELL_HEXA;
      }
      else if (n_node == n_node_pyramid) {
        type = FVMC_CELL_PYRAM;
      }
      else {
        type = FVMC_CELL_PRISM;
      }

      FVMC_ho_basis (type, order, n_node, 1, _uvw_point_fvm, _weights);

    }
  }
}


void LocationToLocalMesh::compute3DMeanValuesPoly(const double point_coords[],
                                                  const int    n_poly_faces,
                                                  const int    n_poly_vertex,
                                                  const int    faceDirection[],
                                                  const int    faceToVertexIdx[],
                                                  const int    faceToVertex[],
                                                  const double vertex_coords[],
                                                  const double characteristicLength,
                                                  const float  distElt,
                                                  double       distBarCoords[])
{

  //
  // Polyhedron
  //


  // Mise a jour des tableaux locaux

  std::vector <double> coo_som_face(3 * n_poly_vertex);
  std::vector <double> dist(n_poly_vertex);
  std::vector <double> s(3 * n_poly_vertex);
  std::vector <double> angle(3); //angle[  v(i) v v(i+1) ]
  std::vector <double> normale(9); //normale

  double sigma = 0;

  int isOnFace = 0;

  double eps_loc = geometricEpsilon(characteristicLength, GEOM_EPS_VOL);

  /**** Inialisation du tableau des coordonnees temporaires a 0 ****/

  for (int isom = 0; isom < n_poly_vertex; isom++)
    distBarCoords[isom] = 0.;

  for (int isom = 0; isom < n_poly_vertex; isom++) {
    s[3 * isom    ] = vertex_coords[3 * isom    ] - point_coords[0];
    s[3 * isom + 1] = vertex_coords[3 * isom + 1] - point_coords[1];
    s[3 * isom + 2] = vertex_coords[3 * isom + 2] - point_coords[2];
  }

  if (distElt > 1.) {

    //
    // Search clostest face
    //

    int n_vtx_max = 0;
    for (int i = 0; i < n_poly_faces; i++) {
      const int n_vtx = faceToVertexIdx[i+1] - faceToVertexIdx[i];
      if (n_vtx > n_vtx_max) {
        n_vtx_max = n_vtx;
      }
    }

    double * poly_vertex = (double *) malloc (3 * sizeof(double) * n_vtx_max);

    double dist_min = DBL_MAX;
    int iface = -1;

    for (int i = 0; i < n_poly_faces; i++) {
      const int n_vtx = faceToVertexIdx[i+1] - faceToVertexIdx[i];
      int k = 0;

      for (int j = faceToVertexIdx[i]; j < faceToVertexIdx[i+1]; j++) {
        int i_vertex = faceToVertex[j] - 1;
        poly_vertex[k++] = vertex_coords[3 * i_vertex];
        poly_vertex[k++] = vertex_coords[3 * i_vertex + 1];
        poly_vertex[k++] = vertex_coords[3 * i_vertex + 2];

      }

      double closest[3];
      double pcoords[3];
      double dist_face;

      fvmc_polygon_evaluate_Position(const_cast<double *> (point_coords),
                                                  n_vtx, poly_vertex, closest,
                                                  pcoords, &dist_face);

      if (dist_face < dist_min) {
        dist_min = dist_face;
        iface = i;
      }

    }

    free (poly_vertex);

    //
    // Compute mean value for this face if point is on face
    //

    if (iface >= 0) {

      const int n_face_vertex = faceToVertexIdx[iface + 1]
                              - faceToVertexIdx[iface];

      //
      // Copy

      std::vector<int> distBarCoordsFaceIdx(2);
      distBarCoordsFaceIdx[0] = 0;
      distBarCoordsFaceIdx[1] = n_face_vertex;
      std::vector<double> distBarCoordsFace(n_face_vertex);
      int face_location = iface + 1;

      computePolygonMeanValues(1,
                               &face_location,
                               point_coords,
                               faceToVertexIdx,
                               faceToVertex,
                               vertex_coords,
                               distBarCoordsFaceIdx,
                               distBarCoordsFace);

      for (int j = 0; j < n_poly_vertex; j++)
        distBarCoords[j] = 0.;

      for (int j = 0; j < n_face_vertex; j++) {
        int vertex = faceToVertex[faceToVertexIdx[iface]+j] - 1;
        distBarCoords[vertex] = distBarCoordsFace[j];
      }

      distBarCoordsFaceIdx.clear();
      distBarCoordsFace.clear();

    }

    else {

      //
      // If point is not in a face, get the closest vertex
      //

      double normS = sqrt(s[3*0]*s[3*0] + s[3*0+1]*s[3*0+1] + s[3*0+2]*s[3*0+2]);
      int closestVertex = 0;
      for (int isom = 1; isom < n_poly_vertex; isom++) {

        double nextNormS = sqrt(  s[3*isom]*s[3*isom]
                                  + s[3*isom+1]*s[3*isom+1]
                                  + s[3*isom+2]*s[3*isom+2]);
        if (nextNormS < normS) {
          closestVertex = isom;
          normS = nextNormS;
        }
      }

      distBarCoords[closestVertex] = 1;

    }

  }

  else if (distElt > 0 ) {

    //
    // Check if point is on a face
    //

    const int   n_points  = 1;
    fvmc_lnum_t point_ids = 0;

    //
    // Search clostest face
    //

    fvmc_lnum_t face_location = -1;
    float face_distance = -1.;

    fvmc_point_dist_closest_polygon(3,
                                    n_poly_faces,
                                    faceToVertexIdx,
                                    faceToVertex,
                                    vertex_coords,
                                    n_points,
                                    &point_ids,
                                    point_coords,
                                    &face_location,
                                    &face_distance);

    //
    // Compute mean value for this face if point is on face
    //

    if (face_location > 0 && face_distance < eps_loc) {

      isOnFace = 1;

      const int face_location_idx = face_location - 1;
      const int n_face_vertex = faceToVertexIdx[face_location_idx+1]
                              - faceToVertexIdx[face_location_idx];

      //
      // Copy

      std::vector<int> distBarCoordsFaceIdx(2);
      distBarCoordsFaceIdx[0] = 0;
      distBarCoordsFaceIdx[1] = n_face_vertex;
      std::vector<double> distBarCoordsFace(n_face_vertex);

      computePolygonMeanValues(1,
                               &face_location,
                               point_coords,
                               faceToVertexIdx,
                               faceToVertex,
                               vertex_coords,
                               distBarCoordsFaceIdx,
                               distBarCoordsFace);

      for (int j = 0; j < n_poly_vertex; j++)
        distBarCoords[j] = 0.;

      for (int j = 0; j < n_face_vertex; j++) {
        int vertex = faceToVertex[faceToVertexIdx[face_location_idx]+j] - 1;

        distBarCoords[vertex] = distBarCoordsFace[j];

        distBarCoordsFaceIdx.clear();
        distBarCoordsFace.clear();

      }

    }

    //
    // General alogorithm for point in polyhedron
    //

    if (!isOnFace) {

      for (int isom = 0; isom < n_poly_vertex; isom++) {

        dist[isom] = sqrt(s[3*isom    ] * s[3*isom    ]
                        + s[3*isom + 1] * s[3*isom + 1]
                        + s[3*isom + 2] * s[3*isom + 2]);

        s[3*isom]     /= dist[isom];
        s[3*isom + 1] /= dist[isom];
        s[3*isom + 2] /= dist[isom];

      }

      //
      // Second loop on faces to commpute barycentric coordinates
      //


      std::vector <int> triangle_vertices(9); //Nombre de sommets apres decoupage en triangle

      for (int iface = 0; iface < n_poly_faces; iface++) {

        const int n_vertex_fac = faceToVertexIdx[iface + 1]
                               - faceToVertexIdx[iface];

        const int ind_fac_som = faceToVertexIdx[iface];

        fvmc_triangulate_state_t *fvmc_triangulate = fvmc_triangulate_state_create(n_vertex_fac);

        int triangle_vertice_size = (n_vertex_fac-2) * 3;

        if (triangle_vertices.size() < triangle_vertice_size)
          triangle_vertices.resize(triangle_vertice_size);

        //
        // Face triangulation
        //

        int n_triangles;

        if (n_vertex_fac == 4) {

          n_triangles = fvmc_triangulate_quadrangle(3,
                                                    vertex_coords,
                                                    NULL,
                                                    faceToVertex + ind_fac_som,
                                                    &(triangle_vertices[0]));

        }

        else if (n_vertex_fac > 4) {

          n_triangles = fvmc_triangulate_polygon(3,
                                                 n_vertex_fac,
                                                 vertex_coords,
                                                 NULL,
                                                 faceToVertex + ind_fac_som,
                                                 FVMC_TRIANGULATE_MESH_DEF,
                                                 &(triangle_vertices[0]),
                                                 fvmc_triangulate);
        }

        else {
          n_triangles = 1;
          for (int i = 0; i < 3; i++) {
            triangle_vertices[i] = faceToVertex[ind_fac_som + i];
          }
        }

        //
        // Loop on triangles
        //

        for (int itri = 0; itri < n_triangles; itri++) {

          //
          // Check triangle surface
          //

          const int i = triangle_vertices[3*itri    ] - 1;
          int j, k;
          if (faceDirection[iface] < 0) {
            j = triangle_vertices[3*itri + 1] - 1;
            k = triangle_vertices[3*itri + 2] - 1;
          }
          else {
            j = triangle_vertices[3*itri + 2] - 1;
            k = triangle_vertices[3*itri + 1] - 1;
          }

          const double coo_ijx = vertex_coords[3*j]   - vertex_coords[3*i];
          const double coo_ijy = vertex_coords[3*j+1] - vertex_coords[3*i+1];
          const double coo_ijz = vertex_coords[3*j+2] - vertex_coords[3*i+2];
          const double coo_ikx = vertex_coords[3*k]   - vertex_coords[3*i];
          const double coo_iky = vertex_coords[3*k+1] - vertex_coords[3*i+1];
          const double coo_ikz = vertex_coords[3*k+2] - vertex_coords[3*i+2];

          const double areaTri_ijk = 0.5 * sqrt((coo_ijy * coo_ikz - coo_ijz * coo_iky)
                                                * (coo_ijy * coo_ikz - coo_ijz * coo_iky)
                                                + (coo_ijz * coo_ikx - coo_ijx * coo_ikz)
                                                * (coo_ijz * coo_ikx - coo_ijx * coo_ikz)
                                                + (coo_ijx * coo_iky - coo_ijy * coo_ikx)
                                                * (coo_ijx * coo_iky - coo_ijy * coo_ikx));

          double eps_face = geometricEpsilon(characteristicLength, GEOM_EPS_SURF);

          if (fabs(areaTri_ijk) > eps_face) {

            for (int isom = 0; isom < 3; isom++) {

              int isuiv;
              int iprec;
              double prod_scal;
              double mod;

              if (faceDirection[iface] < 0) {
                iprec = triangle_vertices[3*itri + (isom + 2) % 3] - 1;
                isuiv = triangle_vertices[3*itri + (isom + 1) % 3] - 1;
              }
              else {
                iprec = triangle_vertices[3*itri + (isom + 1) % 3] - 1;
                isuiv = triangle_vertices[3*itri + (isom + 2) % 3] - 1;
              }

              prod_scal = s[3*iprec    ] * s[3*isuiv    ]
                        + s[3*iprec + 1] * s[3*isuiv + 1]
                        + s[3*iprec + 2] * s[3*isuiv + 2];

              angle[isom] = acos(prod_scal); //s

              normale[3 * isom    ] =  s[3*iprec + 1] * s[3*isuiv + 2]
                                     - s[3*iprec + 2] * s[3*isuiv + 1];
              normale[3 * isom + 1] =  s[3*iprec + 2] * s[3*isuiv    ]
                                     - s[3*iprec    ] * s[3*isuiv + 2];
              normale[3 * isom + 2] =  s[3*iprec    ] * s[3*isuiv + 1]
                                     - s[3*iprec + 1] * s[3*isuiv    ];

              /// verifier norm

              mod = sqrt(normale[3*isom    ] * normale[3*isom    ]
                         + normale[3*isom + 1] * normale[3*isom + 1]
                         + normale[3*isom + 2] * normale[3*isom + 2]);

              if (mod <  eps_face) {
                normale[3*isom    ] = 0.;
                normale[3*isom + 1] = 0.;
                normale[3*isom + 2] = 0.;
              }

              else {

                normale[3*isom    ] /= mod;
                normale[3*isom + 1] /= mod;
                normale[3*isom + 2] /= mod;
              }

            }

            for (int isom = 0; isom < 3; isom++) {

              double ps_nij_njk; //a ameliorer
              double ps_nki_njk; //a ameliorer
              double ps_ei_njk;  //a ameliorer

              const int iprec = (isom + 2) % 3;
              const int isuiv = (isom + 1) % 3;

              ps_nij_njk = normale[3 * isom    ] * normale[3 * isuiv    ]
                + normale[3 * isom + 1] * normale[3 * isuiv + 1]
                + normale[3 * isom + 2] * normale[3 * isuiv + 2];

              ps_nki_njk = normale[3 * isom    ] * normale[3 * iprec    ]
                + normale[3 * isom + 1] * normale[3 * iprec + 1]
                + normale[3 * isom + 2] * normale[3 * iprec + 2];

              // ps_ei_njk --> sur la face


              const int ivertex_tri = triangle_vertices[3*itri + isom] - 1;
              ps_ei_njk = s[3*ivertex_tri    ] * normale[3*isom]
                        + s[3*ivertex_tri + 1] * normale[3*isom + 1]
                        + s[3*ivertex_tri + 2] * normale[3*isom + 2];

              // vérifier ps_ei_njk

              if (fabs(ps_ei_njk) >  eps_face) {
                distBarCoords[ivertex_tri] +=
                  (angle[isom] + angle[isuiv] * ps_nij_njk + angle[iprec] * ps_nki_njk)
                  / (2 * ps_ei_njk);
               }

            } // Loop en vertices

          } // Good triangle

        } // Loop on triangles

        fvmc_triangulate_state_destroy(fvmc_triangulate);

      } // Loop on faces

      for (int isom = 0; isom < n_poly_vertex; isom++) {

        distBarCoords[isom] /= dist[isom];
        sigma += distBarCoords[isom];

      }

      for (int isom = 0; isom < n_poly_vertex; isom++)
        distBarCoords[isom] = distBarCoords[isom] / sigma;

    } // End of general algorithm (if (!isonface))

    //
    // Output results
    //


    if (0 == 1) {

      std::vector <double> test(3);

      for (int i = 0; i < 3; i++)
        test[i] = 0;

      for (int isom = 0; isom < n_poly_vertex; isom++){

        test[0] += distBarCoords[isom] * vertex_coords[3*isom];
        test[1] += distBarCoords[isom] * vertex_coords[3*isom + 1];
        test[2] += distBarCoords[isom] * vertex_coords[3*isom + 2];

      }

      bftc_printf("point distant | verification \n");

      double dd = 0;
      for (int i = 0; i < 3; i++) {
        bftc_printf("  %f       |    %f \n",point_coords[i],test[i]);
        dd += (point_coords[i] - test[i]) * (point_coords[i] - test[i]);
      }

      if (sqrt(dd) > 1e-3)
        bftc_printf(" !!!! Erreur sur les coordonnees baryc directionf: %12.5e %i !!!!\n",sqrt(dd), isOnFace);
      else
        bftc_printf(" ++++ ok                                         : %12.5e %i ++++\n",sqrt(dd), isOnFace);

      bftc_printf("coord :");
      bftc_printf(" %12.5e %12.5e %12.5e", point_coords[0],
                  point_coords[1],
                  point_coords[2] );
      bftc_printf("\n");

      bftc_printf("coo b :");
      for (int isom = 0; isom < n_poly_vertex; isom++)
        bftc_printf(" %f", distBarCoords[isom]);

      bftc_printf("\n");

    }
  }

}

///
/// \brief Compute 3d shape functions and their derivatives given element
/// parametric coordinates.
///
/// This function is adapted from the CGNS interpolation tool.
///
///   @param [in]    elt_type    Type of element
///   @param [in]    uvw[]       Parametric coordinates
///   @param [out]   shapef[]    Barycenter's coordinates
///   @param [out]   deriv [][]  Derivative of shape function
///

void
LocationToLocalMesh::compute_shapef_3d(const cwipi_element_t elt_type,
                                       const double          uvw[3],
                                       double                shapef[8],
                                       double                deriv[8][3])

{
  switch (elt_type) {

  case CWIPI_CELL_TETRA4:

    shapef[0] = 1. - uvw[0] - uvw[1] - uvw[2];
    shapef[1] =      uvw[0];
    shapef[2] =               uvw[1];
    shapef[3] =                        uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -1.0;
      deriv[0][1] = -1.0;
      deriv[0][2] = -1.0;
      deriv[1][0] =  1.0;
      deriv[1][1] =  0.0;
      deriv[1][2] =  0.0;
      deriv[2][0] =  0.0;
      deriv[2][1] =  1.0;
      deriv[2][2] =  0.0;
      deriv[3][0] =  0.0;
      deriv[3][1] =  0.0;
      deriv[3][2] =  1.0;
    }

    break;

  case CWIPI_CELL_HEXA8:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
    shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
    shapef[6] = uvw[0] * uvw[1] * uvw[2];
    shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] = -(1.0 - uvw[1]) * uvw[2];
      deriv[4][1] = -(1.0 - uvw[0]) * uvw[2];
      deriv[4][2] =  (1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[5][0] =  (1.0 - uvw[1]) * uvw[2];
      deriv[5][1] = -uvw[0] * uvw[2];
      deriv[5][2] =  uvw[0] * (1.0 - uvw[1]);
      deriv[6][0] =  uvw[1] * uvw[2];
      deriv[6][1] =  uvw[0] * uvw[2];
      deriv[6][2] =  uvw[0] * uvw[1];
      deriv[7][0] = -uvw[1] * uvw[2];
      deriv[7][1] =  (1.0 - uvw[0]) * uvw[2];
      deriv[7][2] =  (1.0 - uvw[0]) * uvw[1];
    }

    break;

  case CWIPI_CELL_PRISM6:

    shapef[0] = (1.0 - uvw[0] - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[2]);
    shapef[2] = uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0] - uvw[1]) * uvw[2];
    shapef[4] = uvw[0] * uvw[2];
    shapef[5] = uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0] - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[2]);
      deriv[1][1] =  0.0;
      deriv[1][2] = -uvw[0];
      deriv[2][0] =  0.0;
      deriv[2][1] =  (1.0 - uvw[2]);
      deriv[2][2] = -uvw[1];
      deriv[3][0] = -uvw[2];
      deriv[3][1] = -uvw[2];
      deriv[3][2] =  (1.0 - uvw[0] - uvw[1]);
      deriv[4][0] =  uvw[2];
      deriv[4][1] =  0.0;
      deriv[4][2] =  uvw[0];
      deriv[5][0] =  0.0;
      deriv[5][1] =  uvw[2];
      deriv[5][2] =  uvw[1];
    }

    break;

  case CWIPI_CELL_PYRAM5:

    if (0 == 1) {

      // Alternative de coordonnees barycentriques

      shapef[0] = 0.125 * (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
      shapef[1] = 0.125 * (1.0 + uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
      shapef[2] = 0.125 * (1.0 + uvw[0]) * (1.0 + uvw[1]) * (1.0 - uvw[2]);
      shapef[3] = 0.125 * (1.0 - uvw[0]) * (1.0 + uvw[1]) * (1.0 - uvw[2]);
      shapef[4] = 0.5   *                                   (1.0 + uvw[2]);

      if (deriv != NULL) {
        deriv[0][0] = - 0.125 * (1.0 - uvw[1]) * (1.0 - uvw[2]);
        deriv[0][1] = - 0.125 * (1.0 - uvw[0]) * (1.0 - uvw[2]);
        deriv[0][2] = - 0.125 * (1.0 - uvw[0]) * (1.0 - uvw[1]);

        deriv[1][0] =   0.125 * (1.0 - uvw[1]) * (1.0 - uvw[2]);
        deriv[1][1] = - 0.125 * (1.0 + uvw[0]) * (1.0 - uvw[2]);
        deriv[1][2] = - 0.125 * (1.0 + uvw[0]) * (1.0 - uvw[1]);

        deriv[2][0] =   0.125 * (1.0 + uvw[1]) * (1.0 - uvw[2]);
        deriv[2][1] =   0.125 * (1.0 + uvw[0]) * (1.0 - uvw[2]);
        deriv[2][2] = - 0.125 * (1.0 + uvw[0]) * (1.0 + uvw[1]);

        deriv[3][0] = - 0.125 * (1.0 + uvw[1]) * (1.0 - uvw[2]);
        deriv[3][1] =   0.125 * (1.0 - uvw[0]) * (1.0 - uvw[2]);
        deriv[3][2] = - 0.125 * (1.0 - uvw[0]) * (1.0 + uvw[1]);

        deriv[4][0] = 0;
        deriv[4][1] = 0;
        deriv[4][2] = 0.5;
      }
    }

    else {

      shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
      shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
      shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
      shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
      shapef[4] = uvw[2];

      if (deriv != NULL) {
        deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
        deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
        deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
        deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
        deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
        deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
        deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
        deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
        deriv[2][2] = -uvw[0] * uvw[1];
        deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
        deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
        deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
        deriv[4][0] =  0.0;
        deriv[4][1] =  0.0;
        deriv[4][2] =  1.0;
      }
    }

    break;

  default:
    bftc_error(__FILE__, __LINE__, 0,
              "_compute_shapef: unhandled element type %s\n",
              fvmc_element_type_name[elt_type]);

  }

}

///
/// \brief Compute tetrzhedron, hexahedron, pyramid, or prism parametric coordinates
/// for a given point.
///
/// This function is adapted from the CGNS interpolation tool.
///
///   @param [in]    elt_type        Type of element
///   @param [in]    point_coords    Point coordinates
///   @param [in]    vertex_coords[] Pointer to element vertex coordinates
///   @param [in]    tolerance       Location tolerance factor
///   @param [out]   uvw[]           Parametric coordinates of point in element
///
///   @return                        Return 1 if uvw are computed, 0 otherwise
///

int
LocationToLocalMesh::compute_uvw(const cwipi_element_t elt_type,
                                 const double          point_coords[],
                                 const double          vertex_coords[8][3],
                                 const double          tolerance,
                                 double                uvw[3])

{
  int i, j, n_elt_vertices, iter;
  int max_iter = 20;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  switch (elt_type) {

  case CWIPI_CELL_TETRA4:
    n_elt_vertices = 4;
    break;

  case CWIPI_CELL_HEXA8:
    n_elt_vertices = 8;
    break;

  case CWIPI_CELL_PRISM6:
    n_elt_vertices = 6;
    break;

  case CWIPI_CELL_PYRAM5:
    n_elt_vertices = 5;
    break;

  default:
    break;

  }

  if (elt_type == CWIPI_CELL_TETRA4) {

    double vol6;

    double t00, t10, t20, t01, t02, t03, t11, t12, t13, t21, t22, t23;
    double v01[3], v02[3], v03[3];
    double v12[3], v13[3];
    double v23[3];
    double n_v[6];

    for (i = 0; i < 3; i++) {
      v01[i] = vertex_coords[1][i] - vertex_coords[0][i];
      v02[i] = vertex_coords[2][i] - vertex_coords[0][i];
      v03[i] = vertex_coords[3][i] - vertex_coords[0][i];
      v12[i] = vertex_coords[2][i] - vertex_coords[1][i];
      v13[i] = vertex_coords[3][i] - vertex_coords[1][i];
      v23[i] = vertex_coords[3][i] - vertex_coords[2][i];
    }

    n_v[0] = norm(v01);
    n_v[1] = norm(v02);
    n_v[2] = norm(v03);
    n_v[3] = norm(v12);
    n_v[4] = norm(v13);
    n_v[5] = norm(v23);

    double characteristicLength = std::min(n_v[0], n_v[1]);
    for (i = 2; i < 6; i++) {
      characteristicLength = std::min(characteristicLength, n_v[i]);
    }

    vol6 = fabs(  v01[0] * (v02[1]*v03[2] - v02[2]*v03[1])
                  - v02[0] * (v01[1]*v03[2] - v01[2]*v03[1])
                  + v03[0] * (v01[1]*v02[2] - v01[2]*v02[1]));

    double epsilon_denom  = geometricEpsilon(characteristicLength,
                                             GEOM_EPS_VOL);

    if (vol6 < epsilon_denom){
      bftc_error(__FILE__, __LINE__, 0,
                 "compute_uvw : degenerated tetra, volume = %12.5e\n", vol6);
    }

    t00  =   point_coords[0] - vertex_coords[0][0];
    t10  =   point_coords[1] - vertex_coords[0][1];
    t20  =   point_coords[2] - vertex_coords[0][2];

    t01  = - vertex_coords[0][0] + vertex_coords[1][0];
    t02  = - vertex_coords[0][0] + vertex_coords[2][0];
    t03  = - vertex_coords[0][0] + vertex_coords[3][0];

    t11  = - vertex_coords[0][1] + vertex_coords[1][1];
    t12  = - vertex_coords[0][1] + vertex_coords[2][1];
    t13  = - vertex_coords[0][1] + vertex_coords[3][1];

    t21  = - vertex_coords[0][2] + vertex_coords[1][2];
    t22  = - vertex_coords[0][2] + vertex_coords[2][2];
    t23  = - vertex_coords[0][2] + vertex_coords[3][2];

    uvw[0] = (  t00 * (t12*t23 - t13*t22)
              - t10 * (t02*t23 - t22*t03)
              + t20 * (t02*t13 - t12*t03)) / vol6;
    uvw[1] = (- t00 * (t11*t23 - t13*t21)
              + t10 * (t01*t23 - t21*t03)
              - t20 * (t01*t13 - t03*t11)) / vol6;
    uvw[2] = (  t00 * (t11*t22 - t21*t12)
              - t10 * (t01*t22 - t21*t02)
              + t20 * (t01*t12 - t11*t02)) / vol6;
  }

  else if (   elt_type == CWIPI_CELL_HEXA8
              || elt_type == CWIPI_CELL_PRISM6
              || elt_type == CWIPI_CELL_PYRAM5) {

    /* Use Newton-method to determine parametric coordinates and shape function */

    for (i = 0; i < 3; i++)
      uvw[i] = 0.5;

    for (iter = 0; iter < max_iter; iter++) {

      compute_shapef_3d(elt_type, uvw, shapef, dw);

      b[0] = - point_coords[0];
      b[1] = - point_coords[1];
      b[2] = - point_coords[2];

      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
          a[i][j] = 0.0;
      }

      for (i = 0; i < n_elt_vertices; i++) {

        b[0] += (shapef[i] * vertex_coords[i][0]);
        b[1] += (shapef[i] * vertex_coords[i][1]);
        b[2] += (shapef[i] * vertex_coords[i][2]);

        for (j = 0; j < 3; j++) {
          a[0][j] -= (dw[i][j] * vertex_coords[i][0]);
          a[1][j] -= (dw[i][j] * vertex_coords[i][1]);
          a[2][j] -= (dw[i][j] * vertex_coords[i][2]);
        }

      }

      if (inverse_3x3(a, b, x)) {
        bftc_printf("error compute_uvw : matrice non iversible\n");
        return 0;
      }

      dist = 0.0;

      for (i = 0; i < 3; i++) {
        dist += x[i] * x[i];
        uvw[i] += x[i];
      }

      if (dist <= (tolerance * tolerance))
        return 1;

    }

    bftc_printf("Warning compute_uvw : no convergence for the point (%12.5e, %12.5e, %12.5e), dist = %12.5e\n",
                point_coords[0], point_coords[1], point_coords[2], dist);
    bftc_printf("                      in the element : (%12.5e, %12.5e, %12.5e)\n", vertex_coords[0][0],
                vertex_coords[0][1],
                vertex_coords[0][2]);
    for (i = 1; i < n_elt_vertices; i++) {
      bftc_printf("                                       (%12.5e, %12.5e, %12.5e)\n", vertex_coords[i][0],
                  vertex_coords[i][1],
                  vertex_coords[i][2]);
    }
    bftc_printf("                            shapef :");
    for (i = 1; i < n_elt_vertices; i++) {
      bftc_printf(" %12.5e", shapef[i]);
    }
    bftc_printf("\n");

    return 1;

  }

  else {
    bftc_error(__FILE__, __LINE__, 0,
               "compute_uvw : unhandled element type\n");
  }

  return 1;
}


} // Namespace CWIPI

#undef _DOT_PRODUCT
#undef _MODULE
#undef _CROSS_PRODUCT
