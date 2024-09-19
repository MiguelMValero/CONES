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
/*
 * DistantLocation.cxx
 *
 *  Created on: Oct 16, 2009
 *      Author: equemera
 */


#include <mpi.h>
#include <cstring>
#include <fvmc_locator.h>

#include "locationToDistantMesh.hxx"

namespace cwipi
{

///
/// \brief Default constructor
///

LocationToDistantMesh::LocationToDistantMesh(const bool isCoupledRank,
                                             const cwipi_coupling_type_t couplingType,
                                             const ApplicationProperties& localApplicationProperties)
: _isCoupledRank(isCoupledRank), _couplingType(couplingType),
  _localApplicationProperties(localApplicationProperties)
{
  _coordsPointsToLocate = NULL;
  _nPointsToLocate = 0;
  _nUnlocatedPoint = 0 ;
  _nLocatedPoint = 0;
  _unlocatedPoint = NULL;
  _locatedPoint = NULL;
  _elementContaining = NULL;
  _elementContainingBarycentricCoordinates = NULL;
  _elementContainingMPIrankContaining = NULL;
  _elementContainingNVertex = NULL;
  _elementContainingVertex = NULL;
  _elementContainingVertexCoords = NULL;
  _locationInfo = CWIPI_BASIC_INFO;
  _toLocate = true;

}

///
/// \brief Default destructor
///

LocationToDistantMesh::~LocationToDistantMesh()
{
  clear();
}

///
/// \brief Set points to locate
///
///   @param nPointsToLocate      Number of points to locate
///   @param coordsPointsToLocate Coordinates of points to locate
///

void LocationToDistantMesh::setpointsToLocate(int nPointsToLocate, double *coordsPointsToLocate)
{
  _toLocate = true;
  if (_isCoupledRank) {
    _nPointsToLocate = nPointsToLocate;
    _coordsPointsToLocate = coordsPointsToLocate;
  }
}

///
/// \brief Set points to locate
///
///   @param nPointsToLocate      Number of points to locate
///   @param coordsPointsToLocate Coordinates of points to locate
///

void LocationToDistantMesh::synchronize()
{
  const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

  int rootRank;

  if (_isCoupledRank) {
    if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      MPI_Comm_rank(localComm, &rootRank);

      assert(rootRank == 0);

      MPI_Bcast(&_nPointsToLocate, 1, MPI_INT, rootRank, localComm );
      MPI_Bcast(&_nUnlocatedPoint, 1, MPI_INT, rootRank, localComm );
      MPI_Bcast(&_nLocatedPoint, 1, MPI_INT, rootRank, localComm );

      MPI_Bcast(_locatedPoint,
                _nLocatedPoint,
                MPI_INT,
                rootRank,
                localComm );

      MPI_Bcast(_unlocatedPoint, _nUnlocatedPoint, MPI_INT, rootRank, localComm );

      _toLocate = false;

      if (_locationInfo == CWIPI_DISTANT_MESH_INFO) {


      }
    }
  }

  else {
    if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      rootRank = 0;

      _toLocate = false;

      MPI_Bcast(&_nPointsToLocate, 1, MPI_INT, rootRank, localComm );

      MPI_Bcast(&_nUnlocatedPoint, 1, MPI_INT, rootRank, localComm );
      MPI_Bcast(&_nLocatedPoint, 1, MPI_INT, rootRank, localComm );

      if (_locatedPoint != NULL)
        free ( _locatedPoint);

      if (_unlocatedPoint != NULL)
        free ( _unlocatedPoint );

      _locatedPoint =  (int *) malloc (sizeof(int) * (_nLocatedPoint));
      _unlocatedPoint =  (int *) malloc (sizeof(int) * (_nUnlocatedPoint));

      MPI_Bcast(_locatedPoint,
                _nLocatedPoint,
                MPI_INT,
                rootRank,
                localComm);
      MPI_Bcast( _unlocatedPoint, _nUnlocatedPoint, MPI_INT, rootRank, localComm );

      if (_locationInfo == CWIPI_DISTANT_MESH_INFO) {

        if (_elementContaining != NULL)
          free ( _elementContaining);

        if (_elementContainingNVertex != NULL)
          free ( _elementContainingNVertex );

        if (_elementContainingBarycentricCoordinates != NULL)
          free ( _elementContainingBarycentricCoordinates );

        if (_elementContainingVertex != NULL)
          free ( _elementContainingVertex );

        if (_elementContainingVertexCoords != NULL)
          free ( _elementContainingVertexCoords );

        _elementContaining =  (int *) malloc (sizeof(int) * (_nLocatedPoint));
        _elementContainingNVertex =  (int *) malloc (sizeof(int) * (_nLocatedPoint+1));

        MPI_Bcast(_elementContaining,
                  _nLocatedPoint,
                  MPI_INT,
                  rootRank,
                  localComm);

        MPI_Bcast(_elementContainingNVertex,
                  _nLocatedPoint + 1,
                  MPI_INT,
                  rootRank,
                  localComm);

        _elementContainingBarycentricCoordinates =  (double *) malloc (sizeof(double) * (_elementContainingNVertex[_nLocatedPoint]));
        _elementContainingVertex =  (int *) malloc (sizeof(int) * (_elementContainingNVertex[_nLocatedPoint]));
        _elementContainingVertexCoords =  (double *) malloc (sizeof(double) * (3 * _elementContainingNVertex[_nLocatedPoint]));

        MPI_Bcast(_elementContainingBarycentricCoordinates,
                  _elementContainingNVertex[_nLocatedPoint],
                  MPI_DOUBLE,
                  rootRank,
                  localComm);

        MPI_Bcast(_elementContainingVertex,
                  _elementContainingNVertex[_nLocatedPoint],
                  MPI_INT,
                  rootRank,
                  localComm);

        MPI_Bcast(_elementContainingVertexCoords,
                  3 *_elementContainingNVertex[_nLocatedPoint],
                  MPI_DOUBLE,
                  rootRank,
                  localComm);

      }
    }
  }
}

///
/// 
///
size_t LocationToDistantMesh::locationSize()
{
  size_t il_size = 0;
  il_size += 3 * sizeof(int);

  il_size += sizeof(int);
  if (_elementContainingNVertex != NULL) 
    il_size += (_nLocatedPoint+1) * sizeof(int);

  il_size += sizeof(int);
  if (_elementContainingVertex != NULL) 
    il_size += _elementContainingNVertex[_nLocatedPoint] * sizeof(int);

  il_size += sizeof(int);
  if (_elementContainingBarycentricCoordinates != NULL) 
    il_size += _elementContainingNVertex[_nLocatedPoint] * sizeof(double);
  
  il_size += sizeof(int);
  if(_elementContainingMPIrankContaining != NULL) 
    il_size += _nLocatedPoint *sizeof(int);

  il_size += sizeof(int);
  if (_elementContainingVertexCoords != NULL) 
    il_size += (3 * _elementContainingNVertex[_nLocatedPoint]) * sizeof(double);

  il_size += sizeof(int);
  if (_elementContaining != NULL) 
    il_size += _nLocatedPoint * sizeof(int);
  
  return il_size;
  
}

void LocationToDistantMesh::packLocation(unsigned char *buff)
{
  int s;
  void *p;
  size_t s_pack;

  p = (void *)buff;

  s_pack = sizeof(int);
  memcpy(p,(void *)&_nPointsToLocate, s_pack);
  p = (void *) ((char *) p + s_pack);

  s_pack = sizeof(int);
  memcpy(p,(void *)&_nLocatedPoint, s_pack);
  p = (void *) ((char *) p + s_pack);

  s_pack = sizeof(int);
  memcpy(p,(void *)&_nUnlocatedPoint, s_pack);
  p = (void *) ((char *) p + s_pack);
  
  // pour chaque tableau, on commence par stoker sa taille 
  // pour pouvoir l'allouer si nécessaire à la lecture
  
  if (_elementContainingNVertex != NULL) {
    s = _nLocatedPoint+1;

    s_pack = sizeof(int);
    memcpy(p,(void *)&s, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = s*sizeof(int);
    memcpy(p,(void *) _elementContainingNVertex, s_pack);
    p = (void *) ((char *) p + s_pack);
  } 
  else {
    s = 0;

    s_pack = sizeof(int);
    memcpy(p,(void *)&s,  s_pack);
    p = (void *) ((char *) p + s_pack);
  }

  if (_elementContainingVertex != NULL) {
    s = _elementContainingNVertex[_nLocatedPoint];

    s_pack = sizeof(int);
    memcpy(p,(void *)&s,  s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = s*sizeof(int);
    memcpy(p,(void *) _elementContainingVertex, s_pack);
    p = (void *) ((char *) p + s_pack);
  } 
  else {
    s = 0;
    s_pack = sizeof(int);
    memcpy(p,(void *)&s,  s_pack);
    p = (void *) ((char *) p + s_pack);
  }

  if (_elementContainingBarycentricCoordinates != NULL) {
    s = _elementContainingNVertex[_nLocatedPoint];

    s_pack = sizeof(int);
    memcpy(p,(void *)&s,  s_pack);
    p = (void *) ((char *) p + s_pack);
    
    s_pack = s*sizeof(double);
    memcpy(p,(void *) _elementContainingBarycentricCoordinates, s_pack);
    p = (void *) ((char *) p + s_pack);
  } 
  else {
    s = 0;

    s_pack = sizeof(int);
    memcpy(p,(void *)&s,  s_pack);
    p = (void *) ((char *) p + s_pack);
  }

  // attention : alloué dans locationToLocalMesh
  if(_elementContainingMPIrankContaining != NULL) {
    s = _nLocatedPoint;

    s_pack = sizeof(int);
    memcpy(p,(void *)&s,  s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = s*sizeof(int);
    memcpy(p,(void *) _elementContainingMPIrankContaining, s_pack);
    p = (void *) ((char *) p + s_pack);
  } 
  else {
    s = 0;

    s_pack = sizeof(int);
    memcpy(p,(void *)&s, s_pack);
    p = (void *) ((char *) p + s_pack);
  }
  
  if (_elementContainingVertexCoords != NULL) {
    s = (3 * _elementContainingNVertex[_nLocatedPoint]);

    s_pack = sizeof(int);
    memcpy(p,(void *)&s, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = s*sizeof(double);
    memcpy(p,(void *)_elementContainingVertexCoords, s_pack);
    p = (void *) ((char *) p + s_pack);
  } 
  else {
    s = 0;

    s_pack = sizeof(int);
    memcpy(p,(void *)&s, s_pack);
    p = (void *) ((char *) p + s_pack);
  }

  if (_elementContaining != NULL) {
    s = _nLocatedPoint;

    s_pack = sizeof(int);
    memcpy(p,(void *)&s, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = _nLocatedPoint*sizeof(int);
    memcpy(p,(void *)_elementContaining, s_pack);
    p = (void *) ((char *) p + s_pack);
  } 
  else {
    s = 0;

    s_pack = sizeof(int);
    memcpy(p,(void *)&s, s_pack);
    p = (void *) ((char *) p + s_pack);
  }
}

void LocationToDistantMesh::unpackLocation(unsigned char *buff)
{
  int s;
  size_t cur_pos;
  
  cur_pos = 0;

  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_nPointsToLocate,sizeof(int));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_nLocatedPoint,sizeof(int));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&_nUnlocatedPoint,sizeof(int));
 
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContainingNVertex != NULL) free ( _elementContainingNVertex);
    _elementContainingNVertex =  (int *) malloc (sizeof(int) * (s));
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)_elementContainingNVertex,s*sizeof(int));
  }
    
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContainingVertex != NULL) free ( _elementContainingVertex);
    _elementContainingVertex =  (int *) malloc (sizeof(int) * (s));
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *) _elementContainingVertex,s*sizeof(int));
  }

  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContainingBarycentricCoordinates != NULL) free ( _elementContainingBarycentricCoordinates);
    _elementContainingBarycentricCoordinates =  (double *) malloc (sizeof(double) * (s));
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *) _elementContainingBarycentricCoordinates,s*sizeof(double));
  } 

  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if(_elementContainingMPIrankContaining != NULL) free ( _elementContainingMPIrankContaining);
    _elementContainingMPIrankContaining =  (int *) malloc (sizeof(int) * (s));
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *) _elementContainingMPIrankContaining,s*sizeof(int));
  } 
 
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContainingVertexCoords != NULL) free ( _elementContainingVertexCoords);
    _elementContainingVertexCoords = (double *) malloc (sizeof(double) * s) ;
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)_elementContainingVertexCoords,s*sizeof(double));
  } 
  
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)&s, sizeof(int));
  if (s != 0) {
    if (_elementContaining != NULL) free ( _elementContaining);
    _elementContaining =  (int *) malloc (sizeof(int) * (s));
    cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos],(void *)_elementContaining,_nLocatedPoint*sizeof(int));
  } 
  _toLocate = false;
}




///
/// \brief Clear location
///

void LocationToDistantMesh::clear()
{
  if (!_isCoupledRank && _couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {
    if (_locatedPoint != NULL) {
      free (  _locatedPoint);
      _locatedPoint = NULL;
    }
    

    if (_unlocatedPoint != NULL) {
      free (  _unlocatedPoint);
      _unlocatedPoint = NULL;
    }

  }

  if (_elementContainingBarycentricCoordinates != NULL)
    free ( _elementContainingBarycentricCoordinates);
  if (_elementContainingMPIrankContaining != NULL)
    free ( _elementContainingMPIrankContaining);
  if (_elementContainingNVertex != NULL)
    free ( _elementContainingNVertex);
  if (_elementContainingVertex != NULL)
    free ( _elementContainingVertex);
  if (_elementContainingVertexCoords != NULL)
    free ( _elementContainingVertexCoords);
  if (_elementContaining != NULL)
    free ( _elementContaining);

  _elementContainingBarycentricCoordinates = NULL;
  _elementContainingMPIrankContaining = NULL;
  _elementContainingNVertex = NULL;
  _elementContainingVertex = NULL;
  _elementContainingVertexCoords = NULL;
  _elementContaining = NULL;

}

}
