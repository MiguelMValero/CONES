#ifndef LOCATIONTODISTANTMESH_HXX_
#define LOCATIONTODISTANTMESH_HXX_
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

#include <string>
#include <vector>

#include "cwipi.h"
#include "applicationProperties.hxx"

namespace cwipi
{

///
/// \brief This Class realize the location of point in a distant mesh
///
/// This class must be used with symetric class LocationToLocalMesh
///

class LocationToLocalMesh;

class LocationToDistantMesh
{

  friend class LocationToLocalMesh;

public:

  ///
  /// \brief Constructor
  ///
  ///   @param isCoupledrank  indicate if current MPI rank is a coupled rank
  ///

  LocationToDistantMesh(const bool isCoupledrank,
                        const cwipi_coupling_type_t _couplingType,
                        const ApplicationProperties& localApplicationProperties);

  ///
  /// \brief Default destructor
  ///

  virtual ~LocationToDistantMesh();

  ///
  /// \brief Set the type of information to be exchanged at the location
  ///
  ///   @param information
  ///

  inline void setInfo(cwipi_located_point_info_t info);

  ///
  /// \brief Set points to locate
  ///
  ///   @param nPointsToLocate      Number of points to locate
  ///   @param coordsPointsToLocate Coordinates of points to locate
  ///

  void setpointsToLocate(int nPointsToLocate, double *coordsPointsToLocate);

  ///
  /// 
  ///
  size_t  locationSize();

  void packLocation(unsigned char *buff);
  void unpackLocation(unsigned char *buff);

  ///
  /// \brief Points location
  ///
  void locate();

  ///
  /// \brief Get the number of points to locate
  ///

  inline int getNpointsToLocate() const;

  ///
  /// \brief Get the coordinates of points to locate
  ///

  inline const double & getCoordsPointsToLocate() const;
  inline const double * getCoordsPointsToLocate();

  ///
  /// \brief Get the number of unlocated points
  ///

  inline int getNUnlocatedPoint() const;

  ///
  /// \brief Get list of unlocated points
  ///

  inline const int *getUnlocatedPoint() const;

  ///
  /// \brief Get the number of located points
  ///

  inline int getNLocatedPoint() const;

  ///
  /// \brief Get list of located points
  ///

  inline const int *getLocatedPoint() const;

  ///
  /// \brief Get list of number of vertices of containing element
  ///

  inline const int *getElementContainingNVertex() const;

  ///
  /// \brief Get connectivity of element that contains each located point
  ///

  inline const int *getElementContainingVertex() const;

  ///
  /// \brief Get vertices coordinates of the element that contains each located point
  ///

  inline const double *getElementContainingVertexCoords() const;

  ///
  /// \brief Get barycentric coordinates for the element that contains each located point
  ///

  inline const double *getElementContainingBarycentricCoordinates() const;

  ///
  /// \brief Get MPI rank that contains each located point
  ///

  inline const int *getElementContainingMPIrank() const;

  ///
  /// \brief ///< Distant element that contains each located points (size : nLocatedpoint)
  ///

  inline const int *getElementContaining() const;

  ///
  /// \brief Get "tolocate" status
  ///

  inline bool getToLocateStatus() const;

  ///
  /// \brief Set information to exchange
  ///

  void synchronize();

  ///
  /// \brief Clear location
  ///

  void clear();

private :
  LocationToDistantMesh();
  LocationToDistantMesh(const LocationToDistantMesh &other);
  LocationToDistantMesh &operator=(const LocationToDistantMesh &other);
private :
  const bool _isCoupledRank;                       ///< Status indicating whether the current MPI rank is coupled
  const cwipi_coupling_type_t _couplingType;       ///< Coupling type
  const ApplicationProperties& _localApplicationProperties; ///< Local application properties
  int     _nPointsToLocate;                          ///< Number of points to locate
  double *_coordsPointsToLocate;                    ///< Coordinates of points to locate
  bool    _toLocate;                                 ///< Status to generate location

  int     _nUnlocatedPoint;                         ///< Number of unlocated points
  int     _nLocatedPoint;                           ///< Number of located points
  int    *_unlocatedPoint;                          ///< Unlocated points
  int    *_locatedPoint;                            ///< Located points
  int    *_elementContaining;                       ///< Distant element that contains each located points (size : nLocatedpoint)
  int    *_elementContainingNVertex;                ///< Number of vertex of element containing (size : nLocatedpoint + 1)
  int    *_elementContainingVertex;                 ///< Vertices of element containing (size : _nElementContainingVertex[nLocatedpoint])
  double *_elementContainingVertexCoords;           ///< Vertex coordinates of element containing (size : 3*_nElementContainingVertex[nLocatedpoint])
  double *_elementContainingBarycentricCoordinates; ///< Barycentric coordinates (size : _nElementContainingVertex[nLocatedpoint])
  int    *_elementContainingMPIrankContaining;      ///< MPI rank containing (size : nLocatedpoint)
  cwipi_located_point_info_t _locationInfo;         ///< received information for located points
};

///
/// \brief Get the number of points to locate
///

int LocationToDistantMesh::getNpointsToLocate() const
{
  return _nPointsToLocate;
}

///
/// \brief Get the coordinates of points to locate
///

const double & LocationToDistantMesh::getCoordsPointsToLocate() const
{
  return *_coordsPointsToLocate;
}

const double * LocationToDistantMesh::getCoordsPointsToLocate()
{
  return _coordsPointsToLocate;
}

///
/// \brief Get the number of unlocated points
///

int LocationToDistantMesh::getNUnlocatedPoint() const
{
  return _nUnlocatedPoint;
}

///
/// \brief Get list of unlocated points
///

const int *LocationToDistantMesh::getUnlocatedPoint() const
{
  return _unlocatedPoint;
}

///
/// \brief Get the number of located points
///

int LocationToDistantMesh::getNLocatedPoint() const
{
  return _nPointsToLocate - _nUnlocatedPoint;
}

///
/// \brief Get list of located points
///

const int *LocationToDistantMesh::getLocatedPoint() const
{
  return _locatedPoint;
}

///
/// \brief Get list of number of vertices of containing element
///

const int *LocationToDistantMesh::getElementContainingNVertex() const
{
  if (_locationInfo != CWIPI_DISTANT_MESH_INFO)
    bftc_error(__FILE__, __LINE__, 0, "This field is unavailable : Activate with setInfo method" );
  return _elementContainingNVertex;
}

///
/// \brief Get connectivity of element that contains each located point
///

const int *LocationToDistantMesh::getElementContainingVertex() const
{
  if (_locationInfo != CWIPI_DISTANT_MESH_INFO)
    bftc_error(__FILE__, __LINE__, 0, "This field is unavailable : Activate with setInfo method" );
  return _elementContainingVertex;
}

///
/// \brief Get vertices coordinates of the element that contains each located point
///

const double *LocationToDistantMesh::getElementContainingVertexCoords() const
{
  if (_locationInfo != CWIPI_DISTANT_MESH_INFO)
    bftc_error(__FILE__, __LINE__, 0, "This field is unavailable : Activate with setInfo method" );
  return _elementContainingVertexCoords;
}

///
/// \brief Get barycentric coordinates for the element that contains each located point
///

const double *LocationToDistantMesh::getElementContainingBarycentricCoordinates() const
{
  if (_locationInfo != CWIPI_DISTANT_MESH_INFO)
    bftc_error(__FILE__, __LINE__, 0, "This field is unavailable : Activate with setInfo method" );
  return _elementContainingBarycentricCoordinates;
}

///
/// \brief Get MPI rank that contains each located point
///

const int *LocationToDistantMesh::getElementContainingMPIrank() const
{
  if (_locationInfo != CWIPI_DISTANT_MESH_INFO)
    bftc_error(__FILE__, __LINE__, 0, "This field is unavailable : Activate with setInfo method" );
  return _elementContainingMPIrankContaining;
}

///
/// \brief Get "tolocate" status
///

bool LocationToDistantMesh::getToLocateStatus() const
{
  return _toLocate;
}

///
/// \brief Set information type to receive with location result
///
///   @param information type to receive with location result
///

void LocationToDistantMesh::setInfo(cwipi_located_point_info_t info)
{
  _locationInfo = info;
}

///
/// \brief ///< Distant element that contains each located points (size : nLocatedpoint)
///

inline const int *LocationToDistantMesh::getElementContaining() const
{
  if (_locationInfo != CWIPI_DISTANT_MESH_INFO)
    bftc_error(__FILE__, __LINE__, 0, "This field is unavailable : Activate with setInfo method" );
  return _elementContaining;
}

}

#endif /* LOCATIONTODISTANTMESH_HXX_ */
