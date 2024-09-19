#ifndef LOCATIONTOLOCALMESH_HXX_
#define LOCATIONTOLOCALMESH_HXX_
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
#include <cmath>

#include <fvmc_locator.h>
#include <fvmc_nodal.h>
#include <bftc_error.h>
#include <fvmc_triangulate.h>

#include "oldMesh.hxx"
#include "cwipi.h"

namespace cwipi
{

class LocationToDistantMesh;
class ApplicationProperties;

enum {
  X,
  Y,
  Z
} ;

class LocationToLocalMesh
{
public:

  ///
  /// \brief Constructor
  ///
  ///   @param supportMesh                                 Mesh where distant points are localized
  ///   @param solverType                                  Solver type of current application
  ///   @param optBboxStep                                 Discritization step for high order element bounding boxes
  ///   @param tolerance                                   Geometric tolerance for localization
  ///   @param couplingComm                                Coupling MPI communicator
  ///   @param coupledApplicationNRankCouplingComm         Rank number of coupled application
  ///   @param coupledApplicationBeginningRankCouplingComm Beginning rank of coupled application
  ///   @param isCoupledRank                               Indicate if current MPI rank is a coupled rank
  ///   @param entitiesDim                                 Geometric entities dimension
  ///   @param localApplicationProperties                  Local application properties
  ///   @param locationToDistantMesh                       Information about local points location in the distant mesh
  ///

  LocationToLocalMesh(
                      const cwipi_solver_type_t  &solverType,
                      const int optBboxStep,
                      const double &tolerance,
                      const MPI_Comm& couplingComm,
                      const int &coupledApplicationNRankCouplingComm,
                      const int &coupledApplicationBeginningRankCouplingComm,
                      const bool isCoupledRank,
                      const int entitiesDim,
                      const ApplicationProperties& localApplicationProperties,
                      LocationToDistantMesh &locationToDistantMesh);

  ///
  /// \brief Default destructor
  ///

  virtual ~LocationToLocalMesh();

  void packLocation(unsigned char *buff);
  void unpackLocation(unsigned char *buff);
  size_t locationSize();


  ///
  /// \brief distant points location in the local mesh
  ///

  void locate();

  ///
  /// \brief Get barycentric coordinates index  of distant points in the local mesh (size = nDistantPoint + 1)
  ///

  //inline const std::vector <int> & getBarycentricCoordinatesIndex() const;
  inline const int *getBarycentricCoordinatesIndex() const;

  ///
  /// \brief Get barycentric coordinates of distant points in the local mesh (size = BarycentricCoordinatesIndex[nDistantPoint])
  ///

  //inline const std::vector <double> & getBarycentricCoordinates() const;
  inline const double *getBarycentricCoordinates() const;

  ///
  /// \brief Get uvw (size = 3 * nDistantPoint])
  ///

  //inline const std::vector <double> & getBarycentricCoordinates() const;
  inline const double *getUvw() const;

  ///
  /// \brief Return location result (size = nDistantpoint)
  ///

  inline const int *getLocation() const;

  ///
  /// \brief Return distance to location element result (size = nDistantpoint)
  ///

  inline const float *getDistance() const;

  ///
  /// \brief Return number of located distant point
  ///

  inline int getNLocatedDistantPoint() const;

  ///
  /// \brief Return number of distant ranks
  ///

  inline int getNDistantRank() const;

  ///
  /// \brief Return located points distribution on distant ranks (size = n_distant_rank + 1)
  ///

  inline const int *getLocatedPointsDistribution() const;

  ///
  /// \brief Return distant point distribution on distant ranks (size = n_distant_rank + 1)
  ///

  inline const int *getDistantDistribution() const;

  ///
  /// \brief Return number coordinates of located distant point
  ///

  inline const double *getPointCoordinates() const;

  ///
  /// \brief Return fvm locator
  ///

  inline fvmc_locator_t *getFVMLocator() const;

  ///
  /// \brief Exchange field on vertices of cells that contain each located points
  ///
  ///   @param [in]  sendingField    Field defined on local mesh vertices
  ///   @param [out] receivingField  Field defined on vertices of distant
  ///                                elements that contain each located point
  ///   @param [in]  stride          Number of field component
  ///

  void exchangeCellVertexFieldOfElementContaining(double *sendingField, 
                                                  double *receivingField, 
                                                  const int stride);

  ///
  /// \brief Exchange field on cells that contain each located points
  ///
  ///   @param [in]  sendingField    Field defined on local mesh vertices
  ///   @param [out] receivingField  Field defined on vertices of distant
  ///                                elements that contain each located point
  ///   @param [in]  stride          Number of field component
  ///

  void exchangeCellCenterFieldOfElementContaining(double *sendingField,
                                                  double *receivingField, 
                                                  const int stride);

  ///
  /// \brief Set support mesh
  ///
  ///   @param [in]      supportMesh location support mesh
  ///

  inline void setSupportMesh(oldMesh *supportMesh, bool lb_tolocate);


  ///
  /// \brief Set the step for the computaion of the ho optimized element bounding box 
  ///
  ///   @param [in]  step
  ///

  inline void optBboxStep(const int step);

  
private :


  void  _computeBary (int numPts, double *pts, double bary[3]);

  
  void _project_point2(double x[3], double pt_plan[3],
                       double normal[3], double xproj[3]);

  void _computeNormal (int numPts, double *pts, double n[3]);

  ///
  /// \brief Compute Mean Values 2D
  ///
  ///

  void compute2DMeanValues();

///
/// \brief Compute Polygo Mean Values
///

  void computePolygonMeanValues(const int n_dist_points,
                                const fvmc_lnum_t *dist_locations,
                                const fvmc_coord_t *dist_coords,
                                const int *meshConnectivityIndex,
                                const int *meshConnectivity,
                                const double *meshVertexCoords,
                                const std::vector <int>& nDistBarCoords,
                                std::vector <double>& distBarCoords);
  
  ///
  /// \brief Compute Mean Values 3D
  ///
  ///

  void compute3DMeanValues();

  ///
  /// \brief Compute Polyhedra Mean Values 3D
  ///
  ///

  void compute3DMeanValuesPoly(const double point_coords[],
                               const int    n_poly_faces,
                               const int    n_poly_vertex,
                               const int    faceDirection[],
                               const int    faceToVertexIdx[],
                               const int    faceToVertex[],
                               const double vertex_coords[],
                               const double characteristicLength,
                               const float  distElt,
                               double       distBarCoords[]);


  ///
  /// \brief Projection to the midplane
  ///
  ///   @param [in]         nbr_som_fac location support mesh
  ///   @param [inout]      coo_som_fac coordinates of face vertices
  ///   @param [inout]      coo_point_dist coordinates of distant points
  ///

  void midplaneProjection
  (
   const int     nbr_som_fac,
   double *const coo_som_fac,
   double *const coo_point_dist
   );

  ///
  /// \brief Compute tetrzhedron, hexahedron, pyramid, or prism parametric coordinates for a given point.
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
  compute_uvw(const cwipi_element_t elt_type,
              const double          point_coords[],
              const double          vertex_coords[8][3],
              const double          tolerance,
              double                uvw[3]);

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
  compute_shapef_3d(const cwipi_element_t elt_type,
                    const double          uvw[3],
                    double                shapef[8],
                    double                deriv[8][3]);

private :

  oldMesh                       *_supportMesh;                                 ///< Mesh where distant points are localized
  const cwipi_solver_type_t  &_solverType;                                  ///< Solver type of current application
  int                      _optBboxStep;                                 ///< Discitization step for the computation
                                                                            //   of the optimized high order element bounding box 
  const double               &_tolerance;                                   ///< Geometric tolerance for localization
  const MPI_Comm             &_couplingComm;                                ///< Coupling MPI communicator
  const int                  &_coupledApplicationNRankCouplingComm;         ///< Rank number of coupled application
  const int                  &_coupledApplicationBeginningRankCouplingComm; ///< Beginning rank of coupled application
  const bool                  _isCoupledRank;                               ///< Indicate if current MPI rank is a coupled rank
  const int                   _entitiesDim;                                 ///< Geometric entities dimension
  const ApplicationProperties&_localApplicationProperties;                  ///< Application properties
  LocationToDistantMesh      &_locationToDistantMesh;                       ///< Information about local points location in the distant mesh

  fvmc_locator_t             *_fvmLocator;                                  ///< fvm structure that build the location
  std::vector <int>          *_barycentricCoordinatesIndex;                 ///< Barycentric coordinates for each
  std::vector <double>       *_barycentricCoordinates;                      ///< Barycentric coordinates associated to the eleme
  std::vector <double>        *_uvw;                                         ///< uvw coordinates in the element (only for ho elements)
  int                         _nDistantPoint;                               ///< Number of distant points located in the local mesh
  int                        *_location;                                    ///< Local elements that contain distant points
  int                        *_locatedPointsDistribution;                   ///< Located points distribution on distant ranks
  int                        *_distantDistribution;                         ///< Distant point distribution on distant ranks
  float                      *_distance;                                    ///< distance to Local elements that contain distant points
  std::vector <int>          *_nVertex;                                     ///< Vertices number of local elements that contain distant points
  bool                        _toLocate;                                    ///< Status to activate location
  int                         _maxElementContainingNVertex;                 ///< Maximum number of vertices of elements that contain located point
};

///
/// \brief Get barycentric coordinates index  of distant points in the local mesh (size = nDistantPoint + 1)
///

//   const std::vector <int> & LocationToLocalMesh::getBarycentricCoordinatesIndex() const
//   {
//     if (_toLocate)
//       bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
//     return *_barycentricCoordinatesIndex;
//   }

  const int* LocationToLocalMesh::getBarycentricCoordinatesIndex() const
  {
    if (_toLocate)
      bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    if (_barycentricCoordinatesIndex == NULL)
      return NULL;
    else 
      return &(*_barycentricCoordinatesIndex)[0];
  }

///
/// \brief Get barycentric coordinates of distant points in the local mesh (size = BarycentricCoordinatesIndex[nDistantPoint])
///

//   const std::vector <double> & LocationToLocalMesh::getBarycentricCoordinates() const
//   {
//     if (_toLocate)
//       bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
//     return *_barycentricCoordinates;
//   }

  const double *LocationToLocalMesh::getBarycentricCoordinates() const
  {
    if (_toLocate)
      bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    if (_barycentricCoordinates == NULL)
      return NULL;
    else
      return &(*_barycentricCoordinates)[0];
  }

  const double *LocationToLocalMesh::getUvw() const
  {
    if (_toLocate)
      bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
    if (_uvw == NULL)
      return NULL;
    else
      return &(*_uvw)[0];
  }

///
/// \brief Return location result (size = nDistantpoint)
///

const int *LocationToLocalMesh::getLocation() const
{
  if (_toLocate)
    bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _location;
}

///
/// \brief Return number of distant ranks
///

int LocationToLocalMesh::getNDistantRank() const
{
  if (_toLocate)
    bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _coupledApplicationNRankCouplingComm;
}

///
/// \brief Return located points distribution on distant ranks (size = n_distant_rank + 1)
///

const int *LocationToLocalMesh::getLocatedPointsDistribution() const
{
  if (_toLocate)
    bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _locatedPointsDistribution;
}

///
/// \brief Return distant point distribution on distant ranks (size = n_distant_rank + 1)
///

const int *LocationToLocalMesh::getDistantDistribution() const
{
  if (_toLocate)
    bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _distantDistribution;
}

///
/// \brief Return distance to location cell result (size = nDistantpoint)
///

const float *LocationToLocalMesh::getDistance() const
{
  if (_toLocate)
    bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _distance;
}

///
/// \brief Return number of located distant point
///

int LocationToLocalMesh::getNLocatedDistantPoint() const
{
  if (_toLocate)
    bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _nDistantPoint;
}

///
/// \brief Return number coordinates of located distant point
///

const double *LocationToLocalMesh::getPointCoordinates() const
{
  if (_toLocate)
    bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return fvmc_locator_get_dist_coords(_fvmLocator);
}

///
/// \brief Return fvm locator
///

fvmc_locator_t *LocationToLocalMesh::getFVMLocator() const
{
  if (_toLocate)
    bftc_error(__FILE__, __LINE__, 0,"Call 'locate' before this call !\n");
  return _fvmLocator;
}

///
/// \brief Set support mesh
///
///   @param [in]      supportMesh  location support mesh
///

void LocationToLocalMesh::setSupportMesh(oldMesh *supportMesh, bool lb_tolocate = true)
{
  _toLocate = lb_tolocate;
  _supportMesh = supportMesh;
}

///
/// \brief Set the step for the computaion of the ho optimized element bounding box 
///
///   @param [in]  step
///

void LocationToLocalMesh::optBboxStep(const int step)
{
  _optBboxStep = step;
}

} // Namespace cwipi

#endif /* LOCALLOCATION_HXX_ */
