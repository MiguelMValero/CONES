#ifndef __COUPLING_PROPERTIES_I_H__
#define __COUPLING_PROPERTIES_I_H__
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

#include <cstring>

#include "locationToLocalMesh.hxx"
#include "locationToDistantMesh.hxx"
#include "bftc_error.h"


namespace cwipi {


  const float *oldCoupling::distance() const
  {
    return &(_distance[0]);
  }

  
  void oldCoupling::set_data_user(void *data)
  {
    _data_user = data;
  }

  void* oldCoupling::get_data_user(void)
  {
    return _data_user;
  }
  
  void  oldCoupling::set_interpolation_function(cwipi_interpolation_fct_t fct)
  {
    _interpolationFct = fct;
  }


  void  oldCoupling::hoOptionsSet(const char* option, const char* value)
  {
    if (!strcmp(option, "opt_bbox_step")) {
      sscanf(value, "%d", &_optBboxStep);
      for (size_t i = 0; i<_tablelocationToLocalMesh.size(); i++) {
        _tablelocationToLocalMesh[i]->optBboxStep(_optBboxStep);
      }
    }
    else {
      bftc_error(__FILE__, __LINE__, 0, "'%s' is not a high order option\n",
                 option);
    }
  }

  
  void  oldCoupling::set_interpolation_function_f(void * fct)
  {
    _interpolationFct_f = fct;
  }
  
  void  oldCoupling::set_ho_interpolation_function(cwipi_user_interp_ho_fct_t fct)
  {
    _ho_interpolationFct = fct;
  }

  void  oldCoupling::set_ho_interpolation_function_f(void * fct)
  {
    _ho_interpolationFct_f = fct;
  }

  const int * oldCoupling::getDistantLocation() const
  {
    return _locationToLocalMesh->getLocation();
  }

  const float * oldCoupling::getDistantDistance() const
  {
    return _locationToLocalMesh->getDistance();
  }

  int oldCoupling::getNNotlocatedPoint() const
  {
    return _locationToDistantMesh->getNUnlocatedPoint();
  }

  const int *oldCoupling::getNotlocatedPoint() const
  {
    return _locationToDistantMesh->getUnlocatedPoint();
  }

  int oldCoupling::getNLocatedPoint() const
  {
    return _locationToDistantMesh->getNLocatedPoint();
  }

  const int *oldCoupling::getLocatedPoint() const
  {
    return _locationToDistantMesh->getLocatedPoint();
  }

  const int *oldCoupling::getDistantBarycentricCoordinatesIndex() const
  {
    return &_locationToLocalMesh->getBarycentricCoordinatesIndex()[0];
  }

  const double *oldCoupling::getDistantBarycentricCoordinates() const
  {
    return &_locationToLocalMesh->getBarycentricCoordinates()[0];
  }

  int oldCoupling::getNDistantPoint() const
  {
    return _locationToLocalMesh->getNLocatedDistantPoint();
  }

  const double *oldCoupling::getDistantPointCoordinates() const
  {
    return &_locationToLocalMesh->getPointCoordinates()[0];
  }

  int oldCoupling::getNDistantRank() const
  {
    return _locationToLocalMesh->getNDistantRank();
  }

  const int *oldCoupling::getLocatedPointsDistribution() const
  {
    return _locationToLocalMesh->getLocatedPointsDistribution();
  }

  const int *oldCoupling::getDistantDistribution() const
  {
    return _locationToLocalMesh->getDistantDistribution();
  }

  ///
  /// \brief Set the type of information to be exchanged at the location
  ///
  ///   @param information
  ///

  void oldCoupling::setInfo(cwipi_located_point_info_t info)
  {
      return _locationToDistantMesh->setInfo(info);
  }

  ///
  /// \brief Get list of number of vertices of containing element
  ///

  const int *oldCoupling::getElementContainingNVertex() const
  {
    return _locationToDistantMesh->getElementContainingNVertex();
  }

  ///
  /// \brief Get connectivity of element that contains each located point
  ///

  const int *oldCoupling::getElementContainingVertex() const
  {
    return _locationToDistantMesh->getElementContainingVertex();
  }

  ///
  /// \brief Get vertices coordinates of the element that contains each located point
  ///

  const double *oldCoupling::getElementContainingVertexCoords() const
  {
    return _locationToDistantMesh->getElementContainingVertexCoords();
  }

  ///
  /// \brief Get barycentric coordinates for the element that contains each located point
  ///

  const double *oldCoupling::getElementContainingBarycentricCoordinates() const
  {
    return _locationToDistantMesh->getElementContainingBarycentricCoordinates();
  }

  ///
  /// \brief Get MPI rank that contains each located point
  ///

  const int *oldCoupling::getElementContainingMPIrank() const
  {
    return _locationToDistantMesh->getElementContainingMPIrank();
  }


  ///
  /// \brief Get distant element that contain located point
  ///

  const int *oldCoupling::getElementContaining() const
  {
    return _locationToDistantMesh->getElementContaining();
  }

} // name space cwipi


#endif //__COUPLING_PROPERTIES_I_H__
