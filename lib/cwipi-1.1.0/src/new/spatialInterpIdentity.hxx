#ifndef SPATIALINTERPIDENTITY_H_
#define SPATIALINTERPIDENTITY_H_
/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

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

#include "spatialInterp.hxx"

/**
 * \cond
 */

namespace cwipi {
  class SpatialInterpIdentity : public SpatialInterp {
  public:

    SpatialInterpIdentity();

    ~SpatialInterpIdentity() override;

    void weightsCompute() override;

    private:
        void interpolate(Field *referenceField, double **buffer) override;

        void init (
         Coupling                   *coupling,
         CWP_Dof_location_t         localCodeDofLOcation,
         CWP_Dof_location_t         cplCodeDofLOcation,
         SpatialInterpExchDirection exchDirection) override;

        virtual void clear() override;

//    SpatialInterpIdentity *_spatial_interp_cpl;

    int **_src_to_tgt_idx;

  }; //end SpatialInterpIdentity

/**
 * \endcond
 */

}
#endif // SPATIALINTERPIDENTITY_H_
