#ifndef CWP_SPATIALINTERPLOCATIONOCTREE_H
#define CWP_SPATIALINTERPLOCATIONOCTREE_H
/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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

#include "spatialInterpLocation.hxx"
#include <pdm_mesh_location.h>

namespace cwipi {
    const double CWP_MESH_LOCATION_BBOX_TOLERANCE = 1.e-3;

    class SpatialInterpLocationMeshLocation : public SpatialInterpLocation {
    public:
        SpatialInterpLocationMeshLocation() = default;

        virtual ~SpatialInterpLocationMeshLocation();

    protected:

        void localization_init() override;

        void localization_points_cloud_setting() override;

        void localization_surface_setting() override;

        void localization_compute() override;

        void localization_get() override;

        void localization_free() override;

        PDM_mesh_location_t *_id_pdm;

        PDM_mesh_location_method_t _location_method;

    };

    class SpatialInterpLocationMeshLocationOctree : public SpatialInterpLocationMeshLocation {
    public:
        SpatialInterpLocationMeshLocationOctree() {
            _location_method = PDM_MESH_LOCATION_OCTREE;
        };

        ~SpatialInterpLocationMeshLocationOctree() override = default;
    };

    class SpatialInterpLocationMeshLocationDbbtree : public SpatialInterpLocationMeshLocation {
    public:
        SpatialInterpLocationMeshLocationDbbtree() {
            _location_method = PDM_MESH_LOCATION_DBBTREE;
        };

        ~SpatialInterpLocationMeshLocationDbbtree() override = default;
    };

    class SpatialInterpLocationMeshLocationLocateAllTgt : public SpatialInterpLocationMeshLocation {
    public:
        SpatialInterpLocationMeshLocationLocateAllTgt() {
            _location_method = PDM_MESH_LOCATION_LOCATE_ALL_TGT;
        };

        ~SpatialInterpLocationMeshLocationLocateAllTgt() override = default;
    };

/**
 * \endcond
 */

}
#endif //CWP_SPATIALINTERPLOCATIONOCTREE_H
