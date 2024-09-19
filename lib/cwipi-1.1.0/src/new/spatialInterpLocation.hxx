#ifndef SPATIALINTERPLOCATION_H_
#define SPATIALINTERPLOCATION_H_
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

#include "spatialInterp.hxx"
#include "pdm_gnum_location.h"
#include "pdm_part_mesh_nodal.h"

/**
 * \cond
 */

namespace cwipi {
  class SpatialInterpLocation : public SpatialInterp {
  public:

    /**
      *
      * \brief SpatialInterp location constructor.
      *
      */
    SpatialInterpLocation();


    /**
      *
      * \brief SpatialInterp location Init.
      *
      */

    virtual void 
    init (
      Coupling           *coupling, 
      CWP_Dof_location_t localCodeDofLOcation,
      CWP_Dof_location_t cplCodeDofLOcation,
      SpatialInterpExchDirection exchDirection ) override;

    virtual void clear() override;


    /**
      *
      * \brief SpatialInterp location destructor.
      *
      */
    virtual ~SpatialInterpLocation();

    /**
      *
      * \brief Compute of the spatial interpolation weights. Localization and communication
      *        tree building.
      *
      * \param [in] Texch_t    Type of exchange (sending or reception).
      *
      */
    void weightsCompute() override;

    // Get point_*

    inline double **
    points_coords_get()
    {
      return _points_coords;
    }

    inline double **
    points_uvw_get()
    {
      return _points_uvw;
    }

    inline double **
    points_dist2_get()
    {
      return _points_dist2;
    }

    inline double **
    points_projected_coords_get()
    {
      return _points_projected_coords;
    }

    // Get internal cell_vtx ordering

    inline int **
    cell_vtx_idx_get()
    {
      return _cell_vtx_idx;
    }

    inline int **
    cell_vtx_get()
    {
      return _cell_vtx;
    }

    // A mettre en privé !!!

    SpatialInterpLocation *_spatial_interp_cpl;

    // Localization data

  protected:

    PDM_part_mesh_nodal_t* _pdm_CplNodal;

      //
    // Target properties
    
    double **_tgt_distance;                 // Distance to the closest source element surface by partition
    double **_tgt_projected;                // Projected point coordinates (on the closest source element surface)
    CWP_g_num_t **_tgt_closest_elt_gnum;    // Closest source element global numbering

    //
    // Source properties
    int         **_elt_pts_inside_idx;
    PDM_g_num_t **_points_gnum;
    double      **_points_coords;
    double      **_points_uvw;
    double      **_points_dist2;
    double      **_points_projected_coords;

    int         **_cell_vtx_idx;
    int         **_cell_vtx;

    int fake_idx[1] = {0};

  protected:
    /**
      *
      * \brief Interpolation of a point cloud on a reference field.
      *
      * \param [in]   referenceField   Reference field pointer
      *
      */
    void interpolate(Field *referenceField, double **buffer) override;


    /***********************************************************
     ***********************************************************
     **                                                       **
     **             object functions              **
     **                                                       **
     ***********************************************************
     ***********************************************************/
    /**
       *
       * \brief Setting of the points cloud for localization.
       *
       * \param [out] id_dist   Localization object identifier.
       *
       */
    virtual void localization_init();

    /**
       *
       * \brief Setting of the points cloud for localization.
       *
       * \param [out] id_dist   Localization object identifier.
       *
       */
    virtual void localization_points_cloud_setting();

    /**
      *
      * \brief Setting of the surface mesh and cloud points at
      * null for the localization object in a case of receving
      *    * code i.e. code which provides cloud points for interpolation.
      *
      * \param [out] id_dist   Localization object identifier.
      *
      */
    virtual void localization_surface_setting();

    /**
      *
      * \brief Compute of localization of a points cloud on a surface
      *  mesh through the localization object.
      *
      * \param [int] id_dist   Localization object identifier.
      *
      */
    virtual void localization_compute();

    /**
      *
      * \brief Get localization results from localization object.
      *
      * \param [int] id_dist   Localization object identifier.
      *
      */
    virtual void localization_get();

    /**
      *
      * \brief Get localization results from localization object
      * from the coupled code in the case where the both codes are on
      * the same process.
      *
      * \param [int] id_dist   Localization object identifier.
      *
      */

    virtual void localization_free();

  
  }; //end SpatialInterpLocation

/**
 * \endcond
 */

}
#endif // SPATIALINTERPLOCATION_H_
