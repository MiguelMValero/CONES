#ifndef __SPATIAL_INTERP_I_H__
#define __SPATIAL_INTERP_I_H__
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

#include "coupling.hxx"
#include "coupling_i.hxx"

namespace cwipi {

  /**
   *
   * \brief Return the number of uncomputed targets
   *
   * \return                Number of uncomputed targets
   *
   */

  int
  SpatialInterp::nUncomputedTargetsGet(int i_part) const

  {
    return _n_uncomputed_tgt[i_part];
  }

  /**
   *
   * \brief Return uncomputed targets
   *
   * \return                Uncomputed targets
   *
   */

  const int *
  SpatialInterp::uncomputedTargetsGet(int i_part) const
  {
    return _uncomputed_tgt[i_part];
  }

  /**
   *
   * \brief Return the number of computed targets
   *
   * \return                Number of computed targets
   */

  int
  SpatialInterp::nComputedTargetsGet(int i_part) const
  {
    return _n_computed_tgt[i_part];
  }

  /**
   *
   * \brief Return computed targets
   *
   *
   * \return                Computed targets
   *
   */

  const int *
  SpatialInterp::computedTargetsGet(int i_part) const
  {
    return _computed_tgt[i_part];
  }

  /**
   *
   * \brief Return the number of involved source elements
   *
   * \return                Number of source elements
   */

  int
  SpatialInterp::nInvolvedSourcesGet(int i_part) const
  {
    return _n_involved_sources_tgt[i_part];
  }

  /**
   *
   * \brief Return involved source elements
   *
   * \return                Source elements
   */

  const int *
  SpatialInterp::involvedSourcesGet(int i_part) const
  {
    return _involved_sources_tgt[i_part];
  }

}

#endif //__SPATIAL_INTERP_H__
