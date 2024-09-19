/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

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

#include <map>
#include <sstream>

#include <globalData.hxx>
#include "cwp.h"
#include "cwp_priv.h"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  GlobalData::GlobalData(std::string     global_data_id,
                         size_t          s_entity,
                         int             stride,
                         int             n_entity):
  _global_data_id(global_data_id),
  _s_entity(s_entity),
  _stride(stride),
  _n_entity(n_entity),
  _send_data(NULL),
  _recv_data(NULL)
  {
  }

  GlobalData::~GlobalData()
  {}

}

/**
 * \endcond
 */
