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

#include <map>
#include <vector>
#include <sstream>
#include <cassert>

#include <partData.hxx>
#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_part_to_part.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_error.h"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  /**
    * \brief Constructor
    *
    */
  PartData::PartData(std::string           part_data_id,
                     CWP_PartData_exch_t   exch_type,
                     CWP_g_num_t         **gnum_elt,
                     int                  *n_elt,
                     int                   n_part):
  _part_data_id(part_data_id),
  _exch_type(exch_type),
  _ptp(NULL)
  {
    _gnum_elt = gnum_elt;
    _n_elt    = n_elt;
    _n_part   = n_part;
  }


  /**
    * \brief Destructor
    *
    */
  PartData::~PartData()
  {}


  void
  PartData::data_set(int    exch_id,
                     void **data)
  {
    map<int, void **>::iterator it = _data.find(exch_id);

    if (it != _data.end()) {
      PDM_error(__FILE__, __LINE__, 0, "PartData '%s': Data pointer with exch_type %d already set for exch_id %d\n",
                _part_data_id.c_str(), _exch_type, exch_id);
    }

    pair<int, void **> new_pair(exch_id, data);
    _data.insert(new_pair);
  };


  void
  PartData::request_set(int exch_id,
                        int request)
  {
    map<int, int>::iterator it = _request.find(exch_id);

    if (it != _request.end()) {
      PDM_error(__FILE__, __LINE__, 0, "PartData '%s': Request pointer with exch_type %d already set for exch_id %d\n",
                _part_data_id.c_str(), _exch_type, exch_id);
    }

    pair<int, int> new_pair(exch_id, request);
    _request.insert(new_pair);
  };


  int
  PartData::request_get(int exch_id)
  {
    map<int, int>::iterator it = _request.find(exch_id);

    if (it == _request.end()) {
      PDM_error(__FILE__, __LINE__, 0, "PartData '%s': Request for exch_type %d undiefined set for exch_id %d\n",
                _part_data_id.c_str(), _exch_type, exch_id);
    }

    return it->second;
  };


  /* In-place filtering of data coming from multiple origins */
  void
  PartData::recv_data_filter(int exch_id)
  {
    assert(_exch_type == CWP_PARTDATA_RECV);

    map<int, void **>::iterator it = _data.find(exch_id);

    if (it == _data.end()) {
      PDM_error(__FILE__, __LINE__, 0, "PartData '%s': Recv data pointer with was not set for exch_id %d\n",
                _part_data_id.c_str(), exch_id);
    }

    map<int, int>::iterator it_s_unit = _s_unit.find(exch_id);
    int s_unit = it_s_unit->second;

    int         **come_from_idx = NULL;
    PDM_g_num_t **come_from     = NULL;
    PDM_part_to_part_gnum1_come_from_get(_ptp,
                                         &come_from_idx,
                                         &come_from);

    int n_part1 = 0;
    int n_part2 = 0;
    PDM_part_to_part_n_part_get(_ptp,
                                &n_part1,
                                &n_part2);


    unsigned char **data = (unsigned char **) it->second;

    for (int ipart = 0; ipart < n_part2; ipart++) {
      int  n_ref = 0;
      int *ref   = NULL;
      PDM_part_to_part_ref_lnum2_single_part_get(_ptp,
                                                 ipart,
                                                 &n_ref,
                                                 &ref);
      for (int i = 0; i < n_ref; i++) {
        int idx = come_from_idx[ipart][i];
        for (int j = 0; j < s_unit; j++) {
          data[ipart][s_unit*i + j] = data[ipart][s_unit*idx + j];
        }
      }
    }
  };


  void
  PartData::exch_clear(int exch_id)
  {
    _data   .erase(exch_id);
    _request.erase(exch_id);

    _s_unit.erase(exch_id);
  }


}

/**
 * \endcond
 */
