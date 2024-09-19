#ifndef __PARTDATA_H__
#define __PARTDATA_H__

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
#include <iostream>

#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_part_to_part.h"

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class PartData partData.hxx "partData.hxx"
   * \brief Abstract partionned data
   *
   *  This class is partionned data abstract interface
   *
   */


  class PartData {

  public:

    /**
     * \brief Constructor
     *
     */
    PartData(std::string           part_data_id,
             CWP_PartData_exch_t   exch_type,
             CWP_g_num_t         **gnum_elt,
             int                  *n_elt,
             int                   n_part);

    /**
     * \brief Destructor
     *
     */
    ~PartData();

    void
    data_set(int    exch_id,
             void **data);

    void
    request_set(int exch_id,
                int request);

    int
    request_get(int exch_id);

    void
    recv_data_filter(int exch_id);

    void
    exch_clear(int exch_id);

    inline PDM_part_to_part_t *
    ptp_get()
    {
      return _ptp;
    }

    inline void
    ptp_set(PDM_part_to_part_t *ptp)
    {
      _ptp = ptp;
    }

    inline CWP_g_num_t **
    gnum_elt_get()
    {
      return _gnum_elt;
    }

    inline int *
    n_elt_get()
    {
      return _n_elt;
    }

    inline int
    n_part_get()
    {
      return _n_part;
    }

    inline CWP_PartData_exch_t
    exch_type_get()
    {
      return _exch_type;
    }

  private:

    std::string              _part_data_id;
    CWP_PartData_exch_t      _exch_type;
    PDM_part_to_part_t      *_ptp;

    CWP_g_num_t            **_gnum_elt;
    int                     *_n_elt;
    int                      _n_part;
    std::map<int, void **>   _data;

    /* Internal */
    std::map<int, int>       _s_unit;
    std::map<int, int>       _request;
  };

}

/**
 * \endcond
 */

#endif //__PARTDATA_H__
