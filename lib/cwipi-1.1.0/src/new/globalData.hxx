#ifndef __GLOBALDATA_H__
#define __GLOBALDATA_H__
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
#include <iostream>

#include "cwp.h"
#include "cwp_priv.h"

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class GlobalData globalData.hxx "globalData.hxx"
   * \brief Abstract global data
   *
   *  This class is global data abstract interface
   *
   */


  class GlobalData {

  public:

    /**
     * \brief Constructor
     *
     */
    GlobalData() {}

    GlobalData(std::string     global_data_id,
               size_t          s_entity,
               int             stride,
               int             n_entity);

    /**
     * \brief Destructor
     *
     */
    ~GlobalData();

    inline size_t&
    s_entity_get()
    {
      return _s_entity;
    }

    inline int&
    stride_get()
    {
      return _stride;
    }

    inline int&
    n_entity_get()
    {
      return _n_entity;
    }

    inline void *
    send_data_get()
    {
      return _send_data;
    }

    inline void *
    recv_data_get()
    {
      return _recv_data;
    }

    inline void
    send_data_set
    (
     void *send_data
     )
    {
      _send_data = send_data;
    }

    inline void
    recv_data_set
    (
     void *recv_data
     )
    {
      _recv_data = recv_data;
    }

    inline MPI_Request
    data_send_request_get()
    {
      return _data_send_request;
    }

    inline MPI_Request
    data_recv_request_get()
    {
      return _data_recv_request;
    }

    inline void
    data_send_request_set
    (
     MPI_Request req
    )
    {
      _data_send_request = req;
    }

    inline void
    data_recv_request_set
    (
     MPI_Request req
    )
    {
      _data_recv_request = req;
    }

  private:

    std::string     _global_data_id;

    MPI_Request     _data_send_request;
    MPI_Request     _data_recv_request;

    size_t          _s_entity;
    int             _stride;
    int             _n_entity;

    void           *_send_data;

    void           *_recv_data;

  };

}

/**
 * \endcond
 */

#endif //__GLOBALDATA_H__
