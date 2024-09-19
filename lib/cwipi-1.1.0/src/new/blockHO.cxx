

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

#include "blockHO.hxx"
#include "cwp.h"
#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <pdm_error.h>
#include <pdm_printf.h>
#include <pdm_logging.h>
#include <pdm_ho_ordering.h>
#include <map>
#include <vector>

#include <mesh.hxx>
#include <string.h>

/**
 * \cond
 */

namespace cwipi {
  BlockHO::BlockHO()
     :BlockStd::BlockStd()
  {

    _ho_ordering = NULL;
  }


  BlockHO::~BlockHO()
  {
    _connec.clear();
    for (size_t i = 0; i < _connec_ijk.size(); i++) {
      if (_connec_ijk[i] != NULL) {
        free(_connec_ijk[i]);
      }
    }
    _connec_ijk.clear();
    if (_ho_ordering != NULL) {
      free(_ho_ordering);
    }
  }


  void BlockHO::BlockAdd(CWP_Block_t blockType, Mesh* mesh)
  {
    BlockStd::BlockAdd(blockType, mesh);
    _connec.resize(_n_part, NULL); 
    _connec_ijk.resize(_n_part, NULL);
  }


  void BlockHO::blockSet(int i_part,int n_elt,int* connec,CWP_g_num_t* mesh_global_num, const int order, const char *ho_ordering){


    _global_num [i_part] = mesh_global_num;

    _n_elt[i_part] = n_elt;
    _connec[i_part] = connec;

    PDM_Mesh_nodal_elt_t elt_type = (PDM_Mesh_nodal_elt_t) _blockType;
    int elt_node_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    _connec_ijk[i_part] = (int *) malloc(sizeof(int) * n_elt * elt_node_n);
    memcpy(_connec_ijk[i_part], _connec[i_part], sizeof(int) * n_elt * elt_node_n);

    _order = order;
    if (_ho_ordering != NULL) {
      free(_ho_ordering);
      _ho_ordering = NULL;
    }
    if (ho_ordering != NULL) {
      _ho_ordering = (char *) malloc(sizeof(char) * (strlen(ho_ordering) + 1));
      strcpy(_ho_ordering, ho_ordering);
    }
  }


  void BlockHO::blockGet(int i_part,int *n_elts,int** connec,CWP_g_num_t** global_num, int *order, char **ho_ordering) {
    *global_num      = _global_num [i_part];
    *n_elts          = _n_elt[i_part];
    *connec          = _connec[i_part];
    *order           = _order;
    *ho_ordering     = _ho_ordering;
  }


  void BlockHO::geomFinalize(){
    reorder();
  }

  void BlockHO::reorder() {

    if (PDM_ho_ordering_id_get(_ho_ordering) < 0) {
      PDM_error(__FILE__, __LINE__, 0,
                "CWP_Mesh_interf_ho_ordering_from_IJK_set must be called to provide CWIPI with the user's high-order numbering\n");
    }

    PDM_Mesh_nodal_elt_t elt_type = (PDM_Mesh_nodal_elt_t) _blockType;

    for (int ipart = 0; ipart < _n_part; ipart++) {
      PDM_Mesh_nodal_reorder_elt_vtx(elt_type,
                                     _order,
                                     _ho_ordering,
                                     NULL,
                                     _n_elt     [ipart],
                                     _connec    [ipart],
                                     _connec_ijk[ipart]);
    }
  }


}


/**
 * \endcond
 */
