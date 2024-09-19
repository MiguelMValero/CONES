

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

#include "blockStd.hxx"
#include "cwp.h"
#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <pdm_error.h>
#include <pdm_printf.h>
#include <map>
#include <vector>

#include <mesh.hxx>

/**
 * \cond
 */

namespace cwipi {
  BlockStd::BlockStd()
     :Block::Block()
  {

  }


  BlockStd::~BlockStd()
  {
    _connec.clear();
  }


  void BlockStd::BlockAdd(CWP_Block_t blockType, Mesh* mesh)
  {
    Block::BlockAdd(blockType, mesh);
    _connec.resize(_n_part, NULL); 
  }


  void BlockStd::blockSet(int i_part,int n_elt,int* connec,CWP_g_num_t* mesh_global_num){


     _global_num [i_part] = mesh_global_num;

     _n_elt[i_part] = n_elt;
     _connec[i_part] = connec;
  }


  void BlockStd::blockGet(int i_part,int *n_elts,int** connec,CWP_g_num_t** global_num){
    *global_num      = _global_num [i_part];
    *n_elts          = _n_elt[i_part];
    *connec          = _connec[i_part];
  }


  void BlockStd::geomFinalize(){

    // for(int i_part = 0; i_part<_n_part; i_part++){

    //   Visu* visu = ((Mesh*)_mesh)->getVisu();
    //   if(visu->isCreated() && ((Mesh*)_mesh)->getDisplacement() == CWP_DYNAMIC_MESH_STATIC) {
    //     visu->GeomBlockStdSet ( ((Mesh*)_mesh)->getIdVisu( _block_id_cwipi ),
    //                               i_part,
    //                               _n_elt[i_part] ,
    //                               _connec[i_part],
    //                               _global_num [i_part]);
    //   }
    // }
  }


}


/**
 * \endcond
 */
