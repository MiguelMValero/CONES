

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

#include "blockCP.hxx"
#include "cwp.h"
#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <map>
#include <vector>
#include <mesh.hxx>
#include <pdm_priv.h>

namespace cwipi {
  BlockCP::BlockCP()
     :Block::Block()
  {
  }


  BlockCP::~BlockCP()
  {
  
  }


  void BlockCP::BlockAdd(CWP_Block_t blockType, Mesh* mesh)
  {
    Block::BlockAdd(blockType, mesh);
    _n_faces.resize(_n_part);
    _connec_cells_idx.resize(_n_part,NULL);
    _connec_cells.resize(_n_part,NULL);
    _connec_faces_idx.resize(_n_part,NULL);
    _connec_faces.resize(_n_part,NULL);
    _face_ln_to_gn.resize(_n_part,NULL);
  }

  void BlockCP::blockSet(int          i_part,
                         int          n_elts,
                         int          n_faces,
                         int*         connec_faces_idx,
                         int*         connec_faces,
                         int*         connec_cells_idx,
                         int*         connec_cells,
                         CWP_g_num_t* global_num) 
  {


    _global_num[i_part] = global_num;

    _n_elt[i_part]            = n_elts;
    _n_faces[i_part]          = n_faces;
    _connec_cells_idx[i_part] = connec_cells_idx;  
    _connec_cells[i_part]     = connec_cells;   
    _connec_faces_idx[i_part] = connec_faces_idx;  
    _connec_faces[i_part]     = connec_faces;
  
  }



   /**
   * \brief Get a CWIPI block in a partition
   * 
   * \param [in]  i_part     Partition identifier
   * \param [out]  n_elts     Number of elements of the block in the partition
   * \param [out]  n_faces    Number of faces of the block in the partition
   * \param [out]  connec_faces_idx Vertices to faces connectivity index
   * \param [out]  connec_faces     Vertices to faces connectivity
   * \param [out]  connec_cells_idx Faces to cells connectivity index
   * \param [out]  connec_cells     Faces to cells connectivity
   * \param [out]  global_num Mesh  Global numbering of the block
   *
   */ 
         
   void BlockCP::blockGet(int         i_part,
                         int         *n_elts,
                         int         *n_faces,
                         int         **connec_faces_idx, 
                         int         **connec_faces,
                         int         **connec_cells_idx,
                         int         **connec_cells,
                         CWP_g_num_t **global_num)
   {

    *global_num       = _global_num[i_part];
    *n_elts           = _n_elt[i_part];
    *n_faces          = _n_faces[i_part];
    *connec_cells_idx = _connec_cells_idx[i_part];  
    *connec_cells     = _connec_cells[i_part];   
    *connec_faces_idx = _connec_faces_idx[i_part];  
    *connec_faces     = _connec_faces[i_part];
  
   }


  void
  BlockCP::FacesGNumFree(int i_part) {
    // if (_owner_gnum == PDM_OWNERSHIP_KEEP) {
      if (_face_ln_to_gn[i_part] != NULL) {
        free(_face_ln_to_gn[i_part]);
        _face_ln_to_gn[i_part] = NULL;
      }
    // }
  }



  void BlockCP::geomFinalize(){

    // for(int i_part = 0; i_part<_n_part; i_part++){

    //   Visu* visu = ((Mesh*)_mesh)->getVisu();
    //   if(visu->isCreated() && ((Mesh*)_mesh)->getDisplacement() == CWP_DYNAMIC_MESH_STATIC) {
    //     visu->GeomBlockPoly3D( ((Mesh*)_mesh)->getIdVisu( _block_id_cwipi ),
    //                               i_part,
    //                               _n_elt[i_part] ,
    //                               _n_faces[i_part],
    //                               _connec_faces_idx [i_part],
    //                               _connec_faces     [i_part],
    //                               _connec_cells_idx[i_part],
    //                               _connec_cells    [i_part],
    //                               _global_num [i_part]);


    //   }
    // } //end i_part


  }
}

