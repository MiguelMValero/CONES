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

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_gnum.h"
#include "block.hxx"
#include <map>
#include <vector>
#include "mesh.hxx"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  Block::Block()
         :_owner_gnum(PDM_OWNERSHIP_USER)
  {

  }


  void Block::BlockAdd(CWP_Block_t blockType, Mesh* mesh)
  {
     _mesh        = mesh;
     _localComm   = const_cast<MPI_Comm*>(static_cast<Mesh*>(mesh)->getMPICommP());

     _blockType   = blockType;
     _n_part      = static_cast<Mesh*>(mesh)->getNPart();

     _global_num  .resize(_n_part,NULL);
     _n_elt       .resize(_n_part);
     _cells_center.resize(_n_part, NULL);
     _parent_num  .resize(_n_part, NULL);

  }


  Block::~Block(){
   for (int i = 0; i < _n_part; i++) {
      if (_owner_gnum == PDM_OWNERSHIP_KEEP && _global_num[i] != NULL) {
        free (_global_num[i]);
      }

    }

  }


  const double* Block::eltCentersGet(int i_part) {

    if (_cells_center[i_part] == NULL) {
      PDM_part_mesh_nodal_section_elt_center_compute(_mesh->getPdmNodalIndex(), _block_id_pdm, i_part, PDM_OWNERSHIP_KEEP);
    }
    return PDM_part_mesh_nodal_section_elt_center_get(_mesh->getPdmNodalIndex(), _block_id_pdm, i_part, PDM_OWNERSHIP_KEEP);
  }

  CWP_g_num_t*
  Block::GNumMeshGet(int i_part) {
    return _global_num[i_part];
  }

  void
  Block::GNumMeshSet(int i_part,CWP_g_num_t* gnum, PDM_ownership_t owner) {
    _global_num[i_part] = gnum;

    _owner_gnum = owner;
  }

  void
  Block::GNumMeshFree(int i_part) {
    if (_owner_gnum == PDM_OWNERSHIP_KEEP) {
      if (_global_num[i_part] != NULL) {
        free(_global_num[i_part]);
        _global_num[i_part] = NULL;
      }
    }
  }


  void Block::blockSetParentNum(int i_part, int* parent_num, PDM_ownership_t owner){
    _parent_num[i_part] = parent_num;
    _owner_parent_num = owner;
  }


  void Block::ParentNumFree(int i_part) {
    if (_owner_parent_num == PDM_OWNERSHIP_KEEP) {
      if (_parent_num[i_part] != NULL) {
        free(_parent_num[i_part]);
        _parent_num[i_part] = NULL;
      }
    }
  }

}

/**
 * \endcond
 */
