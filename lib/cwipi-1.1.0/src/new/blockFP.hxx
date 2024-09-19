#ifndef __BlockFP_H__
#define __BlockFP_H__
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

#include "block.hxx"
#include "cwp.h"
#include <map>
#include <vector>


/**
 * \cond
 */

namespace cwipi {

  /**
   * \class BlockFP blockFP.hxx "blockFP.hxx"
   * \brief Polygon 2D Block Mesh definition
   *
   *  This class defines the polygon 2D Block objects.
   *  It inherits from the Block class through a Factory.
   *
   */

  class BlockFP :
    public Block {

    public:


     /**
      * \brief Destructor.
      *
      */

      BlockFP();


     /**
      * \brief Destructor.
      *
      */

      virtual ~BlockFP();

     /**
       *
       * \brief Block addition
       *
       * Add a block to the mesh.
       *
       * \param [in] blockType              Type of the block
       * \param [in] mesh                   The Mesh object owning the block
       */

      virtual void BlockAdd(CWP_Block_t blockType, Mesh* mesh);

      /**
       * \brief Set a CWIPI block in a partition
       *
       * \param [in]  i_part     Partition identifier
       * \param [in]  n_elts     Number of elements of the block in the partition.
       * \param [in]  connec_idx Elements connectivity index
       * \param [in]  connec     Elements connectivity
       * \param [in]  global_num Mesh global numbering of the block
       *
       */

      virtual void blockSet(int i_part,int n_elts,
                            int* connec_idx,
                            int* connec,
                            CWP_g_num_t* mesh_global_num);

      /**
       * \brief Get a CWIPI block in a partition
       *
       * \param [in]  i_part     Partition identifier
       * \param [in]  n_elts     Number of elements of the block in the partition.
       * \param [in]  connec_idx Elements connectivity index
       * \param [in]  connec     Elements connectivity
       * \param [in]  global_num Mesh global numbering of the block
       *
       */

      virtual void blockGet(int          i_part,
                            int         *n_elts,
                            int         **connec_idx,
                            int         **connec,
                            CWP_g_num_t **mesh_global_num);

       /**
        *
        * \brief return the element connectivity (Standard or Face_Poly_2D CWP_Block_t) or cells-faces connectivity (Cells_POLY_3D)
        * for each partition.
        *
        *
        */

        inline virtual std::vector<int*>  ConnecGet();

       /**
        *
        * \brief return the element connectivity index (Face_Poly_2D CWP_Block_t) or cells-faces connectivity index (Cells_POLY_3D)
        * for each partition.
        *
        *
        */

        inline virtual std::vector<int*>  ConnecIDXGet();
        inline virtual int* ConnecIDXGet(int i_part);

        inline int*  ConnecGet(int i_part);

        bool  gnumRequired(){
           for(int i_part = 0; i_part<_n_part; i_part++){
             if(_global_num [i_part] == NULL )
               return true;
           }

           return false;
        }


        void GNumBlockSet(int i_part, CWP_g_num_t* global_num){
           _global_num [i_part] = global_num;
        }

        void geomFinalize();



    private:
      std::vector<int*>          _connec_idx;          /*!< Connectivity Index for each partition */
      std::vector<int*>          _connec;              /*!< Connectivity for each partition */

  };


  std::vector<int*>  BlockFP::ConnecGet() {

    return _connec;
  }

  std::vector<int*>  BlockFP::ConnecIDXGet() {
    return _connec_idx;
  }

  int*  BlockFP::ConnecIDXGet(int i_part) {
     return _connec_idx[i_part];
  }



  int*  BlockFP::ConnecGet(int i_part) {

    return _connec[i_part];
  }

}


/**
 * \endcond
 */

#endif //BlockFP
