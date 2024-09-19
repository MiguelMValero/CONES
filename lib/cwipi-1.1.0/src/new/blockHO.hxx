#ifndef __BlockHO_H__
#define __BlockHO_H__
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
#include "blockStd.hxx"
#include "cwp.h"
#include <map>
#include <vector>
#include <cstdlib>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \class BlockStd blockStd.hxx "blockStd.hxx"
   * \brief Standard Block Mesh definition
   *
   *  This class defines the standard Block objects.
   *  It inherits from the Block class through a Factory.
   *
   */

  class BlockHO :
    public BlockStd {

    public:


     /**
      * \brief Constructor.
      *
      */

      BlockHO();

     /**
      * \brief Destructor.
      *
      */

      ~BlockHO();

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
       * \param [in]  connec     Elements connectivity.
       * \param [in]  global_num Mesh global numbering of the block
       *
       */

       void blockSet(int i_part,int n_elts,int* connec,CWP_g_num_t* global_num,const int order, const char *ho_ordering);


      /**
       * \brief Get a CWIPI block
       *
       * \param [in]  i_part     Partition identifier
       * \param [out]  n_elts     Number of elements of the block in the partition.
       * \param [out]  connec     Elements connectivity.
       * \param [out]  global_num Mesh global numbering of the block
       *
       */

       void blockGet(int i_part,int *n_elts,int** connec,CWP_g_num_t** global_num, int *order, char **ho_ordering);

       /**
        *
        * \brief return the element connectivity (Standard or Face_Poly_2D CWP_Block_t) or cells-faces connectivity (Cells_POLY_3D)
        * for each partition.
        *
        *
        */

        inline std::vector <int*>  ConnecGet();

        inline std::vector <int*>  ConnecIJKGet();

        inline int*  ConnecGet(int i_part);

        inline int*  ConnecIJKGet(int i_part);

        inline int  OrderGet();

        inline char * HOOrderingGet();

        bool  gnumRequired(){
           for(int i_part = 0; i_part<_n_part; i_part++){
             if(_global_num [i_part] == NULL )
               return true;
           }

           return false;
        }

        void geomFinalize();

        void reorder();

    private:

      int   _order;
      char *_ho_ordering;
      std::vector<int*> _connec;     /*!< Connectivity in user order for each partition */
      std::vector<int*> _connec_ijk; /*!< Connectivity in Lagrange (IJK) order for each partition */

  };



  int*  BlockHO::ConnecGet(int i_part) {

    return _connec[i_part];
  }


  std::vector<int*>  BlockHO::ConnecGet() {

    return _connec;
  }

  int*  BlockHO::ConnecIJKGet(int i_part) {

    return _connec_ijk[i_part];
  }

  std::vector<int*>  BlockHO::ConnecIJKGet() {

    return _connec_ijk;
  }

  int  BlockHO::OrderGet() {

    return _order;
  }

  char * BlockHO::HOOrderingGet() {

    return _ho_ordering;
  }

}

/**
 * \endcond
 */


#endif //BlockHO
