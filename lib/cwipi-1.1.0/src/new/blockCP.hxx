#ifndef __blockCP_H__
#define __blockCP_H__
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

namespace cwipi {

  /** 
   * \class BlockCP blockCP.hxx "blockCP.hxx"
   * \brief Polyhedron 3D Block Mesh definition
   *
   *  This class defines the polyhedron 3D Block objects.
   *  It inherits from the Block class through a Factory.
   * 
   */

  class BlockCP :
    public Block {
    
    public:
        
     /**
      * \brief Destructor.
      *
      */
    
      BlockCP();


     /**
      * \brief Destructor.
      *
      */
    
      virtual ~BlockCP();
    
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
       * \param [in]  n_elts     Number of elements of the block in the partition
       * \param [in]  n_faces    Number of faces of the block in the partition
       * \param [in]  connec_faces_idx Vertices to faces connectivity index
       * \param [in]  connec_faces     Vertices to faces connectivity
       * \param [in]  connec_cells_idx Faces to cells connectivity index
       * \param [in]  connec_cells     Faces to cells connectivity
       * \param [in]  global_num Mesh  Global numbering of the block
       *
       */ 
             
       virtual void blockSet(int i_part,int n_elts,int n_faces,
                             int* connec_faces_idx, 
                             int* connec_faces,
                             int* connec_cells_idx,
                             int* connec_cells,
                             CWP_g_num_t* global_num);


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
             
       virtual void blockGet(int         i_part,
                             int         *n_elts,
                             int         *n_faces,
                             int         **connec_faces_idx, 
                             int         **connec_faces,
                             int         **connec_cells_idx,
                             int         **connec_cells,
                             CWP_g_num_t **global_num);


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


       /**
        *
        * \brief return the element faces connectivity (Cells_POLY_3D)
        * for each partition.
        * 
        *
        */
     
        inline virtual std::vector<int*>  ConnecFacesGet();


       /**
        *
        * \brief return the element faces connectivity index (Cells_POLY_3D)
        * for each partition.
        * 
        *
        */
     
        inline virtual std::vector<int*>  ConnecFacesIDXGet();

       /**
        *
        * \brief return the number of faces for each partition (Cells_POLY_3D)
        * 
        *
        */
     
        inline virtual std::vector<int >  NFacesGet();

        inline virtual std::vector<CWP_g_num_t *>  FacesGNumGet();

        inline virtual void FacesGNumSet(std::vector<CWP_g_num_t *> face_ln_to_gn);


        void geomFinalize();

        void FacesGNumFree(int i_part);

    private:
      std::vector<int >          _n_faces;             /*!< Number of faces for each partition */
      std::vector<int*>          _connec_faces_idx;    /*!< Faces connectivity Index for each partition */
      std::vector<int*>          _connec_faces;        /*!< Faces connectivity for each partition */
      std::vector<int*>          _connec_cells_idx;    /*!< Cells onnectivity Index for each partition */
      std::vector<int*>          _connec_cells;        /*!< Cells connectivity for each partition */    
      std::vector<CWP_g_num_t*>  _face_ln_to_gn;       /*!< Face global ids */

  }; //BlockCP Class


  std::vector<int*>  BlockCP::ConnecGet() {
    return _connec_cells;
  }

  std::vector<int*>  BlockCP::ConnecIDXGet() {
  
    return _connec_cells_idx;
  }

  std::vector<int*>  BlockCP::ConnecFacesGet() {
    return _connec_faces;
  }

  std::vector<int*>  BlockCP::ConnecFacesIDXGet() {
    return _connec_faces_idx;
  }

  std::vector<int >  BlockCP::NFacesGet() {
    return _n_faces;
  }

  std::vector<CWP_g_num_t *>  BlockCP::FacesGNumGet() {
    return _face_ln_to_gn;
  }

  void BlockCP::FacesGNumSet(std::vector<CWP_g_num_t *> face_ln_to_gn) {
    _face_ln_to_gn = face_ln_to_gn;
  }


} //namespace cwipi

#endif //blockCP

