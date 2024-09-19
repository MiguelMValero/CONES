#ifndef __BLOCK_H__
#define __BLOCK_H__
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

#include "cwp_priv.h"
CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wcast-function-type")
#include "cwp.h"
CWP_GCC_SUPPRESS_WARNING_POP


#include "pdm_mesh_nodal.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include <map>
#include <vector>
#include "pdm_error.h"


/**
 * \cond
 */

namespace cwipi {

  class Mesh;

  /**
   * \class Block block.hxx "block.hxx"
   * \brief Mesh mother Block definition
   *
   *  This class defines the Block objects.
   *  Blocks contain elements, are a part of
   *  a Mesh and are defined on Mesh partitions.
   *
   */

  class Block {

  public :

    /**
     *
     * \brief Constructor
     *
     */

    Block();

    /**
     * \brief Destructor.
     *
     */

    virtual ~Block();

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
     *
     * \brief Return the CWIPI block type.
     *
     * \return The CWIPI block type.
     */

    inline CWP_Block_t blockTypeGet();


    /**
     *
     * \brief Return the block identifier.
     *
     */

    inline int blockIDGet();
    inline int blockIDCWIPIGet();

    inline void blockIDCWIPISet(int block_id);

    inline void blockIDPDMSet(int block_id);


    inline int blockIDPDMGet();
    /**
     *
     * \brief return the number of element for each partition
     *
     *
     */

     inline std::vector <int>  NEltsGet();

     inline int  NEltsGet(int i_part);

    /**
     *
     * \brief return the element global numbering for each partition
     *
     *
     */

    CWP_g_num_t* GNumMeshGet(int i_part);


    void GNumMeshSet(int i_part, CWP_g_num_t* gnum, PDM_ownership_t owner);

    void GNumMeshFree(int i_part);


    const double* eltCentersGet(int i_part);

    /**
     *
     * \brief return the number of partition
     *
     *
     */

     inline int NPartGet();


    /**
     *
     * \brief return true is the global numbering (inside the mesh) is defined
     *
     *
     */

    inline bool globNumDefined();


    /**
     *
     * \brief Return true if the Block is already in Paradigm block database and false otherwise.
     *
     */

    inline bool inPDMDB();


    /**
     *
     * \brief Set the _inPDMDB variable to true, and so indicate that the block is already in the PDM database.
     *
     */

    inline void SetinPDMDB();




    void blockSetParentNum(int i_part, int *parent_num, PDM_ownership_t owner);
    inline int *ParentNumGet(int i_part);
    void ParentNumFree(int i_part);



    virtual void geomFinalize() = 0;

  private :

    /**
     *
     * \brief Assigment operator
     *
     */

    Block
    &operator=
    (const Block &other);

    /**
     *
     * \brief Copy constructor
     *
     */

    Block (const Block& other);

  protected:

    CWP_Block_t                _blockType;              /*!< Block Type */
    PDM_Mesh_nodal_t          *_pdmNodal_handle_index;  /*!< PDM Nodal Index */
    MPI_Comm                  *_localComm;              /*!< Communicator */
    int                        _n_part;                 /*!< Number of partitions */
    std::vector <int>          _n_elt;                  /*!< Number of elements for each partition */
    std::vector <double*>      _cells_center;           /*!< Cell centers */
    int                        _block_id_pdm;           /*!< Block identifier */
    int                        _block_id_cwipi;         /*!< Block identifier */
    std::vector <CWP_g_num_t*> _global_num;             /*!< Global numbering in the Mesh  */
    std::vector <int        *> _parent_num;
    Mesh                      *_mesh;                   /*!< Pointer to the mesh object owning the block */
    PDM_ownership_t            _owner_gnum;             /*!< Owner of global numbers */                
    PDM_ownership_t            _owner_parent_num;       /*!< Owner of parent numbers */

  };

  std::vector <int>
  Block::NEltsGet() {
    return _n_elt;
  }


  int
  Block::NEltsGet(int i_part) {
    return _n_elt[i_part];
  }

  int
  Block::NPartGet() {
    return _n_part;
  }

  CWP_Block_t
  Block::blockTypeGet(){
    return _blockType;
  }


  int
  Block::blockIDGet(){
    return _block_id_pdm;
  }

  int
  Block::blockIDCWIPIGet(){
    return _block_id_cwipi;
  }

  void
  Block::blockIDCWIPISet(int block_id){
    _block_id_cwipi = block_id;
  }

  void
  Block::blockIDPDMSet(int block_id){
    _block_id_pdm = block_id;
  }

  int 
  Block::blockIDPDMGet(){
    return _block_id_pdm;
  }

  int *
  Block::ParentNumGet(int i_part){
    return _parent_num[i_part];
  }

}
#endif //__BLOCK_H__

/**
 * \endcond
 */
