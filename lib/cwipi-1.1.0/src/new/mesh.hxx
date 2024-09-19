#ifndef __MESH_H__
#define __MESH_H__
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

#include <vector>
#include <map>

#include <mpi.h>

#include <pdm_mesh_nodal.h>
#include <pdm_part_mesh_nodal.h>
#include <pdm_printf.h>
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "block.hxx"
#include "blockFP.hxx"
#include "blockCP.hxx"
#include "blockStd.hxx"
#include "blockHO.hxx"
#include "cwp.h"
// #include "visualization.hxx"


/**
 * \cond
 */



namespace cwipi {


  /**
   * \class Mesh mesh.hxx "mesh.hxx"
   * \brief Interface mesh
   *
   *  This class defines the interface mesh objects.
   *
   */

  class Coupling;
  // class Visu;
  class Mesh {

  public:

    /**
     * \brief Mesh constructor
     *
     * Construct the CWIPI mesh by using paradigm nodal methods.
     *
     * \param [in] localComm Coupling Communicator.
     * \param [in] npart     Number of mesh partitions.
     * \param [in] visu      Pointer to the Visu object
     *
     */

    Mesh
    (
      const MPI_Comm    &localComm,
      int                npart,
      CWP_Dynamic_mesh_t displacement,
      Coupling           *cpl
    );


    /**
     * \brief Mesh destructor
     *
     * Destroy the Mesh object.
     *
     */

    virtual ~Mesh();

    /**
     * \brief Set the mesh coordinates
     *
     * \param [in] i_part       Index of the mesh partition
     * \param [in] n_pts        Number of partition vertices
     * \param [in] coords       Array containing the verticies coordinates
     * \param [in] global_num   Global numbering of the vertices or NULL
     *
     */

    void 
    coordSet 
    (
      const int   i_part,
      const int   n_pts,
      double      coords[],
      CWP_g_num_t global_num[]
    );


    /**
     * \brief Mesh deletion and free memory
     *
     */

    void 
    meshDel();

    /**
     * \brief Addition of a block (set of cells) to the mesh partition
     *
     * This function add a block to the geometric mesh.
     *
     * \param [in] block_type  Type of the block addition
     *
     * \return block_id  Block Identifier
     */

    int 
    blockAdd 
    (
      const CWP_Block_t  block_type
    );


    /**
     * \brief Set a standard block to the interface mesh
     *
     * \param [in] i_part      Partition identifier
     * \param [in] block_id    Block identifier
     * \param [in] n_elts      Number of block elements
     * \param [in] connec      Vertices to elements connectivity
     * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
     *
     */

    void
    stdBlockSet
    (
      const int    i_part,
      const int    block_id,
      const int    n_elts,
      int          connec[],
      CWP_g_num_t  global_num[]
    );


    /**
     * \brief Set a high-order block to the interface mesh
     *
     * \param [in] i_part      Partition identifier
     * \param [in] block_id    Block identifier
     * \param [in] n_elts      Number of block elements
     * \param [in] connec      Vertices to elements connectivity
     * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
     * \param [in] order       Element order
     * \param [in] ho_ordering HO ordering name
     *
     */

    void
    HOBlockSet
    (
      const int    i_part,
      const int    block_id,
      const int    n_elts,
      int          connec[],
      CWP_g_num_t  global_num[],
      const int    order,
      const char  *ho_ordering
    );

    /**
     * \brief Get a standard block from the interface mesh
     *
     * \param [in]  i_part      Partition identifier
     * \param [in]  block_id    Block identifier
     * \param [out] n_elts      Number of block elements
     * \param [out] connec      Vertices to elements connectivity
     * \param [out] global_num  Global numbering of the vertices in the block (or NULL)
     *
     */

    void
    stdBlockGet
    (
      const int     i_part,
      const int     block_id,
      int          *n_elts,
      int         **connec,
      CWP_g_num_t **global_num
    );

    /**
     * \brief Get a high-order block from the interface mesh
     *
     * \param [in]  i_part      Partition identifier
     * \param [in]  block_id    Block identifier
     * \param [out] n_elts      Number of block elements
     * \param [out] order       Element order
     * \param [out] connec      Vertices to elements connectivity
     * \param [out] global_num  Global numbering of the vertices in the block (or NULL)
     *
     */

    void
    HOBlockGet
    (
      const int     i_part,
      const int     block_id,
      int          *n_elts,
      int          *order,
      int         **connec,
      CWP_g_num_t **global_num
    );

    /**
     * \brief Get the standard block type
     *
     * \param [in] block_id    Block identifier
     *
     * \return block type
     *
     */

    CWP_Block_t
    stdBlockTypeGet
    (
     const int              block_id
    );

    /**
     * \brief Set a face polygon block to the interface mesh
     *
     * \param [in] i_part      Partition identifier
     * \param [in] block_id    Block identifier
     * \param [in] n_elts      Number of block elements
     * \param [in] connec_idx  Vertices to elements connectivity index
     * \param [in] connec      Vertices to elements connectivity
     * \param [in] global_num  Global numbering of the vertices in the block (or NULL)
     *
     */

    void 
    poly2DBlockSet
    ( 
      const int    i_part,
      const int    block_id,
      const int    n_elts,
      int          connec_idx[],
      int          connec[],
      CWP_g_num_t  global_num[]
    );

    /**
     * \brief Get a face polygon block from the interface mesh
     *
     * \param [in] i_part      Partition identifier
     * \param [in] block_id    Block identifier
     * \param [out] n_elts      Number of block elements
     * \param [out] connec_idx  Vertices to elements connectivity index
     * \param [out] connec      Vertices to elements connectivity
     * \param [out] global_num  Global numbering of the vertices in the block (or NULL)
     *
     */

    void 
    poly2DBlockGet
    ( 
      const int    i_part,
      const int    block_id,
      int         *n_elts,
      int         **connec_idx,
      int         **connec,
      CWP_g_num_t **global_num
    );


    /**
     * \brief Set a face polhedron block to the interface mesh
     *
     * \param [in] i_part            Partition identifier
     * \param [in] block_id          Block identifier
     * \param [in] n_elts            Number of block elements
     * \param [in] n_faces           Number of faces elements (or NULL)
     * \param [in] connec_faces_idx  Vertices to faces connectivity index
     * \param [in] connec_faces      Vertices to faces connectivity
     * \param [in] connec_cells_idx  Faces to cells connectivity index
     * \param [in] connec_cells      Faces to cells connectivity
     * \param [in] global_num        Global numbering of the vertices in the block (or NULL)
     *
     */

    void 
    poly3DBlockSet
    ( 
      const int   i_part,
      const int   block_id,
      const int   n_elts,
      const int   n_faces,
      int         connec_faces_idx[],
      int         connec_faces[],
      int         connec_cells_idx[],
      int         connec_cells[],
      CWP_g_num_t global_num[]
    );


    /**
     * \brief get a face polhedron block 
     *
     * \param [in] i_part            Partition identifier
     * \param [in] block_id          Block identifier
     * \param [out] n_elts            Number of block elements
     * \param [out] n_faces           Number of faces elements (or NULL)
     * \param [out] connec_faces_idx  Vertices to faces connectivity index
     * \param [out] connec_faces      Vertices to faces connectivity
     * \param [out] connec_cells_idx  Faces to cells connectivity index
     * \param [out] connec_cells      Faces to cells connectivity
     * \param [out] global_num        Global numbering of the vertices in the block (or NULL)
     *
     */

    void 
    poly3DBlockGet
    ( 
      const int    i_part,
      const int    block_id,
      int         *n_elts,
      int         *n_faces,
      int         **connec_faces_idx,
      int         **connec_faces,
      int         **connec_cells_idx,
      int         **connec_cells,
      CWP_g_num_t **global_num
    );


    /**
     * \brief Adding a polyhedron block to the geometric mesh from
     * a face-to-cell connectivity and a vertices-to-faces connectivity.
     *
     * This function add a polyhedron 3D block to the geometric mesh from
     * a face-to-cell connectivity and a vertices-to-faces connectivity.
     *
     * \param [in]  i_part            Current partition
     * \param [in]  n_cells           Number of elements
     * \param [in]  cell_face_idx     Polyhedron to face index
     *                                (src_poly_cell_face_idx[0] = 0 and
     *                                 size = n_elts + 1)
     * \param [in]  cell_face         Polyhedron to face connectivity
     *                                (size = cell_face_idx[n_elts])
     * \param [in]  n_faces           Number of faces
     * \param [in]  face_vtx_idx      Polyhedron vertices to faces index
     *                                (face_vtx_idx[0] = 0 and
     *                                 size_idx = max(face_vtx) + 1)
     * \param [in]  face_vtx          Polyhedron vertices to faces connectivity
     *                                (size = face_vtx_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     *
     */

     void fromCellFaceSet(const int   i_part,
                          const int   n_cells,
                          int         cell_face_idx[],
                          int         cell_face[],
                          int         n_faces,
                          int         face_vtx_idx[],
                          int         face_vtx[],
                          CWP_g_num_t parent_num[]);

    /**
     * \brief Adding a polygon 2D block to the geometric mesh from
     * a vertices-to-faces connectivity and a edge-to-face connectivity.
     *
     * This function add a polygon 2D block to the geometric mesh from
     * a vertices-to-faces connectivity and a edge-to-face connectivity.
     *
     * \param [in]  i_part            Current partition
     * \param [in]  n_faces           Number of faces
     * \param [in]  face_edge_idx     Polygon vertices to faces index
     *                                (face_edge_idx[0] = 0 and
     *                                size_idx = max(face_edge) + 1)
     * \param [in]  face_edge         Polyhegon vertices to face connectivity
     *                                (size = face_edge_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     * \param [in]  n_edges           Number of edges
     * \param [in]  edge_vtx          Polygon vertices to edges connectivity
     *                                (size = edge_vtx_idx[size_idx - 1])
     * \param [in]  parent_num        Pointer to parent element number (or NULL)
     *
     */

    void fromFacesEdgeSet(const int   i_part,
                          const int   n_faces,
                          int         face_edge_idx[],
                          int         face_edge[],
                          const int   n_edges,
                          int         edge_vtx[],
                          CWP_g_num_t parent_num[]);

    void fromFacesVtxSet(const int   i_part,
                         const int   n_faces,
                         int         face_vtx_idx[],
                         int         face_vtx[],
                         CWP_g_num_t global_num[]);


    /**
    * \brief Update the block database.
    *
    * This function updates the block database .
    *
    *
    */

    void updateBlockDB();

    /**
    * \brief Get the vertices number of the partition i_part
    *
    * This function gets the vertices number of the partition i_part.
    *
    * \param [in] i_part   Partition identifier
    * \return Vertex number of i_part partition.
    *
    */

    inline int getPartNVertex(int i_part) const;

    /**
    * \brief Get the vertices coordinates of the i_part partition
    *
    * This function gets the vertices number of the i_part partition.
    *
    * \param [in] i_part   Partition identifier
    * \return Vertices coordinates of the i_part partition.
    *
    */

    inline double* getVertexCoords(int i_part);


    inline CWP_g_num_t* getVertexGNum(int i_part);


    /**
    * \brief Get the paradigm nodal object associated with the mesh.
    *
    * This function gets the paradigm nodal object associated with the mesh.
    *
    * \return Paradigm nodal object associated with the mesh.
    *
    */

    inline PDM_part_mesh_nodal_t& getPdmNodal();

    /**
    * \brief Get the number of elements of the id_block block.
    *
    * This function gets the number of elements of the id_block block.
    *
    * \param [in] id_block   Block identifier
    * \return                Number of elements of the id_block block
    *
    */

    inline int getBlockNElts(int id_block, int i_part);

    /**
    * \brief Get the number of elements of the id_part partition
    *
    * This function gets the number of elements of the id_part partition.
    *
    * \param [in] id_part    Partition identifier
    * \return                Number of elements of the id_part partition
    *
    */

    int getPartNElts(int id_part) const;

    /**
    * \brief Get the number of polyhedra in a partition
    *
    * This function gets the number of polyhedra in a partition.
    *
    * \param [in] i_part   Partition identifier
    * \return              Number of polydra of the i_part partition
    *
    */

    inline int getPartNPolyhedra(int i_part) const;

    /**
    * \brief Get a block element connectivity index
    *
    * This function gets the element connectivity index of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Element connectivity index of the id_block block
    *
    */

    inline int* getPoly2DConnectivityIndex(int id_block,int i_part);

    /**
    * \brief Get a block element connectivity
    *
    * This function gets the element connectivity of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Element connectivity of the id_block block
    *
    */

    inline int* getPoly2DConnectivity     (int id_block,int i_part);


    /**
    * \brief Get a block element connectivity
    *
    * This function gets the element connectivity of the id_block block.
    *
    * \param [in] id_block Block identifier
    * \return              Element connectivity of the id_block block
    *
    */

    inline int* getStdConnectivity     (int id_block,int i_part);

    /**
    * \brief True if coordinates are defined on all partitions False otherwise.
    *
    *
    */

    inline bool coordsDefined ();

   // inline const std::vector<double>& getVolume();

   // inline const std::vector<double>& getCellCenterCoords();

   /**
    * \brief Return pdmNodal Handler identifier
    *
    *
    */
    inline  PDM_part_mesh_nodal_t * getPdmNodalIndex();

   /**
    * \brief Return Coordinates of the i_part partition.
    *
    * \param [in] i_part partition identifier.
    * \return Coordinates of the i_part partition.
    *
    */

   inline double* getCoordinates(int i_part);


   /**
    * \brief Return MPI communicator
    *
    *
    */

   inline MPI_Comm getMPIComm();

   /**
    * \brief Return MPI communicator pointer
    *
    *
    */

   inline const MPI_Comm* getMPICommP();

   /**
    * \brief Return number of partitions.
    *
    *
    */

   inline int getNPart();

    //-->>
    inline int getNCell(int id_part) const;

    inline int getNFace(int id_part) const;

    inline int getNEdge(int id_part) const;

    inline int *getCellFaceIndex(int i_part);

    inline int *getCellFace(int i_part);

    inline int *getFaceEdgeIndex(int i_part);

    inline int *getFaceEdge(int i_part);

    inline int *getFaceVtxIndex(int i_part);

    inline int *getFaceVtx(int i_part);

    inline int *getEdgeVtxIndex(int i_part);

    inline int *getEdgeVtx(int i_part);
    //<<--

   void geomFinalize();

   inline bool gnumVtxRequired ();

   inline CWP_Dynamic_mesh_t getDisplacement();
   inline int getIdVisu(int block_id);

   /**
    * \brief Set the Visu pointer object
    *
    * \param [in] visu Pointer to the Visu object
    *
    */

   // inline void setVisu(Visu* visu);

   CWP_g_num_t* globalNumGet(int id_block,int i_part) {
      CWP_g_num_t* gnum = _blockDB[id_block] -> GNumMeshGet(i_part);
      return gnum;
   }

   // inline Visu* getVisu();


   CWP_g_num_t* GNumEltsGet(int i_part);
   double* eltCentersGet(int i_part);
   void eltCentersCompute(int i_part);

   int* blockDBGet() {
     return _blocks_id;
   }

   int nBlockGet() {
     return _nBlocks;
   }

   CWP_Block_t blockTypeGet(int id_block) {
     return _blockDB[id_block] -> blockTypeGet();
   }

   CWP_Block_t* blocksTypeGet() {
     return _blocksType;
   }


   int* eltIdBlockGet(int i_part) {
     return _elt_id_block[i_part];
   }


   int* eltInBlockGet(int i_part) {
     return _elt_in_block[i_part];
   }


   CWP_g_num_t* gnumMeshBlockGet(int id_block,int i_part) {
     return _blockDB[id_block] -> GNumMeshGet(i_part);
   }

   PDM_geometry_kind_t geomKindGet() {
    return _geom_kind;
   }

   PDM_MPI_Comm _pdm_localComm;
  private:

    int _Mesh_nodal_block_std_type_size_get(CWP_Block_t type);

    const MPI_Comm                          &_localComm;              /*!< Communicator */

    int                                     _npart;                  /*!< Number of partition  */
    int                                     _nBlocks;                /*!< Number of blocks of the mesh */
    int*                                    _blocks_id;              /*!< List of block identifiers */
    CWP_Block_t*                            _blocksType;              /*!< Blocks type  */
    std::vector<cwipi::Block*>              _blockDB;                /*!< Blocks database  */

    std::vector<int>                        _nVertex;                /*!< Number of vertices for each partition  */
    std::vector<double*>                    _coords;                 /*!< Vertices coordinate for each partition  */

    std::vector<int>                        _nElts;                  /*!< Number of elements for each partition  */
    std::vector<CWP_g_num_t*>               _gnum_elt;
    std::vector<double*>                    _elt_centers;
    std::vector<int*>                       _elt_id_block;            /*!< Number of elements for each partition  */
    std::vector<int*>                       _elt_in_block;           /*!< Number of elements for each partition  */

    std::vector <CWP_g_num_t*>              _global_num_vtx;         /*!< Global vertices numbering for each partition  */
    std::vector <CWP_g_num_t*>              _global_num_elt;         /*!< Global elements numbering for each partition  */
    std::vector <CWP_g_num_t*>              _global_num_face;        /*!< Global faces numbering for each partition  */

    std::vector<int*>                       _connec_idx;
    std::vector<int*>                       _connec;

    // Visu                                   *_visu;                   /*!< Pointer to the Visu object */
    std::map<int,int>                       _id_visu;                /*!< Map of the PDM_Writer block identifier */
    CWP_Dynamic_mesh_t                      _displacement;           /*!< Type of mesh displacement */
    Coupling                               *_cpl;

    std::vector<int>                        _nCells;        
    std::vector<int*>                       _cellFaceIdx;
    std::vector<int*>                       _cellFace;

    std::vector<int>                        _nFace;
    std::vector<int*>                       _faceEdgeIdx;
    std::vector<int*>                       _faceEdge;
    std::vector<int*>                       _faceVtxIdx;
    std::vector<int*>                       _faceVtx;

    std::vector<int>                        _nEdge;
    std::vector<int*>                       _edgeVtx;

    std::vector<CWP_g_num_t*>               _faceLNToGN;
    std::vector<CWP_g_num_t*>               _cellLNToGN;

    int                                     _faceEdgeMethod;
    int                                     _faceVtxMethod;
    int                                     _cellFaceMethod;

    PDM_part_mesh_nodal_t                  *_pdmNodal_handle_index;  /*!< Mesh (nodal) index for paradigm handler */
    PDM_geometry_kind_t                     _geom_kind;

    bool                                    _isVtxGnumComputed;
    bool                                    _isEltGnumComputed;
    bool                                    _isFaceGnumComputed;



  //   Mesh &operator=(const Mesh &other);  /*!< Assigment operator not available */
  //   Mesh (const Mesh& other);            /*!< Copy constructor not available */

  };




  CWP_Dynamic_mesh_t Mesh::getDisplacement() {
    return _displacement;
  }

  // void Mesh::setVisu(Visu* visu) {
  //   // _visu = visu;
  // }

  // Visu* Mesh::getVisu() {
  //   // return _visu;
  // }

  int  Mesh::getIdVisu(int block_id) {
    return _id_visu[block_id];
  }

  int Mesh::getNPart() {
    return _npart;
  }

  double* Mesh::getCoordinates(int i_part) {
    return _coords[i_part];
  }

  PDM_part_mesh_nodal_t *Mesh::getPdmNodalIndex() {
    return _pdmNodal_handle_index;
  }

  MPI_Comm Mesh::getMPIComm() {
    return _localComm;
  }

  const MPI_Comm* Mesh::getMPICommP() {
    return &_localComm;
  }


  bool Mesh::coordsDefined () {
    for(int i=0; i<_npart;i++){
      if(_coords[i] == NULL)
        return false;
    }
    return true;
  }


  bool Mesh::gnumVtxRequired () {
    for(int i=0; i<_npart;i++){
      if(_global_num_vtx[i] == NULL && _nVertex[i] > 0)
        return true;
    }
    return false;
  }


  int Mesh::getPartNVertex(int i_part) const
  {
    return _nVertex[i_part];
  }

  double* Mesh::getVertexCoords(int i_part)
  {
    return _coords[i_part];
  }

  CWP_g_num_t* Mesh::getVertexGNum(int i_part)
  {
    return _global_num_vtx[i_part];
  }

  int Mesh::getBlockNElts(int id_block,int i_part)
  {
    return _blockDB[id_block] -> NEltsGet()[i_part];
  }

  /*int Mesh::getPartNPolyhedra(int i_part) const
  {
    return _nPolyhedra[i_part];
  }
  */

  int* Mesh::getPoly2DConnectivityIndex(int id_block,int i_part)
  {
    BlockFP * block = dynamic_cast <BlockFP *> (_blockDB[id_block]);
    return block-> ConnecIDXGet()[i_part];
  }

  int* Mesh::getPoly2DConnectivity(int id_block,int i_part)
  {
    BlockFP * block = dynamic_cast <BlockFP *> (_blockDB[id_block]);
    return block -> ConnecGet()[i_part];
  }

  int* Mesh::getStdConnectivity(int id_block,int i_part)
  {
    BlockStd * block = dynamic_cast <BlockStd *> (_blockDB[id_block]);
    return block-> ConnecGet()[i_part];
  }

  //-->>
  int Mesh::getNCell(int id_part) const {
      return _nCells[id_part];
  }

  int Mesh::getNFace(int id_part) const {
      return _nFace[id_part];
  }

  int Mesh::getNEdge(int id_part) const {
      return _nEdge[id_part];
  }

  int *Mesh::getCellFaceIndex(int i_part) {
      return _cellFaceIdx[i_part];
  }

  int *Mesh::getCellFace(int i_part) {
      return _cellFace[i_part];
  }

  int *Mesh::getFaceEdgeIndex(int i_part) {
      return _faceEdgeIdx[i_part];
  }

  int *Mesh::getFaceEdge(int i_part) {
      return _faceEdge[i_part];
  }

  int *Mesh::getFaceVtxIndex(int i_part) {
      return _faceVtxIdx[i_part];
  }

  int *Mesh::getFaceVtx(int i_part) {
      return _faceVtx[i_part];
  }

  int *Mesh::getEdgeVtx(int i_part) {
      return _edgeVtx[i_part];
  }
  //<<--

}

/**
 * \endcond
 */

#endif //__Mesh_H__
