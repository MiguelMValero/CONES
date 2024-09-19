#ifndef __COUPLING_I_H__
#define __COUPLING_I_H__
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
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "mesh.hxx"
#include "field.hxx"

namespace cwipi {


  /*----------------------------------------------------------------------------*
   * Methods about exchange frequency                                           *
   *----------------------------------------------------------------------------*/



    /**
     * \brief Return kind of displacement
     *  
     * 
     * \return kind of displacement
     *
     */

    CWP_Dynamic_mesh_t
    Coupling::DisplacementGet (
    )
    {
      return _displacement;
    }


  /**
   * \brief Setting receiving frequency.
   *
   * This function set receiving frequency. It must be used when
   * the type of receiving frequency is \ref CWP_TIME_EXCH_N_TIME_STEP
   *
   * \param [in]  n_step     Frequency in steps number
   *
   */

  void
  Coupling::recvFreqSet (
    int n_step
  )
  {
    PDM_UNUSED (n_step);
    PDM_error(__FILE__, __LINE__, 0, "recvFreqSet not implemented yet\n");
  }



  /*----------------------------------------------------------------------------*
   * Methods about visualization                                                *
   *----------------------------------------------------------------------------*/


  /**
   *
   * \brief Return the Visu object associated to this coupling
   *
   * \return Visu object pointer
   *
   */

  // Visu* Coupling::visuGet() {
  //    // return &_visu;
  // }



  /**
   *
   * \brief Return the PDM_writer object associated to this coupling
   *
   * \return Visu object pointer
   *
   */

  PDM_writer_t* 
  Coupling::writerGet ()
  {
    return _writer;
  }



  /**
   *
   * \brief Return the PDM_writer object associated to this coupling
   *
   * \return Visu object pointer
   *
   */

  int 
  Coupling::freqWriterGet ()
  {
    return _freq_writer;
  }


  /**
   * \brief Curent number of coupling step
   *  
   * 
   * \return the current number of coupling step
   *
   */

  int 
  Coupling::NStepGet (
  )
  {
    return _n_step;
  }

  /*----------------------------------------------------------------------------*
   * Methods about mesh                                                     *
   *----------------------------------------------------------------------------*/

  bool
  Coupling::has_mesh
  (
   )
  {
    int i_rank;
    MPI_Comm_rank(_communication.unionCommGet(), &i_rank);

    return ((commTypeGet() == CWP_COMM_PAR_WITH_PART                ) ||
            (commTypeGet() == CWP_COMM_PAR_WITHOUT_PART &&
             i_rank == _communication.unionCommLocCodeRootRanksGet()));
  }

  /**
   * \brief Setting vertices
   *
   * This method set partition vertices
   *
   * \param [in]  i_part      Current partition
   * \param [in]  n_pts       Number of points
   * \param [in]  coord       Coordinates (size = 3 * n_pts)
   * \param [in]  global_num  Pointer to global element number (or NULL)
   *
   */

  void 
  Coupling::meshVtcsSet (
    const int          i_part,
    const int          n_pts,
    double             coords[],
    CWP_g_num_t        global_num[]
  )
  {
    if (has_mesh()) {
      _mesh.coordSet(i_part,
                     n_pts,
                     coords,
                     global_num);
    }
    else {
      _mesh.coordSet(i_part,
                     0,
                     NULL,
                     NULL);
    }
  }


 /**
  * \brief Add a block to the interface mesh.
  *
  *
  * \param [in]  block_type       Block type
  *
  * \return block identifier
  */

  int 
  Coupling::meshBlockAdd (
    const CWP_Block_t     block_type
  )
  {
    return _mesh.blockAdd(block_type);
  }


  /**
   * \brief Set a standard block to the interface mesh
   *
   * This function adds a connectivity block to the geometric support.
   *
   *  Definition of element connectivity is :
   *
   *  - edge (\ref CWP_BLOCK_EDGE2) :
   *
   *   \code
   *       1 x-------x 2
   *   \endcode
   *
   *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
   *
   *   \code
   *       1 x-------x 3
   *          \     /
   *           \   /
   *            \ /
   *             x 2
   *   \endcode
   *
   *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
   *
   *   \code
   *          4 x-------x 3
   *           /       /
   *          /       /
   *       1 x-------x2
   *   \endcode
   *
   *   - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) :
   *
   *   \code
   *             x 4
   *            /|\
   *           / | \
   *          /  |  \
   *       1 x- -|- -x 3
   *          \  |  /
   *           \ | /
   *            \|/
   *             x 2
   *   \endcode
   *
   *   - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
   *
   *   \code
   *              5 x
   *               /|\
   *              //| \
   *             // |  \
   *          4 x/--|---x 3
   *           //   |  /
   *          //    | /
   *       1 x-------x 2
   *   \endcode
   *
   *  - prism (\ref CWP_BLOCK_CELL_PRISM6) :
   *
   *   \code
   *       4 x-------x 6
   *         |\     /|
   *         | \   / |
   *       1 x- \-/ -x 3
   *          \ 5x  /
   *           \ | /
   *            \|/
   *             x 2
   *   \endcode
   *
   *  -  hexaedron (\ref CWP_BLOCK_CELL_HEXA8) :
   *
   *   \code
   *          8 x-------x 7
   *           /|      /|
   *          / |     / |
   *       5 x-------x6 |
   *         | 4x----|--x 3
   *         | /     | /
   *         |/      |/
   *       1 x-------x 2
   *   \endcode
   *
   * \param [in]  i_part      Current partition
   * \param [in]  block_id    Block identifier
   * \param [in]  n_elts      Number of elements
   * \param [in]  connec      Connectivity (size = n_vertex_elt * n_elts)
   * \param [in]  global_num  Pointer to global element numbering (or NULL)
   *
   */

  void 
  Coupling::meshStdBlockSet (
    const int           i_part,
    const int           block_id,
    const int           n_elts,
    int                 connec[],
    CWP_g_num_t         global_num[]
  )
  {
    _mesh.stdBlockSet (i_part,
                       block_id,
                       n_elts,
                       connec,
                       global_num);
  }


  void
  Coupling::meshHOBlockSet (
    const int           i_part,
    const int           block_id,
    const int           n_elts,
    int                 connec[],
    CWP_g_num_t         global_num[],
    const int           order,
    const char         *ho_ordering
  )
  {
    _mesh.HOBlockSet(i_part,
                     block_id,
                     n_elts,
                     connec,
                     global_num,
                     order,
                     ho_ordering);
  }

  /**
   * \brief Get a standard block to the interface mesh
   *
   *
   * \param [in]   i_part      Current partition
   * \param [in]   block_id    Block identifier
   * \param [out]  n_elts      Number of elements
   * \param [out]  connec      Connectivity
   * \param [out]  global_num  Pointer to global element number (or NULL)
   *
   */

  void
  Coupling::meshStdBlockGet (
    const int    i_part,
    const int    block_id,
    int         *n_elts,
    int         **connec,
    CWP_g_num_t **global_num
  )
  {
    _mesh.stdBlockGet (i_part,
                       block_id,
                       n_elts,
                       connec,
                       global_num);
  }


  /**
   * \brief Get a standard block to the interface mesh
   *
   *
   * \param [in]   i_part      Current partition
   * \param [in]   block_id    Block identifier
   * \param [out]  n_elts      Number of elements
   * \param [out]  order       Element order
   * \param [out]  connec      Connectivity
   * \param [out]  global_num  Pointer to global element number (or NULL)
   *
   */

  void
  Coupling::meshHOBlockGet (
    const int     i_part,
    const int     block_id,
    int          *n_elts,
    int          *order,
    int         **connec,
    CWP_g_num_t **global_num
  )
  {
    _mesh.HOBlockGet(i_part,
                     block_id,
                     n_elts,
                     order,
                     connec,
                     global_num);
  }

  /**
    * \brief Get the standard block type
    *
    * \param [in]  block_id    Block identifier
    *
    * \return block type
    *
    */

  CWP_Block_t
  Coupling::meshStdBlockTypeGet (
    const int           block_id
  )
  {
    CWP_Block_t block_type = _mesh.stdBlockTypeGet(block_id);

    return block_type;
  }

  /**
   * \brief Set a generic high order block to the interface mesh
   *
   *
   * \param [in]  i_part      Partition identifier
   * \param [in]  block_id    Block identifier
   * \param [in]  n_elts      Number of elements
   * \param [in]  order       SpatialInterp order
   * \param [in]  connec      Connectivity (size = n_vertex_elt * n_elts)
   * \param [in]  global_num  Pointer to global element number (or NULL)
   *
   */

  void
  Coupling::meshHighOrderBlockSet (
    const int           i_part,
    const int           block_id,
    const int           n_elts,
    const int           order,
    int                 connec[],
    CWP_g_num_t         global_num[]
  )
  {
    PDM_UNUSED (i_part);
    PDM_UNUSED (block_id);
    PDM_UNUSED (n_elts);
    PDM_UNUSED (order);
    PDM_UNUSED (connec);
    PDM_UNUSED (global_num);

    PDM_error(__FILE__, __LINE__, 0, "meshHighOrderBlockSet not implemented yet\n"); 
  }


  /**
   * \brief Set the connectivity of a polygon block in a mesh interface partition.
   *
   *
   * \param [in]  i_part      Current partition
   * \param [in]  block_id    Block identifier
   * \param [in]  n_elts      Number of elements
   * \param [in]  connec_idx  Connectivity index (connec_idx[0] = 0 and
   *                          size = n_elts + 1)
   * \param [in]  connec      Connectivity (size = connec_idx[n_elts] * n_elts)
   * \param [in]  global_num  Pointer to global element number (or NULL)
   *
   */

  void 
  Coupling::meshFPolyBlockSet (
    const int            i_part,
    const int            block_id,
    const int            n_elts,
    int                  connec_idx[],
    int                  connec[],
    CWP_g_num_t          global_num[]
  )
  {
    if (has_mesh()) {
      _mesh.poly2DBlockSet(i_part,
                           block_id,
                           n_elts,
                           connec_idx,
                           connec,
                           global_num);
    }
    else {
      _mesh.poly2DBlockSet(i_part,
                           block_id,
                           0,
                           NULL,
                           NULL,
                           NULL);
    }
  }


  /**
   * \brief Set the connectivity of a polygon block in a mesh interface partition.
   *
   *
   * \param [in]   i_part      Current partition
   * \param [in]   block_id    Block identifier
   * \param [out]  n_elts      Number of elements
   * \param [out]  connec_idx  Connectivity index (connec_idx[0] = 0 and
   *                          size = n_elts + 1)
   * \param [out]  connec      Connectivity (size = connec_idx[n_elts] * n_elts)
   * \param [out]  global_num  Pointer to global element number (or NULL)
   *
   */

  void
  Coupling::meshFPolyBlockGet (
    const int            i_part,
    const int            block_id,
    int                 *n_elts,
    int                **connec_idx,
    int                **connec,
    CWP_g_num_t        **global_num
  )
  {
    _mesh.poly2DBlockGet(i_part,
                         block_id,
                         n_elts,
                         connec_idx,
                         connec,
                         global_num);  
  }


  /**
   * \brief Set the connectivity of a polyhedron block in a mesh interface partition.
   *
   * Definition of element connectivity is :
   *
   * \param [in]  i_part            Current partition
   * \param [in]  block_id          Block identifier
   * \param [in]  n_elts            Number of elements
   * \param [in]  connec_cells_idx  Polyhedron to face index
   *                                (src_poly_cell_face_idx[0] = 0 and
   *                                 size = n_elts + 1)
   * \param [in]  connec_cells      Polyhedron to face connectivity
   *                                (size = cell_face_idx[n_elts])
   * \param [in]  n_faces           Number of faces
   * \param [in]  connec_faces_idx  Polyhedron face to vertex index
   *                                (face_vertex_idx[0] = 0 and
   *                                size_idx = max(cell_face_connec) + 1)
   * \param [in]  connec_faces      Polyhedron face to vertex connectivity
   *                                (size = face_vertex_idx[size_idx - 1])
   * \param [in]  global_num        Pointer to global element number (or NULL)
   *
   */

  void 
  Coupling::meshCPolyBlockSet (
    const int           i_part,
    const int           block_id,
    const int           n_elts,
    const int           n_faces,
          int           connec_faces_idx[],
          int           connec_faces[],
          int           connec_cells_idx[],
          int           connec_cells[],
          CWP_g_num_t   global_num[]
  )
  {
    if (has_mesh()) {
      _mesh.poly3DBlockSet(i_part,
                           block_id,
                           n_elts,
                           n_faces,
                           connec_faces_idx,
                           connec_faces,
                           connec_cells_idx,
                           connec_cells,
                           global_num);
    }
    else {
      _mesh.poly3DBlockSet(i_part,
                           block_id,
                           0,
                           0,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL);
    }
  }


  /**
   * \brief Set the connectivity of a polyhedron block in a mesh interface partition.
   *
   * Definition of element connectivity is :
   *
   * \param [in]  i_part            Current partition
   * \param [in]  block_id          Block identifier
   * \param [out]  n_elts            Number of elements
   * \param [out]  connec_cells_idx  Polyhedron to face index
   *                                (src_poly_cell_face_idx[0] = 0 and
   *                                 size = n_elts + 1)
   * \param [out]  connec_cells      Polyhedron to face connectivity
   *                                (size = cell_face_idx[n_elts])
   * \param [out]  n_faces           Number of faces
   * \param [out]  connec_faces_idx  Polyhedron face to vertex index
   *                                (face_vertex_idx[0] = 0 and
   *                                size_idx = max(cell_face_connec) + 1)
   * \param [out]  connec_faces      Polyhedron face to vertex connectivity
   *                                (size = face_vertex_idx[size_idx - 1])
   * \param [out]  global_num        Pointer to global element number (or NULL)
   *
   */

  inline void
  Coupling::meshCPolyBlockGet (
    const int           i_part,
    const int           block_id,
    int                *n_elts,
    int                *n_faces,
    int               **connec_faces_idx,
    int               **connec_faces,
    int               **connec_cells_idx,
    int               **connec_cells,
    CWP_g_num_t       **global_num
  )
  {
    _mesh.poly3DBlockGet(i_part,
                         block_id,
                         n_elts,
                         n_faces,
                         connec_faces_idx,
                         connec_faces,
                         connec_cells_idx,
                         connec_cells,
                         global_num);  
  }


  /**
   * \brief Adding a polyhedron block to the mesh from
   * a face-to-cell connectivity and a vertices-to-faces connectivity.
   *
   * This function adds a polyhedron 3D block to the mesh from
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

  void 
  Coupling::meshFromCellFaceSet(
    const int         i_part,
    const int         n_cells,
          int         cell_face_idx[],
          int         cell_face[],
          int         n_faces,
          int         face_vtx_idx[],
          int         face_vtx[],
          CWP_g_num_t parent_num[]
  )
  {
    if (has_mesh()) {
      _mesh.fromCellFaceSet(i_part,
                            n_cells,
                            cell_face_idx,
                            cell_face,
                            n_faces,
                            face_vtx_idx,
                            face_vtx,
                            parent_num);
    }
    else {
      _mesh.fromCellFaceSet(i_part,
                            0,
                            NULL,
                            NULL,
                            0,
                            NULL,
                            NULL,
                            NULL);
    }
  }


  /**
   * \brief Adding a polygon 2D block to the mesh from
   * a vertices-to-faces connectivity and a edge-to-face connectivity.
   *
   * This function add a polygon 2D block to the mesh from
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
   * \param [in]  edge_vtx_idx      Vertices to edges connectivity index
   *                                (edge_vtx_idx[0] = 0 and
   *                                size_idx = max(edge_vtx) + 1)
   * \param [in]  edge_vtx          Polygon vertices to edges connectivity
   *                                (size = edge_vtx_idx[size_idx - 1])
   * \param [in]  parent_num        Pointer to parent element number (or NULL)
   *
   */

  void 
  Coupling::meshFromFacesEdgeSet (
    const int         i_part,
    const int         n_faces,
          int         face_edge_idx[],
          int         face_edge[],
    const int         n_edges,
          int         edge_vtx[],
          CWP_g_num_t parent_num[]
  )
  {
    if (has_mesh()) {
      _mesh.fromFacesEdgeSet (i_part,
                              n_faces,
                              face_edge_idx,
                              face_edge,
                              n_edges,
                              edge_vtx,
                              parent_num);
    }
    else {
      _mesh.fromFacesEdgeSet (i_part,
                              0,
                              NULL,
                              NULL,
                              0,
                              NULL,
                              NULL);
    }
  }


  void
  Coupling::meshFromFacesVtxSet (
    const int         i_part,
    const int         n_faces,
          int         face_vtx_idx[],
          int         face_vtx[],
          CWP_g_num_t global_num[]
  )
  {
    if (has_mesh()) {
      _mesh.fromFacesVtxSet (i_part,
                             n_faces,
                             face_vtx_idx,
                             face_vtx,
                             global_num);
    }
    else {
      _mesh.fromFacesVtxSet (i_part,
                             0,
                             NULL,
                             NULL,
                             NULL);
    }
  }


  /**
   * \brief SpatialInterp mesh removal
   *
   * This function delete the  mesh
   *
   */

  void 
  Coupling::meshDel()
  {
    if (_writer != NULL) {
      // free geometric variable data
      PDM_writer_geom_data_free(_writer,
                                _id_geom_writer);
      if (_id_field_partitioning_writer >= 0) PDM_writer_var_data_free(_writer,
                                                                        _id_field_partitioning_writer);
      if (_id_field_ranking_writer >= 0) PDM_writer_var_data_free(_writer,
                                                                   _id_field_ranking_writer);
      // free user target geometric variable data
      if (_userTargetN != NULL) {
        PDM_writer_geom_data_free(_writer,
                                  _id_user_tgt_geom_writer);
        PDM_writer_var_data_free(_writer,
                                  _id_user_tgt_field_partitioning_writer);
        PDM_writer_var_data_free(_writer,
                                  _id_user_tgt_field_ranking_writer);
      }
      // free field variables data
      std::map < string, Field * >::iterator itf = _fields.begin();
      while (itf != _fields.end()) {

        if (itf->second->visuStatusGet() == CWP_STATUS_ON) {

          for (int i = 0; i < itf->second->nComponentGet(); i++) {
            int id_send = itf->second->_id_writer_var_send_get()[i];
            int id_recv = itf->second->_id_writer_var_recv_get()[i];

            if (id_send >= 0) PDM_writer_var_data_free(_writer, id_send);
            if (id_recv >= 0) PDM_writer_var_data_free(_writer, id_recv);
          }

          int id_send_status = itf->second->_id_writer_var_send_status_get();
          if (id_send_status >= 0) PDM_writer_var_data_free(_writer, id_send_status);
          int id_recv_status = itf->second->_id_writer_var_recv_status_get();
          if (id_recv_status >= 0) PDM_writer_var_data_free(_writer, id_recv_status);

        } // end if field is visualized

        itf++;
      } // end iterate over fields
    } // end if write is done

    _mesh.meshDel();

    _is_mesh_finalized = 0;
  }


  /**
   *
   * \brief Finalize mesh description (Computation of entities global numbers if not given by the user)
   *
   */

  void Coupling::meshFinalize() 
  {
    _is_mesh_finalized = 1;
    _mesh.geomFinalize();
  }


  /**
   *
   * \brief Return the mesh object associated to this coupling
   *
   * \return Mesh object pointer
   *
   */

  Mesh* Coupling::meshGet() {
     return &_mesh;
  }


  /*----------------------------------------------------------------------------*
   * Methods about field                                                        *
   *----------------------------------------------------------------------------*/


  /**
   *
   * \brief Get nunmber of field components
   *
   * \param [in]   field_id       Field identifier
   *
   * \return                      number of field components
   *
   */

  int
  Coupling::fieldNComponentGet (
    const string &field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
       PDM_error(__FILE__, __LINE__, 0, "'%s' not existing field\n", field_id.c_str());
    }
    return It->second->nComponentGet();
  }


  /**
   *
   * \brief Get location of the degrees of freedom
   *
   * \param [in]   field_id       Field identifier
   *
   * \return                      Field data type
   *
   */

  CWP_Dof_location_t Coupling::fieldDofLOcationGet
  (
    const string &field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                 "'%s' not existing field\n", field_id.c_str());
    }
    return It->second->locationGet();

  }


  /**
   *
   * \brief Get field storage type
   *
   * \param [in]   field_id       Field identifier
   *
   */

  CWP_Field_storage_t
  Coupling::fieldStorageGet (
    const string &field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                 "'%s' not existing field\n", field_id.c_str());
    }

    return It->second->storageTypeGet();

  }


  /**
   *
   * \brief Removing a field
   *
   * \param [in]  field_id       Field identifier
   *
   */

  void
  Coupling::fieldDel (
    const string &field_id
  )
  {

    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      delete It->second;
      It->second = NULL;
      _fields.erase(It->first);
    }

  }

  /*----------------------------------------------------------------------------*
   * Methods about exchange                                                     *
   *----------------------------------------------------------------------------*/

  void
  Coupling::computedTargetsBcastEnable
  (
    const string &field_id
  )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    if (it == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Error computedTargetsBcastEnable : '%s' not existing field\n", field_id.c_str());
    }

    Field* field = it->second;

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_SEND) {
      PDM_error(__FILE__, __LINE__, 0, "Error computedTargetsBcastEnable : '%s' does not receive data\n", field_id.c_str());
    }

    field->computedTgtBcastEnable();
  }

  /**
   *
   * \brief Return the number of uncomputed targets
   *
   * \return                Number of uncomputed targets
   */

  int 
  Coupling::nUncomputedTargetsGet
  (
    const string &field_id,
    const int  i_part
  )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    if (it == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Error nUncomputedTargetsGet : '%s' not existing field\n", field_id.c_str());
    }

    Field* field = it->second;

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_SEND) {
      PDM_error(__FILE__, __LINE__, 0, "Error nUncomputedTargetsGet : '%s' does not receive data\n", field_id.c_str());     
    }

    if (!has_mesh() && !field->computedTgtBcastIsEnabled()) {
      PDM_error(__FILE__, __LINE__, 0, "Error nUncomputedTargetsGet : CWP_Computed_tgts_bcast_enable must be called for field '%s'\n", field_id.c_str());
    }

    return _spatial_interp_recv[make_pair(field->locationGet(), field->linkedFieldLocationGet())]->nUncomputedTargetsGet(i_part);

  }



  /**
   *
   * \brief Return uncomputed targets
   *
   * \return                Uncomputed targets
   */

  const int *
  Coupling::uncomputedTargetsGet (
    const string &field_id,
    const int  i_part
  )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    if (it == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Error uncomputedTargetsGet : '%s' not existing field\n", field_id.c_str());
    }

    Field* field = it->second;

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_SEND) {
      PDM_error(__FILE__, __LINE__, 0, "Error unncomputedTargetsGet : '%s' does not receive data\n", field_id.c_str());     
    }

    if (!has_mesh() && !field->computedTgtBcastIsEnabled()) {
      PDM_error(__FILE__, __LINE__, 0, "Error uncomputedTargetsGet : CWP_Computed_tgts_bcast_enable must be called for field '%s'\n", field_id.c_str());
    }

    return _spatial_interp_recv[make_pair(field->locationGet(), field->linkedFieldLocationGet())]->uncomputedTargetsGet(i_part);
  }

  /**
   *
   * \brief Return the number of computed targets
   *
   * \return                Number of computed targets
   */

  int
  Coupling::nComputedTargetsGet (
    const string &field_id,
    const int  i_part
  )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    if (it == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Error nComputedTargetsGet : '%s' not existing field\n", field_id.c_str());
    }

    Field* field = it->second;

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_SEND) {
      PDM_error(__FILE__, __LINE__, 0, "Error nComputedTargetsGet : '%s' does not receive data\n", field_id.c_str());     
    }

    if (!has_mesh() && !field->computedTgtBcastIsEnabled()) {
      PDM_error(__FILE__, __LINE__, 0, "Error nComputedTargetsGet : CWP_Computed_tgts_bcast_enable must be called for field '%s'\n", field_id.c_str());
    }

    return _spatial_interp_recv[make_pair(field->locationGet(), field->linkedFieldLocationGet())]->nComputedTargetsGet(i_part);
  }

  /**
   *
   * \brief Return computed targets
   *
   * \return                Computed targets
   */

  const int *
  Coupling::computedTargetsGet (
    const string &field_id,
    const int  i_part
  )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    if (it == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Error computedTargetsGet : '%s' not existing field\n", field_id.c_str());
    }

    Field* field = it->second;

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_SEND) {
      PDM_error(__FILE__, __LINE__, 0, "Error computedTargetsGet : '%s' does not receive data\n", field_id.c_str());     
    }

    if (!has_mesh() && !field->computedTgtBcastIsEnabled()) {
      PDM_error(__FILE__, __LINE__, 0, "Error computedTargetsGet : CWP_Computed_tgts_bcast_enable must be called for field '%s'\n", field_id.c_str());
    }

    return _spatial_interp_recv[make_pair(field->locationGet(), field->linkedFieldLocationGet())]->computedTargetsGet(i_part);
  }



  void
  Coupling::involvedSourcesBcastEnable
  (
    const string &field_id
  )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    if (it == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Error involvedSourcesBcastEnable : '%s' not existing field\n", field_id.c_str());
    }

    Field* field = it->second;

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_RECV) {
      PDM_error(__FILE__, __LINE__, 0, "Error involvedSourcesBcastEnable : '%s' does not receive data\n", field_id.c_str());
    }

    field->computedTgtBcastEnable();
  }


  int
  Coupling::nInvolvedSourcesGet(
    const string &field_id,
    const int  i_part
  )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    if (it == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Error nInvolvedSourcesGet : '%s' not existing field\n", field_id.c_str());
    }

    Field* field = it->second;

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_RECV) {
      PDM_error(__FILE__, __LINE__, 0, "Error nInvolvedSourcesGet : '%s' does not send data\n", field_id.c_str());
    }

    if (!has_mesh() && !field->involvedSrcBcastIsEnabled()) {
      PDM_error(__FILE__, __LINE__, 0, "Error nInvolvedSourcesGet : CWP_Involved_srcs_bcast_enable must be called for field '%s'\n", field_id.c_str());
    }

    return _spatial_interp_send[make_pair(field->locationGet(), field->linkedFieldLocationGet())]->nInvolvedSourcesGet(i_part);
  }

  const int *
  Coupling::involvedSourcesGet(
    const string &field_id,
    const int  i_part
  )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    if (it == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Error involvedSourcesGet : '%s' not existing field\n", field_id.c_str());
    }

    Field* field = it->second;

    if (field->exchangeTypeGet() == CWP_FIELD_EXCH_RECV) {
      PDM_error(__FILE__, __LINE__, 0, "Error involvedSourcesGet : '%s' does not send data\n", field_id.c_str());
    }

    if (!has_mesh() && !field->involvedSrcBcastIsEnabled()) {
      PDM_error(__FILE__, __LINE__, 0, "Error involvedSourcesGet : CWP_Involved_srcs_bcast_enable must be called for field '%s'\n", field_id.c_str());
    }

    return _spatial_interp_send[make_pair(field->locationGet(), field->linkedFieldLocationGet())]->involvedSourcesGet(i_part);
  }

  /*----------------------------------------------------------------------------*
   * methods about user interpolation                                           *
   *----------------------------------------------------------------------------*/

  /**
   *
   * \brief Setting of an user interpolation function.
   *
   * This function takes into account an user interpolation function written with
   *  void (* \ref CWP_Interp_function_t) interface.
   *
   * \param [in] fct        Function
   *
   */

  void 
  Coupling::interpFunctionSet (
    const string                     field_id,
          CWP_Interp_function_t      fct
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      It -> second -> interpFunctionSet(fct);
    }
  }

  void
  Coupling::interpFunctionFSet (
    const string                     field_id,
          CWP_Interp_function_t      fct
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      It -> second -> interpFunctionFSet(fct);
    }
  }

  void
  Coupling::interpFunctionPSet (
    const string                     field_id,
          CWP_Interp_function_p_t    fct
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      It -> second -> interpFunctionPSet(fct);
    }
  }

  void
  Coupling::interpFunctionPUnset (
    const string                     field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      It -> second -> interpFunctionPUnset();
    }
  }


  /**
   *
   * \brief Unsetting of an user interpolation function.
   *
   */

  inline void
  Coupling::interpFunctionUnSet (
    const string field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      It -> second -> interpFunctionUnSet();
    }
  }

  inline void
  Coupling::interpFunctionFUnSet (
    const string field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      It -> second -> interpFunctionFUnSet();
    }
  }

  /*----------------------------------------------------------------------------*
   * methods about attributes                                                   *
   *----------------------------------------------------------------------------*/

  /**
   *
   * \brief Return communication type
   *
   * \return CWP_Comm_t Communication Type
   *
   */

  CWP_Comm_t
  Coupling::commTypeGet (
  )
  {
    return _commType;
  }


  /**
   *
   * \brief Return the dimesnion of geometric entities of this coupling
   *
   * \return Entities dimension 
   *
   */

  CWP_Interface_t 
  Coupling::entitiesDimGet() 
  {
     return _entities_dim;
  }


  /**
   *
   * \brief Return the local code fields defined for this coupling
   *
   * \return Fields 
   *
   */

  std::map < string, Field * >* 
  Coupling::fieldsGet() 
  {
     return &_fields;
  }


  /**
   *
   * \brief Return the spatial interpolation objects according to direction of the exchange
   * 
   * \param [in] exchDirection     Direction of the exchange
   *
   * \return Local  spatial interpolation objects
   *
   */

  std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > ,SpatialInterp*>* 
  Coupling::spatialInterpGet(
    CWP_Field_exch_t exchDirection
  )
  {
    if (exchDirection == CWP_FIELD_EXCH_SEND) {
      return &_spatial_interp_send;
    }
    else if (exchDirection == CWP_FIELD_EXCH_RECV) {
      return &_spatial_interp_recv;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "SpatialInterp not found.\n");
      return nullptr;     
    } 

  }

  /**
   *
   * \brief Return the spatial interpolation according to parameters
   * 
   * \param [in] localLocation     Local location of degrees of freedom
   * \param [in] cplLocation       Coupled location of degrees of freedom
   * \param [in] exchDirection     Direction of the exchange
   *
   * \return Local code properties
   *
   */

  SpatialInterp* 
  Coupling::spatialInterpGet (
    CWP_Dof_location_t localLocation, 
    CWP_Dof_location_t cplLocation, 
    CWP_Field_exch_t exch_t) 
  {

    if (exch_t == CWP_FIELD_EXCH_SEND) {
      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > ,SpatialInterp*> ::iterator p;
      p = _spatial_interp_send.find(make_pair(localLocation, cplLocation));
      if (p == _spatial_interp_send.end())
        PDM_error(__FILE__, __LINE__, 0, "SpatialInterp not found.\n");
      return p->second;
    }

    else if (exch_t == CWP_FIELD_EXCH_RECV) {
      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > ,SpatialInterp*> ::iterator p;
      p = _spatial_interp_recv.find(make_pair(localLocation, cplLocation));
      if (p == _spatial_interp_recv.end())
        PDM_error(__FILE__, __LINE__, 0, "SpatialInterp not found.\n");
      return p->second;
    }

    else {
      PDM_error(__FILE__, __LINE__, 0, "SpatialInterp not found.\n");
      return nullptr;     
    }

  }


  /**
   *
   * \brief Return the local code properties
   *
   * \return Local code properties
   *
   */


  CodeProperties* 
  Coupling::localCodePropertiesGet() 
  {
    return const_cast<CodeProperties*>(&_localCodeProperties);
  }


  /**
   *
   * \brief Return the coupled code properties
   *
   * \return Local code properties
   *
   */

  CodeProperties* 
  Coupling::coupledCodePropertiesGet() 
  {
    return const_cast<CodeProperties*>(&_coupledCodeProperties);
  }


  /**
   *
   * \brief Return the communication way 
   *
   * \return Local code properties
   *
   */

  Communication* 
  Coupling::communicationGet() {
    return const_cast<Communication*>(&_communication);
  }


  /**
   *
   * \brief Return the couplings data base
   *
   * \return Coupling data base
   *  
   */

  CouplingDB* 
  Coupling::couplingDBGet() {
    return &_cplDB;
  }

    
  /**
   *
   * \brief Return the coupling id
   *
   * \return id
   *  
   */

  string 
  Coupling::IdGet(){
    return _cplId;
  }



  /**
   *
   * \brief Return sendSpatial map
   *
   * \return id
   *  
   */

  std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &    
  Coupling::sendSpatialInterpGet() {
    return _spatial_interp_send;
  }


  /**
   *
   * \brief Return recvSpatial map
   *
   * \return id
   *  
   */

  std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &    
  Coupling::recvSpatialInterpGet(){
    return _spatial_interp_recv;
  }

  /*----------------------------------------------------------------------------*
   * methods about user target                                                  *
   *----------------------------------------------------------------------------*/

  
  /**
   * \brief Setting user target points
   *
   * This function must be called if the nature of receiving fieldsDouble
   * is \ref CWP_DOF_LOCATION_USER
   *
   * \param [in]  i_part  Current partition
   * \param [in]  n_pts   Number of points
   * \param [in]  coord   Coordinates (size = 3 * n_pts)
   * \param [in]  g_num   global number or NUL (size = n_pts)
   *
   */

  void
  Coupling::userTargetSet (
    const int         iPart,
    const int         nPts,
    const double      coords[],
    const CWP_g_num_t gNum[]
  )
  {
    if (_userTargetN == nullptr) {
      _userTargetN =  (int *) malloc (sizeof(int) * (_nPart));
      for (int iPart1 = 0; iPart1 < _nPart; iPart1++) {
        _userTargetN[iPart1] = 0;
      }
      if (gNum != nullptr) {
        _userTargetGnum = (const PDM_g_num_t**) malloc (sizeof(PDM_g_num_t *) * _nPart); 
        for (int iPart1 = 0; iPart1 < _nPart; iPart1++) {
          _userTargetGnum[iPart1] = nullptr;
        }
      }
      _userTargetCoord = (const double**) malloc (sizeof(double *) * _nPart); 
      for (int iPart1 = 0; iPart1 < _nPart; iPart1++) {
        _userTargetCoord[iPart1] = nullptr;
      }
    }

    _userTargetN[iPart]     = nPts;
    _userTargetCoord[iPart] = coords;

    if (gNum != nullptr) {
      _userTargetGnum[iPart] = gNum;
    }
  }


  /**
   * \brief Return global number of user targets 
   * 
   * \param [in]  iPart     Current partition
   * 
   * 
   * \return Global Number 
   *
   */

  const CWP_g_num_t *
  Coupling::userTargetGNumGet (
    const int         iPart
  )
  {
    if (! _userTargetGnum) {
      userTargetGnumCompute();
    }
    return _userTargetGnum[iPart];
  }


  /**
   * \brief Return coords of user targets 
   * 
   * \param [in]  iPart     Current partition
   * 
   * 
   * \return Coordinates 
   *
   */

  const double *
  Coupling::userTargetCoordsGet (
    const int         iPart
  ) const
  {
    assert (_userTargetCoord != nullptr);
    return _userTargetCoord[iPart]; 
  }


  /**
   * \brief Return number of user targets 
   * 
   * \param [in]  iPart     Current partition
   * 
   * 
   * \return Number of user targets
   *
   */

  int
  Coupling::userTargetNGet (
    const int         iPart
  )
  {
    assert (_userTargetN != nullptr);
    return _userTargetN[iPart];   
  }

  /**
   * \brief Return number of partition 
   *  
   * 
   * \return Number of partition
   *
   */
  
  int
  Coupling::nPartGet (
  ) const
  {
    return _nPart;
  }


  /**
   * \brief Return number of partition of coupled code
   *  
   * 
   * \return Number of partition
   *
   */
  
  int
  Coupling::cplNPartGet (
  )
  {
    return _cplNPart;
  }


  /**
   * \brief Return the user spatial interpolation properties of type double
   *
   * \return Map storing the user spatial interpolation properties of type double
   *
   */
  std::map <std::string, double> &
  Coupling::SpatialInterpPropertiesDoubleGet (
  )
  {
    return _spatial_interp_properties_double;
  }

  /**
   * \brief Return the user spatial interpolation properties of type int
   *
   * \return Map storing the user spatial interpolation properties of type int
   *
   */
  std::map <std::string, int> &
  Coupling::SpatialInterpPropertiesIntGet (
  )
  {
    return _spatial_interp_properties_int;
  }


  int
  Coupling::idGeomWriterGet
  (
   CWP_Dof_location_t dof_location
  )
  {
    if (dof_location == CWP_DOF_LOCATION_USER) {
      return _id_user_tgt_geom_writer;
    } else {
      return _id_geom_writer;
    }
  }

} // name space cwipi


#endif //__COUPLING_PROPERTIES_I_H__
