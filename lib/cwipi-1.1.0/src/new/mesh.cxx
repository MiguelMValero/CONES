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
#include <mesh.hxx>
#include <mpi.h>
#include <cstdlib>

#include <pdm_mesh_nodal.h>
#include <pdm_gnum.h>
#include <pdm_error.h>
#include <pdm_printf.h>
#include <pdm_logging.h>
#include <pdm_array.h>
#include "cwp.h"
#include "factory.hpp"
#include "block.hxx"
#include "blockStd.hxx"
#include "blockCP.hxx"
#include "blockFP.hxx"
#include "blockHO.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

/**
 * \cond
 */

namespace cwipi {

  /**
   * \typedef FB
   *
   * \brief Block Factory
   *
   *  A block \ref Factory wich makes \ref Block
   *  class objects.
   *  The type of block objects build depends on the
   *  block type \ref CWP_Block_t .
   *
   */

  typedef Factory<Block, CWP_Block_t> FB;





  static void _gen_face_gnum
  (
   bool                       with_edges,
   PDM_MPI_Comm               comm,
   int                        n_part,
   std::vector<int>     const &n_face,
   std::vector<int*>    const &face_vtx_idx,
   std::vector<int*>    const &face_vtx,
   std::vector<int*>    const &face_edge,
   std::vector<int*>    const &edge_vtx,
   std::vector<double*> const &vtx_coord,
   std::vector<CWP_g_num_t*>  &face_ln_to_gn
   )
  {
    PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,
                                               n_part,
                                               PDM_FALSE,
                                               1e-3,
                                               comm,
                                               PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

    /* Compute face centers */
    double **face_center = (double **) malloc(sizeof(double *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {
      face_center[i_part] = PDM_array_const_double(n_face[i_part] * 3, 0.);

      if (with_edges) {
        for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
          for (int i_edge = face_vtx_idx[i_part][i_face]; i_edge < face_vtx_idx[i_part][i_face+1]; i_edge++) {
            int edge_id = std::abs(face_edge[i_part][i_edge]) - 1;
            for (int i_vtx = 2*edge_id; i_vtx < 2*(edge_id+1); i_vtx++) {
              int vtx_id = edge_vtx[i_part][i_vtx] - 1;
              for (int j = 0; j < 3; j++) {
                face_center[i_part][3*i_face+j] += vtx_coord[i_part][3*vtx_id+j];
              }
            }
          }

          double normalization = 0.5/(double) (face_vtx_idx[i_part][i_face+1] - face_vtx_idx[i_part][i_face]);
          for (int j = 0; j < 3; j++) {
            face_center[i_part][3*i_face+j] *= normalization;
          }
        }
      }
      else {
        for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
          for (int i_vtx = face_vtx_idx[i_part][i_face]; i_vtx < face_vtx_idx[i_part][i_face+1]; i_vtx++) {
            int vtx_id = face_vtx[i_part][i_vtx] - 1;
            for (int j = 0; j < 3; j++) {
              face_center[i_part][3*i_face+j] += vtx_coord[i_part][3*vtx_id+j];
            }
          }

          double normalization = 1./(double) (face_vtx_idx[i_part][i_face+1] - face_vtx_idx[i_part][i_face]);
          for (int j = 0; j < 3; j++) {
            face_center[i_part][3*i_face+j] *= normalization;
          }
        }
      }
    }

    /* Compute global numbering */
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_gnum_set_from_coords(gen_gnum, i_part, n_face[i_part], face_center[i_part], NULL);
    }

    PDM_gnum_compute(gen_gnum);

    for (int i_part = 0; i_part < n_part; i_part++) {
      face_ln_to_gn[i_part] = const_cast<CWP_g_num_t*>(PDM_gnum_get(gen_gnum, i_part));
    }

    PDM_gnum_free(gen_gnum);

    for (int i_part = 0; i_part < n_part; i_part++) {
      free(face_center[i_part]);
    }
    free(face_center);
  }

  /**
    * \brief Mesh constructor
    *
    * Construct the CWIPI mesh by using paradigm nodal methods.
    *
    * \param [in] npart Number of mesh partitions.
    *
    */

  Mesh::Mesh
  (
    const MPI_Comm &localComm,
    const int npart,
    const CWP_Dynamic_mesh_t displacement,
    Coupling *cpl
  )
  :
    _localComm(localComm),
    _nBlocks(0),
    _blocksType(nullptr),  
                  //_hoOrdering (NULL),
    // _visu(visu),
    _displacement(displacement),
    _cpl(cpl),
    _faceEdgeMethod(0),
    _faceVtxMethod(0),
    _cellFaceMethod(0),
    _pdmNodal_handle_index(),
    _isVtxGnumComputed(false),
    _isEltGnumComputed(false),
    _isFaceGnumComputed(false)
  {

    _pdm_localComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&localComm));
    // pdm_nodal building

    _npart                 = npart;
    _nVertex   .resize(npart,0);
    _nElts     .resize(npart,0);
    _connec_idx.resize(npart,NULL);
    _connec    .resize(npart,NULL);
    _gnum_elt  .resize(npart,NULL);
    _elt_centers  .resize(npart,NULL);
    _elt_id_block  .resize(npart,NULL);
    _elt_in_block  .resize(npart,NULL);

    _coords .resize(npart,NULL);
    _global_num_vtx .resize(npart,NULL);
    _global_num_elt .resize(npart,NULL);
    _global_num_face.resize(npart,NULL);

    _nCells      .resize(npart, 0)  ;
    _cellFaceIdx .resize(npart, NULL);
    _cellFace    .resize(npart, NULL);
    _nFace       .resize(npart,0)   ;
    _faceEdgeIdx .resize(npart,NULL);
    _faceEdge    .resize(npart,NULL);
    _faceVtxIdx  .resize(npart,NULL);
    _faceVtx     .resize(npart,NULL);
    _nEdge       .resize(npart,0)   ;
    _edgeVtx     .resize(npart,NULL);

    _faceLNToGN  .resize(npart,NULL);
    _cellLNToGN  .resize(npart,NULL);

    switch (_cpl->entitiesDimGet()) {
      case CWP_INTERFACE_POINT:
        _geom_kind = PDM_GEOMETRY_KIND_CORNER;
        break;
      case CWP_INTERFACE_LINEAR:
        _geom_kind = PDM_GEOMETRY_KIND_RIDGE;
        break;
      case CWP_INTERFACE_SURFACE:
        _geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
        break;
      case CWP_INTERFACE_VOLUME:
        _geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
        break;
      default:
        PDM_error(__FILE__, __LINE__, 0, "Invalid entity dimension %d\n", _cpl->entitiesDimGet());
    }

  }

  Mesh::~Mesh()
  {
    for (int i = 0; i < _npart; i++) {
      free (_connec_idx[i]);
      free (_connec[i]);
      free (_gnum_elt[i]);
      if (_elt_centers[i] != NULL) {
        free (_elt_centers[i]);
      }
      if (_elt_id_block[i] != NULL) {
        free (_elt_id_block[i]);
      }
      if (_elt_in_block[i] != NULL) {
        free (_elt_in_block[i]);
      }

      if (_isFaceGnumComputed) {
        free (_faceLNToGN[i]);
        // free (_global_num_face[i]);
      }

      if (_isEltGnumComputed) {
        free (_global_num_elt[i]); 
      }

      if (_isVtxGnumComputed) {
        free (_global_num_vtx[i]); 
      }
    }

    for (int i = 0; i < (int) _blockDB.size(); i++) {
      delete _blockDB[i];
    }

    if (_blocksType != NULL) {
      free ( _blocksType);
    }

  }



  int
  Mesh::_Mesh_nodal_block_std_type_size_get
  (
    CWP_Block_t type
  )
  {

    switch (type) {

      case CWP_BLOCK_NODE: return 1;
      break;

      case CWP_BLOCK_EDGE2: return 2;
      break;

      case CWP_BLOCK_FACE_TRIA3: return 3;
      break;

      case CWP_BLOCK_FACE_QUAD4: return 4;
      break;

      case CWP_BLOCK_CELL_TETRA4: return 4;
      break;

      case CWP_BLOCK_CELL_HEXA8: return 8;
      break;

      case CWP_BLOCK_CELL_PYRAM5: return 5;
      break;

      case CWP_BLOCK_CELL_PRISM6: return 6;
      break;

      default: return -1;
      PDM_error(__FILE__, __LINE__, 0, "This CWP_Block_t is not available as function argument.\n");

    }
  }

  void
  Mesh::eltCentersCompute
  (
    int i_part
  )
  {

    int n_elt_part = getPartNElts(i_part);

    if (_elt_centers[i_part] == NULL) {
      _elt_centers[i_part] = (double*)malloc(sizeof(double)*3* n_elt_part);
    }


    int ind_idx=0;
    for(int i=0;i<_nBlocks;i++){
      int n_elt = _blockDB[i]->NEltsGet()[i_part];

      PDM_part_mesh_nodal_section_elt_center_reset(_pdmNodal_handle_index,
                                                   _blockDB[i]->blockIDPDMGet(),
                                                   i_part);
      const double* elt_centers_block = _blockDB[i]->eltCentersGet(i_part);

      for(int j=0;j<n_elt;j++){
        for(int k=0;k<3;k++){
          _elt_centers[i_part][3*ind_idx+k]=elt_centers_block[3*j+k];
        }
        ind_idx++;
      }//end loop j
    }//end loop on block
  }


  CWP_g_num_t*
  Mesh::GNumEltsGet
  (
    int i_part
  )
  {
    // if (_connec_idx[i_part]==NULL) {
    //   connecCompute(i_part);
    // }//end if NULL

    return PDM_part_mesh_nodal_g_num_get_from_part (_pdmNodal_handle_index, _geom_kind, i_part, PDM_OWNERSHIP_KEEP);
  }


  double*
  Mesh::eltCentersGet
  (
    int i_part
  )
  {
    if(_elt_centers[i_part]==NULL || _displacement != CWP_DYNAMIC_MESH_STATIC) {
      eltCentersCompute(i_part);
    }//end if NULL

    return _elt_centers[i_part];
  }


  void
  Mesh::coordSet
  (
    const int   i_part,
    const int   n_vtx,
    double      coords[],
    CWP_g_num_t global_num[]
  )
  {
    _coords[i_part] = coords;
    _nVertex[i_part]  = n_vtx;
    _global_num_vtx[i_part] = global_num;
  }



  /**************************************************************/


  void
  Mesh::updateBlockDB()
  {
    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(_pdmNodal_handle_index,
                                                                                      _geom_kind);


    int  n_block     = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
    int *section_ids = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

    for (int i_block = 0; i_block < n_block; i_block++) {

      CWP_Block_t t_block = (CWP_Block_t) PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                          section_ids[i_block]);

      int block_id = blockAdd((CWP_Block_t) t_block);

      for (int i_part = 0; i_part < _npart; i_part++) {
        int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                                section_ids[i_block],
                                                                i_part);

        if (t_block == CWP_BLOCK_FACE_POLY) {
          int         *connec     = NULL;
          int         *connec_idx = NULL;
          CWP_g_num_t *g_num      = NULL;
          PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                       section_ids[i_block],
                                                       i_part,
                                                       &connec_idx,
                                                       &connec,
                                                       PDM_OWNERSHIP_KEEP);

          g_num = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                      section_ids[i_block],
                                                      i_part,
                                                      PDM_OWNERSHIP_KEEP);

          poly2DBlockSet(i_part,
                         block_id,
                         n_elt,
                         connec_idx,
                         connec,
                         g_num);
        }
        else if (t_block == CWP_BLOCK_CELL_POLY) {
          int          n_face              = 0;
          CWP_g_num_t *face_g_num          = NULL;
          int         *face_vtx_idx        = NULL;
          int         *face_vtx            = NULL;
          CWP_g_num_t *cell_g_num          = NULL;
          int         *cell_face_idx       = NULL;
          int         *cell_face           = NULL;
          int         *_parent_num         = NULL;
          CWP_g_num_t *parent_entity_g_num = NULL;
          PDM_part_mesh_nodal_elmts_section_poly3d_get(pmne,
                                                       section_ids[i_block],
                                                       i_part,
                                                       &n_face,
                                                       &face_g_num,
                                                       &face_vtx_idx,
                                                       &face_vtx,
                                                       &cell_g_num,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       &_parent_num,
                                                       &parent_entity_g_num,
                                                       PDM_OWNERSHIP_KEEP);

          poly3DBlockSet(i_part,
                         block_id,
                         n_elt,
                         n_face,
                         face_vtx_idx,
                         face_vtx,
                         cell_face_idx,
                         cell_face,
                         cell_g_num);
        }
        else {
          int         *connec              = NULL;
          CWP_g_num_t *g_num               = NULL;
          int         *_parent_num         = NULL;
          CWP_g_num_t *parent_entity_g_num = NULL;
          PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                    section_ids[i_block],
                                                    i_part,
                                                    &connec,
                                                    &g_num,
                                                    &_parent_num,
                                                    &parent_entity_g_num,
                                                    PDM_OWNERSHIP_KEEP);

          stdBlockSet(i_part,
                      block_id,
                      n_elt,
                      connec,
                      g_num);
        }



        /* Parent numbering */
        int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                   section_ids[i_block],
                                                                   i_part,
                                                                   PDM_OWNERSHIP_KEEP);

        _blockDB[block_id]->blockSetParentNum(i_part,
                                              parent_num,
                                              PDM_OWNERSHIP_USER); // hack : owner is pmne
      }

      int i_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(_pdmNodal_handle_index,
                                                                        _geom_kind,
                                                                        section_ids[i_block]);
      _blockDB[block_id]->blockIDPDMSet(i_section);
    }
  }



  /**********************************************************************/

  void
  Mesh::geomFinalize
  (
  )
  {
    int unionRank;
    MPI_Comm_rank(_cpl->communicationGet()->unionCommGet(),&unionRank);

    int globalRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);

    int mesh_dimension;
    switch (_cpl->entitiesDimGet()) {
      case CWP_INTERFACE_POINT:
        mesh_dimension = 0;
        break;
      case CWP_INTERFACE_LINEAR:
        mesh_dimension = 1;
        break;
      case CWP_INTERFACE_SURFACE:
        mesh_dimension = 2;
        break;
      case CWP_INTERFACE_VOLUME:
        mesh_dimension = 3;
        break;
      default:
        PDM_error(__FILE__, __LINE__, 0, "Invalid entity dimension %d\n", _cpl->entitiesDimGet());
    }

    _pdmNodal_handle_index = PDM_part_mesh_nodal_create(mesh_dimension, _npart,_pdm_localComm);

    if (gnumVtxRequired()) {
      _isVtxGnumComputed = true;

      PDM_gen_gnum_t *pdmGNum_handle_index = PDM_gnum_create(3, _npart, PDM_FALSE, 1e-3, _pdm_localComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

      for (int i_part = 0; i_part < _npart ; i_part++) {
        PDM_gnum_set_from_coords(pdmGNum_handle_index, i_part, _nVertex[i_part], _coords[i_part], NULL);
      }

      PDM_gnum_compute(pdmGNum_handle_index);

      for (int i_part=0; i_part < _npart; i_part++) {
        _global_num_vtx[i_part] = const_cast<CWP_g_num_t*>(PDM_gnum_get(pdmGNum_handle_index, i_part));
      }

      PDM_gnum_free(pdmGNum_handle_index);
    }

    for (int i_part = 0; i_part < _npart; i_part++) {

      PDM_part_mesh_nodal_coord_set(_pdmNodal_handle_index ,
                                    i_part                 ,
                                    _nVertex       [i_part],
                                    _coords        [i_part],
                                    _global_num_vtx[i_part],
                                    PDM_OWNERSHIP_USER);

      // if(_visu->isCreated() && _displacement == CWP_DYNAMIC_MESH_STATIC) {
      //   _visu->GeomCoordSet(i_part,
      //   _nVertex       [i_part],
      //   _coords        [i_part],
      //   _global_num_vtx[i_part]);
      // }

    }//loop i_part


    if (_faceEdgeMethod == 1) {

      int compute_gnum = 0;
      for (int i_part = 0; i_part < _npart; i_part++) {
        if (_faceLNToGN[i_part] == NULL && _nFace[i_part] > 0) {
          compute_gnum = 1;
          break;
        }
      }

      if(compute_gnum) {
        _isEltGnumComputed = true;

        _gen_face_gnum(true,
                       _pdm_localComm,
                       _npart,
                       _nFace,
                       _faceVtxIdx,
                       _faceVtx,
                       _faceEdge,
                       _edgeVtx,
                       _coords,
                       _faceLNToGN);
      }

      for (int i_part = 0; i_part < _npart; i_part++) {

        PDM_part_mesh_nodal_face2d_faceedge_add(_pdmNodal_handle_index,
                                                i_part,
                                                _nFace      [i_part],
                                                _nEdge      [i_part],
                                                _edgeVtx    [i_part],
                                                _faceEdgeIdx[i_part],
                                                _faceEdge   [i_part],
                                                _faceLNToGN [i_part],
                                                PDM_OWNERSHIP_KEEP);

      }//end i_part loop

      updateBlockDB();

    }

    else if (_faceVtxMethod == 1) {

      int compute_gnum = 0;
      for (int i_part = 0; i_part < _npart; i_part++) {
        if (_faceLNToGN[i_part] == NULL && _nFace[i_part] > 0) {
          compute_gnum = 1;
          break;
        }
      }

      if (compute_gnum) {
        _isEltGnumComputed = true;

        _gen_face_gnum(false,
                       _pdm_localComm,
                       _npart,
                       _nFace,
                       _faceVtxIdx,
                       _faceVtx,
                       _faceEdge,
                       _edgeVtx,
                       _coords,
                       _faceLNToGN);
      }

      for (int i_part = 0; i_part < _npart; i_part++) {
        PDM_part_mesh_nodal_faces_facevtx_add(_pdmNodal_handle_index,
                                              i_part,
                                              _nFace     [i_part],
                                              _faceVtxIdx[i_part],
                                              _faceVtx   [i_part],
                                              _faceLNToGN[i_part],
                                              PDM_OWNERSHIP_KEEP);
      }

      updateBlockDB();
    }

    else if (_cellFaceMethod == 1) {

      int compute_face_gnum = 0;
      for (int i_part = 0; i_part < _npart; i_part++) {
        if (_faceLNToGN[i_part] == NULL && _nCells[i_part] > 0) {
          compute_face_gnum = 1;
          break;
        }
      }

      if (compute_face_gnum) {

        _isFaceGnumComputed = true;

        _gen_face_gnum(false,
                       _pdm_localComm,
                       _npart,
                       _nFace,
                       _faceVtxIdx,
                       _faceVtx,
                       _faceEdge,
                       _edgeVtx,
                       _coords,
                       _faceLNToGN);
      }


      int compute_gnum = 0;
      for (int i_part = 0; i_part < _npart; i_part++) {
        if (_cellLNToGN[i_part] == NULL) {
          compute_gnum = 1;
          break;
        }
      }

      if (compute_gnum) {

        _isEltGnumComputed = true;

        PDM_gen_gnum_t *pdmGNum_handle_index = PDM_gnum_create(3, _npart, PDM_FALSE, 1e-3, _pdm_localComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        double ** cell_center =  (double **) malloc (sizeof(double *) * (_npart));

        for (int i_part = 0; i_part < _npart; i_part++) {
          cell_center[i_part] =  (double *) malloc (sizeof(double) * (3*_nCells[i_part]));
          for (int j = 0; j < 3*_nCells[i_part]; j++) {
            cell_center[i_part][j] = 0.;
          }

          for (int j = 0; j < _nCells[i_part]; j++) {

            int weights = 0;
            int idx_face = _cellFaceIdx[i_part][j];
            int nb_face = _cellFaceIdx[i_part][j] - idx_face;

            for (int j1 = idx_face; j1 < idx_face + nb_face; j1++) {

              int i_face = std::abs(_cellFace[i_part][j1]) - 1;

              int idx = _faceVtxIdx[i_part][i_face];
              int nb = _faceVtxIdx[i_part][i_face] - idx;

              for (int k = idx; k < idx + nb; k++) {

                int i_vtx = _faceVtx[i_part][k] - 1;

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*j+k2] += _coords[i_part][3*i_vtx+k2];
                }
                weights++;
              }
            }

            for (int k2 = 0; k2 < 3; k2++) {
              cell_center[i_part][3*j+k2] /= weights;
            }
          }
        }

        for (int i_part = 0; i_part < _npart; i_part++) {
          PDM_gnum_set_from_coords (pdmGNum_handle_index, i_part, _nCells[i_part], cell_center[i_part], NULL);
        }

        PDM_gnum_compute (pdmGNum_handle_index);

        for(int i_part=0;i_part<_npart;i_part++) {
          _cellLNToGN[i_part] = const_cast<CWP_g_num_t*>(PDM_gnum_get (pdmGNum_handle_index, i_part));
        }

        PDM_gnum_free (pdmGNum_handle_index);

        for (int i_part = 0; i_part < _npart; i_part++) {
          free ( cell_center[i_part]);
        }

        free ( cell_center);

      }

      for(int i_part=0;i_part<_npart;i_part++){
        PDM_part_mesh_nodal_cell3d_cellface_add(_pdmNodal_handle_index,
                                                i_part,
                                                _nCells[i_part],
                                                _nFace[i_part],
                                                _faceVtxIdx[i_part],
                                                _faceVtx[i_part],
                                                _faceLNToGN[i_part],
                                                _cellFaceIdx[i_part],
                                                _cellFace[i_part],
                                                _cellLNToGN[i_part],
                                                PDM_OWNERSHIP_KEEP);


      }//end i_part loop

      updateBlockDB();

    }

    else {

      for(int i_block=0;i_block<_nBlocks;i_block++){
        CWP_Block_t block_type = _blockDB[i_block]->blockTypeGet();
        PDM_Mesh_nodal_elt_t pdm_block_type = (PDM_Mesh_nodal_elt_t) block_type;

        _blockDB[i_block]->blockIDPDMSet(PDM_part_mesh_nodal_section_add(_pdmNodal_handle_index, pdm_block_type));

      } //end loop on block

      // Create face_g_num for CELL_POLY
      for (int i_block = 0; i_block < _nBlocks; i_block++) {
        CWP_Block_t block_type = _blockDB[i_block]->blockTypeGet();

        if (block_type == CWP_BLOCK_CELL_POLY) {
          BlockCP *block = dynamic_cast<BlockCP *>(_blockDB[i_block]);

          std::vector<int>           n_face       = block->NFacesGet();
          std::vector<int *>         face_vtx_idx = block->ConnecFacesIDXGet();
          std::vector<int *>         face_vtx     = block->ConnecFacesGet();
          std::vector<CWP_g_num_t *> face_ln_to_gn(_npart, NULL);

          _gen_face_gnum(false,
                         _pdm_localComm,
                         _npart,
                         n_face,
                         face_vtx_idx,
                         face_vtx,
                         _faceEdge, // unused
                         _edgeVtx,  // unused
                         _coords,
                         face_ln_to_gn);

          block->FacesGNumSet(face_ln_to_gn);
        }
      }

      int compute_gnum = 0;

      for(int i_part = 0; i_part < _npart; i_part++){

        for(int i_block = 0; i_block < _nBlocks; i_block++){

          if (_blockDB[i_block]->GNumMeshGet(i_part) == NULL && _blockDB[i_block]->NEltsGet(i_part) > 0) {
            compute_gnum = 1;
            break;
          }
        }
      }

      if (compute_gnum) {

        _isEltGnumComputed = true;

        PDM_gen_gnum_t *pdmGNum_handle_index = PDM_gnum_create(3, _npart, PDM_FALSE, 1e-3, _pdm_localComm, PDM_OWNERSHIP_KEEP);

        double **cell_center =  (double **) malloc (sizeof(double *) * (_npart));

        for(int i_part = 0; i_part < _npart; i_part++){
          cell_center[i_part] =  (double *) malloc (sizeof(double) * (3*_nElts[i_part]));

          for(int i = 0; i < 3*_nElts[i_part]; i++){
            cell_center[i_part][i] = 0.;
          }

          int ielt = 0;
          for(int i_block = 0; i_block < _nBlocks; i_block++){

            int n_elt = _blockDB[i_block]->NEltsGet(i_part);

            CWP_Block_t block_type = _blockDB[i_block]->blockTypeGet();

            if (block_type == CWP_BLOCK_FACE_POLY) {
              BlockFP *block = dynamic_cast<BlockFP *>(_blockDB[i_block]);

              for (int j = 0; j < n_elt; j++) {
                int idx = block->ConnecIDXGet()[i_part][j];
                int nb  = block->ConnecIDXGet()[i_part][j+1] - idx;

                for (int k = idx; k < idx + nb; k++) {
                  int ivtx = block->ConnecGet()[i_part][k] - 1;
                  for (int k2 = 0; k2 < 3; k2++) {
                    cell_center[i_part][3*ielt+k2] += _coords[i_part][3*ivtx+k2];
                  }
                }

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*ielt+k2] /= nb;
                }

                ielt++;
              }

            }

            else if (block_type == CWP_BLOCK_CELL_POLY) {
              BlockCP *block = dynamic_cast<BlockCP *>(_blockDB[i_block]);

              for (int j = 0; j < n_elt; j++) {

                int weights = 0;
                int idx_face = block->ConnecIDXGet()[i_part][j];
                int nb_face = block->ConnecIDXGet()[i_part][j+1] - idx_face;

                for (int j1 = idx_face; j1 < idx_face + nb_face; j1++) {

                  int i_face = std::abs(block->ConnecGet()[i_part][j1]) - 1;

                  int idx = block->ConnecFacesIDXGet()[i_part][i_face];
                  int nb  = block->ConnecFacesIDXGet()[i_part][i_face+1] - idx;

                  for (int k = idx; k < idx + nb; k++) {

                    int i_vtx = block->ConnecFacesGet()[i_part][k] - 1;

                    for (int k2 = 0; k2 < 3; k2++) {
                      cell_center[i_part][3*ielt+k2] += _coords[i_part][3*i_vtx+k2];
                    }
                    weights++;
                  }
                }

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*ielt+k2] /= weights;
                }

                ielt++;
              }
            }

            else if (block_type == CWP_BLOCK_NODE        ||
                     block_type == CWP_BLOCK_EDGE2       ||
                     block_type == CWP_BLOCK_FACE_TRIA3  ||
                     block_type == CWP_BLOCK_FACE_QUAD4  ||
                     block_type == CWP_BLOCK_CELL_TETRA4 ||
                     block_type == CWP_BLOCK_CELL_HEXA8  ||
                     block_type == CWP_BLOCK_CELL_PRISM6 ||
                     block_type == CWP_BLOCK_CELL_PYRAM5) {
              BlockStd *block = dynamic_cast<BlockStd *>(_blockDB[i_block]);

              PDM_Mesh_nodal_elt_t elt_type = (PDM_Mesh_nodal_elt_t) block_type;

              int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, 1);

              for (int j = 0; j < n_elt; j++) {

                for (int k = 0; k < n_vtx_elt; k++) {

                  int i_vtx = block->ConnecGet()[i_part][j*n_vtx_elt + k] - 1;

                  for (int k2 = 0; k2 < 3; k2++) {
                    cell_center[i_part][3*ielt+k2] += _coords[i_part][3*i_vtx+k2];
                  }  
                }

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*ielt+k2] /= n_vtx_elt;
                }
                ielt++;
              }

              // block->ConnecGet()[i_part];

            }

            else {
              BlockHO *block = dynamic_cast<BlockHO *>(_blockDB[i_block]);

              PDM_Mesh_nodal_elt_t elt_type = (PDM_Mesh_nodal_elt_t) block_type;

              int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, block->OrderGet());

              for (int j = 0; j < n_elt; j++) {

                for (int k = 0; k < n_vtx_elt; k++) {

                  int i_vtx = block->ConnecGet()[i_part][j*n_vtx_elt + k] - 1;

                  for (int k2 = 0; k2 < 3; k2++) {
                    cell_center[i_part][3*ielt+k2] += _coords[i_part][3*i_vtx+k2];
                  }
                }

                for (int k2 = 0; k2 < 3; k2++) {
                  cell_center[i_part][3*ielt+k2] /= n_vtx_elt;
                }
                ielt++;
              }

            }
          }

          PDM_gnum_set_from_coords (pdmGNum_handle_index, i_part, _nElts[i_part], cell_center[i_part], NULL);

        }

        PDM_gnum_compute (pdmGNum_handle_index);

        for (int i_part = 0; i_part < _npart; i_part++) {

          PDM_g_num_t *gnum = const_cast<PDM_g_num_t*>(PDM_gnum_get (pdmGNum_handle_index, i_part));

          for(int i_block = 0; i_block < _nBlocks; i_block++){

            int n_elt = _blockDB[i_block]->NEltsGet(i_part);
            CWP_g_num_t *block_gnum = (CWP_g_num_t *) malloc (sizeof(CWP_g_num_t) * n_elt);

            memcpy(block_gnum, gnum, n_elt * sizeof(CWP_g_num_t));

            _blockDB[i_block]->GNumMeshSet(i_part, block_gnum, PDM_OWNERSHIP_KEEP);

          }


        }

        PDM_gnum_free (pdmGNum_handle_index);

        for (int i_part = 0; i_part < _npart; i_part++) {
          free ( cell_center[i_part]);
        }

        free ( cell_center);

      }

      for(int i_part = 0; i_part < _npart; i_part++){

        int parent_id = 0;

        for(int i_block = 0; i_block < _nBlocks; i_block++){
          int n_elt = _blockDB[i_block]->NEltsGet(i_part);

          int *parent_num = (int *) malloc(sizeof(int) * n_elt);
          for (int i = 0; i < n_elt; i++) {
            parent_num[i] = parent_id++;
          }

          _blockDB[i_block]->blockSetParentNum(i_part,
                                               parent_num,
                                               PDM_OWNERSHIP_KEEP);

          CWP_Block_t block_type = _blockDB[i_block]->blockTypeGet();

          int pdm_id_block = _blockDB[i_block]->blockIDPDMGet();

          if (block_type == CWP_BLOCK_FACE_POLY) {
            BlockFP *block = dynamic_cast<BlockFP *>(_blockDB[i_block]);

            PDM_part_mesh_nodal_section_poly2d_set(_pdmNodal_handle_index,
                                                   pdm_id_block,
                                                   i_part,
                                                   n_elt,
                                                   block->ConnecIDXGet()[i_part],
                                                   block->ConnecGet()[i_part],
                                                   _blockDB[i_block]->GNumMeshGet(i_part),
                                                   _blockDB[i_block]->ParentNumGet(i_part),//NULL,
                                                   PDM_OWNERSHIP_USER);

          }

          else if (block_type == CWP_BLOCK_CELL_POLY) {
            BlockCP *block = dynamic_cast<BlockCP *>(_blockDB[i_block]);

            PDM_part_mesh_nodal_section_poly3d_set(_pdmNodal_handle_index,
                                                   pdm_id_block,
                                                   i_part,
                                                   n_elt,
                                                   block->NFacesGet()[i_part],
                                                   block->ConnecFacesIDXGet()[i_part],
                                                   block->ConnecFacesGet()[i_part],
                                                   block->FacesGNumGet()[i_part],
                                                   block->ConnecIDXGet()[i_part],
                                                   block->ConnecGet()[i_part],
                                                   block->GNumMeshGet(i_part),
                                                   _blockDB[i_block]->ParentNumGet(i_part),//NULL,
                                                   NULL,
                                                   PDM_OWNERSHIP_USER);

          }

          else if (block_type == CWP_BLOCK_NODE        ||
                   block_type == CWP_BLOCK_EDGE2       ||
                   block_type == CWP_BLOCK_FACE_TRIA3  ||
                   block_type == CWP_BLOCK_FACE_QUAD4  ||
                   block_type == CWP_BLOCK_CELL_TETRA4 ||
                   block_type == CWP_BLOCK_CELL_HEXA8  ||
                   block_type == CWP_BLOCK_CELL_PRISM6 ||
                   block_type == CWP_BLOCK_CELL_PYRAM5) {
            BlockStd *block = dynamic_cast<BlockStd *>(_blockDB[i_block]);

            PDM_part_mesh_nodal_section_std_set(_pdmNodal_handle_index,
                                                pdm_id_block,
                                                i_part,
                                                n_elt,
                                                block->ConnecGet()[i_part],
                                                block->GNumMeshGet(i_part),
                                                _blockDB[i_block]->ParentNumGet(i_part),//NULL,
                                                NULL,
                                                PDM_OWNERSHIP_USER);
          }

          else {
            BlockHO *block = dynamic_cast<BlockHO *>(_blockDB[i_block]);

            PDM_part_mesh_nodal_section_std_ho_set(_pdmNodal_handle_index,
                                                   pdm_id_block,
                                                   i_part,
                                                   n_elt,
                                                   block->ConnecIJKGet()[i_part],
                                                   block->GNumMeshGet(i_part),
                                                   _blockDB[i_block]->ParentNumGet(i_part),//NULL,
                                                   NULL,
                                                   block->OrderGet(),
                                                   NULL,
                                                   PDM_OWNERSHIP_USER);
          }
        }
      }   //end loop on block
    }

 
    _blocksType = (CWP_Block_t *) malloc (sizeof(CWP_Block_t) * _nBlocks);

    for (int i_block = 0 ; i_block < _nBlocks ; i_block++) {
      _blockDB[i_block]->geomFinalize();
      _blocksType[i_block] = _blockDB[i_block]->blockTypeGet(); 
    } //Loop on blockDB

    for (int i_part = 0; i_part < _npart; i_part++) {

      _elt_id_block[i_part] = (int *) malloc (sizeof(int) * _nElts[i_part]);
      _elt_in_block[i_part] = (int *) malloc (sizeof(int) * _nElts[i_part]);

      for(int i_block = 0; i_block < _nBlocks; i_block++){
        int n_elt = _blockDB[i_block]->NEltsGet(i_part);

        int pdm_id_block = _blockDB[i_block]->blockIDPDMGet();

        int *parent_num = PDM_part_mesh_nodal_section_parent_num_get(_pdmNodal_handle_index, pdm_id_block, i_part, PDM_OWNERSHIP_KEEP);

        if (parent_num != NULL) {
          for(int i_elt = 0; i_elt < n_elt; i_elt++){
            (_elt_id_block[i_part])[parent_num[i_elt]] = i_block;
            (_elt_in_block[i_part])[parent_num[i_elt]] = i_elt + 1;
          }
        }
        else {
          for(int i_elt = 0; i_elt < n_elt; i_elt++){
            _elt_id_block[i_part][i_elt] = i_block;
            _elt_in_block[i_part][i_elt] = i_elt + 1;
          }        
        }

      }

    }

    _nBlocks     = PDM_part_mesh_nodal_n_section_in_geom_kind_get (_pdmNodal_handle_index, _geom_kind);
    _blocks_id   = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(_pdmNodal_handle_index, _geom_kind);

    // if(_visu->isCreated() && _displacement == CWP_DYNAMIC_MESH_STATIC ) {
    //   _visu->GeomWrite(this);
    // }

  }


  void
  Mesh::stdBlockSet
  (
    const int              i_part,
    const int              block_id,
    const int              n_elts,
    int                    connec[],
    CWP_g_num_t            global_num[]
  )
  {

    BlockStd *block = dynamic_cast <BlockStd *> (_blockDB [block_id]);
    block->blockSet(i_part,n_elts,connec,global_num);

    _nElts[i_part]  += n_elts;

  }

  void
  Mesh::HOBlockSet
  (
   const int    i_part,
   const int    block_id,
   const int    n_elts,
   int          connec[],
   CWP_g_num_t  global_num[],
   const int    order,
   const char  *ho_ordering
   )
  {

    BlockHO *block = dynamic_cast <BlockHO *> (_blockDB [block_id]);
    block->blockSet(i_part,n_elts,connec,global_num, order, ho_ordering);

    _nElts[i_part]  += n_elts;

  }

  CWP_Block_t
  Mesh::stdBlockTypeGet
  (
    const int              block_id
  )
  {
    BlockStd *block = dynamic_cast <BlockStd *> (_blockDB [block_id]);
    CWP_Block_t block_type = block->blockTypeGet();

    return block_type;
  }

  /*************************************************/

  void
  Mesh::stdBlockGet
  (
    const int              i_part,
    const int              block_id,
    int                   *n_elts,
    int                  **connec,
    CWP_g_num_t          **global_num
  )
  {

    BlockStd *block = dynamic_cast <BlockStd *> (_blockDB [block_id]);

    block->blockGet(i_part,
                    n_elts,
                    connec,
                    global_num);

  }

  /*************************************************/

  void
  Mesh::HOBlockGet
  (
    const int              i_part,
    const int              block_id,
    int                   *n_elts,
    int                   *order,
    int                  **connec,
    CWP_g_num_t          **global_num
  )
  {

    BlockHO *block = dynamic_cast <BlockHO *> (_blockDB [block_id]);

    char *ho_ordering;

    block->blockGet(i_part,
                    n_elts,
                    connec,
                    global_num,
                    order,
                    &ho_ordering);

  }

  /*************************************************/

  void
  Mesh::poly2DBlockSet
  (
    const int              i_part,
    const int              block_id,
    const int              n_elts,
    int                    connec_idx[],
    int                    connec[],
    CWP_g_num_t            global_num[]
  )
  {

    BlockFP *block = dynamic_cast <BlockFP *> (_blockDB [block_id]);

    block->blockSet(i_part,
                    n_elts,
                    connec_idx,
                    connec,
                    global_num);
    _nElts[i_part]  += n_elts;

  }

  /*************************************************/

  void
  Mesh::poly2DBlockGet
  (
    const int              i_part,
    const int              block_id,
    int                   *n_elts,
    int                  **connec_idx,
    int                  **connec,
    CWP_g_num_t          **global_num
  )
  {

    BlockFP *block = dynamic_cast <BlockFP *> (_blockDB [block_id]);

    block->blockGet(i_part,
                    n_elts,
                    connec_idx,
                    connec,
                    global_num);
  
  }

  /**********************************************************************/


  void
  Mesh::poly3DBlockSet
  (
    const int              i_part,
    const int              block_id,
    const int              n_elts,
    const int              n_faces,
    int                    connec_faces_idx[],
    int                    connec_faces[],
    int                    connec_cells_idx[],
    int                    connec_cells[],
    CWP_g_num_t            global_num[]
  )
  {

    BlockCP *block = dynamic_cast <BlockCP *> (_blockDB [block_id]);
    block->blockSet(i_part,
                    n_elts,
                    n_faces,
                    connec_faces_idx,
                    connec_faces,
                    connec_cells_idx,
                    connec_cells,
                    global_num);

    _nElts[i_part]  +=  n_elts;

  }


  void 
  Mesh::poly3DBlockGet
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
  )
  {
    BlockCP *block = dynamic_cast <BlockCP *> (_blockDB [block_id]);
    block->blockGet(i_part,
                    n_elts,
                    n_faces,
                    connec_faces_idx,
                    connec_faces,
                    connec_cells_idx,
                    connec_cells,
                    global_num);
  }

  void
  Mesh::meshDel()
  {
    if (_pdmNodal_handle_index != NULL) {
      PDM_part_mesh_nodal_free(_pdmNodal_handle_index);
    }
    _pdmNodal_handle_index = NULL;


    for (int i = 0; i < _npart; i++) {
      if (_elt_id_block[i] != NULL) {
        free(_elt_id_block[i]);
        _elt_id_block[i] = NULL;
      }

      if (_elt_in_block[i] != NULL) {
        free(_elt_in_block[i]);
        _elt_in_block[i] = NULL;
      }

      for (int i_block = 0; i_block < _nBlocks; i_block++) {
        _blockDB[i_block]->ParentNumFree(i);
        _blockDB[i_block]->GNumMeshFree(i);

        CWP_Block_t block_type = _blockDB[i_block]->blockTypeGet();
        if (block_type == CWP_BLOCK_CELL_POLY) {
          BlockCP *block = dynamic_cast<BlockCP *>(_blockDB[i_block]);
          block->FacesGNumFree(i);
        }
      }

      if (_isVtxGnumComputed) {
        if (_global_num_vtx[i] != NULL) {
          free(_global_num_vtx[i]);
          _global_num_vtx[i] = NULL;
        }
      }

      if (_isEltGnumComputed) {
        if (_faceEdgeMethod == 1 || _faceVtxMethod == 1) {
          if (_faceLNToGN[i] != NULL) {
            free(_faceLNToGN[i]);
            _faceLNToGN[i] = NULL;
          }
        }
        else if (_cellFaceMethod == 1) {
          if (_cellLNToGN[i] != NULL) {
            free(_cellLNToGN[i]);
            _cellLNToGN[i] = NULL;
          }
        }
      }
    }

    free(_blocksType);
    _blocksType = NULL;
    _nBlocks = 0;
    for (int i = 0; i < (int) _blockDB.size(); i++) {
      delete _blockDB[i];
    }
    _blockDB.clear();
  }


  int
  Mesh::blockAdd
  (
    const CWP_Block_t      block_type
  )
  {
    Block *myBlock = FB::getInstance().CreateObject(block_type);
    myBlock->BlockAdd(block_type,this);
    int block_id_cwipi = _nBlocks;
    _blockDB.push_back (myBlock);
    // _blockDB.resize(_nBlocks+1);
    // _blockDB[block_id_cwipi] = myBlock;
    myBlock->blockIDCWIPISet(block_id_cwipi);

    _nBlocks = _blockDB.size();
    // _nBlocks++;

    return myBlock->blockIDCWIPIGet();

  }

  void
  Mesh::fromCellFaceSet
  (
    const int   i_part,
    const int   n_cells,
    int         cell_face_idx[],
    int         cell_face[],
    int         n_faces,
    int         face_vtx_idx[],
    int         face_vtx[],
    CWP_g_num_t global_num[]
  )
  {
    _cellFaceMethod = 1;
    _cellLNToGN[i_part] = global_num;

    _faceVtxIdx[i_part] = face_vtx_idx;
    _cellFaceIdx[i_part] = cell_face_idx;
    _faceVtx[i_part] = face_vtx;
    _cellFace[i_part] = cell_face;

    _nFace[i_part] = n_faces;
    _nCells[i_part] = n_cells;

  }


  void
  Mesh::fromFacesEdgeSet
  (
    const int   i_part,
    const int   n_faces,
    int         face_edge_idx[],
    int         face_edge[],
    const int   n_edges,
    int         edge_vtx[],
    CWP_g_num_t global_num[]
  )
  {
    _faceEdgeMethod = 1;

    _faceLNToGN[i_part]  = global_num;
    _faceEdgeIdx[i_part] = face_edge_idx;
    _edgeVtx[i_part]     = edge_vtx;
    _faceEdge[i_part]    = face_edge;
    _nEdge[i_part]       = n_edges;
    _nFace[i_part]       = n_faces;
  }


  void
  Mesh::fromFacesVtxSet
  (
    const int   i_part,
    const int   n_faces,
    int         face_vtx_idx[],
    int         face_vtx[],
    CWP_g_num_t global_num[]
  )
  {
    _faceVtxMethod = 1;

    _faceLNToGN[i_part] = global_num;
    _faceVtxIdx[i_part] = face_vtx_idx;
    _faceVtx[i_part]    = face_vtx;
    _nFace[i_part]      = n_faces;
  }


  int Mesh::getPartNElts(int id_part) const
  {
    CWP_Interface_t ed = _cpl->entitiesDimGet();

    if (ed == CWP_INTERFACE_POINT) {
      return PDM_part_mesh_nodal_n_elmts_get(_pdmNodal_handle_index, PDM_GEOMETRY_KIND_CORNER, id_part);
    }
    else if (ed == CWP_INTERFACE_LINEAR) {
      return PDM_part_mesh_nodal_n_elmts_get(_pdmNodal_handle_index, PDM_GEOMETRY_KIND_RIDGE, id_part);
    }
    else if (ed == CWP_INTERFACE_SURFACE) {
      return PDM_part_mesh_nodal_n_elmts_get(_pdmNodal_handle_index, PDM_GEOMETRY_KIND_SURFACIC, id_part);
    }
    else if (ed == CWP_INTERFACE_VOLUME) {
      return PDM_part_mesh_nodal_n_elmts_get(_pdmNodal_handle_index, PDM_GEOMETRY_KIND_VOLUMIC, id_part);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Interface type %d not supported\n", (int) ed);
    }

    return -1;
  }

}

/**
 * \endcond
 */
