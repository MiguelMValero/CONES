/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

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
#include <cmath>
#include "pdm_part_to_part.h"
#include "pdm_array.h"
#include "pdm_logging.h"

#include "spatialInterpIdentity.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

/**
 * \cond
 */

/**
 * NOTES
 * Use cpl comm or union comm?
 * + : use fields utilities (visu, interlaced/interleaved, ...)
 * - : require "fake" mesh or user_tgt_cloud (either way, only one per coupling)
 * - : limited to 'double'
 */

namespace cwipi {
  SpatialInterpIdentity::SpatialInterpIdentity() = default;

  SpatialInterpIdentity::~SpatialInterpIdentity
  (
   )
  {
    if (_src_to_tgt_idx != NULL) {
      for (int i = 0; i < _nPart; i++) {
        if (_src_to_tgt_idx[i] != NULL) {
          free(_src_to_tgt_idx[i]);
        }
      }
      free(_src_to_tgt_idx);
      _src_to_tgt_idx = NULL;
    }
  }

  /**
    *
    * \brief SpatialInterp location Init.
    *
    */

  void
  SpatialInterpIdentity::init
  (
    Coupling                   *coupling,
    CWP_Dof_location_t          localCodeDofLocation,
    CWP_Dof_location_t          cplCodeDofLocation,
    SpatialInterpExchDirection  exchDirection
  )
  {
    SpatialInterp::init(coupling,
                        localCodeDofLocation,
                        cplCodeDofLocation,
                        exchDirection);

    _interpolation_time = CWP_SPATIAL_INTERP_AT_RECV;

    //
    // Data for PDM_part_to_part_t
    _src_to_tgt_idx = NULL;

    if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
      _src_to_tgt_idx = (int **) malloc (sizeof(int *) * _nPart);
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
          _src_gnum  [i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet(i_part);
          _src_n_gnum[i_part] = _mesh->getPartNElts(i_part);
        }
        else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
          _src_gnum  [i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum(i_part);
          _src_n_gnum[i_part] = _mesh->getPartNVertex(i_part);
        }
        else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
          _src_gnum  [i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet(i_part);
          _src_n_gnum[i_part] = _cpl->userTargetNGet(i_part);
        }

        _src_to_tgt_idx[i_part] = PDM_array_new_idx_from_const_stride_int(1,
                                                                          _src_n_gnum[i_part]);
      }
    }
    else {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
          _tgt_gnum  [i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet(i_part);
          _tgt_n_gnum[i_part] = _mesh->getPartNElts(i_part);
        }
        else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
          _tgt_gnum  [i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum(i_part);
          _tgt_n_gnum[i_part] = _mesh->getPartNVertex(i_part);
        }
        else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
          _tgt_gnum  [i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet(i_part);
          _tgt_n_gnum[i_part] = _cpl->userTargetNGet(i_part);
        }
      }
    }
  }


  void SpatialInterpIdentity::clear()
  {
    SpatialInterp::clear();

    if (_src_to_tgt_idx != NULL) {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_src_to_tgt_idx[i_part] != NULL) {
          free(_src_to_tgt_idx[i_part]);
        }
      }

      free(_src_to_tgt_idx);
      _src_to_tgt_idx = NULL;
    }
  }


  void SpatialInterpIdentity::weightsCompute()
  {

    if (!_coupledCodeProperties->localCodeIs()) {
      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
        for (int i_part = 0; i_part < _nPart; i_part++) {
          _n_involved_sources_tgt[i_part] = _src_n_gnum[i_part];
          _involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * _n_involved_sources_tgt[i_part]);
          for (int i = 0 ; i < _src_n_gnum[i_part] ; ++i) {
            _involved_sources_tgt[i_part][i] = i + 1;
          }
        }

        _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) _src_gnum,
                                        (const int          *) _src_n_gnum,
                                        _nPart,
                                        (const PDM_g_num_t **) NULL,
                                        (const int          *) NULL,
                                        0,
                                        (const int         **) _src_to_tgt_idx,
                                        (const PDM_g_num_t **) _src_gnum,
                                        _pdmCplComm);
      }
      else {
        for (int i_part = 0; i_part < _nPart; i_part++) {
          _n_computed_tgt  [i_part] = _tgt_n_gnum[i_part];
          _n_uncomputed_tgt[i_part] = 0;
          _computed_tgt[i_part] = (int *) malloc(sizeof(int) * _n_computed_tgt[i_part]);
          for (int i = 0; i < _n_computed_tgt[i_part]; i++) {
            _computed_tgt[i_part][i] = i + 1;
          }
        }

        _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) NULL,
                                        (const int          *) NULL,
                                        0,
                                        (const PDM_g_num_t **) _tgt_gnum,
                                        (const int          *) _tgt_n_gnum,
                                        _nPart,
                                        (const int         **) NULL,
                                        (const PDM_g_num_t **) NULL,
                                        _pdmCplComm);
      }
    }
    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        SpatialInterpIdentity *cpl_spatial_interp;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
          cpl_spatial_interp =
          dynamic_cast <SpatialInterpIdentity *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }
        else {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
          cpl_spatial_interp =
          dynamic_cast <SpatialInterpIdentity *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            _n_involved_sources_tgt[i_part] = _src_n_gnum[i_part];
            _involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * _n_involved_sources_tgt[i_part]);
            for (int i = 0; i < _src_n_gnum[i_part]; ++i) {
              _involved_sources_tgt[i_part][i] = i + 1;
            }
          }

          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            cpl_spatial_interp->_n_computed_tgt  [i_part] = cpl_spatial_interp->_tgt_n_gnum[i_part];
            cpl_spatial_interp->_n_uncomputed_tgt[i_part] = 0;
            cpl_spatial_interp->_computed_tgt[i_part] = (int *) malloc(sizeof(int) * cpl_spatial_interp->_n_computed_tgt[i_part]);
            for (int i = 0; i < cpl_spatial_interp->_n_computed_tgt[i_part]; i++) {
              cpl_spatial_interp->_computed_tgt[i_part][i] = i + 1;
            }
          }

          _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) _src_gnum,
                                          (const int          *) _src_n_gnum,
                                          _nPart,
                                          (const PDM_g_num_t **) cpl_spatial_interp->_tgt_gnum,
                                          (const int          *) cpl_spatial_interp->_tgt_n_gnum,
                                          _cplNPart,
                                          (const int         **) _src_to_tgt_idx,
                                          (const PDM_g_num_t **) _src_gnum,
                                          _pdmCplComm);
        }
        else {
          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            cpl_spatial_interp->_n_involved_sources_tgt[i_part] = cpl_spatial_interp->_src_n_gnum[i_part];
            cpl_spatial_interp->_involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * cpl_spatial_interp->_n_involved_sources_tgt[i_part]);
            for (int i = 0; i < cpl_spatial_interp->_src_n_gnum[i_part]; ++i) {
              cpl_spatial_interp->_involved_sources_tgt[i_part][i] = i + 1;
            }
          }

          for (int i_part = 0; i_part < _nPart; i_part++) {
            _n_computed_tgt  [i_part] = _tgt_n_gnum[i_part];
            _n_uncomputed_tgt[i_part] = 0;
            _computed_tgt[i_part] = (int *) malloc(sizeof(int) * _n_computed_tgt[i_part]);
            for (int i = 0; i < _n_computed_tgt[i_part]; i++) {
              _computed_tgt[i_part][i] = i + 1;
            }
          }

          _ptsp = PDM_part_to_part_create((const PDM_g_num_t **) cpl_spatial_interp->_src_gnum,
                                          (const int          *) cpl_spatial_interp->_src_n_gnum,
                                          _cplNPart,
                                          (const PDM_g_num_t **) _tgt_gnum,
                                          (const int          *) _tgt_n_gnum,
                                          _nPart,
                                          (const int         **) cpl_spatial_interp->_src_to_tgt_idx,
                                          (const PDM_g_num_t **) cpl_spatial_interp->_src_gnum,
                                          _pdmCplComm);
        }

        cpl_spatial_interp->_ptsp = _ptsp;
      }
    }
  }


  void SpatialInterpIdentity::interpolate
  (
   Field   *referenceField,
   double **buffer
   )
  {
    int nComponent = referenceField->nComponentGet();
    const CWP_Field_storage_t storage = referenceField->storageTypeGet();

    int  *n_ref = NULL;
    int **ref   = NULL;
    PDM_part_to_part_ref_lnum2_get(_ptsp,
                                   &n_ref,
                                   &ref);

    int         **come_from_idx = NULL;
    PDM_g_num_t **come_from     = NULL;
    PDM_part_to_part_gnum1_come_from_get(_ptsp,
                                         &come_from_idx,
                                         &come_from);

    for (int i_part = 0; i_part < _nPart; i_part++) {

      // if (n_ref[i_part] != _tgt_n_gnum[i_part]) {
      //   PDM_error(__FILE__, __LINE__, 0, "Invalid Part Data (some tgt are not referenced)\n");
      // }

      double *referenceData = (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_TARGET);

      /*
       in case of multiple 'come from', what do we do?
       default : first?, mean?, ...
       user 'interpolation'?
      */
      if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
        for (int j = 0; j < nComponent; j++) {
          for (int i = 0; i < n_ref[i_part]; i++) {
            int idx = come_from_idx[i_part][i];
            referenceData[n_ref[i_part]*j + i] = buffer[i_part][come_from_idx[i_part][i]*j + idx];
          }
        }
      }
      else { // if (storage == CWP_FIELD_STORAGE_INTERLACED) {
        for (int i = 0; i < n_ref[i_part]; i++) {
          int idx = come_from_idx[i_part][i];
          for (int j = 0; j < nComponent; j++) {
            referenceData[nComponent*i + j] = buffer[i_part][nComponent*idx + j];
          }
        }
      }

    }
  }


}

/**
 * \endcond
 */
