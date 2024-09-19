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
#include "spatialInterp.hxx"
#include "spatialInterp_i.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include <mpi.h>

#include "pdm_mesh_nodal.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_gnum.h"
#include "pdm_timer.h"
#include "pdm_gnum_location.h"
#include "pdm_binary_search.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "cwp.h"
#include "cwp_priv.h"
#include <limits>
#include <vector>
#include <algorithm>
#include <cmath>


/**
 * \cond
 */

namespace cwipi {

  SpatialInterp::SpatialInterp()
  {
  }


  SpatialInterp::~SpatialInterp()
  {
    if (!_coupledCodeProperties->localCodeIs()) {
      if (_ptsp != nullptr) {
        PDM_part_to_part_free (_ptsp);
        _ptsp = nullptr;
      }
    }
    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        PDM_part_to_part_free (_ptsp);
        _ptsp = nullptr;
      }
    }

    free ( _n_elt_weights);
    free ( _n_computed_tgt);
    free ( _n_uncomputed_tgt);
    free ( _n_involved_sources_tgt);
    free ( _src_n_gnum);
    free ( _tgt_n_gnum);

    for (int i = 0; i < _nPart; i++) {
      if (_weights_idx[i] != NULL) {
        free (_weights_idx[i]);
      }
      if (_weights[i] != NULL) {
        free (_weights[i]);
      }
      if (_computed_tgt[i] != NULL) {
        free (_computed_tgt[i]);
      }
      if (_involved_sources_tgt[i] != NULL) {
        free(_involved_sources_tgt[i]);
      }
      if (_uncomputed_tgt[i] != NULL) {
        free (_uncomputed_tgt[i]);
      }
    }
    free ( _weights_idx);
    free ( _weights);

    free ( _computed_tgt);
    free ( _involved_sources_tgt);

    free ( _uncomputed_tgt);

    free ( _src_gnum);
    free ( _tgt_gnum);
  }


  void 
  SpatialInterp::init(
    Coupling                  *coupling, 
    CWP_Dof_location_t         localCodeDofLocation,
    CWP_Dof_location_t         cplCodeDofLocation,
    SpatialInterpExchDirection exchDirection 
  )
  {
    _cpl                    = coupling;
    // _visu                   = coupling->visuGet();
    _mesh                   = coupling->meshGet();

    _localCodeDofLocation   = localCodeDofLocation;
    _coupledCodeDofLocation = cplCodeDofLocation;

    _localCodeProperties    = _cpl->localCodePropertiesGet();
    _coupledCodeProperties  = _cpl->coupledCodePropertiesGet();

    _cplComm                = _cpl->communicationGet()->cplCommGet();
    _pdmCplComm             = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_cplComm));

    _unionComm              = _cpl->communicationGet()->unionCommGet();
    _pdmUnionComm           = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_unionComm));

    _nPart                  = _mesh->getNPart();

    _rootRankUnionComm      = _cpl->communicationGet()->unionCommLocCodeRootRanksGet();
    _cplRootRankUnionComm   = _cpl->communicationGet()->unionCommCplCodeRootRanksGet();

    _rootRankCplComm        = _cpl->communicationGet()->cplCommLocCodeRootRanksGet();
    _cplRootRankCplComm     = _cpl->communicationGet()->cplCommCplCodeRootRanksGet();

    _exchDirection          = exchDirection;

    //_connectableRanks_cpl = _cpl->communicationGet()->cplCommCplRanksGet();
    //_connectableRanks     = _cpl->communicationGet()->cplCommLocRanksGet();

    //_id     = _localCodeProperties  ->idGet();
    //_id_cpl = _coupledCodeProperties->idGet();

    _send_buffer.resize(_cpl->fieldsGet()->size());  /*!< Send buffer (size = n_field) */
    _recv_buffer.resize(_cpl->fieldsGet()->size());  /*!< Recv buffer (size = n_field) */

    _send_request.resize(_cpl->fieldsGet()->size()); /*!< Send request (size = n_field) */
    _recv_request.resize(_cpl->fieldsGet()->size()); /*!< Recv request (size = n_field) */


    for (size_t i = 0; i < _cpl->fieldsGet()->size(); i++) {
      _send_buffer[i]  = NULL;
      _recv_buffer[i]  = NULL;
      _send_request[i] = -1;
      _recv_request[i] = -1;
    }

    _cplNPart = _cpl->cplNPartGet();
    _nPart    = _cpl->nPartGet();

    _ptsp = NULL;

    if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
      _src_n_gnum =  (int *) malloc (sizeof(int) * (_nPart));
      _src_gnum =  (const PDM_g_num_t **) malloc (sizeof(PDM_g_num_t) * _nPart);

      _tgt_n_gnum =  (int *) malloc (sizeof(int) * (_cplNPart));
      _tgt_gnum = (const PDM_g_num_t **) malloc (sizeof(PDM_g_num_t) * _cplNPart);

      for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
        _src_n_gnum[i_part] = 0;
        _src_gnum[i_part] = nullptr;
      }

      for (int i_part = 0 ; i_part < _cplNPart ; i_part++) { 
        _tgt_n_gnum[i_part] = 0;
        _tgt_gnum[i_part] = nullptr;
      }
    }
    else {
      _src_n_gnum =  (int *) malloc (sizeof(int) * (_cplNPart));
      _src_gnum = (const PDM_g_num_t **) malloc (sizeof(PDM_g_num_t) * _cplNPart);

      _tgt_n_gnum =  (int *) malloc (sizeof(int) * (_nPart));
      _tgt_gnum = (const PDM_g_num_t **) malloc (sizeof(PDM_g_num_t) * _nPart);

      for (int i_part = 0 ; i_part < _cplNPart ; i_part++) { 
        _src_n_gnum[i_part] = 0;
        _src_gnum[i_part] = nullptr;
      }

      for (int i_part = 0 ; i_part < _nPart ; i_part++) { 
        _tgt_n_gnum[i_part] = 0;
        _tgt_gnum[i_part] = nullptr;
      }
    }


    _n_elt_weights =  (int *) malloc (sizeof(int) * (_nPart));
    _weights_idx =  (int **) malloc (sizeof(int *) * (_nPart));
    _weights =  (double **) malloc (sizeof(double *) * (_nPart));

    _n_computed_tgt =  (int *) malloc (sizeof(int) * (_nPart));
    _computed_tgt =  (int **) malloc (sizeof(int *) * (_nPart));

    _n_involved_sources_tgt =  (int *) malloc (sizeof(int) * (_nPart));
    _involved_sources_tgt =  (int **) malloc (sizeof(int *) * (_nPart));

    _n_uncomputed_tgt =  (int *) malloc (sizeof(int) * (_nPart));
    _uncomputed_tgt =  (int **) malloc (sizeof(int *) * (_nPart));

    for (int i = 0; i < _nPart; i++) {
      _weights_idx[i] = NULL;
      _weights[i] = NULL;
      _computed_tgt[i] = NULL;
      _uncomputed_tgt[i] = NULL;
      _involved_sources_tgt[i] = NULL;

      _n_uncomputed_tgt[i] = 0;
      _n_computed_tgt[i] = 0;
      _n_involved_sources_tgt[i] = 0;
      _n_elt_weights[i] = 0;
    }

  }


  void
  SpatialInterp::clear()
  {
    if (_src_n_gnum != NULL) {
      free(_src_n_gnum);
      _src_n_gnum = NULL;
    }

    if (_src_gnum != NULL) {
      free(_src_gnum);
      _src_gnum = NULL;
    }

    if (_tgt_n_gnum != NULL) {
      free(_tgt_n_gnum);
      _tgt_n_gnum = NULL;
    }

    if (_tgt_gnum != NULL) {
      free(_tgt_gnum);
      _tgt_gnum = NULL;
    }

    if (_n_elt_weights != NULL) {
      free(_n_elt_weights);
      _n_elt_weights = NULL;
    }

    if (_weights_idx != NULL) {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_weights_idx[i_part] != NULL) {
          free(_weights_idx[i_part]);
          _weights_idx[i_part] = NULL;
        }
      }
      free(_weights_idx);
      _weights_idx = NULL;
    }
    if (_weights != NULL) {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_weights[i_part] != NULL) {
          free(_weights[i_part]);
          _weights[i_part] = NULL;
        }
      }
      free(_weights);
      _weights = NULL;
    }

    if (_n_computed_tgt != NULL) {
      free(_n_computed_tgt);
      _n_computed_tgt = NULL;
    }
    if (_computed_tgt != NULL) {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_computed_tgt[i_part] != NULL) {
          free(_computed_tgt[i_part]);
          _computed_tgt[i_part] = NULL;
        }
      }
      free(_computed_tgt);
      _computed_tgt = NULL;
    }

    if (_n_involved_sources_tgt != NULL) {
      free(_n_involved_sources_tgt);
      _n_involved_sources_tgt = NULL;
    }
    if (_involved_sources_tgt != NULL) {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_involved_sources_tgt[i_part] != NULL) {
          free(_involved_sources_tgt[i_part]);
          _involved_sources_tgt[i_part] = NULL;
        }
      }
      free(_involved_sources_tgt);
      _involved_sources_tgt = NULL;
    }

    if (_n_uncomputed_tgt != NULL) {
      free(_n_uncomputed_tgt);
      _n_uncomputed_tgt = NULL;
    }
    if (_uncomputed_tgt != NULL) {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_uncomputed_tgt[i_part] != NULL) {
          free(_uncomputed_tgt[i_part]);
          _uncomputed_tgt[i_part] = NULL;
        }
      }
      free(_uncomputed_tgt);
      _uncomputed_tgt = NULL;
    }


    if (_ptsp != NULL) {
      PDM_part_to_part_free(_ptsp);
      _ptsp = NULL;
    }

    if (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
      SpatialInterp *cpl_spatial_interp;

      cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
        cpl_spatial_interp = cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];
      }
      else {
        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
        cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];
      }

      cpl_spatial_interp->_ptsp = NULL;
    }
  }


/***************************************************************************/
/***************************************************************************/

  void SpatialInterp::issend(Field* referenceField) {

    if (referenceField->is_send_yet_get() == 1) {
      PDM_error(__FILE__, __LINE__, 0,
                "The field has already been exchanged for the current time step. "
                "Each field can only be exchanged once within a time step "
                "(framed by calls to CWP_time_step_beg and CWP_time_step_end functions)\n");
    }

    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId            = referenceField->fieldIDIntGet();
      const CWP_Type_t data_type = referenceField->dataTypeGet();
      CWP_UNUSED (data_type);
      const size_t s_data        = sizeof(double);
      const int stride           = referenceField->nComponentGet();
      const CWP_Field_storage_t storage = referenceField->storageTypeGet();
      PDM_stride_t pdm_storage;
      if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
        pdm_storage = PDM_STRIDE_CST_INTERLEAVED;
      } else {
        pdm_storage = PDM_STRIDE_CST_INTERLACED;
      }

      int           *n_elt1;
      int          **selected_part2_idx;
      PDM_g_num_t  **selected_part2;

      PDM_part_to_part_part1_to_part2_get (_ptsp,
                                           &n_elt1,
                                           &selected_part2_idx,
                                           &selected_part2);

      _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);

      // A n_elt1

      for (int i = 0; i < _nPart; i++) {
        _send_buffer[intId][i] = 
          (double *) malloc(sizeof(double) * stride * selected_part2_idx[i][n_elt1[i]]);
      }

      if (_interpolation_time == CWP_SPATIAL_INTERP_AT_SEND) {
        interpolate (referenceField, _send_buffer[intId]);
      }
      else {
        // send source field (copy?)
        for (int i = 0; i < _nPart; i++) {
          double *referenceData = (double *) referenceField->dataGet(i, CWP_FIELD_MAP_SOURCE);

          if (storage == CWP_FIELD_STORAGE_INTERLACED) {
            for (int j = 0; j < n_elt1[i]; j++) {
              for (int k = selected_part2_idx[i][j]; k < selected_part2_idx[i][j+1]; k++) {
                for (int l = 0; l < stride; l++) {
                  _send_buffer[intId][i][stride*k + l] = referenceData[stride*j + l];
                }
              }
            }
          }
          else {
            for (int l = 0; l < stride; l++) {
              for (int j = 0; j < n_elt1[i]; j++) {
                for (int k = selected_part2_idx[i][j]; k < selected_part2_idx[i][j+1]; k++) {
                  _send_buffer[intId][i][selected_part2_idx[i][n_elt1[i]]*l + k] = referenceData[n_elt1[i]*l + j];
                }
              }
            }
          }
        }
      }


      // Fake receive
      PDM_part_to_part_iexch(_ptsp,
                             PDM_MPI_COMM_KIND_P2P,
                             pdm_storage,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                             stride,
                             s_data,
                             NULL,
             (const void **) _send_buffer[intId],
                             NULL,
                  (void ***) &_recv_buffer[intId],
                             &(_send_request[intId]));
      _recv_request[intId] = _send_request[intId];


    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId            = referenceField->fieldIDIntGet();
        const CWP_Type_t data_type = referenceField->dataTypeGet();
        CWP_UNUSED(data_type);
        const size_t s_data        = sizeof(double);
        const int stride           = referenceField->nComponentGet();
        const CWP_Field_storage_t storage = referenceField->storageTypeGet();
        PDM_stride_t pdm_storage;
        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          pdm_storage = PDM_STRIDE_CST_INTERLEAVED;
        } else {
          pdm_storage = PDM_STRIDE_CST_INTERLACED;
        }

        int           *n_elt1;
        int          **selected_part2_idx;
        PDM_g_num_t  **selected_part2;

        PDM_part_to_part_part1_to_part2_get (_ptsp,
                                             &n_elt1,
                                             &selected_part2_idx,
                                             &selected_part2);

        _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);

        for (int i = 0; i < _nPart; i++) {
          _send_buffer[intId][i] = (double *) malloc(sizeof(double) * stride  * selected_part2_idx[i][n_elt1[i]]);
        }

        if (_interpolation_time == CWP_SPATIAL_INTERP_AT_SEND) {
          interpolate (referenceField, _send_buffer[intId]);
        }
        else {
          // send source field (copy?)
          for (int i = 0; i < _nPart; i++) {
            double *referenceData = (double *) referenceField->dataGet(i, CWP_FIELD_MAP_SOURCE);

            if (storage == CWP_FIELD_STORAGE_INTERLACED) {
              for (int j = 0; j < n_elt1[i]; j++) {
                for (int k = selected_part2_idx[i][j]; k < selected_part2_idx[i][j+1]; k++) {
                  for (int l = 0; l < stride; l++) {
                    _send_buffer[intId][i][stride*k + l] = referenceData[stride*j + l];
                  }
                }
              }
            }
            else {
              for (int l = 0; l < stride; l++) {
                for (int j = 0; j < n_elt1[i]; j++) {
                  for (int k = selected_part2_idx[i][j]; k < selected_part2_idx[i][j+1]; k++) {
                    _send_buffer[intId][i][selected_part2_idx[i][n_elt1[i]]*l + k] = referenceData[n_elt1[i]*l + j];
                  }
                }
              }
            }
          }
        }


        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 
        SpatialInterp *cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();
        cpl_spatial_interp->_send_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);
        // cpl_spatial_interp->_recv_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);

        int  *n_ref_gnum2;
        int **ref_gnum2;

        PDM_part_to_part_ref_lnum2_get (_ptsp,
                                        &n_ref_gnum2,
                                        &ref_gnum2);

        int          **gnum1_come_from_idx;
        PDM_g_num_t  **gnum1_come_from;

        PDM_part_to_part_gnum1_come_from_get (_ptsp,
                                              &gnum1_come_from_idx,
                                              &gnum1_come_from);

        for (int i = 0; i < _cplNPart; i++) {
          cpl_spatial_interp->_send_buffer[cpl_intId][i] = nullptr;
        }
        PDM_part_to_part_iexch(_ptsp,
                               PDM_MPI_COMM_KIND_P2P,
                               pdm_storage,
                               PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                               stride,
                               s_data,
                               NULL,
               (const void **) _send_buffer[intId],
                               NULL,
                    (void ***) &cpl_spatial_interp->_recv_buffer[cpl_intId],
                               &(_send_request[intId]));
        cpl_spatial_interp->_recv_request[cpl_intId] = _send_request[intId];
      }
    }
  }


  void SpatialInterp::waitIssend(Field* referenceField) {

    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId            = referenceField->fieldIDIntGet();

      PDM_part_to_part_iexch_wait (_ptsp, _send_request[intId]);

      int                _n_part1;
      int                _n_part2;
      int               *_n_elt1;
      int               *_n_elt2;
      PDM_part_to_part_n_part_and_n_elt_get(_ptsp,
                                            &_n_part1,
                                            &_n_part2,
                                            &_n_elt1,
                                            &_n_elt2);

      if (_cpl->commTypeGet() == CWP_COMM_PAR_WITHOUT_PART &&
          referenceField->involvedSrcBcastIsEnabled()) {
        const int localRootRank = _localCodeProperties->rootRankGet();
        const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();
        int globalRank;
        MPI_Comm_rank(globalComm, &globalRank);

        MPI_Group globalGroup;
        MPI_Comm_group(_localCodeProperties->globalCommGet(), &globalGroup);
        int localRootRankInIntraConnectableComm;
        MPI_Group_translate_ranks(globalGroup, 1, &localRootRank,
                                  _localCodeProperties->connectableGroupGet(), &localRootRankInIntraConnectableComm);

        int rankInIntraConnectableComm;
        MPI_Comm_rank(_localCodeProperties->connectableCommGet(), &rankInIntraConnectableComm);

        for (int i = 0; i < _nPart; i++) {
          // Broadcast involved src
          MPI_Bcast((void *) &_n_involved_sources_tgt[i],
                    1,
                    MPI_INT,
                    localRootRankInIntraConnectableComm,
                    _localCodeProperties->connectableCommGet());

          if (rankInIntraConnectableComm != localRootRankInIntraConnectableComm) {
            if (_involved_sources_tgt[i] != NULL) {
              free(_involved_sources_tgt[i]);
            }
            _involved_sources_tgt[i] = (int *) malloc(sizeof(int) * _n_involved_sources_tgt[i]);
          }
          MPI_Bcast((void *) _involved_sources_tgt[i],
                    _n_involved_sources_tgt[i],
                    MPI_INT,
                    localRootRankInIntraConnectableComm,
                    _localCodeProperties->connectableCommGet());
        }
      }

      if (_send_buffer[intId] != NULL) {
        for (int i = 0; i < _n_part1; i++) {
          if (_send_buffer[intId][i] != NULL) {
            free (_send_buffer[intId][i]);
            _send_buffer[intId][i] = NULL;
          }
        }
        free (_send_buffer[intId]);
        _send_buffer[intId] = NULL;
      }

      if (_recv_buffer[intId] != NULL) {
        for (int i = 0; i < _n_part2; i++) {
          if (_recv_buffer[intId][i] != NULL) {
            free (_recv_buffer[intId][i]);
            _recv_buffer[intId][i] = NULL;
          }
        }
        free (_recv_buffer[intId]);
        _recv_buffer[intId] = NULL;
      }

      PDM_writer_t* writer =  _cpl->writerGet();

      if (writer != nullptr) {
        if ((_cpl->NStepGet() % _cpl->freqWriterGet()) == 0) {
          referenceField->write(CWP_FIELD_EXCH_SEND);
        }
      }

    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId            = referenceField->fieldIDIntGet();

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet(); 
        SpatialInterp *cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();

        PDM_part_to_part_iexch_wait(_ptsp, cpl_spatial_interp->_recv_request[cpl_intId]);

        if (_cpl->commTypeGet() == CWP_COMM_PAR_WITHOUT_PART &&
          referenceField->involvedSrcBcastIsEnabled()) {
          const int localRootRank = _localCodeProperties->rootRankGet();
          const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();
          int globalRank;
          MPI_Comm_rank(globalComm, &globalRank);

          MPI_Group globalGroup;
          MPI_Comm_group(_localCodeProperties->globalCommGet(), &globalGroup);
          int localRootRankInIntraConnectableComm;
          MPI_Group_translate_ranks(globalGroup, 1, &localRootRank,
                                    _localCodeProperties->connectableGroupGet(), &localRootRankInIntraConnectableComm);

          int rankInIntraConnectableComm;
          MPI_Comm_rank(_localCodeProperties->connectableCommGet(), &rankInIntraConnectableComm);

          for (int i = 0; i < _nPart; i++) {
            // Broadcast involved src
            MPI_Bcast((void *) &_n_involved_sources_tgt[i],
                      1,
                      MPI_INT,
                      localRootRankInIntraConnectableComm,
                      _localCodeProperties->connectableCommGet());

            if (rankInIntraConnectableComm != localRootRankInIntraConnectableComm) {
              if (_involved_sources_tgt[i] != NULL) {
                free(_involved_sources_tgt[i]);
              }
              _involved_sources_tgt[i] = (int *) malloc(sizeof(int) * _n_involved_sources_tgt[i]);
            }
            MPI_Bcast((void *) _involved_sources_tgt[i],
                      _n_involved_sources_tgt[i],
                      MPI_INT,
                      localRootRankInIntraConnectableComm,
                      _localCodeProperties->connectableCommGet());
          }
        }

        if (_send_buffer[intId] != NULL) {
          for (int i = 0; i < _nPart; i++) {
            if (_send_buffer[intId][i] != NULL) {
              free (_send_buffer[intId][i]);
              _send_buffer[intId][i] = NULL;
            }
          }
          free (_send_buffer[intId]);
          _send_buffer[intId] = NULL;
        }
        if (_recv_buffer[intId] != NULL) {
          for (int i = 0; i < _nPart; i++) {
            if (_recv_buffer[intId][i] != NULL) {
              free (_recv_buffer[intId][i]);
              _recv_buffer[intId][i] = NULL;
            }
          }
          free (_recv_buffer[intId]);
          _recv_buffer[intId] = NULL;
        }

        PDM_writer_t* writer =  _cpl->writerGet();

        if (writer != nullptr) {
          if ((_cpl->NStepGet() % _cpl->freqWriterGet()) == 0) {
            referenceField->write(CWP_FIELD_EXCH_SEND);
          }
        }

        const int cplRootRank = _coupledCodeProperties->rootRankGet();
        const MPI_Comm& globalComm = _coupledCodeProperties->globalCommGet();
        int globalRank;
        MPI_Comm_rank(globalComm, &globalRank);

        if (cpl_cpl.commTypeGet() == CWP_COMM_PAR_WITH_PART || globalRank == cplRootRank) {
          if (_interpolation_time == CWP_SPATIAL_INTERP_AT_RECV) {
            cpl_spatial_interp->interpolate (cpl_referenceField, cpl_spatial_interp->_recv_buffer[cpl_intId]);
          }
          else {
            int nComponent                        = cpl_referenceField->nComponentGet();
            int dataTypeSize                      = cpl_referenceField->dataTypeSizeGet();

            int  *ptp2_n_ref_gnum2;
            int **ptp2_ref_gnum2;
            PDM_part_to_part_ref_lnum2_get (_ptsp,
                                            &ptp2_n_ref_gnum2,
                                            &ptp2_ref_gnum2);

            int  *ptp2_n_unref_gnum2;
            int **ptp2_unref_gnum2;
            PDM_part_to_part_unref_lnum2_get (_ptsp,
                                              &ptp2_n_unref_gnum2,
                                              &ptp2_unref_gnum2);


            int         **ptp2_gnum1_come_from_idx;
            PDM_g_num_t **ptp2_gnum1_come_from;
            PDM_part_to_part_gnum1_come_from_get (_ptsp,
                                                  &ptp2_gnum1_come_from_idx,
                                                  &ptp2_gnum1_come_from);


            for (int i = 0; i < _cplNPart; i++) {
              double *referenceData  = (double *) cpl_referenceField->dataGet(i, CWP_FIELD_MAP_TARGET);
              for (int j = 0; j < ptp2_n_ref_gnum2[i]; j++) {
                assert ((ptp2_gnum1_come_from_idx[i][j+1] - ptp2_gnum1_come_from_idx[i][j]) == 1);
              }
              memcpy(referenceData, cpl_spatial_interp->_recv_buffer[cpl_intId][i], dataTypeSize * nComponent * ptp2_n_ref_gnum2[i]);
            }
          }
        }

        if (cpl_cpl.commTypeGet() == CWP_COMM_PAR_WITHOUT_PART) {
          /* Root rank broadcasts the referenceField inside the intraConnectableComm */
          MPI_Group globalGroup;
          MPI_Comm_group(_coupledCodeProperties->globalCommGet(), &globalGroup);
          int cplRootRankInIntraConnectableComm;
          MPI_Group_translate_ranks(globalGroup, 1, &cplRootRank,
                                    _coupledCodeProperties->connectableGroupGet(), &cplRootRankInIntraConnectableComm);

          int rankInIntraConnectableComm;
          MPI_Comm_rank(_coupledCodeProperties->connectableCommGet(), &rankInIntraConnectableComm);

          // Easier way to get 'n_elt'?
          int nComponent = cpl_referenceField->nComponentGet();
          int  n_elt;
          int  n_part1 = 0;
          int  n_part2 = 0;
          int *n_elt1  = NULL;
          int *n_elt2  = NULL;
          PDM_part_to_part_n_part_and_n_elt_get(_ptsp,
                                                &n_part1,
                                                &n_part2,
                                                &n_elt1,
                                                &n_elt2);
          for (int i = 0; i < _nPart; i++) {
            if (globalRank == cplRootRank) {
              n_elt = n_elt2[i];
            }
            MPI_Bcast((void *) &n_elt,
                      1,
                      MPI_INT,
                      cplRootRankInIntraConnectableComm,
                      _coupledCodeProperties->connectableCommGet());

            double *cpl_referenceData = (double *) cpl_referenceField->dataGet(i, CWP_FIELD_MAP_TARGET);

            MPI_Bcast((void *) cpl_referenceData,
                      n_elt * nComponent,
                      MPI_DOUBLE,
                      cplRootRankInIntraConnectableComm,
                      _coupledCodeProperties->connectableCommGet());

            if (cpl_referenceField->computedTgtBcastIsEnabled()) {
              // Broadcast (un)computed tgt as well
              MPI_Bcast((void *) &cpl_spatial_interp->_n_computed_tgt[i],
                        1,
                        MPI_INT,
                        cplRootRankInIntraConnectableComm,
                        _coupledCodeProperties->connectableCommGet());

              if (rankInIntraConnectableComm != cplRootRankInIntraConnectableComm) {
                if (cpl_spatial_interp->_computed_tgt[i] != NULL) {
                  free(cpl_spatial_interp->_computed_tgt[i]);
                }
                cpl_spatial_interp->_computed_tgt[i] = (int *) malloc(sizeof(int) * cpl_spatial_interp->_n_computed_tgt[i]);
              }
              MPI_Bcast((void *) cpl_spatial_interp->_computed_tgt[i],
                        cpl_spatial_interp->_n_computed_tgt[i],
                        MPI_INT,
                        cplRootRankInIntraConnectableComm,
                        _coupledCodeProperties->connectableCommGet());
            }
          }
        }

        if (cpl_spatial_interp->_send_buffer[cpl_intId] != NULL) {
          for (int i = 0; i < _cplNPart; i++) {
            if (cpl_spatial_interp->_send_buffer[cpl_intId][i] != NULL) {
              free (cpl_spatial_interp->_send_buffer[cpl_intId][i]);
              cpl_spatial_interp->_send_buffer[cpl_intId][i] = NULL;
            }
          }
          free (cpl_spatial_interp->_send_buffer[cpl_intId]);
          cpl_spatial_interp->_send_buffer[cpl_intId] = NULL;
        }

        if (cpl_spatial_interp->_recv_buffer[cpl_intId] != NULL) {
          for (int i = 0; i < _cplNPart; i++) {
            if (cpl_spatial_interp->_recv_buffer[cpl_intId][i] != NULL) {
              free (cpl_spatial_interp->_recv_buffer[cpl_intId][i]);
              cpl_spatial_interp->_recv_buffer[cpl_intId][i] = NULL;
            }
          }
          free (cpl_spatial_interp->_recv_buffer[cpl_intId]);
          cpl_spatial_interp->_recv_buffer[cpl_intId] = NULL;
        }

        if (writer != nullptr) {
          if ((_cpl->NStepGet() % _cpl->freqWriterGet()) == 0) {
            cpl_referenceField->write(CWP_FIELD_EXCH_RECV);
          }
        }
      } // end if local code works
    } // end if joint

    // set field to sent
    referenceField->is_send_yet_set(1);
  }


  void SpatialInterp::irecv(Field *referenceField) {

    if (referenceField->is_recv_yet_get()) {
      PDM_error(__FILE__, __LINE__, 0,
                "The field has already been exchanged for the current time step. "
                "Each field can only be exchanged once within a time step "
                "(framed by calls to CWP_time_step_beg and CWP_time_step_end functions)\n");
    }

    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId            = referenceField->fieldIDIntGet();
      const CWP_Type_t data_type = referenceField->dataTypeGet();
      CWP_UNUSED(data_type);
      const size_t s_data        = sizeof(double);
      const int stride           = referenceField->nComponentGet();
      const CWP_Field_storage_t storage = referenceField->storageTypeGet();
      PDM_stride_t pdm_storage;
      if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
        pdm_storage = PDM_STRIDE_CST_INTERLEAVED;
      } else {
        pdm_storage = PDM_STRIDE_CST_INTERLACED;
      }

      int  *n_ref_gnum2;
      int **ref_gnum2;

      PDM_part_to_part_ref_lnum2_get (_ptsp,
                                      &n_ref_gnum2,
                                      &ref_gnum2);

      int          **gnum1_come_from_idx;
      PDM_g_num_t  **gnum1_come_from;

      PDM_part_to_part_gnum1_come_from_get (_ptsp,
                                            &gnum1_come_from_idx,
                                            &gnum1_come_from);

      _send_buffer[intId] = (double **) malloc(sizeof(double *) * _nPart);

      for (int i = 0; i < _nPart; i++) {
        _send_buffer[intId][i] = nullptr;
      }

      // Fake receive
      PDM_part_to_part_iexch(_ptsp,
                             PDM_MPI_COMM_KIND_P2P,
                             pdm_storage,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                             stride,
                             s_data,
                             NULL,
             (const void **) _send_buffer[intId],
                             NULL,
                  (void ***) &_recv_buffer[intId],
                             &(_send_request[intId]));
      _recv_request[intId] = _send_request[intId];
    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId            = referenceField->fieldIDIntGet();
        const CWP_Type_t data_type = referenceField->dataTypeGet();
        CWP_UNUSED(data_type);
        const size_t s_data        = sizeof(double);
        const int stride           = referenceField->nComponentGet();
        const CWP_Field_storage_t storage = referenceField->storageTypeGet();
        PDM_stride_t pdm_storage;
        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          pdm_storage = PDM_STRIDE_CST_INTERLEAVED;
        } else {
          pdm_storage = PDM_STRIDE_CST_INTERLACED;
        }

        int  *n_ref_gnum2;
        int **ref_gnum2;

        PDM_part_to_part_ref_lnum2_get (_ptsp,
                                        &n_ref_gnum2,
                                        &ref_gnum2);

        int          **gnum1_come_from_idx;
        PDM_g_num_t  **gnum1_come_from;
       
        PDM_part_to_part_gnum1_come_from_get (_ptsp,
                                              &gnum1_come_from_idx,
                                              &gnum1_come_from);


        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet(); 
        SpatialInterp *cpl_spatial_interp = cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();
        cpl_spatial_interp->_send_buffer[cpl_intId] = (double **) malloc(sizeof(double *) * _cplNPart);

        int           *cpl_n_elt1;
        int          **cpl_selected_part2_idx;
        PDM_g_num_t  **cpl_selected_part2;

        PDM_part_to_part_part1_to_part2_get (_ptsp,
                                             &cpl_n_elt1,
                                             &cpl_selected_part2_idx,
                                             &cpl_selected_part2);

        for (int i = 0; i < _cplNPart; i++) {
          cpl_spatial_interp->_send_buffer[cpl_intId][i] = (double *) malloc(sizeof(double) * stride * cpl_selected_part2_idx[i][cpl_n_elt1[i]]);
        }

        if (_interpolation_time == CWP_SPATIAL_INTERP_AT_SEND) {
          cpl_spatial_interp->interpolate (cpl_referenceField, cpl_spatial_interp->_send_buffer[cpl_intId]);
        }
        else {
          // send source field (copy?)
          for (int i = 0; i < _cplNPart; i++) {
            double *cpl_referenceData = (double *) cpl_referenceField->dataGet(i, CWP_FIELD_MAP_SOURCE);

            if (storage == CWP_FIELD_STORAGE_INTERLACED) {
              for (int j = 0; j < cpl_n_elt1[i]; j++) {
                for (int k = cpl_selected_part2_idx[i][j]; k < cpl_selected_part2_idx[i][j+1]; k++) {
                  for (int l = 0; l < stride; l++) {
                    cpl_spatial_interp->_send_buffer[cpl_intId][i][stride*k + l] = cpl_referenceData[stride*j + l];
                  }
                }
              }
            }
            else {
              for (int l = 0; l < stride; l++) {
                for (int j = 0; j < cpl_n_elt1[i]; j++) {
                  for (int k = cpl_selected_part2_idx[i][j]; k < cpl_selected_part2_idx[i][j+1]; k++) {
                    cpl_spatial_interp->_send_buffer[cpl_intId][i][cpl_selected_part2_idx[i][cpl_n_elt1[i]]*l + k] = cpl_referenceData[cpl_n_elt1[i]*l + j];
                  }
                }
              }
            }
          }
        }

        PDM_part_to_part_iexch(_ptsp,
                               PDM_MPI_COMM_KIND_P2P,
                               pdm_storage,
                               PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                               stride,
                               s_data,
                               NULL,
               (const void **) cpl_spatial_interp->_send_buffer[cpl_intId],
                               NULL,
                    (void ***) &_recv_buffer[intId],
                               &(cpl_spatial_interp->_send_request[cpl_intId]));
        _recv_request[intId] = cpl_spatial_interp->_send_request[cpl_intId];
      }
    }
  }


  void SpatialInterp::waitIrecv(Field* referenceField) {
    if (!_coupledCodeProperties->localCodeIs()) {

      const int intId = referenceField->fieldIDIntGet();


      PDM_part_to_part_iexch_wait (_ptsp, _send_request[intId]);

      const int localRootRank = _localCodeProperties->rootRankGet();
      const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();
      int globalRank;
      MPI_Comm_rank(globalComm, &globalRank);

      if (_cpl->commTypeGet() == CWP_COMM_PAR_WITH_PART || globalRank == localRootRank) {
        if (_interpolation_time == CWP_SPATIAL_INTERP_AT_RECV ) {
          interpolate (referenceField, _recv_buffer[intId]);
        }
        else {
          int nComponent                        = referenceField->nComponentGet();
          int dataTypeSize                      = referenceField->dataTypeSizeGet();

          int  *ptp2_n_ref_gnum2;
          int **ptp2_ref_gnum2;
          PDM_part_to_part_ref_lnum2_get (_ptsp,
                                          &ptp2_n_ref_gnum2,
                                          &ptp2_ref_gnum2);

          int  *ptp2_n_unref_gnum2;
          int **ptp2_unref_gnum2;
          PDM_part_to_part_unref_lnum2_get (_ptsp,
                                            &ptp2_n_unref_gnum2,
                                            &ptp2_unref_gnum2);


          int         **ptp2_gnum1_come_from_idx;
          PDM_g_num_t **ptp2_gnum1_come_from;
          PDM_part_to_part_gnum1_come_from_get (_ptsp,
                                                &ptp2_gnum1_come_from_idx,
                                                &ptp2_gnum1_come_from);


          for (int i = 0; i < _nPart; i++) {
            double *referenceData  = (double *) referenceField->dataGet(i, CWP_FIELD_MAP_TARGET);
            for (int j = 0; j < ptp2_n_ref_gnum2[i]; j++) {
              assert ((ptp2_gnum1_come_from_idx[i][j+1] - ptp2_gnum1_come_from_idx[i][j]) <= 1);
            }
            memcpy(referenceData, _recv_buffer[intId][i], dataTypeSize * nComponent * ptp2_n_ref_gnum2[i]);
          }
        }
      }

      if (_cpl->commTypeGet() == CWP_COMM_PAR_WITHOUT_PART) {
        /* Root rank broadcasts the referenceField inside the intraConnectableComm */
        MPI_Group globalGroup;
        MPI_Comm_group(_localCodeProperties->globalCommGet(), &globalGroup);
        int localRootRankInIntraConnectableComm;
        MPI_Group_translate_ranks(globalGroup, 1, &localRootRank,
                                  _localCodeProperties->connectableGroupGet(), &localRootRankInIntraConnectableComm);

        int rankInIntraConnectableComm;
        MPI_Comm_rank(_localCodeProperties->connectableCommGet(), &rankInIntraConnectableComm);

        // Easier way to get 'n_elt'?
        int nComponent = referenceField->nComponentGet();
        int  n_elt;
        int  n_part1 = 0;
        int  n_part2 = 0;
        int *n_elt1  = NULL;
        int *n_elt2  = NULL;
        PDM_part_to_part_n_part_and_n_elt_get(_ptsp,
                                              &n_part1,
                                              &n_part2,
                                              &n_elt1,
                                              &n_elt2);
        for (int i = 0; i < _nPart; i++) {
          if (globalRank == localRootRank) {
            n_elt = n_elt2[i];
          }
          MPI_Bcast((void *) &n_elt,
                    1,
                    MPI_INT,
                    localRootRankInIntraConnectableComm,
                    _localCodeProperties->connectableCommGet());

          double *referenceData = (double *) referenceField->dataGet(i, CWP_FIELD_MAP_TARGET);

          MPI_Bcast((void *) referenceData,
                    n_elt * nComponent,
                    MPI_DOUBLE,
                    localRootRankInIntraConnectableComm,
                    _localCodeProperties->connectableCommGet());

          if (referenceField->computedTgtBcastIsEnabled()) {
            // Broadcast (un)computed tgt as well
            MPI_Bcast((void *) &_n_computed_tgt[i],
                      1,
                      MPI_INT,
                      localRootRankInIntraConnectableComm,
                      _localCodeProperties->connectableCommGet());

            if (rankInIntraConnectableComm != localRootRankInIntraConnectableComm) {
              if (_computed_tgt[i] != NULL) {
                free(_computed_tgt[i]);
              }
              _computed_tgt[i] = (int *) malloc(sizeof(int) * _n_computed_tgt[i]);
            }
            MPI_Bcast((void *) _computed_tgt[i],
                      _n_computed_tgt[i],
                      MPI_INT,
                      localRootRankInIntraConnectableComm,
                      _localCodeProperties->connectableCommGet());
          }
        }
      }

      if (_send_buffer[intId] != NULL) {
        for (int i = 0; i < _nPart; i++) {
          if (_send_buffer[intId][i] != NULL) {
            free (_send_buffer[intId][i]);
            _send_buffer[intId][i] = NULL;
          }
        }
        free (_send_buffer[intId]);
        _send_buffer[intId] = NULL;
      }

      if (_recv_buffer[intId] != NULL) {
        for (int i = 0; i < _nPart; i++) {
          if (_recv_buffer[intId][i] != NULL) {
            free (_recv_buffer[intId][i]);
            _recv_buffer[intId][i] = NULL;
          }
        }
        free (_recv_buffer[intId]);
        _recv_buffer[intId] = NULL;
      }

      PDM_writer_t* writer =  _cpl->writerGet();

      if (writer != nullptr) {
        if ((_cpl->NStepGet() % _cpl->freqWriterGet()) == 0) {
          referenceField->write(CWP_FIELD_EXCH_RECV);
        }
      }

    } // end if disjoint

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        const int intId = referenceField->fieldIDIntGet();

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.sendSpatialInterpGet(); 
        SpatialInterp *cpl_spatial_interp = cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)];

        Field* cpl_referenceField = (*cpl_cpl.fieldsGet())[referenceField->fieldIDGet()];

        const int cpl_intId = cpl_referenceField->fieldIDIntGet();

        PDM_part_to_part_iexch_wait (_ptsp, cpl_spatial_interp->_send_request[cpl_intId]);

        const int localRootRank = _localCodeProperties->rootRankGet();
        const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();
        int globalRank;
        MPI_Comm_rank(globalComm, &globalRank);

        if (_cpl->commTypeGet() == CWP_COMM_PAR_WITH_PART || globalRank == localRootRank) {
          if (_interpolation_time == CWP_SPATIAL_INTERP_AT_RECV) {
            interpolate (referenceField, _recv_buffer[intId]);
          }
          else {
            int nComponent                        = referenceField->nComponentGet();
            int dataTypeSize                      = referenceField->dataTypeSizeGet();

            int  *ptp2_n_ref_gnum2;
            int **ptp2_ref_gnum2;
            PDM_part_to_part_ref_lnum2_get (_ptsp,
                                            &ptp2_n_ref_gnum2,
                                            &ptp2_ref_gnum2);

            int  *ptp2_n_unref_gnum2;
            int **ptp2_unref_gnum2;
            PDM_part_to_part_unref_lnum2_get (_ptsp,
                                              &ptp2_n_unref_gnum2,
                                              &ptp2_unref_gnum2);


            int         **ptp2_gnum1_come_from_idx;
            PDM_g_num_t **ptp2_gnum1_come_from;
            PDM_part_to_part_gnum1_come_from_get (_ptsp,
                                                  &ptp2_gnum1_come_from_idx,
                                                  &ptp2_gnum1_come_from);


            for (int i = 0; i < _nPart; i++) {
              double *referenceData  = (double *) referenceField->dataGet(i, CWP_FIELD_MAP_TARGET);
              for (int j = 0; j < ptp2_n_ref_gnum2[i]; j++) {
                assert ((ptp2_gnum1_come_from_idx[i][j+1] - ptp2_gnum1_come_from_idx[i][j]) == 1);
              }
              memcpy(referenceData, _recv_buffer[intId][i], dataTypeSize * nComponent * ptp2_n_ref_gnum2[i]);
            }
          }
        }

        if (_cpl->commTypeGet() == CWP_COMM_PAR_WITHOUT_PART) {
          /* Root rank broadcasts the referenceField inside the intraConnectableComm */
          MPI_Group globalGroup;
          MPI_Comm_group(_localCodeProperties->globalCommGet(), &globalGroup);
          int localRootRankInIntraConnectableComm;
          MPI_Group_translate_ranks(globalGroup, 1, &localRootRank,
                                    _localCodeProperties->connectableGroupGet(), &localRootRankInIntraConnectableComm);

          int rankInIntraConnectableComm;
          MPI_Comm_rank(_localCodeProperties->connectableCommGet(), &rankInIntraConnectableComm);

          // Easier way to get 'n_elt'?
          int nComponent = referenceField->nComponentGet();
          int  n_elt;
          int  n_part1 = 0;
          int  n_part2 = 0;
          int *n_elt1  = NULL;
          int *n_elt2  = NULL;
          PDM_part_to_part_n_part_and_n_elt_get(_ptsp,
                                                &n_part1,
                                                &n_part2,
                                                &n_elt1,
                                                &n_elt2);
          for (int i = 0; i < _nPart; i++) {
            if (globalRank == localRootRank) {
              n_elt = n_elt2[i];
            }
            MPI_Bcast((void *) &n_elt,
                      1,
                      MPI_INT,
                      localRootRankInIntraConnectableComm,
                      _localCodeProperties->connectableCommGet());

            double *referenceData = (double *) referenceField->dataGet(i, CWP_FIELD_MAP_TARGET);

            MPI_Bcast((void *) referenceData,
                      n_elt * nComponent,
                      MPI_DOUBLE,
                      localRootRankInIntraConnectableComm,
                      _localCodeProperties->connectableCommGet());

            if (referenceField->computedTgtBcastIsEnabled()) {
              // Broadcast (un)computed tgt as well
              MPI_Bcast((void *) &_n_computed_tgt[i],
                        1,
                        MPI_INT,
                        localRootRankInIntraConnectableComm,
                        _localCodeProperties->connectableCommGet());

              if (rankInIntraConnectableComm != localRootRankInIntraConnectableComm) {
                if (_computed_tgt[i] != NULL) {
                  free(_computed_tgt[i]);
                }
                _computed_tgt[i] = (int *) malloc(sizeof(int) * _n_computed_tgt[i]);
              }
              MPI_Bcast((void *) _computed_tgt[i],
                        _n_computed_tgt[i],
                        MPI_INT,
                        localRootRankInIntraConnectableComm,
                        _localCodeProperties->connectableCommGet());
            }
          }
        }


        if (cpl_spatial_interp->_send_buffer[cpl_intId] != NULL) {
          for (int i = 0; i < _cplNPart; i++) {
            if (cpl_spatial_interp->_send_buffer[cpl_intId][i] != NULL) {
              free (cpl_spatial_interp->_send_buffer[cpl_intId][i]);
              cpl_spatial_interp->_send_buffer[cpl_intId][i] = NULL;
            }
          }
          free (cpl_spatial_interp->_send_buffer[cpl_intId]);
          cpl_spatial_interp->_send_buffer[cpl_intId] = NULL;
        }
        if (cpl_spatial_interp->_recv_buffer[cpl_intId] != NULL) {
          for (int i = 0; i < _cplNPart; i++) {
            if (cpl_spatial_interp->_recv_buffer[cpl_intId][i] != NULL) {
              free (cpl_spatial_interp->_recv_buffer[cpl_intId][i]);
              cpl_spatial_interp->_recv_buffer[cpl_intId][i] = NULL;
            }
          }
          free (cpl_spatial_interp->_recv_buffer[cpl_intId]);
          cpl_spatial_interp->_recv_buffer[cpl_intId] = NULL;
        }

        if (_send_buffer[intId] != NULL) {
          for (int i = 0; i < _nPart; i++) {
            if (_send_buffer[intId][i] != NULL) {
              free (_send_buffer[intId][i]);
              _send_buffer[intId][i] = NULL;
            }
          }
          free (_send_buffer[intId]);
          _send_buffer[intId] = NULL;
        }
        if (_recv_buffer[intId] != NULL) {
          for (int i = 0; i < _nPart; i++) {
            if (_recv_buffer[intId][i] != NULL) {
              free (_recv_buffer[intId][i]);
              _recv_buffer[intId][i] = NULL;
            }
          }
          free (_recv_buffer[intId]);
          _recv_buffer[intId] = NULL;
        }

        PDM_writer_t* writer =  _cpl->writerGet();

        if (writer != nullptr) {
          if ((_cpl->NStepGet() % _cpl->freqWriterGet()) == 0) {
            cpl_referenceField->write(CWP_FIELD_EXCH_SEND);
          }
        }

        if (writer != nullptr) {
          if ((_cpl->NStepGet() % _cpl->freqWriterGet()) == 0) {
            referenceField->write(CWP_FIELD_EXCH_RECV);
          }
        }

      } // if local rank has to work
    } // end if joint

    // set field to received
    referenceField->is_recv_yet_set(1);
  }



/***************************************************************************/
/***************************************************************************/



} // end namespace cwipi

/**
 * \endcond
 */
