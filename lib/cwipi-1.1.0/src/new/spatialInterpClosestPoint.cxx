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

#include "pdm_closest_points.h"
#include "spatialInterpClosestPoint.hxx"
#include "cwp_priv.h"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_linear_algebra.h"
#include "pdm_vtk.h"

CWP_CLANG_SUPPRESS_WARNING("-Wunused-private-field")


namespace cwipi {
    //CWIPI_CLANG_SUPPRESS_WARNING("-Wunused-private-field")
    SpatialInterpClosestPoint::SpatialInterpClosestPoint() = default;

    SpatialInterpClosestPoint::~SpatialInterpClosestPoint
    (
     )
    {
      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_closest_src_gnum[i_part] != NULL) {
          free(_closest_src_gnum[i_part]);
        }
        // if (_weights[i_part] != NULL) {
        //   free(_weights[i_part]);
        // }
      }

      free ( _closest_src_gnum);
      // free ( _weights);

      for (int i_part = 0; i_part < _nPart; i_part++) {
        if (_tgt_in_src_idx[i_part] != NULL) {
          free(_tgt_in_src_idx[i_part]);
        }
        if (_tgt_in_src_gnum[i_part] != NULL) {
          free(_tgt_in_src_gnum[i_part]);
        }
        if (_tgt_in_src_dist[i_part] != NULL) {
          free(_tgt_in_src_dist[i_part]);
        }
      }

      free ( _tgt_in_src_idx);
      free ( _tgt_in_src_gnum);
      free ( _tgt_in_src_dist);

      if (_send_coord != NULL) {
        free(_send_coord);
        _send_coord = NULL;
      }

      if (_recv_coord != NULL) {
        int n_part_tgt = 0;
        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          n_part_tgt = _nPart;
        }
        else {
          n_part_tgt  = _cplNPart;
        }

        for (int i_part = 0; i_part < n_part_tgt; i_part++) {
          if (_recv_coord[i_part] != NULL) {
            free(_recv_coord[i_part]);
          }
        }
        free(_recv_coord);
        _recv_coord = NULL;
      }
    }


    /**
     *
     * \brief SpatialInterpClosestPoint Init.
     *
     */

    void
    SpatialInterpClosestPoint::init
    (
     Coupling                   *coupling,
     CWP_Dof_location_t          localCodeDofLocation,
     CWP_Dof_location_t          cplCodeDofLocation,
     SpatialInterpExchDirection  exchDirection
     )
    {
      SpatialInterp::init (coupling,
                           localCodeDofLocation,
                           cplCodeDofLocation,
                           exchDirection);

      _coordinates_exchanged = COORD_NOT_YET_EXCHANGED;

      _interpolation_time = CWP_SPATIAL_INTERP_AT_RECV;

      //
      // Data for PDM_part_to_part_t
      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
        for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            _src_gnum  [i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
            _src_n_gnum[i_part] = _mesh->getPartNElts (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            _src_gnum  [i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
            _src_n_gnum[i_part] = _mesh->getPartNVertex (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            _src_gnum  [i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
            _src_n_gnum[i_part] = _cpl->userTargetNGet (i_part);
          }
        }
      }
      else {
        for (int i_part = 0 ; i_part < _nPart ; i_part++) {
          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
            _tgt_n_gnum[i_part] = _mesh->getPartNElts (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
            _tgt_n_gnum[i_part] = _mesh->getPartNVertex (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            _tgt_gnum  [i_part] = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
            _tgt_n_gnum[i_part] = _cpl->userTargetNGet (i_part);
          }
        }
      }

      //
      // Target properties
      _closest_src_gnum = (PDM_g_num_t**) malloc (sizeof(PDM_g_num_t*) * _nPart);

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _closest_src_gnum[i_part] = NULL;
      }

      //
      // Source properties
      _tgt_in_src_idx  =  (int **) malloc (sizeof(int *) * (_nPart));
      _tgt_in_src_gnum = (PDM_g_num_t**) malloc (sizeof(PDM_g_num_t*) * _nPart);
      
      _tgt_in_src_dist =  (double **) malloc (sizeof(double *) * (_nPart));

      for (int i_part = 0; i_part < _nPart; i_part++) {
        _tgt_in_src_idx [i_part] = NULL;
        _tgt_in_src_gnum[i_part] = NULL;
        _tgt_in_src_dist[i_part] = NULL;
      }

      _send_coord = NULL;
      _recv_coord = NULL;
    }




    void SpatialInterpClosestPoint::weightsCompute() {

      /* Set source and target point clouds */
      set_PDM_object();


      /* Compute */
      if (!_coupledCodeProperties->localCodeIs() ||
          (_coupledCodeProperties->localCodeIs() && _localCodeProperties->idGet() < _coupledCodeProperties->idGet())) {
        PDM_closest_points_compute(_id_pdm);

        int dump_times = 0;
        char *env_var = NULL;
        env_var = getenv ("CWP_DUMP_TIMES");
        if (env_var != NULL) {
          dump_times = atoi(env_var);
        }
        if (dump_times) {
          PDM_closest_points_dump_times(_id_pdm);
        }
      }

      /* Get PDM part_to_part object */
      if (_id_pdm != NULL) {
        if (!_reverse) {
          PDM_closest_points_part_to_part_get(_id_pdm,
                                              &_ptsp,
                                              PDM_OWNERSHIP_USER);
        }
      }

      int n_neighbors = CWP_CLOSEST_POINTS_N_NEIGHBORS;
      std::map<std::string, int> prop = _cpl->SpatialInterpPropertiesIntGet();
      std::map<std::string, int>::iterator it;

      it = prop.find("n_neighbors");
      if (it != prop.end()) {
        n_neighbors = it->second;
      }

      /* Get PDM results */
      if (!_coupledCodeProperties->localCodeIs()) {
        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

            if (_reverse) {
              _tgt_in_src_idx[i_part] = PDM_array_new_idx_from_const_stride_int(n_neighbors,
                                                                                _src_n_gnum[i_part]);
              PDM_closest_points_get(_id_pdm,
                                     i_part,
                                     &(_tgt_in_src_gnum[i_part]),
                                     &(_tgt_in_src_dist[i_part]));
            }
            else {
              PDM_closest_points_tgt_in_src_get(_id_pdm,
                                                i_part,
                                                &(_tgt_in_src_idx [i_part]),
                                                &(_tgt_in_src_gnum[i_part]));
              PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                     i_part,
                                                     &(_tgt_in_src_idx [i_part]),
                                                     &(_tgt_in_src_dist[i_part]));
            }

          }
          else {
            _tgt_in_src_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
            _tgt_in_src_idx[i_part][0] = 0;

            if (_reverse) {
              int *tgt_to_src_idx = NULL;
              PDM_closest_points_tgt_in_src_get(_id_pdm,
                                                i_part,
                                                &tgt_to_src_idx,
                                                &(_closest_src_gnum[i_part]));
              free(tgt_to_src_idx);
              PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                     i_part,
                                                     &tgt_to_src_idx,
                                                     &(_weights[i_part]));
            }
            else {
              PDM_closest_points_get(_id_pdm,
                                     i_part,
                                     &(_closest_src_gnum[i_part]),
                                     &(_weights[i_part]));

            }

          }

        }
      }
      else {
        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          SpatialInterpClosestPoint *cpl_spatial_interp;

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          else {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }


          if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
            for (int i_part = 0; i_part < _nPart; i_part++) {
              if (_reverse) {
                _tgt_in_src_idx[i_part] = PDM_array_new_idx_from_const_stride_int(n_neighbors,
                                                                                  _src_n_gnum[i_part]);
                PDM_closest_points_get(_id_pdm,
                                       i_part,
                                       &(_tgt_in_src_gnum[i_part]),
                                       &(_tgt_in_src_dist[i_part]));
              }
              else {
                PDM_closest_points_tgt_in_src_get(_id_pdm,
                                                  i_part,
                                                  &(_tgt_in_src_idx [i_part]),
                                                  &(_tgt_in_src_gnum[i_part]));
                PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                       i_part,
                                                       &(_tgt_in_src_idx [i_part]),
                                                       &(_tgt_in_src_dist[i_part]));
              }
            }

            for (int i_part = 0; i_part < _cplNPart; i_part++) {
              cpl_spatial_interp->_tgt_in_src_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
              cpl_spatial_interp->_tgt_in_src_idx[i_part][0] = 0;

              if (_reverse) {
                int *tgt_to_src_idx = NULL;
                PDM_closest_points_tgt_in_src_get(_id_pdm,
                                                  i_part,
                                                  &tgt_to_src_idx,
                                                  &(cpl_spatial_interp->_closest_src_gnum[i_part]));
                free(tgt_to_src_idx);
                PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                       i_part,
                                                       &tgt_to_src_idx,
                                                       &(cpl_spatial_interp->_weights[i_part]));
              }
              else {
                PDM_closest_points_get(_id_pdm,
                                       i_part,
                                       &(cpl_spatial_interp->_closest_src_gnum[i_part]),
                                       &(cpl_spatial_interp->_weights[i_part]));
              }
            }
          }
          else {
            for (int i_part = 0; i_part < _nPart; i_part++) {
              _tgt_in_src_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
              _tgt_in_src_idx[i_part][0] = 0;

              if (_reverse) {
                int *tgt_to_src_idx = NULL;
                PDM_closest_points_tgt_in_src_get(_id_pdm,
                                                  i_part,
                                                  &tgt_to_src_idx,
                                                  &(_closest_src_gnum[i_part]));
                free(tgt_to_src_idx);
                PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                       i_part,
                                                       &tgt_to_src_idx,
                                                       &(_weights[i_part]));
              }
              else {
                PDM_closest_points_get(_id_pdm,
                                       i_part,
                                       &(_closest_src_gnum[i_part]),
                                       &(_weights[i_part]));
              }
            }

            for (int i_part = 0; i_part < _cplNPart; i_part++) {
              if (_reverse) {
                cpl_spatial_interp->_tgt_in_src_idx[i_part] = PDM_array_new_idx_from_const_stride_int(n_neighbors,
                                                                                                      cpl_spatial_interp->_src_n_gnum[i_part]);

                PDM_closest_points_get(_id_pdm,
                                       i_part,
                                       &(cpl_spatial_interp->_tgt_in_src_gnum[i_part]),
                                       &(cpl_spatial_interp->_tgt_in_src_dist[i_part]));
              }
              else {
                PDM_closest_points_tgt_in_src_get(_id_pdm,
                                                  i_part,
                                                  &(cpl_spatial_interp->_tgt_in_src_idx [i_part]),
                                                  &(cpl_spatial_interp->_tgt_in_src_gnum[i_part]));
                PDM_closest_points_tgt_in_src_dist_get(_id_pdm,
                                                       i_part,
                                                       &(cpl_spatial_interp->_tgt_in_src_idx [i_part]),
                                                       &(cpl_spatial_interp->_tgt_in_src_dist[i_part]));
              }

            }
          }
        }
      }

      /* Free PDM object */
      if (!_coupledCodeProperties->localCodeIs()) {
        PDM_closest_points_free(_id_pdm);
        _id_pdm = nullptr;
      }
      else {
        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          PDM_closest_points_free(_id_pdm);
          _id_pdm = nullptr;

          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          SpatialInterpClosestPoint *cpl_spatial_interp;

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          else {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          cpl_spatial_interp->_id_pdm = nullptr;
        }
      }

      /* Create part_to_part object if null */

      if (!_coupledCodeProperties->localCodeIs()) {
        if (_ptsp == NULL) {
          _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **) _src_gnum,
                                           (const int          *) _src_n_gnum,
                                           _nPart,
                                           (const PDM_g_num_t **) _tgt_gnum,
                                           (const int          *) _tgt_n_gnum,
                                           _nPart,
                                           (const int         **) _tgt_in_src_idx,
                                           (const PDM_g_num_t **) _tgt_in_src_gnum,
                                           _pdmUnionComm);
        }
      }
      else {
        if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
          SpatialInterpClosestPoint *cpl_spatial_interp;

          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }
          else {
            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
          }

          if (_ptsp == NULL) {
            if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
              _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **) _src_gnum,
                                               (const int          *) _src_n_gnum,
                                               _nPart,
                                               (const PDM_g_num_t **) cpl_spatial_interp->_tgt_gnum,
                                               (const int          *) cpl_spatial_interp->_tgt_n_gnum,
                                               _cplNPart,
                                               (const int         **) _tgt_in_src_idx,
                                               (const PDM_g_num_t **) _tgt_in_src_gnum,
                                               _pdmUnionComm);
            }
            else {
              _ptsp = PDM_part_to_part_create ((const PDM_g_num_t **) cpl_spatial_interp->_src_gnum,
                                               (const int          *) cpl_spatial_interp->_src_n_gnum,
                                               _cplNPart,
                                               (const PDM_g_num_t **) _tgt_gnum,
                                               (const int          *) _tgt_n_gnum,
                                               _nPart,
                                               (const int         **) cpl_spatial_interp->_tgt_in_src_idx,
                                               (const PDM_g_num_t **) cpl_spatial_interp->_tgt_in_src_gnum,
                                               _pdmUnionComm);
            }
          }

          cpl_spatial_interp->_ptsp = _ptsp;

        }
      }


      if (_ptsp != NULL) {
        int  *n_ref_tgt = NULL;
        int **ref_tgt   = NULL;
        PDM_part_to_part_ref_lnum2_get(_ptsp, &n_ref_tgt, &ref_tgt);

        int          *n_src          = NULL;
        int         **src_to_tgt_idx = NULL;
        PDM_g_num_t **src_to_tgt     = NULL;
        PDM_part_to_part_part1_to_part2_get(_ptsp,
                                            &n_src,
                                            &src_to_tgt_idx,
                                            &src_to_tgt);

        if (!_coupledCodeProperties->localCodeIs()) {

          if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
            for (int i_part = 0; i_part < _nPart; i_part++) {
              _n_computed_tgt[i_part] = n_ref_tgt[i_part];
              _computed_tgt[i_part] = (int *) malloc(sizeof(int) * _n_computed_tgt[i_part]);
              memcpy(_computed_tgt[i_part], ref_tgt[i_part], sizeof(int) * _n_computed_tgt[i_part]);
            }
          }
          else {
            for (int i_part = 0; i_part < _nPart; i_part++) {
              _n_involved_sources_tgt[i_part] = 0;
              _involved_sources_tgt  [i_part] = (int *) malloc(sizeof(int) * n_src[i_part]);
              for (int i = 0; i < n_src[i_part]; i++) {
                if (src_to_tgt_idx[i_part][i+1] > src_to_tgt_idx[i_part][i]) {
                  _involved_sources_tgt[i_part][_n_involved_sources_tgt[i_part]++] = i+1;
                }
              }
              _involved_sources_tgt[i_part] = (int *) realloc(_involved_sources_tgt[i_part],
                                                              sizeof(int) * _n_involved_sources_tgt[i_part]);
            }
          }

        }
        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            SpatialInterpClosestPoint *cpl_spatial_interp;

            if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
              std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
              cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

              for (int i_part = 0; i_part < _cplNPart; i_part++) {
                cpl_spatial_interp->_n_computed_tgt[i_part] = n_ref_tgt[i_part];
                cpl_spatial_interp->_computed_tgt[i_part] = (int *) malloc(sizeof(int) * cpl_spatial_interp->_n_computed_tgt[i_part]);
                memcpy(cpl_spatial_interp->_computed_tgt[i_part], ref_tgt[i_part], sizeof(int) * cpl_spatial_interp->_n_computed_tgt[i_part]);
              }

              for (int i_part = 0; i_part < _nPart; i_part++) {
                _n_involved_sources_tgt[i_part] = 0;
                _involved_sources_tgt  [i_part] = (int *) malloc(sizeof(int) * n_src[i_part]);
                for (int i = 0; i < n_src[i_part]; i++) {
                  if (src_to_tgt_idx[i_part][i+1] > src_to_tgt_idx[i_part][i]) {
                    _involved_sources_tgt[i_part][_n_involved_sources_tgt[i_part]++] = i+1;
                  }
                }
                _involved_sources_tgt[i_part] = (int *) realloc(_involved_sources_tgt[i_part],
                                                                sizeof(int) * _n_involved_sources_tgt[i_part]);
              }
            }
            else {
              std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
              cpl_spatial_interp = dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

              for (int i_part = 0; i_part < _nPart; i_part++) {
                _n_computed_tgt[i_part] = n_ref_tgt[i_part];
                _computed_tgt[i_part] = (int *) malloc(sizeof(int) * _n_computed_tgt[i_part]);
                memcpy(_computed_tgt[i_part], ref_tgt[i_part], sizeof(int) * _n_computed_tgt[i_part]);
              }

              for (int i_part = 0; i_part < _cplNPart; i_part++) {
                cpl_spatial_interp->_n_involved_sources_tgt[i_part] = 0;
                cpl_spatial_interp->_involved_sources_tgt  [i_part] = (int *) malloc(sizeof(int) * n_src[i_part]);
                for (int i = 0; i < n_src[i_part]; i++) {
                  if (src_to_tgt_idx[i_part][i+1] > src_to_tgt_idx[i_part][i]) {
                    cpl_spatial_interp->_involved_sources_tgt[i_part][cpl_spatial_interp->_n_involved_sources_tgt[i_part]++] = i+1;
                  }
                }
                cpl_spatial_interp->_involved_sources_tgt[i_part] = (int *) realloc(cpl_spatial_interp->_involved_sources_tgt[i_part],
                                                                                    sizeof(int) * cpl_spatial_interp->_n_involved_sources_tgt[i_part]);
              }
            }
          }

        }
      }

    }


    void SpatialInterpClosestPoint::issend(Field *referenceField) {

      if (_coordinates_exchanged == COORD_NOT_YET_EXCHANGED) {
        /* Send source points coordinates */
        if (!_coupledCodeProperties->localCodeIs()) {

          PDM_part_to_part_iexch(_ptsp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                                 1,
                                 3*sizeof(double),
                                 NULL,
                (const void  **) _send_coord,
                                 NULL,
                (      void ***) &_recv_coord,
                                 &(_send_coord_request));
          _recv_coord_request = _send_coord_request;
        }

        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            SpatialInterpClosestPoint *cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            PDM_part_to_part_iexch(_ptsp,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                                   1,
                                   3*sizeof(double),
                                   NULL,
                  (const void  **) _send_coord,
                                   NULL,
                  (      void ***) &cpl_spatial_interp->_recv_coord,
                                   &(_send_coord_request));
            cpl_spatial_interp->_recv_coord_request = _send_coord_request;
          }
        }

        _coordinates_exchanged = COORD_EXCHANGE_INITIALIZED;
      }

      SpatialInterp::issend(referenceField);
    }


    void SpatialInterpClosestPoint::waitIssend(Field *referenceField) {

      assert(_coordinates_exchanged != COORD_NOT_YET_EXCHANGED);

      if (_coordinates_exchanged == COORD_EXCHANGE_INITIALIZED) {

        if (!_coupledCodeProperties->localCodeIs()) {
          PDM_part_to_part_iexch_wait(_ptsp, _send_coord_request);
        }

        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
            SpatialInterpClosestPoint *cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            PDM_part_to_part_iexch_wait(_ptsp, cpl_spatial_interp->_recv_coord_request);
          }
        }

        _coordinates_exchanged = COORD_EXCHANGE_FINALIZED;
      }

      SpatialInterp::waitIssend(referenceField);
    }


    void SpatialInterpClosestPoint::irecv(Field *referenceField) {

      /* Receive source points coordinates */
      if (_coordinates_exchanged == COORD_NOT_YET_EXCHANGED) {

        if (!_coupledCodeProperties->localCodeIs()) {

          PDM_part_to_part_iexch(_ptsp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                                 1,
                                 3*sizeof(double),
                                 NULL,
                (const void  **) _send_coord,
                                 NULL,
                (      void ***) &_recv_coord,
                                 &(_send_coord_request));
          _recv_coord_request = _send_coord_request;
        }

        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            SpatialInterpClosestPoint *cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            PDM_part_to_part_iexch(_ptsp,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                                   1,
                                   3*sizeof(double),
                                   NULL,
                  (const void  **) cpl_spatial_interp->_send_coord,
                                   NULL,
                  (      void ***) &_recv_coord,
                                   &(cpl_spatial_interp->_send_coord_request));
            _recv_coord_request = cpl_spatial_interp->_send_coord_request;
          }
        }

        _coordinates_exchanged = COORD_EXCHANGE_INITIALIZED;
      }

      SpatialInterp::irecv(referenceField);
    }


    void SpatialInterpClosestPoint::waitIrecv(Field *referenceField) {

      assert(_coordinates_exchanged != COORD_NOT_YET_EXCHANGED);

      if (_coordinates_exchanged == COORD_EXCHANGE_INITIALIZED) {

        if (!_coupledCodeProperties->localCodeIs()) {
          PDM_part_to_part_iexch_wait(_ptsp, _send_coord_request);
        }

        else {
          if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

            std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
            SpatialInterpClosestPoint *cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

            PDM_part_to_part_iexch_wait(_ptsp, cpl_spatial_interp->_send_coord_request);
          }
        }

        _coordinates_exchanged = COORD_EXCHANGE_FINALIZED;
      }

      SpatialInterp::waitIrecv(referenceField);
    }





    static void _interp_idw
    (
      const int     n_neighbors,
      const int     stride,
      const double *src_value,
      const double *src_dist2,
            double *tgt_value
     )
    {
      const double eps_dist2 = 1e-30;
      for (int j = 0; j < stride; j++) {
        tgt_value[j] = 0.;
      }

      double sum_w = 0.;
      for (int i = 0; i < n_neighbors; i++) {
        double w = 1./std::max(eps_dist2, src_dist2[i]);
        sum_w += w;
        for (int j = 0; j < stride; j++) {
          tgt_value[j] += w*src_value[stride*i + j];
        }
      }

      sum_w = 1./sum_w;
      for (int j = 0; j < stride; j++) {
        tgt_value[j] *= sum_w;
      }
    }


    static void _interp_least_squares
    (
      const int     n_neighbors,
      const int     stride,
      const double *src_value,
      const double *src_coord,
      const double *src_dist2,
      const double *tgt_coord,
            double *tgt_value
     )
    {
      double A[4*4] = {0.};

      double *b = (double *) malloc(sizeof(double) * 4 *stride);
      double *coeff = (double *) malloc(sizeof(double) * 4 *stride);

      for (int i = 0; i < n_neighbors; i++) {
        b[i] = 0.;
        coeff[i] = 0.;
      }

      for (int i = 0; i < n_neighbors; i++) {
        // log_trace("i = %d / %d\n", i, n_neighbors);
        double x = src_coord[3*i    ];
        double y = src_coord[3*i + 1];
        double z = src_coord[3*i + 2];

        A[4*0+0] += x * x;
        A[4*0+1] += x * y;
        A[4*0+2] += x * z;
        A[4*0+3] += x;

        A[4*1+1] += y * y;
        A[4*1+2] += y * z;
        A[4*1+3] += y;

        A[4*2+2] += z * z;
        A[4*2+3] += z;

        A[4*3+3] += 1.;

        for (int j = 0; j < stride; j++) {
          double f = src_value[stride*i + j];
          b[stride*0 + j] += x * f;
          b[stride*1 + j] += y * f;
          b[stride*2 + j] += z * f;
          b[stride*3 + j] += f;
        }
      }

      /* Symmetrize */
      A[4*1+0] = A[4*0+1];
      A[4*2+0] = A[4*0+2];
      A[4*3+0] = A[4*0+3];

      A[4*2+1] = A[4*1+2];
      A[4*3+1] = A[4*1+3];

      A[4*3+2] = A[4*2+3];

      int stat = PDM_linear_algebra_linsolve_svd(4,
                                                 4,
                                                 stride,
                                                 0.,
                                      (double *) A,
                                      (double *) b,
                                                 coeff);
      if (stat == 0) {
        for (int j = 0; j < stride; j++) {
          tgt_value[j] =
          coeff[stride*0 + j] * tgt_coord[0] +
          coeff[stride*1 + j] * tgt_coord[1] +
          coeff[stride*2 + j] * tgt_coord[2] +
          coeff[stride*3 + j];
        }
      }
      else {
        _interp_idw(n_neighbors,
                    stride,
                    src_value,
                    src_dist2,
                    tgt_value);
        // for (int j = 0; j < stride; j++) {
        //   tgt_value[j] = (j+1)*1000;
        // }
      }

      free(b);
      free(coeff);
    }


    static inline void _basis_vector
    (
     const int     dim,
     const int     degree,
     const double *x,
     const double *x0,
           double *b
     )
    {
      b[0] = 1.;

      for (int j = 0; j < degree; j++) {
        for (int i = 0; i < dim; i++) {
          double _x = x[i];
          if (x0 != NULL) {
            _x -= x0[i];
          }
          b[1+dim*j+i] = 1.;
          for (int k = 0; k <= j; k++) {
            b[1+dim*j+i] *= _x;
          }
        }
      }
    }


    static inline double _weight_function
    (
     const double d2
     )
    {
      const double eps2 = 1e-12;

      return 1. / (d2 + eps2);
    }

    static void _interp_weighted_least_squares
    (
      const int     degree,
      const int     n_neighbors,
      const int     stride,
      const double *src_value,
      const double *src_coord,
      const double *src_dist2,
      const double *tgt_coord,
            double *tgt_value
     )
    {
      int stat = 1;
      int siz = (1 + degree*3);

      double *c = (double *) malloc(sizeof(double) * siz *stride);

      if (n_neighbors >= siz) {

        /* Special case : coincident points */
        const double eps_dist2 = 1e-24;

        for (int i = 0; i < n_neighbors; i++) {
          if (src_dist2[i] < eps_dist2) {
            for (int j = 0; j < stride; j++) {
              tgt_value[j] = src_value[stride*i+j];
            }
            free (c);
            return;
          }
        }


        double *A = (double *) malloc(sizeof(double) * siz * siz);
        double *rhs = (double *) malloc(sizeof(double) * siz *stride);
        double *b = (double *) malloc(sizeof(double) * siz );

        for (int i = 0; i < siz * siz; i++) {
          A[i] = 0;
        }

        for (int i = 0; i < siz * stride; i++) {
          rhs[i] = 0;
        }

        for (int i = 0; i < siz; i++) {
          b[i] = 0;
        }

        for (int i = 0; i < n_neighbors; i++) {
          double wi = _weight_function(src_dist2[i]);

          _basis_vector(3,
                        degree,
                        src_coord + 3*i,
                        tgt_coord,
                        b);

          for (int j = 0; j < siz; j++) {
            for (int k = 0; k < stride; k++) {
              rhs[stride*j+k] += wi * b[j] * src_value[stride*i+k];
            }

            for (int k = j; k < siz; k++) {
              A[siz*j+k] += wi * b[j] * b[k];
            }
            for (int k = 0; k < j; k++) {
              A[siz*j+k] = A[siz*k+j];
            }
          }
        }

        stat = PDM_linear_algebra_linsolve_svd(siz,
                                               siz,
                                               stride,
                                               0.,
                                    (double *) A,
                                    (double *) rhs,
                                               c);

        free (A);
        free (rhs);
        free (b);
      }

      if (stat == 0) {
        for (int j = 0; j < stride; j++) {
          tgt_value[j] = c[j];
        }
      }
      else {
        // log_trace("singular matrix! ");
        // PDM_log_trace_array_double(tgt_coord, 3, "tgt_coord : ");
        _interp_idw(n_neighbors,
                    stride,
                    src_value,
                    src_dist2,
                    tgt_value);
      }

      free (c);

    }




    void SpatialInterpClosestPoint::interpolate(Field *referenceField, double **buffer) {

      int nComponent = referenceField->nComponentGet();
      // CWP_Dof_location_t referenceFieldType = referenceField->locationGet();
      CWP_Interpolation_t interpolationType = referenceField->interpolationTypeGet();
      const CWP_Field_storage_t storage = referenceField->storageTypeGet();

      if (interpolationType == CWP_INTERPOLATION_USER) {
        CWP_Interp_function_t   interpolationFunction  = referenceField->interpolationFunctionGet();
        CWP_Interp_function_t   interpolationFunctionF = referenceField->interpFunctionFGet();
        CWP_Interp_function_p_t interpolationFunctionP = referenceField->interpFunctionPGet();

        CWP_Interp_function_t _interpolationFunction = NULL;
        if (interpolationFunction != NULL) {
          _interpolationFunction = interpolationFunction;
        }
        else if (interpolationFunctionF != NULL) {
          _interpolationFunction = interpolationFunctionF;
        }

        if (_interpolationFunction == NULL &&
            interpolationFunctionP == NULL) {
          PDM_error(__FILE__, __LINE__, 0, "Undefined user interpolation function\n");
        }

        for (int i_part = 0 ; i_part < _nPart ; i_part++) {

          if (0) {
            int  *n_ref = NULL;
            int **ref   = NULL;
            PDM_part_to_part_ref_lnum2_get(_ptsp,
                                           &n_ref,
                                           &ref);

            int  *n_unref = NULL;
            int **unref   = NULL;
            PDM_part_to_part_unref_lnum2_get(_ptsp,
                                             &n_unref,
                                             &unref);

            int         **come_from_idx = NULL;
            PDM_g_num_t **come_from     = NULL;
            PDM_part_to_part_gnum1_come_from_get(_ptsp,
                                                 &come_from_idx,
                                                 &come_from);

            int n_tgt = n_ref[i_part] + n_unref[i_part];

            double *src_coord = _recv_coord[i_part];
            const double *tgt_coord = NULL;
            if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              tgt_coord = _mesh->eltCentersGet(i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              tgt_coord = _mesh->getVertexCoords(i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
              tgt_coord = _cpl->userTargetCoordsGet(i_part);
            }

            int i_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
            char filename[999];
            sprintf(filename, "mapping_%d.vtk", i_rank);

            int n_pts = n_tgt + come_from_idx[i_part][n_ref[i_part]];
            double *coord = (double *) malloc(sizeof(double) * n_pts * 3);
            PDM_g_num_t *gnum = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_pts);

            memcpy(coord, tgt_coord, sizeof(double) * n_tgt * 3);
            memcpy(coord + n_tgt*3, src_coord, sizeof(double) * come_from_idx[i_part][n_ref[i_part]] * 3);

            memcpy(gnum, _tgt_gnum[i_part], sizeof(PDM_g_num_t) * n_tgt);
            for (int i = 0; i < come_from_idx[i_part][n_ref[i_part]]; i++) {
              gnum[n_tgt+i] = -_closest_src_gnum[i_part][i];
            }

            int *connec = (int *) malloc(sizeof(int) * come_from_idx[i_part][n_ref[i_part]] * 2);
            for (int i = 0; i < n_ref[i_part]; i++) {
              for (int j = come_from_idx[i_part][i]; j < come_from_idx[i_part][i+1]; j++) {
                connec[2*j  ] = ref[i_part][i];
                connec[2*j+1] = n_tgt + j + 1;
              }
            }

            PDM_vtk_write_std_elements(filename,
                                       n_pts,
                                       coord,
                                       gnum,
                                       PDM_MESH_NODAL_BAR2,
                                       come_from_idx[i_part][n_ref[i_part]],
                                       connec,
                                       NULL,
                                       0,
                                       NULL,
                                       NULL);
            free(coord);
            free(gnum);
            free(connec);
          }

          if (interpolationFunctionP != NULL) {
            (*interpolationFunctionP)(referenceField->pythonObjectGet(),
                                      i_part,
                                      (double *) buffer[i_part],
                                      (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_TARGET));
          }
          else {
            (*_interpolationFunction)(_localCodeProperties->nameGet().c_str(),
                                      _cpl->IdGet().c_str(),
                                      referenceField->fieldIDGet().c_str(),
                                      i_part,
                                      (double *) buffer[i_part],
                                      (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_TARGET));
          }

        }

      }

      else {

        std::map<std::string, int> prop = _cpl->SpatialInterpPropertiesIntGet();
        std::map<std::string, int>::iterator it;

        int  *n_ref = NULL;
        int **ref   = NULL;
        PDM_part_to_part_ref_lnum2_get(_ptsp,
                                       &n_ref,
                                       &ref);

        int  *n_unref = NULL;
        int **unref   = NULL;
        PDM_part_to_part_unref_lnum2_get(_ptsp,
                                         &n_unref,
                                         &unref);

        int         **come_from_idx = NULL;
        PDM_g_num_t **come_from     = NULL;
        PDM_part_to_part_gnum1_come_from_get(_ptsp,
                                             &come_from_idx,
                                             &come_from);

        int max_n_neighbors = 0;
        for (int i_part = 0; i_part < _nPart; i_part++) {
          for (int i = 0; i < n_ref[i_part]; i++) {
            max_n_neighbors = std::max(max_n_neighbors,
                                         come_from_idx[i_part][i+1] - come_from_idx[i_part][i]);
          }
        }

        int polyfit_degree = CWP_CLOSEST_POINTS_POLYFIT_DEGREE;
        it = prop.find("polyfit_degree");
        if (it != prop.end()) {
          polyfit_degree = it->second;
        }

        int use_idw_interpolation = 0;
        double       *src_coord = NULL;
        const double *tgt_coord = NULL;

        double *src_value = NULL;
        double *tgt_value = NULL;

        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          src_value = (double *) malloc(sizeof(double) * nComponent * max_n_neighbors);
          tgt_value = (double *) malloc(sizeof(double) * nComponent);
        }

        for (int i_part = 0; i_part < _nPart; i_part++) {

          double *referenceData = (double *) referenceField->dataGet(i_part, CWP_FIELD_MAP_TARGET);

          int n_tgt = n_ref[i_part] + n_unref[i_part];

          /* Set values for uncomputed targets */
          for (int i = 0; i < n_unref[i_part]; i++) {
            int id = unref[i_part][i] - 1;
            if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
              for (int j = 0; j < nComponent; j++) {
                referenceData[n_tgt*j + id] = 0.;
              }
            }
            else {
              for (int j = 0; j < nComponent; j++) {
                referenceData[nComponent*id + j] = 0.;
              }
            }
          }

          if (!use_idw_interpolation) {
            if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              tgt_coord = _mesh->eltCentersGet(i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              tgt_coord = _mesh->getVertexCoords(i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
              tgt_coord = _cpl->userTargetCoordsGet(i_part);
            }
          }

          src_coord = _recv_coord[i_part];




          for (int i = 0; i < n_ref[i_part]; i++) {
            int id = ref[i_part][i] - 1;
            // log_trace("tgt pt %d / %d\n", i, n_pts);

            int n_neighbors = come_from_idx[i_part][i+1] - come_from_idx[i_part][i];

            if (0) {
              for (int k = come_from_idx[i_part][i]; k < come_from_idx[i_part][i+1]; k++) {
                log_trace("src point " PDM_FMT_G_NUM " : %f %f %f\n",
                          _closest_src_gnum[i_part][k],
                          src_coord[3*k  ],
                          src_coord[3*k+1],
                          src_coord[3*k+2]);
              }
            }

            if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
              // interlace src_value
              for (int k = 0; k < n_neighbors; k++) {
                for (int j = 0; j < nComponent; j++) {
                  src_value[nComponent*k + j] = buffer[i_part][j*come_from_idx[i_part][n_ref[i_part]] + come_from_idx[i_part][i] + k];
                }
              }
            }
            else {
              src_value = buffer[i_part] + nComponent*come_from_idx[i_part][i];
              tgt_value = referenceData + nComponent*id;
            }

            if (use_idw_interpolation) {
              _interp_idw(n_neighbors,
                          nComponent,
                          src_value,
                          _weights[i_part] + come_from_idx[i_part][i],
                          tgt_value);
            }
            else {
              if (1) {
                _interp_weighted_least_squares(polyfit_degree,
                                               n_neighbors,
                                               nComponent,
                                               src_value,
                                               src_coord + 3*come_from_idx[i_part][i],
                                               _weights[i_part] + come_from_idx[i_part][i],
                                               tgt_coord + 3*id,
                                               tgt_value);
              }
              else {
                _interp_least_squares(n_neighbors,
                                      nComponent,
                                      src_value,
                                      src_coord + 3*come_from_idx[i_part][i],
                                      _weights[i_part] + come_from_idx[i_part][i],
                                      tgt_coord + 3*id,
                                      tgt_value);
              }


              if (0) {
                log_trace("recv :\n");
                for (int k = 0; k < n_neighbors; k++) {
                  log_trace("  from " PDM_FMT_G_NUM " : ", _closest_src_gnum[i_part][come_from_idx[i_part][i] + k]);
                  PDM_log_trace_array_double(src_value + nComponent*k,
                                             nComponent,
                                             "");
                }
                PDM_log_trace_array_double(tgt_value,
                                           nComponent,
                                           "interpolated : ");
                PDM_log_trace_array_double(tgt_coord + 3*id,
                                           3,
                                           "tgt_coord    : ");
              }
            }

            if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
              // de-interlace tgt_value
              for (int j = 0; j < nComponent; j++) {
                referenceData[n_tgt*j + id] = tgt_value[j];
              }
            }

          }
        }

        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          free(src_value);
          free(tgt_value);
        }


      }


    }


    double **
    SpatialInterpClosestPoint::closest_src_coord_get()
    {
      return _recv_coord;
    }





    void SpatialInterpClosestPoint::clear() {

      SpatialInterp::clear();

      _coordinates_exchanged = COORD_NOT_YET_EXCHANGED;


      int cond1 = !_coupledCodeProperties->localCodeIs();
      int cond2 = !cond1 && (_localCodeProperties->idGet() < _coupledCodeProperties->idGet());

      if (!cond1 && !cond2) {
        return;
      }

      if (_tgt_in_src_idx != NULL) {
        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (_tgt_in_src_idx[i_part] != NULL) {
            free(_tgt_in_src_idx [i_part]);
            free(_tgt_in_src_gnum[i_part]);
            free(_tgt_in_src_dist[i_part]);
          }
        }

        free(_tgt_in_src_idx );
        free(_tgt_in_src_gnum);
        free(_tgt_in_src_dist);
        _tgt_in_src_idx  = NULL;
        _tgt_in_src_gnum = NULL;
        _tgt_in_src_dist = NULL;
      }

      if (_closest_src_gnum != NULL) {
        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (_closest_src_gnum[i_part] != NULL) {
            free(_closest_src_gnum [i_part]);
          }
        }

        free(_closest_src_gnum );
        _closest_src_gnum  = NULL;
      }


      if (_send_coord != NULL) {
        free(_send_coord);
        _send_coord = NULL;
      }

      if (_recv_coord != NULL) {
        for (int i_part = 0; i_part < _nPart; i_part++) {
          if (_recv_coord[i_part] != NULL) {
            free(_recv_coord[i_part]);
          }
        }
        free(_recv_coord);
        _recv_coord = NULL;
      }



      if (cond2) {
        SpatialInterpClosestPoint *cpl_spatial_interp;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
          cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }
        else {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
          cpl_spatial_interp =
            dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        if (cpl_spatial_interp->_tgt_in_src_idx != NULL) {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            if (cpl_spatial_interp->_tgt_in_src_idx[i_part] != NULL) {
              free(cpl_spatial_interp->_tgt_in_src_idx [i_part]);
              free(cpl_spatial_interp->_tgt_in_src_gnum[i_part]);
              free(cpl_spatial_interp->_tgt_in_src_dist[i_part]);
            }
          }

          free(cpl_spatial_interp->_tgt_in_src_idx );
          free(cpl_spatial_interp->_tgt_in_src_gnum);
          free(cpl_spatial_interp->_tgt_in_src_dist);
          cpl_spatial_interp->_tgt_in_src_idx  = NULL;
          cpl_spatial_interp->_tgt_in_src_gnum = NULL;
          cpl_spatial_interp->_tgt_in_src_dist = NULL;
        }


        if (cpl_spatial_interp->_closest_src_gnum != NULL) {
          for (int i_part = 0; i_part < _nPart; i_part++) {
            if (cpl_spatial_interp->_closest_src_gnum[i_part] != NULL) {
              free(cpl_spatial_interp->_closest_src_gnum [i_part]);
            }
          }

          free(cpl_spatial_interp->_closest_src_gnum);
          cpl_spatial_interp->_closest_src_gnum = NULL;
        }

        if (cpl_spatial_interp->_send_coord != NULL) {
          free(cpl_spatial_interp->_send_coord);
          cpl_spatial_interp->_send_coord = NULL;
        }

        if (cpl_spatial_interp->_recv_coord != NULL) {
          for (int i_part = 0; i_part < _cplNPart; i_part++) {
            if (cpl_spatial_interp->_recv_coord[i_part] != NULL) {
              free(cpl_spatial_interp->_recv_coord[i_part]);
            }
          }
          free(cpl_spatial_interp->_recv_coord);
          cpl_spatial_interp->_recv_coord = NULL;
        }

      }
    }


    void SpatialInterpClosestPoint::set_PDM_object() {

      int cond1 = !_coupledCodeProperties->localCodeIs();
      int cond2 = !cond1 && (_localCodeProperties->idGet() < _coupledCodeProperties->idGet());

      if (!cond1 && !cond2) {
        return;
      }

      SpatialInterpClosestPoint *cpl_spatial_interp = NULL;
      cwipi::Mesh *cpl_mesh = NULL;

      if (cond2) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();
          cpl_spatial_interp =
          dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }
        else {
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();
          cpl_spatial_interp =
          dynamic_cast <SpatialInterpClosestPoint *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        cpl_spatial_interp->_id_pdm = _id_pdm;
        cpl_mesh = cpl_cpl.meshGet();
      }


      _send_coord = (const double **) malloc(sizeof(double *) * _nPart);
      if (cond2) {
        cpl_spatial_interp->_send_coord = (const double **) malloc(sizeof(double *) * _cplNPart);
      }
      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
        for (int i_part = 0 ; i_part < _nPart; i_part++) {
          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            _send_coord[i_part] = _mesh->eltCentersGet(i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            _send_coord[i_part] = _mesh->getVertexCoords(i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            _send_coord[i_part] = _cpl->userTargetCoordsGet(i_part);
          }
        }
      }
      else {
        // TODO K: why do we do this?
        for (int i_part = 0 ; i_part < _nPart; i_part++) {
          _send_coord[i_part] = NULL;
        }

        if (cond2) {
          cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
          for (int i_part = 0 ; i_part < cpl_spatial_interp->_nPart; i_part++) {
            if (_coupledCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              cpl_spatial_interp->_send_coord[i_part] = cpl_mesh->eltCentersGet(i_part);
            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              cpl_spatial_interp->_send_coord[i_part] = cpl_mesh->getVertexCoords(i_part);
            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_USER) {
              cpl_spatial_interp->_send_coord[i_part] = cpl_cpl.userTargetCoordsGet(i_part);
            }
          }
        }
      }




      int n_neighbors = CWP_CLOSEST_POINTS_N_NEIGHBORS;
      std::map<std::string, int> prop = _cpl->SpatialInterpPropertiesIntGet();
      std::map<std::string, int>::iterator it;

      it = prop.find("n_neighbors");
      if (it != prop.end()) {
        n_neighbors = it->second;
      }

      _id_pdm = PDM_closest_points_create(_pdmUnionComm,
                                          n_neighbors,
                                          PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);


      int local_is_src =
      (_exchDirection == SPATIAL_INTERP_EXCH_SEND && !_reverse) ||
      (_exchDirection == SPATIAL_INTERP_EXCH_RECV &&  _reverse);


      if (local_is_src) {
        PDM_closest_points_n_part_cloud_set(_id_pdm,
                                            _nPart,
                                            _cplNPart);

        /* Source cloud (local) */
        for (int i_part = 0; i_part < _nPart; i_part++) {
          const double      *src_coord = NULL;
          const PDM_g_num_t *src_g_num = NULL;
          int n_src = 0;

          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            src_g_num = (const PDM_g_num_t *) _mesh->GNumEltsGet  (i_part);
            src_coord =                       _mesh->eltCentersGet(i_part);
            n_src     =                       _mesh->getPartNElts (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            src_g_num = (const PDM_g_num_t *) _mesh->getVertexGNum  (i_part);
            src_coord =                       _mesh->getVertexCoords(i_part);
            n_src     =                       _mesh->getPartNVertex (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            src_g_num = (const PDM_g_num_t *) _cpl->userTargetGNumGet  (i_part);
            src_coord =                       _cpl->userTargetCoordsGet(i_part);
            n_src     =                       _cpl->userTargetNGet     (i_part);
          }

          PDM_closest_points_src_cloud_set(_id_pdm,
                                           i_part,
                                           n_src,
                           (double      *) src_coord,
                           (PDM_g_num_t *) src_g_num);
        }

        /* Target cloud (distant) */
        for (int i_part = 0; i_part < _cplNPart; i_part++) {
          const double      *tgt_coord = NULL;
          const PDM_g_num_t *tgt_g_num = NULL;
          int n_tgt = 0;

          if (cond2) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
            if (_coupledCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              tgt_g_num = (const PDM_g_num_t *) cpl_mesh->GNumEltsGet  (i_part);
              tgt_coord =                       cpl_mesh->eltCentersGet(i_part);
              n_tgt     =                       cpl_mesh->getPartNElts (i_part);
            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              tgt_g_num = (const PDM_g_num_t *) cpl_mesh->getVertexGNum  (i_part);
              tgt_coord =                       cpl_mesh->getVertexCoords(i_part);
              n_tgt     =                       cpl_mesh->getPartNVertex (i_part);
            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_USER) {
              tgt_g_num = (const PDM_g_num_t *) cpl_cpl.userTargetGNumGet  (i_part);
              tgt_coord =                       cpl_cpl.userTargetCoordsGet(i_part);
              n_tgt     =                       cpl_cpl.userTargetNGet     (i_part);
            }
          }

          PDM_closest_points_tgt_cloud_set(_id_pdm,
                                           i_part,
                                           n_tgt,
                           (double      *) tgt_coord,
                           (PDM_g_num_t *) tgt_g_num);
        }
      }
      else {
        PDM_closest_points_n_part_cloud_set(_id_pdm,
                                            _cplNPart,
                                            _nPart);

        /* Source cloud (distant) */
        for (int i_part = 0; i_part < _cplNPart; i_part++) {
          const double      *src_coord = NULL;
          const PDM_g_num_t *src_g_num = NULL;
          int n_src = 0;

          if (cond2) {
            cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
            if (_coupledCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              src_g_num = (const PDM_g_num_t *) cpl_mesh->GNumEltsGet  (i_part);
              src_coord =                       cpl_mesh->eltCentersGet(i_part);
              n_src     =                       cpl_mesh->getPartNElts (i_part);
            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              src_g_num = (const PDM_g_num_t *) cpl_mesh->getVertexGNum  (i_part);
              src_coord =                       cpl_mesh->getVertexCoords(i_part);
              n_src     =                       cpl_mesh->getPartNVertex (i_part);
            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_USER) {
              src_g_num = (const PDM_g_num_t *) cpl_cpl.userTargetGNumGet  (i_part);
              src_coord =                       cpl_cpl.userTargetCoordsGet(i_part);
              n_src     =                       cpl_cpl.userTargetNGet     (i_part);
            }
          }

          PDM_closest_points_src_cloud_set(_id_pdm,
                                           i_part,
                                           n_src,
                           (double      *) src_coord,
                           (PDM_g_num_t *) src_g_num);
        }

        /* Target cloud (local) */
        for (int i_part = 0; i_part < _nPart; i_part++) {
          const double      *tgt_coord = NULL;
          const PDM_g_num_t *tgt_g_num = NULL;
          int n_tgt = 0;

          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            tgt_g_num = (const PDM_g_num_t *) _mesh->GNumEltsGet  (i_part);
            tgt_coord =                       _mesh->eltCentersGet(i_part);
            n_tgt     =                       _mesh->getPartNElts (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            tgt_g_num = (const PDM_g_num_t *) _mesh->getVertexGNum  (i_part);
            tgt_coord =                       _mesh->getVertexCoords(i_part);
            n_tgt     =                       _mesh->getPartNVertex (i_part);
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            tgt_g_num = (const PDM_g_num_t *) _cpl->userTargetGNumGet  (i_part);
            tgt_coord =                       _cpl->userTargetCoordsGet(i_part);
            n_tgt     =                       _cpl->userTargetNGet     (i_part);
          }

          PDM_closest_points_tgt_cloud_set(_id_pdm,
                                           i_part,
                                           n_tgt,
                           (double      *) tgt_coord,
                           (PDM_g_num_t *) tgt_g_num);
        }
      }
    }




};
