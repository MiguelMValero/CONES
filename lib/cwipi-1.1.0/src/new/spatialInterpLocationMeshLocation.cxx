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

#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "spatialInterpLocationMeshLocation.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "pdm_logging.h"

/**
 * \cond
 */

namespace cwipi {


  void
  SpatialInterpLocationMeshLocation::localization_init
  (
  )
  {

    if (!_coupledCodeProperties->localCodeIs()) {

      _id_pdm = PDM_mesh_location_create(1, _pdmUnionComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

      // PDM_mesh_location_reverse_results_enable (_id_pdm);

      PDM_mesh_location_method_set(_id_pdm, _location_method);

      double tolerance = CWP_MESH_LOCATION_BBOX_TOLERANCE;
      std::map<std::string, double> prop = _cpl->SpatialInterpPropertiesDoubleGet();
      std::map<std::string, double>::iterator it;

      it = prop.find("tolerance");
      if (it != prop.end()) {
        tolerance = it->second;
      }

      PDM_mesh_location_tolerance_set(_id_pdm, tolerance);

      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

        PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart);

      }

      else {

        PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _cplNPart);

      }

    }

    else {

      // Attention :
      //     - creation d'un objet unique pdm sur le plus petit des 2 id + Copie de l'id dans l'objet du code couple
      //     - rien a faire pour id le plus grand

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        _id_pdm = PDM_mesh_location_create(1, _pdmUnionComm, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        // PDM_mesh_location_reverse_results_enable (_id_pdm);

        PDM_mesh_location_method_set(_id_pdm, _location_method);

        double tolerance = CWP_MESH_LOCATION_BBOX_TOLERANCE;
        std::map<std::string, double> prop = _cpl->SpatialInterpPropertiesDoubleGet();
        std::map<std::string, double>::iterator it;

        it = prop.find("tolerance");
        if (it != prop.end()) {
          tolerance = it->second;
        }

        PDM_mesh_location_tolerance_set(_id_pdm, tolerance);

        SpatialInterpLocationMeshLocation *cpl_spatial_interp;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();

          cpl_spatial_interp =
            dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        else {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();

          cpl_spatial_interp =
            dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);
        }

        cpl_spatial_interp->_id_pdm = _id_pdm;

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _nPart);

        }

        else {

          PDM_mesh_location_n_part_cloud_set(_id_pdm, 0, _cplNPart);

        }
      }
    }
  }


  SpatialInterpLocationMeshLocation::~SpatialInterpLocationMeshLocation
  (
  )
  {

  }




  void
  SpatialInterpLocationMeshLocation::localization_points_cloud_setting
  (
  )
  {

    if (!_coupledCodeProperties->localCodeIs()) {

      if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

        for (int i_part = 0 ; i_part < _nPart ; i_part++) {

          const double      *part_coord = NULL;
          const PDM_g_num_t *part_gnum = NULL;
          int n_elt_size = 0;

          if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
            part_gnum  = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
            part_coord = _mesh->eltCentersGet (i_part);
            n_elt_size = _mesh->getPartNElts (i_part);
            // PDM_log_trace_array_double (part_coord, 3*n_elt_size, "cloud point cell center : ");
            // PDM_log_trace_array_double (_mesh->getVertexCoords (i_part), 3*_mesh->getPartNVertex (i_part), "cloud point cell center : ");
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
            part_gnum  = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
            part_coord = _mesh->getVertexCoords (i_part);
            n_elt_size = _mesh->getPartNVertex (i_part);
            // PDM_log_trace_array_double (part_coord, 3*n_elt_size, "cloud point node : ");
          }
          else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
            part_gnum  = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
            part_coord = _cpl->userTargetCoordsGet (i_part);
            n_elt_size = _cpl->userTargetNGet (i_part);
            // PDM_log_trace_array_double (part_coord, 3*n_elt_size, "cloud point user : ");
          }

          // printf("n_elt_size = %i\n", n_elt_size);

          PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, n_elt_size, (double *) part_coord, (PDM_g_num_t*) part_gnum);

        }
      }

      else {
        for (int i_part = 0 ; i_part < _cplNPart ; i_part++) {
          PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, 0, NULL, NULL);
        }
      }
    }

    else {

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          for (int i_part = 0 ; i_part < _nPart ; i_part++) {

            const double *part_coord = NULL;
            const PDM_g_num_t *part_gnum = NULL;
            int n_elt_size = 0;

            if (_localCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              part_gnum  = (const PDM_g_num_t *) _mesh->GNumEltsGet (i_part);
              part_coord = _mesh->eltCentersGet (i_part);
              n_elt_size = _mesh->getPartNElts (i_part);

            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              part_gnum  = (const PDM_g_num_t *) _mesh->getVertexGNum (i_part);
              part_coord = _mesh->getVertexCoords (i_part);
              n_elt_size = _mesh->getPartNVertex (i_part);
            }
            else if (_localCodeDofLocation == CWP_DOF_LOCATION_USER) {
              part_gnum  = (const PDM_g_num_t *) _cpl->userTargetGNumGet (i_part);
              part_coord = _cpl->userTargetCoordsGet (i_part);
              n_elt_size = _cpl->userTargetNGet (i_part);
            }

            PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, n_elt_size, (double *) part_coord, (PDM_g_num_t*) part_gnum);

          }

        }

        else {

          cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

          for (int i_part = 0 ; i_part < _cplNPart ; i_part++) {

            const double *part_coord = NULL;
            const PDM_g_num_t *part_gnum = NULL;
            int n_elt_size = 0;

            if (_coupledCodeDofLocation == CWP_DOF_LOCATION_CELL_CENTER) {
              part_gnum  = (const PDM_g_num_t *) cpl_mesh->GNumEltsGet (i_part);
              part_coord = cpl_mesh->eltCentersGet (i_part);
              n_elt_size = cpl_mesh->getPartNElts (i_part);

            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_NODE) {
              part_gnum  = (const PDM_g_num_t *) cpl_mesh->getVertexGNum (i_part);
              part_coord = cpl_mesh->getVertexCoords (i_part);
              n_elt_size = cpl_mesh->getPartNVertex (i_part);
            }
            else if (_coupledCodeDofLocation == CWP_DOF_LOCATION_USER) {
              part_gnum  = (const PDM_g_num_t *) cpl_cpl.userTargetGNumGet (i_part);
              part_coord = cpl_cpl.userTargetCoordsGet (i_part);
              n_elt_size = cpl_cpl.userTargetNGet (i_part);
            }

            PDM_mesh_location_cloud_set(_id_pdm, 0, i_part, n_elt_size, (double *) part_coord, (PDM_g_num_t*) part_gnum);
          }
        }
      }
    }
  }


  void
  SpatialInterpLocationMeshLocation::localization_surface_setting
  (
  )
  {

    if (!_coupledCodeProperties->localCodeIs()) {
      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {
        PDM_mesh_location_shared_nodal_mesh_set(_id_pdm, _pdm_CplNodal);
      }
    }

    else {
      // CWP_Interface_t interf_dim = _cpl->entitiesDimGet();

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

          PDM_mesh_location_shared_nodal_mesh_set(_id_pdm, _pdm_CplNodal);
        }

        else {
          cwipi::Coupling &cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());
          cwipi::Mesh *cpl_mesh = cpl_cpl.meshGet();

          PDM_mesh_location_shared_nodal_mesh_set(_id_pdm, cpl_mesh->getPdmNodalIndex());

        }
      }
    }
  }


  void
  SpatialInterpLocationMeshLocation::localization_compute
  (
  )
  {


    if (!_coupledCodeProperties->localCodeIs()) {

      PDM_mesh_location_compute(_id_pdm);

      int dump_times = 0;
      char *env_var = NULL;
      env_var = getenv ("CWP_DUMP_TIMES");
      if (env_var != NULL) {
        dump_times = atoi(env_var);
      }
      if (dump_times) {
        PDM_mesh_location_dump_times(_id_pdm);
      }

    }

    else {
      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        PDM_mesh_location_compute(_id_pdm);

        int dump_times = 0;
        char *env_var = NULL;
        env_var = getenv ("CWP_DUMP_TIMES");
        if (env_var != NULL) {
          dump_times = atoi(env_var);
        }
        if (dump_times) {
          PDM_mesh_location_dump_times(_id_pdm);
        }

      }
    }
  }


  void
  SpatialInterpLocationMeshLocation::localization_get
  (
  )
  {

    if (!_coupledCodeProperties->localCodeIs()) {

      if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

        for (int i_part = 0 ; i_part < _nPart ; i_part++) {


          PDM_mesh_location_cell_vertex_get (_id_pdm,
                                             i_part,
                                             &(_cell_vtx_idx[i_part]),
                                             &(_cell_vtx[i_part]));

          PDM_mesh_location_points_in_elt_get (_id_pdm,
                                               i_part,
                                               0,
                                               &(_elt_pts_inside_idx[i_part]),
                                               &(_points_gnum[i_part]),
                                               &(_points_coords[i_part]),
                                               &(_points_uvw[i_part]),
                                               &(_weights_idx[i_part]),
                                               &(_weights[i_part]),
                                               &(_points_dist2[i_part]),
                                               &(_points_projected_coords[i_part]));

          _n_involved_sources_tgt[i_part] = _elt_pts_inside_idx[i_part][_src_n_gnum[i_part]];
          _involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * _n_involved_sources_tgt[i_part]);

          int count = 0;
          for (int i = 0 ; i < _src_n_gnum[i_part] ; ++i) {
            if (_elt_pts_inside_idx[i_part][i + 1] > _elt_pts_inside_idx[i_part][i]) {
              _involved_sources_tgt[i_part][count] = i + 1;
              ++count;
            }
          }

          _n_involved_sources_tgt[i_part] = count;
          _involved_sources_tgt[i_part] = (int*) realloc(_involved_sources_tgt[i_part], sizeof(int) * count);

          int n_elt = 0;
          CWP_Interface_t interf_dim = _cpl->entitiesDimGet();

          if (interf_dim == CWP_INTERFACE_SURFACE) {
            n_elt = _mesh->getNFace(i_part);
          }
          else if (interf_dim == CWP_INTERFACE_VOLUME) {
            n_elt = _mesh->getNCell(i_part);
          }


          // for (int i = 0; i < n_elt; i++) {

          //   printf ("elt cwp %d :", i+1);
          //   for (int j = _cell_vtx_idx[i_part][i]; j < _cell_vtx_idx[i_part][i+1]; j++) {
          //     printf(" %d", _cell_vtx[i_part][j]);
          //   }
          //   printf("\n");

          //   for (int j = _elt_pts_inside_idx[i_part][i]; j < _elt_pts_inside_idx[i_part][i+1]; j++) {
          //     printf("pts inside :");
          //     for (int k = 0; k < 3; k++) {
          //       printf (" %12.5e", _points_coords[i_part][3*j+k]);
          //     }
          //     printf("\n");
          //     printf("weights :");
          //     for (int k = _weights_idx[i_part][j]; k < _weights_idx[i_part][j+1]; k++) {
          //       printf(" %12.5e", _weights[i_part][k]);
          //     }
          //     printf("\n");
          //   }
          // }


          _n_elt_weights[i_part] = _elt_pts_inside_idx[i_part][n_elt];

        }
      }

      else {

        for (int i_part = 0; i_part < _nPart; i_part++) {

          _n_computed_tgt[i_part] = PDM_mesh_location_n_located_get (_id_pdm,
                                                                     0,
                                                                     i_part);

          _n_uncomputed_tgt[i_part] = PDM_mesh_location_n_unlocated_get (_id_pdm,
                                                                         0,
                                                                         i_part);

          _computed_tgt[i_part] = PDM_mesh_location_located_get (_id_pdm,
                                                                 0,
                                                                 i_part);

          _uncomputed_tgt[i_part] = PDM_mesh_location_unlocated_get (_id_pdm,
                                                                     0,
                                                                     i_part);

          PDM_mesh_location_point_location_get (_id_pdm,
                                                0,
                                                i_part,
                                                &(_tgt_closest_elt_gnum[i_part]),
                                                &(_tgt_distance[i_part]),
                                                &(_tgt_projected[i_part]));

          // For pdm_part1_to_selected_part2
          _elt_pts_inside_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
          _elt_pts_inside_idx[i_part][0] = 0;

        }
      }
    }

    else {

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {
        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        SpatialInterpLocationMeshLocation *cpl_spatial_interp;

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();

          cpl_spatial_interp = dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        else {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();

          cpl_spatial_interp = dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        if (_exchDirection == SPATIAL_INTERP_EXCH_SEND) {

          for (int i_part = 0 ; i_part < _nPart ; i_part++) {

            PDM_mesh_location_cell_vertex_get (_id_pdm,
                                               i_part,
                                               &(_cell_vtx_idx[i_part]),
                                               &(_cell_vtx[i_part]));

            PDM_mesh_location_points_in_elt_get (_id_pdm,
                                                 i_part,
                                                 0,
                                                 &(_elt_pts_inside_idx[i_part]),
                                                 &(_points_gnum[i_part]),
                                                 &(_points_coords[i_part]),
                                                 &(_points_uvw[i_part]),
                                                 &(_weights_idx[i_part]),
                                                 &(_weights[i_part]),
                                                 &(_points_dist2[i_part]),
                                                 &(_points_projected_coords[i_part]));

            _n_involved_sources_tgt[i_part] = _elt_pts_inside_idx[i_part][_src_n_gnum[i_part]];
            _involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * _n_involved_sources_tgt[i_part]);

            int count = 0;
            for (int i = 0 ; i < _src_n_gnum[i_part] ; ++i) {
              if (_elt_pts_inside_idx[i_part][i + 1] > _elt_pts_inside_idx[i_part][i]) {
                _involved_sources_tgt[i_part][count] = i + 1;
                ++count;
              }
            }

            _n_involved_sources_tgt[i_part] = count;
            _involved_sources_tgt[i_part] = (int*) realloc(_involved_sources_tgt[i_part], sizeof(int) * count);

            int n_elt = 0;
            CWP_Interface_t interf_dim = _cpl->entitiesDimGet();

            if (interf_dim == CWP_INTERFACE_SURFACE) {
              n_elt = _mesh->getNFace(i_part);
            }
            else if (interf_dim == CWP_INTERFACE_VOLUME) {
              n_elt = _mesh->getNCell(i_part);
            }

            _n_elt_weights[i_part] = _elt_pts_inside_idx[i_part][n_elt];

          }

          for (int i_part = 0; i_part < _cplNPart; i_part++) {

            cpl_spatial_interp->_n_computed_tgt[i_part] = PDM_mesh_location_n_located_get (_id_pdm,
                                                                                           0,
                                                                                           i_part);


            cpl_spatial_interp->_n_uncomputed_tgt[i_part] = PDM_mesh_location_n_unlocated_get (_id_pdm,
                                                                                               0,
                                                                                               i_part);

            cpl_spatial_interp->_computed_tgt[i_part] = PDM_mesh_location_located_get (_id_pdm,
                                                                                       0,
                                                                                       i_part);

            cpl_spatial_interp->_uncomputed_tgt[i_part] = PDM_mesh_location_unlocated_get (_id_pdm,
                                                                                           0,
                                                                                           i_part);

            PDM_mesh_location_point_location_get (_id_pdm,
                                                  0,
                                                  i_part,
                                                  &(cpl_spatial_interp->_tgt_closest_elt_gnum[i_part]),
                                                  &(cpl_spatial_interp->_tgt_distance[i_part]),
                                                  &(cpl_spatial_interp->_tgt_projected[i_part]));

            // For pdm_part1_to_selected_part2
            cpl_spatial_interp->_elt_pts_inside_idx[i_part] = (int*) malloc (sizeof(int)); // Use malloc not new [] !
            cpl_spatial_interp->_elt_pts_inside_idx[i_part][0] = 0;

          }

        }

        else {

          for (int i_part = 0; i_part < _nPart; i_part++) {

            _n_computed_tgt[i_part] = PDM_mesh_location_n_located_get (_id_pdm,
                                                                       0,
                                                                       i_part);

            _n_uncomputed_tgt[i_part] = PDM_mesh_location_n_unlocated_get (_id_pdm,
                                                                           0,
                                                                           i_part);

            _computed_tgt[i_part] = PDM_mesh_location_located_get (_id_pdm,
                                                                   0,
                                                                   i_part);

            _uncomputed_tgt[i_part] = PDM_mesh_location_unlocated_get (_id_pdm,
                                                                       0,
                                                                       i_part);

            PDM_mesh_location_point_location_get (_id_pdm,
                                                  0,
                                                  i_part,
                                                  &(_tgt_closest_elt_gnum[i_part]),
                                                  &(_tgt_distance[i_part]),
                                                  &(_tgt_projected[i_part]));

            // For pdm_part1_to_selected_part2
            _elt_pts_inside_idx[i_part] = (int*)  malloc(sizeof(int)); // Use malloc not new [] !
            _elt_pts_inside_idx[i_part][0] = 0;


          }

          for (int i_part = 0 ; i_part < _cplNPart ; i_part++) {

            PDM_mesh_location_cell_vertex_get (_id_pdm,
                                               i_part,
                                               &(cpl_spatial_interp->_cell_vtx_idx[i_part]),
                                               &(cpl_spatial_interp->_cell_vtx[i_part]));

            PDM_mesh_location_points_in_elt_get (_id_pdm,
                                                 i_part,
                                                 0,
                                                 &(cpl_spatial_interp->_elt_pts_inside_idx[i_part]),
                                                 &(cpl_spatial_interp->_points_gnum[i_part]),
                                                 &(cpl_spatial_interp->_points_coords[i_part]),
                                                 &(cpl_spatial_interp->_points_uvw[i_part]),
                                                 &(cpl_spatial_interp->_weights_idx[i_part]),
                                                 &(cpl_spatial_interp->_weights[i_part]),
                                                 &(cpl_spatial_interp->_points_dist2[i_part]),
                                                 &(cpl_spatial_interp->_points_projected_coords[i_part]));

            cpl_spatial_interp->_n_involved_sources_tgt[i_part] = cpl_spatial_interp->_elt_pts_inside_idx[i_part][cpl_spatial_interp->_src_n_gnum[i_part]];
            cpl_spatial_interp->_involved_sources_tgt[i_part] = (int*) malloc(sizeof(int) * cpl_spatial_interp->_n_involved_sources_tgt[i_part]);

            int count = 0;
            for (int i = 0 ; i < cpl_spatial_interp->_src_n_gnum[i_part] ; ++i) {
              if (cpl_spatial_interp->_elt_pts_inside_idx[i_part][i + 1] > cpl_spatial_interp->_elt_pts_inside_idx[i_part][i]) {
                cpl_spatial_interp->_involved_sources_tgt[i_part][count] = i + 1;
                ++count;
              }
            }

            cpl_spatial_interp->_n_involved_sources_tgt[i_part] = count;
            cpl_spatial_interp->_involved_sources_tgt[i_part] = (int*) realloc(cpl_spatial_interp->_involved_sources_tgt[i_part], sizeof(int) * count);

            int n_elt = 0;
            CWP_Interface_t interf_dim = cpl_cpl.entitiesDimGet();

            if (interf_dim == CWP_INTERFACE_SURFACE) {
              n_elt = cpl_spatial_interp->_mesh->getNFace(i_part);
            }
            else if (interf_dim == CWP_INTERFACE_VOLUME) {
              n_elt = cpl_spatial_interp->_mesh->getNCell(i_part);
            }

            cpl_spatial_interp->_n_elt_weights[i_part] = cpl_spatial_interp->_elt_pts_inside_idx[i_part][n_elt];

          }

        }
      }
    }

    if (_id_pdm != NULL) {
      PDM_mesh_location_part_to_part_get(_id_pdm,
                                         0,
                                         &_ptsp,
                                         PDM_OWNERSHIP_USER);
    }
  }


  void SpatialInterpLocationMeshLocation::localization_free() {

    if (!_coupledCodeProperties->localCodeIs()) {
      PDM_mesh_location_free(_id_pdm);
      _id_pdm = nullptr;
    }

    else {

      if (_localCodeProperties->idGet() < _coupledCodeProperties->idGet()) {

        PDM_mesh_location_free(_id_pdm);
        _id_pdm = nullptr;

        cwipi::Coupling& cpl_cpl = _cpl->couplingDBGet()->couplingGet(*_coupledCodeProperties, _cpl->IdGet());

        SpatialInterpLocationMeshLocation *cpl_spatial_interp;

        if (_exchDirection == SPATIAL_INTERP_EXCH_RECV) {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_send_map = cpl_cpl.sendSpatialInterpGet();

          cpl_spatial_interp = dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_send_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }

        else {

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> &cpl_spatial_interp_recv_map = cpl_cpl.recvSpatialInterpGet();

          cpl_spatial_interp = dynamic_cast <SpatialInterpLocationMeshLocation *> (cpl_spatial_interp_recv_map[make_pair(_coupledCodeDofLocation, _localCodeDofLocation)]);

        }
        cpl_spatial_interp->_id_pdm = nullptr;
      }
    }
  }
}

/**
 * \endcond
 */
