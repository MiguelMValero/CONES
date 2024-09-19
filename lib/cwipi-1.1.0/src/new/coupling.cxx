/*  This file is part of the CWIPI library.

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

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstring>
#include <sstream>

#include <pdm_error.h>

#include "coupling.hxx"
#include "coupling_i.hxx"

#include "mesh.hxx"
#include "codeProperties.hxx"

#include "solve_ax_b_4.h"
#include "quickSort.h"
#include "cwp.h"

#include "factory.hpp"
#include "field.hxx"
#include "globalData.hxx"
#include "partData.hxx"
#include "pdm_part_to_part.h"
#include "pdm_array.h"

#include "communication.hxx"
// #include "visualization.hxx"
#include "pdm_writer.h"
#include "pdm_logging.h"

#include "spatialInterp.hxx"
#include "spatialInterpLocation.hxx"
#include "spatialInterpClosestPoint.hxx"
#include "spatialInterpIntersection.hxx"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

/**
 * \cond
 */

#if !defined (__hpux) &&  !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

using namespace std;

namespace cwipi {

  /**
   * \typedef FC
   *
   * \brief Communication Factory
   *
   *  A communication \ref Factory wich makes \ref Communication
   *  class objects.
   *  The type of communication objects build depends on the
   *  communication type \ref CWP_Comm_t .
   *
   */

  typedef Factory<Communication, CWP_Comm_t> FC;

   /**
   * \typedef FG
   *
   * \brief Spatial Interpolation Factory
   *
   *  A Spatial Interpolation \ref Factory wich makes  \ref SpatialInterp class objects.
   *  The type of \ref SpatialInterp objects build depends on the
   *  spatial interpolation algorithm type \ref CWP_Spatial_interp_t; .
   *
   */

  typedef Factory<SpatialInterp, CWP_Spatial_interp_t> FG;


  /**
   * \brief Constructor.
   *
   * This function creates a coupling object and defines its properties.
   *
   * \param [in]  cplId                        Coupling identifier
   * \param [in]  commType                     Communication type
   * \param [in]  localCodeProperties          Local code properties
   * \param [in]  coupledCodeProperties        Coupled code properties
   * \param [in]  spatialInterpAlgo                     SpatialInterp algorithm
   * \param [in]  nPart                        Number of interface partitions
   * \param [in]  movingStatus                 Mesh moving status
   * \param [in]  recvFreqType                 Type of receiving frequency
   * \param [in]  cplDB                        Coupling data base where it coupling is stored
   *
   */

  Coupling::Coupling
  (
   const string               &cplId,
   const CWP_Comm_t           cplType,
         CodeProperties       &localCodeProperties,
         CodeProperties       &coupledCodeProperties,
   const CWP_Interface_t      entities_dim,
   const CWP_Spatial_interp_t spatialInterpAlgo,
   const int                  nPart,
   const CWP_Dynamic_mesh_t   displacement,
   const CWP_Time_exch_t      recvFreqType,
   CouplingDB                 &cplDB
   )
  :_cplId(cplId),
   _commType(cplType),
   _communication(*(FC::getInstance().CreateObject(cplType))),
   _localCodeProperties(localCodeProperties),
   _coupledCodeProperties(coupledCodeProperties),
   _entities_dim(entities_dim),
   _mesh(*new Mesh(localCodeProperties.connectableCommGet(),nPart,displacement,this)),
   _recvFreqType (recvFreqType),
   _id_geom_writer(-1),
   _id_field_partitioning_writer(-1),
   _id_field_ranking_writer(-1),
   _id_user_tgt_geom_writer(-1),
   _id_user_tgt_field_partitioning_writer(-1),
   _id_user_tgt_field_ranking_writer(-1),
   _freq_writer(-1),
   _writer(nullptr),
   _fields(*(new map < string, Field * >())),
   _globalData(*(new map < string, GlobalData >())),
   _partData(*(new map < string, PartData >())),
   _cplDB(cplDB),
   _displacement(displacement),
   _spatialInterpAlgo(spatialInterpAlgo),
   _nPart(nPart),
   _cplNPart(-1),
   _userTargetN(nullptr),
   _userTargetGnum(nullptr),
   _localUserTargetGnum(nullptr),
   _userTargetCoord(nullptr),
   _spatial_interp_send(*new std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>()),
   _spatial_interp_recv(*new std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>()),
   _spatial_interp_properties_double(*new std::map<std::string, double>),
   _spatial_interp_properties_int(*new std::map<std::string, int>),
   _is_mesh_finalized(0),
   _is_first_field_created(0),
   _n_step(0)

   {

     //In case where the both codes are on the same MPI process.
    if (coupledCodeProperties.localCodeIs()) {
      if (cplDB.couplingIs(coupledCodeProperties, cplId) ) {
        const int localRootRank = localCodeProperties.rootRankGet();
        const int cplRootRank = coupledCodeProperties.rootRankGet();

        const MPI_Comm& globalComm = localCodeProperties.globalCommGet();

        int globalRank;
        MPI_Comm_rank(globalComm, &globalRank);

        Coupling &distCpl = cplDB.couplingGet(coupledCodeProperties, cplId);

        if ((cplRootRank == globalRank) && (cplRootRank != localRootRank)) {
          distCpl._communication.init(_coupledCodeProperties, _localCodeProperties, cplId, cplDB);
          _communication.init(distCpl._communication);
        }
        else {
          _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);
          distCpl._communication.init(_communication);
        }

      }
    } // if (coupledCodeProperties.localCodeIs())

    else {
      //Communication initialization, MPI communicator creation ...
      _communication.init(_localCodeProperties, _coupledCodeProperties, cplId, cplDB);

    } // end else

    //entitiesDimGet();

  }


  /**
   * \brief Destructor.
   *
   */

  Coupling::~Coupling()
  {
    delete &_spatial_interp_properties_double;
    delete &_spatial_interp_properties_int;

    std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t > , SpatialInterp*>::iterator it = _spatial_interp_send.begin();
    while (it != _spatial_interp_send.end()) {
        delete it->second;
        it++;
    }

    delete &_spatial_interp_send;

    it = _spatial_interp_recv.begin();
    while (it != _spatial_interp_recv.end()) {
       delete it->second;
        it++;
    }

    delete &_spatial_interp_recv;


    std::map < string, Field * >::iterator itf = _fields.begin();
    while (itf != _fields.end()) {
      if (itf->second != NULL) delete itf->second;
      itf++;
    }

    delete &_fields;

    delete &_globalData;

    delete &_partData;

    delete &_communication;

    delete &_mesh;

    if (_userTargetN != nullptr) {
      if (_localUserTargetGnum != nullptr) {
        for (int iPart = 0; iPart < _nPart; iPart++) {
          if (_localUserTargetGnum[iPart] != nullptr) free (_localUserTargetGnum[iPart]);
        }
      }
      free ( _userTargetN);
      free ( _userTargetGnum);
      free ( _userTargetCoord);
    }

    #if defined(DEBUG) && 0
    cout << "destroying '" << _name << "' coupling : TODO" << endl;
    #endif
  }

  /*----------------------------------------------------------------------------*
   * Methods about part data                                                    *
   *----------------------------------------------------------------------------*/

  /**
   * \brief Check if object already exists
   *
   * \param [in] part_data_id
   *
   */

  bool
  Coupling::partDataIs (
   const string &part_data_id
  )
  {
    map<string,PartData>::iterator it = _partData.find(part_data_id.c_str());
    return (it != _partData.end());
  }

  /**
   * \brief Create partitionned data exchange object
   *
   * \param [in] part_data_id
   * \param [in] exch_type
   * \param [in] gnum_elt
   * \param [in] n_elt
   * \param [in] n_part
   *
   */

  void
  Coupling::partDataCreate
  (
   const string          &part_data_id,
   CWP_PartData_exch_t   exch_type,
   CWP_g_num_t         **gnum_elt,
   int                  *n_elt,
   int                   n_part
  )
  {
    map<string, PartData>::iterator it = _partData.find(part_data_id.c_str());
    if (it != _partData.end()) {
      PDM_error(__FILE__, __LINE__, 0, "'%s' partitionned data exchange object already exists\n", part_data_id.c_str());
    }

    // Create object
    cwipi::PartData newPartData(part_data_id,
                                exch_type,
                                gnum_elt,
                                n_elt,
                                n_part);


    // Create ptp
    MPI_Comm unionComm = _communication.unionCommGet();
    PDM_part_to_part_t  *ptp                = NULL;
    int                **part1_to_part2_idx = NULL;

    /* Overlapping intracomms */
    if (_coupledCodeProperties.localCodeIs()) {
      cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);
      map<string, PartData>::iterator cpl_it = cpl_cpl._partData.find(part_data_id.c_str());

      if (cpl_it != cpl_cpl._partData.end()) {
        if (exch_type == CWP_PARTDATA_SEND) {
          CWP_g_num_t **gnum_elt2 = cpl_it->second.gnum_elt_get();
          int          *n_elt2    = cpl_it->second.n_elt_get   ();
          int           n_part2   = cpl_it->second.n_part_get  ();

          part1_to_part2_idx = (int **) malloc(sizeof(int *) * n_part);
          for (int i_part = 0; i_part < n_part; i_part++) {
            part1_to_part2_idx[i_part] = PDM_array_new_idx_from_const_stride_int(1,
                                                                                 n_elt[i_part]);
          }

          ptp = PDM_part_to_part_create((const PDM_g_num_t **) gnum_elt,
                                        n_elt,
                                        n_part,
                                        (const PDM_g_num_t **) gnum_elt2,
                                        n_elt2,
                                        n_part2,
                                        (const int         **) part1_to_part2_idx,
                                        (const PDM_g_num_t **) gnum_elt,
                                        PDM_MPI_mpi_2_pdm_mpi_comm(&unionComm));

          /* part1_to_part2_idx gets copied in PDM_part_to_part_create so we don't need it anymore */
          for (int i_part = 0; i_part < n_part; i_part++) {
            free(part1_to_part2_idx[i_part]);
          }
          free(part1_to_part2_idx);

        }

        else {
          CWP_g_num_t **gnum_elt1 = cpl_it->second.gnum_elt_get();
          int          *n_elt1    = cpl_it->second.n_elt_get   ();
          int           n_part1   = cpl_it->second.n_part_get  ();

          part1_to_part2_idx = (int **) malloc(sizeof(int *) * n_part1);
          for (int i_part = 0; i_part < n_part1; i_part++) {
            part1_to_part2_idx[i_part] = PDM_array_new_idx_from_const_stride_int(1,
                                                                                 n_elt1[i_part]);
          }

          ptp = PDM_part_to_part_create((const PDM_g_num_t **) gnum_elt1,
                                        n_elt1,
                                        n_part1,
                                        (const PDM_g_num_t **) gnum_elt,
                                        n_elt,
                                        n_part,
                                        (const int         **) part1_to_part2_idx,
                                        (const PDM_g_num_t **) gnum_elt1,
                                        PDM_MPI_mpi_2_pdm_mpi_comm(&unionComm));

          /* part1_to_part2_idx gets copied in PDM_part_to_part_create so we don't need it anymore */
          for (int i_part = 0; i_part < n_part1; i_part++) {
            free(part1_to_part2_idx[i_part]);
          }
          free(part1_to_part2_idx);

          int  *n_ref;
          int **ref;
          for (int i_part = 0; i_part < n_part; i_part++) {
            PDM_part_to_part_ref_lnum2_get(ptp,
                                           &n_ref,
                                           &ref);
            if (n_ref[i_part] != n_elt[i_part]) {
              PDM_error(__FILE__, __LINE__, 0, "Error in input data of CWP_Part_data_create : make sure the global numbering is continous\n");
            }
          }

        }

        cpl_it->second.ptp_set(ptp);
      }
    }

    /* Disjoint intracomms */
    else {

      if (exch_type == CWP_PARTDATA_SEND) {
        part1_to_part2_idx = (int **) malloc(sizeof(int *) * n_part);
        for (int i_part = 0; i_part < n_part; i_part++) {
          part1_to_part2_idx[i_part] = PDM_array_new_idx_from_const_stride_int(1,
                                                                               n_elt[i_part]);
        }

        ptp = PDM_part_to_part_create((const PDM_g_num_t **) gnum_elt,
                                      n_elt,
                                      n_part,
                                      NULL,
                                      NULL,
                                      0,
                                      (const int         **) part1_to_part2_idx,
                                      (const PDM_g_num_t **) gnum_elt,
                                      PDM_MPI_mpi_2_pdm_mpi_comm(&unionComm));

        /* part1_to_part2_idx gets copied in PDM_part_to_part_create so we don't need it anymore */
        for (int i_part = 0; i_part < n_part; i_part++) {
          free(part1_to_part2_idx[i_part]);
        }
        free(part1_to_part2_idx);
      }

      else {
        ptp = PDM_part_to_part_create(NULL,
                                      NULL,
                                      0,
                                      (const PDM_g_num_t **) gnum_elt,
                                      n_elt,
                                      n_part,
                                      NULL,
                                      NULL,
                                      PDM_MPI_mpi_2_pdm_mpi_comm(&unionComm));

        int  *n_ref;
        int **ref;
        PDM_part_to_part_ref_lnum2_get(ptp,
                                       &n_ref,
                                       &ref);
        for (int i_part = 0; i_part < n_part; i_part++) {
          if (n_ref[i_part] != n_elt[i_part]) {
            PDM_error(__FILE__, __LINE__, 0, "Error in input data of CWP_Part_data_create : make sure the global numbering is continous\n");
          }
        }
      }

    }



    newPartData.ptp_set(ptp);


    pair<string, PartData> newPair(part_data_id, newPartData);
    _partData.insert(newPair);
  }




  /**
   * \brief Delete partitionned data exchange object
   *
   * \param [in] part_data_id
   *
   */

  void
  Coupling::partDataDel
  (
   const string          &part_data_id
  )
  {
    map<string,PartData>::iterator it = _partData.find(part_data_id.c_str());
    if (it == _partData.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing partitionned data exchange object\n", part_data_id.c_str());
    }

    // free ptp
    PDM_part_to_part_t *ptp = it->second.ptp_get();

    if (!_coupledCodeProperties.localCodeIs() ||
        _localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {
      PDM_part_to_part_free(ptp);
    }

    // remove from map
    _partData.erase(part_data_id.c_str());

  }



  /**
   * \brief Issend partitionned data
   *
   * \param [in] part_data_id
   * \param [in] exch_id
   * \param [in] s_data
   * \param [in] n_components
   * \param [in] send_data
   *
   */

  void
  Coupling::partDataIssend
  (
   const string  &part_data_id,
   const int      exch_id,
         size_t   s_data,
         int      n_components,
         void   **send_data
  )
  {
    map<string, PartData>::iterator it = _partData.find(part_data_id.c_str());
    if (it == _partData.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Invalid PartData '%s'\n", part_data_id.c_str());
    }

    if (it->second.exch_type_get() != CWP_PARTDATA_SEND) {
      PDM_error(__FILE__, __LINE__, 0,
                "PartData '%s' is not in SEND mode for %s\n",
                part_data_id.c_str(),
                _localCodeProperties.nameGet().c_str());
    }

    MPI_Comm unionComm = _communication.unionCommGet();
    PDM_part_to_part_t *ptp = it->second.ptp_get();

    int mpi_tag = _communication._get_tag(part_data_id,
                                          unionComm,
                                          exch_id);

    int send_request = -1;

    PDM_part_to_part_issend(ptp,
                            s_data,
                            n_components,
            (const void **) send_data,
                            mpi_tag,
                            &send_request);

    it->second.data_set(exch_id,
                        send_data);

    it->second.request_set(exch_id,
                           send_request);
  }

  /**
   * \brief Irecv partitionned data
   *
   * \param [in] part_data_id
   * \param [in] exch_id
   * \param [in] s_data
   * \param [in] n_components
   * \param [in] recv_data
   *
   */

  void
  Coupling::partDataIrecv
  (
   const string  &part_data_id,
   const int      exch_id,
         size_t   s_data,
         int      n_components,
         void   **recv_data
  )
  {
    map<string, PartData>::iterator it = _partData.find(part_data_id.c_str());
    if (it == _partData.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Invalid PartData '%s'\n", part_data_id.c_str());
    }

    if (it->second.exch_type_get() != CWP_PARTDATA_RECV) {
      PDM_error(__FILE__, __LINE__, 0,
                "PartData '%s' is not in RECV mode for %s\n",
                part_data_id.c_str(),
                _localCodeProperties.nameGet().c_str());
    }


    MPI_Comm unionComm = _communication.unionCommGet();
    PDM_part_to_part_t *ptp = it->second.ptp_get();

    int mpi_tag = _communication._get_tag(part_data_id,
                                          unionComm,
                                          exch_id);

    int recv_request = -1;

    PDM_part_to_part_irecv(ptp,
                           s_data,
                           n_components,
                           recv_data,
                           mpi_tag,
                           &recv_request);

    it->second.data_set(exch_id,
                        recv_data);

    it->second.request_set(exch_id,
                           recv_request);

  }

  /**
   * \brief Wait issend partitionned data
   *
   * \param [in] part_data_id
   * \param [in] exch_id
   *
   */

  void
  Coupling::partDataWaitIssend
  (
   const string   &part_data_id,
   const int       exch_id
   )
  {
    map<string, PartData>::iterator it = _partData.find(part_data_id.c_str());
    if (it == _partData.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Invalid PartData '%s'\n", part_data_id.c_str());
    }

    if (it->second.exch_type_get() != CWP_PARTDATA_SEND) {
      PDM_error(__FILE__, __LINE__, 0,
                "PartData '%s' is not in SEND mode for %s\n",
                part_data_id.c_str(),
                _localCodeProperties.nameGet().c_str());
    }

    PDM_part_to_part_t *ptp = it->second.ptp_get();

    int request = it->second.request_get(exch_id);

    PDM_part_to_part_issend_wait(ptp,
                                 request);

    it->second.exch_clear(exch_id);
  }

  /**
   * \brief Wait irecv partitionned data
   *
   * \param [in] part_data_id
   * \param [in] exch_id
   *
   */

   void
  Coupling::partDataWaitIrecv
  (
   const string   &part_data_id,
   const int       exch_id
   )
  {
    map<string, PartData>::iterator it = _partData.find(part_data_id.c_str());
    if (it == _partData.end()) {
      PDM_error(__FILE__, __LINE__, 0, "Invalid PartData '%s'\n", part_data_id.c_str());
    }

    if (it->second.exch_type_get() != CWP_PARTDATA_RECV) {
      PDM_error(__FILE__, __LINE__, 0,
                "PartData '%s' is not in RECV mode for %s\n",
                part_data_id.c_str(),
                _localCodeProperties.nameGet().c_str());
    }

    PDM_part_to_part_t *ptp = it->second.ptp_get();

    int request = it->second.request_get(exch_id);

    PDM_part_to_part_irecv_wait(ptp,
                                request);

    it->second.recv_data_filter(exch_id);

    it->second.exch_clear(exch_id);
  }


  int
  Coupling::partDataNPartGet
  (
   const string   &part_data_id
   )
  {
    map<string,PartData>::iterator it = _partData.find(part_data_id.c_str());
    assert(it != _partData.end());

    return it->second.n_part_get();
  }

  int
  Coupling::partDataNRefGet
  (
   const string   &part_data_id,
   const int       i_part
   )
  {
    map<string,PartData>::iterator it = _partData.find(part_data_id.c_str());
    assert(it != _partData.end());

    assert(i_part < it->second.n_part_get());

    PDM_part_to_part_t *ptp = it->second.ptp_get();

    int  *n_ref;
    int **ref;
    PDM_part_to_part_ref_lnum2_get(ptp,
                                   &n_ref,
                                   &ref);

    return n_ref[i_part];
  }

  /*----------------------------------------------------------------------------*
   * Methods about global data                                                  *
   *----------------------------------------------------------------------------*/

  /**
   * \brief Send a data array.
   */

  void
  Coupling::globalDataIssend
  (
   const string    &global_data_id,
   size_t          s_send_entity,
   int             send_stride,
   int             n_send_entity,
   void           *send_data
   )
  {
    // Create an instance of GlobalData
    map<string,GlobalData>::iterator it = _globalData.find(global_data_id.c_str());
    if (it == _globalData.end()) {
      cwipi::GlobalData newGlobalData(global_data_id,
                                      s_send_entity,
                                      send_stride,
                                      n_send_entity);
      newGlobalData.send_data_set(send_data);

      pair<string, GlobalData > newPair(global_data_id, newGlobalData);
      _globalData.insert(newPair);
    } // end if does not exist
    else {
      PDM_error(__FILE__, __LINE__, 0, "GlobalData %s is already in use\n", global_data_id.c_str());
    }
    it = _globalData.find(global_data_id.c_str());

    assert(it != _globalData.end());

    MPI_Request data_send_request;
    _communication.issendGlobalDataBetweenCodesThroughUnionCom(global_data_id,
                                                               &data_send_request,
                                                               s_send_entity,
                                                               send_stride,
                                                               n_send_entity,
                                                               send_data);
    it->second.data_send_request_set(data_send_request);
  }

  /**
   * \brief Receive a data array.
   */

  void
  Coupling::globalDataIrecv
  (
   const string    &global_data_id,
   size_t          s_recv_entity,
   int             recv_stride,
   int             n_recv_entity,
   void           *recv_data
  )
  {
    // Create an instance of GlobalData
    map<string,GlobalData>::iterator it = _globalData.find(global_data_id.c_str());
    if (it == _globalData.end()) {
      cwipi::GlobalData newGlobalData(global_data_id,
                                      s_recv_entity,
                                      recv_stride,
                                      n_recv_entity);
      newGlobalData.recv_data_set(recv_data);

      pair<string, GlobalData > newPair(global_data_id, newGlobalData);
      _globalData.insert(newPair);
    } // end if does not exist
    else {
      PDM_error(__FILE__, __LINE__, 0, "GlobalData %s is already in use\n");
    }
    it = _globalData.find(global_data_id.c_str());

    assert(it != _globalData.end());

    MPI_Request data_recv_request;
    _communication.irecvGlobalDataBetweenCodesThroughUnionCom(global_data_id,
                                                              &data_recv_request,
                                                              s_recv_entity,
                                                              recv_stride,
                                                              n_recv_entity,
                                                              recv_data);
    it->second.data_recv_request_set(data_recv_request);
  }

  /**
   * \brief Wait of send a data array.
   */

  void
  Coupling::globalDataWaitIssend
  (
   const string    &global_data_id
  )
  {
    // Get local
    map<string,GlobalData>::iterator it = _globalData.find(global_data_id.c_str());
    assert(it != _globalData.end());
    MPI_Request  data_send_request = it->second.data_send_request_get();
    size_t       s_entity          = it->second.s_entity_get();
    int          stride            = it->second.stride_get();
    int          n_entity          = it->second.n_entity_get();
    void        *send_data         = it->second.send_data_get();

    void *cpl_recv_data = NULL;
    if (_coupledCodeProperties.localCodeIs()) {
      cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);
      map<string,GlobalData>::iterator cpl_it = cpl_cpl._globalData.find(global_data_id.c_str());
      assert(cpl_it != cpl_cpl._globalData.end());
      cpl_recv_data = cpl_it->second.recv_data_get();
    }

    // Get coupled
    _communication.waitIssendGlobalDataBetweenCodesThroughUnionCom(data_send_request,
                                                                   s_entity,
                                                                   stride,
                                                                   n_entity,
                                                                   send_data,
                                                                   cpl_recv_data);
  }

  /**
   * \brief Wait of receive a data array.
   */

  void
  Coupling::globalDataWaitIrecv
  (
   const string    &global_data_id
   )
  {
    // Get local
    map<string,GlobalData>::iterator it = _globalData.find(global_data_id.c_str());
    assert(it != _globalData.end());
    MPI_Request  data_recv_request = it->second.data_recv_request_get();
    size_t       s_entity          = it->second.s_entity_get();
    int          stride            = it->second.stride_get();
    int          n_entity          = it->second.n_entity_get();
    void        *recv_data         = it->second.recv_data_get();

    void *cpl_send_data = NULL;
    if (_coupledCodeProperties.localCodeIs()) {
      cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);
      map<string,GlobalData>::iterator cpl_it = cpl_cpl._globalData.find(global_data_id.c_str());
      assert(cpl_it != cpl_cpl._globalData.end());
      cpl_send_data = cpl_it->second.send_data_get();
    }

    // Get coupled
    _communication.waitIrecvGlobalDataBetweenCodesThroughUnionCom(data_recv_request,
                                                                  s_entity,
                                                                  stride,
                                                                  n_entity,
                                                                  cpl_send_data,
                                                                  recv_data);
  }

  /*----------------------------------------------------------------------------*
   * Methods about communicators                                                *
   *----------------------------------------------------------------------------*/

 /**
  * \brief Get coupling communicator and coupling ranks.
  *
  * \param [out] cpl_comm             Coupling communicator
  * \param [out] cpl_ranks            Coupling ranks
  *
  * \return Size of \ref cpl_ranks vector
  *
  */

  int
  Coupling::commGet (
    MPI_Comm  *cpl_comm,
    int      **cpl_ranks
  )
  {
    *cpl_comm  = _communication.cplCommGet();

    std::vector<int>* vect_cpl_ranks = _communication.cplCommCplRanksGet();
    *cpl_ranks = vect_cpl_ranks->data();

    return vect_cpl_ranks->size();
  }


  /*----------------------------------------------------------------------------*
   * Methods about exchange frequency                                           *
   *----------------------------------------------------------------------------*/




  /*----------------------------------------------------------------------------*
   * Methods about spatial interpolation                                        *
   *----------------------------------------------------------------------------*/


 /**
  * \brief Set a spatial interpolation property of type double.
  *
  * \param [in]   name       Name of the property
  * \param [in]   value      Value of the property
  *
  */
 void
  Coupling::spatialInterpPropertyDoubleSet (
    std::string name,
    double      value
  )
{
  std::map<std::string, double>::iterator it;
  it = _spatial_interp_properties_double.find(name);

  if (it == _spatial_interp_properties_double.end()) {
    pair<string, double> newPair(name, value);
    _spatial_interp_properties_double.insert(newPair);
  } else {
    it->second = value;
  }
}


/**
  * \brief Set a spatial interpolation property of type int.
  *
  * \param [in]   name       Name of the property
  * \param [in]   value      Value of the property
  *
  */
 void
  Coupling::spatialInterpPropertyIntSet (
    std::string name,
    int         value
  )
  {
    std::map<std::string, int>::iterator it;
    it = _spatial_interp_properties_int.find(name);

    if (it == _spatial_interp_properties_int.end()) {
      pair<string, int> newPair(name, value);
      _spatial_interp_properties_int.insert(newPair);
    } else {
      it->second = value;
    }
  }

  /**
   * \brief Get spatial interpolation algorithm enum.
   */

  CWP_Spatial_interp_t
  Coupling::spatialInterpAlgoGet()
  {
    return _spatialInterpAlgo;
  }

  /**
   * \brief Getters for callbacks.
   */

  // SpatialInterp

  // Get weights
  void
  Coupling::weight_get
  (
    std::string    name,
    int         ***weights_idx,
    double      ***weights
  )
  {
    map <string, Field *>::iterator it  = _fields.find(name);
    if (it != _fields.end()) {
      CWP_Dof_location_t localFieldLocation = it->second->locationGet();
      CWP_Dof_location_t cplFieldLocation = it->second->linkedFieldLocationGet();

      std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*> spatial_interp_map;

      // send or recv?
      CWP_Spatial_interp_t spatial_interp_algo = spatialInterpAlgoGet();
      if (spatial_interp_algo == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT ||
          spatial_interp_algo == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE         ||
          spatial_interp_algo == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE) {
        spatial_interp_map = _spatial_interp_send;
      }
      else if (spatial_interp_algo == CWP_SPATIAL_INTERP_FROM_INTERSECTION                  ||
               spatial_interp_algo == CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES ||
               spatial_interp_algo == CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES) {
        spatial_interp_map = _spatial_interp_recv;
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "Invalid spatial_interp_algo %d\n", spatial_interp_algo);
      }

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2 = spatial_interp_map.find(newKey);

      if (it2 != spatial_interp_map.end()) {
        *weights_idx = it2->second->weights_idx_get();
        *weights     = it2->second->weights_get();
      }
    }
  }

  // Get source ptp data
  void
  Coupling::src_data_get
  (
    std::string    name,
    int           *n_part_src,
    int          **n_elt_src,
    int         ***src_to_tgt_idx,
    CWP_g_num_t ***src_to_tgt_gnum
  )
  {
    map <string, Field *>::iterator it  = _fields.find(name);
    if (it != _fields.end()) {
      CWP_Dof_location_t localFieldLocation = it->second->locationGet();
      CWP_Dof_location_t cplFieldLocation = it->second->linkedFieldLocationGet();

      std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2 = _spatial_interp_send.find(newKey);

      if (it2 != _spatial_interp_send.end()) {
        PDM_part_to_part_t *ptp = it2->second->ptp_get();

        int  n_part_tgt = 0;
        int *n_elt_tgt  = NULL;
        PDM_part_to_part_n_part_and_n_elt_get(ptp,
                                              n_part_src,
                                              &n_part_tgt,
                                              n_elt_src,
                                              &n_elt_tgt);

        PDM_part_to_part_part1_to_part2_get(ptp,
                                            n_elt_src,
                                            src_to_tgt_idx,
                          (PDM_g_num_t ***) src_to_tgt_gnum);
      }
    }
  }

  // Get Target
  void
  Coupling::tgt_data_get
  (
    std::string     name,
    int            *n_part_tgt,
    int           **n_elt_tgt,
    int           **n_referenced_tgt,
    int          ***referenced_tgt,
    int          ***tgt_come_from_src_idx,
    CWP_g_num_t ***tgt_come_from_src
  )
  {
    map <string, Field *>::iterator it  = _fields.find(name);
    if (it != _fields.end()) {
      CWP_Dof_location_t localFieldLocation = it->second->locationGet();
      CWP_Dof_location_t cplFieldLocation = it->second->linkedFieldLocationGet();

      std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2 = _spatial_interp_recv.find(newKey);

      if (it2 != _spatial_interp_recv.end()) {
        PDM_part_to_part_t *ptp = it2->second->ptp_get();

        int  n_part_src = 0;
        int *n_elt_src  = NULL;
        PDM_part_to_part_n_part_and_n_elt_get(ptp,
                                              &n_part_src,
                                              n_part_tgt,
                                              &n_elt_src,
                                              n_elt_tgt);

        PDM_part_to_part_ref_lnum2_get(ptp,
                                       n_referenced_tgt,
                                       referenced_tgt);

        PDM_part_to_part_gnum1_come_from_get(ptp,
                                             tgt_come_from_src_idx,
                           (PDM_g_num_t ***) tgt_come_from_src);

      }
    }
  }

  // SpatialInterpLocation

  // Get point_*
  void
  Coupling::location_point_data_get
  (
   std::string      name,
   double       ***points_coords,
   double       ***points_uvw,
   double       ***points_dist2,
   double       ***points_projected_coords
  )
  {
    map <string, Field *>::iterator it  = _fields.find(name);
    if (it != _fields.end()) {
      CWP_Dof_location_t localFieldLocation = it->second->locationGet();
      CWP_Dof_location_t cplFieldLocation = it->second->linkedFieldLocationGet();

      std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2 = _spatial_interp_send.find(newKey);

      if (it2 != _spatial_interp_send.end()) {
        SpatialInterpLocation* sil = dynamic_cast <SpatialInterpLocation*> (it2->second);
        *points_coords           = sil->points_coords_get();
        *points_uvw              = sil->points_uvw_get();
        *points_dist2            = sil->points_dist2_get();
        *points_projected_coords = sil->points_projected_coords_get();
      }
    }
  }

  // Get internal cell_vtx ordering
  void
  Coupling::location_internal_cell_vtx_get
  (
   std::string    name,
   int         ***cell_vtx_idx,
   int         ***cell_vtx
  )
  {
    map <string, Field *>::iterator it  = _fields.find(name);
    if (it != _fields.end()) {
      CWP_Dof_location_t localFieldLocation = it->second->locationGet();
      CWP_Dof_location_t cplFieldLocation = it->second->linkedFieldLocationGet();

      std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2 = _spatial_interp_send.find(newKey);

      if (it2 != _spatial_interp_send.end()) {
        SpatialInterpLocation* sil = dynamic_cast <SpatialInterpLocation*> (it2->second);
        *cell_vtx_idx = sil->cell_vtx_idx_get();
        *cell_vtx     = sil->cell_vtx_get();
      }
    }
  }

  // Get local target elt volumes
  void
  Coupling::intersection_tgt_elt_volumes_get
  (
   std::string    name,
   double      ***tgt_elt_volumes
  )
  {
    map <string, Field *>::iterator it  = _fields.find(name);
    if (it != _fields.end()) {
      CWP_Dof_location_t localFieldLocation = it->second->locationGet();
      CWP_Dof_location_t cplFieldLocation = it->second->linkedFieldLocationGet();

      std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2 = _spatial_interp_recv.find(newKey);

      if (it2 != _spatial_interp_recv.end()) {
        SpatialInterpIntersection* sil = dynamic_cast <SpatialInterpIntersection*> (it2->second);
        *tgt_elt_volumes = sil->tgt_elt_volumes_get();
      }
    }
  }

  // Get closest src coord
  void
  Coupling::closest_point_src_coord_get
  (
   std::string    name,
   double      ***closest_src_coord
   )
  {
    map <string, Field *>::iterator it  = _fields.find(name);
    if (it != _fields.end()) {
      CWP_Dof_location_t localFieldLocation = it->second->locationGet();
      CWP_Dof_location_t cplFieldLocation = it->second->linkedFieldLocationGet();

      std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

      std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2 = _spatial_interp_recv.find(newKey);

      if (it2 != _spatial_interp_recv.end()) {
        SpatialInterpClosestPoint* sil = dynamic_cast <SpatialInterpClosestPoint*> (it2->second);
        *closest_src_coord = sil->closest_src_coord_get();
      }
    }
  }

  /**
   *
   * \brief Export mesh to Ensight format
   *
   */

  void 
  Coupling::exportMesh(Coupling &cpl)
  {

    if (cpl._writer != NULL) {

      /* First, create geometry and variables if necessary */
      if (cpl._n_step == 0) {

        // Geometry
        cpl._id_geom_writer =PDM_writer_geom_create_from_mesh_nodal (cpl._writer,
                                                                     "geom",
                                                                     cpl._mesh.getPdmNodalIndex());

        if (cpl._userTargetN != nullptr) {

          cpl._id_user_tgt_geom_writer = PDM_writer_geom_create(cpl._writer,
                                                                "user_tgt_geom",
                                                                cpl._nPart);

          int block_id = PDM_writer_geom_bloc_add(cpl._writer,
                                                  cpl._id_user_tgt_geom_writer,
                                                  PDM_WRITER_POINT,
                                                  PDM_OWNERSHIP_KEEP);

          for (int i_part = 0; i_part < cpl._nPart; i_part++) {

            int n_pts = cpl._userTargetN[i_part];
            CWP_g_num_t *tmp_point_gnum = (CWP_g_num_t *) cpl.userTargetGNumGet(i_part);

            CWP_g_num_t *point_gnum = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
            memcpy(point_gnum, tmp_point_gnum, sizeof(CWP_g_num_t) * n_pts);

            PDM_writer_geom_coord_set(cpl._writer,
                                      cpl._id_user_tgt_geom_writer,
                                      i_part,
                                      n_pts,
                                      cpl._userTargetCoord[i_part],
                                      tmp_point_gnum,
                                      PDM_OWNERSHIP_USER);

            int *point_connec = (int *) malloc(sizeof(int) * n_pts);
            for (int i = 0; i < n_pts; i++) {
              point_connec[i] = i+1;
            }

            PDM_writer_geom_bloc_std_set(cpl._writer,
                                         cpl._id_user_tgt_geom_writer,
                                         block_id,
                                         i_part,
                                         n_pts,
                                         point_connec,
                                         point_gnum);

          }

        } // end if there is a user target


        // Variables to show user targets partitionning
        PDM_writer_var_loc_t PDMfieldType = PDM_WRITER_VAR_ELEMENTS;

        PDM_writer_var_dim_t PDMfieldComp = PDM_WRITER_VAR_SCALAR;

        PDM_writer_status_t  st_dep_tps = PDM_WRITER_ON;

        if (cpl._displacement == CWP_DYNAMIC_MESH_STATIC) {
          st_dep_tps = PDM_WRITER_OFF;
        }

        cpl._id_field_partitioning_writer = PDM_writer_var_create(cpl._writer,
                                                                  st_dep_tps,
                                                                  PDMfieldComp,
                                                                  PDMfieldType,
                                                                  "partitioning");

        cpl._id_field_ranking_writer = PDM_writer_var_create(cpl._writer,
                                                             st_dep_tps,
                                                             PDMfieldComp,
                                                             PDMfieldType,
                                                             "ranking");

        if (cpl._userTargetN != nullptr) {
          cpl._id_user_tgt_field_partitioning_writer = PDM_writer_var_create(cpl._writer,
                                                                             st_dep_tps,
                                                                             PDMfieldComp,
                                                                             PDMfieldType,
                                                                             "user_tgt_partitioning");

          cpl._id_user_tgt_field_ranking_writer = PDM_writer_var_create(cpl._writer,
                                                                        st_dep_tps,
                                                                        PDMfieldComp,
                                                                        PDMfieldType,
                                                                        "user_tgt_ranking");
        } // end if there is a user target
      }

      /* Then begin a new time step if none is currently open */
      if (!PDM_writer_is_open_step (cpl._writer) && ((cpl._n_step % cpl._freq_writer) == 0)) {

        double current_time;

        cpl._localCodeProperties.ctrlParamGet("time", &current_time);

        PDM_writer_step_beg (cpl._writer, current_time);

      }


      /* Finally write geometry and variables */
      if (cpl._n_step == 0) {

        // Geometry
        PDM_writer_geom_write(cpl._writer, cpl._id_geom_writer);

        if (cpl._userTargetN != nullptr) {

          PDM_writer_geom_write(cpl._writer,
                              cpl._id_user_tgt_geom_writer);

        } // end if there is a user target

        // Variables to show user targets partitionning
        std::vector <double *> partitioning_field_data(cpl._mesh.getNPart());
        std::vector <double *> ranking_field_data(cpl._mesh.getNPart());
        std::vector <double *> user_tgt_partitioning_field_data(cpl._nPart);
        std::vector <double *> user_tgt_ranking_field_data(cpl._nPart);

        assert(cpl._mesh.getNPart() == cpl._nPart);

        int i_rank;
        MPI_Comm_rank (cpl._localCodeProperties.connectableCommGet(), &i_rank);

        int n_part = cpl._mesh.getNPart();
        int g_n_part = 0;

        MPI_Scan (&n_part, &g_n_part, 1, MPI_INT, MPI_SUM, cpl._localCodeProperties.connectableCommGet());

        g_n_part += -n_part;        

        for(int i_part= 0 ; i_part < cpl._mesh.getNPart(); i_part++){
          partitioning_field_data[i_part] = (double*) malloc(cpl._mesh.getPartNElts(i_part) * sizeof(double) );
          ranking_field_data     [i_part] = (double*) malloc(cpl._mesh.getPartNElts(i_part) * sizeof(double) );

          for(int i_elt = 0; i_elt < cpl._mesh.getPartNElts(i_part); i_elt++){
            partitioning_field_data[i_part][i_elt] = (double) (i_part + g_n_part);
            ranking_field_data[i_part][i_elt] = (double) i_rank;
          }

          if (cpl._userTargetN != nullptr) {
            user_tgt_partitioning_field_data[i_part] = (double*) malloc(cpl._userTargetN[i_part] * sizeof(double) );
            user_tgt_ranking_field_data     [i_part] = (double*) malloc(cpl._userTargetN[i_part] * sizeof(double) );
            for(int i = 0; i < cpl._userTargetN[i_part]; i++) {
              user_tgt_partitioning_field_data[i_part][i] = (double) (i_part + g_n_part);
              user_tgt_ranking_field_data     [i_part][i] = (double) i_rank;
            }
          } // end if there is a user target

          PDM_writer_var_set(cpl._writer, cpl._id_field_partitioning_writer, cpl._id_geom_writer, i_part, (double *) partitioning_field_data[i_part]);
          PDM_writer_var_set(cpl._writer, cpl._id_field_ranking_writer     , cpl._id_geom_writer, i_part, (double *) ranking_field_data     [i_part]);

          if (cpl._userTargetN != nullptr) {
            PDM_writer_var_set(cpl._writer, cpl._id_user_tgt_field_partitioning_writer, cpl._id_user_tgt_geom_writer, i_part, (double *) user_tgt_partitioning_field_data[i_part]);
            PDM_writer_var_set(cpl._writer, cpl._id_user_tgt_field_ranking_writer     , cpl._id_user_tgt_geom_writer, i_part, (double *) user_tgt_ranking_field_data     [i_part]);
          } // end if there is a user target
        }

        PDM_writer_var_write(cpl._writer, cpl._id_field_partitioning_writer);
        PDM_writer_var_write(cpl._writer, cpl._id_field_ranking_writer);

        PDM_writer_var_data_free(cpl._writer, cpl._id_field_partitioning_writer);
        PDM_writer_var_data_free(cpl._writer, cpl._id_field_ranking_writer);

        if (cpl._userTargetN != nullptr) {
          PDM_writer_var_write(cpl._writer, cpl._id_user_tgt_field_partitioning_writer);
          PDM_writer_var_write(cpl._writer, cpl._id_user_tgt_field_ranking_writer);

          PDM_writer_var_data_free(cpl._writer, cpl._id_user_tgt_field_partitioning_writer);
          PDM_writer_var_data_free(cpl._writer, cpl._id_user_tgt_field_ranking_writer);
        } // end if there is a user target

        for(int i_part= 0 ; i_part < cpl._mesh.getNPart(); i_part++){
          free (partitioning_field_data[i_part]);
          free (ranking_field_data     [i_part]);
          partitioning_field_data[i_part] = nullptr;
          ranking_field_data     [i_part] = nullptr;

          if (cpl._userTargetN != nullptr) {
            free (user_tgt_partitioning_field_data[i_part]);
            free (user_tgt_ranking_field_data     [i_part]);
            user_tgt_partitioning_field_data[i_part] = nullptr;
            user_tgt_ranking_field_data     [i_part] = nullptr;
          } // end if there is a user target
        }
      }

      else if (((cpl._n_step % cpl._freq_writer) == 0) && (cpl._displacement != CWP_DYNAMIC_MESH_STATIC)) {

        // Geometry
        if (cpl._displacement == CWP_DYNAMIC_MESH_VARIABLE) {
        

          PDM_writer_geom_set_from_mesh_nodal (cpl._writer, 
                                               cpl._id_geom_writer,
                                               cpl._mesh.getPdmNodalIndex());

          if (cpl._userTargetN != nullptr) {

            int block_id = PDM_writer_geom_bloc_add(cpl._writer,
                                                    cpl._id_user_tgt_geom_writer,
                                                    PDM_WRITER_POINT,
                                                    PDM_OWNERSHIP_KEEP);

            for (int i_part = 0; i_part < cpl._nPart; i_part++) {

              int n_pts = cpl._userTargetN[i_part];
              CWP_g_num_t *tmp_point_gnum = (CWP_g_num_t *) cpl.userTargetGNumGet(i_part);

              CWP_g_num_t *point_gnum = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
              memcpy(point_gnum, tmp_point_gnum, sizeof(CWP_g_num_t) * n_pts);

              PDM_writer_geom_coord_set(cpl._writer,
                                        cpl._id_user_tgt_geom_writer,
                                        i_part,
                                        n_pts,
                                        cpl._userTargetCoord[i_part],
                                        tmp_point_gnum,
                                        PDM_OWNERSHIP_USER);

              int *point_connec = (int *) malloc(sizeof(int) * n_pts);
              for (int i = 0; i < n_pts; i++) {
                point_connec[i] = i+1;
              }

              PDM_writer_geom_bloc_std_set(cpl._writer,
                                           cpl._id_user_tgt_geom_writer,
                                           block_id,
                                           i_part,
                                           n_pts,
                                           point_connec,
                                           point_gnum);

            }

          } // end if there is a user target

        }

        PDM_writer_geom_write(cpl._writer, cpl._id_geom_writer);

        if (cpl._userTargetN != nullptr) {

          PDM_writer_geom_write(cpl._writer, cpl._id_user_tgt_geom_writer);

        } // end if there is a user target

        // Variables to show user targets partitionning
        std::vector <double *> partitioning_field_data(cpl._mesh.getNPart());
        std::vector <double *> ranking_field_data(cpl._mesh.getNPart());
        std::vector <double *> user_tgt_partitioning_field_data(cpl._nPart);
        std::vector <double *> user_tgt_ranking_field_data(cpl._nPart);

        assert(cpl._mesh.getNPart() == cpl._nPart);

        int i_rank;
        MPI_Comm_rank (cpl._localCodeProperties.connectableCommGet(), &i_rank);

        int n_part = cpl._mesh.getNPart();
        int g_n_part = 0;

        MPI_Scan (&n_part, &g_n_part, 1, MPI_INT, MPI_SUM, cpl._localCodeProperties.connectableCommGet());

        g_n_part += -n_part;

        for(int i_part= 0 ; i_part < cpl._mesh.getNPart(); i_part++){
          partitioning_field_data[i_part] = (double*) malloc(cpl._mesh.getPartNElts(i_part) * sizeof(double) );
          ranking_field_data     [i_part] = (double*) malloc(cpl._mesh.getPartNElts(i_part) * sizeof(double) );

          for(int i_elt = 0; i_elt < cpl._mesh.getPartNElts(i_part); i_elt++){
            partitioning_field_data[i_part][i_elt] = (double) (i_part + g_n_part);
            ranking_field_data     [i_part][i_elt] = (double) i_rank;
          }

          if (cpl._userTargetN != nullptr) {
            user_tgt_partitioning_field_data[i_part] = (double*) malloc(cpl._userTargetN[i_part] * sizeof(double) );
            user_tgt_ranking_field_data     [i_part] = (double*) malloc(cpl._userTargetN[i_part] * sizeof(double) );
            for(int i = 0; i < cpl._userTargetN[i_part]; i++) {
              user_tgt_partitioning_field_data[i_part][i] = (double) (i_part + g_n_part);
              user_tgt_ranking_field_data     [i_part][i] = (double) i_rank;
            }
          } // end if there is a user target

          PDM_writer_var_set(cpl._writer, cpl._id_field_partitioning_writer, cpl._id_geom_writer, i_part, (double *) partitioning_field_data[i_part]);
          PDM_writer_var_set(cpl._writer, cpl._id_field_ranking_writer     , cpl._id_geom_writer, i_part, (double *) ranking_field_data     [i_part]);

          if (cpl._userTargetN != nullptr) {
            PDM_writer_var_set(cpl._writer, cpl._id_user_tgt_field_partitioning_writer, cpl._id_user_tgt_geom_writer, i_part, (double *) user_tgt_partitioning_field_data[i_part]);
            PDM_writer_var_set(cpl._writer, cpl._id_user_tgt_field_ranking_writer     , cpl._id_user_tgt_geom_writer, i_part, (double *) user_tgt_ranking_field_data     [i_part]);
          } // end if there is a user target
        }

        PDM_writer_var_write(cpl._writer, cpl._id_field_partitioning_writer);
        PDM_writer_var_write(cpl._writer, cpl._id_field_ranking_writer);

        PDM_writer_var_data_free(cpl._writer, cpl._id_field_partitioning_writer);
        PDM_writer_var_data_free(cpl._writer, cpl._id_field_ranking_writer);

        if (cpl._userTargetN != nullptr) {
          PDM_writer_var_write(cpl._writer, cpl._id_user_tgt_field_partitioning_writer);
          PDM_writer_var_write(cpl._writer, cpl._id_user_tgt_field_ranking_writer);

          PDM_writer_var_data_free(cpl._writer, cpl._id_user_tgt_field_partitioning_writer);
          PDM_writer_var_data_free(cpl._writer, cpl._id_user_tgt_field_ranking_writer);
        } // end if there is a user target

        for(int i_part= 0 ; i_part < cpl._mesh.getNPart(); i_part++){
          free (partitioning_field_data[i_part]);
          free (ranking_field_data     [i_part]);
          partitioning_field_data[i_part] = nullptr;
          ranking_field_data     [i_part] = nullptr;

          if (cpl._userTargetN != nullptr) {
            free (user_tgt_partitioning_field_data[i_part]);
            free (user_tgt_ranking_field_data     [i_part]);
            user_tgt_partitioning_field_data[i_part] = nullptr;
            user_tgt_ranking_field_data     [i_part] = nullptr;
          } // end if there is a user target
        }

      }
    }
  }

  /**
   * \brief Computation spatial interpolation weights
   *
   * This function compute spatial interpolation weights
   *
   */

  void
  Coupling::spatialInterpWeightsCompute ()
  {

    /////////////////////////////////////////////////////////////////////////////
    //                                                                         //
    // Export mesh and associted fields                                        //
    //                                                                         //
    /////////////////////////////////////////////////////////////////////////////



    if (!_coupledCodeProperties.localCodeIs()) {

      exportMesh(*this);

    }

    else {
      int codeID    = localCodePropertiesGet()->idGet();
      int cplCodeID = coupledCodePropertiesGet()->idGet();

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        if (codeID < cplCodeID) {

          exportMesh(*this);
          exportMesh(cpl_cpl);   

        }
      }
    }

    /////////////////////////////////////////////////////////////////////////////
    //                                                                         //
    // Exchange fields properties to obtain spatial intepolation to build      //
    //                                                                         //
    /////////////////////////////////////////////////////////////////////////////

    // - Store Data to send

    if (!_coupledCodeProperties.localCodeIs()) {

      int codeID    = localCodePropertiesGet()->idGet();
      int cplCodeID = coupledCodePropertiesGet()->idGet();


      if (_n_step == 0) {

        std::string localFieldsName="";
        vector<int> localFieldsNameIdx;

        int localNbField = (int) _fields.size();

        localFieldsNameIdx.reserve(localNbField + 1);
        localFieldsNameIdx.push_back(0);

        std::vector<CWP_Field_exch_t> localFieldsExch;
        localFieldsExch.reserve(localNbField);

        std::vector<CWP_Dof_location_t> localFieldLocationV;
        localFieldLocationV.reserve(localNbField);

        std::map <std::string, cwipi::Field *>::iterator it = _fields.begin();

        while(it != _fields.end()){
          cwipi::Field* field = it->second;

          localFieldsName += it->first;
          localFieldsNameIdx.push_back((int) (localFieldsNameIdx[localFieldsNameIdx.size()-1]+it->first.size()));
          localFieldLocationV.push_back(field->locationGet());
          localFieldsExch.push_back(field->exchangeTypeGet());

          it++;
        }

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   1,
                                                                   (void *) &_nPart,
                                                                   -1,
                                                                   NULL,
                                                                   1,
                                                                   (void *) &_cplNPart,
                                                                   -1,
                                                                   NULL);
        // - Exchange number of fields

        int nSendData   = 1;
        int nRecvData   = 1;
        int cplNbField = 0;

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   nSendData,
                                                                   (void *) &localNbField,
                                                                   -1,
                                                                   NULL,
                                                                   nRecvData,
                                                                   (void *) &cplNbField,
                                                                   -1,
                                                                   NULL);

        // - Allocate memory to receive data

        vector<int               > cplFieldNameIdx (cplNbField + 1, 0);
        vector<CWP_Field_exch_t  > cplFieldExch (cplNbField);
        vector<CWP_Dof_location_t> cplFieldLocationV (cplNbField);
        string                     cplFieldName;

        // - Transfer memory to receive data

        nSendData   = localNbField + 1;
        nRecvData   = cplNbField + 1;
        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   nSendData,
                                                                   (void *) &(localFieldsNameIdx[0]),
                                                                   -1,
                                                                   NULL,
                                                                   nRecvData,
                                                                   (void *) &(cplFieldNameIdx[0]),
                                                                   -1,
                                                                   NULL);

        cplFieldName.resize(cplFieldNameIdx[cplNbField]);
        nSendData   = localFieldsNameIdx[localNbField];
        nRecvData   = cplFieldNameIdx[cplNbField];

        localFieldsNameIdx.clear();


        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(char),
                                                                   nSendData,
                                                                   (void *) localFieldsName.c_str(),
                                                                   -1,
                                                                   NULL,
                                                                   nRecvData,
                                                                   (void *) cplFieldName.c_str(),
                                                                   -1,
                                                                   NULL);

        localFieldsName.clear();

        nSendData   = localNbField;
        nRecvData   = cplNbField;
        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Field_exch_t),
                                                                   nSendData,
                                                                   (void *) &(localFieldsExch[0]),
                                                                   -1,
                                                                   NULL,
                                                                   nRecvData,
                                                                   (void *) &(cplFieldExch[0]),
                                                                   -1,
                                                                   NULL);
        localFieldsExch.clear();

        nSendData   = localNbField;
        nRecvData   = cplNbField;
        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                   nSendData,
                                                                   (void *) &(localFieldLocationV[0]),
                                                                   -1,
                                                                   NULL,
                                                                   nRecvData,
                                                                   (void *) &(cplFieldLocationV[0]),
                                                                   -1,
                                                                   NULL);
        localFieldLocationV.clear();

        // Create spatial interpolation objects

        it = _fields.begin();
        while(it != _fields.end()) {
          cwipi::Field* field = it->second;
          string localFieldName = it->first;
          CWP_Field_exch_t   localFieldExch     = field->exchangeTypeGet();
          CWP_Dof_location_t localFieldLocation = field->locationGet();

          for (int j = 0; j < cplNbField; j++) {
            string _cplFieldName = cplFieldName.substr(cplFieldNameIdx[j], cplFieldNameIdx[j+1]-cplFieldNameIdx[j]);

            if (_cplFieldName == localFieldName) {
              field->linkedFieldLocationSet(cplFieldLocationV[j]);

              if (  (localFieldExch == CWP_FIELD_EXCH_SENDRECV && cplFieldExch[j] == CWP_FIELD_EXCH_SENDRECV)
                  ||(localFieldExch == CWP_FIELD_EXCH_SEND     && cplFieldExch[j] == CWP_FIELD_EXCH_RECV)
                  ||(localFieldExch == CWP_FIELD_EXCH_RECV     && cplFieldExch[j] == CWP_FIELD_EXCH_SEND)) {

                if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_SEND) {
                  std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocationV[j]);
                  if (_spatial_interp_send.find(newKey) == _spatial_interp_send.end()) {

                    _spatial_interp_send.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));

                    _spatial_interp_send[newKey]->init(this, localFieldLocation, cplFieldLocationV[j], SPATIAL_INTERP_EXCH_SEND);
                  }
                }

                if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_RECV) {
                  std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocationV[j]);
                  if (_spatial_interp_recv.find(newKey) == _spatial_interp_recv.end()) {

                    _spatial_interp_recv.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));

                    _spatial_interp_recv[newKey]->init(this, localFieldLocation, cplFieldLocationV[j], SPATIAL_INTERP_EXCH_RECV);

                  }
                }
              }
            }
          }
          it++;
        }

        int sis_s = (int) _spatial_interp_send.size();

        vector<CWP_Dof_location_t> sis_loc;
        sis_loc.reserve(2*sis_s);

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          sis_loc.push_back((sis_it->first).first);
          sis_loc.push_back((sis_it->first).second);
          sis_it++;
        }

        int sis_r;

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   1,
                                                                   (void *) &sis_s,
                                                                   -1,
                                                                   NULL,
                                                                   1,
                                                                   (void *) &sis_r,
                                                                   -1,
                                                                   NULL);

        int sir_s = (int) _spatial_interp_recv.size();

        vector<CWP_Dof_location_t> sir_loc;
        sir_loc.reserve(2*sir_s);

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sir_it = _spatial_interp_recv.begin();
        while(sir_it != _spatial_interp_recv.end()) {
          sir_loc.push_back((sir_it->first).first);
          sir_loc.push_back((sir_it->first).second);
          sir_it++;
        }

        int sir_r;

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                   1,
                                                                   (void *) &sir_s,
                                                                   -1,
                                                                   NULL,
                                                                   1,
                                                                   (void *) &sir_r,
                                                                   -1,
                                                                   NULL);


        assert(sir_r == sis_s);

        _sis_loc_r.resize(2*sir_s);

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                   2*sis_s,
                                                                   (void *) &(sis_loc[0]),
                                                                   -1,
                                                                   NULL,
                                                                   2*sis_r,
                                                                   (void *) &(_sis_loc_r[0]),
                                                                   -1,
                                                                   NULL);

        vector<CWP_Dof_location_t> sir_loc_r;
        sir_loc_r.resize(2*sis_s);

        _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                   2*sir_s,
                                                                   (void *) &(sir_loc[0]),
                                                                   -1,
                                                                   NULL,
                                                                   2*sir_r,
                                                                   (void *) &(sir_loc_r[0]),
                                                                   -1,
                                                                   NULL);
      }


      int clear_data = (_n_step > 0);//(_displacement == CWP_DYNAMIC_MESH_VARIABLE && _n_step > 0);

      if (codeID < cplCodeID) {

        // spatial_interp send

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          if (clear_data) {
            sis_it->second->clear();
            sis_it->second->init(this,
                                 sis_it->first.first,
                                 sis_it->first.second,
                                 SPATIAL_INTERP_EXCH_SEND);
          }
          sis_it->second->weightsCompute();
          sis_it++;
        }

        // spatial_interp recv

        for (int i = 0; i < (int) (_sis_loc_r.size()/2); i++) {
          if (clear_data) {
            _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])]->clear();
            _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])]->init(this,
                                                                                      _sis_loc_r[2*i+1],
                                                                                      _sis_loc_r[2*i  ],
                                                                                      SPATIAL_INTERP_EXCH_RECV);
          }
          _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])]->weightsCompute();
        }

      }

      else {

        // spatial_interp recv

        for (int i = 0; i < (int) (_sis_loc_r.size()/2); i++) {
          if (clear_data) {
            _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])]->clear();
            _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])]->init(this,
                                                                                      _sis_loc_r[2*i+1],
                                                                                      _sis_loc_r[2*i  ],
                                                                                      SPATIAL_INTERP_EXCH_RECV);
          }
          _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])]->weightsCompute();
        }

        // spatial_interp send

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
        while(sis_it != _spatial_interp_send.end()) {
          if (clear_data) {
            sis_it->second->clear();
            sis_it->second->init(this,
                                 sis_it->first.first,
                                 sis_it->first.second,
                                 SPATIAL_INTERP_EXCH_SEND);
          }
          sis_it->second->weightsCompute();
          sis_it++;
        }

      }

    }

    // if an instance of each code on the same rank. All work is made in the call from the smallest code id

    else {

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        int codeID    = localCodePropertiesGet()->idGet();
        int cplCodeID = coupledCodePropertiesGet()->idGet();

        if (_n_step == 0) {

          std::string localFieldsName="";
          vector<int> localFieldsNameIdx;

          int localNbField = (int) _fields.size();

          localFieldsNameIdx.reserve(localNbField + 1);
          localFieldsNameIdx.push_back(0);

          std::vector<CWP_Field_exch_t> localFieldsExch;
          localFieldsExch.reserve(localNbField);

          std::vector<CWP_Dof_location_t> localFieldLocationV;
          localFieldLocationV.reserve(localNbField);

          std::map <std::string, cwipi::Field *>::iterator it = _fields.begin();

          while(it != _fields.end()){
            cwipi::Field* field = it->second;

            localFieldsName += it->first;
            localFieldsNameIdx.push_back((int) (localFieldsNameIdx[localFieldsNameIdx.size()-1]+it->first.size()));
            localFieldLocationV.push_back(field->locationGet());
            localFieldsExch.push_back(field->exchangeTypeGet());

            it++;
          }

          std::string cpl_localFieldsName="";
          vector<int> cpl_localFieldsNameIdx;

          int cpl_localNbField = (int) cpl_cpl._fields.size();

          cpl_localFieldsNameIdx.reserve(cpl_localNbField + 1);
          cpl_localFieldsNameIdx.push_back(0);

          std::vector<CWP_Field_exch_t> cpl_localFieldsExch;
          cpl_localFieldsExch.reserve(cpl_localNbField);

          std::vector<CWP_Dof_location_t> cpl_localFieldLocationV;
          cpl_localFieldLocationV.reserve(cpl_localNbField);

          it = cpl_cpl._fields.begin();

          while(it != cpl_cpl._fields.end()){
            cwipi::Field* field = it->second;

            cpl_localFieldsName += it->first;
            cpl_localFieldsNameIdx.push_back((int) (cpl_localFieldsNameIdx[cpl_localFieldsNameIdx.size()-1]+it->first.size()));
            cpl_localFieldLocationV.push_back(field->locationGet());
            cpl_localFieldsExch.push_back(field->exchangeTypeGet());

            it++;
          }

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                     1,
                                                                     (void *) &_nPart,
                                                                     1,
                                                                     (void *) &cpl_cpl._nPart,
                                                                     1,
                                                                     (void *) &_cplNPart,
                                                                     1,
                                                                     (void *) &cpl_cpl._cplNPart);

          // - Exchange number of fields

          int nSendData       = 1;
          int nRecvData       = 1;
          int cplNbField      = 0;

          int cpl_nSendData       = 1;
          int cpl_nRecvData       = 1;
          int cpl_cplNbField      = 0;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                     nSendData,
                                                                     (void *) &localNbField,
                                                                     cpl_nSendData,
                                                                     (void *) &cpl_localNbField,
                                                                     nRecvData,
                                                                     (void *) &cplNbField,
                                                                     cpl_nRecvData,
                                                                     (void *) &cpl_cplNbField);


          // - Allocate memory to receive data

          vector<int               > cplFieldsNameIdx (cplNbField + 1, 0);
          vector<CWP_Field_exch_t  > cplFieldsExch (cplNbField);
          vector<CWP_Dof_location_t> cplFieldsLocationV (cplNbField);
          string                     cplFieldsName;

          vector<int               > cpl_cplFieldsNameIdx (cpl_cplNbField + 1, 0);
          vector<CWP_Field_exch_t  > cpl_cplFieldsExch (cpl_cplNbField);
          vector<CWP_Dof_location_t> cpl_cplFieldsLocationV (cpl_cplNbField);
          string                     cpl_cplFieldsName;

          // - Transfer memory to receive data

          nSendData   = localNbField + 1;
          nRecvData   = cplNbField + 1;

          cpl_nSendData   = cpl_localNbField + 1;
          cpl_nRecvData   = cpl_cplNbField + 1;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                     nSendData,
                                                                     (void *) &(localFieldsNameIdx[0]),
                                                                     cpl_nSendData,
                                                                     (void *) &(cpl_localFieldsNameIdx[0]),
                                                                     nRecvData,
                                                                     (void *) &(cplFieldsNameIdx[0]),
                                                                     cpl_nRecvData,
                                                                     (void *) &(cpl_cplFieldsNameIdx[0]));

          cplFieldsName.resize(cplFieldsNameIdx[cplNbField]);
          nSendData   = localFieldsNameIdx[localNbField];
          nRecvData   = cplFieldsNameIdx[cplNbField];

          cpl_cplFieldsName.resize(cpl_cplFieldsNameIdx[cpl_cplNbField]);
          cpl_nSendData   = cpl_localFieldsNameIdx[cpl_localNbField];
          cpl_nRecvData   = cpl_cplFieldsNameIdx[cpl_cplNbField];

          localFieldsNameIdx.clear();
          cpl_localFieldsNameIdx.clear();


          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(char),
                                                                     nSendData,
                                                                     (void *) localFieldsName.c_str(),
                                                                     cpl_nSendData,
                                                                     (void *) cpl_localFieldsName.c_str(),
                                                                     nRecvData,
                                                                     (void *) cplFieldsName.c_str(),
                                                                     cpl_nRecvData,
                                                                     (void *) cpl_cplFieldsName.c_str());

          localFieldsName.clear();
          cpl_localFieldsName.clear();

          nSendData   = localNbField;
          nRecvData   = cplNbField;

          cpl_nSendData   = cpl_localNbField;
          cpl_nRecvData   = cpl_cplNbField;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Field_exch_t),
                                                                     nSendData,
                                                                     (void *) &(localFieldsExch[0]),
                                                                     cpl_nSendData,
                                                                     (void *) &(cpl_localFieldsExch[0]),
                                                                     nRecvData,
                                                                     (void *) &(cplFieldsExch[0]),
                                                                     cpl_nRecvData,
                                                                     (void *) &(cpl_cplFieldsExch[0]));

          localFieldsExch.clear();
          cpl_localFieldsExch.clear();

          nSendData   = localNbField;
          nRecvData   = cplNbField;

          cpl_nSendData   = cpl_localNbField;
          cpl_nRecvData   = cpl_cplNbField;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                     nSendData,
                                                                     (void *) &(localFieldLocationV[0]),
                                                                     cpl_nSendData,
                                                                     (void *) &(cpl_localFieldLocationV[0]),
                                                                     nRecvData,
                                                                     (void *) &(cplFieldsLocationV[0]),
                                                                     cpl_nRecvData,
                                                                     (void *) &(cpl_cplFieldsLocationV[0]));

          localFieldLocationV.clear();
          cpl_localFieldLocationV.clear();

          // Create spatial interpolation objects

          it = _fields.begin();
          while(it != _fields.end()) {
            cwipi::Field* field = it->second;
            string localFieldName = it->first;
            CWP_Field_exch_t   localFieldExch     = field->exchangeTypeGet();
            CWP_Dof_location_t localFieldLocation = field->locationGet();

            for (int j = 0; j < cplNbField; j++) {

              string cplFieldName = cplFieldsName.substr( cplFieldsNameIdx[j], cplFieldsNameIdx[j+1]-cplFieldsNameIdx[j] );

              if (cplFieldName == localFieldName) {
                field->linkedFieldLocationSet(cplFieldsLocationV[j]);
                if (cpl_cpl._fields.find(cplFieldName) == cpl_cpl._fields.end()) {
                  PDM_error(__FILE__, __LINE__, 0, "'%s' Field not found\n",
                            cplFieldName.c_str());
                }

                cpl_cpl._fields[cplFieldName]->linkedFieldLocationSet(localFieldLocation);

                if (  (localFieldExch == CWP_FIELD_EXCH_SENDRECV && cplFieldsExch[j] == CWP_FIELD_EXCH_SENDRECV)
                    ||(localFieldExch == CWP_FIELD_EXCH_SEND     && cplFieldsExch[j] == CWP_FIELD_EXCH_RECV)
                    ||(localFieldExch == CWP_FIELD_EXCH_RECV     && cplFieldsExch[j] == CWP_FIELD_EXCH_SEND)) {

                  if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_SEND) {
                    std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldsLocationV[j]);
                    std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldsLocationV[j], localFieldLocation);

                    if (_spatial_interp_send.find(newKey) == _spatial_interp_send.end()) {

                      _spatial_interp_send.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));
                      cpl_cpl._spatial_interp_recv.insert(make_pair(cpl_newKey, FG::getInstance().CreateObject(cpl_cpl._spatialInterpAlgo)));
              
                      _spatial_interp_send[newKey]->init(this, localFieldLocation, cplFieldsLocationV[j], SPATIAL_INTERP_EXCH_SEND);
                      cpl_cpl._spatial_interp_recv[cpl_newKey]->init(&cpl_cpl, cplFieldsLocationV[j], localFieldLocation, SPATIAL_INTERP_EXCH_RECV);        

                    }
                  }

                  if (localFieldExch == CWP_FIELD_EXCH_SENDRECV || localFieldExch == CWP_FIELD_EXCH_RECV) {

                    std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldsLocationV[j]);
                    std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldsLocationV[j], localFieldLocation);

                    if (_spatial_interp_recv.find(newKey) == _spatial_interp_recv.end()) {

                      _spatial_interp_recv.insert(make_pair(newKey, FG::getInstance().CreateObject(_spatialInterpAlgo)));
                      cpl_cpl._spatial_interp_send.insert(make_pair(cpl_newKey, FG::getInstance().CreateObject(cpl_cpl._spatialInterpAlgo)));

                      _spatial_interp_recv[newKey]->init(this, localFieldLocation, cplFieldsLocationV[j], SPATIAL_INTERP_EXCH_RECV);
                      cpl_cpl._spatial_interp_send[cpl_newKey]->init(&cpl_cpl, cplFieldsLocationV[j], localFieldLocation, SPATIAL_INTERP_EXCH_SEND);

                    }
                  }
                }
              }
            }
            it++;
          }

          int sis_s = (int) _spatial_interp_send.size();

          vector<CWP_Dof_location_t> sis_loc;
          sis_loc.reserve(2*sis_s);

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
          while(sis_it != _spatial_interp_send.end()) {
            sis_loc.push_back((sis_it->first).first);
            sis_loc.push_back((sis_it->first).second);
            sis_it++;
          }

          int cpl_sis_s = (int) cpl_cpl._spatial_interp_send.size();

          vector<CWP_Dof_location_t> cpl_sis_loc;
          cpl_sis_loc.reserve(2*cpl_sis_s);

          sis_it = cpl_cpl._spatial_interp_send.begin();
          while(sis_it != cpl_cpl._spatial_interp_send.end()) {
            cpl_sis_loc.push_back((sis_it->first).first);
            cpl_sis_loc.push_back((sis_it->first).second);
            sis_it++;
          }

          int sis_r;
          int cpl_sis_r;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                     1,
                                                                     (void *) &sis_s,
                                                                     1,
                                                                     (void *) &cpl_sis_s,
                                                                     1,
                                                                     (void *) &sis_r,
                                                                     1,
                                                                     (void *) &cpl_sis_r);

          int sir_s = (int) _spatial_interp_recv.size();
    
          vector<CWP_Dof_location_t> sir_loc;
          sir_loc.reserve(2*sis_r);
    
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sir_it = _spatial_interp_recv.begin();
          while(sir_it != _spatial_interp_recv.end()) {
            sir_loc.push_back((sir_it->first).first);
            sir_loc.push_back((sir_it->first).second);
            sir_it++;
          }
    
          int cpl_sir_s = (int) cpl_cpl._spatial_interp_recv.size();
    
          vector<CWP_Dof_location_t> cpl_sir_loc;
          cpl_sir_loc.reserve(2*cpl_sir_s);
    
          sir_it = cpl_cpl._spatial_interp_recv.begin();
          while(sir_it != cpl_cpl._spatial_interp_recv.end()) {
            cpl_sir_loc.push_back((sir_it->first).first);
            cpl_sir_loc.push_back((sir_it->first).second);
            sir_it++;
          }
    
          int sir_r;
          int cpl_sir_r;

          // sir_s = 5;
          // cpl_sir_s = 7;

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(int),
                                                                     1,
                                                                     (void *) &sir_s,
                                                                     1,
                                                                     (void *) &cpl_sir_s,
                                                                     1,
                                                                     (void *) &sir_r,
                                                                     1,
                                                                     (void *) &cpl_sir_r);


          _sis_loc_r.resize(2*sir_s);

          _cpl_sis_loc_r.resize(2*cpl_sir_s);
    
          assert(sir_r     == sis_s);
          assert(sis_r     == sir_s);
          assert(cpl_sir_r == cpl_sis_s);
          assert(cpl_sis_r == cpl_sir_s);

          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                     2*sis_s,
                                                                     (void *) &(sis_loc[0]),
                                                                     2*cpl_sis_s,
                                                                     (void *) &(cpl_sis_loc[0]),
                                                                     2*sis_r,
                                                                     (void *) &(_sis_loc_r[0]),
                                                                     2*cpl_sis_r,
                                                                     (void *) &(_cpl_sis_loc_r[0]));
    
          vector<CWP_Dof_location_t> sir_loc_r;
          sir_loc_r.resize(2*sis_s);

          vector<CWP_Dof_location_t> cpl_sir_loc_r;
          cpl_sir_loc_r.resize(2*cpl_sis_s);
    
          _communication.iexchGlobalDataBetweenCodesThroughUnionCom (sizeof(CWP_Dof_location_t),
                                                                     2*sir_s,
                                                                     (void *) &(sir_loc[0]),
                                                                     2*cpl_sir_s,
                                                                     (void *) &(cpl_sir_loc[0]),
                                                                     2*sir_r,
                                                                     (void *) &(sir_loc_r[0]),
                                                                     2*cpl_sir_r,
                                                                     (void *) &(cpl_sir_loc_r[0]));
        }

        int clear_data = (_n_step > 0);//(_displacement == CWP_DYNAMIC_MESH_VARIABLE && _n_step > 0);

        if (codeID < cplCodeID) {
  
          int i = 0;
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = _spatial_interp_send.begin();
          while(sis_it != _spatial_interp_send.end()) {
            SpatialInterp *cpl_spatial_interp = cpl_cpl._spatial_interp_recv[make_pair(_cpl_sis_loc_r[2*i+1], _cpl_sis_loc_r[2*i])];
            if (clear_data) {
              sis_it->second->clear();
              sis_it->second->init(this,
                                   sis_it->first.first,
                                   sis_it->first.second,
                                   SPATIAL_INTERP_EXCH_SEND);
              cpl_spatial_interp->clear();
              cpl_spatial_interp->init(&cpl_cpl,
                                       _cpl_sis_loc_r[2*i+1],
                                       _cpl_sis_loc_r[2*i  ],
                                       SPATIAL_INTERP_EXCH_RECV);
            }
            sis_it->second->weightsCompute();
            cpl_spatial_interp->weightsCompute();
            sis_it++;
            i++;
          }
  
          i = 0;
          sis_it = cpl_cpl._spatial_interp_send.begin();
          while(sis_it != cpl_cpl._spatial_interp_send.end()) {
            SpatialInterp *cpl_spatial_interp = _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])];
            if (clear_data) {
              cpl_spatial_interp->clear();
              cpl_spatial_interp->init(&cpl_cpl,
                                       _sis_loc_r[2*i+1],
                                       _sis_loc_r[2*i  ],
                                       SPATIAL_INTERP_EXCH_RECV);
              sis_it->second->clear();
              sis_it->second->init(this,
                                   sis_it->first.first,
                                   sis_it->first.second,
                                   SPATIAL_INTERP_EXCH_SEND);
            }
            cpl_spatial_interp->weightsCompute();
            sis_it->second->weightsCompute();
            sis_it++;
            i++;
          }
        }

        else {
  
          int i = 0;
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator sis_it = cpl_cpl._spatial_interp_send.begin();
          while(sis_it != cpl_cpl._spatial_interp_send.end()) {
            SpatialInterp *cpl_spatial_interp = _spatial_interp_recv[make_pair(_sis_loc_r[2*i+1], _sis_loc_r[2*i])];
            if (clear_data) {
              cpl_spatial_interp->clear();
              cpl_spatial_interp->init(&cpl_cpl,
                                       _sis_loc_r[2*i+1],
                                       _sis_loc_r[2*i  ],
                                       SPATIAL_INTERP_EXCH_RECV);
              sis_it->second->clear();
              sis_it->second->init(this,
                                   sis_it->first.first,
                                   sis_it->first.second,
                                   SPATIAL_INTERP_EXCH_SEND);
            }
            cpl_spatial_interp->weightsCompute();
            sis_it->second->weightsCompute();
            sis_it++;
            i++;
          }
  
          i = 0;
          sis_it = _spatial_interp_send.begin();
          while(sis_it != _spatial_interp_send.end()) {
            SpatialInterp *cpl_spatial_interp = cpl_cpl._spatial_interp_recv[make_pair(_cpl_sis_loc_r[2*i+1], _cpl_sis_loc_r[2*i])];
            if (clear_data) {
              sis_it->second->clear();
              sis_it->second->init(this,
                                   sis_it->first.first,
                                   sis_it->first.second,
                                   SPATIAL_INTERP_EXCH_SEND);
              cpl_spatial_interp->clear();
              cpl_spatial_interp->init(&cpl_cpl,
                                       _cpl_sis_loc_r[2*i+1],
                                       _cpl_sis_loc_r[2*i  ],
                                       SPATIAL_INTERP_EXCH_RECV);
            }
            sis_it->second->weightsCompute();
            cpl_spatial_interp->weightsCompute();
            sis_it++;
            i++;
          }
        }
      }
    }

  }


  /*----------------------------------------------------------------------------*
   * Methods about visualization                                                *
   *----------------------------------------------------------------------------*/

  /**
   * \brief Enable visualization output
   *
   * This function enable visualization output.
   *
   * \param [in]  freq             Output frequency
   * \param [in]  format           Output format to visualize exchanged fieldsDouble
   *                               on the coupled mesh. Choice between :
   *                               - "EnSight Gold"
   *                               - "MED_ficher"
   *                               - "CGNS"
   *                               .
   * \param [in]  format_option   Output options "opt1, opt2, ..." :
   *                         - text               output text files
   *                         - binary             output binary files (default)
   *                         - big_endian         force binary files
   *                                              to big-endian
   *                         - discard_polygons   do not output polygons
   *                                              or related values
   *                         - discard_polyhedra  do not output polyhedra
   *                                              or related values
   *                         - divide_polygons    tesselate polygons
   *                                              with triangles
   *                         - divide_polyhedra   tesselate polyhedra
   *                                              with tetrahedra and pyramids
   *                                              (adding a vertex near
   *                                               each polyhedron's center)
   *                         .
   *
   */

  void
  Coupling::visuSet (
    const int               freq,
    const CWP_Visu_format_t format,
    const char             *format_option
  )
  {
    CWP_UNUSED(format);

    if (_is_first_field_created) {
      PDM_error(__FILE__, __LINE__, 0,
                "Error : CWP_Visu_set has to be called before CWP_Field_create.\n");
      abort();
    }

    if (_is_mesh_finalized) {
      PDM_error(__FILE__, __LINE__, 0,
                "Error : CWP_Visu_set has to be called before CWP_Mesh_interf_finalize.\n");
      abort();
    }

    string CodeName = _localCodeProperties.nameGet();
    string cplCodeName = _coupledCodeProperties.nameGet();
    string cplId = IdGet();

    string visuDir = "cwipi";
    string output_dir = visuDir+"/"+cplId+"_"+CodeName+"_"+cplCodeName;
    string output_dir_tmp = visuDir+"_writer/"+cplId+"_"+CodeName+"_"+cplCodeName;

    int rank;

    MPI_Comm_rank(_communication.unionCommGet(),&rank);

    // if ((commTypeGet() == CWP_COMM_PAR_WITH_PART ||
    //     (commTypeGet() == CWP_COMM_PAR_WITHOUT_PART && rank == _communication.unionCommLocCodeRootRanksGet())) && freq > 0) {
    if (freq > 0) {

      _freq_writer = freq;

      PDM_MPI_Comm pdmComm = PDM_MPI_mpi_2_pdm_mpi_comm(const_cast<MPI_Comm*>(&_localCodeProperties.connectableCommGet()));

      PDM_writer_topology_t pdm_topology = PDM_WRITER_TOPO_CST;

      if     ( _displacement == CWP_DYNAMIC_MESH_STATIC     ) {
        pdm_topology  = PDM_WRITER_TOPO_CST;
      }
      
      else if( _displacement == CWP_DYNAMIC_MESH_DEFORMABLE ) {
        pdm_topology  = PDM_WRITER_TOPO_DEFORMABLE;
      }

      else { 
        pdm_topology  = PDM_WRITER_TOPO_VARIABLE;
      }

      PDM_writer_fmt_fic_t fmt_fic      = PDM_WRITER_FMT_BIN;
      const char* fmt                   = "Ensight";
      PDM_writer_status_t st_reprise    = PDM_WRITER_OFF;
      const char *options_comp          = "";
      //Proportion of working node for file acess
      int working_node = 1;
  
      PDM_io_kind_t acess_type =   PDM_IO_KIND_MPI_SIMPLE;

      std::string str_options = format_option;
      std::string delimiter   = ",";

      std::string chars = "\t\n\v\f\r ";

      size_t pos = 0;
      std::string option;
      do
      {
        pos = str_options.find(delimiter);
        option = str_options.substr(0, pos);
        option.erase(0, option.find_first_not_of(chars));
        option.erase(option.find_last_not_of(chars)+1,option.length());
        str_options.erase(0, pos + delimiter.length());

        if (option == "text") {
          fmt_fic = PDM_WRITER_FMT_ASCII;
        }
        else if (option == "binary") {
          fmt_fic = PDM_WRITER_FMT_BIN;
        }
        else if (option != "") {
          PDM_error(__FILE__, __LINE__, 0,
                    "Not a valid visualization option.\n");
        } 

      }
      while (pos != std::string::npos);

      _writer = PDM_writer_create(fmt,
                                  fmt_fic,
                                  pdm_topology,
                                  st_reprise,
                                  (char *) output_dir_tmp.c_str(),
                                  (char *) string("chr").c_str(),
                                  pdmComm,
                                  acess_type,
                                  working_node,
                                  options_comp);
    }
  }

  /**
   *
   * \brief End visualization output
   *
   */

  void
  Coupling::visuEnd ()
  {
    if (_writer != NULL) {
      PDM_writer_free (_writer);
      _writer = nullptr;
    }
  }


  /**
   *
   * \brief MPI Barrier on the coupling communicator.
   *
   */

  void
  Coupling::barrier()
  {
    MPI_Comm unionComm = _communication.unionCommGet();

    if (!_coupledCodeProperties.localCodeIs()) {
      MPI_Barrier(unionComm);
      return;
    }

    if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {
      MPI_Barrier(unionComm);
      return;
    }
  }

  /*----------------------------------------------------------------------------*
   * Methods  about mesh                                                        *
   *----------------------------------------------------------------------------*/




  /*----------------------------------------------------------------------------*
   * Methods about field                                                        *
   *----------------------------------------------------------------------------*/

  /**
   *
   * \brief Create a new field
   *
   * \param [in]  field_id       Field id
   * \param [in]  data_type      Data type
   * \param [in]  storage        Storage type
   * \param [in]  n_component    Number of componenent
   * \param [in]  nature         Nature
   * \param [in]  exch_type      Exchange type
   * \param [in]  visu_status    Visualization status
   *
   */

  void
  Coupling::fieldCreate
  (
   const string               &field_id,
   const CWP_Type_t           data_type,
   const CWP_Field_storage_t  storage,
   const int                  n_component,
   const CWP_Dof_location_t    fieldType,
   const CWP_Field_exch_t     exch_type,
   const CWP_Status_t         visu_status
  )
  {

    _is_first_field_created = 1;

    if (fieldIs(field_id)) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' existing field\n", field_id.c_str());
    }

    //
    // Create the new field

    cwipi::Field *newField = new cwipi::Field(field_id,
                                              _fields.size(),
                                              data_type,
                                              this,
                                              fieldType,
                                              storage,
                                              n_component,
                                              exch_type,
                                              visu_status);  

    pair<string, Field* > newPair(string(field_id), newField);
    string localName = _localCodeProperties.nameGet();
    _fields.insert(newPair);



    // if (_visu.isCreated() && newField->visuStatusGet() == CWP_STATUS_ON) {
    //   _visu.WriterFieldCreate(newField);
    // }
  }


  void
  Coupling::fieldPythonObjectSet
  (
   const string &field_id,
   void   *p
   )
  {
    map<string,Field*>::iterator it = _fields.find(field_id.c_str());

    assert(it != _fields.end());

    it->second->pythonObjectSet(p);
  }


  /**
   * \brief Return if a field identifier exists
   *
   * \param [in]  field_id         Field identifier
   *
   * \return status
   */


  bool
  Coupling::fieldIs
  (
   const string &field_id
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    return (It != _fields.end());
  }


 /**
  *
  * \brief Set Field data
  *
  * \param [in]  field_id       Field identifier
  * \param [in]  data           Storage array (mapping)
  *
  */

  void
  Coupling::fieldDataSet
  (
    const string &field_id,
    int i_part,
    const CWP_Field_map_t   map_type,
    void* data
  )
  {

    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      It->second->dataSet(i_part, map_type, data);
      // if (_visu.isCreated() && It->second->visuStatusGet() == CWP_STATUS_ON) {
      //   _visu.fieldDataSet(It->second,map_type, i_part);
      // }
    }
  }

  /**
  *
  * \brief Get Field data
  *
  * \param [in]   field_id       Field identifier
  * \param [out]  data           Storage array (mapping)
  *
  */

  void
  Coupling::fieldDataGet
  (
    const string &field_id,
    int i_part,
    const CWP_Field_map_t   map_type,
    void** data
  )
  {

    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' not existing field\n", field_id.c_str());
    }
    else {
      *data = It->second->dataGet(i_part, map_type);
    }
  }

  /**
   *
   * \brief Get nunmber of field degrees of freedom
   *
   * \param [in]   field_id       Field identifier
   * \param [in]   i_part         Partition identifier
   *
   * \return                      Number of field degrees of freedom
   *
   */

  int
  Coupling::fieldNDOFGet
  (
    const string &field_id,
    int          i_part
  )
  {
    map<string,Field*>::iterator It = _fields.find(field_id.c_str());
    if (It == _fields.end()) {
       PDM_error(__FILE__, __LINE__, 0, "'%s' not existing field\n", field_id.c_str());
    }
    return It->second->nDOFGet(i_part);
  }

  /*----------------------------------------------------------------------------*
   * Methods about exchange                                                     *
   *----------------------------------------------------------------------------*/


  /**
   * \brief data exchange <b>(Not implemented yet)</b>
   *
   * Exchange depending on exchange frequency
   *
   */

  void
  Coupling::exchange ()
  {
    PDM_error(__FILE__, __LINE__, 0, "\nexchange not implemented yet\n");
  }


  /**
   * \brief Exchange data field with the coupled code with blocking
   *        communications. <b>(Not implemented yet)</b>
   *
   * This function exchanges interpolated fieldsDouble between coupled codes.
   *
   * \warning  The size of tgt_field_id size is n_computed_tgt.
   *           If \f$ n\_uncomputed\_tgt \ne n\_tgt\_pts \f$,
   *           user himself must set values for uncomputed target points.
   *
   * \param [in]  src_field_id              Source field (NULL->no sending)
   * \param [in]  tgt_field_id              Target field (NULL->no receiving)
   * \param [in]  ptFortranInterpolationFct Fortran user interpolation (or NULL)
   * \param [out] n_uncomputed_tgt          Number of uncomputed target
   *
   * \return                                Exchange status
   *
   */

  void
  Coupling::sendrecv (
    const string &field_id
  )
  {
    PDM_UNUSED (field_id);
    PDM_error(__FILE__, __LINE__, 0, "\nsendrecv not implemented yet\n");
  }


  /**
   *
   * \brief Sending of data field to the coupled code with nonblocking
   *        communications.
   *
   * This function sends interpolated field to the coupled code.
   *
   * \param [in]  src_id                    Source field
   *
   */

  void
  Coupling::issend (
    const string &sendingFieldID
  )
  {

    if (!_coupledCodeProperties.localCodeIs()) {

      map <string, Field *>::iterator it;
      it = _fields.find(sendingFieldID);

      if (it != _fields.end()) {
        Field* sendingField = it->second;

        CWP_Dof_location_t localFieldLocation = sendingField->locationGet();
        CWP_Dof_location_t cplFieldLocation = sendingField->linkedFieldLocationGet();

        std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;

        it2 = _spatial_interp_send.find(newKey);

        if (it2 == _spatial_interp_send.end()) {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
        }
        else {
          it2->second->issend(sendingField);
        }






      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
      }
    }

    else {

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        map <string, Field *>::iterator it;
        it = _fields.find(sendingFieldID);

        map <string, Field *>::iterator cpl_it;
        cpl_it = cpl_cpl._fields.find(sendingFieldID);

        Field* recvField = NULL;
        if (cpl_it != cpl_cpl._fields.end()) {
          recvField = it->second;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }

        if (it != _fields.end()) {
          Field* sendingField = it->second;

          CWP_Dof_location_t localFieldLocation = sendingField->locationGet();
          CWP_Dof_location_t cplFieldLocation = sendingField->linkedFieldLocationGet();

          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);
          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocation, localFieldLocation);

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator cpl_it2;

          it2 = _spatial_interp_send.find(newKey);
          cpl_it2 = cpl_cpl._spatial_interp_recv.find(cpl_newKey);

          if (cpl_it2 == cpl_cpl._spatial_interp_recv.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            cpl_it2->second->irecv(recvField);
          }

          if (it2 == _spatial_interp_send.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            it2->second->issend(sendingField);
          }
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }
      }
    }
  }


  /**
   *
   * \brief Waiting of the end of exchange related to request.
   *
   * This function waits the end of exchange related to request
   * from \ref CWP_Issend
   *
   * \param [in] src_id                    Source field
   *
   */

  void
  Coupling::waitIssend (
    const string &sendingFieldID
  )
  {
    if (!_coupledCodeProperties.localCodeIs()) {
      map <string, Field *>::iterator it;
      it = _fields.find(sendingFieldID);

      if (it != _fields.end()) {
        Field* sendingField = it->second;

        CWP_Dof_location_t localFieldLocation = sendingField->locationGet();
        CWP_Dof_location_t cplFieldLocation = sendingField->linkedFieldLocationGet();

        std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;

        it2 = _spatial_interp_send.find(newKey);

        if (it2 == _spatial_interp_send.end()) {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
        }
        else {
          it2->second->waitIssend(sendingField);
        }
      }
    }
    else {
      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {
        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        map <string, Field *>::iterator it;
        it = _fields.find(sendingFieldID);

        map <string, Field *>::iterator cpl_it;
        cpl_it = cpl_cpl._fields.find(sendingFieldID);

        Field* recvField = NULL;
        if (cpl_it != cpl_cpl._fields.end()) {
          recvField = cpl_it->second;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }

        if (it != _fields.end()) {
          Field* sendingField = it->second;

          CWP_Dof_location_t localFieldLocation = sendingField->locationGet();
          CWP_Dof_location_t cplFieldLocation = sendingField->linkedFieldLocationGet();

          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);
          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocation, localFieldLocation);

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator cpl_it2;

          it2 = _spatial_interp_send.find(newKey);
          cpl_it2 = cpl_cpl._spatial_interp_recv.find(cpl_newKey);

          if (cpl_it2 == cpl_cpl._spatial_interp_recv.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            cpl_it2->second->waitIrecv(recvField);
          }

          if (it2 == _spatial_interp_send.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            it2->second->waitIssend(sendingField);
          }
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }
      }
    }
  }


  /**
   *
   * \brief Receiving of Data field from the coupled code with nonblocking
   *        communications.
   *
   * This function receives interpolated field from the coupled code
   *
   * \param [in]  receving_field_id       Target field ID
   *
   *
   */

  void
  Coupling::irecv
  (
    const string &receivingFieldID
  )
  {
    if (!_coupledCodeProperties.localCodeIs()) {

      map <string, Field *>::iterator it;
      it = _fields.find(receivingFieldID);

      if (it != _fields.end()) {
        Field* receivingField = it->second;

        CWP_Dof_location_t localFieldLocation = receivingField->locationGet();
        CWP_Dof_location_t cplFieldLocation = receivingField->linkedFieldLocationGet();

        std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;

        it2 = _spatial_interp_recv.find(newKey);

        if (it2 == _spatial_interp_recv.end()) {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
        }
        else {
          it2->second->irecv(receivingField);
        }
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
      }
    }

    else {

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        map <string, Field *>::iterator it;
        it = _fields.find(receivingFieldID);

        map <string, Field *>::iterator cpl_it;
        cpl_it = cpl_cpl._fields.find(receivingFieldID);

        Field* sendField = NULL;
        if (cpl_it != cpl_cpl._fields.end()) {
          sendField = it->second;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }

        if (it != _fields.end()) {
          Field* receivingField = it->second;

          CWP_Dof_location_t localFieldLocation = receivingField->locationGet();
          CWP_Dof_location_t cplFieldLocation = receivingField->linkedFieldLocationGet();

          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);
          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocation, localFieldLocation);

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator cpl_it2;

          it2 = _spatial_interp_recv.find(newKey);
          cpl_it2 = cpl_cpl._spatial_interp_send.find(cpl_newKey);

          if (it2 == _spatial_interp_recv.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            it2->second->irecv(receivingField);
          }

          if (cpl_it2 == cpl_cpl._spatial_interp_send.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            cpl_it2->second->issend(sendField);
          }

        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }
      }
    }

  }


  /**
   *
   * \brief Waiting of the end of exchange related to request.
   *
   * This function waits the end of exchange related to request
   * from \ref CWP_Irecv
   *
   * \param [in]  receving_field_id       Target field ID
   *
   */

  void
  Coupling::waitIrecv (
    const string &receivingFieldID
  )
  {
    if (!_coupledCodeProperties.localCodeIs()) {

      map <string, Field *>::iterator it;
      it = _fields.find(receivingFieldID);

      if (it != _fields.end()) {
        Field* receivingField = it->second;

        CWP_Dof_location_t localFieldLocation = receivingField->locationGet();
        CWP_Dof_location_t cplFieldLocation = receivingField->linkedFieldLocationGet();

        std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);

        std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;

        it2 = _spatial_interp_recv.find(newKey);

        if (it2 == _spatial_interp_recv.end()) {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
        }
        else {
          it2->second->waitIrecv(receivingField);
        }
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
      }
    }

    else {

      if (_localCodeProperties.idGet() < _coupledCodeProperties.idGet()) {

        cwipi::Coupling& cpl_cpl = _cplDB.couplingGet (_coupledCodeProperties, _cplId);

        map <string, Field *>::iterator it;
        it = _fields.find(receivingFieldID);

        map <string, Field *>::iterator cpl_it;
        cpl_it = cpl_cpl._fields.find(receivingFieldID);

        Field* sendField = NULL;
        if (cpl_it != cpl_cpl._fields.end()) {
          sendField = cpl_it->second;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }

        if (it != _fields.end()) {
          Field* receivingField = it->second;

          CWP_Dof_location_t localFieldLocation = receivingField->locationGet();
          CWP_Dof_location_t cplFieldLocation = receivingField->linkedFieldLocationGet();

          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > newKey (localFieldLocation, cplFieldLocation);
          std::pair < CWP_Dof_location_t, CWP_Dof_location_t > cpl_newKey (cplFieldLocation, localFieldLocation);

          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator it2;
          std::map < std::pair < CWP_Dof_location_t, CWP_Dof_location_t >, SpatialInterp*>::iterator cpl_it2;

          it2 = _spatial_interp_recv.find(newKey);
          cpl_it2 = cpl_cpl._spatial_interp_send.find(cpl_newKey);

          if (it2 == _spatial_interp_recv.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            it2->second->waitIrecv(receivingField);
          }

          if (cpl_it2 == cpl_cpl._spatial_interp_send.end()) {
            PDM_error(__FILE__, __LINE__, 0, "\nUnknown spatial interpolation\n");
          }
          else {
            cpl_it2->second->waitIssend(sendField);
          }

        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "\nUnknown field\n");
        }
      }
    }
  }


  /*----------------------------------------------------------------------------*
   * methods about user interpolation                                           *
   *----------------------------------------------------------------------------*/


  /*----------------------------------------------------------------------------*
   * Private methods                                                            *
   *----------------------------------------------------------------------------*/

  /**
   *
   * \brief Compute user target global number (if not given by user)
   *
   */

  void
  Coupling::userTargetGnumCompute()
  {
    if (_userTargetN != nullptr) {
      if (_userTargetGnum == nullptr) {

        PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm ((void *) &(_localCodeProperties.intraCommGet()));

        PDM_gen_gnum_t *pgg  = PDM_gnum_create (3, _nPart, PDM_FALSE, 1e-3, comm,
                                                   PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);
        for (int iPart = 0; iPart <_nPart; iPart++){
          PDM_gnum_set_from_coords (pgg, iPart, _userTargetN[iPart], _userTargetCoord[iPart], NULL);
        }

        PDM_gnum_compute (pgg);

        _localUserTargetGnum = (CWP_g_num_t **) malloc (sizeof(CWP_g_num_t *) * _nPart);

        for (int iPart = 0; iPart < _nPart; iPart++){
          _localUserTargetGnum[iPart] = const_cast <CWP_g_num_t*> (PDM_gnum_get (pgg, iPart));
        }

        PDM_gnum_free (pgg);

        _userTargetGnum = const_cast <const CWP_g_num_t**> (_localUserTargetGnum);
      }
    }
  }

  /**
   * \brief Update time.
   *
   * \param [in]  current_time     Current time
   *
   */

  void
  Coupling::timeUpdate (
    double current_time
  )
  {
    _n_step++;

    if (_writer != NULL) {
      PDM_writer_step_end(_writer);
    }

    // if(_visu.isCreated() and _visu.physicalTimeGet() > -1.0) {
    //    _visu.WriterStepEnd();
    // }

    std::map < string, Field * >::iterator itf = _fields.begin();
    while (itf != _fields.end()) {
      // itf->second->currentStepWasExchangedReset(); // TO DO
      itf++;
    }

    if (_writer != NULL && (_n_step % _freq_writer == 0)) {
      PDM_writer_step_beg(_writer, current_time);
    }

    // if(_visu.isCreated()) {
    //    _visu.WriterStepBegin(current_time, &_mesh);
    // }

  }

  // Begin code time step

  void
  Coupling::time_step_beg (
    double current_time
  )
  {
    // if first time step do not PDM_writer_step_beg, has to be done after variable
    // creation in exportMesh
    // open time step if there is a writer, no step is open and it is a writting step
    if (_n_step != 0 && _writer != NULL && PDM_writer_is_open_step(_writer) == 0 && (_n_step % _freq_writer == 0)) {
      PDM_writer_step_beg(_writer, current_time);
    }
  }

  // End code time step

  void
  Coupling::time_step_end ()
  {
    // if there is a writer, a step is open and it is a writting step
    if (_writer != NULL && PDM_writer_is_open_step(_writer) == 1 && (_n_step % _freq_writer == 0)) {

      // some unexchanged visu ON
      std::map < string, Field * >::iterator it_f1 = _fields.begin();
      while (it_f1 != _fields.end()) {
        // field unexchanged visu ON
        if (it_f1->second->visuStatusGet() == CWP_STATUS_ON) {
          CWP_Field_exch_t exchange_type = it_f1->second->exchangeTypeGet();
          if (exchange_type == CWP_FIELD_EXCH_SEND || exchange_type == CWP_FIELD_EXCH_SENDRECV) {
            // has not been sent
            if (it_f1->second->is_send_yet_get() == 0) {
              it_f1->second->is_send_end_step_set(1);
              it_f1->second->write(CWP_FIELD_EXCH_SEND);
              it_f1->second->is_send_end_step_set(0);
            }
          }
          else if (exchange_type == CWP_FIELD_EXCH_RECV || exchange_type == CWP_FIELD_EXCH_SENDRECV) {
            // has not been received
            if (it_f1->second->is_recv_yet_get() == 0) {
              it_f1->second->is_recv_end_step_set(1);
              it_f1->second->write(CWP_FIELD_EXCH_RECV);
              it_f1->second->is_recv_end_step_set(0);
            }
          }
        }
        it_f1++;
      }

      // end step
      PDM_writer_step_end(_writer);
    }

    // increment step
    _n_step++;

    // reset field status to unexchanged
    std::map < string, Field * >::iterator it_f2 = _fields.begin();
    while (it_f2 != _fields.end()) {
      it_f2->second->is_send_yet_set(0);
      it_f2->second->is_recv_yet_set(0);
      it_f2++;
    }
  }

// A supprimer

  CWP_g_num_t*
  Coupling::globalNumGet(int id_block,int i_part)
  {
    return _mesh.globalNumGet(id_block,i_part);
  }

  // int
  // Coupling::isUpToDateGet ()
  // {
  //   return _is_up_to_date;
  // }

  // void
  // Coupling::isUpToDateSet ()
  // {
  //   _is_up_to_date = 1;
  // }



} // namespace cwipi

/**
 * \endcond
 */
