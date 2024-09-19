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

#include "communication.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "pdm_printf.h"
#include "pdm_logging.h"

/**
 * \cond
 */

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */

  Communication::Communication()
    :_localCodeProperties(NULL),
     _cplCodeProperties(NULL),
     _tag(-1),
     _unionGroup(MPI_GROUP_NULL),
     _unionComm(MPI_COMM_NULL),
     _unionCommCplRanks(NULL),
     _unionCommLocRanks(NULL),
     _cplCommCplRanks(NULL),
     _cplCommLocRanks(NULL),
     _cplGroup(MPI_GROUP_NULL),
     _cplComm(MPI_COMM_NULL),
     _locCodeRootRankUnionComm(-1),
     _cplCodeRootRankUnionComm(-1),
     _locCodeRootRankCplComm(-1),
     _cplCodeRootRankCplComm(-1),
     _isCplRank(false)
  {
  }

  /**
   *
   * \brief Destructor.
   *
   */

  Communication::~Communication()
  {

    if (_cplComm != _unionComm) {
      if (!_localCodeProperties->isActiveRank() ||
          _localCodeProperties->idGet() < _cplCodeProperties->idGet()) {
        if (_cplComm != MPI_COMM_NULL) {
          MPI_Comm_free(&_cplComm);
          _cplComm = MPI_COMM_NULL;
        }

        if (_unionComm != MPI_COMM_NULL) {
          MPI_Comm_free(&_unionComm);
        }
      }
    }

    else {
      if (!_localCodeProperties->isActiveRank()) {
        if (_cplComm != MPI_COMM_NULL) {
          MPI_Comm_free(&_cplComm);
          _cplComm = MPI_COMM_NULL;
        }
      }
      else {
        if (_localCodeProperties->idGet() < _cplCodeProperties->idGet()) {
          if (_cplComm != MPI_COMM_NULL) {
            MPI_Comm_free(&_cplComm);
            _cplComm = MPI_COMM_NULL;
          }
        }       
      }
    }
    
//    if (_unionGroup != MPI_GROUP_NULL) {
//      MPI_Group_free(&_unionGroup);
//    }


    if (_unionCommCplRanks != NULL) {
      delete _unionCommCplRanks;
    }

    if (_unionCommLocRanks != NULL) {
      delete _unionCommLocRanks;
    }

    if (_cplCommCplRanks != NULL) {
      delete _cplCommCplRanks;
    }

    if (_cplCommLocRanks != NULL) {
      delete _cplCommLocRanks;
    }


//    if (_fvmComm != MPI_COMM_NULL)
//      MPI_Comm_free(&_fvmComm);
  }

  /**
   *
   * \brief Initialize coupling communicators.
   *
   * \param [in]  localCodeProperties   Local code properties
   * \param [in]  cplCodeProperties     Coupled code properties
   * \param [in]  cplId                 Coupling identifier
   *
   */

  void
  Communication::init
  (
   const CodeProperties &localCodeProperties,
   const CodeProperties &cplCodeProperties,
   const string         &cplId,
   CouplingDB           &cplDB
  )
  {

    _isCplRank = localCodeProperties.isActiveRank();
    if (_cplComm == MPI_COMM_NULL) {
      _cplCodeProperties = &cplCodeProperties;
      _localCodeProperties = &localCodeProperties;

      const int localRootRank = localCodeProperties.rootRankGet();
      const int cplRootRank = cplCodeProperties.rootRankGet();

      const MPI_Comm& globalComm = localCodeProperties.globalCommGet();

      int globalRank;
      MPI_Comm_rank(globalComm, &globalRank);

      if (!localCodeProperties.isActiveRank()) {
        PDM_printf(
           "Warning CWP_Cpl_create : Call CWP_Cpl_create function"
           " on an uncoupled rank (%d) of the '%s' code\n",
            globalRank,
            localCodeProperties.nameGet().c_str());
      }
      else {

        //Build a specific tag through the cplId
        _tag = 0;
        for (size_t i = 0; i < cplId.size(); i++) {
          _tag += cplId[i];
        }

        CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wdiv-by-zero")

        MPI_Aint  *maxTagTmp;
        int flag; 

        MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &maxTagTmp, &flag);
        int maxTag = (int) *maxTagTmp; 

        if (maxTag > 0) {
          _tag = _tag % maxTag;
        }
        CWP_GCC_SUPPRESS_WARNING_POP

        // Build the union communicator between the two coupled codes

        if (localCodeProperties.idGet() < cplCodeProperties.idGet()) {
          MPI_Group_union (localCodeProperties.connectableGroupGet(),
                           cplCodeProperties.connectableGroupGet(),
                           &_unionGroup);
        }
        else {
          MPI_Group_union (cplCodeProperties.connectableGroupGet(),
                           localCodeProperties.connectableGroupGet(),
                           &_unionGroup);
        }

        MPI_Comm_create_group(globalComm, _unionGroup, _tag, &_unionComm);


        int mergeInterCommSize;

        MPI_Comm_size(_unionComm, &mergeInterCommSize);

        CWP_Comm_t commType = commTypeGet();

        CWP_Comm_t cplCommType;

        if (globalRank == localRootRank) {
          if (cplCodeProperties.localCodeIs()) {
            if (localRootRank == cplRootRank) {
              Coupling &distCpl = cplDB.couplingGet(cplCodeProperties, cplId);
              cplCommType = distCpl.commTypeGet() ;
            }
            else {
              MPI_Sendrecv (&commType, 1, MPI_INT,
                          cplRootRank, _tag,
                          &cplCommType, 1, MPI_INT,
                           cplRootRank, _tag,
                           globalComm, MPI_STATUS_IGNORE);
            }
          }

          else {
            MPI_Sendrecv (&commType, 1, MPI_INT,
                          cplRootRank, _tag,
                          &cplCommType, 1, MPI_INT,
                           cplRootRank, _tag,
                           globalComm, MPI_STATUS_IGNORE);
          }
        }

        MPI_Request request1;
        MPI_Request request2;

        MPI_Ibcast(&cplCommType, 1, MPI_INT, 0,
                  localCodeProperties.connectableCommGet(), &request1);
        MPI_Wait(&request1, MPI_STATUS_IGNORE);


        if (cplCodeProperties.localCodeIs()) {
          CWP_Comm_t &cplcommType2 = commType;

          MPI_Ibcast(&cplcommType2, 1, MPI_INT, 0,
                    cplCodeProperties.connectableCommGet(), & request2);
          MPI_Wait(&request2, MPI_STATUS_IGNORE);
        }

        //
        // Build the coupling communicator

        _cplCommCreate(cplCommType);
        MPI_Group globalGroup;
        MPI_Comm_group(globalComm, &globalGroup);

        if (_isCplRank) {
          MPI_Group_translate_ranks (globalGroup, 1, &localRootRank,
                                     _unionGroup, &_locCodeRootRankUnionComm);

          MPI_Group_translate_ranks (globalGroup, 1, &cplRootRank,
                                     _unionGroup, &_cplCodeRootRankUnionComm);
        }

      }
    }
  }


  /**
   *
   * \brief Initialize coupling communicators.
   *
   * \param [in]  cplCodeComm           Coupled code communication
   *
   */

  void
  Communication::init
  (
   Communication &cplCodeComm
  )
  {
    _localCodeProperties    = cplCodeComm._cplCodeProperties;
    _cplCodeProperties      = cplCodeComm._localCodeProperties;

    _tag                    = cplCodeComm._tag;
    _unionGroup             = cplCodeComm._unionGroup;
    _unionComm              = cplCodeComm._unionComm;
    _cplGroup               = cplCodeComm._cplGroup;
    _cplComm                = cplCodeComm._cplComm;
    _locCodeRootRankUnionComm = cplCodeComm._cplCodeRootRankUnionComm;
    _cplCodeRootRankUnionComm = cplCodeComm._locCodeRootRankUnionComm;
    _unionCommCplRanks = new std::vector<int>(*(cplCodeComm._unionCommLocRanks));
    _unionCommLocRanks = new std::vector<int>(*(cplCodeComm._unionCommCplRanks));
    _cplCommCplRanks = new std::vector<int>(*(cplCodeComm._cplCommLocRanks));
    _cplCommLocRanks = new std::vector<int>(*(cplCodeComm._cplCommCplRanks));


    const int localRootRank = _localCodeProperties->rootRankGet();

    const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();

    int currentRank;
    MPI_Comm_rank(globalComm, &currentRank);

    if (this->commTypeGet() != CWP_COMM_PAR_WITH_PART) {
      _isCplRank = localRootRank == currentRank;
    }
    else {
      _isCplRank = true;
    }

  }

  MPI_Comm
  Communication::unionCommGet
  (
  )
  {
    return _unionComm;
  }

  MPI_Comm
  Communication::cplCommGet
  (
  )
  {
    return _cplComm;
  }


  std::vector<int>*
  Communication::unionCommCplRanksGet
  (
  )
  {
    return _unionCommCplRanks;
  }

  std::vector<int>*
  Communication::cplCommCplRanksGet
  (
  )
  {
    return _cplCommCplRanks;
  }

  std::vector<int>*
  Communication::unionCommLocRanksGet
  (
  )
  {
    return _unionCommLocRanks;
  }

  std::vector<int>*
  Communication::cplCommLocRanksGet
  (
  )
  {
    return _cplCommLocRanks;
  }

  int
  Communication::unionCommCplCodeRootRanksGet
  (
  )
  {
     return  _cplCodeRootRankUnionComm;
  }

  int
  Communication::cplCommCplCodeRootRanksGet
  (
  )
  {
     return  _cplCodeRootRankCplComm;
  }

  int
  Communication::unionCommLocCodeRootRanksGet
  (
  )
  {
    return _locCodeRootRankUnionComm;
  }

  int
  Communication::cplCommLocCodeRootRanksGet
  (
  )
  {
    return _locCodeRootRankCplComm;
  }

uint32_t Communication::_adler32
(
 const void *buf,
 size_t buflength
 )
{

  const uint8_t * buffer = (const uint8_t *)buf;

  uint32_t s1 = 1;
  uint32_t s2 = 0;

  for (size_t n = 0; n < buflength; n++) {
    s1 = (s1 + buffer[n]) % 65521;
    s2 = (s2 + s1) % 65521;
  }

  return (s2 << 16) | s1;
}

uint32_t Communication::_get_tag
(
const string    name,
MPI_Comm        comm,
int             offset
)
{
  MPI_Aint  *maxTagTmp;
  int flag;

  MPI_Comm_get_attr(comm, MPI_TAG_UB, &maxTagTmp, &flag);
  int maxTag = (int) *maxTagTmp;

  uint32_t mpi_tag = ((_adler32 (name.c_str(),
        name.size()) + offset)% (maxTag - 1)) + 1;

  return mpi_tag;
}


  /**
   *
   * \brief Exchange Data between two coupled code through Union Comminucator
   *        All code ranks receive data (sendRecv between root rank then Bcast)
   *
   * \param [in]    s_data           Size of a data tot exchange
   * \param [in]    n_send_data      Number of data to send
   * \param [in]    send_data        Array of data to send
   * \param [in]    n_send_data_cpl  Number of data to send (if current rank is shared betweenn code en coupled code)
   * \param [in]    send_data_cpl    Array of data to send (if current rank is shared betweenn code en coupled code)
   * \param [in]    n_recv_data      Number of data to receive
   * \param [inout] recv_data        Array of data to receive
   * \param [in]    n_recv_data_cpl  Number of data to receive (if current rank is shared betweenn code en coupled code)
   * \param [inout] recv_data_cpl    Array of data to receive (if current rank is shared betweenn code en coupled code)
   *
   */

  void
  Communication::iexchGlobalDataBetweenCodesThroughUnionCom
  (
   size_t       s_data,
   int          n_send_data,
   void        *send_data,
   int          n_send_data_cpl,
   void        *send_data_cpl,
   int          n_recv_data,
   void        *recv_data,
   int          n_recv_data_cpl,
   void        *recv_data_cpl
  )
  {
    int debug = 0;

    if (debug) {
      printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 1\n");
      fflush(stdout);
    }

    int unionCommRank;
    MPI_Status status;

    if (_cplCodeProperties->localCodeIs()) {
      if (n_send_data_cpl > 0) {
        assert(send_data_cpl != NULL);
      }
      if (n_recv_data_cpl > 0) {
        assert(recv_data_cpl != NULL);
      }
    }
    else {
      assert(send_data_cpl == NULL);
      assert(recv_data_cpl == NULL);
    }

    MPI_Comm_rank (_unionComm, &unionCommRank);

    if (debug) {
      printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 2\n");
      fflush(stdout);
    }

    int codeID    = _localCodeProperties->idGet();
    int cplCodeID = _cplCodeProperties->idGet();
    const int localRootRank = _localCodeProperties->rootRankGet();

    const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();

    int currentRank;
    MPI_Comm_rank(globalComm, &currentRank);

    if (debug) {
      printf("%s - unionCommRank _locCodeRootRankUnionComm _locCodeRootRankCplComm rankglobalcomm localRootRank %d %d %d %d %d\n", _localCodeProperties->nameGet().c_str(),
       unionCommRank, 
       _locCodeRootRankUnionComm,
       _locCodeRootRankCplComm,
       currentRank,
       localRootRank);
    }

    if (unionCommRank ==  _locCodeRootRankUnionComm) {
      if (_cplCodeProperties->localCodeIs()) {

        if (_locCodeRootRankUnionComm == _cplCodeRootRankUnionComm) {
          if (debug) {
            printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 2.1\n");
            fflush(stdout);
          }

          if (codeID < cplCodeID) {

            assert (n_send_data     == n_recv_data_cpl);
            memcpy (recv_data_cpl, send_data    , s_data * n_send_data);

            assert (n_send_data_cpl  == n_recv_data);
            memcpy (recv_data,      send_data_cpl    , s_data * n_send_data_cpl);
          }

        }
        else {
          if (debug) {
            printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 2.2\n");
            fflush(stdout);
          }
          MPI_Sendrecv (send_data,
                        (int) s_data * n_send_data,
                        MPI_UNSIGNED_CHAR,
                        _cplCodeRootRankUnionComm,
                        _tag,
                        recv_data,
                        (int) s_data * n_recv_data,
                        MPI_UNSIGNED_CHAR,
                        _cplCodeRootRankUnionComm,
                        _tag,
                        _unionComm,
                        &status);
        }
      }
      else {
        if (debug) {
          printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 2.3\n");
          fflush(stdout);
        }
        MPI_Sendrecv (send_data,
                      (int) s_data * n_send_data,
                      MPI_UNSIGNED_CHAR,
                      _cplCodeRootRankUnionComm,
                      _tag,
                      recv_data,
                      (int) s_data * n_recv_data,
                      MPI_UNSIGNED_CHAR,
                      _cplCodeRootRankUnionComm,
                      _tag,
                      _unionComm,
                      &status);
      }
    }
    else {
      if (_cplCodeProperties->localCodeIs()) {

        if (unionCommRank == _cplCodeRootRankUnionComm) {
          if (debug) {
            printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 2.7\n");
            fflush(stdout);
          }

          if (codeID < cplCodeID) {

            MPI_Sendrecv (send_data_cpl,
                          (int) s_data * n_send_data_cpl,
                          MPI_UNSIGNED_CHAR,
                          _locCodeRootRankUnionComm,
                          _tag,
                          recv_data_cpl,
                          (int) s_data * n_recv_data_cpl,
                          MPI_UNSIGNED_CHAR,
                          _locCodeRootRankUnionComm,
                          _tag,
                          _unionComm,
                          &status);
          }
        }
      }

    }

    if (debug) {
      printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 4.1\n");
      fflush(stdout);
    }

    if (!_cplCodeProperties->localCodeIs()) {

      MPI_Bcast (recv_data,
                 (int) s_data * n_recv_data,
                 MPI_UNSIGNED_CHAR,
                 0,
                 _localCodeProperties->connectableCommGet());
      int rank1;
      MPI_Comm_rank(_localCodeProperties->connectableCommGet(), &rank1);
      int s1;
      MPI_Comm_size(_localCodeProperties->connectableCommGet(), &s1);
    }

    else {

      if (codeID < cplCodeID) {
        MPI_Bcast (recv_data,
                   (int) s_data * n_recv_data,
                   MPI_UNSIGNED_CHAR,
                   0,
                   _localCodeProperties->connectableCommGet());
        int rank1;
        MPI_Comm_rank(_localCodeProperties->connectableCommGet(), &rank1);
        int s1;
        MPI_Comm_size(_localCodeProperties->connectableCommGet(), &s1);

        if (debug) {
          printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 4.2\n");
          fflush(stdout);
        }

        MPI_Bcast (recv_data_cpl,
                   (int) s_data * n_recv_data_cpl,
                   MPI_UNSIGNED_CHAR,
                   0,
                   _cplCodeProperties->connectableCommGet());
        int rank2;
        int s2;
        MPI_Comm_rank(_cplCodeProperties->connectableCommGet(), &rank2);
        MPI_Comm_size(_cplCodeProperties->connectableCommGet(), &s2);
      }
      else {
        if (debug) {
          printf("Communication::iexchGlobalDataBetweenCodesThroughUnionCom - 4.2\n");
          fflush(stdout);
        }

        MPI_Bcast (recv_data_cpl,
                   (int) s_data * n_recv_data_cpl,
                   MPI_UNSIGNED_CHAR,
                   0,
                   _cplCodeProperties->connectableCommGet());
        int rank2;
        int s2;
        MPI_Comm_rank(_cplCodeProperties->connectableCommGet(), &rank2);
        MPI_Comm_size(_cplCodeProperties->connectableCommGet(), &s2);
      }
    }
  }

  /**
   *
   * \brief Non blocking send of global data array
   *
   * \param [in]  global_data_id
   * \param [out] data_send_request
   * \param [in]  s_send_entity
   * \param [in]  send_stride
   * \param [in]  n_send_entity
   * \param [in]  send_data
   */

  void
  Communication::issendGlobalDataBetweenCodesThroughUnionCom
  (
   const string    global_data_id,
   MPI_Request    *data_send_request,
   size_t          s_entity,
   int             stride,
   int             n_entity,
   void           *send_data
  )
  {
    assert(send_data != NULL);

    // Get union communicator ie. union of all active ranks of the codes in the coupling
    int unionCommRank;
    MPI_Comm_rank (_unionComm, &unionCommRank);

    // Get couplings id
    int codeID    = _localCodeProperties->idGet();
    int cplCodeID = _cplCodeProperties->idGet();

    // Isend
    if (unionCommRank ==  _locCodeRootRankUnionComm) {
      if (_cplCodeProperties->localCodeIs()) {
        if (_locCodeRootRankUnionComm == _cplCodeRootRankUnionComm) {
          if (codeID < cplCodeID) {

            // Nothing

          } // end code with smallest id does the action of the rank
        } // end if i_rank is root rank of coupled code
        else {
          uint32_t data_tag = _get_tag(global_data_id,
                                       _unionComm,
                                       0);

          MPI_Issend(send_data,
                     (int) s_entity * stride * n_entity,
                     MPI_UNSIGNED_CHAR,
                     _cplCodeRootRankUnionComm,
                     (int) data_tag,
                     _unionComm,
                     data_send_request);
        }
      } // end if i_rank is joint with the coupled code
      else {
        uint32_t data_tag = _get_tag(global_data_id,
                                     _unionComm,
                                     0);

        MPI_Issend(send_data,
                   (int) s_entity * stride * n_entity,
                   MPI_UNSIGNED_CHAR,
                   _cplCodeRootRankUnionComm,
                   (int) data_tag,
                   _unionComm,
                   data_send_request);
      }
    } // end if i_rank is the root rank of the local code
  }


  /**
   *
   * \brief Non blocking receive of global data array
   *
   * \param [in]  global_data_id
   * \param [out] data_recv_request
   * \param [in]  s_entity
   * \param [in]  stride
   * \param [in]  n_entity
   * \param [out] recv_data
   */

  void
  Communication::irecvGlobalDataBetweenCodesThroughUnionCom
  (
   const string    global_data_id,
   MPI_Request    *data_recv_request,
   size_t          s_entity,
   int             stride,
   int             n_entity,
   void           *recv_data
  )
  {
    assert(recv_data != NULL);

    // Get union communicator ie. union of all active ranks of the codes in the coupling
    int unionCommRank;
    MPI_Comm_rank (_unionComm, &unionCommRank);

    // Get couplings id
    int codeID    = _localCodeProperties->idGet();
    int cplCodeID = _cplCodeProperties->idGet();

    // Irecv
    if (unionCommRank ==  _locCodeRootRankUnionComm) {
      if (_cplCodeProperties->localCodeIs()) {
        if (_locCodeRootRankUnionComm == _cplCodeRootRankUnionComm) {
          if (codeID < cplCodeID) {
            // Nothing
          } // end code with smallest id does the action of the rank
        } // end if i_rank is root rank of coupled code
        else {
          uint32_t data_tag = _get_tag(global_data_id,
                                       _unionComm,
                                       0);
          MPI_Irecv(recv_data,
                    (int) s_entity * stride * n_entity,
                    MPI_UNSIGNED_CHAR,
                    _cplCodeRootRankUnionComm,
                    (int) data_tag,
                    _unionComm,
                    data_recv_request);
        }
      } // end if i_rank is joint with the coupled code
      else {
        uint32_t data_tag = _get_tag(global_data_id,
                                     _unionComm,
                                     0);
        MPI_Irecv(recv_data,
                  (int) s_entity * stride * n_entity,
                  MPI_UNSIGNED_CHAR,
                  _cplCodeRootRankUnionComm,
                  (int) data_tag,
                  _unionComm,
                  data_recv_request);
      }
    } // end if i_rank is the root rank of the local code

  }


  /**
   *
   * \brief Non blocking send wait of global data array
   *
   * \param [in]  data_send_request
   * \param [in]  s_entity
   * \param [in]  stride
   * \param [in]  n_entity
   * \param [in]  send_data
   * \param [out] recv_data
   */

  void
  Communication::waitIssendGlobalDataBetweenCodesThroughUnionCom
  (
   MPI_Request     data_send_request,
   size_t          s_entity,
   int             stride,
   int             n_entity,
   void           *send_data,
   void           *recv_data
  )
  {
    // Get union communicator ie. union of all active ranks of the codes in the coupling
    int unionCommRank;
    MPI_Comm_rank (_unionComm, &unionCommRank);

    // Get couplings id
    int codeID    = _localCodeProperties->idGet();
    int cplCodeID = _cplCodeProperties->idGet();

    // waitIsend
    if (unionCommRank == _locCodeRootRankUnionComm) {
      if (_cplCodeProperties->localCodeIs()) {
        if (_locCodeRootRankUnionComm == _cplCodeRootRankUnionComm) {
          if (codeID < cplCodeID) {
            assert(send_data != NULL);
            memcpy(recv_data, send_data, s_entity * stride * n_entity);
          } // end code with smallest id does the action of the rank
        } // end if i_rank is root rank of coupled code
        else {
          MPI_Wait(&data_send_request, MPI_STATUS_IGNORE);
        }
      } // end if i_rank is joint with the coupled code
      else {
        MPI_Wait(&data_send_request, MPI_STATUS_IGNORE);
      }
    } // end if i_rank is the root rank of the local code
  }


  /**
   *
   * \brief Non blocking receive wait of global data array
   *
   * \param [in]  data_recv_request
   * \param [in]  s_entity
   * \param [in]  stride
   * \param [in]  n_entity
   * \param [in]  send_data
   * \param [out] recv_data
   */

  void
  Communication::waitIrecvGlobalDataBetweenCodesThroughUnionCom
  (
   MPI_Request     data_recv_request,
   size_t          s_entity,
   int             stride,
   int             n_entity,
   void           *send_data,
   void           *recv_data
  )
  {
    // Get union communicator ie. union of all active ranks of the codes in the coupling
    int unionCommRank;
    MPI_Comm_rank (_unionComm, &unionCommRank);

    // Get couplings id
    int codeID    = _localCodeProperties->idGet();
    int cplCodeID = _cplCodeProperties->idGet();

    // waitIrecv
    if (unionCommRank ==  _locCodeRootRankUnionComm) {
      if (_cplCodeProperties->localCodeIs()) {
        if (_locCodeRootRankUnionComm == _cplCodeRootRankUnionComm) {
          if (codeID < cplCodeID) {
            assert(send_data != NULL);
            memcpy(recv_data, send_data, s_entity * stride * n_entity);
          } // end code with smallest id does the action of the rank
        } // end if i_rank is root rank of coupled code
        else {
          MPI_Wait(&data_recv_request, MPI_STATUS_IGNORE);
        }
      } // end if i_rank is joint with the coupled code
      else {
        MPI_Wait(&data_recv_request, MPI_STATUS_IGNORE);
      }
    } // end if i_rank is the root rank of the local code

    // Broadcast
    int connectable_i_rank;
    MPI_Comm_rank(_localCodeProperties->connectableCommGet(), &connectable_i_rank);
    MPI_Bcast(recv_data,
              (int) s_entity * stride * n_entity,
              MPI_UNSIGNED_CHAR,
              0,
              _localCodeProperties->connectableCommGet());
  }

}

/**
 * \endcond
 */
