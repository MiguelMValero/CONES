#ifndef __CODE_PROPERTIES_H__
#define __CODE_PROPERTIES_H__
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


/** \file codeProperties.hxx
 *
 *  \brief Header file of codeProperties class
 *
 *   This file is the header file of codeProperties class.
 *
 */

#include <unistd.h>

#include <mpi.h>
#include <cassert>
#include <cstring>
#include <iostream>

#include <map>
#include <vector>
#include <string>
#include <typeinfo>

#include "pdm_error.h"

#include "cwp.h"

/**
 * \cond
 */

using namespace std;



  /**
   * \namespace cwipi
   *
   * \brief CWIPI namespace
   *
   * CWIPI namespace
   *
   */

namespace cwipi {

  /**
   * \class CodeProperties
   *        codeProperties.hxx
   *        "codeProperties.hxx"
   *
   * \brief Code properties management.
   *
   *  This class manages code properties :
   *  - Local control parameters,
   *  - Distant control parameters,
   *  - MPI communicators
   *  - .
   *
   */

  class CodeProperties {

    friend class CodePropertiesDB;

  public:

    /**
     * \brief Return if the current rank is a coupled rank
     *
     * \return  isActiveRank
     *
     */

    inline bool
    isActiveRank() const;


    /**
     * \brief Set if the current rank is coupled
     *
     * \param[in] status   Coupled rank status
     *
     */

    inline void
    isActiveRankset
    (
    bool status
    );

    /**
     * \brief Constructor.
     *
     * \param [in]  name         Current code name
     * \param [in]  id           Identifier
     * \param [in]  rootRank     Root rank in global communicator
     * \param [in]  isLocal      Is a local code
     * \param [in]  globalComm   MPI communicator containing all processes
     *                           of all codes
     * \param [in]  n_param_max  Maximum number of parameters
     * \param [in]  str_size_max Maximum string size
     *
     */

    CodeProperties
    (
     string         &name,
     int            id,
     int            rootRank,
     bool            isLocal,
     const MPI_Comm  globalComm,
     int            n_param_max,
     int            str_size_max
    );

    /**
     * \brief Copy constructor.
     *
     * \param [in]  other        other code properties
     *
     */

    CodeProperties
    (
     const CodeProperties& other
    );

    /**
     * \brief Destructor
     *
     */

    virtual ~CodeProperties();

    /**
     * \brief Get code name
     *
     * \return    Name
     *
     */

    inline const string &
    nameGet() const;

    /**
     * \brief Get root rank in global communicator
     *
     * \return    Root rank
     *
     */

    inline int
    rootRankGet() const;

    /**
     * \brief Get identifier
     *
     * \return    Identifier
     *
     */

    inline int
    idGet() const;

    /**
     * \brief Get global MPI communicator
     *
     * \return    Global MPI communicator
     *
     */

    inline const MPI_Comm &
    globalCommGet() const;

    /**
     * \brief Get MPI intra-communicator
     *
     * \return    MPI intra-communicator
     *
     */

    inline const MPI_Comm &
    intraCommGet() const;

    /**
     * \brief Set MPI intra-communicator
     *
     * \param[in] intraComm    MPI intra-communicator
     *
     */

    inline void
    intraCommSet
    (
     MPI_Comm intraComm
    );

    /**
     * \brief Get MPI Group in global communicator
     *
     * \return    MPI group
     *
     */

    inline const MPI_Group &
    groupGet() const;


    /**
     * \brief Set MPI intra-communicator
     *
     * \param[in] group    MPI intra-communicator group
     *
     */

    inline void
    groupSet
    (
     MPI_Group group
    );

    /**
     * \brief Get coupled MPI Group in global communicator
     *
     * \return    MPI group
     *
     */

    inline const MPI_Group &
    connectableGroupGet() const;


    /**
     * \brief Get connectable ranks in global communicator
     *
     * \return    Connectable ranks
     *
     */

    inline const vector <int> *
    connectableRanksGet() const;

    inline const vector <int> *
    intraRanksGet() const;

    /**
     * \brief Get connectable MPI Communicator
     *
     * \return    Connectable communicator
     *
     */

    inline const MPI_Comm &
    connectableCommGet() const;

    /**
     * \brief Get a integer control parameter value
     *
     * \param[in] name   Control parameter name
     *
     */

    inline void
    ctrlParamGet
    (
     const string &name,
     int          *value
    );

    /**
     * \brief Get a double control parameter value
     *
     * \param[in] name   Control parameter name
     *
     */

    inline void
    ctrlParamGet
    (
     const string &name,
     double       *value
    );

    /**
     * \brief Get a string control parameter value
     *
     * \param[in] name   Control parameter name
     *
     */

    inline void
    ctrlParamGet
    (
     const string &name,
     char        **value
    );

    /**
     * \brief Set an integer control parameter value
     *
     * \param[in] name    Control parameter name
     * \param[in] value   Value
     *
     */

    inline void
    ctrlParamSet
    (
     const string &name,
     const int     value
    );

    /**
     * \brief Set a double control parameter value
     *
     * \param[in] name    Control parameter name
     * \param[in] value   Value
     *
     */

    inline void
    ctrlParamSet
    (
     const string &name,
     const double  value
    );


    /**
     * \brief Set a string control parameter value
     *
     * \param[in] name    Control parameter name
     * \param[in] value   Value
     *
     */

    inline void
    ctrlParamSet
    (
     const string &name,
     const char   *value
    );

    /**
     * \brief Add integer control parameter value
     *
     * \param[in] name    Control parameter name
     * \param[in] value   Value
     *
     */

    inline void
    ctrlParamAdd
    (
     const string &name,
     const int value
    );

    /**
     * \brief Add a double control parameter value
     *
     * \param[in] name    Control parameter name
     * \param[in] value   Value
     *
     */

    inline void
    ctrlParamAdd
    (
     const string &name,
     const double  value
    );

    /**
     * \brief Add a string control parameter value
     *
     * \param[in] name    Control parameter name
     * \param[in] value   Value
     *
     */

    inline void
    ctrlParamAdd
    (
     const string &name,
     const char   *value
    );

    /**
     * \brief Cancel a control parameter value
     *
     * \param[in] name    Control parameter name
     *
     */

    template<typename T>
    void
    ctrlParamCancel
    (
     const string &name
    );


    /**
     * \brief Return number of parameters
     *
     * \return Number of parameters
     *
     */

    template<typename T>
    int
    ctrlParamNGet
    (
    );


  /**
   * \brief Return list of parameters
   *
   * \param [in]  nParam Pointer to the number of parameter
   * \param [in]  names  Pointer to an array containing the parameter names.
   *
   */

    template<typename T>
    void
    ctrlParamListGet
    (
    int  *nParam,
    char ***names
    );


    /**
     * \brief  Is a parameter ?
     *
     * \param[in] name
     *
     * \return  1 : true / 0 : false
     *
     */

    template<typename T>
    int
    ctrlParamIs
    (
     const string &name
    );


    /**
     * \brief  Is a local code ?
     *
     * \return  1 : true / 0 : false
     *
     */

    inline bool
    localCodeIs
    (
    ) const;


    /**
     * \brief Dump properties
     *
     */

    void
    dump();


    /**
     * \brief Dump string of properties
     *
     */

    string
    str_dump();


    /**
     * \brief Lock access to the control parameters
     *
     */

    inline void
    paramLock();


    /**
     * \brief Is locked param
     *
     */

    inline int
    paramIsLocked();


    /**
     * \brief Unlock access to the control parameters
     *
     */

    inline void
    paramUnLock();

    /**
     * \brief set isLocal code
     *
     * \param[in] status Code is local or not
     */

    inline void
    isLocalSet (bool status);

   /**
     * \brief Set the user structure
     *
     */

    inline void
    userStructureSet (void *userStruct);


   /**
     * \brief Get the user structure
     *
     */

    inline void * 
    userStructureGet();


  private:

    /**
     * \brief Default assignment (unavailable)
     *
     */

    CodeProperties &
    operator=
    (
     const CodeProperties &other
     );


    /**
     * \brief Update integer parameter values
     *
     */

    inline void
    _updateIntValues
    (
    );


    /**
     * \brief Update double parameter values
     *
     *
     */

    inline void
    _updateDoubleValues
    (
    );


    /**
     * \brief Update string parameter values
     *
     */

    inline void
    _updateStrValues
    (
    );



  private:
    string    _name;          /*!< Name */
    int       _id;            /*!< Identifier */
    bool      _isLocal;       /*!< Is a local code */
    int       _rootRankInGlobalComm; /*!< Root rank
                                         *   in MPI global communicator */
    MPI_Comm  _globalComm;    /*!< MPI global communicator */
    MPI_Comm  _intraComm;     /*!< MPI intra communicator */
    bool      _isActiveRank;  /*!< Is a coupled rank */

    void     *_userStruct;    /*!< Generic pointer about a user structure associated to the code */

    MPI_Group _intraGroup;     /*!< MPI group in the global communicator */
    vector <int> *_intraRanks;  /*!< Code ranks in global communicator */
    MPI_Group _intraConnectableGroup; /*!< coupled MPI group in
                                       the global communicator */
    MPI_Comm _intraConnectableComm; /*!< coupled MPI intra communicator */

    vector <int> *_connectableRanks;  /*!< Coupled code ranks in global communicator */

    MPI_Win   _winGlob;        /*!< MPI window to store general parameters informations */
    int       _winGlobData[4]; /*!< \ref _winGlob data (defined only on \ref _rootRankInGlobalComm :
                                   *      - Lock Param Status
                                   *      - Number of int parameters
                                   *      - Number of double parameters
                                   *      - Number of string parameters */

    MPI_Win   _winIntParamIdxName; /*!< Window to store indexes of int param names
                                                 * size = Number of int parameters + 1 */
    MPI_Win   _winIntParamName; /*!< Window to store param names
                                                 * size = \ref _winIntParamIdxName[Number of int parameter] */
    MPI_Win   _winIntParamValue; /*!< Window to store int param values
                                                 * size = Number of int parameters */
    int      *_winIntParamIdxNameData;  /*!< Data of \ref _winIntParamIdxName window */
    char     *_winIntParamNameData; /*!< Data of \ref _winIntParamName window */
    int      *_winIntParamValueData; /*!< Data of \ref _winIntParamValue window */

    MPI_Win   _winDoubleParamIdxName;/*!< Window to store indexes of double param names
                                                   * size = Number of double parameters + 1 */
    MPI_Win   _winDoubleParamName;/*!< Window to store param names
                                                * size = /stck/equemera/workspace/cwipi/cwipi/src/commWithPart.cxx(64): error: expression must have class type\ref _winDoubleParamIdxName[Number of double parameter] */
    MPI_Win   _winDoubleParamValue; /*!< Window to store double param values
                                                  * size = Number of int parameters */
    int      *_winDoubleParamIdxNameData; /*!< Data of \ref _winDoubleParamIdxName window */
    char     *_winDoubleParamNameData; /*!< Data of \ref _winDoubleParamName window */
    double   *_winDoubleParamValueData; /*!< Data of \ref _winDoubleParamValue window */

    MPI_Win   _winStrParamIdxName; /*!< Window to store indexes of string param names
                                                   * size = Number of string parameters + 1 */
    MPI_Win   _winStrParamName; /*!< Window to store param names
                                                * size = \ref _winDoubleParamIdxName[Number of double parameter] */
    MPI_Win   _winStrParamIdxValue; /*!< Window to store indexes of string param values
                                                   * size = Number of string parameters + 1 */
    MPI_Win   _winStrParamValue; /*!< Window to store string param values
                                                  * size = Number of string parameters */
    int      *_winStrParamIdxNameData;  /*!< Data of \ref _winStrParamIdxName window */
    char     *_winStrParamNameData; /*!< Data of \ref _winStrParamName window */
    int      *_winStrParamIdxValueData; /*!< Data of \ref _winStrParamIdxValue window */
    char     *_winStrParamValueData; /*!< Data of \ref _winStrParamValue window */

    int      _n_param_max; /*!< Maximum number of parameters */
    size_t   _str_size_max; /*!< Maximum string size */

  };


  /**
   * \brief set isLocal code
   *
   */

  void
  CodeProperties::isLocalSet
  (
  bool status
  )
  {
    _isLocal = status;
  }


 /**
   * \brief Set the user structure
   *
   */

  void
  CodeProperties::userStructureSet
  (
  void *userStruct
  )
  {
    _userStruct = userStruct;
  }


 /**
   * \brief Get the user structure
   *
   */

  void * 
  CodeProperties::userStructureGet
  (
  )
  {
    return _userStruct;
  }


  /**
   * \brief Lock access to the control parameters
   *
   */

  void
  CodeProperties::paramLock()
  {
    _winGlobData[0] = 1;
  }

  /**
   * \brief Is locked param
   *
   */

  int
  CodeProperties::paramIsLocked()
  {
    return (_winGlobData[0] == 1);
  }


  /**
   * \brief Unlock access to the control parameters
   *
   */

  void
  CodeProperties::paramUnLock()
  {
    _winGlobData[0] = 0;
  }


  /**
   * \brief Return if the current rank is a coupled rank
   *
   * \return  isActiveRank
   *
   */

  bool
  CodeProperties::isActiveRank() const
  {
    return _isActiveRank;
  }


  /**
   * \brief Set if the current rank is coupled
   *
   * \param[in] status   Coupled rank status
   *
   */

  void
  CodeProperties::isActiveRankset
  (
  bool status
  )
  {
    _isActiveRank = status;
  }


  /**
   * \brief Get code name
   *
   * \return    Name
   *
   */

  const string &
  CodeProperties::nameGet() const
  {
    return _name;
  }

  /**
   * \brief Get root rank in global communicator
   *
   * \return    Root rank
   *
   */

  int
  CodeProperties::rootRankGet() const
  {
    return _rootRankInGlobalComm;
  }

  /**
   * \brief Get connectable ranks in global communicator
   *
   * \return    Connectable ranks
   *
   */

  const vector <int> *
  CodeProperties::connectableRanksGet() const
  {
    return _connectableRanks;
  }

  const vector <int> *
  CodeProperties::intraRanksGet() const
  {
    return _intraRanks;
  }


  /**
   * \brief Get identifier
   *
   * \return    Identifier
   *
   */

  int
  CodeProperties::idGet() const
  {
    return _id;
  }

  /**
   * \brief Get global MPI communicator
   *
   * \return    Global MPI communicator
   *
   */

  const MPI_Comm &
  CodeProperties::globalCommGet() const
  {
    return _globalComm;
  }

  /**
   * \brief Get MPI intra-communicator
   *
   * \return    MPI intra-communicator
   *
   */

  const MPI_Comm &
  CodeProperties::intraCommGet() const
  {
    return _intraComm;
  }

  /**
   * \brief Set MPI intra-communicator
   *
   * \param[in] intraComm    MPI intra-communicator
   *
   */

  void
  CodeProperties::intraCommSet
  (
   MPI_Comm intraComm
  )
  {
    _intraComm = intraComm;
  }

  /**
   * \brief Get MPI Group in global communicator
   *
   * \return    MPI group
   *
   */

  const MPI_Group &
  CodeProperties::groupGet() const
  {
    return _intraGroup;
  }

  /**
   * \brief Get connectable MPI Group in global communicator
   *
   * \return    MPI group
   *
   */

  const MPI_Group &
  CodeProperties::connectableGroupGet() const
  {
    return _intraConnectableGroup;
  }

  /**
   * \brief Get connectable MPI communicator in global communicator
   *
   * \return    MPI group
   *
   */

  const MPI_Comm &
  CodeProperties::connectableCommGet() const
  {
    return _intraConnectableComm;
  }

  /**
   * \brief Set MPI intra-communicator
   *
   * \param[in] group    MPI intra-communicator group
   *
   */

  void
  CodeProperties::groupSet
  (
   MPI_Group group
  )
  {
    _intraGroup = group;
  }


  /**
   * \brief Update integer parameter values
   *
   *
   */

  void
  CodeProperties::_updateIntValues
  (
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (rank != _rootRankInGlobalComm) {

      int lockStatus = 1;
      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm, 0, 4,
                  MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          usleep(1); // WARNING mandatory to allow ECLUSIVE to access
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus);

      int nIntParam    = _winGlobData[1];

      if (nIntParam > 0) {
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winIntParamIdxName);
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winIntParamName);
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winIntParamValue);

        MPI_Get (_winIntParamIdxNameData, nIntParam + 1,
                 MPI_INT, _rootRankInGlobalComm, 0, nIntParam + 1,
                 MPI_INT, _winIntParamIdxName);
        MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamIdxName);

        MPI_Get (_winIntParamNameData, _winIntParamIdxNameData[nIntParam],
                 MPI_CHAR, _rootRankInGlobalComm, 0, _winIntParamIdxNameData[nIntParam],
                 MPI_CHAR, _winIntParamName);
        MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamName);

        MPI_Get (_winIntParamValueData, nIntParam,
                 MPI_INT, _rootRankInGlobalComm, 0, nIntParam,
                 MPI_INT, _winIntParamValue);
        MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamValue);
      }
    }

    else {
      int lockStatus = _winGlobData[0];
      if (lockStatus) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");
      }
    }
  }


  /**
   * \brief Get an integer control parameter value
   *
   * \param[in]  name    Control parameter name
   * \param[out] value   Control parameter value
   *
   */

  void
  CodeProperties::ctrlParamGet
  (
   const string &name,
   int          *value
  ) 
  {
    // cout << "ctrlParamGet int " << _name  << " " << name << endl;

    // MPI_Barrier (_intraComm);

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    _updateIntValues ();

    int nIntParam = _winGlobData[1];

    int sName = name.size();
    int found = 0;
    int i;
    for (i = 0; i < nIntParam; i++) {
      int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
      if (sName == sParam) {
       found = !strncmp(name.c_str(),
               _winIntParamNameData + _winIntParamIdxNameData[i],
               sName);
      }
      if (found) break;
    }

    if (!found) {
      PDM_error(__FILE__, __LINE__, 0,
                 "'%s' Unknown parameter on '%s' code\n", name.c_str(),
                                                         _name.c_str());
    }

    *value = _winIntParamValueData[i];

    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);

  }


  /**
   * \brief Update double parameter values
   *
   */

  void
  CodeProperties::_updateDoubleValues
  (
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (rank != _rootRankInGlobalComm) {

      int lockStatus = 1;
      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm, 0, 4,
                  MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          usleep(1); // WARNING mandatory to allow ECLUSIVE to access
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus);

      int nDoubleParam    = _winGlobData[2];

      if (nDoubleParam > 0) {
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winDoubleParamIdxName);
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winDoubleParamName);
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winDoubleParamValue);

        MPI_Get (_winDoubleParamIdxNameData, nDoubleParam + 1,
                 MPI_INT, _rootRankInGlobalComm, 0, nDoubleParam + 1,
                 MPI_INT, _winDoubleParamIdxName);
        MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamIdxName);

        MPI_Get (_winDoubleParamNameData, _winDoubleParamIdxNameData[nDoubleParam],
                 MPI_CHAR, _rootRankInGlobalComm, 0, _winDoubleParamIdxNameData[nDoubleParam],
                 MPI_CHAR, _winDoubleParamName);
        MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamName);

        MPI_Get (_winDoubleParamValueData, nDoubleParam,
                 MPI_DOUBLE, _rootRankInGlobalComm, 0, nDoubleParam,
                 MPI_DOUBLE, _winDoubleParamValue);
        MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamValue);
      }
    }

    else {
      int lockStatus = _winGlobData[0];
      if (lockStatus) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");
      }
    }
  }


  /**
   * \brief Get a double control parameter value
   *
   * \param[in]  name    Control parameter name
   * \param[out] value   Pointer to the control parameter value
   *
   */

  void
  CodeProperties::ctrlParamGet
  (
   const string &name,
   double       *value
  ) 
  {
    // cout << "ctrlParamGet double " << _name  << " " << name << endl;

    // MPI_Barrier (_intraComm);

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    _updateDoubleValues ();

    int nDoubleParam    = _winGlobData[2];

    int sName = name.size();
    int found = 0;
    int i;
    for (i = 0; i < nDoubleParam; i++) {
      int sParam = _winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i];
      if (sName == sParam) {
       found = !strncmp(name.c_str(),
               _winDoubleParamNameData + _winDoubleParamIdxNameData[i],
               sName);
      }
      if (found) break;
    }

    if (!found) {
      PDM_error(__FILE__, __LINE__, 0,
                 "'%s' Unknown parameter on '%s' code\n", name.c_str(),
                                                         _name.c_str());
    }

    *value = _winDoubleParamValueData[i];

    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);
  }


  /**
   * \brief Update string parameter values
   *
   */

  void
  CodeProperties::_updateStrValues
  (
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (rank != _rootRankInGlobalComm) {

      int lockStatus = 1;
      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm, 0, 4,
                  MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          usleep(1); // WARNING mandatory to allow ECLUSIVE to access
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus);

      int nStrParam    = _winGlobData[3];

      if (nStrParam > 0) {

        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winStrParamIdxName);
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winStrParamName);
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winStrParamIdxValue);
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winStrParamValue);

        MPI_Get (_winStrParamIdxNameData, nStrParam + 1,
                 MPI_INT, _rootRankInGlobalComm, 0, nStrParam + 1,
                 MPI_INT, _winStrParamIdxName);
        MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxName);

        MPI_Get (_winStrParamNameData, _winStrParamIdxNameData[nStrParam],
                 MPI_CHAR, _rootRankInGlobalComm, 0, _winStrParamIdxNameData[nStrParam],
                 MPI_CHAR, _winStrParamName);
        MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamName);

        MPI_Get (_winStrParamIdxValueData, nStrParam + 1,
                 MPI_INT, _rootRankInGlobalComm, 0, nStrParam + 1,
                 MPI_INT, _winStrParamIdxValue);
        MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxValue);

        MPI_Get (_winStrParamValueData, _winStrParamIdxValueData[nStrParam],
                 MPI_CHAR, _rootRankInGlobalComm, 0, _winStrParamIdxValueData[nStrParam],
                 MPI_CHAR, _winStrParamValue);
        MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamValue);
      }
    }

    else {
      int lockStatus = _winGlobData[0];
      if (lockStatus) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");
      }
    }
  }

  /**
   * \brief Get a string control parameter value
   *
   * \param[in]  name   Control parameter name
   * \param[out] value  Pointer to the control parameter value
   *
   */

  void
  CodeProperties::ctrlParamGet
  (
   const string &name,
   char        **value
  )
  {
    // cout << "ctrlParamGet char " << _name  << " " << name << endl;

    // MPI_Barrier (_intraComm);

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    _updateStrValues();

    int nStrParam    = _winGlobData[3];

    int sName = name.size();
    int found = 0;
    int i;
    for (i = 0; i < nStrParam; i++) {
      int sParam = _winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i];
      if (sName == sParam) {
       found = !strncmp(name.c_str(),
               _winStrParamNameData + _winStrParamIdxNameData[i],
               sName);
      }
      if (found) break;
    }

    if (!found) {
      PDM_error(__FILE__, __LINE__, 0,
                 "'%s' Unknown parameter on '%s' code\n", name.c_str(),
                                                         _name.c_str());
    }

    int sValue = _winStrParamIdxValueData[i+1] - _winStrParamIdxValueData[i];
    *value = (char * ) malloc(sizeof(char) * (sValue + 1));
    strncpy (*value,
             _winStrParamValueData + _winStrParamIdxValueData[i],
             sValue);
    (*value)[sValue] = '\0';

    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);
  }

  /**
   * \brief Set an integer control parameter value
   *
   * \param[in] name    Control parameter name
   * \param[in] value   Value
   *
   */

  void
  CodeProperties::ctrlParamSet
  (
   const string &name,
   const int     value
  )
  {
    if (!_isLocal) {
      PDM_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Set a distant code parameter is not allowed\n",
                                                   _name.c_str());
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winIntParamValue);

      int sName = name.size();
      int found = 0;
      int i;

      int nIntParam    = _winGlobData[1];

      for (i = 0; i < nIntParam; i++) {
        int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
        if (sName == sParam) {
         found = !strncmp(name.c_str(),
                 _winIntParamNameData + _winIntParamIdxNameData[i],
                 sName);
        }
        if (found) break;
      }

      if (!found) {
        PDM_error(__FILE__, __LINE__, 0,
                   "'%s' Unknown parameter on '%s' code\n", name.c_str(),
                                                           _name.c_str());
      }

      _winIntParamValueData[i] = value;

      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
    MPI_Barrier(_intraConnectableComm);
  }

  /**
   * \brief Set a double control parameter value
   *
   * \param[in] name    Control parameter name
   * \param[in] value   Value
   *
   */

  void
  CodeProperties::ctrlParamSet
  (
   const string &name,
   const double  value
  )
  {
    if (!_isLocal) {
      PDM_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Set a distant code parameter is not allowed\n",
                                                   _name.c_str());
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winDoubleParamValue);

      int sName = name.size();
      int found = 0;
      int i;

      int nDoubleParam = _winGlobData[2];

      for (i = 0; i < nDoubleParam; i++) {
        int sParam = _winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i];
        if (sName == sParam) {
         found = !strncmp(name.c_str(),
                 _winDoubleParamNameData + _winDoubleParamIdxNameData[i],
                 sName);
        }
        if (found) break;
      }

      if (!found) {
        PDM_error(__FILE__, __LINE__, 0,
                   "'%s' Unknown parameter on '%s' code\n", name.c_str(),
                                                           _name.c_str());
      }

      _winDoubleParamValueData[i] = value;

      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
    MPI_Barrier(_intraConnectableComm);
  }

  /**
   * \brief Set a string control parameter value
   *
   * \param[in] name    Control parameter name
   * \param[in] value   Value
   *
   */

  void
  CodeProperties::ctrlParamSet
  (
   const string &name,
   const char  *value
  )
  {
    if (!_isLocal) {
      PDM_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Set a distant code parameter is not allowed\n",
                                                   _name.c_str());
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamIdxValue);

      int sName = name.size();
      int found = 0;
      int i;

      int nStrParam    = _winGlobData[3];

      for (i = 0; i < nStrParam; i++) {
        int sParam = _winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i];
        if (sName == sParam) {
         found = !strncmp(name.c_str(),
                 _winStrParamNameData + _winStrParamIdxNameData[i],
                 sName);
        }
        if (found) break;
      }

      if (!found) {
        PDM_error(__FILE__, __LINE__, 0,
                   "'%s' Unknown parameter on '%s' code\n", name.c_str(),
                                                           _name.c_str());
      }

      int sValue = _winStrParamIdxValueData[i+1] - _winStrParamIdxValueData[i];
      int gap = strlen(value) - sValue;

      if (gap != 0) {
        if (gap > 0) {
          for (int i1 = _winStrParamIdxValueData[nStrParam] - 1; i1 >= _winStrParamIdxValueData[i+1]; i1--) {
            _winStrParamValueData[i1+gap] = _winStrParamValueData[i1];
          }
        }
        else {
          for (int i1 = _winStrParamIdxValueData[i+1]; i1 < _winStrParamIdxValueData[nStrParam]; i1++) {
            _winStrParamValueData[i1+gap] = _winStrParamValueData[i1];
          }
        }
        for (int i1 = i+1; i1 < nStrParam; i1++) {
          _winStrParamIdxValueData[i1] += gap;
        }
      }


      strncpy(_winStrParamValueData + _winStrParamIdxValueData[i],
              value, strlen(value));

      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
    MPI_Barrier(_intraConnectableComm);
  }

  /**
   * \brief Add an integer control parameter value
   *
   * \param[in] name    Control parameter name
   * \param[in] value   Value
   *
   */

  void
  CodeProperties::ctrlParamAdd
  (
   const string &name,
   const int     value
  )
  {

    if (!_isLocal) {
      PDM_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Add a distant code parameter is not allowed\n",
                                                   _name.c_str());
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winIntParamValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winIntParamName);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winIntParamIdxName);

      int sName = name.size();
      int found = 0;
      int i;

      int nIntParam    = _winGlobData[1];

      if (nIntParam >= _n_param_max) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Impossible to create the '%s' parameter. \n"
                   "The maximum number of parameters is exceeded. \n"
                   "To increase the maximum number of parameters : "
                   "Define the environment variable 'CWP_N_PARAM_MAX'\n",
                    name.c_str());
      }

      if (name.size() >= _str_size_max)  {
        PDM_error(__FILE__, __LINE__, 0,
                   "Impossible to create the '%s' parameter. \n"
                   "The maximum name size is exceeded. \n"
                   "To allow longer parameter names : "
                   "Define the environment variable 'CWP_STR_SIZE_MAX'\n",
                    name.c_str());
      }

      for (i = 0; i < nIntParam; i++) {
        int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
        if (sName == sParam) {
         found = !strncmp(name.c_str(),
                 _winIntParamNameData + _winIntParamIdxNameData[i],
                 sName);
        }
        if (found) break;
      }

      if (found) {
        PDM_error(__FILE__, __LINE__, 0,
                   "'%s' is already a parameter off '%s' code\n", name.c_str(),
                                                                 _name.c_str());
      }

      nIntParam += 1;

      _winIntParamIdxNameData[nIntParam] = _winIntParamIdxNameData[nIntParam-1] + name.size();

      strncpy(_winIntParamNameData + _winIntParamIdxNameData[nIntParam-1],
              name.c_str(), name.size());

      _winIntParamValueData[nIntParam-1] = value;

      _winGlobData[1] = nIntParam;

      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamIdxName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winIntParamName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
    MPI_Barrier(_intraConnectableComm);
  }


  /**
   * \brief Add a double control parameter value
   *
   * \param[in] name    Control parameter name
   * \param[in] value   Value
   *
   */

  void
  CodeProperties::ctrlParamAdd
  (
   const string &name,
   const double  value
  )
  {
    if (!_isLocal) {
      PDM_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Add a distant code parameter is not allowed\n",
                                                   _name.c_str());
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winDoubleParamValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winDoubleParamName);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winDoubleParamIdxName);

      int sName = name.size();
      int found = 0;
      int i;

      int nDoubleParam = _winGlobData[2];

      if (nDoubleParam >= _n_param_max) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Impossible to create the '%s' parameter. \n"
                   "The maximum number of parameters is exceeded. \n"
                   "To increase the maximum number of parameters : "
                   "Define the environment variable 'CWP_N_PARAM_MAX'\n",
                    name.c_str());
      }

      if (name.size() >= _str_size_max)  {
        PDM_error(__FILE__, __LINE__, 0,
                   "Impossible to create the '%s' parameter. \n"
                   "The maximum name size is exceeded. \n"
                   "To allow longer parameter names : "
                   "Define the environment variable 'CWP_STR_SIZE_MAX'\n",
                    name.c_str());
      }

      for (i = 0; i < nDoubleParam; i++) {
        int sParam = _winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i];
        if (sName == sParam) {
         found = !strncmp(name.c_str(),
                 _winDoubleParamNameData + _winDoubleParamIdxNameData[i],
                 sName);
        }
        if (found) break;
      }

      if (found) {
        PDM_error(__FILE__, __LINE__, 0,
                   "'%s' is already a parameter off '%s' code\n", name.c_str(),
                                                                 _name.c_str());
      }

      nDoubleParam += 1;

      _winDoubleParamIdxNameData[nDoubleParam] = _winDoubleParamIdxNameData[nDoubleParam-1] + name.size();

      strncpy(_winDoubleParamNameData + _winDoubleParamIdxNameData[nDoubleParam-1],
              name.c_str(), name.size());

      _winDoubleParamValueData[nDoubleParam-1] = value;

      _winGlobData[2] = nDoubleParam;

      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamIdxName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winDoubleParamName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }
    MPI_Barrier(_intraConnectableComm);
  }

  /**
   * \brief Add a string control parameter value
   *
   * \param[in] name    Control parameter name
   * \param[in] value   Value
   *
   */

  void
  CodeProperties::ctrlParamAdd
  (
   const string &name,
   const char   *value
  )
  {

    if (!_isLocal) {
      PDM_error(__FILE__, __LINE__, 0,
           "'%s' is a distant code. Add a distant code parameter is not allowed\n",
                                                   _name.c_str());
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamIdxValue);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamName);
      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winStrParamIdxName);

      int sName = name.size();
      int found = 0;
      int i;

      int nStrParam    = _winGlobData[3];

      if (nStrParam >= _n_param_max) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Impossible to create the '%s' parameter. \n"
                   "The maximum number of parameters is exceeded. \n"
                   "To increase the maximum number of parameters : "
                   "Define the environment variable 'CWP_N_PARAM_MAX'\n",
                    name.c_str());
      }

      if (name.size() >= _str_size_max)  {
        PDM_error(__FILE__, __LINE__, 0,
                   "Impossible to create the '%s' parameter. \n"
                   "The maximum name size is exceeded. \n"
                   "To allow longer parameter names : "
                   "Define the environment variable 'CWP_STR_SIZE_MAX'\n",
                    name.c_str());
      }

      if (strlen(value) >= _str_size_max)  {
        PDM_error(__FILE__, __LINE__, 0,
                   "Impossible to create the string '%s' parameter. \n"
                   "The maximum string size is exceeded. \n"
                   "To allow longer parameters : "
                   "Define the environment variable 'CWP_STR_SIZE_MAX'\n",
                    name.c_str());
      }

      for (i = 0; i < nStrParam; i++) {
        int sParam = _winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i];
        if (sName == sParam) {
         found = !strncmp(name.c_str(),
                 _winStrParamNameData + _winStrParamIdxNameData[i],
                 sName);
        }
        if (found) break;
      }

      if (found) {
        PDM_error(__FILE__, __LINE__, 0,
                   "'%s' is already a parameter off '%s' code\n", name.c_str(),
                                                                 _name.c_str());
      }

      nStrParam += 1;

      _winStrParamIdxValueData[nStrParam] = _winStrParamIdxValueData[nStrParam-1] + strlen(value);

      _winStrParamIdxNameData[nStrParam] = _winStrParamIdxNameData[nStrParam-1] + name.size();

      strncpy(_winStrParamNameData + _winStrParamIdxNameData[nStrParam-1],
              name.c_str(), name.size());

      strncpy(_winStrParamValueData + _winStrParamIdxValueData[nStrParam-1],
              value, strlen(value));

      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamName);
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamIdxValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winStrParamValue);
      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);

      _winGlobData[3] = nStrParam;

    }

    MPI_Barrier(_intraConnectableComm);

  }

  /**
   * \brief Cancel a control parameter value
   *
   * \param[in] name    Control parameter name
   *
   */

  template<typename T>
  void
  CodeProperties::ctrlParamCancel
  (
   const string &name
  )
  {
    if (!_isLocal) {
      PDM_error(__FILE__, __LINE__, 0,
                 "'%s' is a distant code. Delete a distant"
                 " code parameter is not allowed\n",
                                                   _name.c_str());
    }

    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    if (_rootRankInGlobalComm == rank) {

      MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0, _winGlob);

      int nIntParam    = _winGlobData[1];
      int nDoubleParam = _winGlobData[2];
      int nStrParam    = _winGlobData[3];

      MPI_Win  *winTypeParamIdxValue = NULL;
      MPI_Win  *winTypeParamValue = NULL;
      MPI_Win  *winTypeParamIdxName = NULL;
      MPI_Win  *winTypeParamName = NULL;

      int nTypeParam = -1;
      int  *winTypeParamIdxValueData = NULL;
      T *winTypeParamValueData = NULL;
      int  *winTypeParamIdxNameData = NULL;
      char *winTypeParamNameData = NULL;

      if (typeid(T) == typeid(string)) {
        nTypeParam               = nStrParam;
        winTypeParamIdxValue     = &_winStrParamIdxValue;
        winTypeParamValue        = &_winStrParamValue;
        winTypeParamIdxName      = &_winStrParamIdxName;
        winTypeParamName         = &_winStrParamName;
        winTypeParamIdxValueData = _winStrParamIdxValueData;
        winTypeParamValueData    = (T *) _winStrParamValueData;
        winTypeParamIdxNameData  = _winStrParamIdxNameData;
        winTypeParamNameData     = _winStrParamNameData;
      }
      else if (typeid(T) == typeid(int)) {
        nTypeParam              = nIntParam;
        winTypeParamValue       = &_winIntParamValue;
        winTypeParamIdxName     = &_winIntParamIdxName;
        winTypeParamName        = &_winIntParamName;
        winTypeParamValueData   = (T *) _winIntParamValueData;
        winTypeParamIdxNameData = _winIntParamIdxNameData;
        winTypeParamNameData    = _winIntParamNameData;
      }
      else if (typeid(T) == typeid(double)) {
        nTypeParam              = nDoubleParam;
        winTypeParamValue       = &_winDoubleParamValue;
        winTypeParamIdxName     = &_winDoubleParamIdxName;
        winTypeParamName        = &_winDoubleParamName;
        winTypeParamValueData   = (T *) _winDoubleParamValueData;
        winTypeParamIdxNameData = _winDoubleParamIdxNameData;
        winTypeParamNameData    = _winDoubleParamNameData;
      }
      else {
        PDM_error(__FILE__, __LINE__, 0,
                  "Type not taken into account \n");
      }

      int i;
      int found;

      for (i = 0; i < nTypeParam; i++) {
        size_t sParam = winTypeParamIdxNameData[i+1] - winTypeParamIdxNameData[i];
        if (name.size() == sParam) {
         found = !strncmp(name.c_str(),
                 winTypeParamNameData + winTypeParamIdxNameData[i],
                 name.size());
        }
        if (found) break;
      }

      if (found) {
        if (winTypeParamIdxValueData != NULL) {
          MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0,
                        *winTypeParamIdxValue);
        }
        MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0,
                      *winTypeParamValue);
        MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0,
                      *winTypeParamName);
        MPI_Win_lock (MPI_LOCK_EXCLUSIVE, _rootRankInGlobalComm, 0,
                      *winTypeParamIdxName);

        if (winTypeParamIdxValueData != NULL) {
          int gap = winTypeParamIdxValueData[i+1] - winTypeParamIdxValueData[i];
          for (int i1 = winTypeParamIdxValueData[i]; i1 < winTypeParamIdxValueData[nTypeParam] - gap; i1++) {
            winTypeParamValueData[i1] = winTypeParamValueData[i1+gap];
          }
          for (int i1 = i; i1 < nTypeParam - 1; i1++) {
            winTypeParamIdxValueData[i1] = winTypeParamIdxValueData[i1+1] - gap;
          }
        }
        else {
          for (int i1 = i; i1 < nTypeParam - 1; i1++) {
            winTypeParamValueData[i1] = winTypeParamValueData[i1+1];
          }
        }

        int gap2 = winTypeParamIdxNameData[i+1] - winTypeParamIdxNameData[i];
        for (int i1 = winTypeParamIdxNameData[i]; i1 < winTypeParamIdxNameData[nTypeParam] - gap2; i1++) {
          winTypeParamNameData[i1] = winTypeParamNameData[i1+gap2];
        }

        for (int i1 = i; i1 < nTypeParam - 1; i1++) {
          winTypeParamIdxNameData[i1] = winTypeParamIdxNameData[i1+1] - gap2;
        }

        nTypeParam += -1;

        winTypeParamIdxNameData[nStrParam+1] = winTypeParamIdxNameData[nStrParam] + name.size();

        if (winTypeParamIdxValueData != NULL) {
          MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxValue);
        }
        MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamValue);
        MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamName);
        MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxName);

      }

      MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
    }

  }

  /**
   * \brief Return number of parameters
   *
   * \return Number of parameters
   *
   */

  template<typename T>
  int
  CodeProperties::ctrlParamNGet
  (
   )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    int oldLockStatus   = _winGlobData[0];
    int oldNIntParam    = _winGlobData[1];
    int oldNDoubleParam = _winGlobData[2];
    int oldNStrParam    = _winGlobData[3];
    int lockStatus      = oldLockStatus;
    int nIntParam       = oldNIntParam;
    int nDoubleParam    = oldNDoubleParam;
    int nStrParam       = oldNStrParam;

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    if (rank != _rootRankInGlobalComm) {

      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm,
                  0, 4, MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          usleep(1); // WARNING mandatory to allow ECLUSIVE to access
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus);

      nIntParam    = _winGlobData[1];
      nDoubleParam = _winGlobData[2];
      nStrParam    = _winGlobData[3];
    }

    else {
      if (lockStatus) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");
      }

    }

    MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);

    int nParam = 0;
    if (typeid(T) == typeid(string)) {
      nParam = nStrParam;
    }
    else if (typeid(T) == typeid(int)) {
      nParam = nIntParam;
    }
    else if (typeid(T) == typeid(double)) {
      nParam = nDoubleParam;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0,
                "Type not taken into account \n");
    }

    return nParam;
  }


  /**
   * \brief Return list of parameters
   *
   * \param [in]  nParam Pointer to the number of parameter
   * \param [in]  names  Pointer to an array containing the parameter names.
   *
   */

  template<typename T>
  void
  CodeProperties::ctrlParamListGet
  (
   int  *nParam,
   char ***names
  )
  {
    int rank;
    MPI_Comm_rank(_globalComm, &rank);

    int oldLockStatus   = _winGlobData[0];
    int oldNIntParam    = _winGlobData[1];
    int oldNDoubleParam = _winGlobData[2];
    int oldNStrParam    = _winGlobData[3];
    int lockStatus      = oldLockStatus;
    int nIntParam       = oldNIntParam;
    int nDoubleParam    = oldNDoubleParam;
    int nStrParam       = oldNStrParam;

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    if (typeid(T) == typeid(string)) {
      _updateStrValues();
    }
    else if (typeid(T) == typeid(int)) {
      _updateIntValues ();
    }
    else if (typeid(T) == typeid(double)) {
      _updateDoubleValues();
    }

    if (rank != _rootRankInGlobalComm) {

      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm,
                  0, 4, MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          usleep(1); // WARNING mandatory to allow ECLUSIVE to access
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus);

      nIntParam    = _winGlobData[1];
      nDoubleParam = _winGlobData[2];
      nStrParam    = _winGlobData[3];

    }

    else {
      if (lockStatus) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");
      }

    }

    MPI_Win  *winTypeParamIdxName = NULL;
    MPI_Win  *winTypeParamName = NULL;

    int nTypeParam = 0;
    int  *winTypeParamIdxNameData = NULL;
    char *winTypeParamNameData = NULL;

    if (typeid(T) == typeid(string)) {
      nTypeParam               = nStrParam;
      winTypeParamIdxName      = &_winStrParamIdxName;
      winTypeParamName         = &_winStrParamName;
      winTypeParamIdxNameData  = _winStrParamIdxNameData;
      winTypeParamNameData     = _winStrParamNameData;
    }
    else if (typeid(T) == typeid(int)) {
      nTypeParam              = nIntParam;
      winTypeParamIdxName     = &_winIntParamIdxName;
      winTypeParamName        = &_winIntParamName;
      winTypeParamIdxNameData = _winIntParamIdxNameData;
      winTypeParamNameData    = _winIntParamNameData;
    }
    else if (typeid(T) == typeid(double)) {
      nTypeParam              = nDoubleParam;
      winTypeParamIdxName     = &_winDoubleParamIdxName;
      winTypeParamName        = &_winDoubleParamName;
      winTypeParamIdxNameData = _winDoubleParamIdxNameData;
      winTypeParamNameData    = _winDoubleParamNameData;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0,
                "Type not taken into account \n");
    }

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0,
                  *winTypeParamName);
    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0,
                  *winTypeParamIdxName);

    *names = (char **) malloc(sizeof(char *) * nTypeParam);

    for (int i = 0; i < nTypeParam; i++) {
      int sParam = winTypeParamIdxNameData[i+1] - winTypeParamIdxNameData[i];
      (*names)[i] = (char *) malloc(sizeof(char) * (sParam + 1));
      char *curName = (*names)[i];

      strncpy(curName,
              winTypeParamNameData + winTypeParamIdxNameData[i],
              sParam);

      curName[sParam] = '\0';

    }

    *nParam = nTypeParam;

    MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamName);
    MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxName);
    MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);

  }


  /**
   * \brief  Is a parameter ?
   *
   * \param[in] name
   *
   * \return  1 : true / 0 : false
   *
   */

  template<typename T>
  int
  CodeProperties::ctrlParamIs
  (
   const string &name
  )
  {
    // Get MPI rank index
    int i_rank;
    MPI_Comm_rank(_globalComm, &i_rank);

    // Lock global window
    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);

    // Variables to template
    int nTypeParam = 0;

    int  *winTypeParamIdxNameData = NULL;
    char *winTypeParamNameData    = NULL;

    MPI_Win  *winTypeParamIdxName = NULL;
    MPI_Win  *winTypeParamName = NULL;

    // if not rank with window
    if (i_rank != _rootRankInGlobalComm) {

      int lockStatus = 1;
      do {
        MPI_Request rq1;
        MPI_Rget ((void *) _winGlobData, 4, MPI_INT, _rootRankInGlobalComm, 0, 4,
                  MPI_INT, _winGlob, &rq1);
        MPI_Wait (&rq1, MPI_STATUS_IGNORE);
        lockStatus = _winGlobData[0];
        if (lockStatus) {
          MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);
          usleep(1); // WARNING mandatory to allow ECLUSIVE to access
          MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);
        }
      }  while (lockStatus);

      if (typeid(T) == typeid(string)) {

        nTypeParam               = _winGlobData[3];
        winTypeParamIdxNameData  = _winStrParamIdxNameData;
        winTypeParamNameData     = _winStrParamNameData;
        winTypeParamIdxName      = &_winStrParamIdxName;
        winTypeParamName         = &_winStrParamName;

      }
      else if (typeid(T) == typeid(int)) {

        nTypeParam               = _winGlobData[1];
        winTypeParamIdxNameData  = _winIntParamIdxNameData;
        winTypeParamNameData     = _winIntParamNameData;
        winTypeParamIdxName     = &_winIntParamIdxName;
        winTypeParamName        = &_winIntParamName;

      }
      else if (typeid(T) == typeid(double)) {

        nTypeParam               = _winGlobData[2];
        winTypeParamIdxNameData  = _winDoubleParamIdxNameData;
        winTypeParamNameData     = _winDoubleParamNameData;
        winTypeParamIdxName     = &_winDoubleParamIdxName;
        winTypeParamName        = &_winDoubleParamName;

      }
      else {
        PDM_error(__FILE__, __LINE__, 0,
                  "Type not taken into account \n");
      }

      if (nTypeParam > 0) {
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0,
                      *winTypeParamName);
        MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0,
                      *winTypeParamIdxName);

        MPI_Get (winTypeParamIdxNameData, nTypeParam + 1,
                 MPI_INT, _rootRankInGlobalComm, 0, nTypeParam + 1,
                 MPI_INT, *winTypeParamIdxName);

        MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxName);

        MPI_Get (winTypeParamNameData, winTypeParamIdxNameData[nTypeParam],
                 MPI_CHAR, _rootRankInGlobalComm, 0, winTypeParamIdxNameData[nTypeParam],
                 MPI_CHAR, *winTypeParamName);

        MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamName);
      }
    }

    else {

      if (typeid(T) == typeid(string)) {

        nTypeParam               = _winGlobData[3];
        winTypeParamIdxNameData  = _winStrParamIdxNameData;
        winTypeParamNameData     = _winStrParamNameData;
        winTypeParamIdxName      = &_winStrParamIdxName;
        winTypeParamName         = &_winStrParamName;

      }
      else if (typeid(T) == typeid(int)) {

        nTypeParam               = _winGlobData[1];
        winTypeParamIdxNameData  = _winIntParamIdxNameData;
        winTypeParamNameData     = _winIntParamNameData;
        winTypeParamIdxName     = &_winIntParamIdxName;
        winTypeParamName        = &_winIntParamName;

      }
      else if (typeid(T) == typeid(double)) {

        nTypeParam               = _winGlobData[2];
        winTypeParamIdxNameData  = _winDoubleParamIdxNameData;
        winTypeParamNameData     = _winDoubleParamNameData;
        winTypeParamIdxName     = &_winDoubleParamIdxName;
        winTypeParamName        = &_winDoubleParamName;

      }
      else {
        PDM_error(__FILE__, __LINE__, 0,
                  "Type not taken into account \n");
      }

      int lockStatus = _winGlobData[0];
      if (lockStatus) {
        PDM_error(__FILE__, __LINE__, 0,
                   "Unlock parameters before read its on the current rank\n");
      }
    }

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0,
                  *winTypeParamName);
    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0,
                  *winTypeParamIdxName);

    // check correct name if there are any parameters
    int found = 0;
    if (nTypeParam > 0) {
      int sName = name.size();
      for (int i = 0; i < nTypeParam; i++) {
        int sParam = winTypeParamIdxNameData[i+1] - winTypeParamIdxNameData[i];
        if (sName == sParam) {
          found = !strncmp(name.c_str(),
                           winTypeParamNameData + winTypeParamIdxNameData[i],
                           sName);
        }
        if (found) break;
      }
    }
    MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamName);
    MPI_Win_unlock (_rootRankInGlobalComm, *winTypeParamIdxName);

    MPI_Win_unlock (_rootRankInGlobalComm, _winGlob);

    return found;
  }

  /**
   * \brief  Is a local code ?
   *
   * \return  1 : true / 0 : false
   *
   */

  bool
  CodeProperties::localCodeIs
  (
  ) const
  {
    return _isLocal;
  }

}


/**
 * \endcond
 */

#endif /* __CODE_PROPERTIES_H__ */
