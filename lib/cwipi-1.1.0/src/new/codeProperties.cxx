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
#include <mpi.h>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "codeProperties.hxx"
#include "pdm_printf.h"


/**
 * \cond
 */

using namespace std;

namespace cwipi
{

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

  CodeProperties::CodeProperties
  (
   string &name,
   int    id,
   int    rootRank,
   bool   isLocal,
   const MPI_Comm globalComm,
   int            n_param_max,
   int            str_size_max
  ): _name(name), _id(id), _isLocal(isLocal),
     _rootRankInGlobalComm(rootRank),
     _globalComm(globalComm),
     _isActiveRank(false),
     _userStruct(NULL),
     _winIntParamIdxName(MPI_WIN_NULL),
     _winIntParamName(MPI_WIN_NULL),
     _winIntParamValue(MPI_WIN_NULL),
     _winIntParamIdxNameData(NULL),
     _winIntParamNameData(NULL),
     _winIntParamValueData(NULL),
     _winDoubleParamIdxName(MPI_WIN_NULL),
     _winDoubleParamName(MPI_WIN_NULL),
     _winDoubleParamValue(MPI_WIN_NULL),
     _winDoubleParamIdxNameData(NULL),
     _winDoubleParamNameData(NULL),
     _winDoubleParamValueData(NULL),
     _winStrParamIdxName(MPI_WIN_NULL),
     _winStrParamName(MPI_WIN_NULL),
     _winStrParamIdxValue(MPI_WIN_NULL),
     _winStrParamValue(MPI_WIN_NULL),
     _winStrParamIdxNameData(NULL),
     _winStrParamNameData(NULL),
     _winStrParamIdxValueData(NULL),
     _winStrParamValueData(NULL),
     _n_param_max(n_param_max),
     _str_size_max(str_size_max)
  {
    _intraComm = MPI_COMM_NULL;
    _intraConnectableGroup = MPI_GROUP_NULL;
    _intraConnectableComm = MPI_COMM_NULL;
    _intraGroup        = MPI_GROUP_NULL;
    _intraRanks        = NULL;
    _connectableRanks  = NULL;

    _winGlobData[0] = 0; // Unlock parameters access
    _winGlobData[1] = 0; // 0 int param
    _winGlobData[2] = 0; // 0 doube param
    _winGlobData[3] = 0; // 0 str param

  }

  /**
   * \brief Copy constructor.
   *
   * \param [in]  other        other code properties
   *
   */

  CodeProperties::CodeProperties
  (
   const CodeProperties& other
  ): _name(other._name), _id(other._id), _isLocal(other._isLocal),
     _rootRankInGlobalComm(other._rootRankInGlobalComm),
     _globalComm(other._globalComm),
     _intraComm(other._intraComm),
     _isActiveRank(other._isActiveRank),
     _intraGroup(other._intraGroup),
     _intraRanks(other._intraRanks),
     _intraConnectableGroup(other._intraConnectableGroup),
     _intraConnectableComm(other._intraConnectableComm),
     _connectableRanks(other._connectableRanks),
     _winIntParamIdxName(other._winIntParamIdxName),
     _winIntParamName(other._winIntParamName),
     _winIntParamValue(other._winIntParamValue),
     _winIntParamIdxNameData(other._winIntParamIdxNameData),
     _winIntParamNameData(other._winIntParamNameData),
     _winIntParamValueData(other._winIntParamValueData),
     _winDoubleParamIdxName(other._winDoubleParamIdxName),
     _winDoubleParamName(other._winDoubleParamName),
     _winDoubleParamValue(other._winDoubleParamValue),
     _winDoubleParamIdxNameData(other._winDoubleParamIdxNameData),
     _winDoubleParamNameData(other._winDoubleParamNameData),
     _winDoubleParamValueData(other._winDoubleParamValueData),
     _winStrParamIdxName(other._winStrParamIdxName),
     _winStrParamName(other._winStrParamName),
     _winStrParamIdxValue(other._winStrParamIdxValue),
     _winStrParamValue(other._winStrParamValue),
     _winStrParamIdxNameData(other._winStrParamIdxNameData),
     _winStrParamNameData(other._winStrParamNameData),
     _winStrParamIdxValueData(other._winStrParamIdxValueData),
     _winStrParamValueData(other._winStrParamValueData),
     _n_param_max(other._n_param_max),
     _str_size_max(other._str_size_max)

  {

    memcpy(_winGlobData, other._winGlobData, 4 *sizeof(int));

  }

  /**
   * \brief Dump properties
   *
   */

  void
  CodeProperties::dump()
  {
    PDM_printf ("'%s' properties\n",_name.c_str());
    PDM_printf ("  - Identifier : %d\n", _id);
    PDM_printf ("  - Root rank in global communicator : %d\n", _rootRankInGlobalComm);
    PDM_printf ("  - Is it a local code : %d\n", _isLocal);
    if (_isLocal) {
      PDM_printf ("  - Is it an active rank? : %d\n", _isActiveRank);
    }
    PDM_printf ("  - Ranks in global communicator :");
    for (size_t i = 0; i < _intraRanks->size(); i++) {
      PDM_printf (" %d", (*_intraRanks)[i]);
    }
    PDM_printf ("\n");

    PDM_printf ("  - Coupled ranks in global communicator :");
    for (size_t i = 0; i < _connectableRanks->size(); i++) {
      PDM_printf (" %d", (*_connectableRanks)[i]);
    }
    PDM_printf ("\n");

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);


    char tmpName [81];
    char tmpValue[81];
    unsigned int sParamMax;

    /* Integer parameters */
    _updateIntValues ();

    PDM_printf ("  - %d integer control parameters \n", _winGlobData[1]);

    sParamMax = 0;
    for (int i = 0; i < _winGlobData[1]; i++) {
      unsigned int sParam = (unsigned int) (_winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i]);
      sParamMax = max(sParam, sParamMax);
    }

    if (sParamMax > 80) sParamMax = 80;

    char fmtIntName[22];
    sprintf(fmtIntName, "     * %%%u.%us : %%d\n",sParamMax, sParamMax);

    for (int i = 0; i < _winGlobData[1]; i++) {
      int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];

      strncpy(tmpName,
              _winIntParamNameData + _winIntParamIdxNameData[i],
              min (sParam, (int) sParamMax));
      tmpName[sParam] = '\0';

      PDM_printf (fmtIntName, tmpName, _winIntParamValueData[i]);
    }




    /* Double parameters */
    _updateDoubleValues ();

    PDM_printf ("  - %d double control parameters \n", _winGlobData[2]);

    sParamMax = 0;
    for (int i = 0; i < _winGlobData[2]; i++) {
      unsigned int sParam = (unsigned int) (_winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i]);
      sParamMax = max(sParam, sParamMax);
    }

    if (sParamMax > 80) sParamMax = 80;

    char fmtDoubleName[26];
    sprintf(fmtDoubleName, "     * %%%u.%us : %%12.5e\n",sParamMax, sParamMax);

    for (int i = 0; i < _winGlobData[2]; i++) {
      int sParam = _winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i];

      strncpy(tmpName,
              _winDoubleParamNameData + _winDoubleParamIdxNameData[i],
              min (sParam, (int) sParamMax));
      tmpName[sParam] = '\0';

      PDM_printf (fmtDoubleName, tmpName, _winDoubleParamValueData[i]);
    }



    /* Char parameters */
    _updateStrValues ();

    PDM_printf ("  - %d string control parameters \n", _winGlobData[3]);

    sParamMax = 0;
    unsigned int sValueMax = 0;
    for (int i = 0; i < _winGlobData[3]; i++) {
      unsigned int sParam = (unsigned int) (_winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i]);
      sParamMax = max(sParam, sParamMax);
      unsigned int sValue = (unsigned int) (_winStrParamIdxValueData[i+1] - _winStrParamIdxValueData[i]);
      sValueMax = max(sValue, sValueMax);
    }

    char fmtStrName[27];
    sprintf(fmtStrName, "     * %%%u.%us : %%%u.%us\n",
            sParamMax, sParamMax,
            sValueMax, sValueMax);

    for (int i = 0; i < _winGlobData[3]; i++) {
      int sParam = _winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i];
      strncpy (tmpName,
               _winStrParamNameData + _winStrParamIdxNameData[i],
               min (sParam, (int) sParamMax));
      tmpName[sParam] = '\0';

      int sValue = _winStrParamIdxValueData[i+1] - _winStrParamIdxValueData[i];
      strncpy (tmpValue,
               _winStrParamValueData + _winStrParamIdxValueData[i],
               min (sValue, (int) sValueMax));
      tmpValue[sValue] = '\0';

      PDM_printf (fmtStrName, tmpName, tmpValue);
    }

    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);

    PDM_printf_flush();
  }


  /**
   * \brief Dump string of properties
   *
   */

  string
  CodeProperties::str_dump()
  {
    char buffer[1080];
    sprintf(buffer, "'%s' properties\n",_name.c_str());
    string properties = buffer;
    sprintf(buffer, "  - Identifier : %d\n", _id);
    properties.append(buffer);
    sprintf(buffer, "  - Root rank in global communicator : %d\n", _rootRankInGlobalComm);
    properties.append(buffer);
    sprintf(buffer, "  - Is it a local code : %d\n", _isLocal);
    properties.append(buffer);
    if (_isLocal) {
      sprintf(buffer, "  - Is it an active rank? : %d\n", _isActiveRank);
      properties.append(buffer);
    }
    sprintf(buffer, "  - Ranks in global communicator :");
    properties.append(buffer);
    for (size_t i = 0; i < _intraRanks->size(); i++) {
      sprintf(buffer, " %d", (*_intraRanks)[i]);
      properties.append(buffer);
    }
    sprintf(buffer, "\n");
    properties.append(buffer);

    sprintf(buffer, "  - Coupled ranks in global communicator :");
    properties.append(buffer);
    for (size_t i = 0; i < _connectableRanks->size(); i++) {
      sprintf(buffer, " %d", (*_connectableRanks)[i]);
      properties.append(buffer);
    }
    sprintf(buffer, "\n");
    properties.append(buffer);

    MPI_Win_lock (MPI_LOCK_SHARED, _rootRankInGlobalComm, 0, _winGlob);


    char tmpName [81];
    char tmpValue[81];


    // Update and print integer
    _updateIntValues ();

    sprintf(buffer, "  - %d integer control parameters \n", _winGlobData[1]);
    properties.append(buffer);

    unsigned int sParamMax = 0;
    for (int i = 0; i < _winGlobData[1]; i++) {
      unsigned int sParam = (unsigned int) (_winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i]);
      sParamMax = max(sParam, sParamMax);
    }

    if (sParamMax > 80) sParamMax = 80;

    char fmtIntName[22];
    sprintf(fmtIntName, "     * %%%u.%us : %%d\n",sParamMax, sParamMax);

    for (int i = 0; i < _winGlobData[1]; i++) {
      int sParam = _winIntParamIdxNameData[i+1] - _winIntParamIdxNameData[i];
      strncpy (tmpName,
               _winIntParamNameData + _winIntParamIdxNameData[i],
               min (sParam, (int) sParamMax));
      tmpName[sParam] = '\0';

      sprintf(buffer, fmtIntName, tmpName, _winIntParamValueData[i]);
      properties.append(buffer);

    }

    // Update and print double
    _updateDoubleValues ();

    sParamMax = 0;
    for (int i = 0; i < _winGlobData[2]; i++) {
      unsigned int sParam = (unsigned int) (_winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i]);
      sParamMax = max(sParam, sParamMax);
    }

    if (sParamMax > 80) sParamMax = 80;

    char fmtDoubleName[26];
    sprintf(fmtDoubleName, "     * %%%u.%us : %%12.5e\n",sParamMax, sParamMax);

    sprintf(buffer, "  - %d double control parameters \n", _winGlobData[2]);
    properties.append(buffer);

    for (int i = 0; i < _winGlobData[2]; i++) {
      int sParam = _winDoubleParamIdxNameData[i+1] - _winDoubleParamIdxNameData[i];
      strncpy (tmpName,
               _winDoubleParamNameData + _winDoubleParamIdxNameData[i],
               min (sParam, (int) sParamMax));
      tmpName[sParam] = '\0';

      sprintf(buffer, fmtDoubleName, tmpName, _winDoubleParamValueData[i]);
      properties.append(buffer);

    }

    // Update and print string
    _updateStrValues ();

    sParamMax = 0;
    unsigned int sValueMax = 0;
    for (int i = 0; i < _winGlobData[3]; i++) {
      unsigned int sParam = (unsigned int) (_winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i]);
      sParamMax = max(sParam, sParamMax);
      unsigned int sValue = (unsigned int) (_winStrParamIdxValueData[i+1] - _winStrParamIdxValueData[i]);
      sValueMax = max(sValue, sValueMax);
    }

    if (sParamMax > 80) sParamMax = 80;
    if (sValueMax > 80) sValueMax = 80;

    char fmtStrName[27];
    sprintf(fmtStrName, "     * %%%u.%us : %%%u.%us\n",
            sParamMax, sParamMax,
            sValueMax, sValueMax);

    sprintf(buffer, "  - %d string control parameters \n", _winGlobData[3]);
    properties.append(buffer);

    for (int i = 0; i < _winGlobData[3]; i++) {
      int sParam = _winStrParamIdxNameData[i+1] - _winStrParamIdxNameData[i];
      strncpy (tmpName,
               _winStrParamNameData + _winStrParamIdxNameData[i],
               min (sParam, (int) sParamMax));
      tmpName[sParam] = '\0';

      int sValue = _winStrParamIdxValueData[i+1] - _winStrParamIdxValueData[i];
      strncpy (tmpValue,
               _winStrParamValueData + _winStrParamIdxValueData[i],
               min (sValue, (int) sValueMax));
      tmpValue[sValue] = '\0';

      sprintf(buffer, fmtStrName, tmpName, tmpValue);
      properties.append(buffer);

    }

    MPI_Win_unlock ( _rootRankInGlobalComm, _winGlob);

    return properties;
  }

  /**
   * \brief Destructor
   *
   */

  CodeProperties::~CodeProperties()
  {

    if (_winGlob != MPI_WIN_NULL) {
      MPI_Win_free(&_winGlob);
    }
    if (_winIntParamIdxName != MPI_WIN_NULL) {
      MPI_Win_free(&_winIntParamIdxName);
    }
    if (_winIntParamName != MPI_WIN_NULL) {
      MPI_Win_free(&_winIntParamName);
    }
    if (_winIntParamValue != MPI_WIN_NULL) {
      MPI_Win_free(&_winIntParamValue);
    }
    if (_winDoubleParamIdxName != MPI_WIN_NULL) {
      MPI_Win_free(&_winDoubleParamIdxName);
    }
    if (_winDoubleParamName != MPI_WIN_NULL) {
      MPI_Win_free(&_winDoubleParamName);
    }
    if (_winDoubleParamValue != MPI_WIN_NULL) {
      MPI_Win_free(&_winDoubleParamValue);
    }
    if (_winStrParamIdxName != MPI_WIN_NULL) {
      MPI_Win_free(&_winStrParamIdxName);
    }
    if (_winStrParamName != MPI_WIN_NULL) {
      MPI_Win_free(&_winStrParamName);
    }
    if (_winStrParamIdxValue != MPI_WIN_NULL) {
      MPI_Win_free(&_winStrParamIdxValue);
    }
    if (_winStrParamValue != MPI_WIN_NULL) {
      MPI_Win_free(&_winStrParamValue);
    }

    if (_winIntParamIdxNameData != NULL) {
      free (_winIntParamIdxNameData);
      _winIntParamIdxNameData = NULL;
    }

    if (_winIntParamNameData != NULL) {
      free (_winIntParamNameData);
      _winIntParamNameData = NULL;
    }

    if (_winIntParamValueData != NULL) {
      free (_winIntParamValueData);
      _winIntParamValueData = NULL;
    }

    if (_winDoubleParamIdxNameData != NULL) {
      free (_winDoubleParamIdxNameData);
      _winDoubleParamIdxNameData = NULL;
    }

    if (_winDoubleParamNameData != NULL) {
      free (_winDoubleParamNameData);
      _winDoubleParamNameData = NULL;
    }

    if (_winDoubleParamValueData != NULL) {
      free (_winDoubleParamValueData);
      _winDoubleParamValueData = NULL;
    }

    if (_winStrParamIdxNameData != NULL) {
      free (_winStrParamIdxNameData);
      _winStrParamIdxNameData = NULL;
    }

    if (_winStrParamNameData != NULL) {
      free (_winStrParamNameData);
      _winStrParamNameData = NULL;
    }

    if (_winStrParamIdxValueData != NULL) {
      free (_winStrParamIdxValueData);
      _winStrParamIdxValueData = NULL;
    }

    if (_winStrParamValueData != NULL) {
      free (_winStrParamValueData);
      _winStrParamValueData = NULL;
    }

    if (_intraComm != MPI_COMM_NULL) {
      MPI_Comm_free(&_intraComm);
    }
    if (_intraGroup != MPI_GROUP_NULL) {
      MPI_Group_free(&_intraGroup);
    }

    if ( _intraConnectableComm != MPI_COMM_NULL) {
      MPI_Comm_free(&_intraConnectableComm);
    }

    if ( _intraConnectableGroup != MPI_GROUP_NULL) {
      MPI_Group_free(&_intraConnectableGroup);
    }

    if (_intraRanks != NULL) {
      delete _intraRanks;
    }

    if (_connectableRanks != NULL) {
      delete _connectableRanks;
    }
  }
}

/**
 * \endcond
 */
