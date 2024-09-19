/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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
///////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009  ONERA
//
// Description :
//    This class describes application properties :
//        - intra mpi communicator (localComm)
//        - beginning and end rank in the inter mpi communicator (globalComm)
//        - control parameters storage
//
///////////////////////////////////////////////////////////////////////////////

#include <mpi.h>

#include "applicationProperties.hxx"

namespace cwipi
{

  ApplicationProperties::ApplicationProperties(std::string &name,
                                               const MPI_Comm globalComm)

    : _name(name), _globalComm(globalComm),
    _intControlParameters(*(new std::map <std::string, int>())),
    _doubleControlParameters(*(new std::map <std::string, double>())),
    _stringControlParameters(*(new std::map <std::string, std::string>()))
  {
    _localComm = MPI_COMM_NULL;
    _beginningRank = -999;
    _endRank = -999;
  }

  ApplicationProperties::ApplicationProperties(const ApplicationProperties& other)
    : _name(other._name), _globalComm(other._globalComm), _localComm(other._localComm),
      _beginningRank(other._beginningRank), _endRank(other._endRank),
      _intControlParameters(other._intControlParameters),
      _doubleControlParameters(other._doubleControlParameters),
      _stringControlParameters(other._stringControlParameters)
  {
  }

  void ApplicationProperties::dump()
  {
    bftc_printf("'%s' properties\n",_name.c_str());
    bftc_printf("  - Ranks in global MPI_comm : %i <= ranks <= %i \n",
               _beginningRank,
               _endRank);
    bftc_printf("  - Int Control Parameter :\n");

    typedef std::map <std::string, int>::iterator Iterator1;
    for (Iterator1 p = _intControlParameters.begin(); p != _intControlParameters.end(); p++)
      bftc_printf("   * '%s' : %i\n", p->first.c_str(), p->second);
    bftc_printf("  - Double Control Parameter :\n");

    typedef std::map <std::string, double>::iterator Iterator2;
    for (Iterator2 p = _doubleControlParameters.begin(); p != _doubleControlParameters.end(); p++)
      bftc_printf("   * '%s' : %12.5e\n", p->first.c_str(), p->second);

    bftc_printf("  - String Control Parameter :\n");

    typedef std::map <std::string, std::string>::iterator Iterator3;
    for (Iterator3 p = _stringControlParameters.begin(); p != _stringControlParameters.end(); p++)
      bftc_printf("   * '%s' : '%s'\n", p->first.c_str(), p->second.c_str());

    bftc_printf_flush();
  }

  ApplicationProperties::~ApplicationProperties()
  {

    if (!_intControlParameters.empty())
      _intControlParameters.clear();
    delete &_intControlParameters;
    if (!_doubleControlParameters.empty())
      _doubleControlParameters.clear();
    delete &_doubleControlParameters;
    if (!_stringControlParameters.empty())
      _stringControlParameters.clear();
    delete &_stringControlParameters;
    if (_localComm != MPI_COMM_NULL)
      MPI_Comm_free(&_localComm);

  }
}
