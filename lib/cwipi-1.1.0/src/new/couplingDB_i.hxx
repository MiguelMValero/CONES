#ifndef __COUPLING_DB_I_H__
#define __COUPLING_DB_I_H__
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

#include <cassert>

#include <pdm_error.h>

#include "codeProperties.hxx"

using namespace std;

namespace cwipi {

  /**
   * \brief Return a coupling object from it identifier
   *
   * \param [in]  cplId              Coupling identifier
   *
   */

  Coupling& 
  CouplingDB::couplingGet
  (
   const CodeProperties &localCodeProperties,
   const string         &cplId
  ) 
  {
    typedef const map < const cwipi::CodeProperties *, map <string, Coupling * > > ::iterator Iterator;
    typedef map <string, Coupling * > ::iterator Iterator2;
    Iterator p = _couplingDB.find(&localCodeProperties);
    Iterator2 p1;
    if (p == _couplingDB.end()) {
      PDM_error(__FILE__, __LINE__, 0, 
                "'%s' coupling not found for '%s' code\n", cplId.c_str(), 
                localCodeProperties.nameGet().c_str());
    }
    else {
      p1 = p->second.find(cplId);
      if (p1 == p->second.end()) {
        PDM_error(__FILE__, __LINE__, 0, 
                    "'%s' coupling not found '%s' code\n", cplId.c_str(),
                   localCodeProperties.nameGet().c_str());
      }

    }
    assert( p1->second != NULL);
    return *p1->second;
  }
    
  /**
   * \brief Return if a coupling identifier exists  
   *
   * \param [in]  localCodeProperties  Source code
   * \param [in]  cplId                Coupling identifier
   *
   * \return status
   */

  bool 
  CouplingDB::couplingIs
  (
   const CodeProperties &localCodeProperties,
   const string &cplId
  )
  {
    bool status = true;
    typedef const map < const cwipi::CodeProperties *, map <string, Coupling * > > ::iterator Iterator;
    typedef map <string, Coupling * > ::iterator Iterator2;
    Iterator p = _couplingDB.find(&localCodeProperties);
    Iterator2 p1;
    if (p == _couplingDB.end()) {
      status = false;
    }
    else {
      p1 = p->second.find(cplId);
      if (p1 == p->second.end()) {
        status = false;
      }

    }
    return status;
  }
  
  
  
}

#endif
