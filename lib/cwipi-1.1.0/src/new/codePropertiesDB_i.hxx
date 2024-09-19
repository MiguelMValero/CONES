#ifndef __CODE_PROPERTIES_DB_I_H__
#define __CODE_PROPERTIES_DB_I_H__
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
#include <cstdarg>
#include <vector>
#include <typeinfo>
#include <type_traits>

#include "codeProperties.hxx"
#include "codePropertiesDB.hxx"

using namespace std;

namespace cwipi {

  /**
   * \brief Taking into account of a proxy function to redirect output
   *
   * \param [in]  proxyFunction  Function
   *
   */

  void 
  CodePropertiesDB::printfProxySet
  (
   PDM_printf_proxy_t *const proxyFunction
  )
  {
    if (proxyFunction != NULL)
      PDM_printf_proxy_set(proxyFunction);
  }

  /**
   * \brief Return local code MPI intra communicator.
   *
   * \param[in]   localCodeName  Local code name
   * 
   * \return  MPI Intra communicator
   *
   */

  const MPI_Comm &
  CodePropertiesDB::intraCommGet(const string & localCodeName) const
  {
    return _locCodePropertiesDB[localCodeName]->intraCommGet();
  }

  const MPI_Comm &
  CodePropertiesDB::connectableCommGet(const string & localCodeName) const
  {
    return _locCodePropertiesDB[localCodeName]->connectableCommGet();
  }

 /**
   * \brief Set the user structure
   *
   */

   void
  CodePropertiesDB::userStructureSet 
  (
    const string & localCodeName,
    void *userStruct
  )
  {
    _locCodePropertiesDB[localCodeName]->userStructureSet(userStruct);
  }


 /**
   * \brief Get the user structure
   *
   */

  void * 
  CodePropertiesDB::userStructureGet
  (
     const string & localCodeName
  ) const
  {
    return _locCodePropertiesDB[localCodeName]->userStructureGet();
  }


  /**
   * \brief Return MPI communicator containing all processes of all codes.
   *
   * \return  Global MPI Communicator
   *
   */

  const MPI_Comm &
  CodePropertiesDB::globalCommGet() const
  {    
    return _globalComm;
  }

    
  /**
   * \brief Return the code properties.
   *
   * \param [in]  codeName  Code name
   *
   * \return      Properties
   *
   */

  inline  CodeProperties &
  CodePropertiesDB::codePropertiesGet
  (
   const string &codeName
  ) const
  {
    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);
    if (p == _codePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' code not found \n", codeName.c_str());
    assert(p->second != NULL);
    return *(p->second);
  }


    
  /**
   * \brief Return the number of codes known to CWIPI
   *
   * \return   Number of codes
   *
   */

  int
  CodePropertiesDB::codesNbGet
  (
  ) const
  {
    return _codePropertiesDB.size();
  }

    
  /**
   * \brief Return the number of localccodes known to CWIPI
   *
   * \return   Number of local codes
   *
   */

  int
  CodePropertiesDB::localCodesNbGet
  (
  ) const
  {
    return _locCodePropertiesDB.size();    
  }

    
  /**
   * \brief Return the number of codes known to CWIPI
   *
   * \return   Number of codes
   *
   */

  const char **
  CodePropertiesDB::codesListGet
  (
  ) const
  {
    const int nCodes = _codePropertiesDB.size();
    const char **list = (const char **) malloc(sizeof(char *) * nCodes);
    
    typedef map <string, CodeProperties * >::iterator CI;
    int i = 0;
    for (CI p = _codePropertiesDB.begin();
         p != _codePropertiesDB.end(); p++) {
      list[i++] = p->first.c_str();
    }
    return list;
  }

    
  /**
   * \brief Return the number of localccodes known to CWIPI
   *
   * \return   Number of local codes
   *
   */

  const char **
  CodePropertiesDB::localCodesListGet
  (
  ) const
  {
    const int nCodes = _locCodePropertiesDB.size();
    const char **list = (const char **) malloc(sizeof(char *) * nCodes);
    
    typedef map <string, CodeProperties * >::iterator CI;
    int i = 0;
    for (CI p = _locCodePropertiesDB.begin();
         p != _locCodePropertiesDB.end(); p++) {
      list[i++] = p->first.c_str();
    }
    return list;
  }

  /**
   * \brief Add a control paramater.
   *
   * \param [in]  localCodeName   Local code name
   * \param [in]  name            Parameter name
   * \param [in]  value           Initial value 
   *
   */

  template < typename T > 
  void 
  CodePropertiesDB::ctrlParamAdd
  (
   const string &localCodeName, 
   const string &name, 
   const T       value
  )
  {
    
    const map <string, CodeProperties * >::iterator p = 
      _locCodePropertiesDB.find(localCodeName);
    if (p == _locCodePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' is not a local code \n", localCodeName.c_str());
    p->second->ctrlParamAdd(name, value);
  }


  /**
   * \brief set a control paramater.
   *
   * \param [in]  localCodeName   Local code name
   * \param [in]  name            Parameter name
   * \param [in]  value           Initial value 
   *
   */

  template < typename T > 
  void 
  CodePropertiesDB::ctrlParamSet
  (
   const string &localCodeName, 
   const string &name, 
   const T       value
  )
  {
    const map <string, CodeProperties * >::iterator p = 
      _locCodePropertiesDB.find(localCodeName);
    if (p == _locCodePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' is not a local code \n", localCodeName.c_str());
    p->second->ctrlParamSet(name, value);
  }

    
  /**
   * \brief Cancel a control paramater.
   *
   * \param [in]  localCodeName   Local code name
   * \param [in]  name            Parameter name
   *
   */

  template < typename T > 
  void 
  CodePropertiesDB::ctrlParamCancel
  (
   const string &localCodeName, 
   const string &name
  )
  {
    const map <string, CodeProperties * >::iterator p = 
      _locCodePropertiesDB.find(localCodeName);
    if (p == _locCodePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' is not a local code \n", localCodeName.c_str());
    p->second->ctrlParamCancel<T>(name);
  }


  /**
   * \brief Return the number of parameters
   *
   * \param [in]  codeName  Code name
   *
   * \return  Number of parameters
   *
   */

  template < typename T > 
  int 
  CodePropertiesDB::ctrlParamNGet
  (
   const string &codeName
  ) const
  {
    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);
    if (p == _codePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' code not found \n", codeName.c_str());
    return p->second->ctrlParamNGet<T>();    
  }


  /**
   * \brief Get the value of a control paramater.
   *
   * \param [in]  codeName  Code name
   * \param [in]  name      Parameter name
   *
   * \return             Value           
   *
   */

  template < typename T > 
  const T
  CodePropertiesDB::ctrlParamGet
  (
    const string &codeName,
    const string &name
  )
  {

    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);
    if (p == _codePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' code not found \n", codeName.c_str());
    T value;
    p->second->ctrlParamGet(name, &value);
    return value;
  }

  
  /**
     * \brief Reduce a parameter through a list of codes. The available processes
     *        are sum, max and min.
     *
     * \param [in]  op          Operator from \ref CWP_Op_t
     * \param [in]  name        Parameter name
     * \param [in]  nCode       Number of code
     * \param       code_names  Code names
     *
     * \return             Operation result
     *
     */

  template < typename T > 
  void
  CodePropertiesDB::ctrlParamReduce
  (
    const CWP_Op_t  op,
    const string    &name,
    T               *res,
    const int        nCode,
    const char     **code_names
  )
  {
    string codeName1 = string((char*) code_names[0]);

    *res = this->ctrlParamGet < T > (codeName1, name);

    for (int k = 1; k < nCode; k++) {

      string codeName = string((char*) code_names[k]);

      const map <string, CodeProperties * >::iterator p = 
        _codePropertiesDB.find(codeName);

      if (p == _codePropertiesDB.end())
        PDM_error(__FILE__, __LINE__, 0,
                   "'%s' code not found \n", codeName.c_str());

      T distParam = this->ctrlParamGet < T > (codeName, name);

      switch (op) {
      case CWP_OP_MAX:
        *res = max(distParam, *res);
        break;
      
      case CWP_OP_MIN:
        *res = min(distParam, *res);
        break;

      case CWP_OP_SUM:
        *res += distParam;
        break;
      }
    }
  }

  /**
   * \brief Lock access to local parameters from a distant code  
   *
   * \param [in]  codeName  Local code name to lock
   *
   */
  
  void 
  CodePropertiesDB::lock
  (
   const string &codeName
  )
  {
    const map <string, CodeProperties * >::iterator p = 
      _locCodePropertiesDB.find(codeName);
    if (p == _locCodePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' is not a local code \n", codeName.c_str());
    p->second->paramLock();
  }

  /**
   * \brief Is locked param  
   *
   * \param [in]  codeName  Local code name to lock
   *
   */
  
  int 
  CodePropertiesDB::isLocked
  (
   const string &codeName
  )
  {
    const map <string, CodeProperties * >::iterator p = 
      _locCodePropertiesDB.find(codeName);
    if (p == _locCodePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' is not a local code \n", codeName.c_str());
    return p->second->paramIsLocked();
  }
  
  /**
   * \brief unlock access to local parameters from a distant code  
   *
   * \param [in]  codeName  Local code name to unlock
   *
   */

  void 
  CodePropertiesDB::unLock
  (
   const string &codeName
  )
  {
    const map <string, CodeProperties * >::iterator p = 
      _locCodePropertiesDB.find(codeName);
    if (p == _locCodePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                "'%s' is not a local code \n", codeName.c_str());
    p->second->paramUnLock();
  }
  
  /**
   * \brief Return of local parameters
   *
   * \param [in]  codeName  Code name
   *
   * \return  List of parameters
   *
   */

  template < typename T > 
  void
  CodePropertiesDB::ctrlParamListGet
  (
   const string &codeName,
   int  *nParam,
   char ***paramNames
  ) const
  {
    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);

    if (p == _codePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                 "'%s' code not found \n", codeName.c_str());

    return p->second->ctrlParamListGet<T>(nParam, paramNames);
  }

  /**
   * \brief Chek name parameter
   *
   * \param [in]  codeName  Code name
   * \param [in]  name  Parameter name to check
   *
   * \return  1 : true / 0 : false
   *
   */

  template < typename T > 
  int
  CodePropertiesDB::ctrlParamIs
  (
   const string &codeName,
   const string &name
   ) const
  {
    const map <string, CodeProperties * >::iterator p = 
      _codePropertiesDB.find(codeName);

    if (p == _codePropertiesDB.end())
      PDM_error(__FILE__, __LINE__, 0,
                 "'%s' code not found \n", codeName.c_str());

    return p->second->ctrlParamIs<T>(name);
  }
  

} // namespace cwipi

#endif /* __CODE_PROPERTIES_DB_H__ */
