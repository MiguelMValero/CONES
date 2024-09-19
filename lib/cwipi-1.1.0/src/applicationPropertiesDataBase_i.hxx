#ifndef __APPLICATION_PROPERTIES_DATA_BASE_I_H__
#define __APPLICATION_PROPERTIES_DATA_BASE_I_H__
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

#include <cassert>

#include "applicationProperties.hxx"


namespace cwipi {

  inline void ApplicationPropertiesDataBase::setPrintfProxy(bftc_printf_proxy_t *const callBackPrintf)
  {
    if (callBackPrintf != NULL)
      bftc_printf_proxy_set(callBackPrintf);
  }

  // Access to local MPI properties (synchronization with a second application)

  inline const MPI_Comm &ApplicationPropertiesDataBase::getLocalComm() const
  {
    return _localApplicationProperties->getLocalComm();
  }

  inline const std::string &ApplicationPropertiesDataBase::getLocalName() const
  {
    return _localApplicationProperties->getName();
  }

  inline const MPI_Comm &ApplicationPropertiesDataBase::getGlobalComm() const
  {
    return _localApplicationProperties->getGlobalComm();
  }

  inline const int &ApplicationPropertiesDataBase::getBeginningRank() const
  {
    return _localApplicationProperties->getBeginningRank();
  }

  inline const int &ApplicationPropertiesDataBase::getEndRank() const
  {
    return _localApplicationProperties->getEndRank();
  }

  // Access to distant MPI properties (synchronization with a second application)

//   inline const MPI_Comm &ApplicationPropertiesDataBase::getDistantLocalComm(const std::string &applicationName)
//   {
//     std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
//     if (p != _distantApplicationPropertiesDataBase.end())
//       return p->second->getLocalComm();
//     else
//       bftc_error(__FILE__, __LINE__, 0,
//                 "'%s' application not found \n", applicationName.c_str());
//   }

  inline const int &ApplicationPropertiesDataBase::getDistantBeginningRank(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getBeginningRank();
  }

  inline const int &ApplicationPropertiesDataBase::getDistantEndRank(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getEndRank();
  }

  // Access to local control parameters

  inline void ApplicationPropertiesDataBase::addLocalIntControlParameter(const std::string &name, const int initialValue)
  {
    _localApplicationProperties->addIntControlParameter(name, initialValue);
  }

  inline void ApplicationPropertiesDataBase::addLocalDoubleControlParameter(const std::string &name, const double initialValue)
  {
    _localApplicationProperties->addDoubleControlParameter(name, initialValue);
  }

  inline void ApplicationPropertiesDataBase::addLocalStringControlParameter(const std::string &name, const std::string initialValue)
  {
    _localApplicationProperties->addStringControlParameter(name, initialValue);
  }

  inline void ApplicationPropertiesDataBase::setLocalIntControlParameter(const std::string &name, const int value)
  {
    _localApplicationProperties->setIntControlParameter(name, value);
  }

  inline void ApplicationPropertiesDataBase::setLocalDoubleControlParameter(const std::string &name, const double value)
  {
    _localApplicationProperties->setDoubleControlParameter(name, value);
  }

  inline void ApplicationPropertiesDataBase::setLocalStringControlParameter(const std::string &name, const std::string value)
  {
    _localApplicationProperties->setStringControlParameter(name, value);
  }

  inline const int &ApplicationPropertiesDataBase::getLocalIntControlParameter(const std::string &name)
  {
    return _localApplicationProperties->getIntControlParameter(name);
  }

  inline const double &ApplicationPropertiesDataBase::getLocalDoubleControlParameter(const std::string &name)
  {
    return _localApplicationProperties->getDoubleControlParameter(name);
  }

  inline const std::string &ApplicationPropertiesDataBase::getLocalStringControlParameter(const std::string &name)
  {
    return _localApplicationProperties->getStringControlParameter(name);
  }

  inline void ApplicationPropertiesDataBase::eraseLocalIntControlParameter(const std::string &name)
  {
    _localApplicationProperties->eraseIntControlParameter(name);
  }

  inline void ApplicationPropertiesDataBase::eraseLocalDoubleControlParameter(const std::string &name)
  {
    _localApplicationProperties->eraseDoubleControlParameter(name);
  }

  inline void ApplicationPropertiesDataBase::eraseLocalStringControlParameter(const std::string &name)
  {
    _localApplicationProperties->eraseStringControlParameter(name);
  }

  // Access to distant control parameters

  inline const int &ApplicationPropertiesDataBase::getDistantIntControlParameter(const std::string &applicationName, const std::string &name) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getIntControlParameter(name);
  }

  inline const ApplicationProperties &ApplicationPropertiesDataBase::getDistantApplicationProperties(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    assert(p->second != NULL);
    return *(p->second);
  }

  inline const ApplicationProperties &ApplicationPropertiesDataBase::getLocalApplicationProperties() const
  {
    assert(_localApplicationProperties != NULL);
    return *_localApplicationProperties;
  }

  inline const double &ApplicationPropertiesDataBase::getDistantDoubleControlParameter(const std::string &applicationName, const std::string &name) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getDoubleControlParameter(name);
  }

  inline const std::string &ApplicationPropertiesDataBase::getDistantStringControlParameter(const std::string &applicationName, const std::string &name) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getStringControlParameter(name);
  }

  inline int ApplicationPropertiesDataBase::hasLocalIntControlParameter(const std::string &name) const
  {
    return  _localApplicationProperties->hasIntControlParameter(name);
  }

  inline int ApplicationPropertiesDataBase::hasLocalDoubleControlParameter(const std::string &name) const
  {
    return  _localApplicationProperties->hasDoubleControlParameter(name);
  }

  inline int ApplicationPropertiesDataBase::hasLocalStringControlParameter(const std::string &name) const
  {
    return  _localApplicationProperties->hasStringControlParameter(name);
  }

  inline int ApplicationPropertiesDataBase::getLocalNIntControlParameter() const
  {
    return  _localApplicationProperties->getNIntControlParameter();
  }

  inline int ApplicationPropertiesDataBase::getLocalNDoubleControlParameter() const
  {
    return  _localApplicationProperties->getNDoubleControlParameter();
  }

  inline int ApplicationPropertiesDataBase::getLocalNStringControlParameter() const
  {
    return  _localApplicationProperties->getNStringControlParameter();
  }

  inline char** ApplicationPropertiesDataBase::getLocalListIntControlParameter() const
  {
    return  _localApplicationProperties->getListIntControlParameter();
  }

  inline char** ApplicationPropertiesDataBase::getLocalListDoubleControlParameter() const
  {
    return  _localApplicationProperties->getListDoubleControlParameter();
  }

  inline char** ApplicationPropertiesDataBase::getLocalListStringControlParameter() const
  {
    return  _localApplicationProperties->getListStringControlParameter();
  }

  inline int ApplicationPropertiesDataBase::hasDistantIntControlParameter(const std::string &applicationName, const std::string &name) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->hasIntControlParameter(name);
  }

  inline int ApplicationPropertiesDataBase::hasDistantDoubleControlParameter(const std::string &applicationName, const std::string &name) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->hasDoubleControlParameter(name);
  }

  inline int ApplicationPropertiesDataBase::hasDistantStringControlParameter(const std::string &applicationName, const std::string &name) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->hasStringControlParameter(name);
  }

  inline int ApplicationPropertiesDataBase::getDistantNIntControlParameter(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getNIntControlParameter();
  }

  inline int ApplicationPropertiesDataBase::getDistantNDoubleControlParameter(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getNDoubleControlParameter();
  }

  inline int ApplicationPropertiesDataBase::getDistantNStringControlParameter(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getNStringControlParameter();
  }

  inline char** ApplicationPropertiesDataBase::getDistantListIntControlParameter(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getListIntControlParameter();
  }

  inline char** ApplicationPropertiesDataBase::getDistantListDoubleControlParameter(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getListDoubleControlParameter();
  }

  inline char** ApplicationPropertiesDataBase::getDistantListStringControlParameter(const std::string &applicationName) const
  {
    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());
    return p->second->getListStringControlParameter();
  }

} // namespace cwipi

#endif /* __APPLICATION_PROPERTIES_DATA_BASE_H__ */
