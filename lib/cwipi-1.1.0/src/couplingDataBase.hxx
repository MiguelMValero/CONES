#ifndef __APPLICATION_COUPLING_DATA_BASE_H__
#define __APPLICATION_COUPLING_DATA_BASE_H__
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

#include <map>
#include <string>

#include <mpi.h>

#include <bftc_printf.h>

#include "singleton.hpp"
#include "cwipi.h"

namespace cwipi {
  class oldCoupling;
  class ApplicationProperties;
  
  class CouplingDataBase : public Singleton <CouplingDataBase>
  {
    friend class Singleton <CouplingDataBase>;
    friend class oldCoupling;
    
  public:
    void createCoupling(const std::string &name, 
                        const cwipi_coupling_type_t couplingType,
                        const ApplicationProperties& localApplicationProperties,
                        const ApplicationProperties& coupledApplicationProperties,
                        const int entitiesDim,
                        const double tolerance,
                        const cwipi_solver_type_t solverType,
                        const int    outputFrequency,
                        const char  *outputFormat,
                        const char  *outputFormatOption,
			const int nbLocations);
    
    void deleteCoupling(const std::string &name);
    
    inline oldCoupling& getCoupling(const std::string &name);
    
    
  private:
    CouplingDataBase();
    CouplingDataBase(const CouplingDataBase &other);
    CouplingDataBase & operator=(const CouplingDataBase &other);
    virtual ~CouplingDataBase();
    
  private:
    std::map <std::string, oldCoupling * > & _couplingDataBase;
    MPI_Comm  _fvmComm;

  };
}

#endif
