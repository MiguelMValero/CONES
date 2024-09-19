#ifndef __CODE_PROPERTIES_DB_H__
#define __CODE_PROPERTIES_DB_H__
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
#include <map>
#include <vector>

#include "pdm_printf.h"

#include "cwp.h"
#include "singleton.hpp"

using namespace std;

namespace cwipi {

  class CodeProperties;

  /** 
   * \class CodePropertiesDB 
   *        codePropertiesDB.hxx 
   *        "codePropertiesDB.hxx"
   *
   * \brief Codes properties management.
   *
   *  This class manages codes properties :
   *  - Local control parameters,
   *  - Distant control parameters,
   *  - MPI communicators
   *  - . 
   * 
   */

  class CodePropertiesDB 
    : public Singleton <CodePropertiesDB>
  {

    friend class Singleton <CodePropertiesDB>;

  public:

    /**
     * \brief MPI Communicator Initialization.
     *
     * This function builds the current code intra-communicator from
     * the current name and the MPI communicator containing all processes of
     * all codes.
     *
     * \param [in]  globalComm      MPI communicator containing all processes 
     *                              of all codes
     * \param [in]  n_codes         Number of codes on the current rank
     * \param [in]  code_names      Codes names on the current rank 
     * \param [in]  is_active_rank  Current rank is active
     * \param [in]  n_param_max     Maximum number of parameters
     * \param [in]  str_size_max    Maximum size for a string
     * \param [out] intra_comms     Current codes intra-communicators 
     *
     */

    void 
    init
    (
     const MPI_Comm     globalComm,
     const int          n_codes,
     const char**       code_names, 
     const CWP_Status_t is_active_rank,
     const int          n_param_max,
     const int          str_size_max,      
     MPI_Comm           *intra_comms
    );

    /**
     * \brief Taking into account of a proxy function to redirect output
     *
     * \param [in]  proxyFunction  Function
     *
     */

    inline void 
    printfProxySet
    (
     PDM_printf_proxy_t *const proxyFunction
    );

    /**
      * \brief Return local code MPI intra communicator.
      *
      * \param[in]   localCodeName  Local code name
      * 
      * \return  MPI Intra communicator
      *
      */

    inline const MPI_Comm &
    intraCommGet
    (
     const string & localCodeName
    ) const;

    inline const MPI_Comm &
    connectableCommGet
    (
     const string & localCodeName
    ) const;


   /**
     * \brief Set the user structure
     *
     */

    inline void
    userStructureSet 
    (
      const string & localCodeName,
      void *userStruct
    );


   /**
     * \brief Get the user structure
     *
     */

    inline void * 
    userStructureGet
    (
       const string & localCodeName
    ) const;

     
    /**
     * \brief Return MPI communicator containing all processes of all codes.
     *
     * \return  Global MPI Communicator
     *
     */

    inline const MPI_Comm &
    globalCommGet() const;

    
    /**
     * \brief Return the code properties.
     *
     * \param [in]  codeName  Code name
     *
     * \return      Properties
     *
     */

    inline CodeProperties &
    codePropertiesGet
    (
     const string &codeName
    ) const;

    
    /**
     * \brief Return the number of codes known to CWIPI
     *
     * \return   Number of codes
     *
     */

    inline int
    codesNbGet
    (
    ) const;

    
    /**
     * \brief Return the number of localccodes known to CWIPI
     *
     * \return   Number of local codes
     *
     */

    inline int
    localCodesNbGet
    (
    ) const;

    
    /**
     * \brief Return the number of codes known to CWIPI
     *
     * \return   Number of codes
     *
     */

    inline const char **
    codesListGet
    (
    ) const;

    
    /**
     * \brief Return the number of localccodes known to CWIPI
     *
     * \return   Number of local codes
     *
     */

    inline const char **
    localCodesListGet
    (
    ) const;


    /**
     * \brief Set a control paramater.
     *
     * \param [in]  localCodeName   Local code name
     * \param [in]  name            Parameter name
     * \param [in]  value           Initial value 
     *
     */
    
    template < typename T > 
    void 
    ctrlParamAdd
    (
     const string &localCodeName, 
     const string &name, 
     const T       value
    );


    /**
     * \brief Set a control paramater.
     *
     * \param [in]  localCodeName   Local code name
     * \param [in]  name            Parameter name
     * \param [in]  value           Initial value 
     *
     */
    
    template < typename T > 
    void 
    ctrlParamSet
    (
     const string &localCodeName, 
     const string &name, 
     const T       value
    );

    
    /**
     * \brief Cancel a control paramater.
     *
     * \param [in]  localCodeName   Local code name
     * \param [in]  name            Parameter name
     *
     */

    template < typename T > 
    void 
    ctrlParamCancel
    (
     const string &localCodeName, 
     const string &name
    );

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
    ctrlParamNGet
    (
     const string &codeName
    ) const;


    /**
     * \brief Return of the control parameter list
     *
     * \param [in]  codeName    Code name
     * \param [in]  nParam      Number of control parameters
     * \param [in]  paramNames   Pointer to the list of parameters names
     *
     */

    template < typename T > 
    void
    ctrlParamListGet
    (
     const string &codeName,
     int  *nParam,
     char ***paramNames
    ) const;

    /**
     * \brief Chek name parameter existence
     *
     * \param [in]  codeName   Code name
     * \param [in]  name       Parameter to check
     *
     * \return  1 : true / 0 : false
     *
     */

    template < typename T > 
    int
    ctrlParamIs
    (
     const string &codeName,
     const string &name
    ) const;


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
    ctrlParamGet
    (
      const string &codeName,
      const string &name
    );

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
    ctrlParamReduce
    (
     const CWP_Op_t  op, 
     const string    &name,
     T               *res,
     const int        nCode,
     const char     **code_names
    );

    /**
     * \brief Dump properties.  
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
     * \brief Lock access to local parameters from a distant code  
     *
     * \param [in]  codeName  Code name to lock
     *
     */

    inline void 
    lock
    (
    const string &codeName
    );


    /**
     * \brief Is locked param  
     *
     * \param [in]  codeName  Local code name to lock
     *
     */
    
    inline int 
    isLocked
    (
     const string &codeName
    );


    /**
     * \brief unlock access to local parameters from a distant code  
     *
     * \param [in]  codeName  Code name to unlock
     *
     */

    inline void 
    unLock
    (
    const string &codeName
    );

  private:

    /**
     * \brief Default constructor.
     *
     */

    CodePropertiesDB();

    /**
     * \brief Copy constructor (unavailable)
     *
     */

    CodePropertiesDB
    (
     const CodePropertiesDB &other
     );

    /**
     * \brief Default assignment (unavailable)
     *
     */

    CodePropertiesDB & 
    operator=
    (
     const CodePropertiesDB &other
     );

    /**
     * \brief Destructor.
     *
     */

    virtual ~CodePropertiesDB();

  private:
    MPI_Comm                          _globalComm;             /*!< Global communicator */  
    map <string, CodeProperties * > & _codePropertiesDB;       /*!< Distant code 
                                                                    properties data base */
    map <string, CodeProperties * > & _locCodePropertiesDB;    /*!< Local code properties */
    
    int                                _n_param_max;           /*!< Maximum number of parameters */  
    int                                _str_size_max;          /*!< Maximum size for a string */

  private:
    static const int _nIssend;                                        /*!< Number of issend 
                                                                           to send parameters */
  };
}

#endif /* __CODE_PROPERTIES_DB_H__ */
