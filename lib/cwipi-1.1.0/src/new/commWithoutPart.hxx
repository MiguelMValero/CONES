#ifndef __COMMWITHOUTPART_H__
#define __COMMWITHOUTPART_H__
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
#include "cwp.h"

namespace cwipi {

  /** 
   * \class CommWithoutPart commWithoutPart.hxx "commWithoutPart.hxx"
   * \brief Communication without support partitionning
   *
   *  This class manages communication when the support is unpartitionned
   *
   */

  class CommWithoutPart :
    public Communication {
    
  public :

    /**
     *
     * \brief Constructor.
     *
     */

    CommWithoutPart();

    /**
     * \brief Destructor.
     *
     */

    virtual ~CommWithoutPart();

    /**
     *
     * \brief Return Communication Type
     *
     * \return   communication type
     *
     */

    inline CWP_Comm_t commTypeGet(); 

    /**
     *
     * \brief Synchronise
     *
     */
/*
    virtual void
    sync
    (
     void *tab, 
     MPI_Datatype mpiType, 
     int tabSize
    );
    
    */

  private :

    /**
     *
     * \brief Assigment operator
     *
     */
    
    CommWithoutPart 
    &operator=
    (const CommWithoutPart &other);

    /**
     *
     * \brief Copy constructor
     *
     */

    CommWithoutPart(const CommWithoutPart& other); 

  private :

    /**
     *
     * \brief Building coupling communicator.
     *
     */

    void 
    _cplCommCreate
    (
     CWP_Comm_t cplCodeCommType
     );

    /**
     *
     * \brief Communication type
     *
     */

    CWP_Comm_t _commType;

  };

  /**
   *
   * \brief Return Communication Type
   *
   * \return   communication type
   *
   */
  
  CWP_Comm_t 
  CommWithoutPart::commTypeGet()
  {
    return _commType;
  }

}
#endif
