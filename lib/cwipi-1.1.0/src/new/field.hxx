#ifndef __FIELD_H__
#define __FIELD_H__
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

#include <sstream>
#include <mesh.hxx>
#include <map>

#include <coupling.hxx>
#include "cwp.h"
#include "cwp_priv.h"
#include "pdm_writer.h"


/**
 * \cond
 */

namespace cwipi {

  /**
   * \class Field field.hxx "field.hxx"
   * \brief Abstract field
   *
   *  This class is field abstract interface
   *
   */
  class Mesh;
  class Field {

  public:

    /**
     * \brief Constructor
     *
     */

    Field() {}

    Field (std::string            field_id    ,
           int                    fieldIDInt,      
           CWP_Type_t             dataType    ,
           Coupling*              cpl         ,
           CWP_Dof_location_t      fieldType   ,
           CWP_Field_storage_t    storage     ,
           int                    nComponent  ,
           CWP_Field_exch_t       exchangeType,
           CWP_Status_t           visuStatus);

    /**
     * \brief Destructor
     *
     */

    ~Field();

    /**
     * \brief set data array
     *
     * \param [in] data   Data array
     *
     */

    void dataSet (int i_part, const CWP_Field_map_t  map_type, void* data);



    /**
     * \brief Write th
     *
     * \param [in] data   Data array
     *
     */

    void write 
    (
      CWP_Field_exch_t exch_type
    );


    /**
     *
     * \brief Get field storage type
     *
     * \return            Field storage type
     *
     */

    inline CWP_Field_storage_t
    storageTypeGet() const
    {
      return _storage;
    }

    /**
     *
     * \brief Get nunmber of field components
     *
     * \return             Number of field components
     *
     */

    inline int
    nComponentGet() const
    {
      return _nComponent;
    }

    /**
     *
     * \brief Get nunmber of field degrees of freedom
     *
     * \param [in]   i_part         Partition identifier
     *
     * \return             Number of field degrees of freedom
     *
     */

    int
    nDOFGet
    (
     int i_part
    );

    inline std::string
    fieldIDGet() const
    {
      return _fieldID;
    }


    inline int
    fieldIDIntGet() const
    {
      return _fieldIDInt;
    }

    inline int *
    _id_writer_var_send_get()
    {
      return _id_writer_var_send;
    }

    inline int
    _id_writer_var_send_status_get()
    {
      return _id_writer_var_send_status;
    }

    inline int *
    _id_writer_var_recv_get()
    {
      return _id_writer_var_recv;
    }

    inline int
    _id_writer_var_recv_status_get()
    {
      return _id_writer_var_recv_status;
    }

    /**
     *
     * \brief Get field nature
     *
     * \return           Field nature
     *
     */

    inline CWP_Dof_location_t
    locationGet() const
    {
      return _fieldLocation;
    }

    CWP_Type_t
    dataTypeGet() const
    {
      return _dataType;
    }

    int
    dataTypeSizeGet() const
    {
      return _dataTypeSize;
    }

    /**
     *
     * \brief Get exchange type
     *
     * \return          Exchange type
     *
     */

    inline CWP_Field_exch_t
    exchangeTypeGet() const
    {
      return _exchangeType;
    }

    /**
     *
     * \brief Get visu status
     *
     * \return          Exchange type
     *
     */

    inline CWP_Status_t
    visuStatusGet()  const
    {
      return _visuStatus;
    }

    /**
     *
     * \brief Get data
     *
     * \return          Data
     *
     */

    void* dataGet(int i_part,  const CWP_Field_map_t  map_type) const
    {
      if (map_type == CWP_FIELD_MAP_SOURCE) {
        return _data_src[i_part];
      }
      else if (map_type == CWP_FIELD_MAP_TARGET) {
        return _data_tgt[i_part];
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "Field::dataSet Error : unknown data type.\n");
        return nullptr;
      }
    }


    void interpFunctionSet(CWP_Interp_function_t fct)
    {
      _interpolationType     = CWP_INTERPOLATION_USER ;
      _interpolationFunction = fct;
    }

    void interpFunctionFSet(CWP_Interp_function_t fct)
    {
      _interpolationType       = CWP_INTERPOLATION_USER ;
      _interpolationFunction_f = fct;
    }

    void interpFunctionPSet(CWP_Interp_function_p_t fct)
    {
      _interpolationType       = CWP_INTERPOLATION_USER ;
      _interpolationFunction_p = fct;
    }

    CWP_Interp_function_t interpFunctionFGet()
    {
      return _interpolationFunction_f;
    }

    CWP_Interp_function_p_t interpFunctionPGet()
    {
      return _interpolationFunction_p;
    }

    void interpFunctionUnSet()
    {
      _interpolationType     = CWP_INTERPOLATION_DEFAULT ;
      _interpolationFunction = NULL;
    }

    void interpFunctionFUnSet()
    {
      _interpolationType       = CWP_INTERPOLATION_DEFAULT ;
      _interpolationFunction_f = NULL;
    }

    void interpFunctionPUnset()
    {
      _interpolationType       = CWP_INTERPOLATION_DEFAULT ;
      _interpolationFunction_p = NULL;
    }

    CWP_Interp_function_t interpolationFunctionGet() {
      return _interpolationFunction;
    }

    CWP_Interpolation_t interpolationTypeGet() {
      return _interpolationType;
    }


    // std::vector<void*> dataGetAll() const
    // {
    //   return _data;
    // }

    // Field has been sent ?
    void is_send_yet_set (int value)
    {
      _is_send_yet = value;
    }

    int is_send_yet_get ()
    {
      return _is_send_yet;
    }

    void is_send_end_step_set (int value)
    {
      _is_send_end_step = value;
    }

    int is_send_end_step_get ()
    {
      return _is_send_end_step;
    }

    // Field has been received ?
    void is_recv_yet_set (int value)
    {
      _is_recv_yet = value;
    }

    int is_recv_yet_get ()
    {
      return _is_recv_yet;
    }

    void is_recv_end_step_set (int value)
    {
      _is_recv_end_step = value;
    }

    int is_recv_end_step_get ()
    {
      return _is_recv_end_step;
    }

    int computedTgtBcastIsEnabled() const
    {
      return _computed_tgt_bcast_enabled;
    }

    void computedTgtBcastEnable()
    {
      _computed_tgt_bcast_enabled = 1;
    }

    int involvedSrcBcastIsEnabled() const
    {
      return _computed_tgt_bcast_enabled;
    }

    void involvedSrcBcastEnable()
    {
      _computed_tgt_bcast_enabled = 1;
    }


  void lastRequestAdd (int i_proc, MPI_Request request) {
    _last_request[i_proc] = request;
  }


  MPI_Request lastRequestGet (int i_proc) {
    return _last_request[i_proc];

  }

  void lastRequestAdd_p2p (int i_proc, std::vector<MPI_Request> request) {
    _last_request_p2p[i_proc] = request;
  }


  std::vector<MPI_Request> lastRequestGet_p2p (int i_proc) {
    return _last_request_p2p[i_proc];
  }


  void* recvBufferGet () {
    return _recvBuffer;
  }

  void* sendBufferGet () {
    return _sendBuffer;
  }

  void sendBufferSet (void* sendBuffer) {
    _sendBuffer = sendBuffer;
  }

  void linkedFieldLocationSet(CWP_Dof_location_t linkedFieldLocation){
    _linkedFieldLocation = linkedFieldLocation;
  }

  CWP_Dof_location_t linkedFieldLocationGet(){
    return _linkedFieldLocation;
  }

  Coupling *couplingGet(){
    return _cpl;
  }

  Mesh *meshGet(){
    return _mesh;
  }

  void pythonObjectSet(void *p) {
    _python_object = p;
  }

  void *pythonObjectGet() {
    return _python_object;
  }

  private:

    CWP_Field_storage_t                      _storage;        /*!< Storage type */
    int                                      _nComponent;     /*!< Number of component */
    CWP_Dof_location_t                       _fieldLocation;  /*!< Value location Interpolation methods for sender and cloud points type for receiver */
    CWP_Dof_location_t                       _linkedFieldLocation; /*!< Value location Interpolation methods for sender and cloud points type for receiver */
    CWP_Field_exch_t                         _exchangeType;   /*!< Exchange type */
    CWP_Status_t                             _visuStatus;     /*!< Visualization status */
    // std::vector<void* >                      _data_src;       /*!< Pointer to data array Add a another data poiter for send fields */
    // std::vector<void* >                      _data_tgt;       /*!< Pointer to data array Add a another data poiter for recv fields */
    void                                   **_data_src;
    void                                   **_data_tgt;
    CWP_Type_t                               _dataType;
    std::string                              _fieldID;
    int                                      _fieldIDInt;
    void                                    *_sendBuffer;
    void                                    *_recvBuffer;
    Coupling                                *_cpl;
    Mesh                                    *_mesh;
    int                                      _n_part;
    PDM_writer_t                            *_writer;
    int                                     *_id_writer_var_send;
    int                                     *_id_writer_var_recv;
    /* status = -1: not exchanged,
                 0: computed and exchanged,
                 1: not computed and exchanged */
    int                                      _id_writer_var_send_status;
    int                                      _id_writer_var_recv_status;

    std::map <int,MPI_Request>               _last_request;
    std::map <int,std::vector<MPI_Request>>  _last_request_p2p;
    int                                      _dataTypeSize;
    CWP_Interp_function_t                    _interpolationFunction;
    CWP_Interpolation_t                      _interpolationType;
    CWP_Interp_function_t                    _interpolationFunction_f;
    CWP_Interp_function_p_t                  _interpolationFunction_p;

    Field &operator=(const Field &other);       /*!< Assignment operator not available */
    Field (const Field& other);                 /*!< Copy constructor not available */

    int                                      _computed_tgt_bcast_enabled;
    int                                      _involved_src_bcast_enabled;

    // writer
    int _is_send_yet;      /*!< Tells if a field has been sent at a given moment */
    int _is_send_end_step; /*!< Tells if a field has been sent at the end of a step */
    int _is_recv_yet;      /*!< Tells if a field has been received at a given moment */
    int _is_recv_end_step; /*!< Tells if a field has been received at the end of a step */

    void *_python_object;
  };

}

/**
 * \endcond
 */

#endif //__FIELD_H__
