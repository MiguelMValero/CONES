

#include "mesh.hxx"
#include "field.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "pdm_writer.h"

/**
 * \cond
 */

namespace cwipi {

Field::Field (std::string            field_id    ,
              int                     fieldIDInt,      
              CWP_Type_t             dataType    ,
              Coupling*              cpl        ,
              CWP_Dof_location_t      fieldType   ,
              CWP_Field_storage_t    storage     ,
              int                    nComponent  ,
              CWP_Field_exch_t       exchangeType,
              CWP_Status_t           visuStatus):
   _storage        (storage)     ,
   _nComponent     (nComponent)  ,
   _fieldLocation  (fieldType)   ,
   _linkedFieldLocation (CWP_DOF_LOCATION_UNDEF),
   _exchangeType   (exchangeType),
   _visuStatus     (visuStatus)  ,
   _dataType       (dataType)    ,
   _fieldID        (field_id)    ,
   _fieldIDInt     (fieldIDInt)  ,
   _cpl            (cpl)         ,
   _mesh (cpl->meshGet()),
   _writer (nullptr),
   _id_writer_var_send(nullptr),
   _id_writer_var_recv(nullptr),
   _id_writer_var_send_status(-1),
   _id_writer_var_recv_status(-1),
   _interpolationFunction  (NULL),
   _interpolationFunction_f(NULL),
   _interpolationFunction_p(NULL),
   _computed_tgt_bcast_enabled(0),
   _involved_src_bcast_enabled(0),
   _is_send_yet(0),
   _is_send_end_step(0),
   _is_recv_yet(0),
   _is_recv_end_step(0),
   _python_object(nullptr)
{
  _mesh = cpl->meshGet();
  _n_part = _mesh->getNPart();
  // _data_tgt.resize(_n_part,NULL);
  // _data_src.resize(_n_part,NULL);
  _data_tgt = (void **) malloc(sizeof(void *) * _n_part);
  _data_src = (void **) malloc(sizeof(void *) * _n_part);
  _sendBuffer = NULL;
  _recvBuffer = NULL;
  _dataTypeSize = 0;

  switch (_dataType) {
    case CWP_DOUBLE:
      _dataTypeSize = sizeof(double);
      break;
    case CWP_INT:
      _dataTypeSize = sizeof(int);
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "CWP_CHAR is not usable.\n");
  }


  if (visuStatus == CWP_STATUS_ON) {
    _writer = cpl->writerGet();
  }

  _interpolationType = CWP_INTERPOLATION_DEFAULT;

  if (_writer != nullptr) {

    _id_writer_var_send = (int *) malloc (sizeof (int) * _nComponent);
    _id_writer_var_recv = (int *) malloc (sizeof (int) * _nComponent);
    for (int i = 0; i < _nComponent; i++) {
      _id_writer_var_send[i] = -1;
      _id_writer_var_recv[i] = -1;
    }

    // Set field location type
    PDM_writer_var_loc_t pdm_field_type = PDM_WRITER_VAR_ELEMENTS;

    if     ( _fieldLocation == CWP_DOF_LOCATION_CELL_CENTER) { 
      pdm_field_type = PDM_WRITER_VAR_ELEMENTS ;
    }
    else if ( _fieldLocation == CWP_DOF_LOCATION_NODE        ) {
      pdm_field_type = PDM_WRITER_VAR_VERTICES;
    }
    else if ( _fieldLocation == CWP_DOF_LOCATION_USER        ) {
      pdm_field_type = PDM_WRITER_VAR_VERTICES; // vertices on other geometry, PDM_WRITER_VAR_PARTICLES not implemented
    }

    // Create send variables
    if (_exchangeType == CWP_FIELD_EXCH_SEND ||
        _exchangeType == CWP_FIELD_EXCH_SENDRECV) {
      std::string prefix = "s";

      for (int i_comp = 0; i_comp < _nComponent; i_comp++) {

        std::ostringstream num;
        num << (i_comp + 1);

        std::string fieldName = prefix + "_" + _fieldID + num.str();

        _id_writer_var_send[i_comp] = PDM_writer_var_create(_writer,
                                                            PDM_WRITER_ON,
                                                            PDM_WRITER_VAR_SCALAR,
                                                            pdm_field_type,
                                                            fieldName.c_str());

        // printf("WriterFieldCreate - send: '%s' %d\n",fieldName.c_str(), _id_writer_var_send[i_comp]);
        // fflush(stdout);
      }

      // Create variable to tag if field has been exchanged (always time-dependent)
      std::string fieldName = prefix + "_" + _fieldID;
      std::string fieldStatusName = fieldName + "_status";

      _id_writer_var_send_status = PDM_writer_var_create(_writer,
                                                         PDM_WRITER_ON,
                                                         PDM_WRITER_VAR_SCALAR,
                                                         pdm_field_type,
                                                         fieldStatusName.c_str());
    }

    // Create receive variables
    if (_exchangeType == CWP_FIELD_EXCH_RECV ||
        _exchangeType == CWP_FIELD_EXCH_SENDRECV) {
      std::string prefix = "r";

      for (int i_comp = 0; i_comp < _nComponent; i_comp++) {

        std::ostringstream num;
        num << (i_comp + 1);

        std::string fieldName = prefix + "_" + _fieldID + num.str();

        _id_writer_var_recv[i_comp] = PDM_writer_var_create(_writer,
                                                            PDM_WRITER_ON,
                                                            PDM_WRITER_VAR_SCALAR,
                                                            pdm_field_type,
                                                            fieldName.c_str());

        // printf("WriterFieldCreate - recv: '%s' %d\n",fieldName.c_str(), _id_writer_var_recv[i_comp]);
        // fflush(stdout);
      }

      // Create variable to tag if field has been computed/exchanged
      std::string fieldName = prefix + "_" + _fieldID;
      std::string fieldStatusName = fieldName + "_status";

      _id_writer_var_recv_status = PDM_writer_var_create(_writer,
                                                         PDM_WRITER_ON,
                                                         PDM_WRITER_VAR_SCALAR,
                                                         pdm_field_type,
                                                         fieldStatusName.c_str());
    }

  }

}



Field::~Field()
{
  // _data_tgt.clear();
  // _data_src.clear();
  free(_data_tgt);
  free(_data_src);

  if (_id_writer_var_send != nullptr) {
    free (_id_writer_var_send);
  }
  if (_id_writer_var_recv != nullptr) {
    free (_id_writer_var_recv);
  }

  if (_sendBuffer != NULL) {
    free(_sendBuffer);
  }
  if (_recvBuffer != NULL) {
    free(_recvBuffer);
  }
}


void Field::dataSet ( int i_part, const CWP_Field_map_t   map_type, void* data)
{
  if (map_type == CWP_FIELD_MAP_SOURCE) {
    _data_src[i_part] = data;
  }
  else if (map_type == CWP_FIELD_MAP_TARGET) {
    _data_tgt[i_part] = data;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Field::dataSet Error : unknoown data type.\n");
  }

}


/**
 * \brief Write the field send or receive variables
 *
 * \param [in] exch_type Exchange type
 *
 */

void 
Field::write 
(
 CWP_Field_exch_t exch_type
)
{
  assert(CWP_FIELD_EXCH_SENDRECV != exch_type);

  if (_writer == nullptr) {
    return;
  }

  if (_cpl->NStepGet()%_cpl->freqWriterGet() != 0) {
    return;
  }

  if (visuStatusGet() == CWP_STATUS_OFF) {
    return;
  }

  // Get data depending on exchange type
  int   *id_writer_var;
  int    id_writer_var_status;
  void **data_var = NULL;
  int    write_end_step = 0;

  if (exch_type == CWP_FIELD_EXCH_SEND) {
    id_writer_var        = _id_writer_var_send;
    id_writer_var_status = _id_writer_var_send_status;
    data_var             = _data_src;
    write_end_step       = _is_send_end_step;
  }
  else {
    id_writer_var        = _id_writer_var_recv;
    id_writer_var_status = _id_writer_var_recv_status;
    data_var             = _data_tgt;
    write_end_step       = _is_recv_end_step;
  }


  // Write
  double *comp_data = NULL;
  int n_dof_max = 0;
  for (int i_part = 0; i_part < _cpl->nPartGet(); i_part++) {
    int n_dof = 0;

    if (_fieldLocation == CWP_DOF_LOCATION_CELL_CENTER) {
      n_dof = _mesh->getPartNElts(i_part);
    }
    else if (_fieldLocation == CWP_DOF_LOCATION_NODE) {
      n_dof = _mesh->getPartNVertex(i_part);
    }
    else if (_fieldLocation == CWP_DOF_LOCATION_USER) {
      n_dof = _cpl->userTargetNGet(i_part);
    }

    n_dof_max = std::max(n_dof_max, n_dof);
  }

  double *_comp_data = (double *) malloc(sizeof(double) * n_dof_max);

  // Fill in array and write it
  for (int i_part = 0; i_part < _cpl->nPartGet(); i_part++) {

    double *_data_var = (double *) data_var[i_part];

    int n_dof = 0;

    if (_fieldLocation == CWP_DOF_LOCATION_CELL_CENTER) {
      n_dof = _mesh->getPartNElts(i_part);
    }
    else if (_fieldLocation == CWP_DOF_LOCATION_NODE) {
      n_dof = _mesh->getPartNVertex(i_part);
    }
    else if (_fieldLocation == CWP_DOF_LOCATION_USER) {
      n_dof = _cpl->userTargetNGet(i_part);
    }

    // Data
    for (int i_comp = 0; i_comp < _nComponent; i_comp++) {
      if (write_end_step || exch_type == CWP_FIELD_EXCH_RECV) {
        for (int i = 0; i < n_dof; i++) {
          _comp_data[i] = 123456789;
        }
        comp_data = _comp_data;
      }

      if (!write_end_step) {
        if (exch_type == CWP_FIELD_EXCH_SEND) {
          if (_storage == CWP_FIELD_STORAGE_INTERLACED && _nComponent > 1) {
            for (int i = 0; i < n_dof; i++) {
              _comp_data[i] = _data_var[_nComponent*i + i_comp];
            }
            comp_data = _comp_data;
          }
          else {
            // avoid an unnecessary copy
            comp_data = _data_var + i_comp*n_dof;
          }
        }
        else {
          const int *computed_target   = NULL;
          int        n_computed_target = 0;
          if (_cpl->has_mesh()) {
            computed_target   = _cpl->computedTargetsGet (_fieldID, i_part);
            n_computed_target = _cpl->nComputedTargetsGet(_fieldID, i_part);
          }

          if (_storage == CWP_FIELD_STORAGE_INTERLACED) {
            for (int i = 0; i < n_computed_target; i++) {
              int j = computed_target[i]-1;
              _comp_data[j] = _data_var[_nComponent*i + i_comp];
            }
          }
          else {
            for (int i = 0; i < n_computed_target; i++) {
              int j = computed_target[i]-1;
              _comp_data[j] = _data_var[i_comp*n_computed_target + i];
            }
          }
        }
      }

      PDM_writer_var_set(_writer,
                         id_writer_var[i_comp],
                         _cpl->idGeomWriterGet(_fieldLocation),
                         i_part,
                         comp_data);

    }


    // Status
    if (write_end_step) {
      // Current field has NOT been exchanged
      for (int i = 0; i < n_dof; i++) {
        _comp_data[i] = -1;
      }
    }
    else {
      // Current field has been exchanged
      if (exch_type == CWP_FIELD_EXCH_SEND) {
        for (int i = 0; i < n_dof; i++) {
          _comp_data[i] = 0;
        }
      }
      else {
        const int *computed_target     = NULL;
        int        n_computed_target   = 0;
        const int *uncomputed_target   = NULL;
        int        n_uncomputed_target = 0;
        if (_cpl->has_mesh()) {
          computed_target     = _cpl->computedTargetsGet (_fieldID, i_part);
          n_computed_target   = _cpl->nComputedTargetsGet(_fieldID, i_part);
          uncomputed_target   = _cpl->uncomputedTargetsGet (_fieldID, i_part);
          n_uncomputed_target = _cpl->nUncomputedTargetsGet(_fieldID, i_part);
        }

        for (int i = 0; i < n_uncomputed_target; i++) {
          _comp_data[uncomputed_target[i]-1] = 1;
        }
        for (int i = 0; i < n_computed_target; i++) {
          _comp_data[computed_target[i]-1] = 0;
        }
      }
    }

    PDM_writer_var_set(_writer,
                       id_writer_var_status,
                       _cpl->idGeomWriterGet(_fieldLocation),
                       i_part,
                       _comp_data);

  }


  for (int i_comp = 0; i_comp < _nComponent; i_comp++) {
    PDM_writer_var_write    (_writer, id_writer_var[i_comp]);
    PDM_writer_var_data_free(_writer, id_writer_var[i_comp]);
  }
  PDM_writer_var_write    (_writer, id_writer_var_status);
  PDM_writer_var_data_free(_writer, id_writer_var_status);

  free(_comp_data);

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
Field::nDOFGet
(
 int i_part
)
{

  int n_dof = 0;

  if (_fieldLocation == CWP_DOF_LOCATION_CELL_CENTER) {
    n_dof = _mesh->getPartNElts(i_part);
  }
  else if (_fieldLocation == CWP_DOF_LOCATION_NODE) {
    n_dof = _mesh->getPartNVertex(i_part);
  }
  else if (_fieldLocation == CWP_DOF_LOCATION_USER) {
    n_dof = _cpl->userTargetNGet(i_part);
  }

  return n_dof;
}


}

/**
 * \endcond
 */
