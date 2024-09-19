#ifndef __SPATIAL_INTERP_H__
#define __SPATIAL_INTERP_H__
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

#include <cmath>

#include "mesh.hxx"
#include "field.hxx"
#include "codeProperties.hxx"
#include "coupling.hxx"
#include "pdm_part_to_part.h"
#include "pdm_writer.h"

/**
 * \cond
 */
namespace cwipi {

// A conserver et à renommer !!!

  typedef enum {
    CWP_SPATIAL_INTERP_AT_SEND,
    CWP_SPATIAL_INTERP_AT_RECV
  } CWP_SpatialInterp_time_t;

// A supprimer !!!

  struct target_data {
    int          lnum    ;
    int          origin_part     ;
    CWP_g_num_t  closest_elt_gnum;
    int          origin_proc     ;
    int          closest_elt_part;
    int          l_num_origin    ;
    double       projectedX       ;
    double       projectedY       ;
    double       projectedZ       ;
    double       distance       ;
  };


  typedef enum {
    SPATIAL_INTERP_EXCH_SEND,
    SPATIAL_INTERP_EXCH_RECV
  } SpatialInterpExchDirection;


  class Mesh;
  class Field;
  class Visu;

  /**
   *
   * \class SpatialInterp spatialInterp.hxx "spatialInterp.hxx"
   * \brief SpatialInterp algorithm
   *
   *  This class computes the spatial interpolation weights of points cloud into a mesh and
   *  builds a communication graph to transmit interpolated fieldsDouble on
   *  points cloud from fieldsDouble defined on the mesh.
   *
   */

  class SpatialInterp {

  public:

    /**
     * \brief Constructor
     *
     */

    SpatialInterp();

    /**
     * \brief Destructor
     *
     */

    virtual ~SpatialInterp();

    virtual void 
    init (
      Coupling           *coupling, 
      CWP_Dof_location_t localCodeDofLOcation,
      CWP_Dof_location_t cplCodeDofLOcation,
      SpatialInterpExchDirection exchDirection );

    virtual void clear();

    virtual void weightsCompute()  = 0; // A changer de nom
    
    virtual void interpolate (Field* referenceField, double **buffer) = 0;

    virtual void issend(Field* referenceField);

    virtual void irecv(Field* recevingField);

    virtual void waitIssend(Field* referenceField);

    virtual void waitIrecv(Field* recevingField);

    /**
     *
     * \brief Return the number of uncomputed targets
     *
     * \return                Number of uncomputed targets
     *
     */

    int
    nUncomputedTargetsGet(int i_part) const;

    /**
     *
     * \brief Return uncomputed targets
     *
     * \return                Uncomputed targets
     *
     */

    const int *
    uncomputedTargetsGet(int i_part) const;

    /**
     *
     * \brief Return the number of computed targets
     *
     * \return                Number of computed targets
     */

    int
    nComputedTargetsGet(int i_part) const;

    /**
     *
     * \brief Return computed targets
     *
     *
     * \return                Computed targets
     *
     */

    const int *
    computedTargetsGet(int i_part) const;

    int
    nInvolvedSourcesGet(int i_part) const;

    const int *
    involvedSourcesGet(int i_part) const;


    // Get part_to_part

    inline PDM_part_to_part_t *
    ptp_get()
    {
      return _ptsp;
    }

    // Get data related to weights

    inline int **
    weights_idx_get()
    {
      return _weights_idx;
    }

    inline double **
    weights_get()
    {
      return _weights;
    }

  protected:

    Coupling                   *_cpl;
    Mesh                       *_mesh;                  /*!< Interface Mesh */

    // Visu                       *_visu;                  /*!< Visualization object */
    CodeProperties             *_localCodeProperties;   
    CodeProperties             *_coupledCodeProperties; 

    CWP_Dof_location_t        _localCodeDofLocation;     /*!< Type of points cloud treated by this mapping instance (cell centers, vertices or user defined) */
    CWP_Dof_location_t        _coupledCodeDofLocation;   /*!< Type of points cloud treated by this mapping instance (cell centers, vertices or user defined) */

    SpatialInterpExchDirection  _exchDirection;  /*!< Spatial interpolation (for both codes are local case) */
   
    PDM_part_to_part_t *_ptsp;  /*!< Exchange protocol between src and target */
    int* _src_n_gnum;                      /*!< Number of source element by part (used by _ptsp) */
    int* _tgt_n_gnum;                      /*!< Number of target element by part (used by _ptsp) */
    const PDM_g_num_t** _src_gnum;                /*!< Global number of source element (used by _ptsp) */
    const PDM_g_num_t** _tgt_gnum;                /*!< Global number of target element (used by _ptsp) */
    
    CWP_SpatialInterp_time_t     _interpolation_time; /*!< Interpolation time : before or after exchange field */

    int  _nPart; /*!< Mesh partition number                                                    */
    int  _cplNPart;  /*!< Coupled mesh partition number                                                    */

    int _rootRankUnionComm   ;
    int _cplRootRankUnionComm;

    int _rootRankCplComm   ;
    int _cplRootRankCplComm;

    MPI_Comm _cplComm;
    PDM_MPI_Comm _pdmCplComm;
    MPI_Comm _unionComm;
    PDM_MPI_Comm _pdmUnionComm;

    MPI_Comm _localComm;          // Processus involved in the coupling for the local code

    int *_n_elt_weights;
    int **_weights_idx;
    double **_weights;

    int *_n_computed_tgt;
    int **_computed_tgt;

    int *_n_uncomputed_tgt;
    int **_uncomputed_tgt;

    int *_n_involved_sources_tgt;
    int **_involved_sources_tgt;

    std::vector <double **> _send_buffer;   /*!< Send buffer (size = n_field) */
    std::vector <double **> _recv_buffer;   /*!< Recv buffer (size = n_field) */
    std::vector <int>       _send_request;  /*!< Send request (size = n_field) */
    std::vector <int>       _recv_request;  /*!< Recv request (size = n_field) */


  // A conserver ou supprimer 
  protected:
    /* code Properties */
    int _id;
    int _id_cpl;
    string coupledName;
    string localName;


    int  _nPart_cpl                           ;  /*!< Coupled code mesh partition number                                       */


    /* informations about MPI process (rank) */
    int cplComm_rank;       // Rank in cplComm
    int cplComm_size;       // Size of cplComm
    int localComm_size;     // Size of localComm
    int localComm_size_cpl; // Size of localComm of the coupled code

  };
}

#endif //__SPATIAL_INTERP_H__
