#ifndef __OLDCOUPLING_H__
#define __OLDCOUPLING_H__
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

#include <string>
#include <map>
#include <vector>

#include <fvmc_locator.h>
#include <fvmc_writer.h>
#include <fvmc_nodal.h>
#include "cwipi.h"

namespace cwipi {

  class ApplicationProperties;

  class oldMesh;

  class LocationToDistantMesh;

  class LocationToLocalMesh;

  class oldCoupling {

  public:

    oldCoupling(const std::string& name,
             const cwipi_coupling_type_t couplingType,
             const ApplicationProperties& localApplicationProperties,
             const ApplicationProperties& coupledApplicationProperties,
             const int entitiesDim,
             const double tolerance,
             const cwipi_solver_type_t solverType,
             const int    outputFrequency,
             const char  *outputFormat,
             const char  *outputFormatOption,
             const int nb_Locations);

    virtual ~oldCoupling();

    void defineMesh(const int nVertex,
                    const int nElement,
                    double coordinates[],
                    int connectivity_index[],
                    int connectivity[],
                    int order = -1);

    void defineMesh(fvmc_nodal_t* fvmc_nodal);

    void hoOrderingSet (const cwipi_element_t t_elt,
                        const int n_nodes,
                        const int *ordering);

    void hoOrderingFromRefEltSet (const cwipi_element_t t_elt,
                                  const int n_nodes,
                                  const double *coords);

    void defineMeshAddPolyhedra(const int n_element,
                                int face_index[],
                                int cell_to_face_connectivity[],
                                const int nFace,
                                int face_connectivity_index[],
                                int face_connectivity[]);

    cwipi_exchange_status_t exchange(const char    *exchange_name,
                                     const int      stride,
                                     const int      time_step,
                                     const double   time_value,
                                     const char    *sending_field_name,
                                     const double  *sending_field,
                                     char          *receiving_field_name,
                                     double        *receiving_field,
                                     void          *ptFortranInterpolationFct);

    void issend(const char    *exchangeName,
                const int     tag,
                const int     stride,
                const int     timeStep,
                const double  timeValue,
                const char    *sendingFieldName,
                const double  *sendingField,
                void          *ptFortranInterpolationFct,
                int           *request);

    void waitIssend(int request);

    void irecv(const char    *exchangeName,
               const int     tag,
               const int     stride,
               const int     timeStep,
               const double  timeValue,
               const char    *receivingFieldName,
               const double  *receivingField,
               int           *request);

    void waitIrecv(int request);    

    void updateLocation();

    void setLocationIndex(const int index);
    void openLocationFile(const char *file, const char *moderwa);
    void closeLocationFile();
    void saveLocation();
    void loadLocation();

    void setPointsToLocate(const int n_points,
                           double coordinate[]);

    inline void set_interpolation_function(cwipi_interpolation_fct_t fct);

    inline void hoOptionsSet(const char* option, const char* value);

    inline void set_data_user(void *);

    inline void* get_data_user(void);

    inline void set_interpolation_function_f(void* fct);

    inline void set_ho_interpolation_function(cwipi_user_interp_ho_fct_t fct);

    inline void set_ho_interpolation_function_f(void* fct);

    inline int getNNotlocatedPoint() const;

    inline const int *getNotlocatedPoint() const;

    inline int getNLocatedPoint() const;

    inline const int *getLocatedPoint() const;

    inline const int *getLocatedPointsDistribution() const;

    inline const int *getDistantLocation() const;

    inline int getNDistantRank() const;

    inline const int *getDistantDistribution() const;

    inline const float *getDistantDistance() const;

    inline const float *distance() const;

    inline int getNDistantPoint() const;

    inline const double *getDistantPointCoordinates() const;

    inline const int *getDistantBarycentricCoordinatesIndex() const;

    inline const double *getDistantBarycentricCoordinates() const;

    ///
    /// \brief Set the type of information to be exchanged at the location
    ///
    ///   @param [in] information
    ///

    inline void setInfo(cwipi_located_point_info_t info);

    void locate();

    ///
    /// \brief Get distant element that contain located point
    ///

    inline const int *getElementContaining() const;

    ///
    /// \brief Get list of number of vertices of containing element
    ///

    inline const int *getElementContainingNVertex() const;

    ///
    /// \brief Get connectivity of element that contains each located point
    ///

    inline const int *getElementContainingVertex() const;

    ///
    /// \brief Get vertices coordinates of the element that contains each located point
    ///

    inline const double *getElementContainingVertexCoords() const;

    ///
    /// \brief Get barycentric coordinates for the element that contains each located point
    ///

    inline const double *getElementContainingBarycentricCoordinates() const;

    ///
    /// \brief Get MPI rank that contains each located point
    ///

    inline const int *getElementContainingMPIrank() const;

    ///
    /// \brief Exchange field on vertices of cells that contain each located points
    ///

    void exchangeCellVertexFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride);

    ///
    /// \brief Exchange field on cells that contain each located points
    ///

    void exchangeCellCenterFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride);

  private:

    oldCoupling();

    oldCoupling &operator=(const oldCoupling &other);

    std::vector<double> &  _extrapolate(double *cellCenterField, const int stride);

    //
    // TODO: Dans l'avenir créer une fabrique abstraite qui permet de definir differentes methodes d'interpolation

    void _interpolate(double *vertexField,
                      std::vector<double>& interpolatedField,
                      const int stride);

    void _interpolate1D(double *vertexField,
                        std::vector<double>& interpolatedField,
                        const int stride);

    void _interpolate2D(double *vertexField,
                        std::vector<double>& interpolatedField,
                        const int stride);

    void _interpolate3D(double *vertexField,
                        std::vector<double>& interpolatedField,
                        const int stride);

    void _fieldsVisualization(const char *exchangeName,
                              const int stride,
                              const int timeStep,
                              const double timeValue,
                              const char  *sendingFieldName,
                              const void *sendingField,
                              const char  *receivingFieldName,
                              const void *receivingField);

    void _initVisualization();

    void _createCouplingComm();

  private:
    const std::string               _name;
    const cwipi_coupling_type_t     _couplingType;
    const ApplicationProperties&    _localApplicationProperties;
    const ApplicationProperties&    _coupledApplicationProperties;
    const int                       _entitiesDim;
    const double                    _tolerance;
    const cwipi_solver_type_t       _solverType;
    const std::string               _outputFormat;
    const std::string               _outputFormatOption;
    const int                       _outputFrequency;

  private:
    fvmc_writer_t        *_fvmWriter;
    MPI_Comm             _couplingComm;
    int                  *_rankList;
    int                  _nRankList;
    MPI_Comm             _mergeComm;
    MPI_Comm             _fvmComm;
    int                  _coupledApplicationNRankCouplingComm;
    int                  _coupledApplicationBeginningRankCouplingComm;
    bool                 _isCoupledRank;

    cwipi_interpolation_fct_t _interpolationFct;
    cwipi_user_interp_ho_fct_t _ho_interpolationFct;
    
    void * _interpolationFct_f;
    void * _ho_interpolationFct_f;
    bool                 _toLocate;
    oldMesh                *_supportMesh;
    MPI_File _locationsFile;
    size_t _locationsFile_position;

    std::vector<LocationToDistantMesh *> &_tablelocationToDistantMesh;
    std::vector<LocationToLocalMesh *> &_tablelocationToLocalMesh;

    LocationToDistantMesh *_locationToDistantMesh;
    LocationToLocalMesh *_locationToLocalMesh;

  private:
    std::vector<double> *_tmpVertexField;  // Evite une allocation a chaque extrapolation
    std::vector<double> *_tmpDistantField; //TODO: Fusionner _tmpDistantField utiliser pour exchange
                                           // et les comm asynchrones
    std::map<int, std::vector<double> * > &_tmpDistantFieldsIssend; //TODO: temporaire : A revoir lors
                                                                    // de la restructuration
    std::map<int, const double * > &_tmpLocalFieldsIrecv;
    std::map<int, std::string > &_tmpExchangeNameIrecv;
    std::map<int, int > &_tmpStrideIrecv;
    std::map<int, int > &_tmpTimeStepIrecv;
    std::map<int, double > &_tmpTimeValueIrecv;
    std::map<int, std::string > &_tmpFieldNameIrecv;
    std::vector<float> _distance;

    void *_data_user;

    int _optBboxStep; 
  };

}

#endif //__COUPLING_PROPERTIES_H__
