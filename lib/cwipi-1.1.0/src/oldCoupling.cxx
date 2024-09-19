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

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstring>
#include <sstream>

#include <bftc_error.h>
#include <bftc_file.h>
#include <bftc_mem.h>

#include "fvmc_parall.h"

#include "oldCoupling.hxx"
#include "oldCoupling_i.hxx"

#include "oldMesh.hxx"
#include "applicationProperties.hxx"

#include "solve_ax_b_4.h"
#include "quickSort.h"
#include "cwipi.h"
#include "cwipi_config.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux) &&  !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
extern "C" {
  void PROCF(callfortinterpfct, CALLFORTINTERPFCT)
  ( int *entities_dim,
    int *n_local_vertex,
    int *n_local_element,
    int *n_local_polhyedra,
    int *n_distant_point,
    double *local_coordinates,
    int *local_connectivity_index,
    int *local_connectivity,
    int *local_polyhedra_face_index,
    int *local_polyhedra_cell_to_face_connectivity,
    int *local_polyhedra_face_connectivity_index,
    int *local_polyhedra_face_connectivity,
    double *distant_points_coordinates,
    int *distant_points_location,
    float *distant_points_distance,
    int *distant_points_barycentric_coordinates_index,
    double *distant_points_barycentric_coordinates,
    int *data_dimension,
    int *solver_type,
    double *local_field,
    double *distant_field,
    void *ptFortranInterpolationFct
  );
}

extern "C" {
  void PROCF(callforthointerpfct, CALLFORTHOINTERPFCT)
  ( int *entities_dim,
    int *order,
    int *n_local_vertex,
    int *n_local_element,
    int *n_local_polhyedra,
    int *n_distant_point,
    double *local_coordinates,
    int *local_connectivity_index,
    int *local_connectivity,
    int *local_polyhedra_face_index,
    int *local_polyhedra_cell_to_face_connectivity,
    int *local_polyhedra_face_connectivity_index,
    int *local_polyhedra_face_connectivity,
    double *distant_points_coordinates,
    int *distant_points_location,
    float *distant_points_distance,
    int *distant_points_weights_index,
    double *distant_points_weights,
    double *distant_points_uvw,
    int *data_dimension,
    int *solver_type,
    double *local_field,
    double *distant_field,
    void *ptFortranInterpolationFct
  );
}
#endif

namespace cwipi {

  //TODO: Faire une factory sur le type de couplage
  //TODO: Voir l'utilite de _fvmComm (A supprimer ?)

  oldCoupling::oldCoupling(const std::string& name,
                     const cwipi_coupling_type_t couplingType,
                     const ApplicationProperties& localApplicationProperties,
                     const ApplicationProperties& coupledApplicationProperties,
                     const int entitiesDim,
                     const double tolerance,
                     const cwipi_solver_type_t solverType,
                     const int    outputFrequency,
                     const char  *outputFormat,
                     const char  *outputFormatOption,
                     const int nbLocations)
 :_name(name), _couplingType(couplingType),
  _localApplicationProperties(localApplicationProperties),
  _coupledApplicationProperties(coupledApplicationProperties),
  _entitiesDim(entitiesDim),_tolerance(tolerance), _solverType(solverType),
  _outputFormat(outputFormat), _outputFormatOption(outputFormatOption),
  _outputFrequency(outputFrequency), _fvmWriter(NULL),
  _tablelocationToDistantMesh(*new std::vector<LocationToDistantMesh *>(nbLocations)),
  _tablelocationToLocalMesh(*new std::vector<LocationToLocalMesh *>(nbLocations)),
  _tmpDistantFieldsIssend(*new std::map<int, std::vector<double> * > ()),
  _tmpLocalFieldsIrecv(*new std::map<int, const double * > ()),
  _tmpExchangeNameIrecv(*new  std::map<int, std::string > ()),
  _tmpStrideIrecv(*new  std::map<int, int > ()),
  _tmpTimeStepIrecv(*new  std::map<int, int > ()),
  _tmpTimeValueIrecv(*new  std::map<int, double > ()),
  _tmpFieldNameIrecv(*new  std::map<int, std::string > ())

  {
    _tmpVertexField = NULL;
    _tmpDistantField = NULL;
    _supportMesh = NULL;
    _couplingComm = MPI_COMM_NULL;
    _rankList = NULL;
    _coupledApplicationNRankCouplingComm = -1;
    _coupledApplicationBeginningRankCouplingComm = -1;
    _mergeComm = MPI_COMM_NULL;
    _fvmComm = MPI_COMM_NULL;
    _interpolationFct = NULL;
    _interpolationFct_f = NULL;
    _ho_interpolationFct = NULL;
    _ho_interpolationFct_f = NULL;
    _toLocate = true;
    _isCoupledRank = false;
    _locationsFile_position = 0;
    _data_user = NULL;
    _optBboxStep = 10;


    //
    // Create coupling comm

    _createCouplingComm();

    //
    //

    for (int i = 0; i<nbLocations; i++) {
      _tablelocationToDistantMesh[i] = new LocationToDistantMesh(_isCoupledRank,
                                                                 _couplingType,
                                                                 _localApplicationProperties);
      _tablelocationToLocalMesh[i] = new LocationToLocalMesh(_solverType,
                                                             _optBboxStep,
                                                             _tolerance,
                                                             _couplingComm,
                                                             _coupledApplicationNRankCouplingComm,
                                                             _coupledApplicationBeginningRankCouplingComm,
                                                             _isCoupledRank,
                                                             _entitiesDim,
                                                             _localApplicationProperties,
                                                             *_tablelocationToDistantMesh[i]);
    }
    _locationToDistantMesh = _tablelocationToDistantMesh[0];
    _locationToLocalMesh = _tablelocationToLocalMesh[0];


#ifndef NAN
    bftc_printf("Warning : NAN macro is undefined -> receiving checking deactivation\n");
#endif

  }

  std::vector<double> &  oldCoupling::_extrapolate(double *cellCenterField, const int stride)
  {
    if (_tmpVertexField == NULL)
      _tmpVertexField = new  std::vector<double>(_supportMesh->getNVertex()*stride,0.);

    std::vector<double> &vertexField = *_tmpVertexField;

    for (size_t i = 0; i < vertexField.size(); i++)
      vertexField[i] = 0.;

    // TODO: Faire l'allocation qu'une fois comme _tmpVertexField

    std::vector<double> volumeVertex(_supportMesh->getNVertex(), 0.);
    std::vector<bool> orphanVertex(_supportMesh->getNVertex(), true);

    assert(_supportMesh != NULL);

    const int nElts = _supportMesh->getNElts();
    const int nPoly = _supportMesh->getNPolyhedra();
    const int nStandardElt = nElts - nPoly;

    const int *eltConnectivityIndex = _supportMesh->getEltConnectivityIndex();
    const int *eltConnectivity      = _supportMesh->getEltConnectivity();
    const std::vector<double>& cellVolume       = _supportMesh->getVolume();

    for (int i = 0; i < nStandardElt; i++) {
      const int nEltVertex = eltConnectivityIndex[i+1] - eltConnectivityIndex[i];
      const int index = eltConnectivityIndex[i];
      for (int j = 0; j < nEltVertex; j++) {
        int vertex = eltConnectivity[index+j] - 1;
        volumeVertex[vertex] += cellVolume[i];
        if (volumeVertex[vertex] < 0.)
          std::cout << "Volume : " << i << " " << cellVolume[i] << " " << volumeVertex[vertex] << std::endl;
        for (int k = 0; k < stride; k++)
          vertexField[stride*vertex+k] += cellCenterField[stride*i+k] * cellVolume[i];

        orphanVertex[vertex] = false;
      }
    }

    if (nPoly > 0) {
      std::vector<int> vertexPoly;
      vertexPoly.reserve(30);

      const int *polyhedraFaceIndex = _supportMesh->getPolyhedraFaceIndex();
      const int *polyhedraCellToFaceConnectivity = _supportMesh->getPolyhedraCellToFaceConnectivity();
      const int *polyhedraFaceConnectivityIndex = _supportMesh->getPolyhedraFaceConnectivityIndex() ;
      const int *polyhedraFaceConnectivity = _supportMesh->getPolyhedraFaceConnectivity();

      for (int i = 0; i < nPoly; i++) {
        int nFacePolyhedra = polyhedraFaceIndex[i+1] - polyhedraFaceIndex[i];
        int faceIndex = polyhedraCellToFaceConnectivity[i];
//        int nVertexFace = 0;
        for (int j = 0; j < nFacePolyhedra; j++) {
          int iface = polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
          int nVertexLocFace = polyhedraFaceConnectivityIndex[iface+1] - polyhedraFaceConnectivityIndex[iface];
//          nVertexFace += nVertexLocFace;
          int vertexIndex = polyhedraFaceConnectivityIndex[iface];
          for (int k = 0; k < nVertexLocFace; k++) {
            if (vertexPoly.capacity() <= vertexPoly.size()) {
              int capacity = 2*vertexPoly.capacity();
              vertexPoly.reserve(capacity);
            }
            vertexPoly.push_back(polyhedraFaceConnectivity[vertexIndex+k]);
          }

          quickSort(&vertexPoly[0], 0, vertexPoly.size()-1, NULL);

          int ivertex = -1;

          for (size_t j1 = 0; j1 < vertexPoly.size(); j1++) {
            if (ivertex < vertexPoly[j1]) {
              ivertex = vertexPoly[j1];
              volumeVertex[ivertex - 1] += cellVolume[i];
              if (volumeVertex[ivertex - 1] < 0.)
                std::cout << "Volume : " << i << " " << cellVolume[i] << " " << volumeVertex[ivertex - 1] << std::endl;
              for (int k = 0; k < stride; k++)
                vertexField[stride*(ivertex - 1)+k]  += cellCenterField[stride*i+k] * cellVolume[i];
              orphanVertex[ivertex - 1] = false;
            }
          }
        }
      }
    }

    for (int i = 0; i < _supportMesh->getNVertex(); i++) {
      if (orphanVertex[i] == true) {
        std::cout << "Vertex : " << i+1 << " is not connected to the connectivity !" << std::endl;
        vertexField[i] = 0.;
      }
      else
        for (int k = 0; k < stride; k++)
          vertexField[stride*i+k] /= volumeVertex[i];
    }

    return vertexField;
  }

  oldCoupling::~oldCoupling()
  {
#if defined(DEBUG) && 0
    std::cout << "destroying '" << _name << "' coupling" << std::endl;
#endif

    typedef std::map <int, std::vector<double> * >::iterator Iterator1;
    for (Iterator1 p = _tmpDistantFieldsIssend.begin();
         p != _tmpDistantFieldsIssend.end(); p++) {
      if (p->second != NULL)
        delete p->second;
    }

    delete [] _rankList;

    _distance.clear();

    _tmpDistantFieldsIssend.clear();
    delete &_tmpDistantFieldsIssend;

    _tmpLocalFieldsIrecv.clear();
    delete &_tmpLocalFieldsIrecv;

    _tmpExchangeNameIrecv.clear();
    delete &_tmpExchangeNameIrecv;

    _tmpFieldNameIrecv.clear();
    delete &_tmpFieldNameIrecv;

    _tmpStrideIrecv.clear();
    delete &_tmpStrideIrecv;

    _tmpTimeStepIrecv.clear();
    delete &_tmpTimeStepIrecv;

    _tmpTimeValueIrecv.clear();
    delete &_tmpTimeValueIrecv;

    if (_isCoupledRank ) {

      delete _tmpVertexField;

      delete _tmpDistantField;

      delete _supportMesh;

      for (size_t i = 0; i < _tablelocationToDistantMesh.size(); i++)
        delete _tablelocationToDistantMesh[i]; // libère chaque _tablelocationToDistantMesh
      _tablelocationToDistantMesh.clear();     // vide le vecteur
      delete &_tablelocationToDistantMesh;     // pour libère le new std::vector du constructeur

      for (size_t i = 0; i < _tablelocationToLocalMesh.size(); i++)
        delete _tablelocationToLocalMesh[i];
      _tablelocationToLocalMesh.clear();
      delete &_tablelocationToLocalMesh;

      MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
      if (oldFVMComm != MPI_COMM_NULL)
        MPI_Barrier(oldFVMComm);
      fvmc_parall_set_mpi_comm(_fvmComm);

      if (_fvmWriter != NULL)
        fvmc_writer_finalize(_fvmWriter);

      if (_fvmComm != MPI_COMM_NULL)
        MPI_Barrier(_fvmComm);
      fvmc_parall_set_mpi_comm(oldFVMComm);

    }

    if (_mergeComm != MPI_COMM_NULL)
      MPI_Comm_free(&_mergeComm);

    if (_couplingComm != MPI_COMM_NULL)
      MPI_Comm_free(&_couplingComm);

    if (_fvmComm != MPI_COMM_NULL)
      MPI_Comm_free(&_fvmComm);

  }

  void oldCoupling::_interpolate(double *referenceField,
                                 std::vector<double>& interpolatedField,
                                 const int stride)
  {

    //
    // For a cell center field : give the value of the located cell

    if (_solverType == CWIPI_SOLVER_CELL_CENTER) {
      const int nDistantPoint      = _locationToLocalMesh->getNLocatedDistantPoint() ;
      const int *distantLocation   = _locationToLocalMesh->getLocation();
      for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
        int iel = distantLocation[ipoint] - 1;
        for (int k = 0; k < stride; k++)
          interpolatedField[stride*ipoint+k] = referenceField[stride*iel+k];
      }
    }

    //
    // For a cell vertex field : interpolate with mean values

    else {

      double *dataField = referenceField;

      switch(_entitiesDim) {

      case 1 :
        _interpolate1D(dataField, interpolatedField, stride);
        break;

      case 2 :
        _interpolate2D(dataField, interpolatedField, stride);
        break;

      case 3 :
        _interpolate3D(dataField, interpolatedField, stride);
        break;

      default:
        bftc_error(__FILE__, __LINE__, 0, "'%i' bad entities dimension\n",_entitiesDim);
      }
    }
  }

  void oldCoupling::_interpolate1D(double *referenceVertexField,
                                   std::vector<double>& interpolatedField,
                                   const int stride)
  {
    const int nDistantPoint      =  _locationToLocalMesh->getNLocatedDistantPoint() ;
    const int *distantLocation   = _locationToLocalMesh->getLocation();
    const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
    const int *eltsConnec        = _supportMesh->getEltConnectivity();
    const double *barycentricCoordinates = _locationToLocalMesh->getBarycentricCoordinates();
    const int *barycentricCoordinatesIndex = _locationToLocalMesh->getBarycentricCoordinatesIndex();

    const int order = _supportMesh->getOrder();

    if (order == -1) {
      for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
        int iel = distantLocation[ipoint] - 1;
        double coef1 = barycentricCoordinates[barycentricCoordinatesIndex[ipoint]];
        double coef2 = barycentricCoordinates[barycentricCoordinatesIndex[ipoint]+1];
        int index = eltsConnecPointer[iel];
        int pt1 = eltsConnec[index] - 1;
        int pt2 = eltsConnec[index+1] - 1;

        for (int k = 0; k < stride; k++)
          interpolatedField[stride*ipoint+k] = coef1 * referenceVertexField[stride*pt1+k] +
            coef2 * referenceVertexField[stride*pt2+k];
      }
    }

    else {

      for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
        int iel = distantLocation[ipoint] - 1;

        const int *intern_connec =
          fvmc_nodal_get_internal_connec_elt (&(_supportMesh->getFvmNodal()), iel+1);

        const int n_node =  eltsConnecPointer[iel+1] - eltsConnecPointer[iel];

        double *target_field =  &(interpolatedField[0]) + stride * ipoint;

        const double *weight = barycentricCoordinates + barycentricCoordinatesIndex[ipoint];

        const double *src_field = referenceVertexField;

        for (int i = 0; i < n_node; i++) {

          int i_node = intern_connec[i] - 1;

          for (int j = 0; j < stride; j++) {
            target_field[j] += weight[i] * src_field[stride * i_node + j];
          }
        }

      }
    }
  }


  void oldCoupling::_interpolate2D (double *vertexField,
                                    std::vector<double>& interpolatedField,
                                    const int stride)
  {

    const int nDistantPoint      =  _locationToLocalMesh->getNLocatedDistantPoint() ;
    const int *distantLocation   = _locationToLocalMesh->getLocation();

    const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
    const int *eltsConnec        = _supportMesh->getEltConnectivity();

    const double *barycentricCoordinates = _locationToLocalMesh->getBarycentricCoordinates();
    const int *barycentricCoordinatesIndex = _locationToLocalMesh->getBarycentricCoordinatesIndex();

    const int order = _supportMesh->getOrder();

    if (order == -1) {

      for (int ipoint = 0; ipoint <nDistantPoint; ipoint++) {

        int iel = distantLocation[ipoint] - 1;
        int index = barycentricCoordinatesIndex[ipoint];
        int nSom = barycentricCoordinatesIndex[ipoint+1] - index;

        for (int k = 0; k < stride; k++)
          interpolatedField[stride*ipoint + k] = 0;

        for (int isom = 0; isom <  nSom; isom++) {
          for (int k = 0; k < stride; k++) {
            interpolatedField[stride*ipoint+k] +=
              vertexField[stride*(eltsConnec[eltsConnecPointer[iel]+isom]-1)+k]
              *barycentricCoordinates[index+isom];
          }
        }
      }
    }

    else {

      for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
        int iel = distantLocation[ipoint] - 1;

        const int *intern_connec =
          fvmc_nodal_get_internal_connec_elt (&(_supportMesh->getFvmNodal()), iel+1);

        const int n_node =  eltsConnecPointer[iel+1] - eltsConnecPointer[iel];

        double *target_field =  &(interpolatedField[0]) + stride * ipoint;

        const double *weight = barycentricCoordinates + barycentricCoordinatesIndex[ipoint];

        const double *src_field = vertexField;

        for (int i = 0; i < n_node; i++) {

          int i_node = intern_connec[i] - 1;

          for (int j = 0; j < stride; j++) {
            target_field[j] += weight[i] * src_field[stride * i_node + j];
          }
        }

      }
    }

  }


  void oldCoupling::_interpolate3D(double *vertexField,
                                   std::vector<double>& interpolatedField,
                                   const int stride)
  {
    const int nDistantPoint      =  _locationToLocalMesh->getNLocatedDistantPoint() ;
    const int *distantLocation   = _locationToLocalMesh->getLocation();
    const double *distantCoords  = _locationToLocalMesh->getPointCoordinates();

    const int *eltsConnecPointer = _supportMesh->getEltConnectivityIndex();
    const int *eltsConnec        = _supportMesh->getEltConnectivity();

    const std::vector<int> &polyEltsConnec =
      _supportMesh->getPolyhedraCellToVertexConnectivity();
    const std::vector<int> &polyEltsConnecPointer =
      _supportMesh->getPolyhedraCellToVertexConnectivityIndex();

    const int nStandardElt       = _supportMesh->getNElts() - _supportMesh->getNPolyhedra();
    const double *coords         = _supportMesh->getVertexCoords();
    double coeff[4];

    const double *barycentricCoordinates =
      _locationToLocalMesh->getBarycentricCoordinates();
    const int *barycentricCoordinatesIndex =
      _locationToLocalMesh->getBarycentricCoordinatesIndex();

    const bool useMeanValues = true;

    if (useMeanValues) {

      const int order = _supportMesh->getOrder();

      if (order == -1) {

        for (int ipoint = 0; ipoint <nDistantPoint; ipoint++) {
          int iel = distantLocation[ipoint] - 1;
          int index = barycentricCoordinatesIndex[ipoint];
          int nSom = barycentricCoordinatesIndex[ipoint+1] - index;

          for (int k = 0; k < stride; k++)
            interpolatedField[stride*ipoint + k] = 0;

          for (int isom = 0; isom <  nSom; isom++) {
            for (int k = 0; k < stride; k++) {
              if (iel < nStandardElt) {
                interpolatedField[stride*ipoint+k] +=
                  vertexField[stride*(eltsConnec[eltsConnecPointer[iel]+isom]-1)+k]
                  *barycentricCoordinates[index+isom];
              }
              else {
                interpolatedField[stride*ipoint+k] +=
                  vertexField[stride*(polyEltsConnec[polyEltsConnecPointer[iel - nStandardElt]+isom]-1)+k]
                  *barycentricCoordinates[index+isom];
              }
            }
          }
        }

      }

      else {

        for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
          int iel = distantLocation[ipoint] - 1;

          const int *intern_connec =
            fvmc_nodal_get_internal_connec_elt (&(_supportMesh->getFvmNodal()), iel+1);

          const int n_node =  eltsConnecPointer[iel+1] - eltsConnecPointer[iel];

          double *target_field =  &(interpolatedField[0]) + stride * ipoint;

          const double *weight = barycentricCoordinates + barycentricCoordinatesIndex[ipoint];

          const double *src_field = vertexField;

          for (int i = 0; i < n_node; i++) {

            int i_node = intern_connec[i] - 1;

            for (int j = 0; j < stride; j++) {
              target_field[j] += weight[i] * src_field[stride * i_node + j];
            }

          }

        }
      }
    }

    else {
      for (int ipoint = 0; ipoint < nDistantPoint; ipoint++) {
        int iel = distantLocation[ipoint] - 1;
        if (iel < nStandardElt) {
          int index = eltsConnecPointer[iel];
          int nVertex = eltsConnecPointer[iel+1]-index;

          for (int k = 0; k < stride; k++) {
            double a[4][4] = {{0., 0., 0., 0.},
                              {0., 0., 0., 0.},
                              {0., 0., 0., 0.},
                              {0., 0., 0., 0.}};
            double b[4] = {0., 0., 0., 0.};
            for (int i = 0; i < nVertex; i++) {
              int iVertex = eltsConnec[index+i]-1;
              double v_x = coords[3*iVertex];
              double v_y = coords[3*iVertex+1];
              double v_z = coords[3*iVertex+2];
              double v_f = vertexField[stride*iVertex+k];

              a[0][0] += v_x * v_x;
              a[0][1] += v_x * v_y;
              a[0][2] += v_x * v_z;
              a[0][3] += v_x;

              a[1][1] += v_y * v_y;
              a[1][2] += v_y * v_z;
              a[1][3] += v_y;

              a[2][2] += v_z * v_z;
              a[2][3] += v_z;

              a[3][3] += 1.;

              b[0] += v_x * v_f;
              b[1] += v_y * v_f;
              b[2] += v_z * v_f;
              b[3] += v_f;
            }
            a[1][0] = a[0][1];
            a[2][0] = a[0][2];
            a[3][0] = a[0][3];

            a[2][1] = a[1][2];
            a[3][1] = a[1][3];

            a[3][2] = a[2][3];

            if (solve_ax_b_4(a, b, coeff) == 0) {
              interpolatedField[stride*ipoint+k] = (coeff[0] * distantCoords[3*ipoint]
                                                    + coeff[1] * distantCoords[3*ipoint+1]
                                                    + coeff[2] * distantCoords[3*ipoint+2]
                                                    + coeff[3]);
            }
            else {
              interpolatedField[stride*ipoint+k] = vertexField[stride*nVertex+k];
              /* last encountered value */
            }
          }
        }
        else {

          const int *polyhedraFaceIndex = _supportMesh->getPolyhedraFaceIndex();
          const int *polyhedraCellToFaceConnectivity =
            _supportMesh->getPolyhedraCellToFaceConnectivity();
          const int *polyhedraFaceConnectivityIndex =
            _supportMesh->getPolyhedraFaceConnectivityIndex() ;
          const int *polyhedraFaceConnectivity =
            _supportMesh->getPolyhedraFaceConnectivity();

          int nFacePolyhedra = polyhedraFaceIndex[iel+1] - polyhedraFaceIndex[iel];
          int faceIndex = polyhedraCellToFaceConnectivity[iel];
          int nVertexFace = 0;

          for (int j = 0; j < nFacePolyhedra; j++) {
            int iface = polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
            nVertexFace += polyhedraFaceConnectivityIndex[iface+1] - polyhedraFaceConnectivityIndex[iface];
          }

          std::vector<int> vertexPoly(nVertexFace);
          for (int j = 0; j < nFacePolyhedra; j++) {
            int iface = polyhedraCellToFaceConnectivity[faceIndex+j] - 1;
            int nVertexLocFace = polyhedraFaceConnectivityIndex[iface+1] - polyhedraFaceConnectivityIndex[iface];
            int vertexIndex = polyhedraFaceConnectivityIndex[iface];
            for (int k = 0; k < nVertexLocFace; k++)
              vertexPoly.push_back(polyhedraFaceConnectivity[vertexIndex+k]);
          }
          quickSort(&vertexPoly[0], 0, vertexPoly.size()-1, NULL);

          int iVertex = -1;
          double v_x;
          double v_y;
          double v_z;
          double v_f;

          for(int k = 0; k < stride; k++) {
            double a[4][4] = {{0., 0., 0., 0.},
                              {0., 0., 0., 0.},
                              {0., 0., 0., 0.},
                              {0., 0., 0., 0.}};
            double b[4] = {0., 0., 0., 0.};

            for (size_t j = 0; j < vertexPoly.size(); j++) {
              if (iVertex < vertexPoly[j]) {
                iVertex = vertexPoly[j];
                v_x = coords[3*iVertex];
                v_y = coords[3*iVertex+1];
                v_z = coords[3*iVertex+2];
                v_f = vertexField[stride*iVertex+k];

                a[0][0] += v_x * v_x;
                a[0][1] += v_x * v_y;
                a[0][2] += v_x * v_z;
                a[0][3] += v_x;

                a[1][1] += v_y * v_y;
                a[1][2] += v_y * v_z;
                a[1][3] += v_y;

                a[2][2] += v_z * v_z;
                a[2][3] += v_z;

                a[3][3] += 1.;

                b[0] += v_x * v_f;
                b[1] += v_y * v_f;
                b[2] += v_z * v_f;
                b[3] += v_f;
              }
              a[1][0] = a[0][1];
              a[2][0] = a[0][2];
              a[3][0] = a[0][3];

              a[2][1] = a[1][2];
              a[3][1] = a[1][3];

              a[3][2] = a[2][3];

              if (solve_ax_b_4(a, b, coeff) == 0) {
                interpolatedField[stride*ipoint+k] = (coeff[0] * distantCoords[3*ipoint]
                                                      + coeff[1] * distantCoords[3*ipoint+1]
                                                      + coeff[2] * distantCoords[3*ipoint+2]
                                                      + coeff[3]);
              }
              else {
                interpolatedField[stride*ipoint+k] = v_f; /* last encountered value */
              }
            }
          }
        }
      }
    }
  }


  void oldCoupling::defineMesh(const int nVertex,
                            const int nElement,
                            double coordinates[],
                            int connectivity_index[],
                            int connectivity[],
                            int order)
  {


    if (_supportMesh  != NULL)
      bftc_error(__FILE__, __LINE__, 0, "coupling mesh is already created\n");

    if (_isCoupledRank)
      _supportMesh = new oldMesh(_fvmComm,
                                 _entitiesDim,
                                 nVertex,
                                 nElement,
                                 coordinates,
                                 connectivity_index,
                                 connectivity,
                                 order);
    else
      bftc_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning,"
                 " the coupling mesh must be defined only by the root rank\n");

  }


  void oldCoupling::defineMesh(fvmc_nodal_t* fvmc_nodal)
  {
    if (_supportMesh  != NULL)
      bftc_error(__FILE__, __LINE__, 0, "coupling mesh is already created\n");

    if (_isCoupledRank)
      _supportMesh = new oldMesh(_fvmComm,
                              fvmc_nodal);
    else
      bftc_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning,"
                 " the coupling mesh must be defined only by the root rank\n");

  }

 void oldCoupling::hoOrderingSet (const cwipi_element_t t_elt,
                                  const int n_nodes,
                                  const int *ordering)
 {
   if (_supportMesh == NULL) {
     bftc_error(__FILE__, __LINE__, 0, "Define a mesh before setting ho ordering\n");
   }

   _supportMesh->hoOrderingSet (t_elt, n_nodes, ordering);

 }

 void oldCoupling::hoOrderingFromRefEltSet (const cwipi_element_t t_elt,
                                            const int n_nodes,
                                            const double *coords)

 {
   if (_supportMesh == NULL) {
     bftc_error(__FILE__, __LINE__, 0, "Define a mesh before setting ho ordering\n");
   }

   _supportMesh->hoOrderingFromRefEltSet (t_elt, n_nodes, coords);

 }



  void oldCoupling::setPointsToLocate(const int    n_points,
                                   double coordinate[])
  {
    _toLocate = true;
    _locationToDistantMesh->setpointsToLocate(n_points, coordinate);
  }


  void oldCoupling::defineMeshAddPolyhedra(const int n_element,
                                        int face_index[],
                                        int cell_to_face_connectivity[],
                                        const int nFace,
                                        int face_connectivity_index[],
                                        int face_connectivity[])

  {
    if (_isCoupledRank) {
      if (_supportMesh  == NULL)
        bftc_error(__FILE__, __LINE__, 0, "No mesh to add elements\n");

      _supportMesh->addPolyhedra(n_element,
                                 face_index,
                                 cell_to_face_connectivity,
                                 nFace,
                                 face_connectivity_index,
                                 face_connectivity);
    }
    else
      bftc_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning,"
                 " the coupling mesh must be defined only by the root rank\n");

  }


  void oldCoupling::updateLocation()
  {
    if (_isCoupledRank) {
      _toLocate = true;
    }
    else
      bftc_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning,"
                 " updateLocation must be called only by the root rank\n");
  }

  void oldCoupling::setLocationIndex(const int index)
  {
    if (_isCoupledRank) {
      _locationToDistantMesh = _tablelocationToDistantMesh[index];
      _locationToLocalMesh = _tablelocationToLocalMesh[index];
    }
    else
      bftc_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning,"
                 " updateLocation must be called only by the root rank\n");
  }


  void oldCoupling::openLocationFile(const char *file, const char *moderwa)
  {

    int mode = MPI_MODE_CREATE+MPI_MODE_WRONLY;
    _locationsFile_position = 0;


    if (moderwa[0] == 'w')
      mode = MPI_MODE_CREATE+MPI_MODE_WRONLY;
    else if (moderwa[0] == 'r')
      mode = MPI_MODE_RDONLY;
    else if (moderwa[0] == 'a')
      mode = MPI_MODE_APPEND;
    else
      bftc_error(__FILE__, __LINE__, 0, "The file mode must be 'w' 'r' or 'a' for write read or append\n");

    MPI_File_open(_couplingComm, file, mode, MPI_INFO_NULL, &_locationsFile);


  }

  void oldCoupling::closeLocationFile()
  {
    MPI_File_close(&_locationsFile);
  }



  void oldCoupling::saveLocation()
  {
      printf("saveLocation\n");
    if (_toLocate)
      bftc_error(__FILE__, __LINE__, 0, "trying to save a not located locator");
    if (_isCoupledRank) {
      int currentRank, nbRank;
      size_t       il_size;
      unsigned char *copybuf;
      size_t *records_size;

      MPI_Status err;
      size_t locationToLocalMesh_size,  locationToDistantMesh_size;
      size_t il_position;

      MPI_Comm_rank(_couplingComm, &currentRank);
      MPI_Comm_size(_couplingComm, &nbRank);

      // calcul de la taille totale d'un enregistrement sur chaque proc
      int nLocPts = getNLocatedPoint();
      il_size = sizeof(int) + nLocPts*sizeof(float); // pour le tableau _distance et sa taille
      //il_size = 0;
      // taille de ce qu'il y a sauvegarder dans locationToLocalMesh
      locationToLocalMesh_size = _locationToLocalMesh->locationSize();
      il_size +=  locationToLocalMesh_size;
      // + taille de ce qu'il y a sauvegarder dans locationToDistantMesh
      locationToDistantMesh_size =  _locationToDistantMesh->locationSize();
      il_size += locationToDistantMesh_size ;


      // records_size de taille nbRank
      // va contenir les tailles des enregistrements sur les differents procs
      // chaque proc en a besoin pour caluler sa position relative dans le fichier

      BFTC_MALLOC(records_size, nbRank, size_t);
      records_size[currentRank] = il_size; // taille sur le proc courant
      MPI_Allgather((void *)&il_size,sizeof(size_t),MPI_CHAR,records_size,sizeof(size_t),MPI_CHAR, _couplingComm);

      // seul le proc 0 ecrit ce tableau records_size
      if (currentRank == 0) {
        MPI_File_write_at(_locationsFile, _locationsFile_position,
                          (void *)records_size,  nbRank * sizeof(size_t), MPI_CHAR, &err);
      }
      // increment du pointeur global de fichier
      _locationsFile_position += nbRank * sizeof(size_t);

      // remplissage du buffer d'ecriture copybuf
      BFTC_MALLOC(copybuf,il_size, unsigned char);

      il_position = 0;
      memcpy((void *)&copybuf[il_position],(void *)&nLocPts,sizeof(int));
      il_position += sizeof(int);

      memcpy((void *)&copybuf[il_position],(void *)&(_distance[0]),nLocPts*sizeof(float) );
      il_position += nLocPts*sizeof(float);

      _locationToLocalMesh->packLocation(&copybuf[il_position]);
      il_position += locationToLocalMesh_size;
      _locationToDistantMesh->packLocation(&copybuf[il_position]);

      // calcul de l'adresse de chaque proc en tenant compte de ce que vont prendre les autres
      il_position = _locationsFile_position;
      for (int i = 0; i < currentRank ; i++)
        il_position +=  records_size[i];

      // ecriture du tampon
      MPI_File_write_at_all(_locationsFile,  il_position, (void *) copybuf, il_size, MPI_CHAR, &err);


      // mise a jour du pointeur global sur le fichier
      // en tenant compte de tous les procs
      for (int i = 0; i<nbRank ; i++)
        _locationsFile_position += records_size[i];

      // liberation memoire
      BFTC_FREE(copybuf);
      BFTC_FREE(records_size);
      BFTC_FREE(copybuf);

    }
    else
      bftc_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning,"
                 " saveLocation must be called only by the root rank\n");
  }



  void oldCoupling::loadLocation()
  {
    if (_isCoupledRank) {
      int currentRank, nbRank;
      size_t il_position;
      unsigned char *copybuf;
      size_t *records_size;
      MPI_Status err;
      size_t il_size;

      MPI_Comm_rank(_couplingComm, &currentRank);
      MPI_Comm_size(_couplingComm, &nbRank);
      BFTC_MALLOC(records_size,nbRank, size_t);

      // le premier enregistrement contient la taille de chaque location
      // sur tous les procs, seul le proc 0 lit cet enregistrement et
      // fait un broacast pour le passer aux autres
      if (currentRank == 0) {
        il_position = _locationsFile_position;
        MPI_File_read_at(_locationsFile, il_position, (void *)records_size, nbRank * sizeof(size_t), MPI_CHAR, &err);
      }
      MPI_Bcast(records_size, nbRank*sizeof(size_t), MPI_CHAR, 0, _couplingComm);
      _locationsFile_position += nbRank * sizeof(size_t);

      // allocation d'un buffer pour la lecture
      il_size = records_size[currentRank];

      printf("records_size[currentRank] %ld\n",  records_size[currentRank]);

      BFTC_MALLOC(copybuf,il_size, unsigned char);

      // lecture du tampon
      il_position = _locationsFile_position; // pointeur global
      for (int i = 0; i < currentRank ; i++)
        il_position +=  records_size[i];     // ajout des tailles jusqu'au proc courant
      MPI_File_read_at_all(_locationsFile,  il_position, (void *) copybuf, il_size, MPI_CHAR, &err);

      // unpack des donnees

      il_position = 0;
      int nLocPts;
      il_position += fvmc_locator_unpack_elem((void *)&copybuf[il_position],(void *)&nLocPts,sizeof(int));
      if ((int) _distance.size() != nLocPts)
        _distance.resize(nLocPts);
      il_position += fvmc_locator_unpack_elem((void *)&copybuf[il_position],(void *)&(_distance[0]),nLocPts*sizeof(float));

      _locationToLocalMesh->unpackLocation(&copybuf[il_position]);
      size_t locationToLocalMesh_size = _locationToLocalMesh->locationSize();
      il_position += locationToLocalMesh_size;
      _locationToDistantMesh->unpackLocation(&copybuf[il_position]);
      // mise a jour du pointeur global sur le fichier
      for (int i = 0; i<nbRank ; i++)
        _locationsFile_position += records_size[i];

      BFTC_FREE(records_size);
      BFTC_FREE(copybuf);
      _locationToLocalMesh->setSupportMesh(_supportMesh,false);
      _toLocate = false;
      _initVisualization();
    }
    else
      bftc_error(__FILE__, __LINE__, 0, "for a coupling without parallel partitionning,"
                 " loadLocation must be called only by the root rank\n");
  }


  cwipi_exchange_status_t oldCoupling::exchange(const char    *exchangeName,
                                             const int     stride,
                                             const int     timeStep,
                                             const double  timeValue,
                                             const char    *sendingFieldName,
                                             const double  *sendingField,
                                             char          *receivingFieldName,
                                             double        *receivingField,
                                             void          *ptFortranInterpolationFct)


  {


    cwipi_exchange_status_t status = CWIPI_EXCHANGE_OK;

    const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

    int rootRank;

    if (!_isCoupledRank && sendingField != NULL )
      bftc_printf("Warning : sendingField != NULL, "
                  " only field defined by the root rank is sent\n");

    if (_isCoupledRank) {

      int currentRank;
      MPI_Comm_rank(_couplingComm, &currentRank);


      int lLocalName = _name.size() + 1;
      int lDistantName = 0;
      MPI_Status MPIStatus;

      if (currentRank == 0 ||
          currentRank == _coupledApplicationBeginningRankCouplingComm +
          _coupledApplicationNRankCouplingComm) {

        //
        // Check coupling name

        MPI_Sendrecv(&lLocalName,   1, MPI_INT,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     &lDistantName, 1, MPI_INT,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     _couplingComm, &MPIStatus);

        char *distantCouplingName =  (char *) malloc (sizeof(char) * (lDistantName));

        MPI_Sendrecv(const_cast <char*>(_name.c_str()),
                     lLocalName, MPI_CHAR,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     distantCouplingName, lDistantName, MPI_CHAR,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     _couplingComm, &MPIStatus);

        if (strcmp(_name.c_str(), distantCouplingName))
          bftc_error(__FILE__, __LINE__, 0, "'%s' '%s' bad synchronization point\n",
                     _name.c_str(),
                     distantCouplingName);
        free ( distantCouplingName);

        //
        // Check exchange name

        lLocalName = strlen(exchangeName)+1;
        MPI_Sendrecv(&lLocalName,   1, MPI_INT,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     &lDistantName, 1, MPI_INT,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     _couplingComm, &MPIStatus);

        char *distantExchangeName =  (char *) malloc (sizeof(char) * (lDistantName));

        MPI_Sendrecv(const_cast <char*>(exchangeName),
                     lLocalName, MPI_CHAR,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     distantExchangeName, lDistantName, MPI_CHAR,
                     _coupledApplicationBeginningRankCouplingComm, 0,
                     _couplingComm, &MPIStatus);

        if (strcmp(exchangeName, distantExchangeName))
          bftc_error(__FILE__, __LINE__, 0, "'%s' '%s' bad synchronization point\n",
                     exchangeName,
                     distantExchangeName);

        free ( distantExchangeName);
      }
    }

    //
    // Locate

    // TM mars 2014 : ajout du test if (_toLocate) avant le locate

    if (_toLocate)
      locate();

    //
    // Prepare data (interpolate, extrapolate...)

    if (_isCoupledRank) {

      const int nVertex                 = _supportMesh->getNVertex();
      const int nElts                   = _supportMesh->getNElts();
      const int nPoly                   = _supportMesh->getNPolyhedra();
      const int *localConnectivityIndex = _supportMesh->getEltConnectivityIndex();
      const int *localConnectivity      = _supportMesh->getEltConnectivity();
      const int *localPolyhedraFaceIndex = _supportMesh->getPolyhedraFaceIndex() ;
      const int *localPolyhedraCellToFaceConnectivity = _supportMesh->getPolyhedraCellToFaceConnectivity();
      const int *localPolyhedraFaceConnectivity_index = _supportMesh->getPolyhedraFaceConnectivityIndex();
      const int *localPolyhedraFaceConnectivity       = _supportMesh->getPolyhedraFaceConnectivity();

      const int nDistantPoint      =  _locationToLocalMesh->getNLocatedDistantPoint() ;
      const int *distantLocation   = _locationToLocalMesh->getLocation();
      const float *distantDistance   = _locationToLocalMesh->getDistance();
      const double *distantCoords  = _locationToLocalMesh->getPointCoordinates();

      const int nInteriorList     = _locationToDistantMesh->getNLocatedPoint();
      const double *barycentricCoordinates = _locationToLocalMesh->getBarycentricCoordinates();
      const int *barycentricCoordinatesIndex = _locationToLocalMesh->getBarycentricCoordinatesIndex();

      int lDistantField = stride * nDistantPoint;

      if (_tmpDistantField == NULL)
        _tmpDistantField = new std::vector<double> (lDistantField);

      std::vector<double>& tmpDistantField = *_tmpDistantField;
      if ((int) tmpDistantField.size() < lDistantField)
        tmpDistantField.resize(lDistantField);

      //
      // Interpolation

      if (sendingField != NULL) {

        assert(!(_interpolationFct != NULL && ptFortranInterpolationFct != NULL));

        //
        // Callback Fortran


        const int order = _supportMesh->getOrder();

        if (order == -1) {

          if (ptFortranInterpolationFct != NULL) {
#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
            PROCF(callfortinterpfct, CALLFORTINTERPFCT) (
                                                         const_cast <int *> (&_entitiesDim),
                                                         const_cast <int *> (&nVertex),
                                                         const_cast <int *> (&nElts),
                                                         const_cast <int *> (&nPoly),
                                                         const_cast <int *> (&nDistantPoint),
                                                         const_cast <double *> (_supportMesh->getVertexCoords()),
                                                         const_cast <int *> (localConnectivityIndex),
                                                         const_cast <int *> (localConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceIndex),
                                                         const_cast <int *> (localPolyhedraCellToFaceConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity_index),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity),
                                                         const_cast <double *> (distantCoords),
                                                         const_cast <int *> (distantLocation),
                                                         const_cast <float *> (distantDistance),
                                                         const_cast <int *> (barycentricCoordinatesIndex),
                                                         const_cast <double *> (barycentricCoordinates),
                                                         const_cast <int *> (&stride),
                                                         const_cast <int *> ((const int *) &_solverType),
                                                         const_cast <double *> (sendingField),
                                                         const_cast <double *> (&tmpDistantField[0]),
                                                         ptFortranInterpolationFct
                                                         );
#else
            fprintf(stderr,"CWIPI Error : Impossible to use a Fortran user interpolation function");
            abort();
#endif
          }

          //
          // Callback C

          else if (_interpolationFct != NULL)
            _interpolationFct(_entitiesDim,
                              nVertex,
                              nElts,
                              nPoly,
                              nDistantPoint,
                              _supportMesh->getVertexCoords(),
                              localConnectivityIndex,
                              localConnectivity,
                              localPolyhedraFaceIndex,
                              localPolyhedraCellToFaceConnectivity,
                              localPolyhedraFaceConnectivity_index,
                              localPolyhedraFaceConnectivity,
                              distantCoords,
                              distantLocation,
                              distantDistance,
                              barycentricCoordinatesIndex,
                              barycentricCoordinates,
                              stride,
                              _solverType,
                              sendingField,
                              &tmpDistantField[0]);

          //
          // Callback Fortran appele en C

          else if (_interpolationFct_f != NULL) {
#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
            PROCF(callfortinterpfct, CALLFORTINTERPFCT) (
                                                         const_cast <int *> (&_entitiesDim),
                                                         const_cast <int *> (&nVertex),
                                                         const_cast <int *> (&nElts),
                                                         const_cast <int *> (&nPoly),
                                                         const_cast <int *> (&nDistantPoint),
                                                         const_cast <double *> (_supportMesh->getVertexCoords()),
                                                         const_cast <int *> (localConnectivityIndex),
                                                         const_cast <int *> (localConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceIndex),
                                                         const_cast <int *> (localPolyhedraCellToFaceConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity_index),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity),
                                                         const_cast <double *> (distantCoords),
                                                         const_cast <int *> (distantLocation),
                                                         const_cast <float *> (distantDistance),
                                                         const_cast <int *> (barycentricCoordinatesIndex),
                                                         const_cast <double *> (barycentricCoordinates),
                                                         const_cast <int *> (&stride),
                                                         const_cast <int *> ((const int *) &_solverType),
                                                         const_cast <double *> (sendingField),
                                                         const_cast <double *> (&tmpDistantField[0]),
                                                         _interpolationFct_f
                                                         );
#else
            fprintf(stderr,"CWIPI Error : Impossible to use a Fortran user interpolation function");
            abort();
#endif
          }
          else
            _interpolate((double* )sendingField,
                         tmpDistantField,
                         stride);
        }
        else {
          const double *dist_uvw = fvmc_locator_get_dist_uvw(_locationToLocalMesh->getFVMLocator());

          if (ptFortranInterpolationFct != NULL) {
#ifndef CWP_HAVE_NOT_FORTRAN_IN_C

            PROCF(callforthointerpfct, CALLFORTHOINTERPFCT) (
                                                         const_cast <int *> (&_entitiesDim),
                                                         const_cast <int *> (&order),
                                                         const_cast <int *> (&nVertex),
                                                         const_cast <int *> (&nElts),
                                                         const_cast <int *> (&nPoly),
                                                         const_cast <int *> (&nDistantPoint),
                                                         const_cast <double *> (_supportMesh->getVertexCoords()),
                                                         const_cast <int *> (localConnectivityIndex),
                                                         const_cast <int *> (localConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceIndex),
                                                         const_cast <int *> (localPolyhedraCellToFaceConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity_index),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity),
                                                         const_cast <double *> (distantCoords),
                                                         const_cast <int *> (distantLocation),
                                                         const_cast <float *> (distantDistance),
                                                         const_cast <int *> (barycentricCoordinatesIndex),
                                                         const_cast <double *> (barycentricCoordinates),
                                                         const_cast <double *> (dist_uvw),
                                                         const_cast <int *> (&stride),
                                                         const_cast <int *> ((const int *) &_solverType),
                                                         const_cast <double *> (sendingField),
                                                         const_cast <double *> (&tmpDistantField[0]),
                                                         ptFortranInterpolationFct
                                                         );
#else
            fprintf(stderr,"CWIPI Error : Impossible to use a Fortran user interpolation function");
            abort();
#endif
          }

          //
          // Callback C

          else if (_interpolationFct != NULL)
            _ho_interpolationFct(_entitiesDim,
                                 order,
                                 nVertex,
                                 nElts,
                                 nPoly,
                                 nDistantPoint,
                                 _supportMesh->getVertexCoords(),
                                 localConnectivityIndex,
                                 localConnectivity,
                                 localPolyhedraFaceIndex,
                                 localPolyhedraCellToFaceConnectivity,
                                 localPolyhedraFaceConnectivity_index,
                                 localPolyhedraFaceConnectivity,
                                 distantCoords,
                                 distantLocation,
                                 distantDistance,
                                 barycentricCoordinatesIndex,
                                 barycentricCoordinates,
                                 dist_uvw,
                                 stride,
                                 _solverType,
                                 sendingField,
                                 &tmpDistantField[0]);

          //
          // Callback Fortran appele en C

          else if (_interpolationFct_f != NULL) {
#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
            PROCF(callforthointerpfct, CALLFORTHOINTERPFCT) (
                                                         const_cast <int *> (&_entitiesDim),
                                                         const_cast <int *> (&order),
                                                         const_cast <int *> (&nVertex),
                                                         const_cast <int *> (&nElts),
                                                         const_cast <int *> (&nPoly),
                                                         const_cast <int *> (&nDistantPoint),
                                                         const_cast <double *> (_supportMesh->getVertexCoords()),
                                                         const_cast <int *> (localConnectivityIndex),
                                                         const_cast <int *> (localConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceIndex),
                                                         const_cast <int *> (localPolyhedraCellToFaceConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity_index),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity),
                                                         const_cast <double *> (distantCoords),
                                                         const_cast <int *> (distantLocation),
                                                         const_cast <float *> (distantDistance),
                                                         const_cast <int *> (barycentricCoordinatesIndex),
                                                         const_cast <double *> (barycentricCoordinates),
                                                         const_cast <double *> (dist_uvw),
                                                         const_cast <int *> (&stride),
                                                         const_cast <int *> ((const int *) &_solverType),
                                                         const_cast <double *> (sendingField),
                                                         const_cast <double *> (&tmpDistantField[0]),
                                                         _ho_interpolationFct_f
                                                         );
#else
            fprintf(stderr,"CWIPI Error : Impossible to use a Fortran user interpolation function");
            abort();
#endif
          }
          else
            _interpolate((double* )sendingField,
                         tmpDistantField,
                         stride);


        }

      }

      //
      // Exchange

      double *ptSending = NULL;

      if (sendingField != NULL)
        ptSending = &tmpDistantField[0];

#ifdef NAN
      if (receivingField != NULL && nInteriorList > 0){
        const int idx = 0;
        receivingField[idx] = NAN;
      }
#endif

      MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
      if (oldFVMComm != MPI_COMM_NULL)
        MPI_Barrier(oldFVMComm);
      fvmc_parall_set_mpi_comm(_fvmComm);

      fvmc_locator_exchange_point_var(_locationToLocalMesh->getFVMLocator(),
                                      (void *) ptSending,
                                      (void *) receivingField,
                                      NULL,
                                      sizeof(double),
                                      stride,
                                      0);
      if (_fvmComm != MPI_COMM_NULL)
        MPI_Barrier(_fvmComm);
      fvmc_parall_set_mpi_comm(oldFVMComm);

      //
      // Check receiving

#ifdef NAN
      if (receivingField != NULL && nInteriorList > 0) {
        const int idx = 0;
        if (std::isnan(receivingField[idx]))
          status = CWIPI_EXCHANGE_BAD_RECEIVING;
      }
#endif

      //
      // Visualization

      if (status == CWIPI_EXCHANGE_OK)

        _fieldsVisualization(exchangeName,
                             stride,
                             timeStep,
                             timeValue,
                             sendingFieldName,
                             sendingField,
                             receivingFieldName,
                             receivingField);

      if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {


        MPI_Comm_rank(localComm, &rootRank);

        assert(rootRank == 0);

        MPI_Bcast(&status,
                  1,
                  MPI_INT,
                  0,
                  localComm);

        if( receivingField != NULL )

          MPI_Bcast(receivingField,
                    stride* _locationToDistantMesh->getNLocatedPoint(),
                    MPI_DOUBLE,
                    0,
                    localComm);

      }
    }

    else if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      MPI_Bcast(&status,
                1,
                MPI_INT,
                0,
                localComm);

      if (receivingField != NULL)

        MPI_Bcast(receivingField,
                  stride* _locationToDistantMesh->getNLocatedPoint(),
                  MPI_DOUBLE,
                  0,
                  localComm );
    }

    return status;
  }


  void oldCoupling::issend(const char    *exchangeName,
                        const int     tag,
                        const int     stride,
                        const int     timeStep,
                        const double  timeValue,
                        const char    *sendingFieldName,
                        const double  *sendingField,
                        void          *ptFortranInterpolationFct,
                        int           *request)


  {

    if (!_isCoupledRank && sendingField != NULL )
      bftc_printf("Warning : sendingField != NULL, "
                  " only field defined by the root rank is sent\n");

    //
    // Locate

    if (_toLocate)
      bftc_error(__FILE__, __LINE__, 0,
                 "Call cwipi_locate before cwipi_isend\n");

    //
    // Prepare data (interpolate, extrapolate...)

    if (_isCoupledRank) {

      const int nVertex                 = _supportMesh->getNVertex();
      const int nElts                   = _supportMesh->getNElts();
      const int nPoly                   = _supportMesh->getNPolyhedra();
      const int *localConnectivityIndex = _supportMesh->getEltConnectivityIndex();
      const int *localConnectivity      = _supportMesh->getEltConnectivity();
      const int *localPolyhedraFaceIndex = _supportMesh->getPolyhedraFaceIndex() ;
      const int *localPolyhedraCellToFaceConnectivity = _supportMesh->getPolyhedraCellToFaceConnectivity();
      const int *localPolyhedraFaceConnectivity_index = _supportMesh->getPolyhedraFaceConnectivityIndex();
      const int *localPolyhedraFaceConnectivity       = _supportMesh->getPolyhedraFaceConnectivity();

      const int nDistantPoint      =  _locationToLocalMesh->getNLocatedDistantPoint() ;
      const int *distantLocation   = _locationToLocalMesh->getLocation();
      const float *distantDistance   = _locationToLocalMesh->getDistance();
      const double *distantCoords   = _locationToLocalMesh->getPointCoordinates();

      const double *barycentricCoordinates = _locationToLocalMesh->getBarycentricCoordinates();
      const int *barycentricCoordinatesIndex = _locationToLocalMesh->getBarycentricCoordinatesIndex();

      int lDistantField = stride * nDistantPoint;

      // plutot mettre request comme clé du map

      std::vector<double>& tmpDistantField = *new std::vector<double> (lDistantField, 0);

      //
      // Interpolation

      if (sendingField != NULL) {

        const int order = _supportMesh->getOrder();

        if (order == -1) {
          assert(!(_interpolationFct != NULL && ptFortranInterpolationFct != NULL));

          //
          // Callback Fortran

          if (ptFortranInterpolationFct != NULL) {
#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
            PROCF(callfortinterpfct, CALLFORTINTERPFCT) (
                                                         const_cast <int *> (&_entitiesDim),
                                                         const_cast <int *> (&nVertex),
                                                         const_cast <int *> (&nElts),
                                                         const_cast <int *> (&nPoly),
                                                         const_cast <int *> (&nDistantPoint),
                                                         const_cast <double *> (_supportMesh->getVertexCoords()),
                                                         const_cast <int *> (localConnectivityIndex),
                                                         const_cast <int *> (localConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceIndex),
                                                         const_cast <int *> (localPolyhedraCellToFaceConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity_index),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity),
                                                         const_cast <double *> (distantCoords),
                                                         const_cast <int *> (distantLocation),
                                                         const_cast <float *> (distantDistance),
                                                         const_cast <int *> (barycentricCoordinatesIndex),
                                                         const_cast <double *> (barycentricCoordinates),
                                                         const_cast <int *> (&stride),
                                                         const_cast <int *> ((const int *) &_solverType),
                                                         const_cast <double *> (sendingField),
                                                         const_cast <double *> (&tmpDistantField[0]),
                                                         ptFortranInterpolationFct
                                                         );

#else
            fprintf(stderr,"CWIPI Error : Impossible to use a Fortran user interpolation function");
            abort();
#endif
        }
        //
        // Callback C

          else if (_interpolationFct != NULL)
            _interpolationFct(_entitiesDim,
                              nVertex,
                              nElts,
                              nPoly,
                              nDistantPoint,
                              _supportMesh->getVertexCoords(),
                              localConnectivityIndex,
                              localConnectivity,
                              localPolyhedraFaceIndex,
                              localPolyhedraCellToFaceConnectivity,
                              localPolyhedraFaceConnectivity_index,
                              localPolyhedraFaceConnectivity,
                              distantCoords,
                              distantLocation,
                              distantDistance,
                              barycentricCoordinatesIndex,
                              barycentricCoordinates,
                              stride,
                              _solverType,
                              sendingField,
                              &tmpDistantField[0]);

          //
          // Callback Fortran appele en C

          else if (_interpolationFct_f != NULL) {
#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
            PROCF(callfortinterpfct, CALLFORTINTERPFCT) (
                                                         const_cast <int *> (&_entitiesDim),
                                                         const_cast <int *> (&nVertex),
                                                         const_cast <int *> (&nElts),
                                                         const_cast <int *> (&nPoly),
                                                         const_cast <int *> (&nDistantPoint),
                                                         const_cast <double *> (_supportMesh->getVertexCoords()),
                                                         const_cast <int *> (localConnectivityIndex),
                                                         const_cast <int *> (localConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceIndex),
                                                         const_cast <int *> (localPolyhedraCellToFaceConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity_index),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity),
                                                         const_cast <double *> (distantCoords),
                                                         const_cast <int *> (distantLocation),
                                                         const_cast <float *> (distantDistance),
                                                         const_cast <int *> (barycentricCoordinatesIndex),
                                                         const_cast <double *> (barycentricCoordinates),
                                                         const_cast <int *> (&stride),
                                                         const_cast <int *> ((const int *) &_solverType),
                                                         const_cast <double *> (sendingField),
                                                         const_cast <double *> (&tmpDistantField[0]),
                                                         _interpolationFct_f
                                                         );

#else
            fprintf(stderr,"CWIPI Error : Impossible to use a Fortran user interpolation function");
            abort();
#endif
          }
          else
            _interpolate((double* )sendingField,
                         tmpDistantField,
                         stride);
        }
        else {
          const double *dist_uvw = fvmc_locator_get_dist_uvw(_locationToLocalMesh->getFVMLocator());

          if (ptFortranInterpolationFct != NULL) {
#ifndef CWP_HAVE_NOT_FORTRAN_IN_C

            PROCF(callforthointerpfct, CALLFORTHOINTERPFCT) (
                                                         const_cast <int *> (&_entitiesDim),
                                                         const_cast <int *> (&order),
                                                         const_cast <int *> (&nVertex),
                                                         const_cast <int *> (&nElts),
                                                         const_cast <int *> (&nPoly),
                                                         const_cast <int *> (&nDistantPoint),
                                                         const_cast <double *> (_supportMesh->getVertexCoords()),
                                                         const_cast <int *> (localConnectivityIndex),
                                                         const_cast <int *> (localConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceIndex),
                                                         const_cast <int *> (localPolyhedraCellToFaceConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity_index),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity),
                                                         const_cast <double *> (distantCoords),
                                                         const_cast <int *> (distantLocation),
                                                         const_cast <float *> (distantDistance),
                                                         const_cast <int *> (barycentricCoordinatesIndex),
                                                         const_cast <double *> (barycentricCoordinates),
                                                         const_cast <double *> (dist_uvw),
                                                         const_cast <int *> (&stride),
                                                         const_cast <int *> ((const int *) &_solverType),
                                                         const_cast <double *> (sendingField),
                                                         const_cast <double *> (&tmpDistantField[0]),
                                                         ptFortranInterpolationFct
                                                         );
#else
            fprintf(stderr,"CWIPI Error : Impossible to use a Fortran user interpolation function");
            abort();
#endif
          }

          //
          // Callback C

          else if (_interpolationFct != NULL)
            _ho_interpolationFct(_entitiesDim,
                                 order,
                                 nVertex,
                                 nElts,
                                 nPoly,
                                 nDistantPoint,
                                 _supportMesh->getVertexCoords(),
                                 localConnectivityIndex,
                                 localConnectivity,
                                 localPolyhedraFaceIndex,
                                 localPolyhedraCellToFaceConnectivity,
                                 localPolyhedraFaceConnectivity_index,
                                 localPolyhedraFaceConnectivity,
                                 distantCoords,
                                 distantLocation,
                                 distantDistance,
                                 barycentricCoordinatesIndex,
                                 barycentricCoordinates,
                                 dist_uvw,
                                 stride,
                                 _solverType,
                                 sendingField,
                                 &tmpDistantField[0]);

          //
          // Callback Fortran appele en C

          else if (_interpolationFct_f != NULL) {
#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
            PROCF(callforthointerpfct, CALLFORTHOINTERPFCT) (
                                                         const_cast <int *> (&_entitiesDim),
                                                         const_cast <int *> (&order),
                                                         const_cast <int *> (&nVertex),
                                                         const_cast <int *> (&nElts),
                                                         const_cast <int *> (&nPoly),
                                                         const_cast <int *> (&nDistantPoint),
                                                         const_cast <double *> (_supportMesh->getVertexCoords()),
                                                         const_cast <int *> (localConnectivityIndex),
                                                         const_cast <int *> (localConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceIndex),
                                                         const_cast <int *> (localPolyhedraCellToFaceConnectivity),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity_index),
                                                         const_cast <int *> (localPolyhedraFaceConnectivity),
                                                         const_cast <double *> (distantCoords),
                                                         const_cast <int *> (distantLocation),
                                                         const_cast <float *> (distantDistance),
                                                         const_cast <int *> (barycentricCoordinatesIndex),
                                                         const_cast <double *> (barycentricCoordinates),
                                                         const_cast <double *> (dist_uvw),
                                                         const_cast <int *> (&stride),
                                                         const_cast <int *> ((const int *) &_solverType),
                                                         const_cast <double *> (sendingField),
                                                         const_cast <double *> (&tmpDistantField[0]),
                                                         _ho_interpolationFct_f
                                                         );
#else
            fprintf(stderr,"CWIPI Error : Impossible to use a Fortran user interpolation function");
            abort();
#endif
          }
          else
            _interpolate((double* )sendingField,
                         tmpDistantField,
                         stride);


        }


      }

      //
      // Send

      double *ptSending = NULL;

      if (sendingField != NULL)
        ptSending = &tmpDistantField[0];


      MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
      if (oldFVMComm != MPI_COMM_NULL)
        MPI_Barrier(oldFVMComm);
      fvmc_parall_set_mpi_comm(_fvmComm);

      fvmc_locator_issend_point_var(_locationToLocalMesh->getFVMLocator(),
                                    (void *) ptSending,
                                    NULL,
                                    sizeof(double),
                                    stride,
                                    0,
                                    tag,
                                    request);

      _tmpDistantFieldsIssend[*request] = &tmpDistantField;

      if (_fvmComm != MPI_COMM_NULL)
        MPI_Barrier(_fvmComm);
      fvmc_parall_set_mpi_comm(oldFVMComm);

      //
      // Visualization

      _fieldsVisualization(exchangeName,
                           stride,
                           timeStep,
                           timeValue,
                           sendingFieldName,
                           sendingField,
                           "",
                           NULL);

    }
  }


  void oldCoupling::waitIssend(int request)
  {
    if (_isCoupledRank) {
      fvmc_locator_issend_wait(_locationToLocalMesh->getFVMLocator(), request);

      delete _tmpDistantFieldsIssend[request];
      _tmpDistantFieldsIssend[request] = NULL;
    }
  }


  void oldCoupling::irecv(const char    *exchangeName,
                       const int     tag,
                       const int     stride,
                       const int     timeStep,
                       const double  timeValue,
                       const char    *receivingFieldName,
                       const double  *receivingField,
                       int           *request)
  {

    if (_isCoupledRank) {

      //
      // irecv

      fvmc_locator_irecv_point_var(_locationToLocalMesh->getFVMLocator(),
                                   (void*) receivingField,
                                   NULL,
                                   sizeof(double),
                                   stride,
                                   0,
                                   tag,
                                   request);


    }

    if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

      MPI_Bcast(request,
                1,
                MPI_INT,
                0,
                localComm);

    }

    _tmpLocalFieldsIrecv[*request] = receivingField;
    _tmpExchangeNameIrecv[*request] = exchangeName;
    _tmpStrideIrecv[*request] = stride;
    _tmpTimeStepIrecv[*request] = timeStep;
    _tmpTimeValueIrecv[*request] = timeValue;
    _tmpFieldNameIrecv[*request] = receivingFieldName;

  }


  void oldCoupling::waitIrecv(int request)
  {

    if (_isCoupledRank) {
      fvmc_locator_irecv_wait(_locationToLocalMesh->getFVMLocator(), request);

      //
      // Visualization

      _fieldsVisualization(_tmpExchangeNameIrecv[request].c_str(),
                           _tmpStrideIrecv[request],
                           _tmpTimeStepIrecv[request],
                           _tmpTimeValueIrecv[request],
                           "",
                           NULL,
                           _tmpFieldNameIrecv[request].c_str(),
                           _tmpLocalFieldsIrecv[request]);
    }

    if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      if(_tmpStrideIrecv[request] != 0 ) {

        const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

        MPI_Bcast((void *)_tmpLocalFieldsIrecv[request],
                  _tmpStrideIrecv[request] * _locationToDistantMesh->getNLocatedPoint(),
                  MPI_DOUBLE,
                  0,
                  localComm);
      }
    }

    _tmpLocalFieldsIrecv[request] = NULL;
    _tmpExchangeNameIrecv[request] = "";
    _tmpStrideIrecv[request] = 0;
    _tmpTimeStepIrecv[request] = 1;
    _tmpTimeValueIrecv[request] = 0;
    _tmpFieldNameIrecv[request] = "";
  }


  void oldCoupling::_initVisualization()
  {
    if (_fvmWriter == NULL && _outputFrequency > 0) {

      bftc_file_mkdir_default("cwipi");

      std::string pathString = "cwipi/"+_name + "_" +
        _localApplicationProperties.getName()+"_"+
        _coupledApplicationProperties.getName();

      MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
      if (oldFVMComm != MPI_COMM_NULL)
        MPI_Barrier(oldFVMComm);
      fvmc_parall_set_mpi_comm(_fvmComm);

      _fvmWriter = fvmc_writer_init("Chr",
                                    pathString.c_str(),
                                    _outputFormat.c_str(),
                                    _outputFormatOption.c_str(),
                                    FVMC_WRITER_FIXED_MESH);

      int localCommSize;
      MPI_Comm_size(_fvmComm, &localCommSize);
      int localRank = 0;
      MPI_Comm_rank(_fvmComm, &localRank);

      fvmc_writer_export_nodal(_fvmWriter, &(_supportMesh->getFvmNodal()));

      // Export sub domain

      if (localCommSize > 1) {

        const int nElts  = _supportMesh->getNElts();

        int *domLoc =  (int *) malloc (sizeof(int) * (nElts));

        for (int i = 0; i < nElts; i++)
          domLoc[i] = localRank+1;

        fvmc_writer_export_field(_fvmWriter,
                                 const_cast<fvmc_nodal_t *> (&_supportMesh->getFvmNodal()),
                                 "partitioning",
                                 FVMC_WRITER_PER_ELEMENT,
                                 1,
                                 FVMC_NO_INTERLACE,
                                 0,
                                 NULL,
                                 FVMC_INT32,
                                 -1,
                                 0.0,
                                 (const void *const *)  &domLoc);
        free ( domLoc);
      }

      // TODO: A deplacer et a recreer en cas de maillage mobile

      int nNotLocatedPointSum;
      int nUnlocatedPoint = _locationToDistantMesh->getNUnlocatedPoint();
      MPI_Allreduce (&nUnlocatedPoint, &nNotLocatedPointSum,
                     1, MPI_INT, MPI_SUM,
                     _fvmComm);

      if (nNotLocatedPointSum != 0 && _locationToDistantMesh->getCoordsPointsToLocate() == NULL) {
        if (_solverType == CWIPI_SOLVER_CELL_CENTER) {
          const int nElts  = _supportMesh->getNElts();

          int *domLoc =  (int *) malloc (sizeof(int) * (nElts));

          for (int i = 0; i < nElts; i++)
            domLoc[i] = 1;

          for (int i = 0; i < _locationToDistantMesh->getNUnlocatedPoint(); i++)
            domLoc[_locationToDistantMesh->getUnlocatedPoint()[i]-1] = 0;

          fvmc_writer_export_field(_fvmWriter,
                                   const_cast<fvmc_nodal_t *> (&_supportMesh->getFvmNodal()),
                                   "location",
                                   FVMC_WRITER_PER_ELEMENT,
                                   1,
                                   FVMC_NO_INTERLACE,
                                   0,
                                   NULL,
                                   FVMC_INT32,
                                   -1,
                                   0.0,
                                   (const void *const *)  &domLoc);
        }
        else {
          const int nVertex = _supportMesh->getNVertex();

          int *domLoc =  (int *) malloc (sizeof(int) * (nVertex));

          for (int i = 0; i < nVertex; i++)
            domLoc[i] = 1;

          for (int i = 0; i < _locationToDistantMesh->getNUnlocatedPoint(); i++)
            domLoc[_locationToDistantMesh->getUnlocatedPoint()[i]-1] = 0;

          fvmc_writer_export_field(_fvmWriter,
                                   const_cast<fvmc_nodal_t *> (&_supportMesh->getFvmNodal()),
                                   "location",
                                   FVMC_WRITER_PER_NODE,
                                   1,
                                   FVMC_NO_INTERLACE,
                                   0,
                                   NULL,
                                   FVMC_INT32,
                                   -1,
                                   0.0,
                                   (const void *const *)  &domLoc);
          free ( domLoc);
        }
      }


      if (_fvmComm != MPI_COMM_NULL)
        MPI_Barrier(_fvmComm);
      fvmc_parall_set_mpi_comm(oldFVMComm);

    }
  }


  void oldCoupling::_fieldsVisualization(const char *exchangeName,
                                      const int stride,
                                      const int timeStep,
                                      const double timeValue,
                                      const char  *sendingFieldName,
                                      const void *sendingField,
                                      const char  *receivingFieldName,
                                      const void *receivingField)
  {

    const double couplingDoubleMax = 0.0;
    //const double couplingDoubleMax = NAN;

    std::string localName;

    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();

    if ((_outputFrequency > 0) && (timeStep % _outputFrequency == 0)) {

      assert(_fvmWriter != NULL);

      if (oldFVMComm != MPI_COMM_NULL)
        MPI_Barrier(oldFVMComm);
      fvmc_parall_set_mpi_comm(_fvmComm);

      fvmc_writer_var_loc_t fvmc_writer_var_loc;
      fvmc_interlace_t fvmc_interlace;
      int dim;
      int fieldSize;

      if (_solverType == CWIPI_SOLVER_CELL_CENTER) {
        fieldSize = _supportMesh->getNElts() + _supportMesh->getNPolyhedra();
        fvmc_writer_var_loc = FVMC_WRITER_PER_ELEMENT;
      }
      else{
        fieldSize = _supportMesh->getNVertex();
        fvmc_writer_var_loc = FVMC_WRITER_PER_NODE;
      }
      //
      // Visualization if stride is scalar or vector

      //if (stride == 1 || stride == 3) {

      //     dim = stride;
      //     if (stride == 1)
      //       fvmc_interlace = FVMC_NO_INTERLACE;
      //     else
      //       fvmc_interlace = FVMC_INTERLACE;

      dim = 1;
      fvmc_interlace = FVMC_NO_INTERLACE;

      //TODO : revoir les types correspondance void / double

      if (sendingFieldName != NULL) {


        double *cpSendingField = NULL;
        void **ptField = NULL;

        for (int kk = 0; kk < stride; kk++) {
          std::ostringstream num;
          num << kk;

          localName = "S_c" +num.str()+ "_" + std::string(exchangeName) +
            "_" + std::string(sendingFieldName)+"_c"+num.str();

          if (sendingField != NULL) {

            ptField = const_cast <void **> (&sendingField);
            if (stride > 1) {
              if (cpSendingField == NULL)
                cpSendingField =  (double *) malloc (sizeof(double) * (fieldSize));
              for (int jj =0; jj < fieldSize; jj++) {
                cpSendingField[jj] = ((double *) sendingField)[jj*stride+kk];
              }
              ptField = (void **) &cpSendingField;
            }


            fvmc_writer_export_field(_fvmWriter,
                                     const_cast<fvmc_nodal_t *> (&_supportMesh->getFvmNodal()),
                                     localName.c_str(),
                                     fvmc_writer_var_loc,
                                     dim,
                                     fvmc_interlace,
                                     0,
                                     NULL,
                                     FVMC_DOUBLE,
                                     timeStep,
                                     timeValue,
                                     (const void *const *) ptField);
          }
        }

        if (stride > 1 && cpSendingField != NULL)
          free ( cpSendingField);

      }

      if (receivingFieldName != NULL) {

        //
        // Not located point treatment

        if (receivingField != NULL && _locationToDistantMesh->getCoordsPointsToLocate() == NULL) {

          int lReceivingField = _locationToDistantMesh->getNpointsToLocate();
          double *cpReceivingField =  (double *) malloc (sizeof(double) * (lReceivingField));
          const int nLocatedPoint = _locationToDistantMesh->getNpointsToLocate() - _locationToDistantMesh->getNUnlocatedPoint();

          for (int j = 0; j < stride; j++) {

            for (int i = 0; i < lReceivingField; i++)
              cpReceivingField[i] = couplingDoubleMax;

            std::ostringstream num;
            num << j;

            localName = "R_c" + num.str() + "_" + std::string(exchangeName) +
              "_" + std::string(receivingFieldName);

            for (int i = 0; i < nLocatedPoint; i++)
              cpReceivingField[(_locationToDistantMesh->getLocatedPoint()[i]-1)] = ((double*) receivingField)[i*stride+j];

            fvmc_writer_export_field(_fvmWriter,
                                     const_cast<fvmc_nodal_t *> (&_supportMesh->getFvmNodal()),
                                     localName.c_str(),
                                     fvmc_writer_var_loc,
                                     dim,
                                     fvmc_interlace,
                                     0,
                                     NULL,
                                     FVMC_DOUBLE,
                                     timeStep,
                                     timeValue,
                                     (const void *const *) &cpReceivingField);
          }
          free ( cpReceivingField);
        }
      }
    }

    if (_fvmComm != MPI_COMM_NULL)
      MPI_Barrier(_fvmComm);
    fvmc_parall_set_mpi_comm(oldFVMComm);
    // } //if stride

  }


  void oldCoupling::locate()
  {
    const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();

    if (_isCoupledRank) {

      if (oldFVMComm != _fvmComm) {
        if (oldFVMComm != MPI_COMM_NULL)
          MPI_Barrier(oldFVMComm);
        fvmc_parall_set_mpi_comm(_fvmComm);
      }
    }

    if (_toLocate) {
      _locationToLocalMesh->setSupportMesh(_supportMesh);
    }

    _locationToLocalMesh->locate();

    /* Ajouter les échanges des distances */

    int nLocPts = getNLocatedPoint();

    if ((int) _distance.size() !=  nLocPts)
      _distance.resize(nLocPts);

    if (_isCoupledRank) {

      float *distantDistance = const_cast <float *> (_locationToLocalMesh->getDistance());

      fvmc_locator_exchange_point_var(_locationToLocalMesh->getFVMLocator(),
                                      (void *) distantDistance,
                                      (void *) &(_distance[0]),
                                      NULL,
                                      sizeof(float),
                                      1,
                                      0);
    }

    if (_couplingType == CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING) {

      MPI_Bcast(&(_distance[0]), nLocPts, MPI_FLOAT, 0, localComm );

    }

    if (_isCoupledRank)
      _initVisualization();

    if (_fvmComm != MPI_COMM_NULL)
      MPI_Barrier(_fvmComm);
    fvmc_parall_set_mpi_comm(oldFVMComm);

    _toLocate = false;

  }


  void oldCoupling::_createCouplingComm()
  {

    const MPI_Comm& localComm = _localApplicationProperties.getLocalComm();

    if( _couplingComm == MPI_COMM_NULL ) {

      const int localBeginningRank = _localApplicationProperties.getBeginningRank();
      const int localEndRank       = _localApplicationProperties.getEndRank();

      const MPI_Comm& globalComm = _localApplicationProperties.getGlobalComm();
      int currentRank;
      MPI_Comm_rank(globalComm, &currentRank);

      const int distantBeginningRank = _coupledApplicationProperties.getBeginningRank();
      const int distantEndRank = _coupledApplicationProperties.getEndRank();

      //
      // Define coupledApplicationComm

      int nLocalRank = localEndRank - localBeginningRank + 1;
      int nDistantRank = distantEndRank - distantBeginningRank + 1;

      assert(localBeginningRank != distantBeginningRank);

      //
      // TODO: Trouver un tag unique pour une paire de codes (Risque de pb dans MPI_Intercomm_create)
      //       Lors de la creation de 2 objets couplages de maniere simultanee

      const int tag = 1;

      MPI_Comm tmpInterComm;

      MPI_Intercomm_create(localComm, 0, globalComm, distantBeginningRank, tag, &tmpInterComm);

      int order;
      if (localBeginningRank < distantBeginningRank) {
        order = 0;
        _coupledApplicationBeginningRankCouplingComm = nLocalRank ;
      }
      else {
        order = 1;
        _coupledApplicationBeginningRankCouplingComm = 0 ;
      }

      MPI_Intercomm_merge(tmpInterComm, order, &_mergeComm);

      MPI_Comm_free(&tmpInterComm);

      int coupledApplicationCommSize;

      MPI_Comm_size(_mergeComm, &coupledApplicationCommSize);

      assert(coupledApplicationCommSize == nLocalRank + nDistantRank);

      //
      // Exchange coupling type

      cwipi_coupling_type_t* couplingTypes =
           (cwipi_coupling_type_t*) malloc (sizeof(cwipi_coupling_type_t) * coupledApplicationCommSize);

      MPI_Allgather((void*)& _couplingType,
                    1,
                    MPI_INT,
                    couplingTypes,
                    1,
                    MPI_INT,
                    _mergeComm);

      //
      // Check coupling type

      int begRank = 0;
      int endRank = 0;

      if (localBeginningRank < distantBeginningRank) {
        begRank = 0;
        endRank = nLocalRank;
      }
      else {
        begRank = nDistantRank;
        endRank = nLocalRank + nDistantRank;
      }

      for (int i = begRank; i < endRank; i ++)
        if (couplingTypes[i] != _couplingType)
          bftc_error(__FILE__, __LINE__, 0,
                     "Two different coupling types for the '%s' application\n",
                     _localApplicationProperties.getName().c_str());


      _coupledApplicationNRankCouplingComm = nDistantRank;

      _isCoupledRank = _localApplicationProperties.getBeginningRank() == currentRank ||
        _couplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING;

      cwipi_coupling_type_t distantCouplingType;

      if (localBeginningRank < distantBeginningRank)
        distantCouplingType = couplingTypes[nLocalRank];
      else
        distantCouplingType = couplingTypes[0];

      free ( couplingTypes);

      //
      // Build coupling communicator

      if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING ||
          distantCouplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {

        _rankList = new int[coupledApplicationCommSize];
        _nRankList = -1;

        if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING &&
            distantCouplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {

          _nRankList = 2;
          _rankList[0] = 0;

          _coupledApplicationNRankCouplingComm = 1;
          if (localBeginningRank < distantBeginningRank) {
            _rankList[1] = nLocalRank;
            _coupledApplicationBeginningRankCouplingComm = 1;
          }
          else {
            _rankList[1] = nDistantRank;
            _coupledApplicationBeginningRankCouplingComm = 0;
          }
        }

        else if (distantCouplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
          _nRankList = 1 + nDistantRank;
          _coupledApplicationNRankCouplingComm = nDistantRank;
          if (localBeginningRank < distantBeginningRank) {
            _rankList[0] = 0;
            _coupledApplicationBeginningRankCouplingComm = 1;
            for (int i = 0; i < nDistantRank; i++)
              _rankList[i+1] = nLocalRank + i;
          }
          else {
            _coupledApplicationBeginningRankCouplingComm = 0;
            for (int i = 0; i < nDistantRank; i++)
              _rankList[i] = i;
            _rankList[nDistantRank] = nDistantRank;
          }
        }

        else if (_couplingType == CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {
          _nRankList = 1 + nLocalRank;
          _coupledApplicationNRankCouplingComm = 1;
          if (localBeginningRank < distantBeginningRank) {
            _coupledApplicationBeginningRankCouplingComm = nLocalRank;
            for (int i = 0; i < nLocalRank; i++)
              _rankList[i] = i;
            _rankList[nLocalRank] = nLocalRank;
          }

          else {
            _rankList[0] = 0;
            _coupledApplicationBeginningRankCouplingComm = 0;
            for (int i = 0; i < nLocalRank; i++)
              _rankList[i+1] = nDistantRank + i;
          }
        }

        else {

          bftc_error(__FILE__, __LINE__, 0,
                     "Error in 'build coupling communicator'\n");
        }

        MPI_Group mergeGroup = MPI_GROUP_NULL;
        MPI_Group couplingGroup = MPI_GROUP_NULL;

        MPI_Comm_group(_mergeComm, &mergeGroup);

        MPI_Group_incl(mergeGroup, _nRankList, _rankList, &couplingGroup);

        MPI_Comm_create(_mergeComm, couplingGroup, &_couplingComm);

        MPI_Group_free(&couplingGroup);
        MPI_Group_free(&mergeGroup);

      }
      else {

        // added for get communicator
        _nRankList = endRank - begRank;
        _rankList = new int[_nRankList];
        for (int i = 0; i < _nRankList; i++) {
          _rankList[i] = begRank + i;
        }
        MPI_Comm_dup(_mergeComm, &_couplingComm);
      }
    }


    if (_fvmComm == MPI_COMM_NULL) {

      if (_couplingType != CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING) {

        int list1 = 0;
        MPI_Group localGroup = MPI_GROUP_NULL;
        MPI_Group fvmGroup = MPI_GROUP_NULL;

        MPI_Comm dupLocalComm = MPI_COMM_NULL;
        MPI_Comm_dup(localComm, &dupLocalComm);

        MPI_Comm_group(dupLocalComm, &localGroup);
        MPI_Group_incl(localGroup, 1, &list1, &fvmGroup);
        MPI_Comm_create(localComm,
                        fvmGroup,
                        &_fvmComm);
        MPI_Group_free(&localGroup);
        MPI_Group_free(&fvmGroup);
        MPI_Comm_free(&dupLocalComm);
      }

      else

        MPI_Comm_dup(localComm, &_fvmComm);

    }
  }

  ///
  /// \brief Exchange field on vertices of cells that contain each located points
  ///

  void oldCoupling::exchangeCellVertexFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride)
  {
    _locationToLocalMesh->exchangeCellVertexFieldOfElementContaining(sendingField, receivingField, stride);
  }

  ///
  /// \brief Exchange field on cells that contain each located points
  ///

  void oldCoupling::exchangeCellCenterFieldOfElementContaining (double *sendingField,  double *receivingField, const int stride)
  {
    _locationToLocalMesh->exchangeCellCenterFieldOfElementContaining(sendingField, receivingField, stride);
  }


} // namespace cwipi
