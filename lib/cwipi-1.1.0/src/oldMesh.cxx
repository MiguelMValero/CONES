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

#include <iostream>

#include <bftc_error.h>
#include <bftc_printf.h>

#include <fvmc_nodal_append.h>
#include <fvmc_nodal_order.h>
#include <fvmc_parall.h>

#include "oldMesh.hxx"
#include "quickSort.h"
#include "vectorUtilities.hxx"
#include "geomUtilities.hxx"

namespace cwipi {

  oldMesh::oldMesh(const MPI_Comm &localComm,
             const int nDim,
             const int nVertex,
             const int nElts,
             double* coords,
             int *eltConnectivityIndex,
             int *eltConnectivity,
             int order
             )
    : _localComm(localComm),
      _nDim(nDim), _order(order), _nVertex(nVertex), _nElts(nElts), _nPolyhedra(0), _coords(coords),
      _polygonIndex(NULL), _eltConnectivityIndex(eltConnectivityIndex), _eltConnectivity(eltConnectivity),
      //_hoOrdering (NULL),
      _polyhedraFaceIndex(NULL), _polyhedraCellToFaceConnectivity(NULL), _polyhedraNFace(0),
      _polyhedraFaceConnectivityIndex(NULL), _polyhedraFaceConnectivity(NULL), 
      _polyhedraCellToVertexConnectivity(NULL),
      _polyhedraCellToVertexConnectivityIndex(NULL), 
      _isNodalFinalized(false), _cellCenterCoords(NULL),
      _cellVolume(NULL),_fvmNodal(NULL), _normalFace(NULL), 
      _characteristicLength(NULL),
      _isDegenerated(NULL)

  {

    

    int n_vtx_tria = 3; 
    int n_vtx_quad = 4; 
    int n_vtx_tetra = 4; 
    int n_vtx_hexa = 8; 
    int n_vtx_prism = 6;
    int n_vtx_pyramid = 5;

    if (order >= 1) { 

      n_vtx_tria = (order+1)*(order+2)/2; 
      n_vtx_quad = (order+1)*(order+1); 
      n_vtx_tetra = (order+1)*(order+2)*(order+3)/6; 
      n_vtx_hexa = (order+1)*(order+1)*(order+1); 
      n_vtx_prism = (order+1)*(order+1)*(order+2)/2; 
      n_vtx_pyramid = (order+1)*(order+2)*(2*order+3)/6;

    }
    
    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvmc_parall_set_mpi_comm(localComm);

    //
    // Check dim

    if (_nDim > 3 || _nDim < 1)
      bftc_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);

    //
    // Check order

    bool sorted = true;

    int nbTriangle   = 0;
    int nbQuadrangle = 0;
    int nbPoly       = 0;

    int nbTetra    = 0;
    int nbPyramid  = 0;
    int nbPrism    = 0;
    int nbHexaedra = 0;

    if (_nDim > 1) {

      if (_nDim == 2) {
        for (int i = 0; i < _nElts; i++) {
          int nCurrentEltVertex = eltConnectivityIndex[i+1] - eltConnectivityIndex[i];
          if (nCurrentEltVertex == n_vtx_tria) {
            if (nbQuadrangle != 0 ||
                nbPoly       != 0)
              sorted = false;
            ++nbTriangle;
          }

          else if (nCurrentEltVertex == n_vtx_quad) {
            if (nbPoly != 0)
              sorted = false;
            ++nbQuadrangle;
          }

          else if (nCurrentEltVertex > n_vtx_quad) {
            if (order == -1) {
              ++nbPoly;
            }
            else {
              bftc_error(__FILE__, __LINE__, 0, "order > 1 for a polygon is not available\n");
            }
          }

          else
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");

        }
      }

      else if (_nDim == 3) {
        
        for (int i = 0; i < _nElts; i++) {
          int nCurrentEltVertex = eltConnectivityIndex[i+1] - eltConnectivityIndex[i];
          if (nCurrentEltVertex == n_vtx_tetra) {
            if (nbPyramid  != 0  ||
                nbPrism    != 0  ||
                nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbTetra;
          }

          else if (nCurrentEltVertex == n_vtx_pyramid) {
            if (nbPrism    != 0  ||
                nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbPyramid;
          }

          else if (nCurrentEltVertex == n_vtx_prism) {
            if (nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbPrism;
          }

          else if (nCurrentEltVertex == n_vtx_hexa) {
            if (nbPoly     != 0)
              sorted = false;
            ++nbHexaedra;
          }

//          else if (nCurrentEltVertex > 8) {
//            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
//            ++nbPoly;
//          }

          else
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }
      }
      else {
        bftc_error(__FILE__, __LINE__, 0, "order > 1 : not implemented yet\n");          
      }
    }

    //
    // Sorting

    if (!sorted) {

      switch (_nDim) {

      case 1 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Bug for edges\n");
        break;

      case 2 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : triangle, quadrangle\n");
        break;

      case 3 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : tetraedra, pyramid, prism, hexaedra\n");
        break;

      default :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "unknown dimension : %i\n", _nDim);
        break;
      }
    }

    //
    // fvmc_nodal building

    _fvmNodal = fvmc_nodal_create("Mesh", 3, order);

    //
    // Sections building

    switch (_nDim) {

    case 1 :

      fvmc_nodal_append_shared(_fvmNodal,
                              _nElts,
                              FVMC_EDGE,
                              NULL,
                              NULL,
                              NULL,
                              _eltConnectivity,
                              NULL);
      
      break;

    case 2 :

      int nTriangleSum;
      int nQuadrangleSum;
      int nPolySum;

      MPI_Allreduce (&nbTriangle, &nTriangleSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbQuadrangle, &nQuadrangleSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPoly, &nPolySum,
                     1, MPI_INT, MPI_SUM,
                     localComm);
      
      if (nbTriangle != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbTriangle,
                                FVMC_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTriangleSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbQuadrangle != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbQuadrangle,
                                FVMC_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 3*nbTriangle,
                                NULL);


      else if (nQuadrangleSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPoly != 0) {
        _polygonIndex =  (int *) malloc (sizeof(int) * (nbPoly+1));
        for(int i = 0; i < nbPoly+1; i++) {
          _polygonIndex[i] = _eltConnectivityIndex[nbTriangle+nbQuadrangle+i]-_eltConnectivityIndex[nbTriangle+nbQuadrangle];
        }

        fvmc_nodal_append_shared(_fvmNodal,
                                nbPoly,
                                FVMC_FACE_POLY,
                                NULL,
                                NULL,
                                _polygonIndex,
                                _eltConnectivity + 3*nbTriangle + 4*nbQuadrangle,
                                NULL);
      }

      else if (nPolySum != 0) {

        //bftc_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for a empty polygon section\n");

        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_POLY,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);
      }
      
    break;

    case 3 :
      int nTetraSum;
      int nPyramidSum;
      int nPrismSum;
      int nHexaedraSum;

      MPI_Allreduce (&nbTetra, &nTetraSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPyramid, &nPyramidSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPrism, &nPrismSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbHexaedra, &nHexaedraSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      if (nbTetra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbTetra,
                                FVMC_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTetraSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPyramid != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbPyramid,
                                FVMC_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra,
                                NULL);

      else if (nPyramidSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPrism != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbPrism,
                                FVMC_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra + 5*nbPyramid,
                                NULL);

      else if (nPrismSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbHexaedra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbHexaedra,
                                FVMC_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity +  4*nbTetra + 5*nbPyramid + 6*nbPrism,
                                NULL);

      else if (nHexaedraSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);


      break;
    }

     MPI_Barrier(localComm);
     fvmc_parall_set_mpi_comm(oldFVMComm);

  }

  void oldMesh::_finalizeNodal()
  {

    if (!_isNodalFinalized) {

      MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
      if (oldFVMComm != MPI_COMM_NULL)
        MPI_Barrier(oldFVMComm);
      fvmc_parall_set_mpi_comm(_localComm);

      _isNodalFinalized = true;

      //
      // Shared vertices
      
      fvmc_nodal_set_shared_vertices(_fvmNodal, _coords);
      
      //
      // Order Fvm_nodal
      
      int localCommSize = 0;
      int localRank = 0;
      
      MPI_Comm_size(_localComm, &localCommSize);
      MPI_Comm_rank(_localComm, &localRank);
      
      int *allNElts     =  (int *) malloc (sizeof(int) * (localCommSize));
      unsigned int *globalEltNum =  (unsigned int *) malloc (sizeof(unsigned int) * (_nElts));
      
      MPI_Allgather((void *) const_cast<int*> (&_nElts),
                    1,
                    MPI_INT,
                    allNElts,
                    1,
                    MPI_INT,
                    _localComm);
      
      int nGlobal = 0;
      for(int i = 0; i < localRank; i++)
        nGlobal += allNElts[i];
      
      for(int i = 0; i < _nElts; i++)
      globalEltNum[i] = nGlobal + i + 1;
      
      switch (_nDim) {
      case 2 :
        fvmc_nodal_order_faces(_fvmNodal, globalEltNum);
        break;
      case 3 :
        fvmc_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
      }

      fvmc_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

      free ( globalEltNum);

      //
      // global vertex num
      
      unsigned int *globalVertexNum =  (unsigned int *) malloc (sizeof(unsigned int) * (_nVertex));
      
      MPI_Allgather((void *) const_cast<int*> (&_nVertex),
                    1,
                    MPI_INT,
                    allNElts,
                    1,
                    MPI_INT,
                    _localComm);
      
      nGlobal = 0;
      for(int i = 0; i < localRank; i++)
        nGlobal += allNElts[i];
      
      for(int i = 0; i < _nVertex; i++)
        globalVertexNum[i] = nGlobal + i + 1;
      
      
      fvmc_nodal_order_vertices(_fvmNodal, globalVertexNum);
      fvmc_nodal_init_io_num(_fvmNodal, globalVertexNum, 0);
      
      free ( globalVertexNum);
      free ( allNElts);
      
#if defined(DEBUG) && 0
      fvmc_nodal_dump(_fvmNodal);
#endif
      
      MPI_Barrier(_localComm);
      fvmc_parall_set_mpi_comm(oldFVMComm);
    }
  }


  /////////
  /////////////
  ///////////
  /////////

  oldMesh::oldMesh(const MPI_Comm &localComm,
             fvmc_nodal_t* fvmc_nodal)
    : _localComm(localComm),
      _nDim(fvmc_nodal_get_dim(fvmc_nodal)), _nVertex(0),
      _nElts(0), _nPolyhedra(0), _coords(NULL),
      _polygonIndex(NULL),    
      _eltConnectivityIndex(NULL), _eltConnectivity(NULL),
      _polyhedraFaceIndex(NULL), _polyhedraCellToFaceConnectivity(NULL),
      _polyhedraNFace(0),
      _polyhedraFaceConnectivityIndex(NULL), _polyhedraFaceConnectivity(NULL), 
      _polyhedraCellToVertexConnectivity(NULL), 
      _polyhedraCellToVertexConnectivityIndex(NULL),
      _isNodalFinalized(false), _cellCenterCoords(NULL),
      _cellVolume(NULL), _fvmNodal(NULL), _normalFace(NULL),
      _characteristicLength(NULL),
      _isDegenerated(NULL)
  {

    int order = -1;
    
    //
    // Copy
    
    // TODO: Attention dans ce cas on alloue les vecteurs en interne sans les liberer
    //       Restructurer en creant des const * pour les valeurs partagees et détruite les valeurs 
    //       creees dans le destructeur



    fvmc_nodal_get_vertex(fvmc_nodal,
                         &_nElts,
                         &_eltConnectivityIndex,
                         &_eltConnectivity);


    fvmc_nodal_get_coords(fvmc_nodal,
                         &_nVertex,
                         &_coords);


    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvmc_nodal_destroy(fvmc_nodal); 

    fvmc_parall_set_mpi_comm(localComm);

    //
    // Check dim

    if (_nDim > 3 || _nDim < 1)
      bftc_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);

    //
    // Check order

    bool sorted = true;

    int nbTriangle   = 0;
    int nbQuadrangle = 0;
    int nbPoly       = 0;

    int nbTetra    = 0;
    int nbPyramid  = 0;
    int nbPrism    = 0;
    int nbHexaedra = 0;

    if (_nDim > 1) {

      if (_nDim == 2) {


        for (int i = 0; i < _nElts; i++) {
          int nCurrentEltVertex = _eltConnectivityIndex[i+1] - _eltConnectivityIndex[i];
          if (nCurrentEltVertex == 3) {
            if (nbQuadrangle != 0 ||
                nbPoly       != 0)
              sorted = false;
            ++nbTriangle;
          }

          else if (nCurrentEltVertex == 4) {
            if (nbPoly != 0)
              sorted = false;
            ++nbQuadrangle;
          }

          else if (nCurrentEltVertex > 4) {
            ++nbPoly;
          }

          else
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }
      }

      else if (_nDim == 3) {

        for (int i = 0; i < _nElts; i++) {
          int nCurrentEltVertex = _eltConnectivityIndex[i+1] - _eltConnectivityIndex[i];
          if (nCurrentEltVertex == 4) {
            if (nbPyramid  != 0  ||
                nbPrism    != 0  ||
                nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbTetra;
          }

          else if (nCurrentEltVertex == 5) {
            if (nbPrism    != 0  ||
                nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbPyramid;
          }

          else if (nCurrentEltVertex == 6) {
            if (nbHexaedra != 0  ||
                nbPoly     != 0)
              sorted = false;
            ++nbPrism;
          }

          else if (nCurrentEltVertex == 8) {
            if (nbPoly     != 0)
              sorted = false;
            ++nbHexaedra;
          }

          else if (nCurrentEltVertex > 8) {
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
            ++nbPoly;
          }

          else
            bftc_error(__FILE__, __LINE__, 0, "Erreur dans l'index de connectivite\n");
        }

      }
    }

    //
    // Sorting

    if (!sorted) {

      switch (_nDim) {

      case 1 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Bug for edges\n");
        break;

      case 2 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : triangle, quadrangle\n");
        break;

      case 3 :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "Specified order : tetraedra, pyramid, prism, hexaedra\n");
        break;

      default :
        bftc_error(__FILE__, __LINE__, 0, "Connectivity is not ordered\n"
                  "unknown dimension : %i\n", _nDim);
        break;
      }
    }

    //
    // fvmc_nodal building

    _fvmNodal = fvmc_nodal_create("Mesh", 3, order);

    //
    // Sections building

    switch (_nDim) {

    case 1 :


      fvmc_nodal_append_shared(_fvmNodal,
                              _nElts,
                              FVMC_EDGE,
                              NULL,
                              NULL,
                              NULL,
                              _eltConnectivity,
                              NULL);
      break;

    case 2 :

      int nTriangleSum;
      int nQuadrangleSum;
      int nPolySum;

      MPI_Allreduce (&nbTriangle, &nTriangleSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbQuadrangle, &nQuadrangleSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPoly, &nPolySum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      if (nbTriangle != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbTriangle,
                                FVMC_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTriangleSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_TRIA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbQuadrangle != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbQuadrangle,
                                FVMC_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 3*nbTriangle,
                                NULL);


      else if (nQuadrangleSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_QUAD,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPoly != 0) {
        _polygonIndex =  (int *) malloc (sizeof(int) * (nbPoly+1));
        for(int i = 0; i < nbPoly+1; i++) {
          _polygonIndex[i] = _eltConnectivityIndex[nbTriangle+nbQuadrangle+i]-_eltConnectivityIndex[nbTriangle+nbQuadrangle];
        }

        fvmc_nodal_append_shared(_fvmNodal,
                                nbPoly,
                                FVMC_FACE_POLY,
                                NULL,
                                NULL,
                                _polygonIndex,
                                _eltConnectivity + 3*nbTriangle + 4*nbQuadrangle,
                                NULL);
      }

      else if (nPolySum != 0) {

        //bftc_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for a empty polygon section\n");

        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_FACE_POLY,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);
      }
    break;

    case 3 :
      int nTetraSum;
      int nPyramidSum;
      int nPrismSum;
      int nHexaedraSum;

      MPI_Allreduce (&nbTetra, &nTetraSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPyramid, &nPyramidSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbPrism, &nPrismSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      MPI_Allreduce (&nbHexaedra, &nHexaedraSum,
                     1, MPI_INT, MPI_SUM,
                     localComm);

      if (nbTetra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbTetra,
                                FVMC_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity,
                                NULL);

      else if (nTetraSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_TETRA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPyramid != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbPyramid,
                                FVMC_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra,
                                NULL);

      else if (nPyramidSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_PYRAM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbPrism != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbPrism,
                                FVMC_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity + 4*nbTetra + 5*nbPyramid,
                                NULL);

      else if (nPrismSum != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_PRISM,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);

      if (nbHexaedra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                nbHexaedra,
                                FVMC_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                _eltConnectivity +  4*nbTetra + 5*nbPyramid + 6*nbPrism,
                                NULL);

      else if (nbHexaedra != 0)
        fvmc_nodal_append_shared(_fvmNodal,
                                0,
                                FVMC_CELL_HEXA,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL);


      break;
    }

    //
    // Shared vertices

    fvmc_nodal_set_shared_vertices(_fvmNodal, _coords);

    //
    // Order Fvm_nodal

    int localCommSize = 0;
    int localRank = 0;

    MPI_Comm_size(_localComm, &localCommSize);
    MPI_Comm_rank(_localComm, &localRank);

    int *allNElts     =  (int *) malloc (sizeof(int) * (localCommSize));
    unsigned int *globalEltNum =  (unsigned int *) malloc (sizeof(unsigned int) * (_nElts));

    MPI_Allgather((void *) const_cast<int*> (&_nElts),
                  1,
                  MPI_INT,
                  allNElts,
                  1,
                  MPI_INT,
                  _localComm);

    int nGlobal = 0;
    for(int i = 0; i < localRank; i++)
      nGlobal += allNElts[i];

    for(int i = 0; i < _nElts; i++)
      globalEltNum[i] = nGlobal + i + 1;

    switch (_nDim) {
      case 2 :
        fvmc_nodal_order_faces(_fvmNodal, globalEltNum);
        break;
      case 3 :
        fvmc_nodal_order_cells(_fvmNodal, globalEltNum);
        break;
    }

    fvmc_nodal_init_io_num(_fvmNodal, globalEltNum, _nDim);

    free ( globalEltNum);

    //
    // global vertex num

    unsigned int *globalVertexNum =  (unsigned int *) malloc (sizeof(unsigned int) * (_nVertex));

    MPI_Allgather((void *) const_cast<int*> (&_nVertex),
                  1,
                  MPI_INT,
                  allNElts,
                  1,
                  MPI_INT,
                  _localComm);

    nGlobal = 0;
    for(int i = 0; i < localRank; i++)
      nGlobal += allNElts[i];

    for(int i = 0; i < _nVertex; i++)
      globalVertexNum[i] = nGlobal + i + 1;


    fvmc_nodal_order_vertices(_fvmNodal, globalVertexNum);
    fvmc_nodal_init_io_num(_fvmNodal, globalVertexNum, 0);

    free ( globalVertexNum);
    free ( allNElts);

    #if defined(DEBUG) && 0
    fvmc_nodal_dump(_fvmNodal);
    #endif

    MPI_Barrier(localComm);
    fvmc_parall_set_mpi_comm(oldFVMComm);

  }

  oldMesh::~oldMesh()
  {
    delete _cellCenterCoords;
    delete _cellVolume;
    delete _normalFace;
    delete _characteristicLength;
    delete _isDegenerated;
    free ( _polygonIndex);
    delete _polyhedraCellToVertexConnectivity;
    delete _polyhedraCellToVertexConnectivityIndex; 

    fvmc_nodal_destroy(_fvmNodal);
  }


  void oldMesh::addPolyhedra(const int nElt,
                          int *faceIndex,
                          int *cellToFaceConnectivity,
                          const int nFace,
                          int *faceConnectivityIndex,
                          int *faceConnectivity)
  {
    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
    if (oldFVMComm != MPI_COMM_NULL)
      MPI_Barrier(oldFVMComm);
    fvmc_parall_set_mpi_comm(_localComm);

    if (_fvmNodal == NULL)
      bftc_error(__FILE__, __LINE__, 0, "No mesh to add element\n");

    _nPolyhedra += nElt;
    _nElts += nElt;

    _polyhedraNFace                  = nFace;
    _polyhedraFaceIndex              = faceIndex;
    _polyhedraCellToFaceConnectivity = cellToFaceConnectivity;
    _polyhedraFaceConnectivityIndex  = faceConnectivityIndex;
    _polyhedraFaceConnectivity       = faceConnectivity;

    if (nElt > 0)

      fvmc_nodal_append_shared(_fvmNodal,
                              nElt, 
                              FVMC_CELL_POLY,
                              faceIndex,
                              cellToFaceConnectivity,
                              faceConnectivityIndex,
                              faceConnectivity,
                              NULL);
    else {
      //bftc_error(__FILE__, __LINE__, 0, "define Mesh : unresolved bug in fvm for an empty polyedron section\n");


      fvmc_nodal_append_shared(_fvmNodal,
                                 0,
                                 FVMC_CELL_POLY,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL);
    }

    if (_cellCenterCoords != NULL || _cellVolume != NULL)
      _computeMeshProperties();

    fvmc_parall_set_mpi_comm(oldFVMComm);
  }

  
  void oldMesh::hoOrderingSet (const cwipi_element_t t_elt,
                               const int n_nodes,
                               const int *ordering)
  {

    fvmc_element_t _t_elt = (fvmc_element_t) 0;
    
    switch(t_elt) {
    case CWIPI_NODE:
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' Element not taking into account by fvm \n");
      break;
    case CWIPI_EDGEHO:  
    case CWIPI_EDGE2:         
      _t_elt = FVMC_EDGE;
      break;
    case CWIPI_FACE_TRIA3:    
    case CWIPI_FACE_TRIAHO:   
      _t_elt = FVMC_FACE_TRIA;
      break;
    case CWIPI_FACE_QUAD4:
    case CWIPI_FACE_QUADHO:
      _t_elt = FVMC_FACE_QUAD;
      break;
    case CWIPI_FACE_POLY:
      _t_elt = FVMC_FACE_POLY;
      break;
    case CWIPI_CELL_TETRA4:
    case CWIPI_CELL_TETRAHO:
      _t_elt = FVMC_CELL_TETRA;
      break;
    case CWIPI_CELL_PYRAM5:
    case CWIPI_CELL_PYRAMHO:
      _t_elt = FVMC_CELL_PYRAM;
      break;
    case CWIPI_CELL_PRISM6:
    case CWIPI_CELL_PRISMHO:
      _t_elt = FVMC_CELL_PRISM;
      break;
    case CWIPI_CELL_HEXA8:
    case CWIPI_CELL_HEXAHO:
      _t_elt = FVMC_CELL_HEXA;
      break;
    case CWIPI_CELL_POLY:
      _t_elt = FVMC_CELL_POLY;
      break;
    default:
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' Element not taking into account by fvm \n");
      
    }
    
    fvmc_nodal_ho_ordering_set (_fvmNodal, _t_elt, n_nodes, ordering);
    
  }
  
  void oldMesh::hoOrderingFromRefEltSet (const cwipi_element_t t_elt,
                                         const int n_nodes,
                                         const double *coords)
  {

    fvmc_element_t _t_elt = (fvmc_element_t) 0;
    
    switch(t_elt) {
    case CWIPI_NODE:
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' Element not taking into account by fvm \n");
      break;
    case CWIPI_EDGEHO:
    case CWIPI_EDGE2:         
      _t_elt = FVMC_EDGE;
      break;
    case CWIPI_FACE_TRIA3:    
    case CWIPI_FACE_TRIAHO:   
      _t_elt = FVMC_FACE_TRIA;
      break;
    case CWIPI_FACE_QUAD4:
    case CWIPI_FACE_QUADHO:
      _t_elt = FVMC_FACE_QUAD;
      break;
    case CWIPI_FACE_POLY:
      _t_elt = FVMC_FACE_POLY;
      break;
    case CWIPI_CELL_TETRA4:
    case CWIPI_CELL_TETRAHO:
      _t_elt = FVMC_CELL_TETRA;
      break;
    case CWIPI_CELL_PYRAM5:
    case CWIPI_CELL_PYRAMHO:
      _t_elt = FVMC_CELL_PYRAM;
      break;
    case CWIPI_CELL_PRISM6:
    case CWIPI_CELL_PRISMHO:
      _t_elt = FVMC_CELL_PRISM;
      break;
    case CWIPI_CELL_HEXA8:
    case CWIPI_CELL_HEXAHO:
      _t_elt = FVMC_CELL_HEXA;
      break;
    case CWIPI_CELL_POLY:
      _t_elt = FVMC_CELL_POLY;
      break;
    default:
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' Element not taking into account by fvm \n");


    }

    fvmc_nodal_ho_ordering_from_ref_elt_set (_fvmNodal, _t_elt, n_nodes, coords);

  }

  
  void oldMesh::update()
  {
    if (_cellCenterCoords != NULL || _cellVolume != NULL)
      _computeMeshProperties();
  }


  void oldMesh::_computeCellCenterCoordsWithVertex(const int i,
                                                const int nCurrentEltVertex,
                                                const int index,
                                                const int *eltConnectivity,
                                                std::vector<double> *cellCenterCoords)
  {
    assert (_cellCenterCoords != NULL);

    std::vector<double> &refCellCenterCoords = *cellCenterCoords;

    refCellCenterCoords[3*i]   = 0.;
    refCellCenterCoords[3*i+1] = 0.;
    refCellCenterCoords[3*i+2] = 0.;

    for (int j = 0; j < nCurrentEltVertex ; j++) {
      int ivertex = eltConnectivity[index+j] - 1;
      refCellCenterCoords[3*i]   += _coords[3*ivertex];
      refCellCenterCoords[3*i+1] += _coords[3*ivertex+1];
      refCellCenterCoords[3*i+2] += _coords[3*ivertex+2];
    }
    refCellCenterCoords[3*i]   /= nCurrentEltVertex;
    refCellCenterCoords[3*i+1] /= nCurrentEltVertex;
    refCellCenterCoords[3*i+2] /= nCurrentEltVertex;
  }

  void oldMesh::_computeMeshProperties1D()
  {
    std::vector<double> &refCellVolume           = *_cellVolume;
    std::vector<double> &refCharacteristicLength = *_characteristicLength;
    std::vector<double> &refCellCenterCoords     = *_cellCenterCoords;
    std::vector<int>   &refIsDegenerated        = *_isDegenerated;

    edgesProperties (_nElts,
                     _eltConnectivity,
                     _nVertex,
                     _coords,
                     &refCellCenterCoords[0],
                     &refCellVolume[0],
                     &refCharacteristicLength[0],
                     &refIsDegenerated[0]);
  }


  void oldMesh::_computeMeshProperties2D()
  {

    std::vector<double> &refCellVolume           = *_cellVolume;
    std::vector<double> &refCharacteristicLength = *_characteristicLength;
    std::vector<int>   &refIsDegenerated        = *_isDegenerated;
    std::vector<double> &refNormalFace           = *_normalFace;
    std::vector<double> &refCellCenterCoords     = *_cellCenterCoords;

    int convergence = 1;

    //TODO : 11/03/12 : Ecriture temporaire en attendant la séparation en blocs suivant le type d'element

    for (int i = 0; i < _nElts ; i++) {

      const int index      = _eltConnectivityIndex[i];
      const int nEltVertex = _eltConnectivityIndex[i+1] - index;
      const int nLocElt    = 1;
      int eltConvergence  = 1;
      
      if (nEltVertex == 3) 
        
        triangleProperties (nLocElt,
                            _eltConnectivity + index,
                            _nVertex,
                            _coords,
                            &refNormalFace[3*i],
                            &refCellCenterCoords[3*i],
                            &refCharacteristicLength[i],
                            &refIsDegenerated[i]);

      else if (nEltVertex == 4) 
        
        eltConvergence =  quadrangleProperties (nLocElt,
                                                _eltConnectivity + index,
                                                _nVertex,
                                                _coords,
                                                &refNormalFace[3*i],
                                                &refCellCenterCoords[3*i],
                                                &refCharacteristicLength[i],
                                                &refIsDegenerated[i]);
      else
        
        eltConvergence = polygonProperties (nLocElt,   
                                            _eltConnectivityIndex + i,
                                            _eltConnectivity,
                                            _nVertex,
                                            _coords,
                                            &refNormalFace[3*i],
                                            &refCellCenterCoords[3*i],
                                            &refCharacteristicLength[i],
                                            &refIsDegenerated[i]);

      refCellVolume[i] = norm(&refNormalFace[3*i]);
      
      if (!eltConvergence)
        convergence = eltConvergence;
          
    }

    if (!convergence) 
      bftc_printf("Warning computeMeshProperties : some faces are not planar !\n");

  }

  void oldMesh::_computeMeshProperties3D()
  {

    // TODO: Refaire completement les calcul des centres cellules et volume (comme dans CEDRE)

    std::vector<double> &refCellVolume           = *_cellVolume;
    std::vector<double> &refCharacteristicLength = *_characteristicLength;
    std::vector<int>   &refIsDegenerated        = *_isDegenerated;
    std::vector<double> &refCellCenterCoords     = *_cellCenterCoords;

    const int nStandardElement = _nElts - _nPolyhedra;

    for (int i = 0; i < nStandardElement; i++) {
      
      const int index      = _eltConnectivityIndex[i];
      const int nEltVertex = _eltConnectivityIndex[i+1] - index;
      const int nLocElt    = 1;
      
      if (nEltVertex == 4) {
        
        tetrahedraProperties (nLocElt,
                              _eltConnectivity + index,
                              _nVertex,
                              _coords,
                              &refCellVolume[i],
                              &refCellCenterCoords[3*i],
                              &refCharacteristicLength[i],
                              &refIsDegenerated[i]);
      }
      
      else if (nEltVertex == 5) {
        
        pyramidProperties (nLocElt,
                           _eltConnectivity + index,
                           _nVertex,
                           _coords,
                           &refCellVolume[i],
                           &refCellCenterCoords[3*i],
                           &refCharacteristicLength[i],
                           &refIsDegenerated[i]);
        
      }
      
      else if (nEltVertex == 6) {
        
        prismProperties (nLocElt,
                         _eltConnectivity + index,
                         _nVertex,
                         _coords,
                         &refCellVolume[i],
                         &refCellCenterCoords[3*i],
                         &refCharacteristicLength[i],
                         &refIsDegenerated[i]);
        
      }
      
      else if (nEltVertex == 8) {
        
        hexahedraProperties (nLocElt,
                             _eltConnectivity + index,
                             _nVertex,
                             _coords,
                             &refCellVolume[i],
                             &refCellCenterCoords[3*i],
                             &refCharacteristicLength[i],
                             &refIsDegenerated[i]);
        
      }
    }
    
    if (_nPolyhedra > 0)
      polyhedraProperties (_nPolyhedra,
                           _polyhedraNFace,
                           _polyhedraFaceConnectivityIndex,
                           _polyhedraFaceConnectivity,
                           _polyhedraFaceIndex,
                           _polyhedraCellToFaceConnectivity,
                           _nVertex,
                           _coords,
                           &refCellVolume[nStandardElement],
                           &refCellCenterCoords[3*nStandardElement],
                           &refCharacteristicLength[nStandardElement],
                           &refIsDegenerated[nStandardElement]);

  }

  void oldMesh::_computeMeshProperties()
  {
    if (_cellCenterCoords == NULL) {
      _cellCenterCoords = new std::vector<double>(3*_nElts);
      _cellVolume = new std::vector<double>(_nElts);
      _characteristicLength = new std::vector<double>(_nElts);
      _isDegenerated = new std::vector<int>(_nElts);
    }
    else {
      _cellCenterCoords->resize(3*_nElts);
      _cellVolume->resize(_nElts);
      _characteristicLength->resize(_nElts);
      _isDegenerated->resize(_nElts);
    }

    if (_nDim == 1)
        _computeMeshProperties1D();
    else if (_nDim == 2) {
      if (_normalFace == NULL)
        _normalFace = new std::vector<double>(3*_nElts);
      else
        _normalFace->resize(3*_nElts);
      _computeMeshProperties2D();
    }
    else if (_nDim == 3)
      _computeMeshProperties3D();

  }
  
  void oldMesh::_computeMeshPolyhedraProperties(){
    
    int nbrVertex = 0;
    int nbrVertexOld = 0;
    int sizeTabLoc = 12;
    int indFaceVertex;
    int indVertex;
    int nbrVertexFace;
    
    std::vector<int> vertexPolyLoc (sizeTabLoc);
    std::vector<int> vertexPolyBool(_nVertex);

    if(_polyhedraCellToVertexConnectivity == NULL)
       _polyhedraCellToVertexConnectivity = new std::vector<int>(0);
    
    if (_polyhedraCellToVertexConnectivityIndex == NULL)
      _polyhedraCellToVertexConnectivityIndex = new std::vector<int>(_nPolyhedra + 1);

    std::vector<int> & refPolyhedraCellToVertexConnectivity = *_polyhedraCellToVertexConnectivity;
    std::vector<int> & refPolyhedraCellToVertexConnectivityIndex = *_polyhedraCellToVertexConnectivityIndex;

    for (int i = 0; i < _nVertex; i++)
      vertexPolyBool[i] = 0;    
    
    /**** Premiere boucle pour connaitre le nombre total de sommets ****/

    for(int iPoly = 0; iPoly < _nPolyhedra ; iPoly++){

       int nFacePolyhedra = _polyhedraFaceIndex[iPoly+1] - _polyhedraFaceIndex[iPoly];

       for(int iFace = 0; iFace < nFacePolyhedra ; iFace++){

         int face = abs(_polyhedraCellToFaceConnectivity[_polyhedraFaceIndex[iPoly] + iFace]);


         nbrVertexFace = _polyhedraFaceConnectivityIndex[face]
                       - _polyhedraFaceConnectivityIndex[face - 1];
         
         indFaceVertex = _polyhedraFaceConnectivityIndex[face - 1];

         for (int iVertex = 0 ; iVertex < nbrVertexFace ; iVertex++){

           indVertex = _polyhedraFaceConnectivity[indFaceVertex + iVertex];

           if(vertexPolyBool[indVertex - 1] == 0){

             if (nbrVertex >= sizeTabLoc){
               sizeTabLoc *=2;
               vertexPolyLoc.resize(sizeTabLoc);
             }

             vertexPolyLoc[nbrVertex - nbrVertexOld] = indVertex;
             nbrVertex++;
             vertexPolyBool[indVertex - 1] = 1;       
           }           
           
         }

       }

       for(int iVertexPoly = 0 ; iVertexPoly < nbrVertex - nbrVertexOld ; iVertexPoly++)
         vertexPolyBool[vertexPolyLoc[iVertexPoly] - 1] = 0;
 
       nbrVertexOld = nbrVertex;
    }

    vertexPolyLoc.clear();

    refPolyhedraCellToVertexConnectivity.resize(nbrVertex);
    refPolyhedraCellToVertexConnectivityIndex.resize(_nPolyhedra+1);
    refPolyhedraCellToVertexConnectivityIndex[0] = 0;

    nbrVertex = 0;
    nbrVertexOld = 0;

    /**** Calcul des tableaux d'index et de connectivite entre poly et sommets ****/

    for(int iPoly = 0; iPoly < _nPolyhedra ; iPoly++){

       int nFacePolyhedra = _polyhedraFaceIndex[iPoly+1] - _polyhedraFaceIndex[iPoly];

       for(int iFace = 0; iFace < nFacePolyhedra ; iFace++){

         int face = abs(_polyhedraCellToFaceConnectivity[_polyhedraFaceIndex[iPoly] + iFace]);

         nbrVertexFace = _polyhedraFaceConnectivityIndex[face] - _polyhedraFaceConnectivityIndex[face - 1];
         
         indFaceVertex = _polyhedraFaceConnectivityIndex[face - 1];

         for (int iVertex = 0 ; iVertex < nbrVertexFace ; iVertex++){

          indVertex = _polyhedraFaceConnectivity[indFaceVertex + iVertex];

          if(vertexPolyBool[indVertex - 1] == 0){
            refPolyhedraCellToVertexConnectivity[nbrVertex] = indVertex;            
            nbrVertex++;
            vertexPolyBool[indVertex - 1] = 1;                   
          }           
          
        }
        
      }
      
      for(int iVertexPoly = nbrVertexOld ; iVertexPoly < nbrVertex ; iVertexPoly++)
        vertexPolyBool[refPolyhedraCellToVertexConnectivity[iVertexPoly] - 1] = 0;
      
      refPolyhedraCellToVertexConnectivityIndex[iPoly+1] = nbrVertex;
      
      nbrVertexOld = nbrVertex;

    } 

    // bftc_printf("vertex polyhedra : \n");
    // for (int iPoly = 0; iPoly < _nPolyhedra ; iPoly++){
    //   for (int j = refPolyhedraCellToVertexConnectivityIndex[iPoly]; j < refPolyhedraCellToVertexConnectivityIndex[iPoly+1]; j++) {
    //     bftc_printf("%i ",  refPolyhedraCellToVertexConnectivity[j]);
    //   }
    //   bftc_printf("\n");
    // }  
  }      

}
