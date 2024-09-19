#ifndef __MESH_H__
#define __MESH_H__
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

#include <algorithm>
#include <vector>
#include <map>

#include <mpi.h>

#include <fvmc_nodal.h>
#include <bftc_error.h>

#include "cwipi.h"

namespace cwipi {

  class oldMesh {

  public:
    oldMesh(const MPI_Comm &localComm,
         const int nDim,
         const int nVertex,
         const int nElts,
         double* coords,
         int *eltConnectivityIndex,
         int *eltConnectivity,
         int order = -1);

    oldMesh(const MPI_Comm &localComm,
         fvmc_nodal_t* fvmc_nodal);

    virtual ~oldMesh();

    void hoOrderingSet (const cwipi_element_t t_elt,
                        const int n_nodes,
                        const int *ordering);
    
    void hoOrderingFromRefEltSet (const cwipi_element_t t_elt,
                                  const int n_nodes,
                                  const double *coords);

    void addPolyhedra(const int nElt,
                      int *faceIndex,
                      int *cellToFaceConnectivity,
                      const int nFace,
                      int *faceConnectivityIndex,
                      int *faceConnectivity);

    inline const int& getNVertex() const;

    inline const int& getOrder() const;

    inline const double* getVertexCoords() const;

    inline fvmc_nodal_t& getFvmNodal() const;

    inline fvmc_nodal_t& getFvmNodal();

    inline const int& getNElts() const;

    inline const int& getNPolyhedra() const;

    inline const int* getEltConnectivityIndex() const;

    inline const int* getEltConnectivity() const;

    inline const std::vector<double>& getVolume();

    inline const std::vector<double>& getCellCenterCoords();

    inline const std::vector<double>& getNormalFace();    

    inline const std::vector<double>& getCharacteristicLength();    

    inline const std::vector<int>& getIsDegenerated();    

    inline const int *getPolyhedraFaceIndex() const;

    inline const int *getPolyhedraCellToFaceConnectivity() const;

    inline const int *getPolyhedraFaceConnectivityIndex() const;

    inline const int *getPolyhedraFaceConnectivity() const;

    inline const std::vector<int>& getPolyhedraCellToVertexConnectivity();

    inline const std::vector<int>& getPolyhedraCellToVertexConnectivityIndex();

    void update();

  private :
    oldMesh();

    oldMesh(const oldMesh&);

    oldMesh& operator=(const oldMesh&);

  protected :

    void _computeCellCenterCoordsWithVertex(const int i,
                                            const int nCurrentEltVertex,
                                            const int index,
                                            const int *eltConnectivity,
                                            std::vector<double> *cellCenterCoords) ;

    void _computeMeshProperties1D() ;

    void _computeMeshProperties2D();

    void _computeMeshProperties3D();

    void _computeMeshProperties();

    void _finalizeNodal();

    void _computeMeshPolyhedraProperties();


  private:
    // TODO: renommer _nDim par entitesDim
    const MPI_Comm & _localComm;
    int            _nDim;
    int            _order;
    int           _nVertex;
    int           _nElts;
    int           _nPolyhedra;
    double       *_coords;
    int          *_polygonIndex;

    int          *_eltConnectivityIndex;
    int          *_eltConnectivity;

    int              *_polyhedraFaceIndex;
    int              *_polyhedraCellToFaceConnectivity;

    int              _polyhedraNFace;
    int              *_polyhedraFaceConnectivityIndex;
    int              *_polyhedraFaceConnectivity;

    std::vector<int> *_polyhedraCellToVertexConnectivity;
    std::vector<int> *_polyhedraCellToVertexConnectivityIndex;

    bool         _isNodalFinalized;

    std::vector<double>  *_cellCenterCoords;
    std::vector<double>  *_cellVolume;
    fvmc_nodal_t *_fvmNodal;
    std::vector<double>  *_normalFace;
    std::vector<double>  *_characteristicLength;
    std::vector<int>     *_isDegenerated;
  };

  const int& oldMesh::getNVertex()  const
  {
    return _nVertex;
  }

  const int& oldMesh::getOrder()  const
  {
    return _order;
  }
  
  const double* oldMesh::getVertexCoords()  const
  {
    return _coords;
  }

  fvmc_nodal_t& oldMesh::getFvmNodal()
  {
    _finalizeNodal();
    return *_fvmNodal;
  }

  fvmc_nodal_t& oldMesh::getFvmNodal() const 
  {
    if (! _isNodalFinalized) {
      bftc_error(__FILE__, __LINE__, 0, "'%i' bad dimension\n", _nDim);
    }
    return *_fvmNodal;
  }

  const int& oldMesh::getNElts() const
  {
    return _nElts;
  }

  const int& oldMesh::getNPolyhedra() const
  {
    return _nPolyhedra;
  }

  const int* oldMesh::getEltConnectivityIndex() const
  {
    return _eltConnectivityIndex;
  }

  const int* oldMesh::getEltConnectivity() const
  {
    return _eltConnectivity;
  }

  const std::vector<double>& oldMesh::getVolume()
  {
    if (_cellVolume == NULL)
      _computeMeshProperties();
    return *_cellVolume;
  }

  const std::vector<double>& oldMesh::getCellCenterCoords()
  {
    if (_cellCenterCoords == NULL)
      _computeMeshProperties();

    return *_cellCenterCoords;
  }

  const std::vector<double>& oldMesh::getNormalFace()
  {
    if (_normalFace == NULL)
      _computeMeshProperties();

    return *_normalFace;
  }

  inline const std::vector<double>& oldMesh::getCharacteristicLength()
  {
    if (_characteristicLength == NULL)
      _computeMeshProperties();

    return *_characteristicLength;
  }

  inline const std::vector<int>& oldMesh::getIsDegenerated()
  {
    if (_isDegenerated == NULL)
      _computeMeshProperties();

    return *_isDegenerated;
  }

  const int *oldMesh::getPolyhedraFaceIndex() const
  {
    return _polyhedraFaceIndex;
  }

  const int *oldMesh::getPolyhedraCellToFaceConnectivity() const
  {
    return _polyhedraCellToFaceConnectivity;
  }

  const int *oldMesh::getPolyhedraFaceConnectivityIndex() const
  {
    return _polyhedraFaceConnectivityIndex;
  }

  const int *oldMesh::getPolyhedraFaceConnectivity() const
  {
    return _polyhedraFaceConnectivity;
  }

  const std::vector<int>& oldMesh::getPolyhedraCellToVertexConnectivity()
  {
    if (_polyhedraCellToVertexConnectivity == NULL)
      _computeMeshPolyhedraProperties();
    
    return *_polyhedraCellToVertexConnectivity;
  }
  
  const std::vector<int>& oldMesh::getPolyhedraCellToVertexConnectivityIndex()
  {
    if (_polyhedraCellToVertexConnectivityIndex == NULL)
      _computeMeshPolyhedraProperties();
    
    return *_polyhedraCellToVertexConnectivityIndex;
  }

}

#endif //__COUPLING_MESH_H__
