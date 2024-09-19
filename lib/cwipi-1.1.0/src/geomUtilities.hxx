#ifndef __CWIPI_GEOMUTILITIES_H__
#define __CWIPI_GEOMUTILITIES_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2012  ONERA

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

#include "cwipi_priv.h"
#include "vectorUtilities.hxx"

#include <algorithm>

namespace cwipi {

  /// Minimum value allowed for geometric computation
  
  const double GEOM_EPS_MIN  = 1e-30; 

  /// Constant value used to compute geomtric epsilon for volume

  const double GEOM_EPS_VOL  = 1e-9;
 
  /// Constant value used to compute geomtric epsilon for surface

  const double GEOM_EPS_SURF = 1e-9;

  /// Minimum distance between two vertices 

  const double GEOM_EPS_DIST = 1e-9; 

  ///
  /// \brief Compute a dynamic geometric epsilon from a characteristic length
  ///
  ///   @param [in]  characteristicLength  Characteristic length
  ///   @param [in]  consEpsilon           Constant part
  ///   @return                            Geometric epsilon
  ///

  inline double geometricEpsilon(const double characteristicLength,
                                 const double constEpsilon)
    
  {
    return std::max(constEpsilon * characteristicLength, GEOM_EPS_MIN);
  }


  ///
  /// \brief Triangle surface vector
  /// 
  /// @param [in]  nTriangle      Number of triangles
  /// @param [in]  connectivity   Connectivity
  /// @param [in]  coords         Vertice coordinates
  /// @param [out] surfaceVector  Surface Vector
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///

  inline void triangleSurfaceVector (const int     nTriangle,
                                     const int    *connectivity,
                                     const double *coords,
                                     double       *surfaceVector,
                                     double       *characteristicLength,
                                     int         *isDegenerated)

  {
    for (int itri = 0; itri < nTriangle; itri++) {

      const int i = connectivity[3*itri    ] - 1;
      const int j = connectivity[3*itri + 1] - 1;
      const int k = connectivity[3*itri + 2] - 1;

      double v1[3];
      double v2[3];
      double *surfaceVectorTri = surfaceVector + 3*itri;

      v1[0] = coords[3*j    ] - coords[3*i    ]; 
      v1[1] = coords[3*j + 1] - coords[3*i + 1]; 
      v1[2] = coords[3*j + 2] - coords[3*i + 2]; 

      v2[0] = coords[3*k    ] - coords[3*i    ]; 
      v2[1] = coords[3*k + 1] - coords[3*i + 1]; 
      v2[2] = coords[3*k + 2] - coords[3*i + 2]; 

      if (characteristicLength != NULL) {

        double v3[3];

        v3[0] = coords[3*k    ] - coords[3*j    ]; 
        v3[1] = coords[3*k + 1] - coords[3*j + 1]; 
        v3[2] = coords[3*k + 2] - coords[3*j + 2]; 

        double normV1 = norm(v1);
        double normV2 = norm(v2);
        double normV3 = norm(v3);

        characteristicLength[itri] = std::min(std::min(normV1, normV2), normV3);
 
      }

      crossProduct(v1, v2, surfaceVectorTri);

      surfaceVectorTri[0] *= 0.5;
      surfaceVectorTri[1] *= 0.5;
      surfaceVectorTri[2] *= 0.5;

      if ((characteristicLength != NULL) && (isDegenerated != NULL)) {

        double normSurfaceVectorTri = norm(surfaceVectorTri);
        double eps_loc = geometricEpsilon(characteristicLength[itri], GEOM_EPS_SURF);
        isDegenerated[itri] = 0;
        if (normSurfaceVectorTri <= eps_loc) 
          isDegenerated[itri] = 1;

      }
    }
  }

  ///
  /// \brief Triangle area
  /// 
  /// @param [in]  nTriangle      Number of triangles
  /// @param [in]  surfaceVector         SurfaceVector vectors
  /// @param [out] area           Area                
  ///

  inline void triangleArea (const int     nTriangle,
                            const double *surfaceVector,
                                  double *area)   

  {
    for (int itri = 0; itri < nTriangle; itri++) {
      const double *surfaceVectorTri = surfaceVector + 3*itri;
      area[itri] = norm(surfaceVectorTri);
    }
  }

  ///
  /// \brief Triangle center
  /// 
  /// @param [in]  nTriangle      Number of triangles
  /// @param [in]  connectivity   Connectivity
  /// @param [in]  coords         Vertice coordinates
  /// @param [out] center         center              
  ///

  inline void triangleCenter (const int     nTriangle,
                              const int    *connectivity,
                              const double *coords,
                                    double *center)
  {

    for (int itri = 0; itri < nTriangle; itri++) {

      const int i = connectivity[3*itri    ] - 1;
      const int j = connectivity[3*itri + 1] - 1;
      const int k = connectivity[3*itri + 2] - 1;

      double *centerTri = center + 3*itri;

      centerTri[0] = (coords[3*i    ] + coords[3*j    ] + coords[3*k    ])/3;
      centerTri[1] = (coords[3*i + 1] + coords[3*j + 1] + coords[3*k + 1])/3;
      centerTri[2] = (coords[3*i + 2] + coords[3*j + 2] + coords[3*k + 2])/3;
    }
  }

  ///
  /// \brief Tetrahedra oriented volume
  /// 
  /// @param [in]  nTetrahedra           Number of tetrahedra
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  coords                Vertice coordinates
  /// @param [out] volume                Volume              
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)
  ///

  inline void tetrahedraOrientedVolume (const int     nTetrahedra,
                                        const int    *connectivity,
                                        const double *coords,
                                        double       *volume,
                                        double       *characteristicLength,
                                        int         *isDegenerated)

  {
    
    for (int itet = 0; itet < nTetrahedra; itet++) {

      const int *connectivityTet = connectivity + 4*itet;

      const int i1 = connectivityTet[0] - 1;
      const int i2 = connectivityTet[1] - 1;
      const int i3 = connectivityTet[2] - 1;
      const int i4 = connectivityTet[3] - 1;

      const int nTriangle = 1;
      double surfaceVector[3];
      double faceCenter[3];
      double fC_i4[3];

      triangleSurfaceVector(nTriangle,
                            connectivityTet,
                            coords,
                            surfaceVector,
                            NULL,
                            NULL);

      if (characteristicLength != NULL) {

        double vectV1V2[3] =
          {coords[3*i2    ] - coords[3*i1    ],
           coords[3*i2 + 1] - coords[3*i1 + 1],
           coords[3*i2 + 2] - coords[3*i1 + 2]};

        double vectV1V3[3] =
          {coords[3*i3    ] - coords[3*i1    ],
           coords[3*i3 + 1] - coords[3*i1 + 1],
           coords[3*i3 + 2] - coords[3*i1 + 2]};

        double vectV2V3[3] =
          {coords[3*i3    ] - coords[3*i2    ],
           coords[3*i3 + 1] - coords[3*i2 + 1],
           coords[3*i3 + 2] - coords[3*i2 + 2]};

        double vectV1V4[3] =
          {coords[3*i4    ] - coords[3*i1    ],
           coords[3*i4 + 1] - coords[3*i1 + 1],
           coords[3*i4 + 2] - coords[3*i1 + 2]};

        double normV1V2 = norm(vectV1V2);
        double normV1V3 = norm(vectV1V3);
        double normV2V3 = norm(vectV2V3);
        double normV1V4 = norm(vectV1V4);

        characteristicLength[itet] = std::min(normV1V2,                    normV1V3);
        characteristicLength[itet] = std::min(characteristicLength[itet],  normV2V3);
        characteristicLength[itet] = std::min(characteristicLength[itet],  normV1V2);
        characteristicLength[itet] = std::min(characteristicLength[itet],  normV1V3);
        characteristicLength[itet] = std::min(characteristicLength[itet],  normV1V4);

      }

      triangleCenter(nTriangle,
                     connectivity,
                     coords,
                     faceCenter);

      fC_i4[0] = coords[3*i4    ] - faceCenter[0];
      fC_i4[1] = coords[3*i4 + 1] - faceCenter[1];
      fC_i4[2] = coords[3*i4 + 2] - faceCenter[2];

      volume[itet] = 1./3. * dotProduct(fC_i4, surfaceVector);

      if ((characteristicLength != NULL) && (isDegenerated != NULL)) {
     
        double eps_loc = geometricEpsilon(characteristicLength[itet], GEOM_EPS_VOL);
        isDegenerated[itet] = 0;
        if (volume[itet]  <= eps_loc) 
          isDegenerated[itet] = 1;

      }
    }
  }
  
  ///
  /// \brief Tetrahedra center
  /// 
  /// @param [in]  nTetrahedra    Number of tetrahedra
  /// @param [in]  connectivity   Connectivity
  /// @param [in]  coords         Vertice coordinates
  /// @param [out] center         center              
  ///

 inline void   tetrahedraCenter (const int     nTetrahedra,
                                 const int    *connectivity,
                                 const double *coords,
                                       double *center)

  {
    for (int itet = 0; itet < nTetrahedra; itet++) {

      const int i1 = connectivity[4*itet    ] - 1;
      const int i2 = connectivity[4*itet + 1] - 1;
      const int i3 = connectivity[4*itet + 2] - 1;
      const int i4 = connectivity[4*itet + 3] - 1;

      double *centerTet = center + 3*itet;

      centerTet[0] = (coords[3*i1    ] + coords[3*i2    ] 
                    + coords[3*i3    ] + coords[3*i4    ])/4;
      centerTet[1] = (coords[3*i1 + 1] + coords[3*i2 + 1] 
                    + coords[3*i3 + 1] + coords[3*i4 + 1])/4;
      centerTet[2] = (coords[3*i1 + 2] + coords[3*i2 + 2] 
                    + coords[3*i3 + 2] + coords[3*i4 + 2])/4;
    }
  }


  ///
  /// \brief Tetrahedra Faces
  /// 
  /// @param [in]  nTetrahedra       Number of tetrahedra
  /// @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
  /// @param [in]  connectivity      Connectivity
  /// @param [out] faceConnectivity  Face connectivity
  ///

  inline void tetrahedraFaces (const int     nTetrahedra,
                               const int     orientation,
                               const int    *connectivity,
                               int          *faceConnectivityIndex,
                               int          *faceConnectivity)
  {
    const int nFace       = 4;         // 4 triangles
    const int nVertex     = nFace * 3; // 4 * 3 vertices
    const int nVertexElt  = 4;         // 4 vertices
        
    faceConnectivityIndex[0] = 0;

    for (int ielt = 0; ielt < nTetrahedra; ielt++) {

      for (int iface = 0; iface < nFace; iface++) 
        faceConnectivityIndex[ielt * nFace + iface + 1] = 
          faceConnectivityIndex[ielt * nFace + iface] + 3;

      if (orientation == 0) {

        faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 2];

        faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 1];
        
        faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 3];
        
        faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 2];

      }
        
      else {

        faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt    ];
        
        faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt    ];
        
        faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt    ];
        
        faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 1];

      }
    }
  }


  ///
  /// \brief HexahedraFaces
  /// 
  /// @param [in]  nHexahedra        Number of hexahedra
  /// @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
  /// @param [in]  connectivity      Connectivity
  /// @param [out] faceConnectivity  Face connectivity
  ///

  inline void hexahedraFaces (const int     nHexahedra,
                              const int     orientation,
                              const int    *connectivity,
                              int          *faceConnectivityIndex,
                              int          *faceConnectivity)
  {
    const int nFace       = 6;         // 6 quadrangles
    const int nVertex     = nFace * 4; // 6 * 4 vertices
    const int nVertexElt  = 8;         // 8 vertices
        
    faceConnectivityIndex[0] = 0;

    for (int ielt = 0; ielt < nHexahedra; ielt++) {

      for (int iface = 0; iface < nFace; iface++) 
        faceConnectivityIndex[ielt * nFace + iface + 1] = 
          faceConnectivityIndex[ielt * nFace + iface] + 4;

      if (orientation == 0) {

        faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 3];
        
        faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 7];        
        faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 6];
        
        faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 3];        
        faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 7];
        faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 4];
        
        faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 3];        
        faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 6];
        faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 7];
    
        faceConnectivity[nVertex * ielt + 16] = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 17] = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 18] = connectivity[nVertexElt * ielt + 6];
        faceConnectivity[nVertex * ielt + 19] = connectivity[nVertexElt * ielt + 2];
        
        faceConnectivity[nVertex * ielt + 20] = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 21] = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 22] = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 23] = connectivity[nVertexElt * ielt + 1];

      }

      else {

        faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt    ];
        
        faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 6];
        faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 7];
        faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 4];        
        faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 5];
        
        faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 7];        
        faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt    ];
        
        faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 7];        
        faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 6];
        faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 3];
    
        faceConnectivity[nVertex * ielt + 16] = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 17] = connectivity[nVertexElt * ielt + 6];
        faceConnectivity[nVertex * ielt + 18] = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 19] = connectivity[nVertexElt * ielt + 1];
        
        faceConnectivity[nVertex * ielt + 20] = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 21] = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 22] = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 23] = connectivity[nVertexElt * ielt + 0];

      }
    }
  }


  ///
  /// \brief Prism Faces
  /// 
  /// @param [in]  nPrism            Number of Prism
  /// @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
  /// @param [in]  connectivity      Connectivity
  /// @param [out] faceConnectivity  Face connectivity
  ///

  inline void prismFaces (const int     nPrism,     
                          const int     orientation,
                          const int    *connectivity,
                          int          *faceConnectivityIndex,
                          int          *faceConnectivity)
  {
    const int nFace       = 5;         // 3 quadrangles + 2 triangles
    const int nVertex     = 3*4 + 2*3; // 
    const int nVertexElt  = 6;         // 6 vertices
        
    faceConnectivityIndex[0] = 0;

    for (int ielt = 0; ielt < nPrism; ielt++) {

      faceConnectivityIndex[ielt * nFace + 1] = faceConnectivityIndex[ielt * nFace    ] + 3;
      faceConnectivityIndex[ielt * nFace + 2] = faceConnectivityIndex[ielt * nFace + 1] + 3;
      faceConnectivityIndex[ielt * nFace + 3] = faceConnectivityIndex[ielt * nFace + 2] + 4;
      faceConnectivityIndex[ielt * nFace + 4] = faceConnectivityIndex[ielt * nFace + 3] + 4;
      faceConnectivityIndex[ielt * nFace + 5] = faceConnectivityIndex[ielt * nFace + 4] + 4;

      if (orientation == 0) {

        faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 2];
        
        faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 3];        
        faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 4];
        
        faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 2];        
        faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 5];
        
        faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 3];        
        faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 4];
        
        faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 16] = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 17] = connectivity[nVertexElt * ielt + 3];

      }

      else {

        faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt    ];
        
        faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 4];        
        faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt + 3];
        
        faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 5];        
        faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 2];
        
        faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt    ];        
        faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 1];
        
        faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 5];
        faceConnectivity[nVertex * ielt + 16] = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 17] = connectivity[nVertexElt * ielt    ];

      }
    }
  }

  ///
  /// \brief Pyramid Faces
  /// 
  /// @param [in]  nPyramid          Number of pyramid  
  /// @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
  /// @param [in]  connectivity      Connectivity
  /// @param [out] faceConnectivity  Face connectivity
  ///

  inline void pyramidFaces (const int     nPyramid,
                            const int     orientation,
                            const int    *connectivity,
                            int          *faceConnectivityIndex,
                            int          *faceConnectivity)
  {
    const int nFace       = 5;         // 1 quadrangle + 4 triangles
    const int nVertex     = 1*4 + 4*3; // 
    const int nVertexElt  = 5;         // 5 vertices
        
    faceConnectivityIndex[0] = 0;

    for (int ielt = 0; ielt < nPyramid; ielt++) {

      faceConnectivityIndex[ielt * nFace + 1] = faceConnectivityIndex[ielt * nFace    ] + 4;
      faceConnectivityIndex[ielt * nFace + 2] = faceConnectivityIndex[ielt * nFace + 1] + 3;
      faceConnectivityIndex[ielt * nFace + 3] = faceConnectivityIndex[ielt * nFace + 2] + 3;
      faceConnectivityIndex[ielt * nFace + 4] = faceConnectivityIndex[ielt * nFace + 3] + 3;
      faceConnectivityIndex[ielt * nFace + 5] = faceConnectivityIndex[ielt * nFace + 4] + 3;

      if (orientation == 0) {

        faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt + 3];
        
        faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 4];
        
        faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 4];
    
        faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 4];
        
        faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt + 4];

      }

      else {

        faceConnectivity[nVertex * ielt + 0]  = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 1]  = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 2]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 3]  = connectivity[nVertexElt * ielt    ];
        
        faceConnectivity[nVertex * ielt + 4]  = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 5]  = connectivity[nVertexElt * ielt    ];
        faceConnectivity[nVertex * ielt + 6]  = connectivity[nVertexElt * ielt + 1];
        
        faceConnectivity[nVertex * ielt + 7]  = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 8]  = connectivity[nVertexElt * ielt + 1];
        faceConnectivity[nVertex * ielt + 9]  = connectivity[nVertexElt * ielt + 2];
    
        faceConnectivity[nVertex * ielt + 10] = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 11] = connectivity[nVertexElt * ielt + 2];
        faceConnectivity[nVertex * ielt + 12] = connectivity[nVertexElt * ielt + 3];
        
        faceConnectivity[nVertex * ielt + 13] = connectivity[nVertexElt * ielt + 4];
        faceConnectivity[nVertex * ielt + 14] = connectivity[nVertexElt * ielt + 3];
        faceConnectivity[nVertex * ielt + 15] = connectivity[nVertexElt * ielt    ];

      }
    }
  }


#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable:869)
#elif defined(__clang__)	
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value" 	
#endif

  ///
  /// \brief Edges properties 
  /// 
  /// @param [in]  nEdges                Number of edges     
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  nVertices             Number of vertices
  /// @param [in]  coords                Vertices coordinates
  /// @param [out] center                Center              
  /// @param [out] length                Length              
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///

  inline void edgesProperties (const int     nEdges,
                               const int    *connectivity,
                               const int     nVertices,
                               const double *coords,
                               double       *center,
                               double       *length,
                               double       *characteristicLength,
                               int         *isDegenerated)

  {

    CWIPI_UNUSED (nVertices);

    for (int iedge = 0; iedge < nEdges; iedge++) {

      const int *connectivityEdge = connectivity + 2*iedge;

      const int i1 = connectivityEdge[0] - 1;
      const int i2 = connectivityEdge[1] - 1;

      double *centerEdge = center + 3*iedge;

      for (int i = 0; i < 3; i++) 
        centerEdge[i] = 0.5 * (coords[3*i1 + i] + coords[3*i2 + i]);

      length[iedge] = sqrt((coords[3*i2 + 0] - coords[3*i1 + 0]) 
                           * (coords[3*i2 + 0] - coords[3*i1 + 0])
                         + (coords[3*i2 + 1] - coords[3*i1 + 1]) 
                           * (coords[3*i2 + 1] - coords[3*i1 + 1])   
                         + (coords[3*i2 + 2] - coords[3*i1 + 2]) 
                           * (coords[3*i2 + 2] - coords[3*i1 + 2]));

      if (characteristicLength != NULL) 
        characteristicLength[iedge] = length[iedge];

      if (isDegenerated != NULL) {
        isDegenerated[iedge] = 0;
        if (length[iedge] < GEOM_EPS_DIST)
          isDegenerated[iedge] = 1;
      }
    }    
  }
  
  ///
  /// \brief Triangle properties
  /// 
  /// @param [in]  nTriangle             Number of triangles
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  nVertices             Number of vertices
  /// @param [in]  coords                Vertices coordinates
  /// @param [out] surfaceVector         Surface vector
  /// @param [out] center                Center              
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///

  inline void triangleProperties (const int     nTriangle,
                                  const int    *connectivity,
                                  const int     nVertices,
                                  const double *coords,
                                  double       *surfaceVector,
                                  double       *center,
                                  double       *characteristicLength,
                                  int         *isDegenerated)

  {
   
    CWIPI_UNUSED(nVertices);

    triangleSurfaceVector (nTriangle,
                           connectivity,
                           coords,
                           surfaceVector,
                           characteristicLength,
                           isDegenerated);

    triangleCenter (nTriangle,
                    connectivity,
                    coords,
                    center);

  }

#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#elif defined(__clang__)	
#pragma clang diagnostic pop
#endif

  ///
  /// \brief Quadrangle properties
  /// 
  /// @param [in]  nTriangle             Number of quadrangles
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  nVertices             Number of vertices
  /// @param [in]  coords                Vertices coordinates
  /// @param [out] surfaceVector         Surface vector
  /// @param [out] center                Center              
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///
  /// @return                     The status of properties computation convergence                
  ///

  int quadrangleProperties (const int     nQuadrangle,
                             const int    *connectivity,
                             const int     nVertices,
                             const double *coords,
                             double       *surfaceVector,
                             double       *center,
                             double       *characteristicLength,
                             int         *isDegenerated);


  ///
  /// \brief Polygon properties
  /// 
  /// @param [in]  nPolygon              Number of polygon
  /// @param [in]  connectivityIndex     Connectivity Index
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  nVertices             Number of vertices
  /// @param [in]  coords                Vertices coordinates
  /// @param [out] surfaceVector         Surface vector
  /// @param [out] center                Center
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///
  /// @return                        The status of properties computation convergence
  ///

  int polygonProperties (const int     nPolygon,   
                          const int    *connectivityIndex,
                          const int    *connectivity,
                          const int     nVertices,
                          const double *coords,
                          double       *surfaceVector,
                          double       *center,
                          double       *characteristicLength,
                          int         *isDegenerated);



#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable:869)
#elif defined(__clang__)	
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value" 	
#endif

  ///
  /// \brief Tetrahedra properties
  /// 
  /// @param [in]  nTetrahedra           Number of tetrahedra 
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  nVertices             Number of vertices
  /// @param [in]  coords                Vertices coordinates
  /// @param [out] volume                Volume
  /// @param [out] center                Center              
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///


  inline void tetrahedraProperties (const int     nTetrahedra,
                                    const int    *connectivity,
                                    const int     nVertices,
                                    const double *coords,
                                    double       *volume,
                                    double       *center,
                                    double       *characteristicLength,
                                    int         *isDegenerated)

  {

    CWIPI_UNUSED(nVertices);

    tetrahedraOrientedVolume (nTetrahedra,
                              connectivity,
                              coords,
                              volume,
                              characteristicLength,
                              isDegenerated);
    
    tetrahedraCenter (nTetrahedra,
                      connectivity,
                      coords,
                      center);
  }
  
#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#elif defined(__clang__)	
#pragma clang diagnostic pop
#endif
  
  ///
  /// \brief Hexahedra properties
  /// 
  /// @param [in]  nHexahedra            Number of hexahedra  
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  nVertices             Number of vertices
  /// @param [in]  coords                Vertices coordinates
  /// @param [out] volume                Volume
  /// @param [out] center                Center              
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///

  void hexahedraProperties (const int     nHexahedra,   
                            const int    *connectivity,
                            const int     nVertices,
                            const double *coords,
                            double       *volume,
                            double       *center,
                            double       *characteristicLength,
                            int         *isDegenerated);


  ///
  /// \brief Prism properties
  /// 
  /// @param [in]  nPrism                Number of prism      
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  nVertices             Number of vertices
  /// @param [in]  coords                Vertices coordinates
  /// @param [out] volume                Volume
  /// @param [out] center                Center              
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///

  void prismProperties (const int     nPrism,     
                        const int    *connectivity,
                        const int     nVertices,
                        const double *coords,
                        double       *volume,
                        double       *center,
                        double       *characteristicLength,
                        int         *isDegenerated);


  ///
  /// \brief Pyramid properties
  /// 
  /// @param [in]  nPyramid              Number of pyramid    
  /// @param [in]  connectivity          Connectivity
  /// @param [in]  nVertices             Number of vertices
  /// @param [in]  coords                Vertices coordinates
  /// @param [out] volume                Volume
  /// @param [out] center                Center              
  /// @param [out] characteristicLength  Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated         Degenerated edge indicator (active if != NULL)             
  ///

  void pyramidProperties (const int     nPyramid,   
                          const int    *connectivity,
                          const int     nVertices,
                          const double *coords,
                          double       *volume,
                          double       *center,
                          double       *characteristicLength,
                          int         *isDegenerated);


  ///
  /// \brief Polyhedra properties
  /// 
  /// @param [in]  nPolyhedra                 Number of polyhedra  
  /// @param [in]  nFace                      Number of faces      
  /// @param [in]  faceConnectivityIdx        Face connectivity index
  /// @param [in]  faceConnectivity           Face connectivity
  /// @param [in]  cellToFaceConnectivityIdx  Cell to face connectivity index
  /// @param [in]  cellToFaceConnectivity     Cell to face connectivity
  /// @param [in]  nVertices                  Number of vertices
  /// @param [in]  coords                     Vertices coordinates
  /// @param [out] volume                     Volume
  /// @param [out] center                     Center              
  /// @param [out] characteristicLength       Characteristic length (active if != NULL)             
  /// @param [out] isDegenerated              Degenerated edge indicator (active if != NULL)
  ///

  void polyhedraProperties (const int     nPolyhedra,
                            const int     nFace,
                            const int    *faceConnectivityIdx,
                            const int    *faceConnectivity,   
                            const int    *cellToFaceConnectivityIdx,
                            const int    *cellToFaceConnectivity,
                            const int     nVertices,
                            const double *coords,
                            double       *volume,
                            double       *center,
                            double       *characteristicLength,
                            int         *isDegenerated);


  //inline double distPointOnSurfaceTriangle();
  
  //inline int isCoplanarPointOnSurfaceTriangle();

  //inline int isCoplanarPointOnSurfacePolygon();

  //inline int isPointInVolumeTetra();


}

#endif //__CWIPI_GEOMUTILITIES_H__

