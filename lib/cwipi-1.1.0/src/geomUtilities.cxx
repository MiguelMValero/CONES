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

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include <bftc_printf.h>

#include <geomUtilities.hxx>

namespace cwipi {

  ///
  /// \brief Quadrangle properties
  /// 
  /// @param [in]  nQuadrangle           Number of quadrangles
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
                             int         *isDegenerated)
  {
    int *connectivityIndex =  (int *) malloc (sizeof(int) * (nQuadrangle + 1));
    int convergence;

    connectivityIndex[0] = 0;
    for (int i = 1; i < nQuadrangle + 1; i++)
      connectivityIndex[i] = connectivityIndex[i-1] + 4; 

    convergence = polygonProperties(nQuadrangle,
                                    connectivityIndex,
                                    connectivity,
                                    nVertices,
                                    coords,
                                    surfaceVector,
                                    center,
                                    characteristicLength,
                                    isDegenerated);

    free ( connectivityIndex);

    return convergence;
  }


#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable:869)
#elif defined(__clang__)	
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value" 	
#endif

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
                         int          *isDegenerated)

  {
    CWIPI_UNUSED (nVertices);

    int convergence = 1;
    
    const double dispMin = 1e-9; // Minimum displacement
    const double big = 1e30;     // Big value;

    const int nIterMax = 100;    // Maximum iteration number 

    for (int ifac = 0; ifac < nPolygon; ifac++) {

      //
      // Define local pointer
      // --------------------

      const int nVerticesFace = connectivityIndex[ifac + 1] 
                              - connectivityIndex[ifac    ];

      const int *connectivityFace = connectivity + connectivityIndex[ifac];

      double *surfaceVectorFace = surfaceVector + 3*ifac; // Face surface vector

      double *centerFace = center + 3*ifac; // Face Center

      if (characteristicLength != NULL) {
        characteristicLength[ifac] = big;
      }

      //
      // Initialization
      // --------------

      int nIter = 0;       // Number of iteration

      for (int i = 0; i < 3; i++) {
        surfaceVectorFace[i] = 0;
        centerFace[i]        = 0;
      }
      
      // Initialize face center to the barycenter

      for (int ivert = 0; ivert < nVerticesFace; ivert++) {
        const int vert = connectivityFace[ivert] - 1;  
        for (int i = 0; i < 3; i++) 
          centerFace[i] += coords[3*vert + i];
      }

      for (int i = 0; i < 3; i++) 
        centerFace[i] /= nVerticesFace;

      //
      // Compute cell center and surface vector
      // --------------------------------------

      while (1) {

        double displacement[3] = {0, 0, 0}; // Cell center displacement

        nIter += 1;

        for (int i = 0; i < 3; i++)
          surfaceVectorFace[i] = 0;

        double areaFace = 0; // Face area

        for (int ivert = 0; ivert < nVerticesFace; ivert++) {

          const int vert1 = connectivityFace[ivert] - 1;  
          const int vert2 = connectivityFace[(ivert + 1) % nVerticesFace] - 1;

          // Edge center

          const double edgeCenter[3] = {0.5 * (coords[3*vert1    ] + coords[3*vert2    ]),
                                        0.5 * (coords[3*vert1 + 1] + coords[3*vert2 + 1]),
                                        0.5 * (coords[3*vert1 + 2] + coords[3*vert2 + 2])};

          // Edge vector

          const double vectV1V2[3] = {coords[3*vert2    ] - coords[3*vert1    ],
                                      coords[3*vert2 + 1] - coords[3*vert1 + 1],
                                      coords[3*vert2 + 2] - coords[3*vert1 + 2]};

          if (characteristicLength != NULL) {
            double norm_V1V2 = norm(vectV1V2);
            characteristicLength[ifac] = std::min(characteristicLength[ifac], norm_V1V2);
          }

          // Vector face center -> edge center 

          const double vectFECenter[3] = {edgeCenter[0] - centerFace[0],
                                          edgeCenter[1] - centerFace[1],
                                          edgeCenter[2] - centerFace[2]};

          // Compute area of the triangle (face center, vert1, vert2) 

          double surfaceVectorTria[3];
          crossProduct(vectFECenter, vectV1V2, surfaceVectorTria);

          for (int i = 0; i < 3; i++) 
             surfaceVectorTria[i] *= 0.5;
          
          const double areaTri = norm(surfaceVectorTria);

          areaFace += areaTri;
          for (int i = 0; i < 3; i++) {
            surfaceVectorFace[i] += surfaceVectorTria[i];
            displacement[i] += areaTri * vectFECenter[i];
          }

          areaFace += areaTri;

        }

        double denomAreaFace = 1. / std::max(fabs(areaFace), GEOM_EPS_MIN);

        for (int i = 0; i < 3; i++) {
          displacement[i] = 2./3. * denomAreaFace * displacement[i];
          centerFace[i] += displacement[i];
        }

        //
        // Check convergence

        const double normDisp = norm(displacement);

        if (normDisp < dispMin) {
          break;
        }

        //
        // Check Number of iteration

        else if (nIterMax < nIter) {
          convergence = false;
          break;
        }
      } // while (1)

      if ((characteristicLength != NULL) && (isDegenerated != NULL)) {
     
        double normSurfaceVector = norm(surfaceVectorFace);
        double eps_loc = geometricEpsilon(characteristicLength[ifac], GEOM_EPS_SURF);
        isDegenerated[ifac] = 0;
        if (normSurfaceVector <= eps_loc) 
          isDegenerated[ifac] = 1;
      }
    } // for (int ifac = 0; ifac < nPolygon; ifac++)
 
    if (0 == 1) {
      bftc_printf("surfacevector : ");
      for (int ifac = 0; ifac < 3*nPolygon; ifac++) {
        bftc_printf("%12.5e ",surfaceVector[ifac]);
      }
      bftc_printf("\n");
      
      bftc_printf("center : ");
      for (int ifac = 0; ifac < 3*nPolygon; ifac++) {
        bftc_printf("%12.5e ",center[ifac]);
      }
      bftc_printf("\n");
    }
    return convergence;
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
                            int         *isDegenerated)

  {

    const int orientation = 1; //  Surface vector oriented towards inside cell outside

    const int nQuadrangle = 6;
    const int nTriangle = 0;
    const int nHexahedraFaces = nQuadrangle + nTriangle;
    const int nFaces = nHexahedraFaces * nHexahedra;

    int *faceConnectivity          =  (int *) malloc (sizeof(int) * ((nQuadrangle*4 + nTriangle*3) * nHexahedra));
    int *faceConnectivityIdx       =  (int *) malloc (sizeof(int) * (nFaces + 1));
    int *cellToFaceConnectivityIdx =  (int *) malloc (sizeof(int) * (nHexahedra + 1));
    int *cellToFaceConnectivity    =  (int *) malloc (sizeof(int) * (nFaces));

    //
    // Get hexahedra faces

    hexahedraFaces (nHexahedra,
                    orientation,
                    connectivity,
                    faceConnectivityIdx,
                    faceConnectivity);

    //
    // Define cell to face connectivity

    for (int i = 0; i < nFaces; i++)
      cellToFaceConnectivity[i] = i + 1;

    cellToFaceConnectivityIdx[0] = 0;
    for (int i = 1; i < nHexahedra + 1; i++)
      cellToFaceConnectivityIdx[i] = cellToFaceConnectivityIdx[i-1] + nHexahedraFaces;
    
    //
    // Compute Volume and center 
    
    polyhedraProperties (nHexahedra,
                         nFaces,
                         faceConnectivityIdx,
                         faceConnectivity,   
                         cellToFaceConnectivityIdx,
                         cellToFaceConnectivity,
                         nVertices,
                         coords,
                         volume,
                         center,
                         characteristicLength,
                         isDegenerated);

    
    //
    // Free
    
    free ( faceConnectivity);
    free ( faceConnectivityIdx);
    free ( cellToFaceConnectivity);
    free ( cellToFaceConnectivityIdx);

  }



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
                        int         *isDegenerated)

  {

    const int orientation = 1; //  Surface vector oriented towards inside cell outside

    const int nQuadrangle = 3;
    const int nTriangle = 2;
    const int nPrismFaces = nQuadrangle + nTriangle;
    const int nFaces = nPrismFaces * nPrism;

    int *faceConnectivity          =  (int *) malloc (sizeof(int) * ((nQuadrangle*4 + nTriangle*3) * nPrism));
    int *faceConnectivityIdx       =  (int *) malloc (sizeof(int) * (nFaces + 1));
    int *cellToFaceConnectivityIdx =  (int *) malloc (sizeof(int) * (nPrism + 1));
    int *cellToFaceConnectivity    =  (int *) malloc (sizeof(int) * (nFaces));

    //
    // Get prism faces

    prismFaces (nPrism,
                orientation,
                connectivity,
                faceConnectivityIdx,
                faceConnectivity);

    //
    // Define cell to face connectivity

    for (int i = 0; i < nFaces; i++)
      cellToFaceConnectivity[i] = i + 1;

    cellToFaceConnectivityIdx[0] = 0;
    for (int i = 1; i < nPrism + 1; i++)
      cellToFaceConnectivityIdx[i] = cellToFaceConnectivityIdx[i-1] + nPrismFaces;
    
    //
    // Compute Volume and center 
    
    polyhedraProperties (nPrism,
                         nFaces,
                         faceConnectivityIdx,
                         faceConnectivity,   
                         cellToFaceConnectivityIdx,
                         cellToFaceConnectivity,
                         nVertices,
                         coords,
                         volume,
                         center,
                         characteristicLength,
                         isDegenerated);

    
    //
    // Free
    
    free ( faceConnectivity);
    free ( faceConnectivityIdx);
    free ( cellToFaceConnectivity);
    free ( cellToFaceConnectivityIdx);
    
  }

  
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
                          int         *isDegenerated)

  {
    const int orientation = 1; //  Surface vector oriented towards inside cell outside

    const int nQuadrangle = 1;
    const int nTriangle = 4;
    const int nPyramidFaces = nQuadrangle + nTriangle;
    const int nFaces = nPyramidFaces * nPyramid;

    int *faceConnectivity          =  (int *) malloc (sizeof(int) * ((nQuadrangle*4 + nTriangle*3) * nPyramid));
    int *faceConnectivityIdx       =  (int *) malloc (sizeof(int) * (nFaces + 1));
    int *cellToFaceConnectivityIdx =  (int *) malloc (sizeof(int) * (nPyramid + 1));
    int *cellToFaceConnectivity    =  (int *) malloc (sizeof(int) * (nFaces));

    //
    // Get pyramid faces

    pyramidFaces (nPyramid,
                  orientation,
                  connectivity,
                  faceConnectivityIdx,
                  faceConnectivity);

    //
    // Define cell to face connectivity

    for (int i = 0; i < nFaces; i++)
      cellToFaceConnectivity[i] = i + 1;

    cellToFaceConnectivityIdx[0] = 0;
    for (int i = 1; i < nPyramid + 1; i++)
      cellToFaceConnectivityIdx[i] = cellToFaceConnectivityIdx[i-1] + nPyramidFaces;
    
    //
    // Compute Volume and center 

    polyhedraProperties (nPyramid,
                         nFaces,
                         faceConnectivityIdx,
                         faceConnectivity,   
                         cellToFaceConnectivityIdx,
                         cellToFaceConnectivity,
                         nVertices,
                         coords,
                         volume,
                         center,
                         characteristicLength,
                         isDegenerated);
    
    //
    // Free
    
    free ( faceConnectivity);
    free ( faceConnectivityIdx);
    free ( cellToFaceConnectivity);
    free ( cellToFaceConnectivityIdx);
    
  }


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
                            int         *isDegenerated)

  {
    const double big = 1e30;
    int convergence = 1;
    int *colorVertice =  (int *) malloc (sizeof(int) * (nVertices));

    int  lPolyhedraVertices = 24;
    std::vector <int>  polyhedraVertices(lPolyhedraVertices); // First

    // int *tmpCellToFaceConnectivity = NULL;
    // if (nPolyhedra > 0 )
    //   tmpCellToFaceConnectivity =  (int *) malloc (sizeof(int) * (cellToFaceConnectivityIdx[nPolyhedra]));

    for (int i = 0; i < nVertices; i++)
      colorVertice[i] = false;

    //
    // Compute face properties

    double *surfaceVector = (double *) malloc (sizeof(double) *(3 * nFace)); 
    double *faceCenter    = (double *) malloc (sizeof(double) *(3 * nFace)); 

    int convergenceFace = polygonProperties (nFace,
                                             faceConnectivityIdx,
                                             faceConnectivity,
                                             nVertices,
                                             coords,
                                             surfaceVector,
                                             faceCenter,
                                             NULL,
                                             NULL);

    if (0 == 1) {

      bftc_printf("faceConnectivity : \n");
      for (int ipoly = 0; ipoly < nFace; ipoly++) {
        bftc_printf("  - face %i : ", ipoly+1);
        for (int j = faceConnectivityIdx[ipoly]; j < faceConnectivityIdx[ipoly+1]; j++) {
          bftc_printf("%i ",faceConnectivity[j]);
        }
        bftc_printf("\n");
      }

      bftc_printf("surfacevector : ");
      for (int ipoly = 0; ipoly < 3 * nFace; ipoly++) {
        bftc_printf("%12.5e ",surfaceVector[ipoly]);
      }
      bftc_printf("\n");
      
      bftc_printf("facecenter : ");
      for (int ipoly = 0; ipoly < 3 * nFace; ipoly++) {
        bftc_printf("%12.5e ",faceCenter[ipoly]);
      }
      bftc_printf("\n");
    }

    //
    // Loop on polyhedra

//    double volume_t =0;
    for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {

      double *polyCenter = center + 3*ipoly;
      double  disp[3];

      const int polyIdx   = cellToFaceConnectivityIdx[ipoly];
      const int nPolyFace = cellToFaceConnectivityIdx[ipoly + 1] - polyIdx;

      if (characteristicLength != NULL)
        characteristicLength[ipoly] = big;

      //
      // Intialize cell center to the barycenter

      // Search polyhedra vertices

      int nPolyhedraVertices = 0;
      for (int iface = 0; iface < nPolyFace; iface++) {

        const int face          = abs(cellToFaceConnectivity[polyIdx + iface]) - 1;
        const int faceIdx       = faceConnectivityIdx[face];
        const int nFaceVertices = faceConnectivityIdx[face+1] - faceIdx;

        for (int ivert = 0; ivert < nFaceVertices; ivert++) {
          const int vertex = faceConnectivity[faceIdx + ivert] - 1;
          if (!colorVertice[vertex]) {
            colorVertice[vertex] = 1;

            if (nPolyhedraVertices >= (int) polyhedraVertices.size())
              polyhedraVertices.resize(2 * polyhedraVertices.size());
            polyhedraVertices[nPolyhedraVertices++] = vertex;
          }
        }

      }

      // Compute barycenter

      for (int i = 0; i < 3; i++)
        polyCenter[i] = 0.;
      
      for (int j = 0; j < nPolyhedraVertices; j++) {
        const int vertex = polyhedraVertices[j];
        colorVertice[vertex] = false;
        for (int i = 0; i < 3; i++)
          polyCenter[i] += coords[3*vertex + i];
      }

      for (int i = 0; i < 3; i++)
        polyCenter[i] /= nPolyhedraVertices;

      nPolyhedraVertices = 0;


      //
      // Intialize volume

      volume[ipoly] = 0.;

      //
      // Intialize cell center displacement

      for (int i = 0; i < 3; i++)
         disp[i] = 0.;
 
      //
      // Loop on faces

      for (int iface = 0; iface < nPolyFace; iface++) {

        const int face          = abs(cellToFaceConnectivity[polyIdx + iface]) - 1;
        const int direction     = (cellToFaceConnectivity[polyIdx + iface] < 0) ? -1 : 1;

        const int faceIdx       = faceConnectivityIdx[face];
        const int nFaceVertices = faceConnectivityIdx[face+1] - faceIdx;

        //
        // Loop on vertices

        for (int ivert = 0; ivert < nFaceVertices; ivert++) {
          
          int vert1;
          int vert2;
          int ivert1 = ivert;

          if (direction > 0) {
            vert1 = faceConnectivity[faceIdx + ivert1] - 1;
            vert2 = faceConnectivity[faceIdx + (ivert1 + 1) % nFaceVertices] - 1;
          }
          else {
            ivert1 = nFaceVertices - 1 - ivert; 
            vert1 = faceConnectivity[faceIdx + (ivert1 + 1) % nFaceVertices] - 1;
            vert2 = faceConnectivity[faceIdx + ivert1] - 1;
          }

          if (characteristicLength != NULL)  {

            // Vector vert1 -> vert2

            const double vectV1V2[3] = 
              {coords[3*vert2    ] - coords[3*vert1    ],
               coords[3*vert2 + 1] - coords[3*vert1 + 1],
               coords[3*vert2 + 2] - coords[3*vert1 + 2]};

            double normV1V2 = norm (vectV1V2);

            characteristicLength[ipoly] = std::min(characteristicLength[ipoly], normV1V2);

          }

          // Vector face center -> vert1

          const double vectFCV1[3] = 
            {coords[3*vert1    ] - faceCenter[3 * face    ],
             coords[3*vert1 + 1] - faceCenter[3 * face + 1],
             coords[3*vert1 + 2] - faceCenter[3 * face + 2]};
            
          // Vector face center -> vert2

          const double vectFCV2[3] = 
            {coords[3*vert2    ] - faceCenter[3 * face    ],
             coords[3*vert2 + 1] - faceCenter[3 * face + 1],
             coords[3*vert2 + 2] - faceCenter[3 * face + 2]};

          // Vector cell center -> face center

          const double vectCCFC[3] = 
            {faceCenter[3 * face    ] - polyCenter[0],
             faceCenter[3 * face + 1] - polyCenter[1],
             faceCenter[3 * face + 2] - polyCenter[2]};

          double surfaceVectorTri[3];
          crossProduct (vectFCV1, vectFCV2, surfaceVectorTri);

          for (int i = 0; i < 3; i++)
            surfaceVectorTri[i] *= 0.5;
          
          // Oriented volume

          double volumeTet = 1./3 * dotProduct (surfaceVectorTri, vectCCFC);

          volume[ipoly] += volumeTet;

          for (int i = 0; i < 3; i++)
            disp[i] = disp[i] + volumeTet * vectCCFC[i];

        }

      }

      int signeVol = (volume[ipoly] < 0.) ? -1 : 1 ;

      if (signeVol == -1) {
        bftc_printf("Warning polyhedraProperties : volume < 0 for polyhedron '%i' : %12.5e\n", 
                    ipoly + 1, volume[ipoly]);
      } 

      double denomVol = 1 / std::max(fabs(volume[ipoly]), GEOM_EPS_MIN);
      
      for (int i = 0; i < 3; i++)
        polyCenter[i] =  polyCenter[i] + signeVol * denomVol * disp[i];

//      volume_t += volume[ipoly];

      //
      // Check convergence

      if (!convergenceFace)
        convergence = false;

      if ((characteristicLength != NULL) && (isDegenerated != NULL)) {
     
        double eps_loc = geometricEpsilon(characteristicLength[ipoly], GEOM_EPS_VOL);
        isDegenerated[ipoly] = 0;
        if (fabs(volume[ipoly]) <= eps_loc) 
          isDegenerated[ipoly] = 1;
      }

    }
    // if (nPolyhedra > 0) {
    //   bftc_printf("connec : ");
    //   for (int i = 0; i < cellToFaceConnectivityIdx[nPolyhedra]; i++)
    //     bftc_printf("%i ", tmpCellToFaceConnectivity[i]);
    //   bftc_printf("\n");
    //   free ( tmpCellToFaceConnectivity);
    // }

    if (0 == 1) {

      bftc_printf("surfacevector : ");
      for (int ipoly = 0; ipoly < 3 * nFace; ipoly++) {
        bftc_printf("%12.5e ",surfaceVector[ipoly]);
      }
      bftc_printf("\n");

      bftc_printf("facecenter : ");
      for (int ipoly = 0; ipoly < 3 * nFace; ipoly++) {
        bftc_printf("%12.5e ",faceCenter[ipoly]);
      }
      bftc_printf("\n");

      bftc_printf("isDegenrated : ");
      for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {
        bftc_printf("%i ",isDegenerated[ipoly]);
      }
      bftc_printf("\n");

      bftc_printf("characteristicLength  : ");
      for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {
        bftc_printf("%12.5e ",characteristicLength[ipoly]);
      }
      bftc_printf("\n");

      bftc_printf("volume  : ");
      for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {
        bftc_printf("%12.5e ",volume[ipoly]);
      }
      bftc_printf("\n");

      bftc_printf("center  : ");
      for (int ipoly = 0; ipoly < 3 * nPolyhedra; ipoly++) {
        bftc_printf("%12.5e ",center[ipoly]);
      }
      bftc_printf("\n");

    }

    free (surfaceVector); 
    free ( faceCenter);
    free ( colorVertice);

    if (!convergence)
      bftc_printf("Warning polyhedraProperties : some polyhedra faces are not planar\n");
    
  }
}
