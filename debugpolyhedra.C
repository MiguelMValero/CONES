#include <stdlib.h>
#include <iostream>

int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{

const int nElemSeg = 2;
const int nVertexSeg = nElemSeg + 1;
const int nElts = nElemSeg * nElemSeg;
const int nVertex = nVertexSeg * nVertexSeg * 2;

// Number of faces
const int nFacesBorderZ = nElemSeg * nElemSeg;
const int nFacesBorderY = nElemSeg;
const int nFacesBorderX = nElemSeg;
const int nFaces = nFacesBorderZ * 2 + nFacesBorderY * nVertexSeg + nFacesBorderX * nVertexSeg;

double *coords = NULL;         // Element coordinates
int *eltsConnecPointer = NULL; // Connectivity index
int *eltsConnec = NULL;        // Connectivity
int *faceIndex = NULL;         // Face index
int *cell2faceConnec = NULL;   // Cell to face connectivity
int *faceConnecIndex = NULL;   // Face connectivity index
int *faceConnec = NULL;        // Face connectivity

coords = (double *) malloc(sizeof(double) * 3 * nVertex); // 3 coordinates (X, Y, Z) for each vertex
eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1)); // always equal to the number of elements + 1
eltsConnec = (int *) malloc(sizeof(int) * 8 * nElts); // 8 because as we are in 3D the elements are hexahedra
faceIndex = (int *) malloc(sizeof(int) * (nElts + 1)); // always equal to the number of elements + 1
cell2faceConnec = (int *) malloc(sizeof(int) * 6 * nElts); // for tetrahedra it is always 6 times the number of elements
faceConnecIndex = (int *) malloc(sizeof(int) * (nFaces+ 1)); // always equal to the number of faces + 1
faceConnec = (int *) malloc(sizeof(int) * 4 * nFaces); // 4 because faces have 4 vertices

  const double xmin = 0;
  const double xmax = 1;
  const double ymin = 0;
  const double ymax = 1;
  const double zmin = 0;
  const double zmax = 0.1;

int l = 0;
  for (int k = 0; k < 2; ++k){
    for (int j = 0; j < nVertexSeg; ++j){
      for (int i = 0; i < nVertexSeg; ++i){
        coords[l] = (xmax - xmin) / nElemSeg * i;
        std::cout << "Coordinate x: "<< coords[l] << "\n";
        coords[l + 1] = (ymax - ymin) / nElemSeg * j;
        std::cout << "Coordinate y: "<< coords[l + 1] << "\n";
        coords[l + 2] = (zmax - zmin) * k;
        std::cout << "Coordinate z: "<< coords[l + 2] << "\n";
        l = l + 3;
      }
    }
  }

eltsConnecPointer[0] = 0;
std::cout << "eltsConnecPointer: "<< eltsConnecPointer[0] << "\n";
int m = 0;
  for (int i(0); i < nElts; ++i){
    eltsConnecPointer[i + 1] = m + 8;
    std::cout << "eltsConnecPointer: "<< eltsConnecPointer[i + 1] << "\n";
    m = m + 8;
  }

int n = 0;
  for (int j(0); j < nElemSeg; ++j){
    for (int i(0); i < nElemSeg; ++i){
      eltsConnec[n] = j * nVertexSeg + nVertexSeg * nVertexSeg + i;
      std::cout<< "connectivity: "<< eltsConnec[n]<< "\n";
      eltsConnec[n + 1] = j * nVertexSeg + nVertexSeg * nVertexSeg + i + 1;
      std::cout<< "connectivity: "<< eltsConnec[n + 1]<< "\n";
      eltsConnec[n + 2] = j * nVertexSeg + i + 1;
      std::cout<< "connectivity: "<< eltsConnec[n + 2]<< "\n";
      eltsConnec[n + 3] = j * nVertexSeg + i;
      std::cout<< "connectivity: "<< eltsConnec[n + 3]<< "\n";
      eltsConnec[n + 4] = (j + 1) * nVertexSeg + nVertexSeg * nVertexSeg + i;
      std::cout<< "connectivity: "<< eltsConnec[n + 4]<< "\n";
      eltsConnec[n + 5] = (j + 1) * nVertexSeg + nVertexSeg * nVertexSeg + i + 1;
      std::cout<< "connectivity: "<< eltsConnec[n + 5]<< "\n";
      eltsConnec[n + 6] = (j + 1) * nVertexSeg + i + 1;
      std::cout<< "connectivity: "<< eltsConnec[n + 6]<< "\n";
      eltsConnec[n + 7] = (j + 1) * nVertexSeg + i;
      std::cout<< "connectivity: "<< eltsConnec[n + 7]<< "\n";
      n = n + 8;
    }
  }

faceIndex[0] = 0;
std::cout<< "faceIndex: "<< faceIndex[0]<< "\n";
int o = 0;
  for (int i(0); i < nElts; ++i){
    faceIndex[i + 1] = o + 6;
    std::cout<< "faceIndex: "<< faceIndex[i + 1]<< "\n";
    o = o + 6;
  }

int p = 0;
  for (int j(0); j < nElemSeg; ++j){
    for (int i(0); i < nElemSeg; ++i){
    cell2faceConnec[p] = j * nElemSeg + i;
    std::cout<< "faceConnec: "<< cell2faceConnec[p]<< "\n";
    cell2faceConnec[p + 1] = j * nElemSeg + i + nElts;
    std::cout<< "faceConnec: "<< cell2faceConnec[p + 1]<< "\n";
    cell2faceConnec[p + 2] = j * nElemSeg + i + 2 * nElts;
    std::cout<< "faceConnec: "<< cell2faceConnec[p + 2]<< "\n";
    cell2faceConnec[p + 3] = j * nElemSeg + i + 2 * nElts + nElemSeg;
    std::cout<< "faceConnec: "<< cell2faceConnec[p + 3]<< "\n";
    cell2faceConnec[p + 4] = j * nElemSeg + i + 3 * nElts + nElemSeg + j;
    std::cout<< "faceConnec: "<< cell2faceConnec[p + 4]<< "\n";
    cell2faceConnec[p + 5] = j * nElemSeg + i + 3 * nElts + nElemSeg + j + 1;
    std::cout<< "faceConnec: "<< cell2faceConnec[p + 5]<< "\n";
    p = p + 6;
    }
  }

faceConnecIndex[0] = 0;
std::cout << "faceConnecIndex: "<< faceConnecIndex[0] << "\n";
int q = 0;
  for (int i(0); i < nFaces; ++i){
    faceConnecIndex[i + 1] = q + 4;
    std::cout << "faceConnecIndex: "<< faceConnecIndex[i + 1] << "\n";
    q = q + 4;
  }

// connectivity faces and points
int r = 0;
  for (int k(0); k < 2; ++k){
    for (int j(0); j < nElemSeg; ++j){
      for (int i(0); i < nElemSeg; ++i){
        faceConnec[r] = i + j * nVertexSeg + k * nVertexSeg * nVertexSeg;
        std::cout << "faceConnec: "<< faceConnec[r] << "\n";
        faceConnec[r + 1] = i + j * nVertexSeg + k * nVertexSeg * nVertexSeg + 1;
        std::cout << "faceConnec: "<< faceConnec[r + 1] << "\n";
        faceConnec[r + 2] = i + j * nVertexSeg + k * nVertexSeg * nVertexSeg + 1 + nVertexSeg;
        std::cout << "faceConnec: "<< faceConnec[r + 2] << "\n";
        faceConnec[r + 3] = i + j * nVertexSeg + k * nVertexSeg * nVertexSeg + nVertexSeg;
        std::cout << "faceConnec: "<< faceConnec[r + 3] << "\n";
        r = r + 4;
      }
    }
  }

  for (int j(0); j < nVertexSeg; ++j){
    for (int i(0); i < nElemSeg; ++i){
      faceConnec[r] = i + j * nVertexSeg + nVertexSeg * nVertexSeg;
      std::cout << "faceConnec: "<< faceConnec[r] << "\n";
      faceConnec[r + 1] = i + j * nVertexSeg + nVertexSeg * nVertexSeg + 1;
      std::cout << "faceConnec: "<< faceConnec[r + 1] << "\n";
      faceConnec[r + 2] = i + j * nVertexSeg + 1;
      std::cout << "faceConnec: "<< faceConnec[r + 2] << "\n";
      faceConnec[r + 3] = i + j * nVertexSeg;
      std::cout << "faceConnec: "<< faceConnec[r + 3] << "\n";
      r = r + 4;
    }
  }

  for (int j(0); j < nVertexSeg; ++j){
    for (int i(0); i < nElemSeg; ++i){
      faceConnec[r] = i + j * nElemSeg + nVertexSeg * nVertexSeg;
      std::cout << "faceConnec: "<< faceConnec[r] << "\n";
      faceConnec[r + 1] = i + j * nElemSeg + nVertexSeg * nVertexSeg + nVertexSeg;
      std::cout << "faceConnec: "<< faceConnec[r + 1] << "\n";
      faceConnec[r + 2] = i + j * nElemSeg + 3;
      std::cout << "faceConnec: "<< faceConnec[r + 2] << "\n";
      faceConnec[r + 3] = i + j * nElemSeg;
      std::cout << "faceConnec: "<< faceConnec[r + 3] << "\n";
      r = r + 4;     
    }
  }
}