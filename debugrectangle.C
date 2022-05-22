#include <stdlib.h>
#include <iostream>

int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{

int nx = 3, ny = 3;
int nvertex;
int nelts;
int coord_id = 1;

nvertex = nx * ny;
nelts = (nx - 1) * (ny - 1);

const double xmin = -0.2;
const double xmax = 3.2;
const double ymin = -0.2;
const double ymax = 3.2;

double *coords = NULL;
coords = (double*) malloc(sizeof(double) * 3 * nvertex);
int *connecindex = NULL;
connecindex = (int*) malloc(sizeof(int) * (nelts + 1));
int *connec = NULL;
connec = (int*) malloc(sizeof(int) * 4 * nelts);

double *send_values = NULL;
send_values = (double*) malloc(sizeof(double) * nvertex);
double *recv_values = NULL;
recv_values = (double*) malloc(sizeof(double) * nvertex);

// Coordinates table
int l = 0;
  for (int j = 0; j < ny; ++j){
    for (int i = 0; i < nx; ++i){
      coords[l] = xmin + (xmax - xmin) / (nx - 1) * i; // X
      std::cout << "coords: " << coords[l] << "\n"; 
      coords[l + 1] = ymin + (ymax - ymin) / (ny - 1) * j; // Y
      std::cout << "coords: " << coords[l + 1] << "\n"; 
      coords[l + 2] = 0; // Z
      std::cout << "coords: " << coords[l + 2] << "\n"; 
      l = l + 3;
    }
  }

// Connectivity
connecindex[0] = 0;
std::cout << "connecindex: " << connecindex[0] << "\n"; 
  for (int i = 0; i < nelts; ++i){
    connecindex[i + 1] = i * 4 + 4;
    std::cout << "connecindex: " << connecindex[i + 1] << "\n"; 
  }

l = 0;
  for (int j = 0; j < (ny - 1); ++j){
    for (int i = 0; i < (nx - 1); ++i){
      connec[l] = i + j * nx;
      std::cout << "connec: " << connec[l] << "\n"; 
      connec[l + 1] = i + 1 + j * nx;
      std::cout << "connec: " << connec[l + 1] << "\n"; 
      connec[l + 2] = i + 1 + (j + 1) * nx;
      std::cout << "connec: " << connec[l + 2] << "\n"; 
      connec[l + 3] = i + (j + 1) * nx;
      std::cout << "connec: " << connec[l + 3] << "\n"; 
    }
  }

for (int i = 0; i < nvertex; ++i){
  send_values[i] = coords[3 * i + coord_id];
  std::cout << "sendValues: " << send_values[i] << "\n";
  recv_values[i] = 0;
  std::cout << "receivedValues: " << recv_values[i] << "\n";
}

// std::cout << "size of coordinates: " << coords.length() << "\n";
// std::cout << "size of connecindex: " << connecindex.length() << "\n";
// std::cout << "size of connec: " << connec.length() << "\n";
}