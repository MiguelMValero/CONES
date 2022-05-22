#include <stdlib.h>
#include <iostream>

int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{

int nvertex = 11;
int nelts = 5;
int coord_id = 0;

double *coords = NULL;
coords = (double*) malloc(sizeof(double) * 3 * nvertex);
int *connecindex = NULL;
connecindex = (int*) malloc(sizeof(int) * (nelts + 1));
int *connec = NULL;
connec = (int*) malloc(sizeof(int) * 21);

double *values = NULL;
values = (double*) malloc(sizeof(double) * nvertex);
double *localvalues = NULL;
localvalues = (double*) malloc(sizeof(double) * nvertex);

coords[0] = 0, coords[1] = 0, coords[2] = 0;
coords[3] = 1, coords[4] = 0, coords[5] = 0;
coords[6] = 2, coords[7] = 0, coords[8] = 0;
coords[9] = 3, coords[10] = 0, coords[11] = 0;
coords[12] = 0, coords[13] = 1, coords[14] = 0;
coords[15] = 2, coords[16] = 1, coords[17] = 0;
coords[18] = 3, coords[19] = 1, coords[20] = 0;
coords[21] = 1, coords[22] = 2, coords[23] = 0;
coords[24] = 0, coords[25] = 3, coords[26] = 0;
coords[27] = 2, coords[28] = 3, coords[29] = 0;
coords[30] = 3, coords[31] = 3, coords[32] = 0;

connecindex[0] = 0;
connecindex[1] = 3;
connecindex[2] = 7;
connecindex[3] = 11;
connecindex[4] = 16;
connecindex[5] = 21;

connec[0] = 1, connec[1] = 2, connec[2] = 5;
connec[3] = 3, connec[4] = 4, connec[5] = 7, connec[6] = 6;
connec[7] = 5, connec[8] = 8, connec[9] = 10, connec[10] = 9;
connec[11] = 5, connec[12] = 2, connec[13] = 3, connec[14] = 6, connec[15] = 8;
connec[16] = 6, connec[17] = 7, connec[18] = 11, connec[19] = 10, connec[20] = 8;

for (int i = 0; i < nvertex; ++i){
  values[i] = coords[3 * i + coord_id];
  std::cout << "values to be sent: " << values[i] << "\n";
}

}