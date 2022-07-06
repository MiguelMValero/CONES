#include "fvCFD.H"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <cwipi.h>

void define_mesh(double* pointsCoords,
				int* connecIdx,
				int* connec,
				const fvMesh& mesh);