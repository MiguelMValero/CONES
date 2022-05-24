double *coords = NULL;
int *connecindex = NULL;
int *connec = NULL;

const label& nPoints = mesh.nPoints();
const label& nCells = mesh.nCells();
const label& nInternalFaces = mesh.nInternalFaces();
const faceList& Faces = mesh.faces(); // Face to node
int sizeFaces = Faces.size();

// Number of vertices in every direction
int nx = 2;
int ny = nInternalFaces + 1;
int nz = 2;

coords = (double*) malloc(sizeof(double) * 3 * nPoints);
connecindex = (int*) malloc(sizeof(int) * (nCells + 1));
connec = (int*) malloc(sizeof(int) * 4 * nCells);

// Mesh definition
if (rank == 0) printf("        Create mesh\n");

// Coordinates
for (int nbpoint(0); nbpoint < nPoints; ++nbpoint){
double x(Points[nbpoint].x());
double y(Points[nbpoint].y());
double z(Points[nbpoint].z());
}

int variable1 = 0;
for (int nbpoint(0); nbpoint < nPoints; ++nbpoint){
    coords[variable1] = x[nbpoint];
    coords[variable1 + 1] = y[nbpoint];
    coords[variable1 + 2] = z[nbpoint];
    variable1 = variable1 + 3;
}

// Connectivity
connecindex[0] = 0;
  for (int nbcells(0); nbcells < nCells; ++nbcells){
    connecindex[nbcells + 1] = nbcells * 4 + 4;
  }

int variable2 = 0;
  for (int ycells(0); ycells < (ny - 1); ++ycells){
    connec[variable2] = (ny + 1) * 2 + ycells * 2;
    connec[variable2 + 1] = (ny + 1) * 2 + 1 + ycells * 2;
    connec[variable2 + 2] = 1 + ycells * 2;
    connec[variable2 + 3] = ycells * 2;
    connec[variable2 + 4] = (ny + 1) * 2 + ycells * 2 + 2;
    connec[variable2 + 5] = (ny + 1) * 2 + 1 + ycells * 2 + 2;
    connec[variable2 + 6] = 1 + ycells * 2 + 2;
    connec[variable2 + 7] = ycells * 2 + 2;
  }

//******************************************************************** To fill
// Define mesh
cwipi_define_mesh(cl_coupling_name,
                  nvertex,
                  nelts,
                  coords,
                  connecindex,
                  connec);

//******************************************************************** End To fill