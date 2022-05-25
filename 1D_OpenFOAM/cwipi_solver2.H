int rank;
int coord_id = 1; // value of the coordinate that we want to send (Y)
int stride;
int nNotLocatedPoints = 0;

double *coords = NULL;
int *connecindex = NULL;
int *connec = NULL;
double *sendvalues = NULL;
double *receivedvalues = NULL;

const label& nPoints = mesh.nPoints();
const label& nCells = mesh.nCells();
const label& nInternalFaces = mesh.nInternalFaces();
const faceList& Faces = mesh.faces(); // Face to node
int sizeFaces = Faces.size();
const pointField& Points = mesh.points(); // Node to coordinates

// Number of vertices in every direction
int ny = nInternalFaces + 1;

coords = (double*) malloc(sizeof(double) * 3 * nPoints);
connecindex = (int*) malloc(sizeof(int) * (nCells + 1));
connec = (int*) malloc(sizeof(int) * 4 * nCells);
sendvalues = (double*) malloc(sizeof(double) * nPoints);
receivedvalues = (double*) malloc(sizeof(double) * nPoints);

// Mesh definition
if (rank == 0) printf("        Create mesh\n");

// Coordinates
int variable1 = 0;
for (int nbpoint(0); nbpoint < nPoints; ++nbpoint){
  double x(Points[nbpoint].x());
  double y(Points[nbpoint].y());
  double z(Points[nbpoint].z());

  coords[variable1] = x;
  coords[variable1 + 1] = y;
  coords[variable1 + 2] = z;

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

scalar t = mesh.time().value(); 

//******************************************************************** To fill
// Define mesh
if (t == runTime.deltaTValue()){
    cwipi_define_mesh(cl_coupling_name,
                  nPoints,
                  nCells,
                  coords,
                  connecindex,
                  connec);
}

else {
    cwipi_update_location(cl_coupling_name);
    }

//******************************************************************** End To fill

// Send receive
if (rank == 0) printf("        Exchange\n");

for (int nbpoint(0); nbpoint < nPoints; ++nbpoint){
  sendvalues[nbpoint] = coords[3 * nbpoint + coord_id];
  receivedvalues[nbpoint] = 0;
}
  
stride = 1;

sprintf(cl_exchange_name, "exch1");
sprintf(cl_sending_field_name, "cooy");
sprintf(cl_receiving_field_name, "recv");

// We define the iterations

cwipi_exchange_status_t status = cwipi_exchange(cl_coupling_name,
                                cl_exchange_name,
                                stride,
                                1,     // n_step
                                t,   // physical_time
                                cl_sending_field_name,
                                sendvalues,
                                cl_receiving_field_name,
                                receivedvalues,
                                &nNotLocatedPoints);

//******************************************************************** End To fill