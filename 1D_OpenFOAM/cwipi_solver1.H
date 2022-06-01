int grank;
int commWorldSize;
int coord_id;
float geom_tol = 0.1;

// MPI Initilization
MPI_Comm localcomm;
char *codeName;
char *codeCoupledName;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &grank);
MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

if (grank < commWorldSize / 2) {
  codeName = "code1";
  coord_id = 0; // coordinate to send (X)
  codeCoupledName = "code2";
}
else {
  codeName = "code2";
  coord_id = 1; // coordinate to send (Y)
  codeCoupledName = "code1";
}

//******************************************************************** To fill
// Initialization of the coupling
cwipi_init(MPI_COMM_WORLD,
          codeName,
          &localcomm);

//******************************************************************** End To fill

MPI_Comm_rank(localcomm, &rank);

//******************************************************************** To fill
// Create coupling
if (rank == 0) printf("        Create coupling\n");
cwipi_solver_type_t solver_type;
solver_type = CWIPI_SOLVER_CELL_VERTEX;
sprintf(output_format,"EnSight Gold");
sprintf(output_format_option,"text");
cwipi_create_coupling(cl_coupling_name,
                      CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                      codeCoupledName,                           // Coupled application id
                      3,                                         // Geometric entities dimension
                      geom_tol,                                  // Geometric tolerance
                      CWIPI_STATIC_MESH,                         // Mesh type
                      solver_type,                               // Solver type
                      1,                                         // Postprocessing frequency
                      output_format,                             // Postprocessing format
                      output_format_option);

//******************************************************************** End To fill

const label& nPoints = mesh.nPoints();
const label& nCells = mesh.nCells();
const label& nInternalFaces = mesh.nInternalFaces();
const pointField& Points = mesh.points(); // Node to coordinates

// Number of vertices in every direction
int ny = nInternalFaces + 1;

coords = (double*) malloc(sizeof(double) * 3 * nPoints);
connecindex = (int*) malloc(sizeof(int) * (nCells + 1));
connec = (int*) malloc(sizeof(int) * 4 * nCells);

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
for (int nbcells; nbcells < nCells; ++nbcells){
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