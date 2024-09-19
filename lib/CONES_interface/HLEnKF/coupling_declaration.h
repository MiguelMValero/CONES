/*--------------------------------*- C++ -*----------------------------------*\
|                         __________  _   _____________                       |
|                        / ____/ __ \/ | / / ____/ ___/                       |
|                       / /   / / / /  |/ / __/  \__ \                        |
|                      / /___/ /_/ / /|  / /___ ___/ /                        |
|                      \____/\____/_/ |_/_____//____/                         |
|                                                                             |
| C.  : Coupling               |    C.O.N.Es. version : 2405                  |
| O.  : OpenFOAM (with)        |    OpenFOAM version  : OpenFOAM-org v9       |
| N.  : Numerical              |    Lucas Villanueva   /   Miguel M.Valero    |
| Es. : EnvironmentS           |    Tom Moussie        /   Sarp Er            |
|      ===============         |    Paolo Errante      /   Marcello Meldi     |
\*---------------------------------------------------------------------------*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**
 * @file coupling_declaration.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Coupling declarations.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
//========== Declaration of variables ==========

  char cl_coupling_name[20];
  char output_format[20], output_format_option[20];
  char codeName[20] = {"KF"};
  char codeCoupledName[20];
  char recv_field_name[20] = {"U,V,W"};
  char indexChar[20];

  static int recvTag = 10;
  static int recvTag_params = 11;

  static int status;
  static int sendTag = 12;
  static int sendTag_params = 13;
  
  static int status2;
  MPI_Status status3;

  int grank;
  int rank;
  int commWorldSize;
  // int nNotLocatedPoints = 0;

  //========== MPI Initilization ==========

  MPI_Comm localcomm;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &grank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  //========== CWIPI Initialization ==========

  cwipi_init(MPI_COMM_WORLD, codeName, &localcomm);
  MPI_Comm_rank(localcomm, &rank);

  //== Retrieve mesh information

  int nb_cells = mesh.nCells();
  if (cwipiVerbose) std::cout << "The number of cells from EnKF is " << mesh.nCells() << std::endl << "\n";

  double *values = (double *)calloc(sizeof(double), 3 * nb_cells);
  double *sendValues = (double *)calloc(sizeof(double), 3 * nb_cells);
  double *paramsValues = (double *)calloc(sizeof(double), cwipiParams);
  double *paramsSendValues = (double *)calloc(sizeof(double), cwipiParams);

  double *pointCoords = new double[3 * mesh.nPoints()];
  int *face_index = new int[mesh.nCells() + 1];
  int *face_connectivity_index = new int[mesh.nFaces() + 1];

  int c2fconnec_size = 0;
  forAll(mesh.cells(), i)
  {
    c2fconnec_size = c2fconnec_size + mesh.cells()[i].size();
  }
  int *cell_to_face_connectivity = new int[c2fconnec_size];

  int fconnec_size = 0;
  forAll(mesh.faces(), i)
  {
    fconnec_size = fconnec_size + mesh.faces()[i].size();
  }
  int *face_connectivity = new int[fconnec_size];

  //========== Create coupling ==========

  if (cwipiVerbose)
    printf("Create coupling\n");
  cwipi_solver_type_t solver_type;
  solver_type = CWIPI_SOLVER_CELL_CENTER;
  sprintf(output_format, "EnSight Gold");
  sprintf(output_format_option, "text");

  //* Coupling declaration in a loop in order to have a coupling with each OF instance *
  for (int j = 1; j < (cwipiMembers + 1); j++)
  {
    sprintf(cl_coupling_name, "cwipiFoamCoupling");
    sprintf(codeCoupledName, "cwipiFoam");
    sprintf(indexChar, "%i", j);
    strcat(cl_coupling_name, indexChar);
    strcat(codeCoupledName, indexChar);

    if (cwipiVerbose) std::cout<< "From the KF side the name of my coupling is " << 
    cl_coupling_name << std::endl;

    cwipi_create_coupling(cl_coupling_name,
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          codeCoupledName,   // Coupled application id
                          3,                 // Geometric entities dimension
                          geometricTolerance,          // Geometric tolerance
                          CWIPI_STATIC_MESH, // Mesh type
                          solver_type,       // Solver type
                          0,                 // Postprocessing frequency
                          output_format,     // Postprocessing format
                          output_format_option);

    cwipi_synchronize_control_parameter(codeCoupledName);

    if (cwipiVerbose) std::cout<< "Coupling called " << cl_coupling_name << " correctly created " << std::endl;

    //========== Mesh definition ============

    if (cwipiVerbose) printf("%s: KF : Create mesh \n", __func__);

    define_mesh(pointCoords, face_index, cell_to_face_connectivity, face_connectivity_index, face_connectivity, c2fconnec_size, fconnec_size, mesh, cl_coupling_name, cwipiVerbose);
  }
