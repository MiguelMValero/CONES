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
 * @file receive_procedure.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief KF receiving procedure.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
for (int j = 1; j < (cwipiMembers + 1); j++)
    {
      sprintf(cl_coupling_name, "cwipiFoamCoupling");
      sprintf(indexChar, "%i", j);
      strcat(cl_coupling_name, indexChar);

      //**** Receive velocity field ****
      if (cwipiVerbose == 1) std::cout << "Before receive state from member " << j << " in EnKF \n";

      cwipi_irecv(cl_coupling_name, "ex1", recvTag, 3, i + 1, time, recv_field_name, values, &status);
      cwipi_wait_irecv(cl_coupling_name, status);

      if (cwipiVerbose == 1) std::cout << "After receive state from member " << j << " in EnKF \n";
      // if (cwipiVerbose == 1) std::cout << "The Length of the Array is : "<< sizeof(values)/sizeof(values[0]) << std::endl; //length

      switch (status)
      {
      case CWIPI_EXCHANGE_OK:
        if (cwipiVerbose)
          std::cout << "Receive state OK from member " << j << " in EnKF \n";
        break;
      case CWIPI_EXCHANGE_BAD_RECEIVING:
        if (cwipiVerbose)
          std::cout << "Bad receiving state from member " << j << " in EnKF \n";
        break;
      default:
        if (cwipiVerbose)
          std::cout << "Error: bad exchange status receiving state from member" << j << " in EnKF \n";
      }

      //**** Receive parameters ****
      if (cwipiVerbose)
        std::cout << "Before receive the parameters from member " << j << " in EnKF \n";
      
      int coming_rank = j*subdomains-(subdomains-1);  //Master rank of each member
      MPI_Recv(paramsValues, cwipiParams, MPI_DOUBLE, coming_rank, recvTag_params, MPI_COMM_WORLD, &status3);

      if (cwipiVerbose){
        std::cout << "After receive the parameters from member " << j << " in EnKF \n";
        std::cout << "==== Parameters well received ==== Parameter 1 = " << paramsValues[0] << std::endl << "\n";
        std::cout << "The number of cells in the side of EnKF is " << nb_cells << std::endl << "\n";
      }

      //======== Configuration as Eigen Matrices =========

      for (int k = 0; k < nb_cells; k++)
      {
        stateMatrix(k, j-1) = values[3*k + 0];
        stateMatrix(k + nb_cells, j-1) = values[3*k + 1];
        stateMatrix(k + 2*nb_cells, j-1) = values[3*k + 2];
      }

      //* The parameters are added at the end of the state vector *
      for (int k = 0; k < cwipiParams; k++)
      {
        stateMatrix(3*nb_cells + k, j-1) = paramsValues[k];
      }
    }