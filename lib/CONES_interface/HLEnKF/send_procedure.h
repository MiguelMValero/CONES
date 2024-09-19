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
 * @file send_procedure.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief KF sending procedure.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */

for (int j = 1; j < cwipiMembers + 1; j++)
    {
      if (j == 1) tic();
      prepare_sendBack(sendValues, paramsSendValues, values, stateMatrixUpt, nb_cells, cwipiParams, j, cwipiVerbose);
      sprintf(cl_coupling_name, "cwipiFoamCoupling");

      sprintf(indexChar, "%i", j);
      strcat(cl_coupling_name, indexChar);

      if (cwipiVerbose){
        std::cout << "Before re-sent the parameters from EnKF to member " << j << "\n";
        std::cout << "The name of the coupling is " << cl_coupling_name << "\n";
      }

      cwipi_issend(cl_coupling_name, "ex2", sendTag, 3, i+1, time, recv_field_name, sendValues, &status2);
      cwipi_wait_issend(cl_coupling_name, status2);

      switch (status2)
      {
        case CWIPI_EXCHANGE_OK:
          if (cwipiVerbose)
            std::cout << "Re-sent Ok from EnKF to member " << j << "\n";
          break;
        case CWIPI_EXCHANGE_BAD_RECEIVING:
          if (cwipiVerbose)
            std::cout << "Bad re-sent from EnKF to member " << j << "\n";
          break;
        default:
          if (cwipiVerbose)
            std::cout << "Error: bad re-sent status from EnKF to member " << j << "\n";
      }
      // if (cwipiVerbose)
      //   printf("After re-send \n");

      int going_rank = j*subdomains-(subdomains-1);
      MPI_Send(paramsSendValues, cwipiParams, MPI_DOUBLE, going_rank, sendTag_params, MPI_COMM_WORLD);
      if (j == 1) toc("send back procedure of one member");
    }