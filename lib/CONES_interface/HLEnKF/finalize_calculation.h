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
 * @file finalize_calculation.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Calculation finalization.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
//=========== Delete coupling after all loops performed ============
if (rank == 0)
printf("Delete Cwipi coupling from c++ file\n");

for (int j = 1; j < cwipiMembers + 1; j++)
{
sprintf(cl_coupling_name, "cwipiFoamCoupling");
sprintf(indexChar, "%i", j);
strcat(cl_coupling_name, indexChar);

cwipi_delete_coupling(cl_coupling_name);
}

if (cwipiVerbose)
std::cout << "After cwipi delete coupling from c++ file" << std::endl;

//========== Freeing memory ==========

delete[] pointCoords;
// delete[] connecIdx;
// delete[] connec;
delete[] face_index;
delete[] cell_to_face_connectivity;
delete[] face_connectivity_index;
delete[] face_connectivity;

free(sendValues);
free(values);
free(paramsValues);
free(paramsSendValues);
// free(configValues);

if (cwipiVerbose)
std::cout << "Arrays deallocated in KF" << std::endl;

//========== Finalize environments ==========
cwipi_finalize();

if (cwipiVerbose)
std::cout << "Cwipi environment finalized from KF " << std::endl;

MPI_Finalize();

if (cwipiVerbose)
std::cout << "MPI environment finalized from KF" << std::endl;