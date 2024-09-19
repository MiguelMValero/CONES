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
 * @file init_Path_Foam.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief OpenFOAM initialization.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
//========== OpenFOAM environment for mesh definition ==========

std::string stringRootPath = std::filesystem::current_path();

Foam::Info << "Le chemin dans CONES est : " << stringRootPath << Foam::endl;

char RunCase[250]= "/";

Foam::Time runTime(Foam::Time::controlDictName, stringRootPath, {RunCase});

// if (cwipiVerbose)
  std::cout << "Runtime created" << std::endl
            << "\n";

Foam::fvMesh mesh(
    Foam::IOobject
    (
      Foam::fvMesh::defaultRegion,
      runTime.timeName(),
      runTime,
      Foam::IOobject::MUST_READ
    )
);