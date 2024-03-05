// ======== We write a file inside the root directory with the following format
//     time    Lift    Clift    Drag    Cdrag

// ========= DRAG AND LIFT INFORMATION ======== //

scalar dragCoeff = -drag / (0.5*rhoref*sqr(Uref)*zCoord*diameter);
scalar liftCoeff = -lift / (0.5*rhoref*sqr(Uref)*zCoord*diameter);

std::string forceCoeffs = writefCoeffPath + "/forceCoeffsIBM.dat";

if (runTime.value() <= runTime.startTime().value()) {
    const char* forceCoeffsChar = forceCoeffs.c_str();
    remove(forceCoeffsChar);
}

std::ofstream forceCoeffsFile;
forceCoeffsFile.open(forceCoeffs, std::ios::out | std::ios::app);
forceCoeffsFile << runTime.value() << " " << -lift << " " << liftCoeff << " " << -drag << " " << dragCoeff << "\n";
forceCoeffsFile.close();

// ========= Print the lift and drag coefficients =========
Info<< nl << "========== Force coefficients =========" << nl << endl;
Info<< "Cl = " << liftCoeff << endl;
Info<< "Cd = " << dragCoeff << nl << endl;
