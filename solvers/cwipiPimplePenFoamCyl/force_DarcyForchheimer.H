//========== We initialise the coefficients D with a value equal to 0 everywhere ==========//
volScalarField viscosity = laminarTransport.nu();
forAll(mesh.cells(), cellID)
{ 
    D[cellID][0] = 0;
    D[cellID][1] = 0;
    D[cellID][2] = 0;
}

//========= Uniform value of D at the cylinder ==========//
Info<< "Defining cells with coefficient different from zero" << endl;
std::vector<int> forces;
for (int i = 0; i < cylinderCellsNumber; ++i){
    
    D[cylinderCellID[i]][0] = penaltyCoeff*viscosity[0];
    D[cylinderCellID[i]][1] = penaltyCoeff*viscosity[0];

    forces.push_back(cylinderCellID[i]);
}

////========= Loop over boundary patches ==========//
//forAll (mesh.boundary(), patch)
//{
//    const word& patchName = mesh.boundary()[patch].name();
//    if (patchName == "baffleCylinderFaces_master" || patchName == "baffleCylinderFaces_slave"){
//        forAll(mesh.boundary()[patch], facei)
//        {
//            label bCell = boundaryMesh[patch].faceCells()[facei];
//            D[bCell][0] = penaltyCoeff*viscosity[0];
//            D[bCell][1] = penaltyCoeff*viscosity[0];
//
//            if (std::find(forces.begin(), forces.end(), bCell) == forces.end()){
//                // If not present add it to the vector
//                forces.push_back(bCell);
//            }
//        }
//    }
//}

Info<< "Coefficients defined with a value equal to " << D[cylinderCellID[0]][0] << endl;
