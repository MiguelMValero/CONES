
if(Pstream::nProcs() > 1){
  Info << "This application is bugged in parallel !" <<endl;
  return 1;
}

Info << "Reading DA settings" << endl;

IOdictionary ensembleDict
(
  IOobject
  (
     "ensembleDict",
     runTime.system(),
     mesh,
     IOobject::NO_READ,
     IOobject::AUTO_WRITE
  )
);



IOdictionary daDict
(
  IOobject
  (
    "daDict",
    runTime.system(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
  )
);


// Read number of Ensemble Members
scalar nEns (readScalar(daDict.lookup("nEns")));
Info << "No. of Ensemble Members = " << nEns << endl;

// Read the state 
dictionary state = daDict.subDict("state");
word stateType = state.get<word>("type");
word patchName = state.get<word>("patchName");
word fieldName = state.get<word>("fieldName");
word fieldType = state.get<word>("fieldType");
word sdType = state.get<word>("sdType");
vector sdValue(state.lookup("sdValue"));



if( stateType == "BC" ){
  Info << "State type: "  << stateType <<endl;
  Info << "Patch name: " << patchName <<endl;
  Info << "Field name: " << fieldName <<endl;
  Info << "Field type: " << fieldType << endl;
  label patchID = mesh.boundaryMesh().findPatchID(patchName);

  // Check if the patch indicated exists
  if(patchID==-1) {
    Pout << "Patch not found: " << patchName << " ->the program will crash" <<  endl;
    return 1;
  }
  
  if(fieldName == "U"){
    if(fieldType == "fixedValue"){
      
      // Check which proc got the boundaryPatch
      vector Ubc = {0, 0, 0};
      bool goodProc = false;
      if(Pstream::myProcNo() == 0){
        Foam::system("rm -r Ensemble*");
        if(mesh.boundaryMesh()[patchID].nFaces() > 0){
          Ubc = U.boundaryField()[patchID][0];
          Pout << "Master has the value" <<endl;
        }else{
          for(int i = 1; i < Pstream::nProcs() ; i++){
            Pout << "Master receive from proc " << i << endl;
            IPstream fromSlave(Pstream::commsTypes::blocking, i);
            fromSlave >> goodProc;
            if(goodProc){
              fromSlave >> Ubc;
              Pout << "received " << Ubc << endl;
            }
          }
        }
      }else{
        if(mesh.boundaryMesh()[patchID].nFaces() > 0){
          Pout << "I'm proc "<< Pstream::myProcNo() <<" with " <<  mesh.boundaryMesh()[patchID].nFaces() << endl;
          Ubc = U.boundaryField()[patchID][0];
          goodProc = true;
          OPstream toMaster(Pstream::commsTypes::blocking, 0);
          toMaster << goodProc;
          if(goodProc){
            toMaster << Ubc;
            Pout << "sent " << Ubc << " to master" << endl;
          }
        }else{
          OPstream toMaster(Pstream::commsTypes::blocking, 0);
          toMaster << goodProc;
        }
      }
      

      volVectorField Ubound
      (
        IOobject
        (
          fieldName,
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh
      );


     
      if(Pstream::myProcNo() == 0){
        Pout << "Ubc = " << Ubc << endl;
        #include "Random.H"
        Random rnd = Random(clock::getTime());
        Info << "Variance " << sdValue << endl;
        scalar rndNorm;
        scalar rndVal;
        scalar rndSigm;
        List<vector> xList(nEns, vector::zero);
        for(int i = 0 ; i < nEns ; i++){
          for(int j = 0; j < 3; j++){
            rndNorm = rnd.GaussNormal<scalar>();
            rndSigm = Ubc[j]*sdValue[j]/100;
            rndVal = rndNorm*rndSigm + Ubc[j];
            xList[i][j] = rndVal;
          }
        }



        
        ensembleDict.add("Ensembles", xList, true);
        ensembleDict.Foam::regIOobject::write();
        List<vector> xState(ensembleDict.lookup("Ensembles")); 
        string clone_case = "foamCloneCase . Ensemble";
        for(int i = 0; i < nEns; i++){
          Pout << "Writing boundary " << xState[i] <<endl;
          forAll(Ubound.boundaryFieldRef()[patchID], iFace){
            Ubound.boundaryFieldRef()[patchID][iFace] = xState[i];
          }

          //Ubound.boundaryFieldRef()[patchID] = xState[i];
          Pout << "Ubound boundary " << Ubound.boundaryFieldRef()[patchID] << endl;
          Ubound.Foam::regIOobject::write();
          clone_case.append(std::to_string(i));
          Foam::system("echo ---");
          Foam::system(clone_case);
          clone_case = "foamCloneCase . Ensemble";
       }
        /*
        forAll(Ubound.boundaryFieldRef()[patchID], iFace){
          Ubound.boundaryFieldRef()[patchID][iFace] = Ubc;
        }
        Ubound.Foam::regIOobject::write();
        */


        /*
        //string decompose_case = "decomposePar -case Ensemble";
        string decompose_case = "cd Ensemble";
        for(int i = 0; i < nEns ; i++){
          decompose_case.append(std::to_string(i));
          decompose_case.append(" decomposePar ; cd ../");
          Foam::system(decompose_case);
          decompose_case = "cd Ensemble";
        }
        */
        Pout << "xstate = " << xState << endl;
      }
      //Pout << "Ubound " << Ubound.boundaryField() << endl;

      }else{
        Info << "Field not found: " << fieldName << " ->the program will crash" << endl;
        return 1;
      }
    }
  }



List<vector> xObs ( daDict.lookup("xObs"));
List<vector> Uobs ( daDict.lookup("Uobs"));
labelList cellObs(5);

int obsCount = 0;
forAll(xObs, x)
{
  if(mesh.findCell(xObs[x]) > 0){
    cellObs[obsCount] = mesh.findCell(xObs[x]);
    //Pout << "Closest" << nl << cellObs[obsCount] << endl;
    obsCount++;
  }
}
Pout << "Got " << obsCount << " observations " << endl;
int rankObs = -1;
if(obsCount > 0){
  rankObs = Pstream::myProcNo();
  forAll(cellObs, x)
  {
    Pout << "x = " << xObs[x] <<  "\t Uobs = " << Uobs[x] << "\t Cell ID " << cellObs[x] << endl;
  }
}
