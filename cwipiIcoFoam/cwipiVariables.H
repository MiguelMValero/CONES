Switch cwipiSwitch = runTime.controlDict().lookup("cwipiSwitch");
int cwipiMembers = readScalar(runTime.controlDict().lookup("cwipiMembers"));
int cwipiObs = readScalar(runTime.controlDict().lookup("cwipiObs"));
int cwipiParams = readScalar(runTime.controlDict().lookup("cwipiParams"));
scalar cwipiStep = readScalar(runTime.controlDict().lookup("cwipiStep"));

IOdictionary decomposeParDict
(
    IOobject
    (
        "decomposeParDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dictionary simpleCoeffs = decomposeParDict.subDict("simpleCoeffs");

int nbParts = readScalar(decomposeParDict.lookup("numberOfSubdomains"));
List<scalar> partsRepart =  simpleCoeffs.lookup("n");

scalar cwipiTimestep = 1;
scalar numberCwipiPhaseDouble = (runTime.endTime().value()/runTime.deltaTValue())/cwipiStep;
int numberCwipiPhase = floor(numberCwipiPhaseDouble);
int cwipiIteration = 1;
scalar cwipiPhaseCheck = 0;