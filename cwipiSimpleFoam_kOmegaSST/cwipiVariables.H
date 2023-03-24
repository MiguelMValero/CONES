Switch cwipiSwitch = runTime.controlDict().lookup("cwipiSwitch");

std::ifstream cFile ("cwipiConfig");

int values[17]={0};
int k = 0;

if (cFile.is_open())
{
    std::string line;
    while(std::getline(cFile, line))
    {
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace),line.end());
        if( line.empty() || line[0] == '#' )
        {
            continue;
        }
        auto delimiterPos = line.find("=");
        // auto name = line.substr(0, delimiterPos);
        auto value = line.substr(delimiterPos + 1);

        values[k] = stod(value);
        // std::cout << values[k] << '\n';
        k = k+1;
    }
}
else 
{
    std::cerr << "Couldn't open config file for reading.\n";
}

int cwipiStep = values[0];
int cwipiMembers = values[1];
int cwipiObsU = values[2];
int cwipiObsp = values[3];
int cwipiParams = values[4];
double geom_tol = values[5];
int cwipiVerbose = values[12];
int cwipiParamsObs = values[14];

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

double* pointCoords = new double[3*mesh.nPoints()];
//int* connecIdx = new int[mesh.nCells()+1];
//int* connec = new int[mesh.nCells()*8];

int* face_index = new int[mesh.nCells()+1];
int* face_connectivity_index = new int[mesh.nFaces()+1];

int c2fconnec_size = 0;
forAll(mesh.cells(),i){
    c2fconnec_size = c2fconnec_size + mesh.cells()[i].size();
}
int* cell_to_face_connectivity = new int[c2fconnec_size];

int fconnec_size = 0;
forAll(mesh.faces(),i){
    fconnec_size = fconnec_size + mesh.faces()[i].size();
}
int* face_connectivity = new int[fconnec_size];

if (cwipiVerbose == 1) Info<< "All the parameters from the configuration file have already been defined" << nl << endl;