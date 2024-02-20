Info<< "Calculating the force F..." << nl << endl;

//========== Reading the cells affected by the IBM forcing ==========
// The volumeVectorField F is 0 in all cells with the exception of the cells
// affected by the IBM


//========== Interior Cells of Square ==========//

std::string interiorCells_file = globalCasePath + "/../interiorCells.txt";
std::ifstream interiorCells;
interiorCells.open(interiorCells_file);
int interiorCellsnumber = 1600;
int interiorCellID[interiorCellsnumber];

Info<< "The path with my file containing the cell IDs of the bottom wall is: " << interiorCells_file << nl << endl;
Info<< "The number of cells at the bottom wall is " << interiorCellsnumber << nl << endl;

if (interiorCells.is_open())
{
    std::string line;
    int count = 0;
    while( std::getline(interiorCells, line) )
    {
        std::stringstream ss(line);
        std::string interiorCellID_str;
        std::getline(ss, interiorCellID_str, ',');

        interiorCellID[count] = std::stoi(interiorCellID_str);

        ++count;
	}
}
else
{
   Info<< "Couldn't find interiorCells.txt" << nl << endl;
}
interiorCells.close();

//========== First Layer Cells ==========//

std::string firstCells_file = globalCasePath + "/../firstCells.txt";
std::ifstream firstCells;
firstCells.open(firstCells_file);
int firstCellsnumber = 164;
int firstCellID[firstCellsnumber];

Info<< "Path of the file containing cell IDs of first layer cells: " << firstCells_file << nl << endl;
Info<< "Number of cells in the first layer : " << firstCellsnumber << nl << endl;

if (firstCells.is_open())
{
    std::string line;
    int count = 0;
    while( std::getline(firstCells, line) )
    {
        std::stringstream ss(line);
        std::string firstCellID_str;
        std::getline(ss, firstCellID_str, ',');

        firstCellID[count] = std::stoi(firstCellID_str);

        ++count;
	}
}
else
{
   Info<< "Couldn't find firstCells.txt" << nl << endl;
}
firstCells.close();

Info<< "Before initialization of the force.." << endl;

//========== We initialise the force F with a value equal to 0 everywhere ==========//
volScalarField viscosity = laminarTransport.nu();
forAll(mesh.cells(), cellID)
{ 
    D[cellID][0] = 0;
    D[cellID][1] = 0;
    D[cellID][2] = 0;
    nu[cellID] = viscosity[0];
}

//========= Specifying D for the cells inside the solid region  ==========//
for (int i = 0; i < interiorCellsnumber; ++i){
    
    D[interiorCellID[i]][0] = 1e5;
    D[interiorCellID[i]][1] = 1e5;
    D[interiorCellID[i]][2] = 1e5;
}

//========= Specifying D for the cells in the first layer adjacent to the solid ==========//                
for (int i = 0; i < firstCellsnumber; ++i){
    
    D[firstCellID[i]][0] = 0.0;
    D[firstCellID[i]][1] = 0.0;
    D[firstCellID[i]][2] = 0.0; 

 }

Info<< "Force defined!" << nl << endl;
