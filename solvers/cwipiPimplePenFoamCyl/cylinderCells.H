Info<< "Setting the coefficients D..." << nl << endl;

const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
const scalarField& V = mesh.V();

//========== We read those cells affected by the force ==========
// The volumeVectorField F is 0 in all cells with the exception of the cells in the
// cylinder and the adjacent cells

//========== CYLINDER CELL SET ==========//

std::string cylinderCells_file = globalRootPath + "/interiorCells.txt";
std::ifstream cylinderCells;
cylinderCells.open(cylinderCells_file);
int cylinderCellsNumber{0};
int *cylinderCellID;

Info<< "The path with my file containing the cell IDs of the cylinder is: " << cylinderCells_file 
<< nl << endl;

if (cylinderCells.is_open())
{
    std::string line;
    int count = 0;
    while( std::getline(cylinderCells, line) )
    {
        std::stringstream ss(line);
        std::string cylinderCellsID_str;
        std::getline(ss, cylinderCellsID_str, ',');

        if (count == 0){
            cylinderCellsNumber = std::stoi(cylinderCellsID_str);
            cylinderCellID = new (std::nothrow) int[cylinderCellsNumber];
        }

        else{
            cylinderCellID[count-1] = std::stoi(cylinderCellsID_str);
        }
        ++count;
	}
}
else
{
   Info<< "Couldn't find cylinderCells.txt" << nl << endl;
}
cylinderCells.close();

Info<< "The number of cells at the cylinder is " << cylinderCellsNumber << nl << endl;
