/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "interpolationCellPointWallModified.H"
#include "cwipiPstream.H"
#include <cwipi.h>

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

static int sendTag = 1;
static int sendTag_params = 3;
static int recvTag = 2;
static int recvTag_params = 4;
static int status;
static int status2;
static int status3;

void addControlParams(int numberCwipiPhase, int cwipiStep, double deltaT, double currentTime, int cwipiMembers, int cwipiObs, int cwipiParams)
{
    cwipi_add_local_int_control_parameter("numberCwipiPhase", numberCwipiPhase);
    cwipi_add_local_int_control_parameter("cwipiStep", cwipiStep);
    cwipi_add_local_double_control_parameter("deltaT", deltaT);
    cwipi_add_local_double_control_parameter("currentTime", currentTime);
    cwipi_add_local_int_control_parameter("cwipiMembers", cwipiMembers);
    cwipi_add_local_int_control_parameter("cwipiObs", cwipiObs);
    cwipi_add_local_int_control_parameter("cwipiParams", cwipiParams);
}

void cwipiCoupling(const fvMesh& mesh, double* pointCoords, int* connecIdx, int* connec)
{
    scalar nCells = mesh.nCells();
    scalar nPoints = mesh.nPoints();
    
    Info << "Here we are 2" << nl << endl;

    // double* pointCoords = new double[3*mesh.nPoints()];
    forAll(mesh.points(),i)
    {
        pointCoords[3*i+0]=mesh.points()[i].x();
        pointCoords[3*i+1]=mesh.points()[i].y();
        pointCoords[3*i+2]=mesh.points()[i].z();
    }
    
    Info << "Here we are 3" << nl << endl;

    // int* connecIdx = new int[mesh.nCells()+1];
    connecIdx[0]=0;
    forAll(mesh.cells(),i)
    {
        connecIdx[i+1]=connecIdx[i]+8;
    }

    Info << "Here we are 4" << nl << endl;

    // int* connec = new int[mesh.nCells()*8];
    forAll(mesh.cells(),i)
    {
        forAll(mesh.cellShapes()[i],j)
        {
            connec[8*i+j]=mesh.cellShapes()[i][j]+1;
        }
    }
    
    Info << "Here we are 5" << nl << endl;

    /*Options :
    CWIPI_COUPLING_SEQUENTIAL
    CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING
    CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING
    */

    cwipi_create_coupling("cwipiFoamCoupling",
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          "FOAM_APE",
                          3,
                          1.0,
                          CWIPI_STATIC_MESH,
                          CWIPI_SOLVER_CELL_CENTER,
                          1,
                          "Ensight Gold",
                          "text");
    Info << "Here we are 8" << nl << endl;

    cwipi_synchronize_control_parameter("FOAM_APE");
    
    cwipi_dump_application_properties();
    cwipi_define_mesh("cwipiFoamCoupling", nPoints, nCells, pointCoords,connecIdx,connec);
    cwipi_locate("cwipiFoamCoupling");
}

void cwipiSend(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration)
{

    Info << "Here we are 9" << nl << endl;
    double* fieldsToSend = new double[3*mesh.nCells()];
    Info << "number of cells: " << mesh.nCells() << endl;
    forAll(mesh.cellCentres(),i)
    {
        fieldsToSend[3*i+0]=vf[i].component(0);
        fieldsToSend[3*i+1]=vf[i].component(1);
        fieldsToSend[3*i+2]=vf[i].component(2);
    }

    Info << "Here we are 10" << nl << endl;
    scalar t = runTime.value();

    Info<< "Before sending to KF" << endl;  
    cwipi_issend("cwipiFoamCoupling","ex1",sendTag,3,cwipiIteration,t,"u0,v0,w0",fieldsToSend,&status);
    cwipi_wait_issend("cwipiFoamCoupling",status);
    Info<< "After sending to KF" << endl;

    switch(status)
    {
        case CWIPI_EXCHANGE_OK :
        Info << "Sent Ok" << endl << "\n";
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        printf("Bad receiving\n");
        break;
        default :
        printf("Error : bad exchange status\n");
    }
    
    delete[] fieldsToSend;
}

void cwipiSendParams(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams)
{

    Info << "Here we are 11" << endl;
    // double* ParamsToSend = new double[mesh.nCells()];
    double* ParamsToSend = new double[cwipiParams];
    // Info << "number of cells: " << mesh.nCells() << endl;

    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (myRank == 4)
    {
        std::cout << "Before sending params to KF" << std::endl; 
        label top = mesh.boundaryMesh().findPatchID("movingWall");
        fvPatchVectorField movingWallU = vf.boundaryField()[top];
        
        ParamsToSend[0]=movingWallU[0].component(0);
 
        MPI_Send(ParamsToSend,1,MPI_DOUBLE,8,sendTag_params,MPI_COMM_WORLD);
        std::cout << "After sending params to KF" << std::endl;
    }

    // std::cout << movingWallU[0].component(0) << "and" << myRank << std::endl << "\n";

    scalar t = runTime.value();

    delete[] ParamsToSend;
}

void cwipiRecv(const fvMesh& mesh, volVectorField& U, const Time& runTime, int cwipiIteration)
{
    double* fieldsToRecv = new double[3 * mesh.nCells()];

    Info << "Before Re-receive" << endl;

    scalar t = runTime.value();
    cwipi_irecv("cwipiFoamCoupling","ex2",recvTag,3,cwipiIteration,t,"u0,v0,w0",fieldsToRecv,&status2);

    cwipi_wait_irecv("cwipiFoamCoupling", status2);

    switch(status2)
    {
        case CWIPI_EXCHANGE_OK :
        Info << "Re-receive Ok" << endl;
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        printf("Bad re-receiving\n");
        break;
        default :
        printf("Error : bad exchange status\n");
    }

    forAll(mesh.cells(), i)
    {
        U[i].component(0) = fieldsToRecv[3*i+0];
        U[i].component(1) = fieldsToRecv[3*i+1];
        U[i].component(2) = fieldsToRecv[3*i+2];
    }
    
    Info << "After Re-receive" << endl << "\n";

    delete[] fieldsToRecv;
}

void cwipiRecvParams(const fvMesh& mesh, volVectorField& U, int cwipiParams)
{
    double* paramsToRecv = new double[cwipiParams];
    
    Info << "Before Re-receive params" << endl;

    MPI_Status status4;
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (myRank > 3)
    {
        MPI_Recv(paramsToRecv,1, MPI_DOUBLE, 8, recvTag_params, MPI_COMM_WORLD, &status4);

        label top = mesh.boundaryMesh().findPatchID("movingWall");
        fvPatchVectorField& movingWallU = U.boundaryFieldRef()[top];

        forAll(movingWallU,faceI)
        {
            Info << movingWallU[faceI].component(0) << endl;

            movingWallU[faceI]=vector(paramsToRecv[0],0,0);
            
            Info << movingWallU[faceI].component(0) << endl;
        }
    }
    
    Info << "After Re-receive params" << endl << "\n";

    delete[] paramsToRecv;
}

void cwipideleteCoupling(double* pointCoords, int* connecIdx, int* connec)
{
    Info << "Delete Cwipi coupling from OF" << endl << "\n";
    cwipi_delete_coupling("cwipiFoamCoupling");
    delete[] pointCoords;
    delete[] connecIdx;
    delete[] connec;
}

void UInterpolation(volVectorField& U, fvMesh& mesh){

const char* proj_file = "/home/miguel/OpenFOAM/miguel-8/run/CWIPI_OpenFOAM/first_test/coord_probes.txt";
const char* UInt = "/home/miguel/OpenFOAM/miguel-8/run/CWIPI_OpenFOAM/first_test/UInt.txt";

std::ifstream file;
// MPI_File file;
// MPI_File_open(MPI_COMM_SELF, proj_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
file.open(proj_file);

// string data = "";
int columns = 0;
int parameters = 3; // 3 coordinates of the probes

interpolationCellPointWallModified<vector> triangulateCellsU(U);
point pointCoord;

// if (Pstream::master()){
// 	if (remove(UInt) != 0)
// 		perror( "Error deleting file with interpolated velocities" );
// 	else
// 		//puts( "File successfully deleted " );
// 		Pout << "File successfully deleted with my processor " << Pstream::master() << nl << endl;
// }
if (remove(UInt) != 0) {
	perror( "Error deleting file with interpolated velocities" );
    // MPI_File_delete(UInt, MPI_INFO_NULL);
}
else
	//puts( "File successfully deleted " );
	Pout << "File successfully deleted with my processor " << Pstream::myProcNo() << nl << endl;
// MPI_File_delete(UInt, MPI_INFO_NULL);

// synchronize the processes
// int dummy;
// Foam::Pstream::scatter(dummy, Foam::Pstream::blocking);

// get the number of processors
// int n = Foam::Pstream::nProcs();

Pout << "Processor number " << Pstream::myProcNo() << " is creating the file with interpolated velocities " <<
nl << endl;

myfile.open(UInt, std::ios::out | std::ios::app);
std::ofstream myfile;
// MPI_File myfile;
// double* buf = (double *)malloc( parameters * sizeof(double) );
// char buf[1000000];
// int count = parameters;
// MPI_Status status5;
// MPI_File_open(MPI_COMM_SELF, UInt,  MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &myfile);

std::ifstream inFile(proj_file);
    if (inFile.is_open())
    {
        std::string line;
        while( std::getline(inFile,line) )
        {
            std::stringstream ss(line);

            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX,',');
			pointCoord.x() = stod(coordX);
            std::getline(ss, coordY,',');
			pointCoord.y() = stod(coordY);
            std::getline(ss, coordZ,','); 
			pointCoord.z() = stod(coordZ);

			label ownCell = mesh.findCell(pointCoord);
			if (ownCell == -1){
				Pout<< "Cell ID not found in my processor" << nl << endl;
				continue;
			}
			vector UatSP(vector::zero);
			UatSP = triangulateCellsU.interpolate(pointCoord, ownCell, -1);

		    /* WE WRITE THE NEW VALUES IN A TXT FILE */
            // buf[0] = UatSP.x();
            // buf[1] = UatSP.y();
            // buf[2] = UatSP.z();
            // buf[3] = pointCoord.x();
            // buf[4] = pointCoord.y();
            // buf[5] = pointCoord.z();
            // MPI_File_write(myfile, buf, count, MPI_DOUBLE_PRECISION, &status5);
			// myfile << UatSP.x() << ",";
			// myfile << UatSP.y() << ",";
			// myfile << UatSP.z() << ",";
			// myfile << pointCoord.x() << ",";
			// myfile << pointCoord.y() << ",";
			// myfile << pointCoord.z() << "\n";

			++columns;
			if (columns == parameters){
				columns = 0;
			}

        }
    }

file.close();
myfile.close();
// MPI_File_close(&file);
// MPI_File_close(&myfile);
}

}

// ************************************************************************* //
