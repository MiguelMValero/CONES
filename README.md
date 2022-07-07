# CWIPI - EnKF in the cavity case of OpenFOAM;
In this git repository it is intended to perform the EnKF in the "cavity" case of OpenFOAM. The communications between the EnKF and OpenFOAM are performed by a coupler called CWIPI. Here we include all the steps to make CWIPI work:

Firstly, you need to create the .so executable with all CWIPI functions. In order to do that, go to the folder "cwipiPstream" and compile it by doing:\
**wclean**\
**wmake**

Secondly, it is required to compile the solver cwipiIcoFoam. Go to the folder "cwipiIcoFoam" and compile it by doing:\
**wclean**\
**wmake**

Next, the EnKF has to be compiled. Go to the folder "first_test" and run:\
**make allclean**\
**make all**

To finish with, run the "cavity" case from the OpenFOAM tutorials by doing the following:\
**mpirun -np 1 ./KF_coupling.exe : -np <number_processors> cwipiIcoFoam -parallel** \
NOTE: even though the case in OpenFOAM is launched with 1 processor, it is necessary to specify the option "-parallel"

For any doubts regarding the code you can contact Miguel MARTINEZ VALERO: miguel.martinez_valero@ensam.eu or Lucas VILLANUEVA: lucas.villanueva@ensma.fr
