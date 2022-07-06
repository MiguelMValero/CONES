<<<<<<< HEAD
# CWIPI - EnKF in the cavity case of OpenFOAM
Here we include all the steps to make CWIPI work:;

Firstly, you need to create the .so executable with all CWIPI functions. In order to do that, go to the folder "cwipiPstream" and compile it by doing:;
wclean;
wmake;

Secondly, it is required to compile the solver cwipiIcoFoam. Go to the folder "cwipiIcoFoam" and compile it by doing:;
wclean;
wmake;

Next, the EnKF has to be compiled. Go to the folder "first_test" and run:;
make allclean;
make all;

To finish with, run the "cavity" case from the OpenFOAM tutorials by doing the following:;
mpirun -np <number_processors> cwipiIcoFoam -parallel : -np 1 ./KF_coupling.exe # NOTE: even though the case in OpenFOAM is launched with 1 processor, it is necessary to specify the option "-parallel";
=======
# CWIPI - EnKF in the cavity case of OpenFOAM;
Here we include all the steps to make CWIPI work:

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
**mpirun -np <number_processors> cwipiIcoFoam -parallel : -np 1 ./KF_coupling.exe**\
NOTE: even though the case in OpenFOAM is launched with 1 processor, it is necessary to specify the option "-parallel"

For any doubts regarding the code you can contact Lucas VILLANUEVA: lucas.villanueva@ensma.fr
>>>>>>> 07e2976861dc457f6d25fa2d7f70d99aa1382789
