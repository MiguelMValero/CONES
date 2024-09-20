# CONES (Coupling OpenFOAM with Numerical EnvironmentS)

This git repository aims to couple the CFD software OpenFOAM with any other kind of open-source codes. It is currently employed to carry out sequential Data Assimilation techniques based on the Ensemble Kalman Filter (EnKF), in particular the application allows for implementing multiple Hyper-Localized EnKF (HLEnKF). The communications between the EnKF and OpenFOAM are performed by a coupler called CWIPI. Documentation can be found in the ```doc``` directory.

Firstly, clone the directory:
```
git clone https://github.com/MiguelMValero/CONES.git
```

## Some requirements

In order to use CONES, you need to have OpenFOAM-v9 (OpenFOAM-v8 if using the branch OpenFOAM8) and cwipi-0.12.0 already compiled and installed. Follow the indications from their official sites. Other versions of these two softwares have not been tested.

**1) OpenFOAM-v9:** https://openfoam.org/download/9-ubuntu/

**2) cwipi-1.1.0:** you can compile and install it from the source files located in /libs/cwipi-1.1.0. To compile CWIPI go to the folder and follow the indications of the corresponding README.md file. You can also follow the instructions that you will find here below (be careful, -DCWP_ENABLE_FORTRAN should only be switched on when employing Ubuntu 22.04 version or higher):
```
mkdir build && cd build
cmake -DCWP_ENABLE_Fortran=ON -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpif90 -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx -DMPI_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=/complete_path/cwipi-1.1.0/build ..
make
make install
```

**3) Eigen:** (version 3.3.7 or higher)
```
sudo apt install libeigen3-dev
```

**4) Python3:** (version 3.8.10 or higher) https://phoenixnap.com/kb/how-to-install-python-3-ubuntu

It is possible that some specific testcases require additional libraries. Contact the developers in case of further problems.


## Procedure to run a case (e.g., cavity test case)

In the cavity case of OpenFOAM the parameter to be optimized is the velocity at the top wall (initially with a **U**(1 0 0)), and the observations correspond to a simulation where the velocity is equal to **U**(5 0 0) m/s at the top wall. Here we include all the steps to make CWIPI work:

To start with, you need to create the .so executable with all CWIPI functions. In order to do that, go to the folder "libs/cwipiPstreamPar" and compile it by doing (change the path of the headers and libraries of CWIPI inside the "libs/cwipiPstreamPar/Make/options" file in advance):
```
wclean
wmake
```

Then, it is required to compile the cwipiHLIcoFoamPar solver. Go to the folder "solvers/HLEnKF/cwipiHLIcoFoamPar" and compile them by doing (change the path of the headers of "cwipiPstreamPar" inside the "solvers/HLEnKF/cwipiHLIcoFoamPar/Make/options" file in advance):
```
wclean
wmake
```

Next, the HLEnKF algorithm has to be compiled. Go to the folder "lib/CONES_interface/HLEnKF" and run:
```
make allclean
make all
```

To finish with, run the "cavity_test" case with the HLEnKF code. Go to "tests/cavity_test" folder and run the bash file by doing:
```
./initCwipiIco
```

NOTE: even though the case in OpenFOAM is launched with 1 processor, it is necessary to specify the option "-parallel"

For any doubts regarding the code you can contact Sarp ER: sarp.er@ensam_eu, Paolo ERRANTE: paolo.errante@ensam.eu, Miguel MARTINEZ VALERO: miguel.martinez_valero@ensam.eu, Tom MOUSSIE: tom.moussie@ensam.eu or Lucas VILLANUEVA: lucas.villanueva@ensma.fr. More information can also be found on the publication "L. Villanueva 2024.pdf".
