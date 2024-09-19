# CONES (Coupling OpenFOAM with Numerical EnvironmentS)

Documentation can be found in the ```doc``` directory.

## Procedure to run a case (e.g., cavity test case)

In the cavity case of OpenFOAM the parameter to be optimized is the velocity at the top wall (initially with a **U**(1 0 0)), and the observations correspond to a simulation where the velocity is equal to **U**(5 0 0) m/s at the top wall. Here we include all the steps to make CWIPI work:

To start with, you need to create the .so executable with all CWIPI functions. In order to do that, go to the folder "libs/cwipiPstream" and compile it by doing (change the path of the headers and libraries of CWIPI inside the "libs/cwipiPstream/Make/options" file in advance):
```
wclean
wmake
```

Then, it is required to compile the different solvers (cwipiIcoFoam, cwipiControlIcoFoamPar and cwipiAncillaryIcoFoamPar). Go to the folder "solvers/cwipiIcoFoam" and so on, and compile them by doing (change the path of the headers of "cwipiPstream" inside the "solver/cwipiIcoFoam/Make/options" file in advance):
```
wclean
wmake
```

Next, the MFEnKF has to be compiled. Go to the folder "lib/CONES_interface" and run:
```
make allclean
make all
```

To finish with, run the "cavity_test" case with the MFEnKF. Go to "tests/cavity_test" folder and run the bash file by doing:
```
./initCwipiIco
```

NOTE: even though the case in OpenFOAM is launched with 1 processor, it is necessary to specify the option "-parallel"

For any doubts regarding the code you can contact Miguel MARTINEZ VALERO: miguel.martinez_valero@ensam.eu or Lucas VILLANUEVA: lucas.villanueva@ensma.fr. More information can also be found on the publication "L. Villanueva 2024.pdf".
