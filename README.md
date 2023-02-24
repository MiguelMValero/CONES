# CONES (Coupling OpenFOAM for Numerical EnvironmentS)

This git repository aims to couple the CFD software OpenFOAM with any other kind of open-source codes. It is currently employed to carry out sequential Data Assimilation techniques and, more specifically, an Ensemble Kalman Filter (EnKF). The communications between the EnKF and OpenFOAM are performed by a coupler called CWIPI.

# Some requirements

In order to use CONES, you need to have OpenFOAM-v8 and cwipi-0.12.0 already compiled and installed. Follow the indications from their official sites. Other versions of these two softwares have not been tested.

**1) OpenFOAM-v8:** https://openfoam.org/download/8-ubuntu/

**2) cwipi-0.12.0:**  https://w3.onera.fr/cwipi/telechargement

**3) Eigen:** (version 3.3.7 or higher)
```
sudo apt install libeigen3-dev
```
**4) Python3:** (version 3.8.10 or higher) https://phoenixnap.com/kb/how-to-install-python-3-ubuntu

It is possible that some specific testcases require additional libraries. Contact the developers in case of further problems.

# Procedure to run a case (e.g. cavity test case)

In the cavity case of OpenFOAM the parameter to be optimized is the velocity at the top wall (initially with a **U**(1 0 0)), and the observations correspond to a simulation where the velocity is equal to **U**(5 0 0) m/s at the top wall. Here we include all the steps to make CWIPI work:

Firstly, clone the directory:
```
git clone https://github.com/MiguelMValero/CWIPI.git
```

Secondly, you need to create the .so executable with all CWIPI functions. In order to do that, go to the folder "cwipiPstream" and compile it by doing (change the path of the headers and libraries of CWIPI inside the "cwipiPstream/Make/options" file in advance):
```
wclean
wmake
```

Then, it is required to compile the solver cwipiIcoFoam. Go to the folder "cwipiIcoFoam" and compile it by doing (change the path of the headers and libraries of "cwipiPstream" inside the "cwipiIcoFoam/Make/options" file in advance):
```
wclean
wmake
```

Next, the EnKF has to be compiled. Go to the folder "cavity_test" and run:
```
make allclean
make all
```

To finish with, run the "cavity_test" case from the OpenFOAM tutorials. Go the "cavity_test" folder and run the bash file by doing:
```
./initCwipiIco
```

NOTE: even though the case in OpenFOAM is launched with 1 processor, it is necessary to specify the option "-parallel"

For any doubts regarding the code you can contact Miguel MARTINEZ VALERO: miguel.martinez_valero@ensam.eu or Lucas VILLANUEVA: lucas.villanueva@ensma.fr
