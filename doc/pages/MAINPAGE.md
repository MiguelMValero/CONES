# CONES (Coupling OpenFOAM with Numerical EnvironmentS)

\image html header_cones_v2.png width=85%

[TOC]

## Introduction

CONES is an application aiming to couple the CFD software OpenFOAM with any other open-source code. 

Three main elements constitutes CONES:

1. Ensemble of \f$m\f$ simulations in OpenFOAM.
2. Code with the Data Assimilation (DA) algorithm.
3. Open source coupler *CWIPI*.

\image html EnKF_Cones-1.png width=75%

## Data Assimilation

### What is Data Assimilation ?
Data Assimilation is an approach / method for combining **high-fidelity observations** with **low-fidelity model output** to improve the latter.

We want to predict the state of the system and its future in the **best possible way**.

\image html ./img/DA_methods-1.png  width=75%

### Sequential Data Assimilation : Kalman Filter

**Initialization** : initial estimate for \f$\mathbf{x}_0^a\f$

At every \f$ k \f$ DA cycle,

1. **Prediction / Forecast step** : the forecast state \f$ \mathbf{x}_k^f \f$ is given by the forward model \f$\mathcal{M}_{k:k-1}\f$ of our system.

2. **Correction / analysis step** : we perform a correction of the state \f$\mathbf{x}_k^a\f$ based on the forecast state \f$\mathbf{x}_k^f\f$ and some high-fidelity observation \f$\mathbf{y}_k\f$ through the estimation of the so-called *Kalman gain* matrix \f$\mathbf{K}_k\f$ that minimizes the error covariance matrix of the updated state \f$\mathbf{P}_k^a\f$.

\image html ./img/KF_algorithm-1.png  width=50%

The inconveniences of the Kalman Filter are:
1. The Kalman Filter - also the Extended Kalman Filter (EKF) - only works for moderate deviations from linearity and Gaussianity:
\f[ \mathcal{M}_{k:k-1} \longrightarrow \textrm{Navier-Stokes equations (non-linear)} \f]
2. **High dimensionality** of the error covariance matrix \f$\mathbf{P}_k\f$, with a number of degree of freedom equal to \f$3 \times n_\textrm{cells}\f$:
\f[ \left[ \mathbf{P}_k \right]_{3 \times n_\textrm{cells}, \, 3 \times n_\textrm{cells}} \f]

The possible alternatives to the KF are :
1. Particle filter: cannot yet to be applied to very high dimensional systems (\f$m\f$ = size of the ensembles)
2. Ensemble Kalman Filter (EnKF):

    a. High dimensional systems \f$\rightarrow\f$ we avoid the explicit definition of \f$\mathbf{P}_k\f$.

    b. Non-linear models \f$\rightarrow\f$ Monte-Carlo realisations.

    c. But underestimation of \f$\mathbf{P}_k^a\f$ \f$\rightarrow\f$ inflation & localisation.

### Ensemble Kalman Filter with extended state

The system state \f$\mathbf{x}_k\f$ is as follows:
\f[
\begin{bmatrix}
    \mathbf{u}_k \\
    \theta_k
\end{bmatrix}_{3 \times n_\textrm{cells}+n_\theta,m}
\f]
where

- \f$\mathbf{u}_k\f$ is the velocity field of the whole domain at the instant \f$k\f$.
- \f$\theta_k\f$ are the coefficients from the model we want to infer at the instant \f$k\f$.

The observation data \f$\mathbf{y}_k\f$:
\f[
\left[ \mathbf{y}_k \right]_{n_o,m} \rightarrow \mathcal{N}\left( y_k,\mathbf{R}_k \right)
\f]

1. Set of probes with local velocity \f$\mathbf{u}\f$ and pressure \f$p\f$ values.
2. Instantaneous global force coeffcients \f$C_D\f$, \f$C_L\f$, \f$C_f\f$, ...

\image html ./img/EnKF_algorithm-1.png  width=50%

The following different matrices are needed:
1. Observation covariance matrix \f$\mathbf{R}_k\f$ : we assume Gaussian, non-correlated uncertainty for the observation (\f$\sigma_{i,k}\f$ is the standard deviation for the variable \f$i\f$ and the instant \f$k\f$)
\f[
    \left[ \mathbf{R}_k \right]_{n_o,\, n_o} = \sigma^2_{i,k} \left[ I \right]_{n_o, \, n_o}
\f]

2. Samplig matrix \f$\mathcal{H} \left( \mathbf{x}_k^f \right)\f$ : it is the projection of the model into the position of the observations.
\f[
    \left[ \mathcal{H} \left( \mathbf{x}_k^f \right) \right]_{n_o, \, m}
\f]

3. Anomaly matrices for the system's state \f$\mathbf{X}_k^f\f$ and sampling \f$\mathbf{S}_k^f\f$ : they measure the deviation of each realisations with respect to the mean.
\f[
    \left[ \mathbf{X}_k^f \right]_{3 \times n_\textrm{cells}+n_\theta,m} , \qquad \left[ \mathbf{S}_k^f \right]_{n_o,\,m}
\f]

4. Kalman gain \f$\mathbf{K}_k\f$ and covariance localisation \f$\mathbf{L}\f$.
\f[
    \left[ \mathbf{K}_k \right]_{3 \times n_\textrm{cells}+n_\theta,n_o} , \qquad \left[ \mathbf{L} \right]_{3 \times n_\textrm{cells}+n_\theta,n_o} 
\f]



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
