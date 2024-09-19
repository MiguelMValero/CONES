# MultiFidelity Ensemble Kalman Filter (MFEnKF)

[TOC]

## Introduction

The MFEnKF is built using 3 sources of information:

1.  The principal members, denoted as \f$\mathbf{X}_i\f$

2.  The control members, denoted as \f$\mathbf{\hat{U}}_i\f$

3.  The ancillary members, denoted as \f$\mathbf{U}_i\f$

Control members \f$\mathbf{\hat{U}}_i\f$ are a projection of the principal members \f$\mathbf{X}_i\f$ using the projection operator \f$\boldsymbol{\Phi}_r^\star\f$:

\f[ \mathbf{\hat{U}}_i = \boldsymbol{\Phi}_r^\star \left( \mathbf{X}_i \right) \f]

The MFEnKF methodology is sketched below, from [@popov_multifidelity_2021]:

\image html MFEnKF.png width=50%

## Kalman Gain matrix

The Kalman gain \f$\mathbf{\tilde{K}}_i\f$ is computed using the 3 ensemble
types.

### Ensemble means

\f[
\begin{aligned}
    \boldsymbol{\tilde{\mu}}_\mathbf{X} &= N_{\mathbf{X}}^{-1} \mathsf{E}_\mathbf{X} \boldsymbol{1}_{N_\mathbf{X}}, \\
    \boldsymbol{\tilde{\mu}}_\mathbf{\hat{U}} &= N_{\mathbf{X}}^{-1} \mathsf{E}_\mathbf{\hat{U}} \boldsymbol{1}_{N_\mathbf{X}}, \\
    \boldsymbol{\tilde{\mu}}_\mathbf{U} &= N_{\mathbf{U}}^{-1} \mathsf{E}_\mathbf{U} \boldsymbol{1}_{N_\mathbf{U}}, \\
    \nonumber \\
    \boldsymbol{\tilde{\mu}}_\mathbf{\mathcal{H}\left(X\right)} &= N_{o}^{-1} \mathsf{E}_{\mathcal{H}\left( \mathbf{X} \right)} \boldsymbol{1}_{N_o }, \\
    \boldsymbol{\tilde{\mu}}_{\mathcal{H} (\mathbf{\hat{U}})} &= N_{o}^{-1} \mathsf{E}_{\mathcal{H} (\mathbf{\hat{U}})} \boldsymbol{1}_{N_o}, \\
    \boldsymbol{\tilde{\mu}}_{\mathcal{H} \left(\mathbf{U} \right)} &= N_{o}^{-1} \mathsf{E}_{\mathcal{H}(\mathbf{U})} \boldsymbol{1}_{N_o}\end{aligned}
\f]

### Ensemble anomaly matrices

\f[
\begin{aligned}
    \mathsf{A}_{\mathbf{X}} &= \left( \mathsf{E}_\mathbf{X} - \boldsymbol{\tilde{\mu}}_\mathbf{X} \mathbf{1}_{N_\mathbf{X}}^\top \right) \left( N_\mathbf{X} - 1 \right)^{-\frac{1}{2}}, \\
    \mathsf{A}_\mathbf{\hat{U}} &= \left( \mathsf{E}_\mathbf{\hat{U}} - \boldsymbol{\tilde{\mu}}_\mathbf{\hat{U}} \mathbf{1}_{N_\mathbf{X}}^\top \right) \left( N_\mathbf{X} - 1 \right)^{-\frac{1}{2}}, \\
    \mathsf{A}_\mathbf{U} &= \left( \mathsf{E}_\mathbf{U} - \boldsymbol{\tilde{\mu}}_\mathbf{U} \mathbf{1}_{N_\mathbf{U}}^\top \right) \left( N_\mathbf{U} - 1 \right)^{-\frac{1}{2}}, \\
    \nonumber \\
    \mathsf{A}_{\mathcal{H}(\mathbf{X})} &= \left( \mathsf{E}_{\mathcal{H}(\mathbf{X})} - \boldsymbol{\tilde{\mu}}_{\mathcal{H}(\mathbf{X})} \mathbf{1}_{N_o}^\top \right) \left( N_\mathbf{X} - 1 \right)^{-\frac{1}{2}}, \\
    \mathsf{A}_{\mathcal{H}(\mathbf{\hat{U}})} &= \left( \mathsf{E}_{\mathcal{H}(\mathbf{\hat{U}})} - \boldsymbol{\tilde{\mu}}_{\mathcal{H}(\mathbf{\hat{U}})} \mathbf{1}_{N_o}^\top \right) \left( N_\mathbf{X} - 1 \right)^{-\frac{1}{2}}, \\
    \mathsf{A}_{\mathcal{H}(\mathbf{U})} &= \left( \mathsf{E}_{\mathcal{H}(\mathbf{U})} - \boldsymbol{\tilde{\mu}}_{\mathcal{H}(\mathbf{U})} \mathbf{1}_{N_o}^\top \right) \left( N_\mathbf{U} - 1 \right)^{-\frac{1}{2}}\end{aligned}
\f]

### Ensemble covariances

\f[
\begin{aligned}
    \boldsymbol{\tilde{\Sigma}}_\mathbf{X,\mathcal{H}(X)} &= \mathsf{A}_\mathbf{X} \mathsf{A}_{\mathcal{H}(\mathbf{X})}^\top \\
    \boldsymbol{\tilde{\Sigma}}_{\mathbf{\hat{U}},\mathcal{H}(\mathbf{\Phi_r \hat{U}})} &= \mathsf{A}_\mathbf{\hat{U}} \mathsf{A}_{\mathcal{H}(\mathbf{\Phi_r \hat{U}})}^\top \\
    \boldsymbol{\tilde{\Sigma}}_{\mathbf{X},\mathcal{H}(\mathbf{\Phi_r \hat{U}})} &= \mathsf{A}_\mathbf{X} \mathsf{A}_{\mathcal{H}(\mathbf{\Phi_r \hat{U}})}^\top \\
    \boldsymbol{\tilde{\Sigma}}_{\mathbf{\hat{U}},\mathcal{H}(\mathbf{X})} &= \mathsf{A}_\mathbf{\hat{U}} \mathsf{A}_{\mathcal{H}(\mathbf{X})}^\top \\
    \boldsymbol{\tilde{\Sigma}}_{\mathbf{{U}},\mathcal{H}(\mathbf{\Phi_r \mathbf{U}})} &= \mathsf{A}_\mathbf{U} \mathsf{A}_{\mathcal{H}(\Phi_r \mathbf{U})}^\top \\ 
    \nonumber \\    
    \boldsymbol{\tilde{\Sigma}}_\mathbf{\mathcal{H}(X),\mathcal{H}(X)} &= \mathsf{A}_{\mathcal{H}(\mathbf{X})} \mathsf{A}_{\mathcal{H}(\mathbf{X})}^\top \\
    \boldsymbol{\tilde{\Sigma}}_{\mathcal{H}(\mathbf{\Phi_r \hat{U}}),\mathcal{H}(\mathbf{\Phi_r \hat{U}})} &= \mathsf{A}_{\mathcal{H}(\mathbf{\Phi_r \hat{U}})} \mathsf{A}_{\mathcal{H}(\mathbf{\Phi_r \hat{U}})}^\top \\
    \boldsymbol{\tilde{\Sigma}}_{\mathcal{H}(\mathbf{X}),\mathcal{H}(\mathbf{\Phi_r \hat{U}})} &= \mathsf{A}_{\mathcal{H}(\mathbf{X})} \mathsf{A}_{\mathcal{H}(\mathbf{\Phi_r \hat{U}})}^\top \\
    \boldsymbol{\tilde{\Sigma}}_{\mathbf{\mathcal{H}(\Phi_r \hat{U})},\mathcal{H}(\mathbf{X})} &= \mathsf{A}_{\mathcal{H}(\Phi_r \mathbf{\hat{U}})} \mathsf{A}_{\mathcal{H}(\mathbf{X})}^\top \\
    \boldsymbol{\tilde{\Sigma}}_{\mathcal{H}(\Phi_r \mathbf{{U}}),\mathcal{H}(\mathbf{\Phi_r \mathbf{U}})} &= \mathsf{A}_{\mathcal{H}(\Phi_r \mathbf{U})} \mathsf{A}_{\mathcal{H}(\Phi_r \mathbf{U})}^\top \end{aligned}
\f]

### Kalman gain intermediate covariances

\f[
\begin{aligned}
    \boldsymbol{\tilde{\Sigma}_{Z, \mathcal{H}(\mathbf{Z})}} &= \boldsymbol{\tilde{\Sigma}}_\mathbf{X,\mathcal{H}(X)} + \frac{1}{4} \mathbf{\Phi_r} \boldsymbol{\tilde{\Sigma}}_{\mathbf{\hat{U}},\mathcal{H}(\mathbf{\Phi_r \hat{U}})} - \frac{1}{2}  \boldsymbol{\tilde{\Sigma}}_{\mathbf{X},\mathcal{H}(\mathbf{\Phi_r \hat{U}})} - \frac{1}{2} \mathbf{\Phi_r} \boldsymbol{\tilde{\Sigma}}_{\mathbf{\hat{U}},\mathcal{H}(\mathbf{X})} + \frac{1}{4} \mathbf{\Phi_r} \boldsymbol{\tilde{\Sigma}}_{\mathbf{{U}},\mathcal{H}(\mathbf{\Phi_r \mathbf{U}})} \\
    \nonumber \\
    \boldsymbol{\tilde{\Sigma}_{\mathcal{H}(\mathbf{Z}), \mathcal{H}(\mathbf{Z})}} &= \boldsymbol{\tilde{\Sigma}}_\mathbf{\mathcal{H}(X),\mathcal{H}(X)} + \frac{1}{4} \boldsymbol{\tilde{\Sigma}}_{\mathcal{H}(\mathbf{\Phi_r \hat{U}}),\mathcal{H}(\mathbf{\Phi_r \hat{U}})} - \frac{1}{2} \boldsymbol{\tilde{\Sigma}}_{\mathcal{H}(\mathbf{X}),\mathcal{H}(\mathbf{\Phi_r \hat{U}})} - \frac{1}{2} \boldsymbol{\tilde{\Sigma}}_{\mathbf{\mathcal{H}(\Phi_r \hat{U})},\mathcal{H}(\mathbf{X})} + \frac{1}{4} \boldsymbol{\tilde{\Sigma}}_{\mathcal{H}(\Phi_r \mathbf{{U}}),\mathcal{H}(\mathbf{\Phi_r \mathbf{U}})}\end{aligned}
\f]

### Computing Kalman Gain

\f[
\begin{aligned}
    \mathbf{\tilde{K}} = \boldsymbol{\tilde{\Sigma}_{Z, \mathcal{H}(\mathbf{Z})}} \left( \boldsymbol{\tilde{\Sigma}_{\mathcal{H}(\mathbf{Z}), \mathcal{H}(\mathbf{Z})}} + \boldsymbol{\tilde{\Sigma}_{\eta,\eta}} \right)^{-1}\end{aligned}
\f]

## Update and total variate

Each variable can be updated: 
\f[
\begin{aligned}
    \mathbf{X}^a &= \mathbf{X}^f + \mathbf{\tilde{K}} \left( \mathbf{y} - \mathcal{H}(\mathbf{X}) \right) \\
    \mathbf{\hat{U}}^a &= \mathbf{\hat{U}}^f + \mathbf{\Phi_r^\star}\mathbf{\tilde{K}} \left( \mathbf{y} - \mathcal{H}(\mathbf{\hat{U}}) \right) \\
    \mathbf{U}^a &= \mathbf{U}^f +  \mathbf{\Phi_r^\star} \mathbf{\tilde{K}} \left( \mathbf{y} - \mathcal{H}(\mathbf{U}) \right) \\\end{aligned}
\f]

The total variate \f$\mathbf{Z}\f$ represents the combined prior and
posterior knowledge through the linear control variate technique:

\f[
\begin{aligned}
    \mathbf{Z}^f = \mathbf{X}^f - \mathbf{S}\left( \mathbf{\hat{U}}^f - \overline{\mathbf{U}^f} \right) \\
    \mathbf{Z}^a = \mathbf{X}^a - \mathbf{S}\left( \mathbf{\hat{U}}^a - \overline{\mathbf{U}^a} \right)\end{aligned}
\f]

where \f$\mathbf{S} = \mathbf{\Phi_r}/2\f$
