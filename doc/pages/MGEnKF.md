# MultiGrid Ensemble Kalman Filter (MGEnKF)

The MGEnKF is built with \f$n_e\f$ coarse grid ensembles and one fine grid ensemble. The objective is to save computational costs.

Fine members are denoted as \f$\textrm{F}\f$ superscript, \f$\textrm{C}\f$ for coarse members.

1. All ensemble members, on fine and coarse grids, are advanced in time.
\f[
\begin{equation}
    \left( \mathbf{x}^\textrm{F}_k \right)^f = \mathcal{M}^\mathrm{F}_{k:k-1} \left( \mathbf{x}_{k-1}^\textrm{F} \right)
\end{equation}
\f]

\f[
\begin{equation}
    \left( \mathbf{x}^\textrm{C}_k \right)^f = \mathcal{M}^\mathrm{C}_{k:k-1} \left( \mathbf{x}_{k-1}^\textrm{C} \right)
\end{equation}
\f]

2. The fine grid member forecast \f$\left( \mathbf{x}^\textrm{F}_k \right)^f\f$ is projected on the coarse grid using the projection operator \f$\Pi\f$.
\f[
\begin{equation}
    \left( \mathbf{x}^\textrm{C}_k \right)^\star = \Pi_{\textrm{F} \rightarrow \textrm{C}} \left( \left( \mathbf{x}^\textrm{C}_k \right)^f \right)
\end{equation}
\f]

3. The Kalman gain matrix \f$\mathbf{K}_k^\textrm{C}\f$ is computed using the data from the coarse grid members. Coarse ensemble forecasts are updated.

4. The fine grid state projected on the coarse grid is updated. A classical Kalman Filter analysis is realized using the previously obtained Kalman gain.
\f[
\begin{equation}
    \left( \mathbf{x}^\textrm{C}_k \right)^\prime = \left( \mathbf{x}^\textrm{C}_k \right)^\star + \mathbf{K}_k^\textrm{C} \left( \mathbf{y}_k^\textrm{C} - \mathcal{H}_k^\textrm{C} \left( \mathbf{x}_k^\textrm{C} \right)^\star  \right)
\end{equation}
\f]

5. The fine grid updated state \f$\left( \mathbf{x}^\textrm{F}_k \right)^\prime\f$ is determined using the results previously obtained on the coarse grid.
\f[
\begin{equation}
    \left( \mathbf{x}^\textrm{F}_k \right)^\prime = \left( \mathbf{x}^\textrm{F}_k \right)^f + \Pi_{\textrm{C} \rightarrow \textrm{F}}\underbrace{\left( \left( \mathbf{x}^\textrm{C}_k \right)^\prime - \left( \mathbf{x}^\textrm{C}_k \right)^\star \right)}_{= \mathbf{K}_k^\textrm{C} \left( \mathbf{y}_k^\textrm{C} - \mathcal{H}_k^\textrm{C} \left( \mathbf{x}_k^\textrm{C} \right)^\star  \right)}
\end{equation}
\f]
