.. _bibliography:

Bibliography
============

*"It's not about aptitude. It's about attitude."* Masayuki Yano

Generating the first mesh
-------------------------

To adapt a mesh, you need a mesh. Here are several suggested tools:

  - vizir4 + EGADS : C. Peyret uses this tool from INRIA-GAMMA but mainly on Mac (Linux and Windows versions are delayed)
  - GMSH : C. Benazet and B. Maugars use the Python API of GMSH (maintaining the group information)

Metric
------

The metric field is a field of (mesh dimension :math:`\times` mesh dimension) matrices defined at each vertex.
This field will allow us to compute the edge length for that given metric and determine whether it is "too long" or "too short".

The metric length in a Riemannian metric space is computed for a segment as :math:`\mathbf{ab} = [\mathbf{a}, \mathbf{b}]` as such :math:`t \in [0,1]`
let :math:`l_{\mathcal{M}}(\mathbf{ab}) = \int_{0}^{1} \sqrt{\mathbf{ab}^{T}\mathcal{M}(\mathbf{a}+t\mathbf{ab})\mathbf{ab}} dt`.

Usual operations on metrics are interpolation (when a new point is inserted) and intersection (several metrics defined at a given point).

Metric interpolation
--------------------

In feflo.a the metric field is interpolated on a background field to avoid losing anisotropy by diffusion.
MMG uses the differential geometry concept of parallel transport to avoid this losing anisotropy but also main direction
deviation.

Metric aligned
--------------

Taking the orthogonal of the metric allows the edges to be aligned with main
direction of the hessian. It is supposed to have the same effect has prism
layers in the boundary layer. According to [1] it is supposed to improve
numerical results but according to J. Vanharen it turns out not to be that useful in practice.

Error Estimate
--------------

On could want the surface triangulation to be an accurate approximation of the smooth boundary of the geometry.
To measure that a geometric error estimates is used. [2] suggested to evaluate the Hausdorff distance between those two surfaces.

To reduce the computational effort and obtain a more accurate result, we want to minimize the difference between
the solution to the PDE (*a priori* estimate) and the finite element approximation (*a posteriori*
estimate).  More details on *a posteriori* error estimates can be found in "Theory and Practice of Finite Elements" in
section 10 on "A Posteriori Error Estimates and Adaptive Meshes".

F. Alauzet provides more details on metric based error estimates aiming to mesure the interpolation error in [3] and [4].
Special attention is given to this topic in section 4 of the PhD thesis of L. Fraza [5].

Scheduling
----------

A. Loseille schedules his operations as:

  - collapse (create a unit mesh)
  - split (create a unit mesh)
  - swap (improve the mesh quality)
  - smooth (improve the mesh quality)

Int his PhD thesis [6], P. Caplan suggests to interleave swaps within
the collapse and split operators to weave out of restrictive (geometry or
visibility-related) topological configurations.

Cavity operator
---------------

Loseille developed a way of abstracting remeshing operations in a unique way : the cavity operator [7].

Swap
----

It is either used when the degree at a vertex is too high (output edges from that vertex) or the element quality too poor.
[8] H. Rakotoarivelo (with F. Ledoux at CEA) decides to focus on lowering the vertex degree since it improves numerical
interpolation (stability and precision) and deals with element quality at the smoothing step.
The following allows to determine the optimal degree : :math:`\min_{T_{h}}R(T_{h}) = \min_{T_{h}} \left\lVert d - d^{*} \right\rVert_{2} = \min_{(P, M, n)} \sum_{i=1}^{n} \left(d[p_{i}] - d^{*}\right)^{\frac{1}{2}}`

According to [8], no method exists to find a global optimum (5-6-7 scheme, puzzle solving). He chooses operate iteratively, by swapping an edge (his tool is 2D only):

  - if the degree of the 4 vertices is in mean value reduced
  - if the quality of the worst element is improved
  - if the deviation of the elements with respect to the tangent planes is no more a given angle

Instead of looping over vertices and handling the linked edges, Rakotoarivelo loops over element pairs (ridges are ignored).

Smoothing
---------

It can either be done by Laplacian smoothing or by optimization.
H. Rakotoarivelo, chooses to do Laplacian smoothing and upon failure solve an optimization problem (minimize the inverse of the element quality).

Projection
----------

Projections upon BREPs created by proprietary CAD software is prone to error since model continuity is only enforced up to a tolerance
often higher than required mesh sizes. :math:`P^{3}` meshes are the first degree for which :math:`G^{1}` continuity at vertices may be enforced.
In [9], two methods to construct a :math:`P^{3}` from a :math:`P^{1}` mesh are studied.
The first is to initialize the Lagrange nodes of each element at the straight position in physical space.
These points are then projected onto the surface using the CAD model.
The second approach evaluates the Lagrange nodes on CAD faces directly.
The last optional step is to apply a Lagrange-to-Bézier transformation since the Bézier representation is a more convenient and generalizable one.
This approach does not guarantee to be :math:`G^{1}` continuous contrarily to the approach chosen in MMG [10], the tangent plane method [11].
Still, [9] consider :math:`G^{1}` continuity not mandatory after benchmark results. Here we talk about :math:`G^{1}` at vertices and not at edges.
Indeed [8], highlights in his PhD thesis that the tangent plane method is compute efficient but it does not guarantee the unicity of tangent planes
of the points on the boundary of the high-order element reconstruction.
To enforce precision and regularity (:math:`G^{1}` continuity at vertices and edges), he uses a quadratic spline patch representation.
To interconnect patches he builds a Gregory patch by blending twists points [12].
This is useful for projections in the smoothing phase when point are not inserted in on a given element but moved around to the optimal position.

Our idea is to use :math:`P^{3}` local reconstruction as a pre-conditioner for the projection direction on a refined background mesh.
This allows to remain accurate to the real world representation if in a given step the mesh is collapsed and later refined again.

Gradation
---------

To ensure gradual evolution of the edges sizes of the mesh, the metric field is smoothed (aka. gradation).
H. Rakotoarivelo proved that the gradation measure choice is arbitrary since they are equivalent (in practice  h-variation is chosen out of simplicity).
From a continuous point of view, the mesh gradation process consists in verifying the uniform continuity of the metric field:

:math:`\forall (x, y) \in \Omega^{2},  \left\lVert M(y) - M(x) \right\rVert \le \left\lVert x - y \right\rVert_{2}`

where :math:`C` is a constant and :math:`\left\lVert . \right\rVert` a matrix norm.
This is an algorithm of quadratic complexity. Alternative less CPU-costly algorithms have been suggested in [13].

Groups
------

The article we relied on to develop the group maintain algorithm is [14].

Graph
-----

An important aspect of parallel mesh adaptation is the scheduling of the remeshing tasks.

Rakotoarivelo [8] uses a graph of tasks from which a maximum stable is extracted. The graph is constructed on mesh entity couples rather than edges.
The data structure induced cache misses when working with edges. An analysis of several approaches revealed that the method of Çatalyurek
proved best. There are no sequential conflit handlings but the are some synchronisation barriers. Morevover, it is an algorithm for
shared-memory parallelism.

Lachat [15] aims to extract independant zones (sub-meshes) on which to call the sequential remesher MMG3D rather than independant tasks.
The zones are not attached to a given partition to avoid marking the output mesh with the partitioning.
His idea is to fuse the graph contraction and seed expansion algorithms. It turns out that a multi-level partitioning
algorithm is more costly than a graph contraction algorithm but less than fusing the above cited algorithms. The
boundary of the extracted zones are not remeshed as is the case with a cavity.

The `SCOTCH_graphColor` algorithm in SCOTCH implements Luby's algorithm. In ParaDiGM we used the Greedy algorithm
for the cavity-cavity graph coloring. The development of a parallel propagation algorithm would be the solution
to the graph coloring problem (as well as creation of a Voronoï diagram or mesh generation).

Tools
=====

The Unstructured Grid Adaptation Working Group is an open gathering of researchers working on adapting simplicial meshes to conform to a metric field.
They have created benchmarks available here: https://github.com/UGAWG.

In this section, we established a list of known mesh adaptation tools. Feel free to add other ones.
Let's start with tool that mention a form of parallelism:

  - https://github.com/hobywan/trinity (C++, 2D surface shared-memory)
  - https://github.com/MmgTools (C, 2D and 3D, ParMMG MPI-partitioned)
  - https://github.com/sandialabs/omega_h/tree/main (C++, optionally MPI, OpenMP, CUDA, 2D, 3D)
  - https://github.com/AMReX-Codes/amrex (massively parallel but block-structured, C++)
  - CDT3D (parallel)
  - https://github.com/nasa/refine (C, MPI)
  - https://gitlab.inria.fr/PaMPA/PaMPA (C using MMG3D and PT-SCOTCH)

Let's move on to other tools:

  - https://github.com/tucanos/tucanos (Airbus, Rust, 2D and 3D)
  - Yams by P. Frey
  - feflo.a by A. Loseille
  - EPIC (Boeing)
  - Pragmatic (Imperial College London)
  - https://github.com/hpc-maths/samurai (C++)
  - http://www.p4est.org/
  - https://optimad.github.io/PABLO/

ParaDiGM's approach
===================

Why do we want a parallel mesh adaptation tool?
-----------------------------------------------

To answer this question, we focus on the analysis by G. Puigt of the CREATE compressor (available at Ecole Centrale Lyon in the LMFA).
The compressor is composed of 4 stator rows (an inlet guide vane and 3 stators) and 3 rotors.
Using periodicity to simplify the geometry is not possible if one wants to simulate a rotating stall.
Simulating coarsely (50000 nodes per blade) the full machine leads to a final mesh of 29.6 million points for the 592 blades.
A simulation on a single core with this coarse grid might not be possible due to insufficient memory.
The number of control volumes proposed by GAMMA in their latest paper on turbomachinery leads to 11 millions cells for a Rotor 37 blade.
With a simple extrapolation that means 6512 billion cells for the CREATE compressor. Taking 1 million grid cells per CPU means we need 7000 CPUs.

Let us underline that the usual approach for mesh adaptation in parallel is by working with on partitions and to
refine the partition boundaries to remove those fake ridges as shown in [7]. They considered themselves this not to be an optimal solution.

What points do we want to work on?
----------------------------------

- implement swap (adapt seed for face-swap) (local and global?)
- implement smoothing (local when staring and global?)
- provide a score to prioritize cavities
- extract area in which mesh adaptation will be done
- cavity prison to avoid blocking cavity (limit cavity growth) or interleave swap operations like P. Caplan suggests
- check background mesh is coherent with the volume mesh
- multi-section background mesh
- Hilbert partitioning
- try not all groups on all procs
- change metric interpolation (with a background mesh for instance)
- change quality tolerance
- asynchronism
- check independent to parallelism
- unitary tests
- projection using :math:`P^3` reconstruction for direction (local or global?)
- does MMG3D do projections ?
- propagation algorithm for graph coloring

References
==========

[1] A. Loseille. “Recent Improvements on Cavity-Based Operators for RANS Mesh Adaptation”. In: (2018).

[2] G. Balarac F. Basile P. Benard F. Bordeu J.-B. Chapelier L. Cirrottola G. Caumon C.Dapogny P. Frey A. Froehly G. Ghigliotti
R. Laraufie G. Lartigue C. Legentil R. Mercier V. Moureau C. Nardoni S. Pertant M. Zakari. “Tetrahedral Remeshing in the Context
of Large-Scale Numerical Simulation and High Performance Computing”. In: MathematicS In Action (2022).

[3] P. Frey F. Alauzet. “Estimateur d’erreur géométrique et métriques anisotropes pour l’adaptation de maillage”. In: (2003).

[4] F. Alauzet. “Metric-Based Anisotropic Mesh Adaptation”. In: (2010).

[5] L. Fraza. “3D anisotropic mesh adaptation for Reynolds Averaged Navier-Stokes simulations”. In: (2020).

[6] P. Caplan. “Four-Dimensional Anisotropic Mesh Adaptation for Spacetime Numerical Simulations”. In: (2012).

[7] V. Menier A. Loseille F. Alauzet. “Unique cavity-based operator and hierarchical domain partitioning for fast parallel generation of
anisotropic meshes”. In: Computer-Aided Design (2017).

[8] H. Rakotoarivelo. “Contribution au co-design de noyaux irréguliers
sur accélérateurs manycore: application au remaillage anisotrope
pour le calcul numérique intensif”. In: (2018).

[9] L. Rochery A. Loseille. “P3 Bézier CAD surrogates for anisotropic
mesh adaptation”. In: Computer-Aided Design (2023).

[10] P. Frey C. Dapogny C. Dobrzynski. “Three-dimensional adaptive
domain remeshing, implicit domain meshing, and applications to
free and moving boundary problems”. In: Journal of Computational
Physics (2014).

[11] Vlachos A. Peters J. Boyd C. Mitchell J.L. “Curved PN triangle”.
In: Proceedings of the 2001 Symposium on Interactive 3D Graphics (2001).

[12] Walton and Meek. “A triangular G1 patch from boundary curves”. In: Computer-Aided Design (1996).

[13] F. Alauzet. “Size gradation control of anisotropic meshes”. In: Finite Elements in Analysis and Design (2009).

[14] D. Marcum A. Loseille R. Löhner. “Robust Boundary Layer Mesh Generation”. In: ().

[15] C. Lachat. "Conception et validation d’algorithmes de remaillage parallèles à mémoire distribuée
basés sur un remailleur séquentiel". In: (2013).
