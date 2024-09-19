.. _general:

General presentation
####################


Key concepts & terminology
==========================

.. _concept_global_id:

* local/global ids
* partitions
* block
* distribution
* stride
* ownership


Unstructured meshes
===================

* ...

Conventions
===========

Standard elements definition
----------------------------

.. list-table:: Numbering convention used in ParaDiGM for standard elements
  :widths: 50 50

  * - .. figure:: ../../../images/pdm_mesh_nodal_point.svg
        :alt: PDM_MESH_NODAL_POINT

        ``PDM_MESH_NODAL_POINT``

    - .. figure:: ../../../images/pdm_mesh_nodal_bar2.svg
        :alt: PDM_MESH_NODAL_BAR2

        ``PDM_MESH_NODAL_BAR2``


  * - .. figure:: ../../../images/pdm_mesh_nodal_tria3.svg
        :alt: PDM_MESH_NODAL_TRIA3

        ``PDM_MESH_NODAL_TRIA3``


    - .. figure:: ../../../images/pdm_mesh_nodal_quad4.svg
        :alt: PDM_MESH_NODAL_QUAD4

        ``PDM_MESH_NODAL_QUAD4``


  * - .. figure:: ../../../images/pdm_mesh_nodal_tetra4.svg
        :alt: PDM_MESH_NODAL_TETRA4

        ``PDM_MESH_NODAL_TETRA4``


    - .. figure:: ../../../images/pdm_mesh_nodal_pyram5.svg
        :alt: PDM_MESH_NODAL_PYRAM5

        ``PDM_MESH_NODAL_PYRAM5``


  * - .. figure:: ../../../images/pdm_mesh_nodal_prism6.svg
        :alt: PDM_MESH_NODAL_PRISM6

        ``PDM_MESH_NODAL_PRISM6``


    - .. figure:: ../../../images/pdm_mesh_nodal_hexa8.svg
        :alt: PDM_MESH_NODAL_HEXA8

        ``PDM_MESH_NODAL_HEXA8``


|

* Orientation/sign
* Connectivities
* ...


Philosophy/design choices
=========================

* C-contiguous arrays
* coords always 3D
