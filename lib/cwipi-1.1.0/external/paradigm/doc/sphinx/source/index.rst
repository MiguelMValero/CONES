************************************
Welcome to ParaDiGM's documentation!
************************************

**ParaDiGM** (*Parallel Distributed General Mesh*) is a parallel computational geometry library developed at ONERA.
It provides a progressive framework, which consists of a set of low-, mid- and high-level services usable by developers of scientific computing software.

.. The library is written in C but also has Fortran and Python/Numpy APIs.


###############
Getting started
###############

.. container:: toc-cards

  .. container:: card

    .. figure:: ../../images/index_installation.png
      :target: getting_started/installation.html

    :ref:`Installation <installation>`
      A guide to install and configure ParaDiGM


  .. container:: card

    :ref:`General presentation <general>`
      ParaDiGM's key concepts, terminology, conventions and philosophy


.. toctree::
   :caption: Getting Started
   :maxdepth: 1
   :hidden:

   getting_started/installation
   getting_started/general



|

###########
User manual
###########

.. container:: toc-cards

  .. container:: card

    .. .. figure:: ../../images/icosphere.png
    ..   :target: user_manual/simple_mesh_gen/index.html

    :ref:`Simple mesh generation <simple_mesh_gen>`
      Automatic generation of meshes with simple parametrable shapes


  .. container:: card

    :ref:`Partitioning <partitioning>`
      Graph partitioning, connectivity reconstruction, partition extension and local renumbering


  .. container:: card

    :ref:`Communication graphs <comm_graph>`
      High-level capabilities for building and operating generic communication graphs


  .. container:: card

    :ref:`Global numbering <gnum>`
      Global ID generation


  .. container:: card

    :ref:`Pre-/co-/post-processing <prepro_algo>`
      Distributed geometric and topological algorithms for pre-/co-/post-processing


  .. container:: card

    I/O
      MPI-IO wrappings to read/write files in parallel


  .. container:: card

    High-order meshes
      Utilities for high-order, curved meshes

      .. * ho ordering
      .. * ...


  .. container:: card

    :ref:`FAQ <faq>`
      Frequently asked questions about ParaDiGM

.. * (Reste?)

..   * cellface_orient
..   * geom_elem
..   * triangulate
..   * ...


.. toctree::
   :caption: User manual
   :maxdepth: 1
   :hidden:

   user_manual/simple_mesh_gen/index
   user_manual/partitioning/index
   user_manual/comm_graph/index
   user_manual/gnum/index
   user_manual/prepro_algo/index
   user_manual/faq

|

################
Developer manual
################


.. container:: toc-cards

  .. .. container:: card
  ..
  ..   :ref:`Low level features`
  ..     * octree, dbbtree
  ..     * ...


  .. container:: card

    :ref:`API reference <api>`
      .. toto


  .. container:: card

    :ref:`Coding rules & guidelines<coding_rules>`
      .. toto


.. toctree::
   :caption: Developer manual
   :maxdepth: 1
   :hidden:

   developer_manual/api/index
   developer_manual/coding_rules
   developer_manual/mesh_adaptation/index





.. toctree::
  :maxdepth: 1
  :caption: Appendix
  :hidden:

  changelog
  license


|
