.. _api_partitioning:

Partitioning
============

.. _api_multipart:

Multipart
---------

.. .. doxygenfile:: pdm_multipart.h
..    :project: paradigm

C API
"""""

Enumerators
~~~~~~~~~~~

.. doxygenenum:: PDM_split_dual_t

.. doxygenenum:: PDM_part_size_t

Initialization
~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_multipart_create


Set inputs
~~~~~~~~~~

.. doxygenfunction:: PDM_multipart_dmesh_nodal_set

.. doxygenfunction:: PDM_multipart_dmesh_set

.. doxygenfunction:: PDM_multipart_block_set


Renumbering options
~~~~~~~~~~~~~~~~~~~

.. todo::
  List available renumbering methods

.. doxygenfunction:: PDM_multipart_set_reordering_options

.. doxygenfunction:: PDM_multipart_set_reordering_options_vtx


Perform partitioning
~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_multipart_compute

.. doxygenfunction:: PDM_multipart_stat_get


Get outputs
~~~~~~~~~~~

.. doxygenfunction:: PDM_multipart_part_n_entity_get

.. doxygenfunction:: PDM_multipart_part_connectivity_get

.. doxygenfunction:: PDM_multipart_part_ln_to_gn_get

.. doxygenfunction:: PDM_multipart_part_vtx_coord_get

.. doxygenfunction:: PDM_multipart_get_part_mesh_nodal

.. doxygenfunction:: PDM_multipart_group_get

.. doxygenfunction:: PDM_multipart_domain_interface_shared_set

.. doxygenfunction:: PDM_multipart_part_ghost_infomation_get

.. doxygenfunction:: PDM_multipart_partition_color_get

.. doxygenfunction:: PDM_multipart_part_hyperplane_color_get

.. doxygenfunction:: PDM_multipart_part_thread_color_get


Finalize
~~~~~~~~

.. doxygenfunction:: PDM_multipart_free



Fortran API
"""""""""""

.. ifconfig:: enable_fortran_doc == 'ON'

  .. todo::
    TODO

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



Python API
""""""""""

.. ifconfig:: enable_python_doc == 'ON'

  Initialization
  ~~~~~~~~~~~~~~

  .. autoclass:: Pypdm.Pypdm.MultiPart


  Set inputs
  ~~~~~~~~~~

  .. autofunction:: Pypdm.Pypdm.MultiPart.dmesh_nodal_set

  .. autofunction:: Pypdm.Pypdm.MultiPart.dmesh_set

  .. autofunction:: Pypdm.Pypdm.MultiPart.joins_set


  Renumbering options
  ~~~~~~~~~~~~~~~~~~~

  .. todo::
    List available renumbering methods

  .. autofunction:: Pypdm.Pypdm.MultiPart.set_reordering

  .. autofunction:: Pypdm.Pypdm.MultiPart.set_reordering_vtx


  Perform partitioning
  ~~~~~~~~~~~~~~~~~~~~

  .. autofunction:: Pypdm.Pypdm.MultiPart.compute


  Get outputs
  ~~~~~~~~~~~

  .. autofunction:: Pypdm.Pypdm.MultiPart.n_entity_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.connectivity_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.ln_to_gn_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.vtx_coord_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.part_mesh_nodal_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.graph_comm_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.ghost_information_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.color_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.hyper_plane_color_get

  .. autofunction:: Pypdm.Pypdm.MultiPart.thread_color_get


.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


.. _api_part_connectivity_transform:


Transformation of partitioned connectivities
--------------------------------------------

.. doxygenfile:: pdm_part_connectivity_transform.h
   :project: paradigm
