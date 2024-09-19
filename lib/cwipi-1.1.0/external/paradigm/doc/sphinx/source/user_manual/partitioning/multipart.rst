.. _multipart:

Multipart
=========

C API
-----

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

.. .. doxygenfunction:: PDM_multipart_domain_interface_shared_set


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

.. _c_multipart_get_part_mesh_nodal:

.. doxygenfunction:: PDM_multipart_get_part_mesh_nodal

.. doxygenfunction:: PDM_multipart_group_get

.. doxygenfunction:: PDM_multipart_part_ghost_infomation_get

.. doxygenfunction:: PDM_multipart_partition_color_get

.. doxygenfunction:: PDM_multipart_part_hyperplane_color_get

.. doxygenfunction:: PDM_multipart_part_thread_color_get

.. doxygenfunction:: PDM_multipart_part_graph_comm_get


Finalize
~~~~~~~~

.. doxygenfunction:: PDM_multipart_free


Partitioned nodal mesh
~~~~~~~~~~~~~~~~~~~~~~

Here we describe the getters of the structure retrieved using :ref:`PDM_multipart_get_part_mesh_nodal <c_multipart_get_part_mesh_nodal>`.
This allows to have the arrays corresponding to the partitioned mesh described in nodal connectivity style.

.. doxygenfunction:: PDM_part_mesh_nodal_section_n_elt_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_std_get

.. doxygenfunction:: PDM_part_mesh_nodal_vtx_g_num_get


Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  ~~~~~~~~~~~~~~

  .. f:autosubroutine:: PDM_multipart_create_

  Set inputs
  ~~~~~~~~~~

  .. f:autosubroutine PDM_multipart_dmesh_nodal_set
  .. f:subroutine:: pdm_multipart_dmesh_nodal_set(multipart, i_domain, dmesh_nodal)

    Set distributed mesh data for the input domain. The mesh is described by nodal connectivity

    :p c_ptr   multipart  [in]: Multipart instance
    :p integer i_domain   [in]: Domain identifier
    :p c_ptr   dmesh_nodal[in]: Dmesh Nodal instance

  .. f:autosubroutine PDM_multipart_dmesh_set
  .. f:subroutine:: pdm_multipart_dmesh_set(multipart, i_domain, dmesh)

    Set distributed mesh data for the input domain

    :p c_ptr   multipart[in]: Multipart instance
    :p integer i_domain [in]: Domain identifier
    :p c_ptr   dmesh    [in]: Dmesh instance

  .. .. f:autosubroutine:: PDM_multipart_block_set_

  .. PDM_multipart_domain_interface_shared_set

  Renumbering options
  ~~~~~~~~~~~~~~~~~~~

  .. f:autosubroutine:: PDM_multipart_set_reordering_options_

  .. f:autosubroutine:: PDM_multipart_set_reordering_options_vtx_

  Perform partitioning
  ~~~~~~~~~~~~~~~~~~~~

  .. f:autosubroutine PDM_multipart_compute
  .. f:subroutine:: pdm_multipart_compute(multipart)

    Construct the partitioned meshes on all domains

    :p c_ptr multipart[in]: Multipart instance

  Get outputs
  ~~~~~~~~~~~

  .. f:autosubroutine:: PDM_multipart_part_connectivity_get_

  .. f:autosubroutine:: PDM_multipart_part_ln_to_gn_get_

  .. f:autosubroutine:: PDM_multipart_part_vtx_coord_get_

  .. f:autosubroutine:: PDM_multipart_get_part_mesh_nodal_

  .. f:autosubroutine:: PDM_multipart_group_get_

  .. f:autosubroutine:: PDM_multipart_partition_color_get_

  .. f:autosubroutine:: PDM_multipart_part_ghost_infomation_get_

  .. f:autosubroutine:: PDM_multipart_part_graph_comm_get_

  Finalize
  ~~~~~~~~

  .. f:autosubroutine PDM_multipart_free
  .. f:subroutine:: pdm_multipart_free(mpart)

    Free a Multipart structure

    :p c_ptr mpart[inout]: Multipart instance

  Partitioned nodal mesh
  ~~~~~~~~~~~~~~~~~~~~~~

  Here we describe the getters of the structure retrieved using ``PDM_multipart_get_part_mesh_nodal``.
  This allows to have the arrays corresponding to the partitioned mesh described in nodal connectivity style.

  .. f:autosubroutine:: PDM_part_mesh_nodal_section_n_elt_get_

  .. f:autosubroutine:: PDM_part_mesh_nodal_section_std_get_

  .. f:autosubroutine:: PDM_part_mesh_nodal_vtx_g_num_get_

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. py:class:: MultiPart

    Python structure to perform multiple domain mesh partitioning. 
    Once initialized, all the following
    methods apply to a :class:`MultiPart` instance.

    .. rubric:: Initialization

    .. autofunction:: Pypdm.Pypdm.MultiPart.__init__

    .. rubric:: Methods summary

    .. autosummary::
      :nosignatures:

      ~Pypdm.Pypdm.MultiPart.dmesh_nodal_set
      ~Pypdm.Pypdm.MultiPart.dmesh_set
      ~Pypdm.Pypdm.MultiPart.reordering_set
      ~Pypdm.Pypdm.MultiPart.reordering_vtx_set
      ~Pypdm.Pypdm.MultiPart.compute
      ~Pypdm.Pypdm.MultiPart.n_entity_get
      ~Pypdm.Pypdm.MultiPart.connectivity_get
      ~Pypdm.Pypdm.MultiPart.ln_to_gn_get
      ~Pypdm.Pypdm.MultiPart.vtx_coord_get
      ~Pypdm.Pypdm.MultiPart.part_mesh_nodal_get
      ~Pypdm.Pypdm.MultiPart.graph_comm_get
      ~Pypdm.Pypdm.MultiPart.ghost_information_get
      ~Pypdm.Pypdm.MultiPart.color_get
      ~Pypdm.Pypdm.MultiPart.hyper_plane_color_get
      ~Pypdm.Pypdm.MultiPart.thread_color_get

    .. rubric:: Set inputs
    
    .. automethod:: Pypdm.Pypdm.MultiPart.dmesh_nodal_set
    .. automethod:: Pypdm.Pypdm.MultiPart.dmesh_set

    .. rubric:: Renumbering options

    .. todo::
      List available renumbering methods

    .. automethod:: Pypdm.Pypdm.MultiPart.reordering_set
    .. automethod:: Pypdm.Pypdm.MultiPart.reordering_vtx_set

    .. rubric:: Perform partitioning

    .. automethod:: Pypdm.Pypdm.MultiPart.compute

    .. rubric:: Get outputs

    .. automethod:: Pypdm.Pypdm.MultiPart.n_entity_get
    .. automethod:: Pypdm.Pypdm.MultiPart.connectivity_get
    .. automethod:: Pypdm.Pypdm.MultiPart.ln_to_gn_get
    .. automethod:: Pypdm.Pypdm.MultiPart.vtx_coord_get
    .. automethod:: Pypdm.Pypdm.MultiPart.part_mesh_nodal_get
    .. automethod:: Pypdm.Pypdm.MultiPart.graph_comm_get
    .. automethod:: Pypdm.Pypdm.MultiPart.ghost_information_get
    .. automethod:: Pypdm.Pypdm.MultiPart.color_get
    .. automethod:: Pypdm.Pypdm.MultiPart.hyper_plane_color_get
    .. automethod:: Pypdm.Pypdm.MultiPart.thread_color_get


  Partitioned nodal mesh
  ~~~~~~~~~~~~~~~~~~~~~~~

  Here we describe the getters of the structure retrieved using :py:func:`~Pypdm.Pypdm.MultiPart.part_mesh_nodal_get`.
  This allows to have the arrays corresponding to the partitioned mesh described in nodal connectivity style.

  .. .. autoclass:: Pypdm.Pypdm.PMeshNodal

  .. py:class:: PartMeshNodalCaspule

    .. automethod:: Pypdm.Pypdm.PartMeshNodalCaspule.get_sections

  .. .. autofunction:: Pypdm.Pypdm.PartMeshNodalCaspule.vtx_g_num_get


.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
