.. _part_extension:

Part extension
==============

C API
-----

Enumerators
~~~~~~~~~~~

.. doxygenenum:: PDM_extend_type_t

Initialization
~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_part_extension_create

Set inputs
~~~~~~~~~~

.. doxygenfunction:: PDM_part_extension_connectivity_set

.. doxygenfunction:: PDM_part_extension_vtx_coord_set

.. doxygenfunction:: PDM_part_extension_ln_to_gn_set

.. doxygenfunction:: PDM_part_extension_part_bound_graph_set

.. doxygenfunction:: PDM_part_extension_group_set

.. .. doxygenfunction:: PDM_part_extension_set_part

.. .. doxygenfunction:: PDM_part_extension_part_domain_interface_shared_set

Perform exchange of extended partition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_part_extension_compute

Get outputs
~~~~~~~~~~~

.. doxygenfunction:: PDM_part_extension_connectivity_get

.. doxygenfunction:: PDM_part_extension_ln_to_gn_get

.. doxygenfunction:: PDM_part_extension_vtx_coord_get

.. doxygenfunction:: PDM_part_extension_group_get

.. .. doxygenfunction:: PDM_part_extension_interface_get

.. .. doxygenfunction:: PDM_part_extension_composed_interface_get


Finalize
~~~~~~~~

.. doxygenfunction:: PDM_part_extension_free

Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  ~~~~~~~~~~~~~~

  .. f:autosubroutine:: PDM_part_extension_create

  Set inputs
  ~~~~~~~~~~

  .. f:autosubroutine PDM_part_extension_set_part

  .. f:autosubroutine:: PDM_part_extension_connectivity_set

  .. f:autosubroutine:: PDM_part_extension_vtx_coord_set

  .. f:autosubroutine:: PDM_part_extension_ln_to_gn_set

  .. f:autosubroutine:: PDM_part_extension_part_bound_graph_set

  .. f:autosubroutine:: PDM_part_extension_group_set

  Perform exchange of extended partition
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  .. f:autosubroutine PDM_part_extension_compute

  .. f:subroutine:: pdm_part_extension_compute(part_ext)

    Compute extended partitions

    :p c_ptr part_ext[in]: Part Extension instance

  Get outputs
  ~~~~~~~~~~~

  .. f:autosubroutine:: PDM_part_extension_connectivity_get

  .. f:autosubroutine:: PDM_part_extension_ln_to_gn_get

  .. f:autosubroutine:: PDM_part_extension_vtx_coord_get

  .. f:autosubroutine:: PDM_part_extension_group_get

  Finalize
  ~~~~~~~~

  .. f:autosubroutine PDM_part_extension_free

  .. f:subroutine:: pdm_part_extension_free(part_ext)

    Free a Part Extension structure

    :p c_ptr part_ext[inout]: Part Extension instance

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)

Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. py:class:: PartExtension

    Python structure to perform partition extension.
    Once initialized, all the following
    methods apply to a :class:`PartExtension` instance.

    .. rubric:: Initialization

    .. automethod:: Pypdm.Pypdm.PartExtension.__init__

    .. rubric:: Methods summary

    .. autosummary::
      :nosignatures:

      ~Pypdm.Pypdm.PartExtension.connectivity_set
      ~Pypdm.Pypdm.PartExtension.vtx_coord_set
      ~Pypdm.Pypdm.PartExtension.ln_to_gn_set
      ~Pypdm.Pypdm.PartExtension.part_bound_graph_set
      ~Pypdm.Pypdm.PartExtension.group_set
      ~Pypdm.Pypdm.PartExtension.compute
      ~Pypdm.Pypdm.PartExtension.connectivity_get
      ~Pypdm.Pypdm.PartExtension.vtx_coord_get
      ~Pypdm.Pypdm.PartExtension.ln_to_gn_get
      ~Pypdm.Pypdm.PartExtension.group_get

    .. rubric:: Set inputs

    .. automethod:: Pypdm.Pypdm.PartExtension.connectivity_set
    .. automethod:: Pypdm.Pypdm.PartExtension.vtx_coord_set
    .. automethod:: Pypdm.Pypdm.PartExtension.ln_to_gn_set
    .. automethod:: Pypdm.Pypdm.PartExtension.part_bound_graph_set
    .. automethod:: Pypdm.Pypdm.PartExtension.group_set

    ..  .. autofunction:: Pypdm.Pypdm.PartExtension.set_part
    ..  .. autofunction:: Pypdm.Pypdm.PartExtension.part_domain_interface_shared_set

    .. rubric:: Perform exchange of extended partition

    .. automethod:: Pypdm.Pypdm.PartExtension.compute

    .. rubric:: Get outputs

    .. automethod:: Pypdm.Pypdm.PartExtension.connectivity_get
    .. automethod:: Pypdm.Pypdm.PartExtension.vtx_coord_get
    .. automethod:: Pypdm.Pypdm.PartExtension.ln_to_gn_get
    .. automethod:: Pypdm.Pypdm.PartExtension.group_get

    .. .. autofunction:: Pypdm.Pypdm.PartExtension.get_interface
    .. .. autofunction:: Pypdm.Pypdm.PartExtension.get_composed_interface

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
