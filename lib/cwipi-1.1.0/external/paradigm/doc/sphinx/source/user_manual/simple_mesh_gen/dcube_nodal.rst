.. _dcube_nodal:

Distributed nodal square/cube mesh
==================================


C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_dcube_nodal_gen_create

Options
"""""""

.. doxygenfunction:: PDM_dcube_nodal_gen_random_factor_set

.. doxygenfunction:: PDM_dcube_nodal_gen_ordering_set


Mesh generation
"""""""""""""""

.. doxygenfunction:: PDM_dcube_nodal_gen_build

.. doxygenfunction:: PDM_dcube_nodal_gen_dmesh_nodal_get


Multi-domain mesh generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_dcube_nodal_cart_topo


Finalization
""""""""""""

.. doxygenfunction:: PDM_dcube_nodal_gen_free




Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  """"""""""""""

  .. f:autosubroutine:: PDM_dcube_nodal_gen_create_

  Options
  """""""

  .. f:autosubroutine PDM_dcube_nodal_gen_random_factor_set

  Mesh generation
  """""""""""""""

  .. f:autosubroutine:: PDM_dcube_nodal_gen_build_

  .. f:autosubroutine:: PDM_dcube_nodal_gen_dmesh_nodal_get_

  Finalization
  """"""""""""

  .. f:autosubroutine PDM_dcube_nodal_gen_free

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)




Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. py:class:: DCubeNodalGenerator

    Python structure to generate a distributed mesh nodal from a cartesian topology.

    .. rubric:: Initialization

    .. autofunction:: Pypdm.Pypdm.DCubeNodalGenerator.__init__

    .. rubric:: Setting options

    .. automethod:: Pypdm.Pypdm.DCubeNodalGenerator.set_random_factor
    .. automethod:: Pypdm.Pypdm.DCubeNodalGenerator.set_ordering

    .. rubric:: Mesh generation

    .. automethod:: Pypdm.Pypdm.DCubeNodalGenerator.compute

    .. automethod:: Pypdm.Pypdm.DCubeNodalGenerator.get_dmesh_nodal
    

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)

