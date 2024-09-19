.. _closest_points:

Closest points
==============

C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_closest_points_create

Point clouds definition
"""""""""""""""""""""""

.. doxygenfunction:: PDM_closest_points_n_part_cloud_set

.. doxygenfunction:: PDM_closest_points_src_cloud_set

.. doxygenfunction:: PDM_closest_points_tgt_cloud_set


Computation
"""""""""""

.. doxygenfunction:: PDM_closest_points_compute

.. doxygenfunction:: PDM_closest_points_dump_times


Results
"""""""

.. doxygenfunction:: PDM_closest_points_part_to_part_get

.. doxygenfunction:: PDM_closest_points_get


Finalization
""""""""""""

.. doxygenfunction:: PDM_closest_points_free



Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  .. todo::
    ...

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)


Python API
----------

.. ifconfig:: enable_python_doc == 'ON'


  .. py:class:: ClosestPoints

      Python structure to look for the closest points of a point cloud
      (target) in an other point cloud (source).

    .. rubric:: Initialization

    .. automethod:: Pypdm.Pypdm.ClosestPoints.__init__

    .. rubric:: Point clouds definition

    .. automethod:: Pypdm.Pypdm.ClosestPoints.n_part_cloud_set
    .. automethod:: Pypdm.Pypdm.ClosestPoints.src_cloud_set
    .. automethod:: Pypdm.Pypdm.ClosestPoints.tgt_cloud_set

    .. rubric:: Computation

    .. automethod:: Pypdm.Pypdm.ClosestPoints.compute
    .. automethod:: Pypdm.Pypdm.ClosestPoints.dump_times

    .. rubric:: Results

    .. automethod:: Pypdm.Pypdm.ClosestPoints.part_to_part_get
    .. automethod:: Pypdm.Pypdm.ClosestPoints.points_get

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
