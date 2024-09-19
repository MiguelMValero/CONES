.. _mesh_location:

Mesh location
=============


C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_mesh_location_create

Source mesh definition
""""""""""""""""""""""

.. doxygenfunction:: PDM_mesh_location_mesh_n_part_set

.. doxygenfunction:: PDM_mesh_location_part_set

.. doxygenfunction:: PDM_mesh_location_nodal_part_set

.. doxygenfunction:: PDM_mesh_location_part_set_2d

.. doxygenfunction:: PDM_mesh_location_nodal_part_set_2d

Target point clouds definition
""""""""""""""""""""""""""""""

.. doxygenfunction:: PDM_mesh_location_n_part_cloud_set

.. doxygenfunction:: PDM_mesh_location_cloud_set

Location computation
""""""""""""""""""""

.. doxygenfunction:: PDM_mesh_location_method_set

.. doxygenfunction:: PDM_mesh_location_tolerance_set

.. doxygenfunction:: PDM_mesh_location_compute

.. doxygenfunction:: PDM_mesh_location_dump_times

Results
"""""""
.. doxygenfunction:: PDM_mesh_location_n_located_get

.. doxygenfunction:: PDM_mesh_location_located_get

.. doxygenfunction:: PDM_mesh_location_n_unlocated_get

.. doxygenfunction:: PDM_mesh_location_unlocated_get

.. doxygenfunction:: PDM_mesh_location_points_in_elt_get

.. doxygenfunction:: PDM_mesh_location_point_location_get

.. doxygenfunction:: PDM_mesh_location_cell_vertex_get

.. doxygenfunction:: PDM_mesh_location_part_to_part_get

Finalization
""""""""""""

.. doxygenfunction:: PDM_mesh_location_free



Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  """"""""""""""

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_create_

  Source mesh definition
  """"""""""""""""""""""

  .. f:subroutine:: pdm_mesh_location_mesh_n_part_set(mloc, n_part)

    Set the number of partitions of the source mesh

    :param c_ptr   mesh_loc [in]: C pointer to PDM_mesh_location_t object
    :param integer n_part   [in]:   Number of partitions


  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_part_set_

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_nodal_part_set

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_part_set_2d_

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_nodal_part_set_2d

  Target point clouds definition
  """"""""""""""""""""""""""""""

  .. f:autosubroutine pdm_mesh_location/pdm_mesh_location_n_part_cloud_set
  .. f:subroutine:: pdm_mesh_location_n_part_cloud_set(mloc, i_point_cloud, n_part)

    Set the number of partitions of a point cloud

    :param c_ptr   mesh_loc      [in]: C pointer to PDM_mesh_location_t object
    :param integer i_point_cloud [in]: Point cloud identifier
    :param integer n_part        [in]: Number of partitions

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_cloud_set_

  Location computation
  """"""""""""""""""""

  .. f:autosubroutine pdm_mesh_location/pdm_mesh_location_method_set
  .. f:subroutine:: pdm_mesh_location_method_set(mesh_loc, method)

    Set the method for computing location (preconditioning stage)

    .. note::
      This is an optional setting

    Admissible values are :
      - ``PDM_MESH_LOCATION_OCTREE``         : Use point octree (default method)
      - ``PDM_MESH_LOCATION_DBBTREE``        : Use bounding-box tree
      - ``PDM_MESH_LOCATION_LOCATE_ALL_TGT`` : All target points are guaranteed to be located

    :p c_ptr mesh_loc[in]: Mesh location instance
    :p integer method[in]: Preconditioning method



  .. f:subroutine:: pdm_mesh_location_tolerance_set(mesh_loc, tol)

    Set the relative tolerance for bounding boxes

    .. note::
      This is an optional setting. By default a relative tolerance equal to 0 is used.

    :p c_ptr mesh_loc[in]: Mesh location instance
    :p real    method[in]: Tolerance



  .. f:autosubroutine pdm_mesh_location/pdm_mesh_location_compute
  .. f:subroutine:: pdm_mesh_location_compute(mesh_loc)

    Compute point location

    :p c_ptr mesh_loc[in]: Mesh location instance


  .. f:subroutine:: pdm_mesh_location_dump_times(mesh_loc)

    Dump elapsed and CPU times

    :p c_ptr mesh_loc[in]: Mesh location instance

  Results
  """""""
  .. f:autosubroutine pdm_mesh_location/pdm_mesh_location_n_located_get
  .. f:function:: pdm_mesh_location_n_located_get(mloc, i_point_cloud, i_part) result(n_located)

    Get the number of located points

    :p c_ptr mesh_loc[in]: Mesh location instance
    :p integer i_point_cloud[in]: Point cloud identifier
    :p integer i_part[in]: Partition identifier
    :p integer n_located[out]: Number of located points

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_located_get_

  .. f:autosubroutine pdm_mesh_location/pdm_mesh_location_n_unlocated_get
  .. f:function:: pdm_mesh_location_n_unlocated_get(mloc, i_point_cloud, i_part) result(n_unlocated)

    Get the number of unlocated points

    :p c_ptr mesh_loc[in]: Mesh location instance
    :p integer i_point_cloud[in]: Point cloud identifier
    :p integer i_part[in]: Partition identifier
    :p integer n_unlocated[out]: Number of unlocated points

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_unlocated_get_

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_points_in_elt_get_

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_point_location_get_

  .. f:autosubroutine pdm_mesh_location/pdm_mesh_location_cell_vertex_get_
  .. f:subroutine:: pdm_mesh_location_cell_vertex_get(mloc, i_part, cell_vtx_idx, cell_vtx)

    Get the cell→vertex connectivity used for internal computations

    .. note::
      For non-standard elements, this connectivity is built by ParaDiGM and is necessary to associate
      the ``points_weights`` array (returned by **pdm_mesh_location_points_in_elt_get**)
      to the appropriate mesh vertices.

    :p c_ptr mesh_loc[in]:           Mesh location instance
    :p integer i_part[in]:           Partition identifier
    :p integer(:) cell_vtx_idx[out]: Index for cell → vertex connectivity
    :p integer(:) cell_vtx[out]:     Cell → vertex connectivity

  .. f:autosubroutine pdm_mesh_location/pdm_mesh_location_part_to_part_get_
  .. f:subroutine:: pdm_mesh_location_part_to_part_get(mesh_loc, icloud, ptp, owner)

    Get Part-to-part instance to exchange data between the source mesh and a target point cloud

    :p c_ptr mesh_loc[in]: Mesh location instance
    :p integer icloud[in]: Point cloud identifier
    :p c_ptr ptp[out]:     Part-to-part instance
    :p integer owner[in]:  Ownership for ``ptp``

  Finalization
  """"""""""""

  .. f:autosubroutine pdm_mesh_location/pdm_mesh_location_free
  .. f:subroutine:: pdm_mesh_location_free(mesh_loc)

    Free a Mesh location structure

    :p c_ptr mesh_loc[inout]: Mesh location instance

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)


Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. py:class:: MeshLocation

    .. we don't use autoclass here, because it does not render good with the autoclass_content = 'both' option.
       If this option is removed, move class description in py file and use autoclass:: Pypdm.Pypdm.MeshLocation

    Python structure to perform mesh location operations. Once initialized, all the following
    methods apply to a :class:`MeshLocation` instance.

    .. rubric:: Initialization

    .. autofunction:: Pypdm.Pypdm.MeshLocation.__init__

    .. rubric:: Instance attributes

    .. autoattribute:: Pypdm.Pypdm.MeshLocation.tolerance
    .. autoattribute:: Pypdm.Pypdm.MeshLocation.method

    .. rubric:: Methods summary

    .. autosummary::
      :nosignatures:

      ~Pypdm.Pypdm.MeshLocation.mesh_n_part_set
      ~Pypdm.Pypdm.MeshLocation.part_set
      ~Pypdm.Pypdm.MeshLocation.nodal_part_set
      ~Pypdm.Pypdm.MeshLocation.part_set_2d
      ~Pypdm.Pypdm.MeshLocation.nodal_part_set_2d
      ~Pypdm.Pypdm.MeshLocation.n_part_cloud_set
      ~Pypdm.Pypdm.MeshLocation.cloud_set
      ~Pypdm.Pypdm.MeshLocation.compute
      ~Pypdm.Pypdm.MeshLocation.dump_times
      ~Pypdm.Pypdm.MeshLocation.located_get
      ~Pypdm.Pypdm.MeshLocation.unlocated_get
      ~Pypdm.Pypdm.MeshLocation.location_get
      ~Pypdm.Pypdm.MeshLocation.points_in_elt_get
      ~Pypdm.Pypdm.MeshLocation.cell_vertex_get
      ~Pypdm.Pypdm.MeshLocation.part_to_part_get


    .. rubric:: Source mesh definition

    .. automethod:: Pypdm.Pypdm.MeshLocation.mesh_n_part_set

    .. automethod:: Pypdm.Pypdm.MeshLocation.part_set

    .. automethod:: Pypdm.Pypdm.MeshLocation.nodal_part_set

    .. automethod:: Pypdm.Pypdm.MeshLocation.part_set_2d

    .. automethod:: Pypdm.Pypdm.MeshLocation.nodal_part_set_2d

    .. rubric:: Target point clouds definition

    .. automethod:: Pypdm.Pypdm.MeshLocation.n_part_cloud_set

    .. automethod:: Pypdm.Pypdm.MeshLocation.cloud_set

    .. rubric:: Location computation

    .. automethod:: Pypdm.Pypdm.MeshLocation.compute

    .. automethod:: Pypdm.Pypdm.MeshLocation.dump_times

    .. rubric:: Results

    .. automethod:: Pypdm.Pypdm.MeshLocation.located_get

    .. automethod:: Pypdm.Pypdm.MeshLocation.unlocated_get

    .. automethod:: Pypdm.Pypdm.MeshLocation.location_get

    .. automethod:: Pypdm.Pypdm.MeshLocation.points_in_elt_get

    .. automethod:: Pypdm.Pypdm.MeshLocation.cell_vertex_get

    .. automethod:: Pypdm.Pypdm.MeshLocation.part_to_part_get

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
