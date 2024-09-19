.. _new_cwipi:

New CWIPI
#########

Released in 2023, this version relies on the parallel computational geometry library **ParaDiGM** which is developed by the same team as CWIPI.

General concepts
================

.. include:: concepts.rst


.. _spatial_interp:

Spatial interpolation methods
=============================

.. include:: spatial_interp.rst

Advanced features
=================


User-defined spatial interpolation
----------------------------------

Users can perform customized spatial interpolation based on the geometric mapping computed by CWIPI.
To do so, one needs to provide CWIPI with a pointer to a user-written interpolation function (via ``CWP_Field_interp_function_set`` in C).
Each Field object can have a specific user-defined spatial interpolation.
The following sections show how to retrieve all the information required to write such a function.

Data getters for user interpolation functions in C
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../../pattern/user_interp_function_c_pattern.txt
   :language: c

Data getters for user interpolation functions in Fortran
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../../pattern/user_interp_function_fortran_pattern.txt
   :language: fortran

Data getters for user interpolation functions in Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../../pattern/user_interp_function_python_pattern.txt
   :language: python


Additional exchange features
----------------------------

CWIPI v.1 offers two additional means to exchange data between the coupled codes.
As opposed to fields, these data need not be related to the interface meshes.

.. include:: global_data.rst

.. include:: part_data.rst



C API documentation
===================

.. doxygenfile:: cwp.h
   :project: cwipi

Fortran API documentation
=========================

:download:`Fortran documentation (pdf) <doc_cwp_fortran.pdf>`

.. .. f:automodule:: cwp

Initialization/Finalization
---------------------------

.. f:autosubroutine:: CWP_Init_

.. .. f:autosubroutine:: CWP_Finalize


Coupling
--------

.. f:autosubroutine:: CWP_Cpl_create_

.. f:autosubroutine:: CWP_Cpl_barrier_

.. f:autosubroutine:: CWP_Visu_set_

.. f:autosubroutine:: CWP_User_tgt_pts_set_

.. f:autosubroutine:: CWP_Cpl_del_


Control parameters
------------------

.. .. f:autosubroutine:: CWP_Param_add

.. .. f:autosubroutine:: CWP_Param_set

.. .. f:autosubroutine:: CWP_Param_get

.. f:autosubroutine:: CWP_Param_del_

.. f:autofunction::   CWP_Param_n_get_

.. f:autosubroutine:: CWP_Param_list_get_

.. f:autofunction::   CWP_Param_is_

.. .. f:autosubroutine:: CWP_Param_reduce

.. f:autosubroutine:: CWP_Param_lock_

.. f:autosubroutine:: CWP_Param_unlock_


Interface mesh
--------------

.. f:autosubroutine:: CWP_Mesh_interf_vtx_set_

.. f:autofunction::   CWP_Mesh_interf_block_add_

.. f:autosubroutine:: CWP_Mesh_interf_block_std_set_

.. f:autosubroutine:: CWP_Mesh_interf_f_poly_block_set_

.. f:autosubroutine:: CWP_Mesh_interf_c_poly_block_set_

.. f:autosubroutine:: CWP_Mesh_interf_from_cellface_set_

.. f:autosubroutine:: CWP_Mesh_interf_from_faceedge_set_

.. f:autosubroutine:: CWP_Mesh_interf_finalize_

.. f:autosubroutine:: CWP_Mesh_interf_del_


Spatial interpolation
---------------------

.. f:autosubroutine:: CWP_Spatial_interp_property_set_

.. f:autosubroutine:: CWP_Spatial_interp_weights_compute_

.. f:autosubroutine:: CWP_Computed_tgts_bcast_enable_

.. f:autosubroutine:: CWP_Involved_srcs_bcast_enable_

.. f:autosubroutine:: CWP_N_computed_tgts_get_

.. f:autosubroutine:: CWP_Computed_tgts_get_

.. f:autosubroutine:: CWP_N_uncomputed_tgts_get_

.. f:autosubroutine:: CWP_Uncomputed_tgts_get_

.. f:autosubroutine:: CWP_N_involved_srcs_get_

.. f:autosubroutine:: CWP_Involved_srcs_get_


Fields
------

.. f:autosubroutine:: CWP_Field_create_

.. f:autosubroutine:: CWP_Field_data_set_

.. f:autosubroutine:: CWP_Field_issend_

.. f:autosubroutine:: CWP_Field_wait_issend_

.. f:autosubroutine:: CWP_Field_irecv_

.. f:autosubroutine:: CWP_Field_wait_irecv_

.. f:autosubroutine:: CWP_Field_del_


User-defined spatial interpolation
----------------------------------

.. f:autosubroutine:: CWP_Field_interp_function_set_

.. f:autosubroutine:: CWP_Field_interp_function_unset_

.. f:autosubroutine:: CWP_Field_n_components_get_

.. f:autosubroutine:: CWP_Field_dof_location_get_

.. f:autosubroutine:: CWP_Field_storage_get_

.. f:autosubroutine:: CWP_Field_n_dof_get_

.. f:autosubroutine:: CWP_Field_src_data_properties_get_

.. f:autosubroutine:: CWP_Field_tgt_data_properties_get_

.. f:autosubroutine:: CWP_Cpl_spatial_interp_algo_get_

Location methods
""""""""""""""""

.. f:autosubroutine:: CWP_Field_location_weights_get_

.. f:autosubroutine:: CWP_Field_location_point_data_get_

.. f:autosubroutine:: CWP_Field_location_point_data_get_

Intersection method
"""""""""""""""""""

.. f:autosubroutine:: CWP_Field_intersection_volumes_get_

.. f:autosubroutine:: CWP_Field_intersection_tgt_elt_volumes_get_

Nearest neighbors methods
"""""""""""""""""""""""""

.. f:autosubroutine:: CWP_Field_nearest_neighbors_distances_get_

.. f:autosubroutine:: CWP_Field_nearest_neighbors_coord_get_


Global data
-----------

.. .. f:autosubroutine:: CWP_Global_data_issend_

.. f:autosubroutine:: CWP_Global_data_wait_issend_

.. .. f:autosubroutine:: CWP_Global_data_irecv_

.. f:autosubroutine:: CWP_Global_data_wait_irecv_


Partitioned data
----------------

.. f:autosubroutine:: CWP_Part_data_create_

.. f:autosubroutine:: CWP_Part_data_issend_

.. f:autosubroutine:: CWP_Part_data_wait_issend_

.. f:autosubroutine:: CWP_Part_data_irecv_

.. f:autosubroutine:: CWP_Part_data_wait_irecv_

.. f:autosubroutine:: CWP_Part_data_del_


Python API documentation (:mod:`pycwp`)
=======================================

.. note::
   All enumerators from the :ref:`C API <C API documentation>` are available in the Python API, except that the "`CWP_`" prefix is omitted (eg ``CWP_STATUS_ON`` becomes ``STATUS_ON``).

.. currentmodule:: pycwp

.. automodule:: pycwp
   :imported-members:
   :members:
   :undoc-members:
   :show-inheritance:


