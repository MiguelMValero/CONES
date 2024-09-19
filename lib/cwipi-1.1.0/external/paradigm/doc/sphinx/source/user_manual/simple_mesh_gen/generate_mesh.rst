.. _generate_mesh:

Generation of simple partitioned meshes
=======================================

C API
-----

.. doxygenfunction:: PDM_generate_mesh_rectangle_ngon


Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  .. f:subroutine:: pdm_generate_mesh_rectangle_ngon(comm, elt_type, xmin, ymin, zmin, lengthx, lengthy, n_x, n_y, n_part, part_method, pn_vtx, pn_edge, pn_face, pvtx_coord, pedge_vtx, pface_edge_idx, pface_edge, pface_vtx, pvtx_ln_to_gn, pedge_ln_to_gn, pface_ln_to_gn [, random_factor])

    Create a partitioned rectangular mesh (2D) with descending connectivities

    Admissible values for ``elt_type``:

      - ``PDM_MESH_NODAL_TRIA3``   : triangles
      - ``PDM_MESH_NODAL_QUAD4``   : quadrangles
      - ``PDM_MESH_NODAL_POLY_2D`` : mixed polygons (triangles, quadrangles and octagons)

    :p integer comm [in]:                        MPI communicator
    :p integer elt_type [in]:                    Element type
    :p real xmin [in]:                           Minimal x-coordinate
    :p real ymin [in]:                           Minimal y-coordinate
    :p real zmin [in]:                           Minimal z-coordinate
    :p real lengthx [in]:                        Length of the rectangle in the x-direction
    :p real lengthy [in]:                        Length of the rectangle in the y-direction
    :p integer(pdm_g_num_s) n_x [in]:            Number of points in the x-direction
    :p integer(pdm_g_num_s) n_y [in]:            Number of points in the y-direction
    :p integer n_part [in]:                      Number of partitions
    :p integer part_method [in]:                 Partitioning method
    :p integer(pdm_l_num_s) pn_vtx(:) [out]:     Number of vertices
    :p integer(pdm_l_num_s) pn_edge(:) [out]:    Number of edges
    :p integer(pdm_l_num_s) pn_face(:) [out]:    Number of faces
    :p pdm_pointer_array_t pvtx_coord [out]:     Vertex coordinates
    :p pdm_pointer_array_t pedge_vtx [out]:      Edge->vertex connectivity
    :p pdm_pointer_array_t pface_edge_idx [out]: Index of face->edge connectivity
    :p pdm_pointer_array_t pface_edge [out]:     Face->edge connectivity
    :p pdm_pointer_array_t pface_vtx [out]:      Face->vertex connectivity
    :p pdm_pointer_array_t pvtx_ln_to_gn [out]:  Vertex global ids
    :p pdm_pointer_array_t pedge_ln_to_gn [out]: Edge global ids
    :p pdm_pointer_array_t pface_ln_to_gn [out]: Face global ids
    :p real random_factor [optional]:            Randomization factor (between 0 and 1)

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)


Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. autofunction:: Pypdm.Pypdm.generate_mesh_rectangle_ngon

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)

