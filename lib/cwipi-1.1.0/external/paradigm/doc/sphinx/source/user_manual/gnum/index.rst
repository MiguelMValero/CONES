.. _gnum:


################
Global numbering
################


.. _gen_gnum:

Global numbering generation
===========================

**ParaDiGM** relies heavily on the notion of :ref:`global numbering <concept_global_id>`.
If your code does not use global ids, these can be generated from geometric data.
Such global numbering is achieved by encoding Cartesian coordinates along the `Morton space-filling curve <https://en.wikipedia.org/wiki/Z-order_curve>`_.

.. Alternatively, from parents & nuplets...
.. Either way, the first step consists in creating an instance of ``PDM_gen_gnum_t`` (or :class:`Pypdm.Pypdm.GlobalNumbering` in Python).


C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_gnum_create

Set inputs
""""""""""

.. doxygenfunction:: PDM_gnum_set_from_coords

.. doxygenfunction:: PDM_gnum_set_from_parents

.. doxygenfunction:: PDM_gnum_set_parents_nuplet

Build global numbering
""""""""""""""""""""""

.. doxygenfunction:: PDM_gnum_compute

Get outputs
"""""""""""

.. doxygenfunction:: PDM_gnum_get

Finalization
""""""""""""

.. doxygenfunction:: PDM_gnum_free


Example
"""""""
The following example shows how to build a global numbering from a set of geometric coordinates (extract from the test case ``pdm_t_gen_gnum.c``).

.. code:: c

  #include "pdm.h"
  #include "pdm_gnum.h"

  // First, create a PDM_gen_gnum_t instance and set some parameters
  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,     // dimension
                                             n_part,
                                             merge,
                                             1.e-3, // tolerance
                                             PDM_MPI_COMM_WORLD,
                                             PDM_OWNERSHIP_USER);

  // Then, provide the coordinates array for each partition
  // (`char_length` can be NULL if `merge` is disabled)
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_coords(gen_gnum,
                             i_part,
                             n_elts[i_part],
                             coords[i_part],
                             char_length);
  }

  // Once all partitions have been set, build the global numbering
  PDM_gnum_compute(gen_gnum);

  // Finally, retrieve the computed global id arrays
  PDM_g_num_t **gnum = malloc(sizeof(PDM_g_num_t) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    gnum[i_part] = PDM_gnum_get(gen_gnum,
                                i_part);

  }

  // Deallocate the PDM_gen_gnum_t instance
  PDM_gnum_free(gen_gnum);



.. The ``dim``, ``merge`` and ``tolerance`` arguments are only relevant if you want the global numbering to be based on geometric data.


Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  .. .. f:automodule:: pdm_gnum

  Initialization
  """"""""""""""

  .. f:autosubroutine:: pdm_gnum_create_

  Set inputs
  """"""""""

  .. f:autosubroutine:: pdm_gnum_set_from_coords_

  .. f:autosubroutine:: pdm_gnum_set_from_parents_

  .. .. f:autosubroutine:: pdm_gnum_set_parents_nuplet_

  Build global numbering
  """"""""""""""""""""""

  .. f:subroutine:: pdm_gnum_compute(gen_gnum)

    Build global numbering

    :p c_ptr gen_gnum [in]:  C pointer to PDM_gen_gnum_t object

  Get outputs
  """""""""""

  .. f:autosubroutine:: pdm_gnum_get_

  Finalization
  """"""""""""

  .. f:subroutine:: pdm_gnum_free(gen_gnum)

    Free the Global Numbering Generation object

    :p c_ptr gen_gnum [in]:  C pointer to PDM_gen_gnum_t object

  Example
  """""""
  The following example shows how to build a global numbering from a set of geometric coordinates.

  .. code:: fortran

    use pdm
    use pdm_gnum
    use iso_c_binding

    type(cptr)                    :: gen_gnum
    double precision,     pointer :: coord => null()
    integer(pdm_g_num_s), pointer :: gnum  => null()
    integer                       :: i_part

    ! First, create a PDM_gen_gnum_t instance and set some parameters
    call pdm_gnum_create(gen_gnum,           &
                         dim,                &
                         n_part,             &
                         merge,              &
                         tolerance,          &
                         MPI_COMM_WORLD,     &
                         PDM_OWNERSHIP_USER)

    ! Then, provide the coordinates array for each partition
    ! (`char_length` can be null() if `merge` is disabled)
    do i_part = 1, n_part
      ! get coordinates pointer for current partition
      coords = my_data_structure(i_part)%coords

      call pdm_gnum_set_from_coords(gen_gnum,       &
                                    i_part,         &
                                    n_elts(i_part), &
                                    coords,         &
                                    null())
    enddo

    ! Once all partitions have been set, build the global numbering
    call pdm_gnum_compute(gen_gnum)

    ! Finally, retrieve the computed global id arrays
    do i_part = 1, n_part
      call pdm_gnum_get(gen_gnum, &
                        i_part,   &
                        gnum)
    enddo

    ! Deallocate gen_gnum
    call pdm_gnum_free(gen_gnum)

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  API reference
  """""""""""""

  .. py:class:: GlobalNumbering

    Python structure to create a global numbering.
    Once initialized, all the following
    methods apply to a :class:`GlobalNumbering` instance.

    .. rubric:: Initialization

    .. automethod:: Pypdm.Pypdm.GlobalNumbering.__init__

    .. rubric:: Set inputs

    .. automethod:: Pypdm.Pypdm.GlobalNumbering.set_from_coords
    .. automethod:: Pypdm.Pypdm.GlobalNumbering.set_from_parent
    .. automethod:: Pypdm.Pypdm.GlobalNumbering.set_parents_nuplet

    .. rubric:: Build global numbering

    .. automethod:: Pypdm.Pypdm.GlobalNumbering.compute

    .. rubric:: Get outputs

    .. automethod:: Pypdm.Pypdm.GlobalNumbering.get

  Example
  """""""
  The following example shows how to build a global numbering from a set of geometric coordinates (extract from the test case ``pdm_t_gnum_p.py``).


  .. literalinclude:: ../../../../../test/pdm_t_gnum_p.py
    :name: python_gen_gnum_ex
    :language: python
    :dedent: 2
    :lines: 6,7,62-83

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)

