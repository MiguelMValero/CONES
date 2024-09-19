.. _ptp:

Part to Part
============

C API
-----

Enumerators
"""""""""""

.. doxygenenum:: PDM_mpi_comm_kind_t

.. doxygenenum:: PDM_part_to_part_data_def_t

Initialization
""""""""""""""

.. doxygenfunction:: PDM_part_to_part_create

.. doxygenfunction:: PDM_part_to_part_create_from_num2_triplet


Information on Part 2 side
""""""""""""""""""""""""""

.. doxygenfunction:: PDM_part_to_part_ref_lnum2_get

.. doxygenfunction:: PDM_part_to_part_unref_lnum2_get

.. doxygenfunction:: PDM_part_to_part_gnum1_come_from_get


Exchange
""""""""

.. doxygenfunction:: PDM_part_to_part_iexch

.. doxygenfunction:: PDM_part_to_part_iexch_wait

.. doxygenfunction:: PDM_part_to_part_reverse_iexch

.. doxygenfunction:: PDM_part_to_part_reverse_iexch_wait

.. .. doxygenfunction:: PDM_part_to_part_issend

.. .. doxygenfunction:: PDM_part_to_part_issend_wait

.. .. doxygenfunction:: PDM_part_to_part_reverse_issend

.. .. doxygenfunction:: PDM_part_to_part_reverse_issend_wait

.. .. doxygenfunction:: PDM_part_to_part_irecv

.. .. doxygenfunction:: PDM_part_to_part_irecv_wait

.. .. doxygenfunction:: PDM_part_to_part_reverse_irecv

.. .. doxygenfunction:: PDM_part_to_part_reverse_irecv_wait


Finalization
""""""""""""

.. doxygenfunction:: PDM_part_to_part_free



Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  """"""""""""""

  .. f:autosubroutine:: PDM_part_to_part_create

  Information on Part 2 side
  """"""""""""""""""""""""""

  .. f:autosubroutine:: PDM_part_to_part_ref_lnum2_get

  .. f:autosubroutine:: PDM_part_to_part_unref_lnum2_get

  .. f:autosubroutine:: PDM_part_to_part_gnum1_come_from_get

  Exchange
  """"""""

  .. f:autosubroutine:: PDM_part_to_part_iexch

  .. f:autosubroutine PDM_part_to_part_iexch_wait

  .. f:subroutine:: pdm_part_to_part_iexch_wait(ptp, request)

    Finalize a non-blocking exchange (Part1→Part2)

    :p c_ptr ptp [in]:       Part-to-Part instance
    :p integer request [in]: Request

  .. f:autosubroutine:: PDM_part_to_part_reverse_iexch

  .. f:autosubroutine PDM_part_to_part_reverse_iexch_wait

  .. f:subroutine:: pdm_part_to_part_reverse_iexch_wait(ptp, request)

    Finalize a non-blocking exchange (Part2→Part1)

    :p c_ptr ptp [in]:       Part-to-Part instance
    :p integer request [in]: Request

  .. f:autosubroutine PDM_part_to_part_issend

  .. f:autosubroutine PDM_part_to_part_issend_wait

  .. f:autosubroutine PDM_part_to_part_irecv_raw

  .. f:autosubroutine PDM_part_to_part_irecv_wait_raw

  Finalization
  """"""""""""

  .. f:autosubroutine PDM_part_to_part_free

  .. f:subroutine:: pdm_part_to_part_free(ptp)

    Free a Part-to-Part structure

    :p c_ptr ptp [inout]: Part-to-part instance


.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. py:class:: PartToPart

    Python structure to perform partition-to-partition data exchanges. Once initialized, all the following
    methods apply to a :class:`PartToPart` instance.

    .. rubric:: Initialization

    .. autoclass:: Pypdm.Pypdm.PartToPart.__init__


    .. rubric:: Methods summary

    .. autosummary::
      :nosignatures:

      ~Pypdm.Pypdm.PartToPart.get_referenced_lnum2
      ~Pypdm.Pypdm.PartToPart.get_unreferenced_lnum2
      ~Pypdm.Pypdm.PartToPart.get_gnum1_come_from
      ~Pypdm.Pypdm.PartToPart.iexch
      ~Pypdm.Pypdm.PartToPart.wait
      ~Pypdm.Pypdm.PartToPart.reverse_iexch
      ~Pypdm.Pypdm.PartToPart.reverse_wait


    .. rubric:: Information on Part 2 side

    .. automethod:: Pypdm.Pypdm.PartToPart.get_referenced_lnum2

    .. automethod:: Pypdm.Pypdm.PartToPart.get_unreferenced_lnum2

    .. automethod:: Pypdm.Pypdm.PartToPart.get_gnum1_come_from


    .. rubric:: Exchange

    .. automethod:: Pypdm.Pypdm.PartToPart.iexch

    .. automethod:: Pypdm.Pypdm.PartToPart.wait

    .. automethod:: Pypdm.Pypdm.PartToPart.reverse_iexch

    .. automethod:: Pypdm.Pypdm.PartToPart.reverse_wait


.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
