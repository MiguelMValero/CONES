.. _ptb:

Part to Block
=============


Initialization
""""""""""""""

.. doxygenfunction:: PDM_part_to_block_create

.. doxygenfunction:: PDM_part_to_block_create_from_distrib

.. doxygenfunction:: PDM_part_to_block_geom_create


Information on block-distributed frame
""""""""""""""""""""""""""""""""""""""

.. doxygenfunction:: PDM_part_to_block_n_elt_block_get

.. doxygenfunction:: PDM_part_to_block_block_gnum_get

.. doxygenfunction:: PDM_part_to_block_block_gnum_count_get

.. doxygenfunction:: PDM_part_to_block_distrib_index_get


Exchange
""""""""

.. doxygenfunction:: PDM_part_to_block_exch

.. doxygenfunction:: PDM_part_to_block_reverse_exch

.. doxygenfunction:: PDM_part_to_block_iexch

.. doxygenfunction:: PDM_part_to_block_iexch_wait

.. doxygenfunction:: PDM_part_to_block_reverse_iexch

.. doxygenfunction:: PDM_part_to_block_reverse_iexch_wait


Finalization
""""""""""""

.. doxygenfunction:: PDM_part_to_block_free
