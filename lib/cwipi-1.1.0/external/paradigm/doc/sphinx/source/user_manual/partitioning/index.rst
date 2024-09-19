.. _partitioning:

############
Partitioning
############

.. container:: toc-cards

  .. container:: card

    :ref:`Multipart <multipart>`
      Multi-domain mesh partitioning


  .. container:: card

    :ref:`Extract part <extract_part>`
      Extraction and redistribution of mesh partitions


  .. container:: card

    :ref:`Connectivity transformation <connec_transform>`
      Utilities for processing unstructured mesh connectivities (partitioned and distributed)

  .. container:: card

    :ref:`Part extension <part_extension>`
      Generate extended mesh partitions

Enumerators
-----------

Here we present enumerators that are usefull for the features detailed in this section.

.. doxygenenum:: PDM_connectivity_type_t

.. doxygenenum:: PDM_mesh_entities_t

* renumbering?
* agglomeration?
* domain interfaces/joints?


.. toctree::
   :caption: Partitioning
   :maxdepth: 1
   :hidden:

   multipart
   extract_part
   connec_transform
   part_extension

