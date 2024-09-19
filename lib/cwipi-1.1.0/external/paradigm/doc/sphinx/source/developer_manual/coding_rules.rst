.. _coding_rules:

#########################
Coding rules & guidelines
#########################

Structure
---------

Most features in ParaDiGM are structured as the following five steps:

  1. **create**  : initialize the structure associated to the feature
  2. **set**     : provide the input data to the structure
  3. **compute** : run the algorithm
  4. **get**     : retrieve the output data
  5. **free**    : free the memory of the feature structure

Some features consist in establishing a communication graph between entities distributed on different processes.
This communication graph is in the form of a :ref:`Part-to-Part <ptp>` object.
This Part-to-Part object holds everything the user needs for addressing arrays both on the sender and receive side (except from the geometric mapping data).
Therefore, Part-to-Part accessors should always be used instead of adding redundant "src_to_tgt(_idx)_get" or "tgt_to_src(_idx)_get" functions.

All features should support meshes with several partitions per MPI rank (``n_part > 1``).
Ask J. Coulet and B. Maugars if you want to know whether you should handle multiple domains.

Indentation
-----------

C, Fortran and Python sources use 2 **spaces** for indentation.


Naming
------

Function, variable and C struct names should follow the ``snake_case`` convention.
Python class names follow the ``CamelCase`` convention.

.. Function names should respect the following templates:
.. ...

Function signatures
-------------------

In functions signatures featuring a connectivity and its index, always put the index argument before the connectivity argument:

.. list-table::
  :widths: 50 50
  :header-rows: 1

  * - |ok| Do
    - |ko| Don't
  * - .. code:: c

        my_function(...,
                    connect_idx,
                    connect,
                    ...);
    - .. code:: c

        my_function(...,
                    connect,
                    connect_idx,
                    ...);


0-based vs 1-based integer arrays
---------------------------------

All arrays with values that are addressable indices should be documented as being 0-based or 1-based in comments as well as in functions signatures.

.. Add with an array is 0 or 1-based in the function signature documentation. If it is 0-based, refer to it as ID rather than number.

