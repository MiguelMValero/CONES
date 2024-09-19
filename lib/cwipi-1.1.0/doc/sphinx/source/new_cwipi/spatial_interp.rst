Since version 1.0, CWIPI offers multiple spatial interpolation methods.
A specific spatial interpolation method can be affected to each Coupling object.

Each method has a set of properties that can be adjusted (with ``CWP_Spatial_interp_property_set`` in C).
Note that some method do not support all types of dof location.
Finally, each method relies on some geometric mapping data to perform spatial interpolation.
This data is computed by CWIPI and can be accessed to design :ref:`customized interpolation functions <Advanced functionalities>`.

.. TODO: maybe mettre ce bloc warning à un autre endroit?

.. warning::
  Some spatial interpolation methods leave some target dofs unmapped if the interface meshes overlap partially.
  The interpolated field for such targets is simply not computed.
  Only the first :math:`n_{mapped} \times n_{comp}` values of the receiving data array are thus filled by CWIPI, where :math:`n_{mapped}` is the number of mapped targets (returned by ``CWP_N_computed_tgts_get``), and :math:`n_{comp}` the number of field components (specified when creating the field and returned by ``CWP_Field_n_components_get``).
  Please refer to section :ref:`Advanced functionalities` for some examples of use.

The spatial interpolation methods can be divided into the following four families.

Interpolation from nearest neighbors
------------------------------------

This family includes two methods:

  * ``CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES``: each target dof is mapped to its ``n_neighbors`` nearest source dofs.

  * ``CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES``: each source dof is mapped to its ``n_neighbors`` nearest target dofs (note that some targets may be mapped to zero source).

.. TODO: schéma?

**Default interpolation**
  Both methods perform a Least Square interpolation weighted by inverse squared distance (see `[Nealen 2004] <http://www.nealen.de/projects/mls/asapmls.pdf>`_ for more detail).

**Properties**

  The following properties can be adjusted by the user:

    * ``n_neighbors``: number of neighbors to interpolate from (`integer`, default value: 4)

    * ``polyfit_degree``: degree of polynomial fit (`integer`, default value: 1)

**Acceptable degrees-of-freedom**

  Since these are point cloud-based methods, they support all kinds of dof locations (mesh nodes, cell centers and user-defined points) for both source and target fields.


Interpolation from mesh intersection
------------------------------------

This family includes only one method: ``CWP_SPATIAL_INTERP_FROM_INTERSECTION`` which consists on computing the overlap fractions between the elements of both interface meshes.
It is suited for conservative interpolation.

**Default interpolation**

  The interpolated field at target :math:`T` is computed as follows

  .. math::

   f(T) = \frac{\sum_{S} w_S f(S)}{\sum_{S} w_S},

  where :math:`S` designates a source element and :math:`w_S` is the volume (or area for surface coupling interfaces) of :math:`T \cap S`.

  Note that for each mapped target element, the interpolation weights sum to one.
  Target elements that intersect no source element are simply not mapped.


**Properties**

  .. _bbox paragraph:

  The intersection algorithm uses efficient bounding box collisions detection as a first step.
  The boxes can be inflated using a relative tolerance.
  This is especially useful for non-planar, surface interface meshes where alignment with cartesian axes might cause some detection misses if the tolerance is set too low.
  The ``tolerance`` property is of type `double` and is set to 0.001 by default, but be can adjusted as well.

**Acceptable degrees-of-freedom**

  Naturally, this method cannot be used to interpolate fields with dofs located at user-defined points.
  Note that in version 1.0 this interpolation method is also not available for node-based source and target fields, but it will be soon.

.. warning::
  This method is only available for *surface* and *volume* interface meshes.


Interpolation from location
---------------------------

This family includes three methods:

  * ``CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE``: only the target dofs contained inside the bounding box of at least one source mesh element are mapped. This candidate-detection step is performed efficiently using a distributed point octree structure. The exact location of these targets are then computed. If a target is associated to multiple source elements, only the nearest one will be kept (usually the one that contains that target).

  * ``CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE``: similar to the previous method except a bounding box tree is used for the first step instead.

  * ``CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT``: similar to the first method except for the candidate-detection step which guarantees all targets are mapped. This method is therefore slower than the two previous ones in cases where the interface meshes overlap only partially.

.. note::
  These first two methods are equivalent to the spatial interpolation method of versions 0.x but have been re-implemented to improve performance.
  Both methods yield the same results but the octree generally performs better.

**Default interpolation**

  If the source field dofs are at cell centers, each mapped target simply gets the field value of its associated source mesh element.
  If the source field dofs are located at nodes, the generalized barycentric coordinates of the target are used as interpolation weights if the source field dofs are located at nodes.


**Properties**

  As for the :ref:`intersection method <bbox paragraph>`, the first step of the localization method uses bounding boxes inflated by a relative geometric tolerance.
  The ``tolerance`` property is of type `double` and is set to 0.001 by default but be can adjusted as well.
  Again, it is recommended to tune this property for non-planar, surface interface meshes when using either of the first two methods, in order to avoid true positive detection errors.
  The third method is however much less sensitive to this property so it should be kept to a low value to avoid detecting too many false positive candidates, which could impair performance.

**Acceptable degrees-of-freedom**

  This interpolation method requires the source field dofs to be located at either mesh nodes or cell centers.
  However, target field dofs can also be located at user-defined points.


Interpolation from identity
---------------------------

This family includes only one method: ``CWP_SPATIAL_INTERP_FROM_IDENTITY`` which consists in mapping each target dof to the source dof with same global id.
It might be useful for conforming interface meshes with different MPI partitioning.

**Default interpolation**
  Each target dof is assigned the field value of its mapped source dof.

  .. Note that if a source dof with same global id

**Properties**

  This method has no properties.

**Acceptable degrees-of-freedom**

  This method supports all kinds of dof location.
