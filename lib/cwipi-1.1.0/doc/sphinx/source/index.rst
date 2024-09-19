*********************************
Welcome to CWIPI's documentation!
*********************************

**CWIPI** (Coupling With Interpolation Parallel Interface) is a C/C++ parallel coupling library under LGPL3 license.
It allows fully parallel data exchanges based on distributed mesh definition. Those meshes are differently partitioned
on several processes. Spatial interpolation is performed between non-coincident meshes on-the-fly. Arbitrarily many codes can be coupled
using this library. The philosophy of CWIPI is to let the parallelism be fully transparent for the user.

The library does not rely on a central coupling instance. The coupling is minimally invasive because it suffices
to call primitives of the library in the codes. Still, such a central supervisor can be set up by writing a Python supervisor script
calling CWIPI's Python API.

####################
General presentation
####################

:ref:`Code philosophy <philosophy>`

:ref:`Installation guide <installation>`

:ref:`License <license>`

###########
User manual
###########

:ref:`Quick Start <quick_start>` references basic information for inexperienced CWIPI users.

:ref:`Old CWIPI <old_cwipi>` is the documentation of the 0.x version of CWIPI. It describes most of the high-level API provided by CWIPI.

:ref:`Old to New <old_to_new>` allows to find equivalences between functions of versions 0.x and 1.x.

:ref:`New CWIPI <new_cwipi>` is the documentation of the 1.x version of CWIPI. It describes most of the high-level API provided by CWIPI.

:ref:`Client-Server <client_server_cwipi>` is the documentation of a client-server mode of CWIPI (based upon the 1.x API).

:ref:`FAQ <faq>` is a compilation of frequently asked questions.

.. toctree::
  :hidden:
  :maxdepth: 1
  :caption: Summary

  philosophy
  installation
  license
  quick_start
  old_cwipi/old_cwipi
  old_to_new/old_to_new
  new_cwipi/new_cwipi
  client_server_cwipi/client_server_cwipi
  faq

##################
Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



