.. _faq:

CWIPI's FAQ
===========

This is a list of Frequently Asked Questions about CWIPI.  Feel free to
suggest new entries!

Why...
------

... does CWIPI fail with the given data entry?
   Before accusing CWIPI have you verified if your data entry complies with the CWIPI
   data entry written in the documentation?

... does CWIPI fail at certain runs?
   This is due to the asynchronous exchanges. It is not possible to get a data that has
   no been set by the other code. A possible fix is to do the get in a while loop to wait
   until the data is available.

... does my Python application behave erratically?
   Python's garbage collector may over zealously decide to reuse some memory areas if no reference to that area is kept explicitly. This results in random behavior of your code from one run to the next, and possibly catastrophic failure. To prevent this, make sure you keep explicit references to all data passed to CWIPI via the Python API. These references can be stored in your code's internal classes or in a dictionary.

