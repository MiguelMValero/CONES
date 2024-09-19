.. _philosophy:

Philosophy
##########

In this section we, the CWIPI developers, would like to raise awareness about philosophy that guides our coding choices.
The library is quite low level. Thus, to ensure it is standalone, we restrict as much as possible the mandatory dependencies.
For that same reason no particular data structure is used. Only simple types and arrays (C-contiguous) are used.
For instance in Python, Numpy arrays are used for all physical and geometrical data, otherwise only basic types (list, int, bool...) are used.
Since the library does not allocate the memory of the data that users can get, users need to manage the memory allocation and deallocation themselves.
Let us note that coordinates arrays always have 3 components (x, y and z), even if the coupling interface is of lower dimension than 3D.
