---
jupytext:
  text_representation:
    extension: '.md'
    format_name: myst
    format_version: '0.7'
    jupytext_version: 1.4.0+dev
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

+++ {"editable": false, "deletable": false}

# Exercise 2 : Coupling with a deformed mesh over time

+++ {"editable": false, "deletable": false}

Now that you know how to set up a basic coupling, let's go further by doing several coupling iterations.
At each iteration, the coupling interface mesh of `code1` is deformed.

+++ {"editable": false, "deletable": false}

## Load magic commands

As usual we start by loading the custom magic commands.

```{code-cell}
---
"editable": false
"deletable": false
---
import os, sys
module_path = os.path.abspath(os.path.join('../../utils'))
if module_path not in sys.path:
    sys.path.append(module_path)
```

```{code-cell}
---
"editable": false
"deletable": false
---
%reload_ext visu_magics
%reload_ext code_magics
%reload_ext figure_magics
```

+++ {"editable": false, "deletable": false}

## Initialization

Since the set up is roughly the same as in the previous exercise, it is not split up in as many code cells. Here we ask you to:
  - Initialize MPI
  - Initialize CWIPI for `code1`
  - Set up the coupling : What value should be chosen for `CWP_Dynamic_mesh_t` since the coupling interface mesh of `code1` is deformed over time?
  - Ask for visualization outputs

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 1

#!/usr/bin/env python

import mpi4py.MPI as MPI

# Initialize MPI
comm = MPI.COMM_WORLD
i_rank = comm.rank
n_rank = comm.size

import numpy as np
import math

# pycwp
try:
    from pycwp import pycwp
    if i_rank == 0:
        print("Yes, we have found the module :D")
except:
    if i_rank == 0:
        print("Oh no, couldn't find the module :'(")

# pypdm
import Pypdm.Pypdm as PDM

n_code = 1
code_name = ["code1"]
is_active_rank = True
intra_comm = pycwp.init() # ??

coupled_code_name = ["code2"]
n_part = 1
oupling_name = "coupling"
recv_freq_type = pycwp.TIME_EXCH_USER_CONTROLLED
cpl = pycwp.Coupling() # ??

format_option = "text"
cpl.visu_set() # ??
```

+++ {"editable": false, "deletable": false}

## What changes when the mesh is deformed over time?

Let's have a look again at the pseudo-code of the introduction.

<!-- ```{prf:algorithm} basic coupling algorithm

**Inputs** Given $code1$ with a mesh $m1$ on which a field that will be sent is defined $sf1$. Given $code2$ with a mesh $m2$ on which a field that will be received is defined $rf2$.

**Output** $rf2$, which is $sf1$ interpolated on $m2$

1. Initialize CWIPI
2. Set coupling between $code1$ and $code2$
3. Describe codes:
   1. $code1$ has a mesh $m1$ on which we define a field $sf1$
   2. $code2$ has a mesh $m2$ and a receive buffer field $rf2$
4. Operate solver iterations:
   1. $code1$ sends $sf1$
   2. $code2$ receives $rf2$
5. Finalize CWIPI
``` -->
**Inputs** Given `code1` with a mesh `m1` on which a field that will be sent is defined `sf1`. Given `code2` with a mesh `m2` on which a field that will be received is defined `rf2`.

**Output** `rf2`, which is `sf1` interpolated on `m2`

1. Initialize CWIPI
2. Set coupling between `code1` and `code2`
3. Describe codes:
   1. `code1` has a mesh `m1` on which we define a field `sf1`
   2. `code2` has a mesh `m2` and a receive buffer field `rf2`
4. Operate solver iterations:
   1. `code1` sends `sf1`
   2. `code2` receives `rf2`
5. Finalize CWIPI


In this exercise `code1` receives a field send by `code2`.
Here we decide to deform `m1` at each `code1` iteration.

*What does that change in our coupling code?
What would happen if `code1` would send `sf1`?*

### Mesh

First we use a simple mesh generation function from ParaDiGM to create our coupling interface mesh : a square **(nothing to do)**.
It is composed of triangle elements (i.e. `BLOCK_FACE_TRIA3`).

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 2

# Create mesh :
mesh = PDM.generate_mesh_rectangle_simplified(intra_comm[0], 10)

n_vtx       = mesh["n_vtx"]
n_elt       = mesh["n_elt"]
coords      = mesh["coords"]
elt_vtx_idx = mesh["elt_vtx_idx"]
elt_vtx     = mesh["elt_vtx"]
```

+++ {"editable": false, "deletable": false}

The mesh will change at each iteration. Since it is deformed, only its coordinates change.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 3

cpl.mesh_interf_vtx_set() # ??

block_id = cpl.mesh_interf_block_add() # ??

cpl.mesh_interf_block_std_set() # ??

# coupling interface finalization ??
```

+++ {"editable": false, "deletable": false}

### Field

In this exercise, `code1` receives a field from `code2`.
The mesh changes at each iteration but here we decided that the field `code2` sends wouldn't.
In a real case application, `code2` would send a different field a each iteration.
Since the mesh topology does not change, the coupling code would be the same since a pointer on the field is provided and there are no internal copies of the provided fields inside CWIPI.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 4

n_components = 1
field_name = "a super fancy field" # ??
field = cpl.field_create() # ??

n_vtx = len(coords)//3
field_data = np.empty(n_vtx*n_components, dtype=np.double)

field.data_set() # ??
```

+++ {"editable": false, "deletable": false}

In the case the topology of your mesh changes at each iteration, the new field array will be set at each iteration.
It is important to know that the field should still be created before starting the first time step.

## Time iterations

At the beginning of each coupling iteration, we begin a new time step using **time_step_beg** which we will terminate
at the end of the iteration with **time_step_end**.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 6

itdeb = 1
itend = 10
ttime = 0.0
dt    = 0.1
alpha = math.degrees(2.0)
sina  = math.sin(alpha)
cosa  = math.cos(alpha)

for it in range(itdeb, itend+1):

  ttime = (it-itdeb)*dt

  pycwp.time_step_beg() # ??
```

+++ {"editable": false, "deletable": false}

Let's rotate the mesh of `code1` with respect to `code2`.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 7

  if (it > itdeb):
    for i in range(n_vtx):
      x = coords[i * 3    ]
      y = coords[i * 3 + 1]
      coords[i * 3    ] = cosa*x - sina*y
      coords[i * 3 + 1] = sina*x + cosa*y
      field_data[i] = coords[3 * i]

```

+++ {"editable": false, "deletable": false}

The aim is to interpolate the field of `code2` onto the mesh of `code1`.
*What happens to the interpolation weights when the mesh of `code1` is deformed?
Thus, what does that induce in your code?*

The chosen tolerance does not change here over time, so we set it before the iteration loop.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 5

property_name = "tolerance"
cpl.spatial_interp_property_set() # ??
```

+++ {"editable": false, "deletable": false}

But the weights need to be computed at each iteration after the mesh has been deformed, so that is done in the iteration loop.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 8

  cpl. # spatial interpolation weights ??

```

+++ {"editable": false, "deletable": false}

Now we receive the field send by `code2`.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 9

  field. # receive field ??

  field. # wait receive field ??
```

+++ {"editable": false, "deletable": false}

Earlier we set a tolerance for the localization algorithm.
To check if that tolerance was large enough, the function **n_uncomputed_tgts_get** can be called to retrieve the number of unlocated vertices of the coupling interface of `code1`.
To know which vertices were unlocated the **uncomputed_tgts_get** is called.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1_code_1 -i 13

  n_uncomputed_tgts = field. # number of uncomputed targets ??
  uncomputed_tgts   = field. # array of uncomputed targets ??
```

+++ {"editable": false, "deletable": false}

Let's have a sneak peek in this algorithm through this animation which will help you understand what we mean by unlocated points.

```{code-cell}
---
"editable": false
"deletable": false
---
%%localization
unlocated
```

+++ {"editable": false, "deletable": false}

*Spoiler : At the end of the exercise you will see that since the coupling interface mesh of `code1` moves
there are unlocated points with the tolerance set to 0.001. Increasing it will eventually let all points be located
but at the cost of the time taken by the algorithm. You call play around with the tolerance once you finish the exercise.*

Let's end the iteration.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 10

  pycwp.time_step_end() # ??
```

+++ {"editable": false, "deletable": false}

## Finalize

Let us finish the coupling by freeing the memory allocated for it and ending this program.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_2_code_1 -i 11

# Delete field

# delete mesh interface ??

# Delete the coupling

# finalize pycwp ??

MPI.Finalize()

```

+++ {"editable": false, "deletable": false}

## Execution and visualization

Run the following cells to execute to program you just wrote and visualize the basic coupling you implemented.

```{code-cell}
---
"deletable": false
---
%merge_code_blocks -l python -p exercise_2_code_1 -n 1 -v -c
```

```{code-cell}
---
"deletable": false
---
%%visualize_dynamic
cwipi_writer/coupling_code1_code2/CHR.case : r_a~super~fancy~field1
cwipi_writer/coupling_code2_code1/CHR.case : s_a~super~fancy~field1
```

+++ {"editable": false, "deletable": false}

<span style="color:red">*You reached the end of this training. Congratulations ! Feel free to give us feedback.*</span>
