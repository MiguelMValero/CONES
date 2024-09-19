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

# Exercise 1 : Mesh partitioning

+++ {"editable": false, "deletable": false}

It's time for some hands on experience with **ParaDiGM**!
Using the API referenced [here](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_formation/user_manual/partitioning/multipart.html#fortran-api),
you will have to fill in the code cells to partition a mesh, i.e. to cut it in subdomains that will be mapped onto the processors of a parallel machine.
In the first section, we generate a block-distributed cube mesh for you. In the next section, you'll start running the partitioning algorithm.
After that, you will be able to retrieve the arrays describing the partitioned mesh.

+++ {"editable": false, "deletable": false}

## Load magic commands
We start by loading the custom magic commands for the proper functioning of the Notebook.

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
```

+++ {"editable": false, "deletable": false}

Fortran forces to define the variables at the top of the program.
In this notebook, we define the variables you need for a given function call in a separate cell above the one you will fill in.

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 1

program pdm_t_mesh_partitioning_f

  use pdm
  use pdm_multipart
  use pdm_vtk
  use pdm_dcube_nodal_gen
  use pdm_dmesh_nodal
  use pdm_part_mesh_nodal
  use pdm_mesh_nodal
  use iso_c_binding
  use pdm_fortran
  use pdm_writer_wrapper
  use pdm_part_connectivity_transform
  use pdm_part_extension

  implicit none

  include "mpif.h"

  !-----------------------------------------------------------
  integer (c_int)    :: i

  ! MPI
  integer            :: code
  integer            :: i_rank
  integer            :: n_rank
  integer, parameter :: comm = MPI_COMM_WORLD

  ! Mesh generation with dcube_nodal_gen
  integer          :: n_x, n_y, n_z
  integer          :: elt_type, order
  double precision :: length
  double precision :: xmin, ymin, zmin
  type (c_ptr)     :: dcube
  type (c_ptr)     :: dmn

```

+++ {"editable": false, "deletable": false}

## Generate the mesh

In this section, **ParaDiGM** tools are used to generate a simple mesh for this exercise: a cube made of tetrahedra.
You have **nothing to do here**. Still if you are curious about this feature, you can have a look [here](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_formation/user_manual/simple_mesh_gen/dcube_nodal.html#fortran-api).

In your numerical simulation software you rarely generate a mesh.
This step actually generates a block-distributed mesh which is distributed in the same way as a mesh you would have **read in parallel**.

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 17

  ! Initialize MPI environment
  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  ! Generate block-distributed parallelepided mesh
  n_x      = 10
  n_y      = 10
  n_z      = 10
  length   = 1.
  xmin     = 0.
  ymin     = 0.
  zmin     = 0.
  elt_type = PDM_MESH_NODAL_TETRA4
  order    = 1
  call PDM_dcube_nodal_gen_create(dcube,     &
                                  comm,      &
                                  n_x,       &
                                  n_y,       &
                                  n_z,       &
                                  length,    &
                                  xmin,      &
                                  ymin,      &
                                  zmin,      &
                                  elt_type,  &
                                  order,     &
                                  PDM_OWNERSHIP_USER)

  call PDM_dcube_nodal_gen_build(dcube, dmn)

  call PDM_dcube_nodal_gen_dmesh_nodal_get(dcube, dmn)

  call PDM_dmesh_nodal_generate_distribution(dmn)

  call PDM_dcube_nodal_gen_free(dcube)

```

Here you can see that the mesh were stored in a Dirstibuted-Nodal-Mesh structure (`dmn`).
This is an internal mesh structure to **ParaDiGM** not for user purpose.
Each feature is made such that you can set the mesh using basic arrays.

+++ {"editable": false, "deletable": false}

Now that we have our mesh, let's partition it !

## Mesh partitioning

For mesh partitioning, as for all other **ParaDiGM** features, there are 5 main steps:
1. **create** the feature structure
2. **set** the data necessary to operate with that feature
3. **compute**, operate the algorithm of the feature
4. **get**, retrieve the output of the algorithm
5. **free** the memory allocated to operate the feature

Following this logic, let's start **creating** (step 1) the mesh partitioning structure for **homogeneously** balanced subdomains.

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 2

  type (c_ptr)                       :: mpart
  integer (c_int)                    :: n_domain = 1
  integer(kind=PDM_l_num_s), pointer :: n_part(:)        => null()  ! since there could be several domains, this is an array
  integer (c_int)                    :: merge_domains, part_method, part_size_method
  double precision,          pointer :: part_fraction(:) => null()
  integer (c_int)                    :: i_domain = 0
  integer (c_int)                    :: i_part   = 0
```

+++ {"editable": false, "deletable": false}

*Remark : since this is a basic example, we ask you to stick with the fixed values for n_domain, n_part, i_domain, i_part and merge_domains.
To get insight about the concepts behind those values you can have a look [here](#Annex-1)*

**ParaDiGM** offers multiple partitioning methods.
Here, we chose to partition the cube with the **Hilbert method**.
This method is favored within the **ParaDiGM** algorithms since it provides quickly a good load balance, though it does not ensure the connectedness of each subdomain.
To ensure the partitions are connected, you should use either
`PDM_SPLIT_DUAL_WITH_PARMETIS` or `PDM_SPLIT_DUAL_WITH_PTSCOTCH` which call the external libraries ParMETIS and PT-Scotch.

*Remark : In this exercise we do not provide weights for the partitioning.*

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 18

  ! Create partitioning structure
  allocate(n_part(n_domain))

  do i = 1, n_domain
    n_part(i) = 1
  end do

  merge_domains    = PDM_FALSE
  part_method      = PDM_SPLIT_DUAL_WITH_HILBERT
  part_size_method = PDM_PART_SIZE_HOMOGENEOUS
  part_fraction    => null() ! unused here since the subdomains are homogeneous
  call PDM_multipart_create() ! ??


```

+++ {"editable": false, "deletable": false}

After mapping the partitioned subdomains on the processors, it is interesting to renumber the entities
of the mesh on each processor for performance through cache blocking but it also provides interesting properties for the application.
This is an **advanced setting we won't be using here**, so we just specify that no renumbering should be performed.
<!-- You can here call the renumbering function but by telling it not to do any renumbering for a start. -->

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 3

  integer(kind=PDM_l_num_s), pointer    :: renum_cell_properties(:) => null()
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 19

  call PDM_multipart_set_reordering_options(mpart,                      &
                                            i_domain,                   &
                                            "PDM_PART_RENUM_CELL_NONE", &
                                            renum_cell_properties,      &
                                            "PDM_PART_RENUM_FACE_NONE")

```

+++ {"editable": false, "deletable": false}

Now that you have created a mesh partitioning structure `mpart`, you can **set** (step 2) the cube mesh to it.
For simplicity of the exercise, we here set the mesh using the **Dirstibuted-Nodal-Mesh** structure (`dmn`).
This is a pratice internal to **ParaDiGM** algorithms. In your software you would just set the mesh using basic arrays.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 20

  call PDM_multipart_dmesh_nodal_set() ! ??

```

+++ {"editable": false, "deletable": false}

At this point you have provided all the information necessary to run the mesh partitioning algorithm. You can call the function to
**compute** (step 3) the subdomains that make up the partitioned cube (i.e.`PDM_multipart_compute`).

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 21

  call ! compute ??
```

+++ {"editable": false, "deletable": false}

## Get the partitioned mesh

You can now **get** (step 4) the output mesh of the partitioning algorithm. Depending on the numerical method, the mesh has to be
described in a different way. For Finite-Element methods a nodal connectivity ([option 1](#Nodal-connectivity-(i.e.-Finite-Element-style)))) usually
suffices while for Finite-Volume methods all descending connectivities ([option 2](#Descending-connectivity-(i.e.-Finite-Volume-style))) are of interest.
Choose which one suits you best and go further in the exercise to the associated section.

### Nodal connectivity (i.e. Finite-Element style)

You choose to get the partitioned mesh in nodal connectivity, i.e. cell->vertex connectivity.

*Remark : The structure in **ParaDiGM** in which partitioned nodal meshes are stored is `part_mesh_nodal`.
Here we get this structure from `mpart` to have a direct access to the arrays we are interested in.
For more information about this structure, have a look [here](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_formation/user_manual/partitioning/multipart.html#id9)*

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 4

  type (c_ptr) :: pmn
```

+++ {"editable": false, "deletable": false}

Let's start with the **vertices** composing the subdomain. How many vertices are there? What are their global ids? What are their coordinates?

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 5

  double precision,      pointer :: coords(:,:)     => null()
  integer(c_int)                 :: n_vtx
  integer (pdm_g_num_s), pointer :: vtx_ln_to_gn(:) => null()
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 22

  ! use PDM_OWNERSHIP_USER
  call PDM_multipart_part_vtx_coord_get() ! ??


  call PDM_multipart_get_part_mesh_nodal(mpart,    &
                                         i_domain, &
                                         pmn,      &
                                         PDM_OWNERSHIP_USER)


  call PDM_part_mesh_nodal_vtx_g_num_get() ! ??

```

+++ {"editable": false, "deletable": false}

Let's move on to the **cells**. How are the vertices connected to form cells? What are their global ids? How many cells are there?

*Remark : since this is a basic example, we ask you to stick with the fixed value for i_section.
To get insight about the concept behind this value you can have a look [here](#Annex-1)*

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 6

  integer (c_int)                     :: n_elt
  integer(kind=PDM_l_num_s), pointer  :: elt_vtx(:)             => null()
  integer (pdm_g_num_s),     pointer  :: elt_ln_to_gn(:)        => null()
  integer(kind=PDM_l_num_s), pointer  :: parent_num(:)          => null()
  integer (pdm_g_num_s),     pointer  :: parent_entity_g_num(:) => null()
  integer (c_int)                     :: i_section = 0
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 23

  call PDM_part_mesh_nodal_section_n_elt_get() ! ??

  ! use PDM_OWNERSHIP_KEEP
  call PDM_part_mesh_nodal_section_std_get() ! elt->vtx ??

```

+++ {"editable": false, "deletable": false}

Now we write the mesh that we just got to be able to visualize it later on **(nothing to do)**.

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 7

  integer(pdm_l_num_s),      pointer :: pn_vtx(:)
  integer(pdm_l_num_s),      pointer :: pn_elt(:)

  type(PDM_pointer_array_t), pointer :: pcoords       => null()
  type(PDM_pointer_array_t), pointer :: pvtx_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer :: pelt_vtx      => null()
  type(PDM_pointer_array_t), pointer :: pelt_ln_to_gn => null()
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 24

  allocate(pn_vtx(1), &
           pn_elt(1))

  pn_vtx(1) = n_vtx
  pn_elt(1) = n_elt

  call PDM_pointer_array_create(pcoords,        1, PDM_TYPE_DOUBLE)
  call PDM_pointer_array_create(pvtx_ln_to_gn,  1, PDM_TYPE_G_NUM)
  call PDM_pointer_array_create(pelt_vtx,       1, PDM_TYPE_INT)
  call PDM_pointer_array_create(pelt_ln_to_gn,  1, PDM_TYPE_G_NUM)

  call PDM_pointer_array_part_set(pcoords,       0, coords)
  call PDM_pointer_array_part_set(pvtx_ln_to_gn, 0, vtx_ln_to_gn)
  call PDM_pointer_array_part_set(pelt_vtx,      0, elt_vtx)
  call PDM_pointer_array_part_set(pelt_ln_to_gn, 0, elt_ln_to_gn)

  call writer_wrapper(comm,          &
                      "visu",        &
                      "pmesh",       &
                      1,             &
                      pn_vtx,        &
                      pcoords,       &
                      pvtx_ln_to_gn, &
                      pn_elt,        &
                      pelt_vtx_idx,  &
                      pelt_vtx,      &
                      pelt_ln_to_gn, &
                      cell_t = elt_type)

  call PDM_pointer_array_free(pcoords)
  call PDM_pointer_array_free(pvtx_ln_to_gn)
  call PDM_pointer_array_free(pelt_vtx)
  call PDM_pointer_array_free(pelt_ln_to_gn)

  deallocate(pn_vtx, &
             pn_elt)

  call PDM_part_mesh_nodal_free(pmn)
```

+++ {"editable": false, "deletable": false}

### Descending connectivity (i.e. Finite-Volume style)

You chose to get the partitioned mesh in descending connectivity, i.e. **cell->face**, **face->vtx** connectivities.
Generic getters have been implemented in **ParaDiGM** for the connectivities and global identifier arrays.
Enumerators allow to specify which data is requested (see details in [documentation](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_formation/user_manual/partitioning/index.html#enumerators)). Here you will need for the mesh connectivity:
- **PDM_CONNECTIVITY_TYPE_CELL_FACE** : cell->face connectivity
- **PDM_CONNECTIVITY_TYPE_FACE_VTX**  : face->vertex connectivity
For the global identifier arrays you will use:
- **PDM_MESH_ENTITY_CELL**  : cell entity
- **PDM_MESH_ENTITY_FACE**  : face entity
- **PDM_MESH_ENTITY_VTX**  : vertex entity

Let's start from the top with **cell** data. How many cells are there? What are their global ids? Which faces compose the cells?

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 8

  integer(kind=PDM_g_num_s), pointer :: cell_ln_to_gn(:) => null()
  integer(c_int)                     :: n_cell
  integer(kind=PDM_l_num_s), pointer :: cell_face(:)     => null()
  integer(kind=PDM_l_num_s), pointer :: cell_face_idx(:) => null()
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 25

  ! use PDM_OWNERSHIP_KEEP
  call PDM_multipart_part_ln_to_gn_get() ! ??

  ! use PDM_OWNERSHIP_KEEP
  call PDM_multipart_part_connectivity_get() ! cell->face ??

```

+++ {"editable": false, "deletable": false}

For the **faces** we proceed in a similar way. How many faces are there? What are their global ids? Which vertices compose the faces?

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 9

  integer(kind=PDM_g_num_s), pointer :: face_ln_to_gn(:) => null()
  integer(c_int)                     :: n_face
  integer(kind=PDM_l_num_s), pointer :: face_vtx(:)      => null()
  integer(kind=PDM_l_num_s), pointer :: face_vtx_idx(:)
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 26

  ! use PDM_OWNERSHIP_KEEP
  call PDM_multipart_part_ln_to_gn_get() ! ??

  ! use PDM_OWNERSHIP_KEEP
  call PDM_multipart_part_connectivity_get() ! face->vtx ??

```

+++ {"editable": false, "deletable": false}

To finish with, we need to have the description of the **vertices**.

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 11

  double precision,      pointer :: coords(:,:)     => null()
  integer(c_int)                 :: n_vtx
  integer (pdm_g_num_s), pointer :: vtx_ln_to_gn(:) => null()
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 28

  ! use PDM_OWNERSHIP_KEEP
  call PDM_multipart_part_ln_to_gn_get() ! ??

  ! use PDM_OWNERSHIP_KEEP
  call PDM_multipart_part_vtx_coord_get() ! ??

```

+++ {"editable": false, "deletable": false}

Now we write the mesh that we just got to be able to visualize it later on **(nothing to do)**.

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 12

  integer(pdm_l_num_s),      pointer :: pn_vtx(:)
  integer(pdm_l_num_s),      pointer :: pn_elt(:)
  integer(pdm_l_num_s),      pointer :: pn_face(:)

  type(PDM_pointer_array_t), pointer :: pcoords        => null()
  type(PDM_pointer_array_t), pointer :: pvtx_ln_to_gn  => null()
  type(PDM_pointer_array_t), pointer :: pelt_vtx_idx   => null()
  type(PDM_pointer_array_t), pointer :: pelt_vtx       => null()
  type(PDM_pointer_array_t), pointer :: pelt_ln_to_gn  => null()
  type(PDM_pointer_array_t), pointer :: pcell_face_idx => null()
  type(PDM_pointer_array_t), pointer :: pcell_face     => null()
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 29

  allocate(pn_vtx(1), &
           pn_elt(1), &
           pn_face(1))

  pn_vtx(1)  = n_vtx
  pn_elt(1)  = n_cell
  pn_face(1) = n_face

  call PDM_pointer_array_create(pcoords,        1, PDM_TYPE_DOUBLE)
  call PDM_pointer_array_create(pvtx_ln_to_gn,  1, PDM_TYPE_G_NUM)
  call PDM_pointer_array_create(pelt_vtx_idx,   1, PDM_TYPE_INT)
  call PDM_pointer_array_create(pelt_vtx,       1, PDM_TYPE_INT)
  call PDM_pointer_array_create(pelt_ln_to_gn,  1, PDM_TYPE_G_NUM)
  call PDM_pointer_array_create(pcell_face_idx, 1, PDM_TYPE_INT)
  call PDM_pointer_array_create(pcell_face,     1, PDM_TYPE_INT)

  call PDM_pointer_array_part_set(pcoords,        0, coords)
  call PDM_pointer_array_part_set(pvtx_ln_to_gn,  0, vtx_ln_to_gn)
  call PDM_pointer_array_part_set(pelt_vtx_idx,   0, face_vtx_idx)
  call PDM_pointer_array_part_set(pelt_vtx,       0, face_vtx)
  call PDM_pointer_array_part_set(pelt_ln_to_gn,  0, cell_ln_to_gn)
  call PDM_pointer_array_part_set(pcell_face_idx, 0, cell_face_idx)
  call PDM_pointer_array_part_set(pcell_face,     0, cell_face)

  call writer_wrapper(comm,                            &
                      "visu",                          &
                      "pmesh",                         &
                      1,                               &
                      pn_vtx,                          &
                      pcoords,                         &
                      pvtx_ln_to_gn,                   &
                      pn_elt,                          &
                      pelt_vtx_idx,                    &
                      pelt_vtx,                        &
                      pelt_ln_to_gn,                   &
                      n_face         = pn_face,        &
                      pcell_face_idx = pcell_face_idx, &
                      pcell_face     = pcell_face)

  call PDM_pointer_array_free(pcoords)
  call PDM_pointer_array_free(pvtx_ln_to_gn)
  call PDM_pointer_array_free(pelt_vtx_idx)
  call PDM_pointer_array_free(pelt_vtx)
  call PDM_pointer_array_free(pelt_ln_to_gn)
  call PDM_pointer_array_free(pcell_face_idx)
  call PDM_pointer_array_free(pcell_face)

  deallocate(pn_vtx, &
             pn_elt, &
             pn_face)
```

+++ {"editable": false, "deletable": false}

## Execution and visualization

First, we finalize the the code you juste wrote by with the last step :  **free** (step 5).

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 30

  deallocate(n_part)
  call PDM_DMesh_nodal_free(dmn)
  call PDM_multipart_free() ! ??

  ! Finalize MPI environment
  call mpi_finalize(code)

  if (i_rank == 0) then
    print *, "End :)"
  endif

end program pdm_t_mesh_partitioning_f

```

+++ {"editable": false, "deletable": false}

Run the following cells to execute the program you just wrote and visualize the output partitioned mesh.

```{code-cell}
---
"deletable": false
---
%merge_code_blocks -l fortran -p exercise_1 -n 2 -v
```

```{code-cell}
---
"deletable": false
---
%%visualize
visu/PMESH.case : i_part
```

+++ {"editable": false, "deletable": false}

## Bonus : Extended partition

If you are reading this, you finished quickly the partitioning exercise. Thus, it means you understood well the **5 step scheme** for using **ParaDiGM** features.

*Remark : To do this bonus you need to have retrieved the mesh in descending connectivity. If you haven't done that yet, please comment your
work on nodal connectivities and get the mesh in descending connectivity first.*

In this bonus, we want to get one layer of extended cells by nodes for our mesh partitions.
This bonus is not guided, so you should have a close look at the [documentation](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_formation/user_manual/partitioning/part_extension.html#fortran-api).

### Step 1

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 13

  type(c_ptr) :: part_ext = C_NULL_PTR
  integer     :: extend_type
  integer     :: depth
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 31

  extend_type = PDM_EXTEND_FROM_VTX
  depth       = 1
  ! use PDM_OWNERSHIP_KEEP
  call PDM_part_extension_create () ! ??
```

+++ {"editable": false, "deletable": false}

### Step 2

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 14

  integer(PDM_l_num_s), pointer   :: vtx_part_bound_proc_idx(:)  => null()
  integer(PDM_l_num_s), pointer   :: vtx_part_bound_part_idx(:)  => null()
  integer(PDM_l_num_s), pointer   :: vtx_part_bound(:)           => null()
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 32

  ! extension by vertex : PDM_MESH_ENTITY_VTX
  ! use PDM_OWNERSHIP_KEEP
  call PDM_multipart_part_graph_comm_get() ! ??

  ! set the above arrays
  call PDM_part_extension_part_bound_graph_set() ! ??

  ! set cell->face connectivity
  call PDM_part_extension_connectivity_set() ! ??

  ! set face->vertex connectivity
  call PDM_part_extension_connectivity_set() ! ??

  call PDM_part_extension_vtx_coord_set() ! ??

  ! set cell global identifier array
  call PDM_part_extension_ln_to_gn_set() ! ??

  ! set face global identifier array
  call PDM_part_extension_ln_to_gn_set() ! ??

  ! set vertex global identifier array
  call PDM_part_extension_ln_to_gn_set() ! ??
```

+++ {"editable": false, "deletable": false}

### Step 3

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 33

call ! compute ??
```

+++ {"editable": false, "deletable": false}

### Step 4

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 15

  integer                         :: n_cell_ext
  integer(pdm_l_num_s), pointer   :: cell_face_ext(:)      => null()
  integer(pdm_l_num_s), pointer   :: cell_face_ext_idx(:)  => null()
  integer(PDM_g_num_s), pointer   :: cell_ln_to_gn_ext(:)  => null()

  integer                         :: n_face_ext
  integer(pdm_l_num_s), pointer   :: face_vtx_ext(:)       => null()
  integer(pdm_l_num_s), pointer   :: face_vtx_ext_idx(:)   => null()
  integer(PDM_g_num_s), pointer   :: face_ln_to_gn_ext(:)  => null()

  integer                         :: n_vtx_ext
  double precision,     pointer   :: vtx_coord_ext(:,:)    => null()
  integer(PDM_g_num_s), pointer   :: vtx_ln_to_gn_ext(:)   => null()
```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 34

  ! Cell
  call PDM_part_extension_ln_to_gn_get () ! ??

  call PDM_part_extension_connectivity_get() ! cell->face ??

  ! Face
  call PDM_part_extension_ln_to_gn_get () ! ??

  call PDM_part_extension_connectivity_get () ! face->vtx ??

  ! Vertices
  call PDM_part_extension_vtx_coord_get () ! ??

  call PDM_part_extension_ln_to_gn_get () ! ??
```

+++ {"editable": false, "deletable": false}

### Step 5

Before handling the allocated memory, we output the mesh partition extension **(nothing to do)**.

```{code-cell}
---
"editable": false
"deletable": false
---
%%code_block -p exercise_1 -i 16

  double precision, pointer          :: total_coords(:,:)      => null()
  integer(c_int)                     :: total_n_vtx

  integer (pdm_g_num_s), pointer     :: total_vtx_ln_to_gn(:)  => null()

  integer(kind=PDM_g_num_s), pointer :: total_face_ln_to_gn(:) => null()
  integer(c_int)                     :: total_n_face
  integer(kind=PDM_l_num_s), pointer :: total_face_vtx(:)      => null()
  integer(kind=PDM_l_num_s), pointer :: total_face_vtx_idx(:)  => null()

  integer(kind=PDM_g_num_s), pointer :: total_cell_ln_to_gn(:) => null()
  integer(c_int)                     :: total_n_cell
  integer(kind=PDM_l_num_s), pointer :: total_cell_face(:)     => null()
  integer(kind=PDM_l_num_s), pointer :: total_cell_face_idx(:) => null()

  type(my_field_t)                   :: elt_fields(1)
  double precision, pointer          :: total_cell_color(:) => null()

```

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 35

  ! Visualisation
  total_n_cell = n_cell + n_cell_ext
  total_n_face = n_face + n_face_ext
  total_n_vtx  = n_vtx  + n_vtx_ext

  allocate(total_vtx_ln_to_gn(total_n_vtx), &
           total_coords(3,total_n_vtx),     &
           total_face_ln_to_gn(total_n_face), &
           total_face_vtx_idx (total_n_face+1), &
           total_cell_ln_to_gn(total_n_cell), &
           total_cell_face_idx(total_n_cell+1))

  ! Cell
  do i = 1, n_cell
    total_cell_ln_to_gn(i) = cell_ln_to_gn(i)
  end do
  do i = 1, n_cell_ext
    total_cell_ln_to_gn(n_cell + i) = cell_ln_to_gn_ext(i)
  end do

  do i = 1, n_cell+1
    total_cell_face_idx(i) = cell_face_idx(i)
  end do
  do i = 1, n_cell_ext+1
    total_cell_face_idx(n_cell + i) = cell_face_idx(n_cell + 1) + cell_face_ext_idx(i)
  end do

  ! Face
  do i = 1, n_face
    total_face_ln_to_gn(i) = face_ln_to_gn(i)
  end do
  do i = 1, n_face_ext
    total_face_ln_to_gn(n_face + i) = face_ln_to_gn_ext(i)
  end do

  do i = 1, n_face+1
    total_face_vtx_idx(i) = face_vtx_idx(i)
  end do
  do i = 1, n_face_ext+1
    total_face_vtx_idx(n_face + i) = face_vtx_idx(n_face + 1) + face_vtx_ext_idx(i)
  end do

  ! Vertex
  do i = 1, n_vtx
    total_vtx_ln_to_gn(i) = vtx_ln_to_gn(i)
  end do
  do i = 1, n_vtx_ext
    total_vtx_ln_to_gn(n_vtx + i) = vtx_ln_to_gn_ext(i)
  end do

  do i = 1, n_vtx
    total_coords(1:3, i) = coords(1:3, i)
  end do
  do i = 1, n_vtx_ext
    total_coords(1:3, n_vtx + i) = vtx_coord_ext(1:3, i)
  end do

  allocate(total_face_vtx (total_face_vtx_idx (total_n_face+1)), &
           total_cell_face(total_cell_face_idx(total_n_cell+1)))

  ! Cell
  do i = 1, cell_face_idx(n_cell+1)
    total_cell_face(i) = cell_face(i)
  end do
  do i = 1, cell_face_ext_idx(n_cell_ext+1)
    total_cell_face(cell_face_idx(n_cell+1) + i) = cell_face_ext(i)
  end do

  ! Face
  do i = 1, face_vtx_idx(n_face+1)
    total_face_vtx(i) = face_vtx(i)
  end do
  do i = 1, face_vtx_ext_idx(n_face_ext+1)
    total_face_vtx(face_vtx_idx(n_face+1) + i) = face_vtx_ext(i)
  end do

  call PDM_pointer_array_create(pcoords,        1, PDM_TYPE_DOUBLE)
  call PDM_pointer_array_create(pvtx_ln_to_gn,  1, PDM_TYPE_G_NUM)
  call PDM_pointer_array_create(pelt_vtx_idx,   1, PDM_TYPE_INT)
  call PDM_pointer_array_create(pelt_vtx,       1, PDM_TYPE_INT)
  call PDM_pointer_array_create(pelt_ln_to_gn,  1, PDM_TYPE_G_NUM)
  call PDM_pointer_array_create(pcell_face_idx, 1, PDM_TYPE_INT)
  call PDM_pointer_array_create(pcell_face,     1, PDM_TYPE_INT)

  call PDM_pointer_array_part_set(pcoords,        0, total_coords)
  call PDM_pointer_array_part_set(pvtx_ln_to_gn,  0, total_vtx_ln_to_gn)
  call PDM_pointer_array_part_set(pelt_vtx_idx,   0, total_face_vtx_idx)
  call PDM_pointer_array_part_set(pelt_vtx,       0, total_face_vtx)
  call PDM_pointer_array_part_set(pelt_ln_to_gn,  0, total_cell_ln_to_gn)
  call PDM_pointer_array_part_set(pcell_face_idx, 0, total_cell_face_idx)
  call PDM_pointer_array_part_set(pcell_face,     0, total_cell_face)

  allocate(pn_vtx(1), &
           pn_elt(1), &
           pn_face(1))

  pn_vtx(1)  = total_n_vtx
  pn_elt(1)  = total_n_cell
  pn_face(1) = total_n_face


  allocate(total_cell_color(total_n_cell))
  total_cell_color(1:n_cell)            = 2*i_rank
  total_cell_color(n_cell:total_n_cell) = 2*i_rank+1

  call pdm_pointer_array_create(elt_fields(1)%pa, &
                                n_part(1),        &
                                PDM_TYPE_DOUBLE)
  call pdm_pointer_array_part_set(elt_fields(1)%pa, &
                                  0,                &
                                  total_cell_color)
  elt_fields(1)%name = "extension"

  call writer_wrapper(comm,           &
                      "visu",         &
                      "pext",         &
                      1,              &
                      pn_vtx,         &
                      pcoords,        &
                      pvtx_ln_to_gn,  &
                      pn_elt,         &
                      pelt_vtx_idx,   &
                      pelt_vtx,       &
                      pelt_ln_to_gn,  &
                      -1,             &
                      pn_face,        &
                      pcell_face_idx, &
                      pcell_face,     &
                      elt_field=elt_fields)

  deallocate(pn_vtx, &
             pn_elt, &
             pn_face)

  deallocate(total_vtx_ln_to_gn,  &
             total_coords,        &
             total_face_ln_to_gn, &
             total_face_vtx_idx,  &
             total_face_vtx,      &
             total_cell_ln_to_gn, &
             total_cell_face_idx, &
             total_cell_face)

```

+++ {"editable": false, "deletable": false}

Now you can do step 5.

```{code-cell}
---
"deletable": false
---
%%code_block -p exercise_1 -i 36

  ! free
  call PDM_part_extension_free () ! ??

  deallocate(n_part, &
             part_fraction)
  call PDM_DMesh_nodal_free(dmn)
  call PDM_multipart_free(mpart)

  ! Finalize MPI environment
  call mpi_finalize(code)

  if (i_rank == 0) then
    print *, "End :)"
  endif

end program pdm_t_mesh_partitioning_f

```

+++ {"editable": false, "deletable": false}

## Execution and visualization

Run the following cells to execute the program you just wrote and visualize the mesh partition extension.

```{code-cell}
---
"deletable": false
---
%merge_code_blocks -l fortran -p exercise_1 -n 2 -c
```

```{code-cell}
---
"deletable": false
---
%%visualize
visu/PEXT.case : extension
```

+++ {"editable": false, "deletable": false}

## Annex 1

In some cases, the mesh is an assembly of several sub-meshes. These are called *domains*.
In the figure bellow, we can see a mesh made of two domains.

<img src="mesh.png" width="180">

Each *domain* is partitioned in subdomains which
are mapped to the processors of the parallel machine. On a processor the subdomain (of a mesh or a domain) can be subdivided in *parts*.
This figure shows rank 0 for the above mesh which has a subdomain of each domain with two parts for the subdomain of domain 1.

<img src="processor.png" width="180">

A mesh can be composed of several element types (tetrahedra, hexahedra, prisms...). In certain settings, the mesh definition for each specific element type
is stored in a separate *section*. So in a *section* one will find data for a specific element type.
Here we can see part 1 of the subdomain on rank 0 of domain 1 which has two sections.

<img src="part.png" width="180">


# Exercise 2

You can now move on to [Exercise 2](./../03_Exercise_2/exercise_2.ipynb).
