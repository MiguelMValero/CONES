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

# Exercise 0 : Get started with the training notebooks

+++ {"editable": false, "deletable": false}

## Notebooks

This training course uses [Jupyter Notebooks](https://jupyter.org/) as a support.

A notebook is simply a document divided into cells of different kinds:
- Text cells (Markdown) : they contain the instructions you should follow to carry out the exercises. These cells are protected against edition but feel free to create new text cells if you wish to take notes.
- Python cells : they contain basic Python instructions (essentially for loading the modules required to run the notebooks) and are protected against edition as well
- `%%code_block` cells : they contain the pieces of code you will write during the exercises.
- `%merge_code_blocks` cells : they let you run the code you wrote piece by piece.
- `%visualize` cells : they let you visualize the fruits of your labor interactively in the notebook.

You will progress through the notebook step by step by running cells one after the other.


+++ {"editable": false, "deletable": false}

## Custom cells

+++ {"editable": false, "deletable": false}

The last three kinds of cells are custom cells developed specifically for this training course.
To enable them, we first need to load the appropriate modules.

```{code-cell} ipython3
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
%reload_ext code_magics
%reload_ext visu_magics
```

+++ {"editable": false, "deletable": false}

## Code cells

+++ {"editable": false, "deletable": false}

Now we will run a simple program in parallel.
The program is cut into several pieces, each one written in a distinct `%%code_block` cell.

Each `%%code_block` cell takes the following mandatory command line arguments:
- `-p` (`--prefix`) `<prefix>`: prefix for the temporary file in which the piece of code of the cell is written (*mandatory*)
- `-i` (`--id`) `<id>`: code cell identifier (this is used to merge the code cells in the right order)

**We strongly advise you not to mess with these arguments, otherwise you may experience trouble when trying to run your code!**

+++ {"editable": false, "deletable": false}

In this first cell, we initialize MPI.

```{code-cell}
%%code_block -p exercise_0 -i 1
program exercise_0

  implicit none
  include "mpif.h"

  integer :: i_rank, err

  call mpi_init(err)

```

+++ {"editable": false, "deletable": false}

In this second cell, each MPI rank will print a message.

```{code-cell}
%%code_block -p exercise_0 -i 2
  call mpi_comm_rank(MPI_COMM_WORLD, i_rank, err)

  print "(a15,1x,i0,a1)", "Hello from rank", i_rank, "!"

```

+++ {"editable": false, "deletable": false}

In this last cell, we finalize MPI.

```{code-cell}
%%code_block -p exercise_0 -i 3
  call mpi_finalize(err)

end program exercise_0
```



## Merge the code pieces, compile and run

+++ {"editable": false, "deletable": false}

Once all the code pieces have been written, it is time to run the program.
The `%merge_code_blocks` line will merge the code pieces, compile and execute the program with `mpirun`.

The line takes the following command line arguments:
- `-l` (`--language`) `<language>`: programming language (c, fortran or python) (*mandatory*)
- `-p` (`--prefix`) `<prefix>`: prefix for the temporary files (*mandatory*)
- `-n` (`--n_rank`) `<n_rank>`: number of MPI ranks (*optional*, the program in run only if a strictly positive number is provided)
- `-v` (`--verbose`): print the merged source code
- `-c` (`--clear`): delete all generated files once execution is complete


```{code-cell}
%merge_code_blocks -l fortran -p exercise_0 -n 2 -c -v
```

# Exercise 1

You can now move on to [Exercise 1](./../02_Exercise_1/exercise_1.ipynb).
