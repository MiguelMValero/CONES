#include "pdm_configf.h"

program testf

use pdm
use iso_c_binding

implicit none

type(c_ptr) :: c = C_NULL_PTR


! integer, pointer :: f(:) => null()
! integer, pointer :: g(:) => null()
integer, allocatable, target :: f(:)
integer, allocatable :: g(:)
integer, pointer     :: p => null()

allocate(f(3))
f(1) = 1
f(2) = 2
f(3) = 3

allocate(g(2))
g(1) = 4
g(2) = 5


! print *, c_associated(c), c_associated(c, c_loc(f)), c_associated(c, c_loc(g))

p => f(1)
c = c_loc(p)
! c = c_loc(f)

! print *, c_associated(c), c_associated(c, c_loc(f)), c_associated(c, c_loc(g))

deallocate(f)

! print *, c_associated(c), c_associated(c, c_loc(f)), c_associated(c, c_loc(g))

! deallocate(g)

! call c_f_pointer(c, g, [3])

! print *, c_associated(c), c_associated(c, c_loc(f)), c_associated(c, c_loc(g))
! print *, g(:)

! deallocate(g)

end program testf
