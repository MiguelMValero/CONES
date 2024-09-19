module pdm
#include "pdmf.h"

integer, parameter :: PDM_FALSE = 0
integer, parameter :: PDM_TRUE  = 1

interface

function PDM_MPI_Comm_f2c (f_comm) &
result (c_comm)                    &
bind (c, name = 'PDM_MPI_Comm_f2c')

  use iso_c_binding
  implicit none

  integer(c_int), value :: f_comm
  integer(c_int)        :: c_comm

end function PDM_MPI_Comm_f2c



function PDM_MPI_Comm_c2f (c_comm) &
result (f_comm)                    &
bind (c, name = 'PDM_MPI_Comm_c2f')

  use iso_c_binding
  implicit none

  integer(c_int), value :: c_comm
  integer(c_int)        :: f_comm

end function PDM_MPI_Comm_c2f

end interface

end module pdm
