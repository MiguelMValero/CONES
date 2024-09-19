program exercise_0_code_2

  implicit none
  include "mpif.h"

  integer :: err

  call mpi_init(err)

  call mpi_finalize(err)

end program exercise_0_code_2
