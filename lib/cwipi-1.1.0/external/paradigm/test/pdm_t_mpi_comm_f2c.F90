!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2020  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

program testf

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_gnum
  use iso_c_binding



  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  integer        :: f_comm = MPI_COMM_WORLD
  integer        :: c_comm
  type(c_ptr)    :: gen_gnum
  integer        :: code


  call mpi_init(code)

  c_comm = PDM_MPI_Comm_f2c(f_comm)

  print *, f_comm, c_comm


  call PDM_gnum_create (gen_gnum, &
                        3,        &
                        1,        &
                        0,        &
                        1.d0,     &
                        f_comm,   &
                        0)

  call PDM_gnum_free(gen_gnum)

  call mpi_finalize(code)

end program testf
