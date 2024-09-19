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

module pdm_fortran

  use pdm

  implicit none

  integer, parameter :: PDM_TYPE_INT    = 0
#ifdef PDM_LONG_G_NUM
  integer, parameter :: PDM_TYPE_G_NUM  = 1
#else
  integer, parameter :: PDM_TYPE_G_NUM  = 0
#endif
  integer, parameter :: PDM_TYPE_DOUBLE   = 2
  integer, parameter :: PDM_TYPE_COMPLEX8 = 3
  integer, parameter :: PDM_TYPE_COMPLEX4 = 4
  integer, parameter :: PDM_TYPE_REAL4    = 5
  integer, parameter :: PDM_TYPE_CPTR     = 6

  interface

    subroutine pdm_fortran_free_c (ptrC) &
      bind (c, name = 'free')

      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptrC


    end subroutine pdm_fortran_free_c

  end interface

end module pdm_fortran
