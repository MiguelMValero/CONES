!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2011-2017  ONERA
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
!
!********************************************************************************
!
! Set output logical unit
!
! parameters
!    output_listing      <-- Output listing file
!
!********************************************************************************
!
subroutine cwipi_set_output_logical_unit (outputUnit)

  use cwipi_printfort

  implicit none

  integer :: outputUnit
  
  ifile = outputUnit

  call cwipi_set_output_listing_f(outputUnit)

end subroutine cwipi_set_output_logical_unit
