!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2011  ONERA
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

module cwp_printfort

  implicit none

  integer, save :: ifile

  interface

    !> cwp_set_output_listing_cf bind with C

    subroutine cwp_set_output_listing_cf() &

        bind(c, name = 'cwp_set_output_listing_cf')

        use, intrinsic :: iso_c_binding
        implicit none

    end subroutine cwp_set_output_listing_cf

  end interface

contains

  !>
  !! \brief Set up the file used for the output listing.
  !!
  !!  \param [in]  outputUnit        Fortran output unit
  !!
  !!

  subroutine cwp_set_output_listing_f (outputUnit)

      use, intrinsic :: iso_c_binding
      implicit none

      integer :: outputUnit

      ifile =  outputUnit

      call cwp_set_output_listing_cf

  end subroutine cwp_set_output_listing_f

end module cwp_printfort

!>
!! \brief Set up the file used for the output listing.
!!
!!  \param [in]  chaine        String from C function
!!  \param [in]  taille        String size
!!
!!

subroutine printfortran (chaine, taille) &

  bind(c, name = 'printfortran')

  use, intrinsic :: iso_c_binding
  use cwp_printfort
  implicit none

  character (kind=c_char, len=1), dimension(*), intent (in) :: chaine
  integer (c_int), intent (inout) :: taille
  character(len = 16384) :: chloc
  integer       ii

  taille = min(taille, 16384 - 1)
  !
  do ii = 1, taille
     chloc(ii:ii) = chaine(ii)
  enddo
  !
  write(ifile, 1000, advance='no') chloc(1:taille)
  !
  return

  1000 format(a)

end subroutine printfortran

subroutine flushfortran () &

  bind(c, name = 'flushfortran')

  use, intrinsic :: iso_c_binding
  use cwp_printfort
  implicit none

  call flush(ifile)

end subroutine flushfortran
