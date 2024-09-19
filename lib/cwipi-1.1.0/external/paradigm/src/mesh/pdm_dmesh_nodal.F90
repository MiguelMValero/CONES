!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2023  ONERA
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

module pdm_dmesh_nodal

  use pdm

  implicit none

  interface

    !>
    !!
    !! \brief Setup global distribution of all elements register in current structure
    !!
    !! \param [inout]   dmn   Pointer to \ref PDM_dmesh_nodal_t object
    !!

    subroutine PDM_dmesh_nodal_generate_distribution (dmn) &
    bind(c, name='PDM_dmesh_nodal_generate_distribution')

      use iso_c_binding
      implicit none

      type (c_ptr), value :: dmn
    end subroutine PDM_dmesh_nodal_generate_distribution

    !>
    !!
    !! \brief Free \ref PDM_dmesh_nodal_t object
    !!
    !! \param [inout]   dmn   Pointer to \ref PDM_dmesh_nodal_t object
    !!

    subroutine PDM_DMesh_nodal_free (dmn) &
    bind(c, name='PDM_DMesh_nodal_free')

      use iso_c_binding
      implicit none

      type (c_ptr), value :: dmn
    end subroutine PDM_DMesh_nodal_free

  end interface

end module pdm_dmesh_nodal
