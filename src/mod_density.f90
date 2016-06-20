#include "real.h"
module mod_density
use mod_grid
use mod_params, only : m_rho, m_K
implicit none

! Parameters.
real(kind=8), private, parameter :: g = 9.81d0

contains

   subroutine density_initialize(dg,nlon,nlat)
      type(data_grids), target, intent(inout) :: dg 
      integer(kind=4), intent(in) :: nlon, nlat

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, dz, m_rhoC
      integer(kind=4) :: i, j

      hz     => dg%wave_field%hz
      dz     => dg%depth_field%dz
      m_rhoC => dg%wave_field%m_rhoC

!$omp parallel do private(i)
      do j = 1, nlat
         do i = 1, nlon
            m_rhoC(i,j) = (1.0d0 + m_rho*dz(i,j)*g/(2.0d0*m_K)) &
                        / (1.0d0 + m_rho*dz(i,j)*g/m_K)
         end do
      end do

      return
   end subroutine density_initialize

end module mod_density
