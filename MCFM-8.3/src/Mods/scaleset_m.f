      module scaleset_m
          use types
          implicit none

          public :: set_scales

          contains

          subroutine set_scales(renscale_in,facscale_in)
              use constants, only: pi
              implicit none

              include 'couple.f'
              include 'qcdcouple.f'
              include 'nlooprun.f'
              include 'scale.f'
              include 'facscale.f'

              real(dp), intent(in) :: renscale_in, facscale_in

              real(dp) :: alphas

              scale = renscale_in
              facscale = facscale_in

              as = alphas(scale,amz,nlooprun)

              ason2pi = as/2/pi
              ason4pi = as/4/pi
              gsq = 4*pi*as
              musq = scale**2
          end subroutine
      end module
