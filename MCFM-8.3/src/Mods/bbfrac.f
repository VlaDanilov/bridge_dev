      module bbfrac_m
        use types
        implicit none

        private
        include 'maxd.f'

        ! 40 = maxd
        real(dp), save, public :: bbfrac(0:maxd)
!$omp threadprivate(bbfrac)

      end module
