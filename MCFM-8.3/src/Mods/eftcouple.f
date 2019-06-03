      ! initialize additional calculated parameters, maybe eft-modified
      module eftcouple
          use types
          implicit none

          public

          include 'ewcouple.f'
          include 'qcdcouple.f'
          include 'masses.f'

          real(dp), public, save, protected :: ecossin
          real(dp), public, save, protected :: eftgw, gb

          public :: eftcouple_init

          contains

          subroutine eftcouple_init()
              implicit none

              gb = sqrt(esq)/sqrt(1-xw)
              ecossin = sqrt(gw**2 + gb**2)

          end subroutine

      end module
