      module anomcoup_tbW
          use types
          implicit none

          complex(dp), public, save, protected :: anomc1, anomc2, anomc3, anomc4, anomc6, anomc7, anomc8, anomc9

          logical, public, save :: disable_sm = .false.
          logical, public, save :: enable_lambda2 = .false.
          logical, public, save :: enable_lambda4 = .false.

          logical, public, save :: mode_anomcoup = .false.

          real(dp), public, save, protected :: lambda
          public :: anomcoup_tbW_set_lambda

          complex(dp), private, save :: priv_c1, priv_c2, priv_c3, priv_c4, priv_c6, priv_c7, priv_c8, priv_c9

          ! QCD gauge invariant subsets, only for debugging
          logical, public, parameter :: enable_resonant = .true.
          logical, public, parameter :: enable_nonresonant = .true.

      contains

        subroutine anomcoup_init()
          implicit none
          include 'ewcouple.f'
          include 'masses.f'
          include 'qcdcouple.f'
          include 'nproc.f'

          ! these mt are real valued normalization factors
          ! while in the amplitudes the complex valued mt as part of the complex mass scheme is used

          ! note that c3,c4,c6,c7 have a minus sign added here, this is the correct normalization
          ! but the studies by Zhang et al. and Protos use a different sign (i.e. no minus sign here)

          anomc1 = priv_c1 * mt**2 ! phi_q (V_L)
          anomc2 = priv_c2 * mt**2 ! phi_phi (V_R)
          anomc3 = -priv_c3 * mt ! uW (g_R)
          anomc4 = -priv_c4 * mt ! dW (g_L)
          anomc6 = -priv_c6 * mt ! t_G
          anomc7 = -priv_c7 * mt ! b_G
          anomc8 = priv_c8 / gw**2 ! 4L      ! undo normalization with gw^2
          anomc9 = priv_c9 / gw**2 ! 4R      ! undo normalization with gw^2

        end subroutine

        subroutine anomcoup_cc()
            implicit none

            anomc2 = conjg(anomc2) 
            anomc3 = conjg(anomc3) 
            anomc4 = conjg(anomc4) 
            anomc6 = conjg(anomc6) 
            anomc7 = conjg(anomc7) 

        end subroutine

        subroutine anomcoup_tbW_set_lambda(l)
            implicit none
            real(dp), intent(in) :: l

            lambda = l
        end subroutine

        subroutine anomcoup_tbW_set_c1(c1)
            implicit none
            complex(dp), intent(in) :: c1

            priv_c1 = c1
            if (c1 /= (0._dp,0._dp)) enable_lambda2 = .true.
        end subroutine

        subroutine anomcoup_tbW_set_c2(c2)
            implicit none
            complex(dp), intent(in) :: c2

            priv_c2 = c2
            if (c2 /= (0._dp,0._dp)) enable_lambda2 = .true.
        end subroutine

        subroutine anomcoup_tbW_set_c3(c3)
            implicit none
            complex(dp), intent(in) :: c3

            priv_c3 = c3
            if (c3 /= (0._dp,0._dp)) enable_lambda2 = .true.
        end subroutine

        subroutine anomcoup_tbW_set_c4(c4)
            implicit none
            complex(dp), intent(in) :: c4

            priv_c4 = c4
            if (c4 /= (0._dp,0._dp)) enable_lambda2 = .true.
        end subroutine

        subroutine anomcoup_tbW_set_c6(c6)
            implicit none
            complex(dp), intent(in) :: c6

            priv_c6 = c6
            if (c6 /= (0._dp,0._dp)) enable_lambda2 = .true.
        end subroutine

        subroutine anomcoup_tbW_set_c7(c7)
            implicit none
            complex(dp), intent(in) :: c7

            priv_c7 = c7
            if (c7 /= (0._dp,0._dp)) enable_lambda2 = .true.
        end subroutine

        subroutine anomcoup_tbW_set_c8(c8)
            implicit none
            complex(dp), intent(in) :: c8

            priv_c8 = c8
            if (c8 /= (0._dp,0._dp)) enable_lambda2 = .true.
        end subroutine

        subroutine anomcoup_tbW_set_c9(c9)
            implicit none
            complex(dp), intent(in) :: c9

            priv_c9 = c9
            if (c9 /= (0._dp,0._dp)) enable_lambda2 = .true.
        end subroutine

      end module

