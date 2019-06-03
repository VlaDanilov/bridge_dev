      module differentiation_m
        use, intrinsic :: iso_fortran_env
        implicit none

        private
        public :: adaptive_deriv_forward

        include 'types.f'

        abstract interface
          function evalf(x, params)
            include 'types.f'
            real(dp), intent(in) :: x
            real(dp), intent(in) :: params(:)
            real(dp) :: evalf
          end function
        end interface

        type, public :: diff_function
          real(dp), allocatable :: params(:)
          procedure(evalf), pointer, nopass :: funPtr
        end type

        contains

        ! based on libgsl deriv/deriv.c,
        ! see there for further explanations of formulas and estimations
        subroutine deriv_forward(func, x, h, ret, abserr_round, abserr_trunc)
          implicit none
          include 'types.f'

          type(diff_function), intent(in) :: func
          real(dp), intent(in) :: x
          real(dp), intent(in) :: h
          real(dp), intent(out) :: ret
          real(dp), intent(out) :: abserr_round, abserr_trunc

          real(dp) :: f1,f2,f3,f4, r2,r4
          real(dp) :: eps
          real(dp) :: e4, dy

          eps = epsilon(1._dp)

          f1 = func%funPtr(x + h / 4._dp, func%params)
          f2 = func%funPtr(x + h / 2._dp, func%params)
          f3 = func%funPtr(x + (3._dp/4._dp)*h, func%params)
          f4 = func%funPtr(x + h, func%params)

          r2 = 2._dp*(f4 - f2)
          r4 = (22._dp / 3._dp) * (f4 - f3) - (62._dp / 3._dp) * (f3 - f2) +
     &          (52._dp / 3._dp) * (f2 - f1)

          ! rounding error for r4
          e4 = 2._dp * 20.67_dp * (abs(f4) + abs(f3) + abs(f2) + abs(f1)) * eps

          ! due to finite precision in x+h = O(eps*x)
          dy = max(abs(r2/h), abs(r4/h)) * abs(x/h) *  eps

          ret = r4 / h
          abserr_trunc = abs((r4-r2)/h)
          abserr_round = abs(e4/h) + dy

        end subroutine

        subroutine adaptive_deriv_forward(func, x, h, ret, sumerr)
          implicit none
          include 'types.f'

          type(diff_function), intent(in) :: func
          real(dp), intent(in) :: x
          real(dp), intent(in) :: h
          real(dp), intent(out) :: ret, sumerr

          real(dp) :: r_0, err_round, err_trunc, err

          call deriv_forward(func, x, h, r_0, err_round, err_trunc)
          err = err_round + err_trunc

          if (err_round < err_trunc .and. (err_round > 0._dp .and. err_trunc > 0._dp)) then
            block
              real(dp) :: r_opt, err_round_opt, err_trunc_opt, err_opt
              real(dp) :: h_opt

              h_opt = h * sqrt(err_round/err_trunc)
              call deriv_forward(func, x, h_opt, r_opt, err_round_opt, err_trunc_opt)
              err_opt = err_round_opt + err_trunc_opt

              if (err_opt < err .and. abs(r_opt - r_0) < 4._dp * err) then
                r_0 = r_opt
                err = err_opt
              endif

            end block
          endif

          ret = r_0
          sumerr = err

        end subroutine

      end module

