      function cln(x,sgn)
        implicit none
        include 'types.f'
        include 'constants.f'

        complex(dp) :: cln
        real(dp), intent(in) :: x, sgn

        cln = log(abs(x))
        if (x < 0._dp) then
          cln = cln + sign(1._dp,sgn)*im*pi
        endif

      end function
