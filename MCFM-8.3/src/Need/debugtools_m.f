      module debugtools_m
        use, intrinsic :: iso_fortran_env
        implicit none

      contains
      
      ! define SR in Mathematica as:
      ! SR[str_] := ToExpression@StringReplace[str, {"E+" -> "*^", "E-" -> "*^-"}];
      subroutine sam_declarespinor(p,n)
        implicit none
        include 'types.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: n
        integer :: j

        do j=1,n
          write (error_unit,'(A,I1,A,ES21.14,A,ES21.14,A,ES21.14,A,ES21.14,A)') 
     &     "DeclareSpinorMomentum[p",j,',{ SR["',
     &     p(j,4), '"], SR["', p(j,1), '"], SR["', p(j,2),
     &       '"], SR["',  p(j,3), '"] }]' 
        enddo

      end subroutine

      subroutine sam_zazb(za,zb,n)
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        integer, intent(in) :: n

        integer :: j1,j2

        do j1=1,n
          do j2=1,n
            write (error_unit,"(A,I1,A,I1,A,ES21.14,A,ES21.14,A)") "za[p",j1,
     &                   ", p",j2,'] = SR["',real(za(j1,j2)),
     &                   " + I*",aimag(za(j1,j2)),'"]'
            write (error_unit,"(A,I1,A,I1,A,ES21.14,A,ES21.14,A)") "zb[p",j1,
     &                   ", p",j2,'] = SR["',real(zb(j1,j2)),
     &                   " + I*",aimag(zb(j1,j2)),'"]'
          enddo
        enddo
      end subroutine

      ! for now we set 2 for + helicity, 1 for - helicity
      ! i1, i3 are the color ordering adjacent fields to the gluon ig
      subroutine check_softfact(amp_single, amp_double, za,zb, i1, ig, i3, hg)
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp), intent(in) :: amp_single, amp_double
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        integer, intent(in) :: i1,ig,i3,hg

        complex(dp) :: soft

        if (hg == 2) then
          soft = za(i1,i3)/za(i1,ig)/za(ig,i3)
        else
          soft = -zb(i1,i3)/zb(i1,ig)/zb(ig,i3)
        endif

        write (error_unit,*) "Soft factorization check"
        write (error_unit,*) "soft*single/double ="//achar(27)//"[31m ", 
     &            amp_single*soft/amp_double, achar(27)//"[0m"
        write (error_unit,*) ""
      end subroutine

      ! again hg1, hg2 are 1 for - helicity, 2 for +
      ! label for all incoming
      subroutine check_collfact_gg(amp_single_plus, amp_single_minus,
     &      amp_double, za, zb, z, ig1, ig2, hg1, hg2 )
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp), intent(in) :: amp_single_plus, amp_single_minus
        complex(dp), intent(in) :: amp_double
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        real(dp), intent(in) :: z
        integer, intent(in) :: ig1, ig2
        integer, intent(in) :: hg1, hg2

        complex(dp) :: split_plus, split_minus

        write (error_unit,*) "gg collinear factorization check"

        ! sum over two intermediate particle helicities

        if (hg1 == 2 .and. hg2 == 2) then
          split_plus = 0._dp
          split_minus = 1/sqrt(z*(1-z))/za(ig1,ig2)
        elseif (hg1 == 1 .and. hg2 == 1) then
          split_plus = -1/sqrt(z*(1-z))/zb(ig1,ig2)
          split_minus = 0._dp
        elseif (hg1 == 2 .and. hg2 == 1) then
          split_plus = (1-z)**2/sqrt(z*(1-z))/za(ig1,ig2)
          split_minus = -z**2/sqrt(z*(1-z))/zb(ig1,ig2)
        elseif (hg1 == 1 .and. hg2 == 2) then
          split_plus = z**2/sqrt(z*(1-z))/za(ig1,ig2)
          split_minus = -(1-z)**2/sqrt(z*(1-z))/zb(ig1,ig2)
        else
          write (error_unit,*) "check_collfact_gg: unknown helicity config"
          call abort
        endif

        write (error_unit,*) "split*single/double = "//achar(27)//"[31m ",
     &      (split_plus*amp_single_minus + split_minus*amp_single_plus)/
     &       amp_double, achar(27)//"[0m"
        write (error_unit,*) ""

      end subroutine

      subroutine check_collfact_qq(amp_single_plus, amp_single_minus,
     &      amp_double, za, zb, z, iq1, iq2, hq1, hq2)
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp), intent(in) :: amp_single_plus, amp_single_minus
        complex(dp), intent(in) :: amp_double
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        real(dp), intent(in) :: z
        integer, intent(in) :: iq1, iq2
        integer, intent(in) :: hq1, hq2

        complex(dp) :: split_plus, split_minus

        write (error_unit,*) "qq collinear factorization check"

        if (hq1 == 1 .and. hq2 == 2) then
          split_minus = -(1-z)/zb(iq1,iq2)
          split_plus = -z/za(iq1,iq2)
        elseif (hq1 == 2 .and. hq2 == 1) then
          split_minus = -z/zb(iq1,iq2)
          split_plus = -(1-z)/za(iq1,iq2)
        else
          write (error_unit,*) "check_collfact_qq: unknown helicity config"
          call abort
        endif

        write (error_unit,*) "split*single/double = "//achar(27)//"[31m ",
     &      (split_plus*amp_single_minus + split_minus*amp_single_plus)/
     &       amp_double, achar(27)//"[0m"
        write (error_unit,*) ""

      end subroutine

      end module
