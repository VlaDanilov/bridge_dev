      subroutine wrecons(pl_in, pnu_in, pnu_out)
          use types
          implicit none
          include 'masses.f'

          real(dp), intent(in) :: pl_in(4)
          real(dp), intent(in) :: pnu_in(4)
          real(dp), intent(out) :: pnu_out(4)

          integer, parameter :: ielectron = 1
          integer, parameter :: ineutrino = 2

          real(dp) :: ps(0:3, 2)
          real(dp) :: pn(0:3)

          real(dp) :: a,b,c,d, fix
          logical ::  good

          integer :: i,j

          ps(0:3,1) = [pl_in(4), pl_in(1), pl_in(2), pl_in(3)]
          ps(0:3,2) = [pnu_in(4), pnu_in(1), pnu_in(2), pnu_in(3)]

          pn(:) = ps(:,2)


          ! from MG 1

          i = ielectron
          j = ineutrino
          a=wmass*wmass*.5d0+ps(1,i)*ps(1,j)+ps(2,i)*ps(2,j)
          b = ps(0,i)
          c = dsqrt(ps(1,j)**2+ps(2,j)**2)
          d = ps(3,i)
          good =(a**2 - b**2*c**2 + c**2*d**2 .gt. 0d0)
          fix=1d0
          do while (.not. good)
             fix=fix*0.9d0
             ps(1,j)=ps(1,j)*.9d0
             ps(2,j)=ps(2,j)*.9d0
             a = wmass*wmass*.5d0+ps(1,i)*ps(1,j)+ps(2,i)*ps(2,j)
             b = ps(0,i)
             c = dsqrt(ps(1,j)**2+ps(2,j)**2)
             d = ps(3,i)
             good =(a**2 - b**2*c**2 + c**2*d**2 .gt. 0d0)
          enddo
          ps(3,j)  = (2d0*a*d/(b**2 - d**2) + 
     &        2d0*b*Sqrt(max(a**2 - b**2*c**2 + c**2*d**2,0d0))/
     &     ((b - d)*(b + d)))/2d0
          
          pn(3)   = (2d0*a*d/(b**2 - d**2) - 
     &        2d0*b*Sqrt(max(a**2 - b**2*c**2 + c**2*d**2,0d0))/
     &     ((b - d)*(b + d)))/2d0
          if (abs(ps(3,j)) .gt. abs(pn(3))) then
             a=ps(3,j)
             ps(3,j)=pn(3)
             pn(3)=a
          endif
          ps(0,j) = dsqrt(ps(1,j)**2+ps(2,j)**2+ps(3,j)**2)
          pn(0) = dsqrt(pn(1)**2+pn(2)**2+pn(3)**2)

          pnu_out(4) = ps(0,j)
          pnu_out(1:3) = ps(1:3,j)


      end subroutine
