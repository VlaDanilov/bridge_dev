! Simple derivative by Taylor expansion about midpoint:
!    f'(x)= 1/(2a)*(f(a+x)-f(x-a)) + O(f''')
      subroutine midpoint_devpdf(ih1,ibeam,x,pdf)
      implicit none
      include 'types.f'
      include 'scale.f'
      include 'facscale.f'
      real(dp) :: x,pdf(-5:5),upx,downx
      real(dp) :: uppdfs(-5:5),downpdfs(-5:5)
      real(dp) :: smallxfac,tinyx
      integer ih1,ibeam,j

      smallxfac=1e-4_dp

      tinyx=smallxfac*x
      upx=x+tinyx
      downx=x-tinyx
      
      call fdist(ih1,upx,facscale,uppdfs)
      call fdist(ih1,downx,facscale,downpdfs)

      pdf(:)=0._dp
      do j=-5,5
         pdf(j)=(uppdfs(j)-downpdfs(j))/(2._dp*tinyx)
      enddo

      return
      end
