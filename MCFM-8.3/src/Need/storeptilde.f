      subroutine storeptilde(nd,p)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'

      integer, intent(in) :: nd
      real(dp), intent(in) :: p(mxpart,4)
      integer:: i,j

      ptilde(nd,:,:) = p(:,:)
      
      end
      
