      subroutine donothing_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k,in
      real(dp), intent(in) :: p(mxpart,4), n(4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

      msq(:,:) = 0._dp
      end


