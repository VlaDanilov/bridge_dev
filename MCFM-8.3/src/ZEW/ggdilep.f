      subroutine ggdilep(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      real(dp) :: p(mxpart,4),msq(fn:nf,fn:nf),t1,t2,fac

      call dotem(3,p,s)

      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)

      msq(:,:)=zip
!      esq=fourpi/137.03599911_dp
      fac=esq**2*q1**4/four !four is spin-ave
      msq(0,0) = fac*eight*(t1/t2+t2/t1)

      return
      end
