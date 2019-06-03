      subroutine qqbggHampsq(i1,i2,i3,i4,p,msq,totsq,ABsq,BAsq,SLsq)
      implicit none
!--   This routine calculates the amplitude squared                                      
!--   for the process q^+(1) qb^-(2) g(3) g(4) H(p1234)
!--   with all momenta outgoing.    
!--   Amplitude squared, summed over colours and helicities
!--   is returned with an overall factor of g^8/16/pi^4/v^2 removed
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'scale.f'
      complex(dp)::Amp(2,2,2,2),AB(2,2,2),BA(2,2,2),sum(2,2,2)
      integer::i1,i2,i3,i4,h12,h3,h4
      real(dp):: p(mxpart,4),msq,totsq,ABsq,BAsq,SLsq
      logical::flip
      logical,save::first=.true.

      call spinoru(4,p,za,zb)

!     1st index Color ordering
!     2nd index helicity line 1
!     3rd index helicity gluon 3
!     4rd index helicity gluon 4
      call qqbggHamp(i1,i2,i3,i4,flip,p,msq,Amp)

      sum(:,:,:)=Amp(1,:,:,:)+Amp(2,:,:,:)

      ABsq=zero
      BAsq=zero
      SLsq=zero
      do h12=1,2
      do h3=1,2
      do h4=1,2
      ABsq=ABsq+cf*xn**2/2d0*Amp(1,h12,h3,h4)*Conjg(Amp(1,h12,h3,h4))
      BAsq=BAsq+cf*xn**2/2d0*Amp(2,h12,h3,h4)*Conjg(Amp(2,h12,h3,h4))
      SLsq=SLsq-cf/2d0*sum(h12,h3,h4)*Conjg(sum(h12,h3,h4))
      enddo
      enddo
      enddo
      totsq=ABsq+BAsq+SLsq
      return
      end

