      subroutine phase_realstop(r,p,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      !include 'limits.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'breit.f'
      include 'x1x2.f'
      include 'phasemin.f'
      include 'energy.f'
      include 'limits.f'

      real(dp), intent(in) :: r(mxdim)
      real(dp), intent(out) :: p(mxpart,4), wt

      real(dp) :: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp) :: p12(4),p126(4),p34(4),p345(4)
      real(dp) :: wt0,wt1,wt2,wt3,wt4
      real(dp) :: smin

      real(dp) :: tau,y,xjac
      
      parameter(wt0=1._dp/twopi**2)

      smin = 0._dp

      tau=exp(log(taumin)*r(12))
      y=0.5_dp*log(tau)*(1._dp-2._dp*r(13))
      xjac=log(taumin)*tau*log(tau)

      xx(1)=sqrt(tau)*exp(+y)
      xx(2)=sqrt(tau)*exp(-y)

      if   ((xx(1) > 1._dp) 
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half

      p12 = -p1-p2
      call phi1_2m(0._dp,r(1),r(2),r(3),smin,p12,p6,p126,wt1,*99)
      call phi1_2m_bw(0._dp,r(4),r(5),r(6),smin,p126,p7,p345,mt,twidth,wt2,*99)
      call phi1_2m_bw(0._dp,r(7),r(8),r(9),smin,p345,p5,p34,wmass,wwidth,wt3,*99)
      call phi3m0(r(10),r(11),p34,p3,p4,wt4,*99)

      wt = xjac*wt0*wt1*wt2*wt3*wt4

      p = 0._dp
      p(1,:) = p1
      p(2,:) = p2
      p(3,:) = p3
      p(4,:) = p4
      p(5,:) = p5
      p(6,:) = p6
      p(7,:) = p7

      return

 99   wt=0._dp
      return 1
      end

