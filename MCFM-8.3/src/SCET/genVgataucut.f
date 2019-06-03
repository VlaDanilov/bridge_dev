      subroutine genVgataucut(r,p,wt,*)
      implicit none
      include 'types.f'
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p34(->p3+p4)+p5+p6
c----
c---- with p6 generated using tau as a variable of integration,
c---- with minimum value taucut
      
      include 'constants.f'
      include 'mxpart.f'
      include 'limits.f'
      include 'vegas_common.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'x1x2.f'
      include 'taucut.f'
      include 'energy.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'debug.f'
      include 'cutoff.f'
      include 'kprocess.f'
      real(dp):: r(mxdim),p(mxpart,4),ph(mxpart,4),Q(4),p3(4),p4(4),p5(4),p34(4)
      real(dp):: wbw,wt,wt34,wt345,Qsq,rtshat,Qsqmin,Qsqmax,t0
      real(dp), parameter:: wt0=one/twopi**2

      p(:,:)=zip

      wt=zip

      Qsqmin=max(wsqmin,one) ! ensure minimum value of m(345)>1 GeV
      Qsqmax=sqrts**2*0.999d0 ! ensure maximum value a bit below s
c--- generate invariant mass of Q=p3+p4+p5
      wbw=one
      call pick(2,Qsq,Qsqmin,Qsqmax,r(1),wbw)
      
      rtshat=sqrt(Qsq)

c--- in the tauboost case, we do not know what value of tau we will eventually end up with,
c--- so just choose it very small;  it is still beneficial to use this routine since it
c--- provides better sampling than the alternatives in the region of very small tau
      if (tauboost) then
        t0=1.e-12_dp
      else
        t0=taucut
      endif

c--- generate pa+pb -> Q
      call genQ(rtshat,r(2),t0,p,wt)
            
c--- generate extra parton
c--- note that r(ndim+1) is a uniform random variable not adapted by VEGAS
      call genparton(2,p,r(3),r(4),r(5),r(ndim+1),t0,ph,wt)

      Q(:)=ph(3,:)
      if (kcase == kWgajet) then
        call phi1_2m_bw(zip,r(6),r(7),r(8),wsqmin,Q,p5,p34,wmass,wwidth,wt345,*999)
      else
        call phi1_2m_bw(zip,r(6),r(7),r(8),wsqmin,Q,p5,p34,zmass,zwidth,wt345,*999)
      endif
      call phi3m0(r(9),r(10),p34,p3,p4,wt34,*999)

c--- translate to momenta to be returned
      p(1,:)=-ph(1,:)
      p(2,:)=-ph(2,:)
      p(3,:)=p3(:)
      p(4,:)=p4(:)
      p(5,:)=p5(:)
      p(6,:)=ph(4,:)
      p(7,:)=zip
      
      xx(1)=-two*p(1,4)/sqrts
      xx(2)=-two*p(2,4)/sqrts
      
      wt=wt0*wbw*wt*wt345*wt34
c--- the factor below will be added later, in lowint
      wt=wt*xx(1)*xx(2)*sqrts**2
      
c      call writeout(p)
c      pause

! trap pathological points      
      if (p(3,4) .ne. p(3,4)) then
!        write(6,*) 'p(3,4) NaN'
        return 1
      endif
      
      if   ((xx(1) > 1._dp) 
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
        if (debug) write(6,*) 'problems with xx(1),xx(2) in gen3taucut',
     & xx(1),xx(2)  
        return 1 
      endif
          
      return

  999 wt=zip
      return 1

      end
