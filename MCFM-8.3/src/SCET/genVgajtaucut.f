      subroutine genVgajtaucut(r,p,wt,*)
      implicit none
      include 'types.f'
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p34(->p3+p4)+p5+p6+p7
c----
c---- with p6 and p7 generated using tau as a variable of integration,
c---- with minimum value t0
      
      include 'constants.f'
      include 'mxpart.f'
      include 'limits.f'
      include 'vegas_common.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'debug.f'
      include 'kprocess.f'
      real(dp):: r(mxdim),p(mxpart,4),ph(mxpart,4),Q(4),p3(4),p4(4),p5(4),p34(4)
      real(dp):: wbw,wt,wt34,wt345,Qsq,rtshat,mass,t0,Qsqmin,Qsqmax
      real(dp), parameter:: wt0=one/twopi**2
      logical pass

      p(:,:)=zip

      wt=zip

      Qsqmin=max(wsqmin,one) ! ensure minimum value of m(345)>1 GeV
      Qsqmax=sqrts**2*0.999d0 ! ensure maximum value a bit below s
c--- generate invariant mass of Q=p3+p4+p5
      wbw=one
      call pick(2,Qsq,Qsqmin,Qsqmax,r(1),wbw)
      
      rtshat=sqrt(Qsq)

      t0=1.e-12_dp

c--- generate pa+pb -> Q
      call genQ(rtshat,r(2),2._dp*t0,p,wt)

c--- generate extra partons
c--- note that r(ndim+1) is a uniform random variable not adapted by VEGAS
      call genparton2(2,p,r(3),r(4),r(5),r(6),r(7),r(8),
     &                r(ndim+1),t0,ph,wt,pass)

      if (pass .eqv. .false.) then
        goto 999
      endif

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
      p(7,:)=ph(5,:)
      p(8,:)=zip

      xx(1)=-two*p(1,4)/sqrts
      xx(2)=-two*p(2,4)/sqrts

      wt=wt0*wbw*wt*wt345*wt34
c--- the factor below will be added later, in realint
      wt=wt*xx(1)*xx(2)*sqrts**2
      
c      call writeout(p)
c      pause
      
      if   ((xx(1) > 1._dp) 
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
        if (debug) write(6,*) 'problems with xx(1),xx(2) in gen4taucut',
     &  xx(1),xx(2)  
        return 1 
      endif
          
          
      if (p(3,4) == p(3,4)) then
        continue
      else
        write(6,*) 'warning: discarding point in gen4taucut'
c        write(6,*) '      vector(1)=',r(1)
c        write(6,*) '      vector(2)=',r(2)
c        write(6,*) '      vector(3)=',r(3)
c        write(6,*) '      vector(4)=',r(4)
c        write(6,*) '      vector(5)=',r(5)
c        write(6,*) '      vector(6)=',r(6)
c        write(6,*) '      vector(7)=',r(7)
c        write(6,*) '      vector(8)=',r(8)
c        write(6,*) '      vector(9)=',r(9)
c        write(6,*) '      vector(10)=',r(10)
        goto 999
      endif
          
      return

  999 wt=zip
      return 1

      end
