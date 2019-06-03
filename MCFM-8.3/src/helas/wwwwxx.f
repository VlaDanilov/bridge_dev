c
c ======================================================================
c
      subroutine wwwwxx(wm1,wp1,wm2,wp2,gwwa,gwwz,zmass,zwidth , vertex)
c
c this subroutine computes an amplitude of the four-point w-/w+         
c coupling, including the contributions of photon and z exchanges.  the 
c photon propagator is given in feynman gauge and the z propagator is   
c given in unitary gauge.                                               
c                                                                       
c input:                                                                
c       complex wm1(0:3)       : first  flow-out w-                  wm1
c       complex wp1(0:3)       : first  flow-out w+                  wp1
c       complex wm2(0:3)       : second flow-out w-                  wm2
c       complex wp2(0:3)       : second flow-out w+                  wp2
c       real    gwwa           : coupling constant of w and a       gwwa
c       real    gwwz           : coupling constant of w and z       gwwz
c       real    zmass          : mass  of z                             
c       real    zwidth         : width of z                             
c                                                                       
c output:                                                               
c       complex vertex         : amplitude        gamma(wm1,wp1,wm2,wp2)
c
      complex*16    wm1(6),wp1(6),wm2(6),wp2(6),vertex
      complex*16 dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),
     &           j12(0:3),j34(0:3),j14(0:3),j32(0:3),dvertx,
     &           sv1,sv2,sv3,sv4,tv1,tv2,tv3,tv4,dzs,dzt,
     &           v12,v13,v14,v23,v24,v34,js12,js34,js14,js32,js,jt
      real*8       pwm1(0:3),pwp1(0:3),pwm2(0:3),pwp2(0:3),
     &           gwwa,gwwz,zmass,zwidth
      real*8     q(0:3),k(0:3),dp1(0:3),dp2(0:3),dp3(0:3),dp4(0:3),
     &           dgwwa2,dgwwz2,dgw2,dmz,dwidth,s,t,das,dat
c
      real*8 r_zero, r_one, r_two
      parameter( r_zero=0.0d0, r_one=1.0d0, r_two=2.0d0 )
c
      pwm1(0)=dble( wm1(5))
      pwm1(1)=dble( wm1(6))
      pwm1(2)=dimag(wm1(6))
      pwm1(3)=dimag(wm1(5))
      pwp1(0)=dble( wp1(5))
      pwp1(1)=dble( wp1(6))
      pwp1(2)=dimag(wp1(6))
      pwp1(3)=dimag(wp1(5))
      pwm2(0)=dble( wm2(5))
      pwm2(1)=dble( wm2(6))
      pwm2(2)=dimag(wm2(6))
      pwm2(3)=dimag(wm2(5))
      pwp2(0)=dble( wp2(5))
      pwp2(1)=dble( wp2(6))
      pwp2(2)=dimag(wp2(6))
      pwp2(3)=dimag(wp2(5))
c
      dv1(0)=dcmplx(wm1(1))
      dv1(1)=dcmplx(wm1(2))
      dv1(2)=dcmplx(wm1(3))
      dv1(3)=dcmplx(wm1(4))
      dp1(0)=dble(pwm1(0))
      dp1(1)=dble(pwm1(1))
      dp1(2)=dble(pwm1(2))
      dp1(3)=dble(pwm1(3))
      dv2(0)=dcmplx(wp1(1))
      dv2(1)=dcmplx(wp1(2))
      dv2(2)=dcmplx(wp1(3))
      dv2(3)=dcmplx(wp1(4))
      dp2(0)=dble(pwp1(0))
      dp2(1)=dble(pwp1(1))
      dp2(2)=dble(pwp1(2))
      dp2(3)=dble(pwp1(3))
      dv3(0)=dcmplx(wm2(1))
      dv3(1)=dcmplx(wm2(2))
      dv3(2)=dcmplx(wm2(3))
      dv3(3)=dcmplx(wm2(4))
      dp3(0)=dble(pwm2(0))
      dp3(1)=dble(pwm2(1))
      dp3(2)=dble(pwm2(2))
      dp3(3)=dble(pwm2(3))
      dv4(0)=dcmplx(wp2(1))
      dv4(1)=dcmplx(wp2(2))
      dv4(2)=dcmplx(wp2(3))
      dv4(3)=dcmplx(wp2(4))
      dp4(0)=dble(pwp2(0))
      dp4(1)=dble(pwp2(1))
      dp4(2)=dble(pwp2(2))
      dp4(3)=dble(pwp2(3))
      dgwwa2=dble(gwwa)**2
      dgwwz2=dble(gwwz)**2
      dgw2  =dgwwa2+dgwwz2
      dmz   =dble(zmass)
      dwidth=dble(zwidth)
c
      v12= dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
      v13= dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
      v14= dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
      v23= dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
      v24= dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
      v34= dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)
c
      q(0)=dp1(0)+dp2(0)
      q(1)=dp1(1)+dp2(1)
      q(2)=dp1(2)+dp2(2)
      q(3)=dp1(3)+dp2(3)
      k(0)=dp1(0)+dp4(0)
      k(1)=dp1(1)+dp4(1)
      k(2)=dp1(2)+dp4(2)
      k(3)=dp1(3)+dp4(3)
c
      s=q(0)**2-q(1)**2-q(2)**2-q(3)**2
      t=k(0)**2-k(1)**2-k(2)**2-k(3)**2
c
      das=-r_one/s
      dat=-r_one/t
      dzs=-r_one/dcmplx( s-dmz**2 , dmax1(dsign(dmz*dwidth,s),r_zero) )
      dzt=-r_one/dcmplx( t-dmz**2 , dmax1(dsign(dmz*dwidth,t),r_zero) )
c
      sv1= (dp2(0)+q(0))*dv1(0) -(dp2(1)+q(1))*dv1(1)
     &    -(dp2(2)+q(2))*dv1(2) -(dp2(3)+q(3))*dv1(3)
      sv2=-(dp1(0)+q(0))*dv2(0) +(dp1(1)+q(1))*dv2(1)
     &    +(dp1(2)+q(2))*dv2(2) +(dp1(3)+q(3))*dv2(3)
      sv3= (dp4(0)-q(0))*dv3(0) -(dp4(1)-q(1))*dv3(1)
     &    -(dp4(2)-q(2))*dv3(2) -(dp4(3)-q(3))*dv3(3)
      sv4=-(dp3(0)-q(0))*dv4(0) +(dp3(1)-q(1))*dv4(1)
     &    +(dp3(2)-q(2))*dv4(2) +(dp3(3)-q(3))*dv4(3)
c
      tv1= (dp4(0)+k(0))*dv1(0) -(dp4(1)+k(1))*dv1(1)
     &    -(dp4(2)+k(2))*dv1(2) -(dp4(3)+k(3))*dv1(3)
      tv2=-(dp3(0)-k(0))*dv2(0) +(dp3(1)-k(1))*dv2(1)
     &    +(dp3(2)-k(2))*dv2(2) +(dp3(3)-k(3))*dv2(3)
      tv3= (dp2(0)-k(0))*dv3(0) -(dp2(1)-k(1))*dv3(1)
     &    -(dp2(2)-k(2))*dv3(2) -(dp2(3)-k(3))*dv3(3)
      tv4=-(dp1(0)+k(0))*dv4(0) +(dp1(1)+k(1))*dv4(1)
     &    +(dp1(2)+k(2))*dv4(2) +(dp1(3)+k(3))*dv4(3)
c
      j12(0)=(dp1(0)-dp2(0))*v12 +sv1*dv2(0) +sv2*dv1(0)
      j12(1)=(dp1(1)-dp2(1))*v12 +sv1*dv2(1) +sv2*dv1(1)
      j12(2)=(dp1(2)-dp2(2))*v12 +sv1*dv2(2) +sv2*dv1(2)
      j12(3)=(dp1(3)-dp2(3))*v12 +sv1*dv2(3) +sv2*dv1(3)
      j34(0)=(dp3(0)-dp4(0))*v34 +sv3*dv4(0) +sv4*dv3(0)
      j34(1)=(dp3(1)-dp4(1))*v34 +sv3*dv4(1) +sv4*dv3(1)
      j34(2)=(dp3(2)-dp4(2))*v34 +sv3*dv4(2) +sv4*dv3(2)
      j34(3)=(dp3(3)-dp4(3))*v34 +sv3*dv4(3) +sv4*dv3(3)
c
      j14(0)=(dp1(0)-dp4(0))*v14 +tv1*dv4(0) +tv4*dv1(0)
      j14(1)=(dp1(1)-dp4(1))*v14 +tv1*dv4(1) +tv4*dv1(1)
      j14(2)=(dp1(2)-dp4(2))*v14 +tv1*dv4(2) +tv4*dv1(2)
      j14(3)=(dp1(3)-dp4(3))*v14 +tv1*dv4(3) +tv4*dv1(3)
      j32(0)=(dp3(0)-dp2(0))*v23 +tv3*dv2(0) +tv2*dv3(0)
      j32(1)=(dp3(1)-dp2(1))*v23 +tv3*dv2(1) +tv2*dv3(1)
      j32(2)=(dp3(2)-dp2(2))*v23 +tv3*dv2(2) +tv2*dv3(2)
      j32(3)=(dp3(3)-dp2(3))*v23 +tv3*dv2(3) +tv2*dv3(3)
c
      js12=q(0)*j12(0)-q(1)*j12(1)-q(2)*j12(2)-q(3)*j12(3)
      js34=q(0)*j34(0)-q(1)*j34(1)-q(2)*j34(2)-q(3)*j34(3)
      js14=k(0)*j14(0)-k(1)*j14(1)-k(2)*j14(2)-k(3)*j14(3)
      js32=k(0)*j32(0)-k(1)*j32(1)-k(2)*j32(2)-k(3)*j32(3)
c
      js=j12(0)*j34(0)-j12(1)*j34(1)-j12(2)*j34(2)-j12(3)*j34(3)
      jt=j14(0)*j32(0)-j14(1)*j32(1)-j14(2)*j32(2)-j14(3)*j32(3)
c
      dvertx = (v12*v34 +v14*v23 -r_two*v13*v24)*dgw2

c     &        +(dzs*dgwwz2+das*dgwwa2)*js -dzs*dgwwz2*js12*js34/dmz**2
c     &        +(dzt*dgwwz2+dat*dgwwa2)*jt -dzt*dgwwz2*js14*js32/dmz**2
c
      vertex = -dcmplx( dvertx )
c
      return
      end
