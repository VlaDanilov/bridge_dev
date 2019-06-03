c
c ----------------------------------------------------------------------
c
      subroutine ggggxx(wm,w31,wp,w32,g, vertex)
c
c this subroutine computes an amplitude of the four-point coupling of   
c the w-, w+ and two w3/z/a.  the amplitude includes the contributions  
c of w exchange diagrams.  the internal w propagator is given in unitary
c gauge.  if one sets wmass=0.0, then the gggg vertex is given (see sect
c 2.9.1 of the manual).
c                                                                       
c input:                                                                
c       complex wm(0:3)        : flow-out w-                         wm 
c       complex w31(0:3)       : first    w3/z/a                     w31
c       complex wp(0:3)        : flow-out w+                         wp 
c       complex w32(0:3)       : second   w3/z/a                     w32
c       real    g              : coupling of w31 with w-/w+             
c                                                  (see the table below)
c                                                                       
c the possible sets of the inputs are as follows:                       
c   -------------------------------------------                         
c   |  wm  |  w31 |  wp  |  w32 |  g31 |  g32 |                         
c   -------------------------------------------                         
c   |  w-  |  w3  |  w+  |  w3  |  gw  |  gw  |                         
c   |  w-  |  w3  |  w+  |  z   |  gw  | gwwz |                         
c   |  w-  |  w3  |  w+  |  a   |  gw  | gwwa |                         
c   |  w-  |  z   |  w+  |  z   | gwwz | gwwz |                         
c   |  w-  |  z   |  w+  |  a   | gwwz | gwwa |                         
c   |  w-  |  a   |  w+  |  a   | gwwa | gwwa |                         
c   -------------------------------------------                         
c where all the bosons are defined by the flowing-out quantum number.   
c                                                                       
c output:                                                               
c       complex vertex         : amplitude          gamma(wm,w31,wp,w32)
c
      implicit none
      complex*16    wm(6),w31(6),wp(6),w32(6),vertex
      complex*16 dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),
     &           dvertx,v12,v13,v14,v23,v24,v34
      real*8       pwm(0:3),pw31(0:3),pwp(0:3),pw32(0:3),g
      real*8     dp1(0:3),dp2(0:3),dp3(0:3),dp4(0:3)
c
      real*8 r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )
c
      pwm(0)=dble( wm(5))
      pwm(1)=dble( wm(6))
      pwm(2)=dimag(wm(6))
      pwm(3)=dimag(wm(5))
      pwp(0)=dble( wp(5))
      pwp(1)=dble( wp(6))
      pwp(2)=dimag(wp(6))
      pwp(3)=dimag(wp(5))
      pw31(0)=dble( w31(5))
      pw31(1)=dble( w31(6))
      pw31(2)=dimag(w31(6))
      pw31(3)=dimag(w31(5))
      pw32(0)=dble( w32(5))
      pw32(1)=dble( w32(6))
      pw32(2)=dimag(w32(6))
      pw32(3)=dimag(w32(5))
c
      dv1(0)=dcmplx(wm(1))
      dv1(1)=dcmplx(wm(2))
      dv1(2)=dcmplx(wm(3))
      dv1(3)=dcmplx(wm(4))
      dp1(0)=dble(pwm(0))
      dp1(1)=dble(pwm(1))
      dp1(2)=dble(pwm(2))
      dp1(3)=dble(pwm(3))
      dv2(0)=dcmplx(w31(1))
      dv2(1)=dcmplx(w31(2))
      dv2(2)=dcmplx(w31(3))
      dv2(3)=dcmplx(w31(4))
      dp2(0)=dble(pw31(0))
      dp2(1)=dble(pw31(1))
      dp2(2)=dble(pw31(2))
      dp2(3)=dble(pw31(3))
      dv3(0)=dcmplx(wp(1))
      dv3(1)=dcmplx(wp(2))
      dv3(2)=dcmplx(wp(3))
      dv3(3)=dcmplx(wp(4))
      dp3(0)=dble(pwp(0))
      dp3(1)=dble(pwp(1))
      dp3(2)=dble(pwp(2))
      dp3(3)=dble(pwp(3))
      dv4(0)=dcmplx(w32(1))
      dv4(1)=dcmplx(w32(2))
      dv4(2)=dcmplx(w32(3))
      dv4(3)=dcmplx(w32(4))
      dp4(0)=dble(pw32(0))
      dp4(1)=dble(pw32(1))
      dp4(2)=dble(pw32(2))
      dp4(3)=dble(pw32(3))
c
      v12= dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
      v13= dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
      v14= dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
      v23= dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
      v24= dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
      v34= dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)

         dvertx = v14*v23 -v13*v24
c
         vertex = dcmplx( dvertx ) * (g*g)
c
      return
      end
