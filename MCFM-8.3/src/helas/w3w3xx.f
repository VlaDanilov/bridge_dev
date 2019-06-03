c
c ----------------------------------------------------------------------
c
      subroutine w3w3xx(wm,w31,wp,w32,g31,g32,wmass,wwidth , vertex)
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
c       real    g31            : coupling of w31 with w-/w+             
c       real    g32            : coupling of w32 with w-/w+             
c                                                  (see the table below)
c       real    wmass          : mass  of w                             
c       real    wwidth         : width of w                             
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
      complex*16    wm(6),w31(6),wp(6),w32(6),vertex
      complex*16 dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),dvertx,
     &           v12,v13,v14,v23,v24,v34
      real*8     g31,g32,wmass,wwidth
c
      real*8 r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )

      dv1(0)=dcmplx(wm(1))
      dv1(1)=dcmplx(wm(2))
      dv1(2)=dcmplx(wm(3))
      dv1(3)=dcmplx(wm(4))
      dv2(0)=dcmplx(w31(1))
      dv2(1)=dcmplx(w31(2))
      dv2(2)=dcmplx(w31(3))
      dv2(3)=dcmplx(w31(4))
      dv3(0)=dcmplx(wp(1))
      dv3(1)=dcmplx(wp(2))
      dv3(2)=dcmplx(wp(3))
      dv3(3)=dcmplx(wp(4))
      dv4(0)=dcmplx(w32(1))
      dv4(1)=dcmplx(w32(2))
      dv4(2)=dcmplx(w32(3))
      dv4(3)=dcmplx(w32(4))
c
      if ( dble(wmass) .ne. r_zero ) then
c         dm2inv = r_one / dmw2
c
         v12= dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
         v13= dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
         v14= dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
         v23= dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
         v24= dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
         v34= dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)
c
         dvertx = v12*v34 +v14*v23 -2.d0*v13*v24
c
         vertex = dcmplx( dvertx ) * (g31*g32)
c
      else
         v12= dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
         v13= dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
         v14= dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
         v23= dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
         v24= dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
         v34= dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)
c

         dvertx = v14*v23 -v13*v24
c
         vertex = dcmplx( dvertx ) * (g31*g32)
      end if
c
      return
      end
