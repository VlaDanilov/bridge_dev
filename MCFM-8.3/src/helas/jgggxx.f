c
c ----------------------------------------------------------------------
c
      subroutine jgggxx(w1,w2,w3,g, jw3w)
c
c this subroutine computes an off-shell w+, w-, w3, z or photon current 
c from the four-point gauge boson coupling, including the contributions 
c of w exchange diagrams.  the vector propagator is given in feynman    
c gauge for a photon and in unitary gauge for w and z bosons.  if one   
c sets wmass=0.0, then the ggg-->g current is given (see sect 2.9.1 of 
c the manual).                                                          
c                                                                       
c input:                                                                
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       complex w3(6)          : third  vector                        w3
c       real    g             : first  coupling constant               
c                                                  (see the table below)
c                                                                       
c output:                                                               
c       complex jw3w(6)        : w current             j^mu(w':w1,w2,w3)
c
      implicit none
      complex*16  w1(6),w2(6),w3(6),jw3w(6)
      complex*16 dw1(0:3),dw2(0:3),dw3(0:3),
     &           jj(0:3),dv,w32,w13
      real*8     p1(0:3),p2(0:3),p3(0:3),q(0:3),g,dg2,q2
c
      real*8 r_zero
      parameter( r_zero=0.0d0 )
c
      jw3w(5) = w1(5)+w2(5)+w3(5)
      jw3w(6) = w1(6)+w2(6)+w3(6)
c
      dw1(0)=dcmplx(w1(1))
      dw1(1)=dcmplx(w1(2))
      dw1(2)=dcmplx(w1(3))
      dw1(3)=dcmplx(w1(4))
      dw2(0)=dcmplx(w2(1))
      dw2(1)=dcmplx(w2(2))
      dw2(2)=dcmplx(w2(3))
      dw2(3)=dcmplx(w2(4))
      dw3(0)=dcmplx(w3(1))
      dw3(1)=dcmplx(w3(2))
      dw3(2)=dcmplx(w3(3))
      dw3(3)=dcmplx(w3(4))
      p1(0)=dble(      w1(5))
      p1(1)=dble(      w1(6))
      p1(2)=dble(dimag(w1(6)))
      p1(3)=dble(dimag(w1(5)))
      p2(0)=dble(      w2(5))
      p2(1)=dble(      w2(6))
      p2(2)=dble(dimag(w2(6)))
      p2(3)=dble(dimag(w2(5)))
      p3(0)=dble(      w3(5))
      p3(1)=dble(      w3(6))
      p3(2)=dble(dimag(w3(6)))
      p3(3)=dble(dimag(w3(5)))
      q(0)=-(p1(0)+p2(0)+p3(0))
      q(1)=-(p1(1)+p2(1)+p3(1))
      q(2)=-(p1(2)+p2(2)+p3(2))
      q(3)=-(p1(3)+p2(3)+p3(3))

      q2 =q(0)**2 -(q(1)**2 +q(2)**2 +q(3)**2)

      dg2=dble(g)*dble(g)
c
      dv = 1.0d0/dcmplx( q2 )

c  for the running width, use below instead of the above dv.
c      dv = 1.0d0/dcmplx( q2 -mv2 , dmax1(dwv*q2/dmv,0.d0) )
c
      w32=dw3(0)*dw2(0)-dw3(1)*dw2(1)-dw3(2)*dw2(2)-dw3(3)*dw2(3)
c
c     
      w13=dw1(0)*dw3(0)-dw1(1)*dw3(1)-dw1(2)*dw3(2)-dw1(3)*dw3(3)
c     
      jj(0)=dg2*( dw1(0)*w32 - dw2(0)*w13 )
      jj(1)=dg2*( dw1(1)*w32 - dw2(1)*w13 )
      jj(2)=dg2*( dw1(2)*w32 - dw2(2)*w13 )
      jj(3)=dg2*( dw1(3)*w32 - dw2(3)*w13 )
c     
      jw3w(1) = dcmplx( jj(0)*dv )
      jw3w(2) = dcmplx( jj(1)*dv )
      jw3w(3) = dcmplx( jj(2)*dv )
      jw3w(4) = dcmplx( jj(3)*dv )
c
      return
      end
