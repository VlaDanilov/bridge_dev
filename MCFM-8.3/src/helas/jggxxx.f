c
c ----------------------------------------------------------------------
c
      subroutine jggxxx(v1,v2,g, jvv)
c
c this subroutine computes an off-shell vector current from the three-  
c point gauge boson coupling.  the vector propagator is given in feynman
c gauge for a massless vector and in unitary gauge for a massive vector.
c                                                                       
c input:                                                                
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       real    g              : coupling constant (see the table below)
c                                                                       
c output:                                                               
c       complex jvv(6)         : vector current            j^mu(v:v1,v2)
c
      complex*16 v1(6),v2(6),jvv(6),j12(0:3),
     &        sv1,sv2,v12
      real*8    p1(0:3),p2(0:3),q(0:3),g,gs,s
c
      real*8 r_zero
      parameter( r_zero=0.0d0 )
c
      jvv(5) = v1(5)+v2(5)
      jvv(6) = v1(6)+v2(6)
c
      p1(0)=dble( v1(5))
      p1(1)=dble( v1(6))
      p1(2)=dimag(v1(6))
      p1(3)=dimag(v1(5))
      p2(0)=dble( v2(5))
      p2(1)=dble( v2(6))
      p2(2)=dimag(v2(6))
      p2(3)=dimag(v2(5))
      q(0)=-dble( jvv(5))
      q(1)=-dble( jvv(6))
      q(2)=-dimag(jvv(6))
      q(3)=-dimag(jvv(5))
      s=q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
c
      v12=v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4)
      sv1= (p2(0)-q(0))*v1(1) -(p2(1)-q(1))*v1(2)
     &    -(p2(2)-q(2))*v1(3) -(p2(3)-q(3))*v1(4)
      sv2=-(p1(0)-q(0))*v2(1) +(p1(1)-q(1))*v2(2)
     &    +(p1(2)-q(2))*v2(3) +(p1(3)-q(3))*v2(4)
      j12(0)=(p1(0)-p2(0))*v12 +sv1*v2(1) +sv2*v1(1)
      j12(1)=(p1(1)-p2(1))*v12 +sv1*v2(2) +sv2*v1(2)
      j12(2)=(p1(2)-p2(2))*v12 +sv1*v2(3) +sv2*v1(3)
      j12(3)=(p1(3)-p2(3))*v12 +sv1*v2(4) +sv2*v1(4)
c
      gs=-g/s
c
      jvv(1) = gs*j12(0)
      jvv(2) = gs*j12(1)
      jvv(3) = gs*j12(2)
      jvv(4) = gs*j12(3)
c
      return
      end
