c
c ----------------------------------------------------------------------
c
      subroutine jvvxxx(v1,v2,g,vmass,vwidth , jvv)
c
c this subroutine computes an off-shell vector current from the three-  
c point gauge boson coupling.  the vector propagator is given in feynman
c gauge for a massless vector and in unitary gauge for a massive vector.
c                                                                       
c input:                                                                
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       real    g              : coupling constant (see the table below)
c       real    vmass          : mass  of output vector v               
c       real    vwidth         : width of output vector v               
c                                                                       
c the possible sets of the inputs are as follows:                       
c    ------------------------------------------------------------------ 
c    |   v1   |   v2   |  jvv   |      g       |   vmass  |  vwidth   | 
c    ------------------------------------------------------------------ 
c    |   w-   |   w+   |  a/z   |  gwwa/gwwz   | 0./zmass | 0./zwidth | 
c    | w3/a/z |   w-   |  w+    | gw/gwwa/gwwz |   wmass  |  wwidth   | 
c    |   w+   | w3/a/z |  w-    | gw/gwwa/gwwz |   wmass  |  wwidth   | 
c    ------------------------------------------------------------------ 
c where all the bosons are defined by the flowing-out quantum number.   
c                                                                       
c output:                                                               
c       complex jvv(6)         : vector current            j^mu(v:v1,v2)
c
      complex*16 v1(6),v2(6),jvv(6),j12(0:3),js,dg,
     &        sv1,sv2,s11,s12,s21,s22,v12
      real*8    p1(0:3),p2(0:3),q(0:3),g,vmass,vwidth,gs,s,vm2,m1,m2
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
      if ( vmass .ne. r_zero ) then
         vm2=vmass**2
         m1=p1(0)**2-(p1(1)**2+p1(2)**2+p1(3)**2)
         m2=p2(0)**2-(p2(1)**2+p2(2)**2+p2(3)**2)
         s11=p1(0)*v1(1)-p1(1)*v1(2)-p1(2)*v1(3)-p1(3)*v1(4)
         s12=p1(0)*v2(1)-p1(1)*v2(2)-p1(2)*v2(3)-p1(3)*v2(4)
         s21=p2(0)*v1(1)-p2(1)*v1(2)-p2(2)*v1(3)-p2(3)*v1(4)
         s22=p2(0)*v2(1)-p2(1)*v2(2)-p2(2)*v2(3)-p2(3)*v2(4)
         js=(v12*(-m1+m2) +s11*s12 -s21*s22)/vm2
c
         dg=-g/dcmplx( s-vm2 , max(sign( vmass*vwidth ,s),r_zero) )
c
c  for the running width, use below instead of the above dg.
c         dg=-g/dcmplx( s-vm2 , max( vwidth*s/vmass ,r_zero) )
c
         jvv(1) = dg*(j12(0)-q(0)*js)
         jvv(2) = dg*(j12(1)-q(1)*js)
         jvv(3) = dg*(j12(2)-q(2)*js)
         jvv(4) = dg*(j12(3)-q(3)*js)
c
      else
         gs=-g/s
c
         jvv(1) = gs*j12(0)
         jvv(2) = gs*j12(1)
         jvv(3) = gs*j12(2)
         jvv(4) = gs*j12(3)
      end if
c
      return
      end
