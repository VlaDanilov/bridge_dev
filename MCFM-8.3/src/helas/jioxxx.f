c
c ----------------------------------------------------------------------
c
      subroutine jioxxx(fi,fo,g,vmass,vwidth , jio)
c
c this subroutine computes an off-shell vector current from an external 
c fermion pair.  the vector boson propagator is given in feynman gauge  
c for a massless vector and in unitary gauge for a massive vector.      
c                                                                       
c input:                                                                
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       real    g(2)           : coupling constants                  gvf
c       real    vmass          : mass  of output vector v               
c       real    vwidth         : width of output vector v               
c                                                                       
c output:                                                               
c       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
c
      complex*16 fi(6),fo(6),jio(6),c0,c1,c2,c3,cs,d
      real*8    g(2),q(0:3),vmass,vwidth,q2,vm2,dd
c
      real*8 r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )
      complex*16 c_imag
      parameter( c_imag=(0d0,1d0) )
c
      jio(5) = fo(5)-fi(5)
      jio(6) = fo(6)-fi(6)
c
      q(0)=dble( jio(5))
      q(1)=dble( jio(6))
      q(2)=dimag(jio(6))
      q(3)=dimag(jio(5))
      q2=q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2=vmass**2
c
      if (vmass.ne.r_zero) then
c
         d=r_one/dcmplx( q2-vm2 , max(sign( vmass*vwidth ,q2),r_zero) )
c  for the running width, use below instead of the above d.
c      d=r_one/dcmplx( q2-vm2 , max( vwidth*q2/vmass ,r_zero) )
c
         if (g(2).ne.r_zero) then
c
            c0=  g(1)*( fo(3)*fi(1)+fo(4)*fi(2))
     &          +g(2)*( fo(1)*fi(3)+fo(2)*fi(4))
            c1= -g(1)*( fo(3)*fi(2)+fo(4)*fi(1))
     &          +g(2)*( fo(1)*fi(4)+fo(2)*fi(3))
            c2=( g(1)*( fo(3)*fi(2)-fo(4)*fi(1)) 
     &          +g(2)*(-fo(1)*fi(4)+fo(2)*fi(3)))*c_imag
            c3=  g(1)*(-fo(3)*fi(1)+fo(4)*fi(2))
     &          +g(2)*( fo(1)*fi(3)-fo(2)*fi(4))
         else
c
            d=d*g(1)
            c0=  fo(3)*fi(1)+fo(4)*fi(2)
            c1= -fo(3)*fi(2)-fo(4)*fi(1)
            c2=( fo(3)*fi(2)-fo(4)*fi(1))*c_imag
            c3= -fo(3)*fi(1)+fo(4)*fi(2)
         end if
c
         cs=(q(0)*c0-q(1)*c1-q(2)*c2-q(3)*c3)/vm2
c
         jio(1) = (c0-cs*q(0))*d
         jio(2) = (c1-cs*q(1))*d
         jio(3) = (c2-cs*q(2))*d
         jio(4) = (c3-cs*q(3))*d
c
      else
         dd=r_one/q2
c
         if (g(2).ne.r_zero) then
            jio(1) = ( g(1)*( fo(3)*fi(1)+fo(4)*fi(2))
     &                +g(2)*( fo(1)*fi(3)+fo(2)*fi(4)) )*dd
            jio(2) = (-g(1)*( fo(3)*fi(2)+fo(4)*fi(1))
     &                +g(2)*( fo(1)*fi(4)+fo(2)*fi(3)) )*dd
            jio(3) = ( g(1)*( fo(3)*fi(2)-fo(4)*fi(1))
     &                +g(2)*(-fo(1)*fi(4)+fo(2)*fi(3)))
     $           *dcmplx(r_zero,dd)
            jio(4) = ( g(1)*(-fo(3)*fi(1)+fo(4)*fi(2))
     &                +g(2)*( fo(1)*fi(3)-fo(2)*fi(4)) )*dd
c
         else
            dd=dd*g(1)
c
            jio(1) =  ( fo(3)*fi(1)+fo(4)*fi(2))*dd
            jio(2) = -( fo(3)*fi(2)+fo(4)*fi(1))*dd
            jio(3) =  ( fo(3)*fi(2)-fo(4)*fi(1))*dcmplx(r_zero,dd)
            jio(4) =  (-fo(3)*fi(1)+fo(4)*fi(2))*dd
         end if
      end if
c
      return
      end
