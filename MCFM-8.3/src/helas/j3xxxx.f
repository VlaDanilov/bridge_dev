c
c ----------------------------------------------------------------------
c
      subroutine j3xxxx(fi,fo,gaf,gzf,zmass,zwidth , j3)
c
c this subroutine computes the sum of photon and z currents with the    
c suitable weights ( j(w3) = cos(theta_w) j(z) + sin(theta_w) j(a) ).   
c the output j3 is useful as an input of vvvxxx, jvvxxx or w3w3xx.      
c the photon propagator is given in feynman gauge, and the z propagator 
c is given in unitary gauge.                                            
c                                                                       
c input:                                                                
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       real    gaf(2)         : fi couplings with a                 gaf
c       real    gzf(2)         : fi couplings with z                 gzf
c       real    zmass          : mass  of z                             
c       real    zwidth         : width of z                             
c                                                                       
c output:                                                               
c       complex j3(6)          : w3 current             j^mu(<fo|w3|fi>)
c
      complex*16 fi(6),fo(6),j3(6),
     &        c0l,c1l,c2l,c3l,csl,c0r,c1r,c2r,c3r,csr,dz,ddif
      real*8    gaf(2),gzf(2),q(0:3),zmass,zwidth,zm2,zmw,q2,da,ww,
     &        cw,sw,gn,gz3l,ga3l
c
      real*8 r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )
      complex*16 c_imag
      parameter( c_imag=(0d0,1d0) )
c
      j3(5) = fo(5)-fi(5)
      j3(6) = fo(6)-fi(6)
c
      q(0)=-dble( j3(5))
      q(1)=-dble( j3(6))
      q(2)=-dimag(j3(6))
      q(3)=-dimag(j3(5))
      q2=q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      zm2=zmass**2
      zmw=zmass*zwidth
c
      da=r_one/q2
      ww=max(dsign( zmw ,q2),r_zero)
      dz=r_one/dcmplx( q2-zm2 , ww )
      ddif=dcmplx( -zm2 , ww )*da*dz
c
c ddif is the difference : ddif=da-dz
c  for the running width, use below instead of the above ww,dz and ddif.
c      ww=max( zwidth*q2/zmass ,r_zero)
c      dz=r_one/dcmplx( q2-zm2 , ww )
c      ddif=dcmplx( -zm2 , ww )*da*dz
c
      cw=r_one/sqrt(r_one+(gzf(2)/gaf(2))**2)
      sw=sqrt((r_one-cw)*(r_one+cw))
      gn=gaf(2)*sw
      gz3l=gzf(1)*cw
      ga3l=gaf(1)*sw
      c0l=  fo(3)*fi(1)+fo(4)*fi(2)
      c0r=  fo(1)*fi(3)+fo(2)*fi(4)
      c1l=-(fo(3)*fi(2)+fo(4)*fi(1))
      c1r=  fo(1)*fi(4)+fo(2)*fi(3)
      c2l= (fo(3)*fi(2)-fo(4)*fi(1))*c_imag
      c2r=(-fo(1)*fi(4)+fo(2)*fi(3))*c_imag
      c3l= -fo(3)*fi(1)+fo(4)*fi(2)
      c3r=  fo(1)*fi(3)-fo(2)*fi(4)
      csl=(q(0)*c0l-q(1)*c1l-q(2)*c2l-q(3)*c3l)/zm2
      csr=(q(0)*c0r-q(1)*c1r-q(2)*c2r-q(3)*c3r)/zm2
c
      j3(1) =  gz3l*dz*(c0l-csl*q(0))+ga3l*c0l*da
     &       + gn*(c0r*ddif-csr*q(0)*dz)
      j3(2) =  gz3l*dz*(c1l-csl*q(1))+ga3l*c1l*da
     &       + gn*(c1r*ddif-csr*q(1)*dz)
      j3(3) =  gz3l*dz*(c2l-csl*q(2))+ga3l*c2l*da
     &       + gn*(c2r*ddif-csr*q(2)*dz)
      j3(4) =  gz3l*dz*(c3l-csl*q(3))+ga3l*c3l*da
     &       + gn*(c3r*ddif-csr*q(3)*dz)
c
      return
      end
