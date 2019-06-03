c
c ----------------------------------------------------------------------
c
      subroutine jwwwxx(w1,w2,w3,gwwa,gwwz,zmass,zwidth,wmass,wwidth ,
     &                  jwww)
c
c this subroutine computes an off-shell w+/w- current from the four-    
c point gauge boson coupling, including the contributions of photon and 
c z exchanges.  the vector propagators for the output w and the internal
c z bosons are given in unitary gauge, and that of the internal photon  
c is given in feynman gauge.                                            
c                                                                       
c input:                                                                
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       complex w3(6)          : third  vector                        w3
c       real    gwwa           : coupling constant of w and a       gwwa
c       real    gwwz           : coupling constant of w and z       gwwz
c       real    zmass          : mass  of internal z                    
c       real    zwidth         : width of internal z                    
c       real    wmass          : mass  of output w                      
c       real    wwidth         : width of output w                      
c                                                                       
c the possible sets of the inputs are as follows:                       
c   ------------------------------------------------------------------- 
c   |  w1  |  w2  |  w3  |gwwa|gwwz|zmass|zwidth|wmass|wwidth || jwww | 
c   ------------------------------------------------------------------- 
c   |  w-  |  w+  |  w-  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  w+  | 
c   |  w+  |  w-  |  w+  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  w-  | 
c   ------------------------------------------------------------------- 
c where all the bosons are defined by the flowing-out quantum number.   
c                                                                       
c output:                                                               
c       complex jwww(6)        : w current             j^mu(w':w1,w2,w3)
c
      complex*16  w1(6),w2(6),w3(6),jwww(6)
      complex*16 dw1(0:3),dw2(0:3),dw3(0:3),
     &           jj(0:3),js(0:3),jt(0:3),j4(0:3),
     &           jt12(0:3),jt32(0:3),j12(0:3),j32(0:3),
     &           dzs,dzt,dw,w12,w32,w13,p1w2,p2w1,p3w2,p2w3,
     &           jk12,jk32,jsw3,jtw1,p3js,ksw3,p1jt,ktw1,jq
      real*8     gwwa,gwwz,zmass,zwidth,wmass,wwidth
      real*8     p1(0:3),p2(0:3),p3(0:3),q(0:3),ks(0:3),kt(0:3),
     &           dgwwa2,dgwwz2,dgw2,dmz,dwz,dmw,dww,mz2,mw2,q2,ks2,kt2,
     &           das,dat
c
      jwww(5) = w1(5)+w2(5)+w3(5)
      jwww(6) = w1(6)+w2(6)+w3(6)
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
      ks(0)=p1(0)+p2(0)
      ks(1)=p1(1)+p2(1)
      ks(2)=p1(2)+p2(2)
      ks(3)=p1(3)+p2(3)
      kt(0)=p2(0)+p3(0)
      kt(1)=p2(1)+p3(1)
      kt(2)=p2(2)+p3(2)
      kt(3)=p2(3)+p3(3)
      q2 =q(0)**2 -(q(1)**2 +q(2)**2 +q(3)**2)
      ks2=ks(0)**2-(ks(1)**2+ks(2)**2+ks(3)**2)
      kt2=kt(0)**2-(kt(1)**2+kt(2)**2+kt(3)**2)
      dgwwa2=dble(gwwa)**2
      dgwwz2=dble(gwwz)**2
      dgw2  =dgwwa2+dgwwz2
      dmz=dble(zmass)
      dwz=dble(zwidth)
      dmw=dble(wmass)
      dww=dble(wwidth)
      mz2=dmz**2
      mw2=dmw**2
c
      das=-dgwwa2/ks2
      dat=-dgwwa2/kt2
      dzs=-dgwwz2/dcmplx( ks2-mz2 , dmax1(dsign(dmz*dwz,ks2),0.d0) )
      dzt=-dgwwz2/dcmplx( kt2-mz2 , dmax1(dsign(dmz*dwz,kt2),0.d0) )
      dw =-1.0d0/dcmplx( q2 -mw2 , dmax1(dsign(dmw*dww,q2 ),0.d0) )
c  for the running width, use below instead of the above dw.
c      dw =-1.0d0/dcmplx( q2 -mw2 , dmax1(dww*q2/dmw,0.d0) )
c
      w12=dw1(0)*dw2(0)-dw1(1)*dw2(1)-dw1(2)*dw2(2)-dw1(3)*dw2(3)
      w32=dw3(0)*dw2(0)-dw3(1)*dw2(1)-dw3(2)*dw2(2)-dw3(3)*dw2(3)
c
      p1w2= (p1(0)+ks(0))*dw2(0)-(p1(1)+ks(1))*dw2(1)
     &     -(p1(2)+ks(2))*dw2(2)-(p1(3)+ks(3))*dw2(3)
      p2w1= (p2(0)+ks(0))*dw1(0)-(p2(1)+ks(1))*dw1(1)
     &     -(p2(2)+ks(2))*dw1(2)-(p2(3)+ks(3))*dw1(3)
      p3w2= (p3(0)+kt(0))*dw2(0)-(p3(1)+kt(1))*dw2(1)
     &     -(p3(2)+kt(2))*dw2(2)-(p3(3)+kt(3))*dw2(3)
      p2w3= (p2(0)+kt(0))*dw3(0)-(p2(1)+kt(1))*dw3(1)
     &     -(p2(2)+kt(2))*dw3(2)-(p2(3)+kt(3))*dw3(3)
c
      jt12(0)= (p1(0)-p2(0))*w12 + p2w1*dw2(0) - p1w2*dw1(0)
      jt12(1)= (p1(1)-p2(1))*w12 + p2w1*dw2(1) - p1w2*dw1(1)
      jt12(2)= (p1(2)-p2(2))*w12 + p2w1*dw2(2) - p1w2*dw1(2)
      jt12(3)= (p1(3)-p2(3))*w12 + p2w1*dw2(3) - p1w2*dw1(3)
      jt32(0)= (p3(0)-p2(0))*w32 + p2w3*dw2(0) - p3w2*dw3(0)
      jt32(1)= (p3(1)-p2(1))*w32 + p2w3*dw2(1) - p3w2*dw3(1)
      jt32(2)= (p3(2)-p2(2))*w32 + p2w3*dw2(2) - p3w2*dw3(2)
      jt32(3)= (p3(3)-p2(3))*w32 + p2w3*dw2(3) - p3w2*dw3(3)
c
      jk12=(jt12(0)*ks(0)-jt12(1)*ks(1)-jt12(2)*ks(2)-jt12(3)*ks(3))/mz2
      jk32=(jt32(0)*kt(0)-jt32(1)*kt(1)-jt32(2)*kt(2)-jt32(3)*kt(3))/mz2
c
      j12(0)=jt12(0)*(das+dzs)-ks(0)*jk12*dzs
      j12(1)=jt12(1)*(das+dzs)-ks(1)*jk12*dzs
      j12(2)=jt12(2)*(das+dzs)-ks(2)*jk12*dzs
      j12(3)=jt12(3)*(das+dzs)-ks(3)*jk12*dzs
      j32(0)=jt32(0)*(dat+dzt)-kt(0)*jk32*dzt
      j32(1)=jt32(1)*(dat+dzt)-kt(1)*jk32*dzt
      j32(2)=jt32(2)*(dat+dzt)-kt(2)*jk32*dzt
      j32(3)=jt32(3)*(dat+dzt)-kt(3)*jk32*dzt
c
      jsw3=j12(0)*dw3(0)-j12(1)*dw3(1)-j12(2)*dw3(2)-j12(3)*dw3(3)
      jtw1=j32(0)*dw1(0)-j32(1)*dw1(1)-j32(2)*dw1(2)-j32(3)*dw1(3)
c
      p3js= (p3(0)-q(0))*j12(0)-(p3(1)-q(1))*j12(1)
     &     -(p3(2)-q(2))*j12(2)-(p3(3)-q(3))*j12(3)
      ksw3= (ks(0)-q(0))*dw3(0)-(ks(1)-q(1))*dw3(1)
     &     -(ks(2)-q(2))*dw3(2)-(ks(3)-q(3))*dw3(3)
      p1jt= (p1(0)-q(0))*j32(0)-(p1(1)-q(1))*j32(1)
     &     -(p1(2)-q(2))*j32(2)-(p1(3)-q(3))*j32(3)
      ktw1= (kt(0)-q(0))*dw1(0)-(kt(1)-q(1))*dw1(1)
     &     -(kt(2)-q(2))*dw1(2)-(kt(3)-q(3))*dw1(3)
c
      js(0)= (ks(0)-p3(0))*jsw3 + p3js*dw3(0) - ksw3*j12(0)
      js(1)= (ks(1)-p3(1))*jsw3 + p3js*dw3(1) - ksw3*j12(1)
      js(2)= (ks(2)-p3(2))*jsw3 + p3js*dw3(2) - ksw3*j12(2)
      js(3)= (ks(3)-p3(3))*jsw3 + p3js*dw3(3) - ksw3*j12(3)
      jt(0)= (kt(0)-p1(0))*jtw1 + p1jt*dw1(0) - ktw1*j32(0)
      jt(1)= (kt(1)-p1(1))*jtw1 + p1jt*dw1(1) - ktw1*j32(1)
      jt(2)= (kt(2)-p1(2))*jtw1 + p1jt*dw1(2) - ktw1*j32(2)
      jt(3)= (kt(3)-p1(3))*jtw1 + p1jt*dw1(3) - ktw1*j32(3)
c
      w13=dw1(0)*dw3(0)-dw1(1)*dw3(1)-dw1(2)*dw3(2)-dw1(3)*dw3(3)
c
      j4(0)=dgw2*( dw1(0)*w32 + dw3(0)*w12 - 2.d0*dw2(0)*w13 )
      j4(1)=dgw2*( dw1(1)*w32 + dw3(1)*w12 - 2.d0*dw2(1)*w13 )
      j4(2)=dgw2*( dw1(2)*w32 + dw3(2)*w12 - 2.d0*dw2(2)*w13 )
      j4(3)=dgw2*( dw1(3)*w32 + dw3(3)*w12 - 2.d0*dw2(3)*w13 )
c
c      jj(0)=js(0)+jt(0)+j4(0)
c      jj(1)=js(1)+jt(1)+j4(1)
c      jj(2)=js(2)+jt(2)+j4(2)
c      jj(3)=js(3)+jt(3)+j4(3)

      jj(0)=j4(0)
      jj(1)=j4(1)
      jj(2)=j4(2)
      jj(3)=j4(3)
c
      jq=(jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/mw2
c

      jwww(1) = dcmplx( (jj(0)-jq*q(0))*dw )
      jwww(2) = dcmplx( (jj(1)-jq*q(1))*dw )
      jwww(3) = dcmplx( (jj(2)-jq*q(2))*dw )
      jwww(4) = dcmplx( (jj(3)-jq*q(3))*dw )
c
      return
      end

