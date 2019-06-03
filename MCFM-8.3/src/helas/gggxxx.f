c
c ======================================================================
c
      subroutine gggxxx(wm,wp,w3,g , vertex)
c
c this subroutine computes an amplitude of the three-point coupling of  
c the gauge bosons.                                                     
c                                                                       
c input:                                                                
c       complex wm(6)          : vector               flow-out w-       
c       complex wp(6)          : vector               flow-out w+       
c       complex w3(6)          : vector               j3 or a    or z   
c       real    g              : coupling constant    gw or gwwa or gwwz
c                                                                       
c output:                                                               
c       complex vertex         : amplitude               gamma(wm,wp,w3)
c
      complex*16 wm(6),wp(6),w3(6),vertex, 
     &        xv1,xv2,xv3,v12,v23,v31,p12,p13,p21,p23,p31,p32
      real*8    pwm(0:3),pwp(0:3),pw3(0:3),g
c
      real*8 r_zero, r_tenth
      parameter( r_zero=0.0d0, r_tenth=0.1d0 )
c
      pwm(0)=dble( wm(5))
      pwm(1)=dble( wm(6))
      pwm(2)=dimag(wm(6))
      pwm(3)=dimag(wm(5))
      pwp(0)=dble( wp(5))
      pwp(1)=dble( wp(6))
      pwp(2)=dimag(wp(6))
      pwp(3)=dimag(wp(5))
      pw3(0)=dble( w3(5))
      pw3(1)=dble( w3(6))
      pw3(2)=dimag(w3(6))
      pw3(3)=dimag(w3(5))
c
      v12=wm(1)*wp(1)-wm(2)*wp(2)-wm(3)*wp(3)-wm(4)*wp(4)
      v23=wp(1)*w3(1)-wp(2)*w3(2)-wp(3)*w3(3)-wp(4)*w3(4)
      v31=w3(1)*wm(1)-w3(2)*wm(2)-w3(3)*wm(3)-w3(4)*wm(4)
      xv1=r_zero
      xv2=r_zero
      xv3=r_zero
      if ( abs(wm(1)) .ne. r_zero ) then
         if (abs(wm(1)).ge.max(abs(wm(2)),abs(wm(3)),abs(wm(4)))
     $        *r_tenth)
     &      xv1=pwm(0)/wm(1)
      endif
      if ( abs(wp(1)) .ne. r_zero) then
         if (abs(wp(1)).ge.max(abs(wp(2)),abs(wp(3)),abs(wp(4)))
     $        *r_tenth)
     &      xv2=pwp(0)/wp(1)
      endif
      if ( abs(w3(1)) .ne. r_zero) then
         if ( abs(w3(1)).ge.max(abs(w3(2)),abs(w3(3)),abs(w3(4)))
     $        *r_tenth)
     &      xv3=pw3(0)/w3(1)
      endif
      p12= (pwm(0)-xv1*wm(1))*wp(1)-(pwm(1)-xv1*wm(2))*wp(2)
     &    -(pwm(2)-xv1*wm(3))*wp(3)-(pwm(3)-xv1*wm(4))*wp(4)
      p13= (pwm(0)-xv1*wm(1))*w3(1)-(pwm(1)-xv1*wm(2))*w3(2)
     &    -(pwm(2)-xv1*wm(3))*w3(3)-(pwm(3)-xv1*wm(4))*w3(4)
      p21= (pwp(0)-xv2*wp(1))*wm(1)-(pwp(1)-xv2*wp(2))*wm(2)
     &    -(pwp(2)-xv2*wp(3))*wm(3)-(pwp(3)-xv2*wp(4))*wm(4)
      p23= (pwp(0)-xv2*wp(1))*w3(1)-(pwp(1)-xv2*wp(2))*w3(2)
     &    -(pwp(2)-xv2*wp(3))*w3(3)-(pwp(3)-xv2*wp(4))*w3(4)
      p31= (pw3(0)-xv3*w3(1))*wm(1)-(pw3(1)-xv3*w3(2))*wm(2)
     &    -(pw3(2)-xv3*w3(3))*wm(3)-(pw3(3)-xv3*w3(4))*wm(4)
      p32= (pw3(0)-xv3*w3(1))*wp(1)-(pw3(1)-xv3*w3(2))*wp(2)
     &    -(pw3(2)-xv3*w3(3))*wp(3)-(pw3(3)-xv3*w3(4))*wp(4)
c
      vertex = -(v12*(p13-p23)+v23*(p21-p31)+v31*(p32-p12))*g
c
      return
      end
