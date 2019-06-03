c
c ----------------------------------------------------------------------
c
      subroutine rotxxx(p,q , prot)
c
c this subroutine performs the spacial rotation of a four-momentum.     
c the momentum p is assumed to be given in the frame where the spacial  
c component of q points the positive z-axis.  prot is the momentum p    
c rotated to the frame where q is given.                                
c                                                                       
c input:                                                                
c       real    p(0:3)         : four-momentum p in q(1)=q(2)=0 frame   
c       real    q(0:3)         : four-momentum q in the rotated frame   
c                                                                       
c output:                                                               
c       real    prot(0:3)      : four-momentum p in the rotated frame   
c
      real*8    p(0:3),q(0:3),prot(0:3),qt2,qt,psgn,qq,p1
c
      real*8 r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )
c
      prot(0) = p(0)
c
      qt2=q(1)**2+q(2)**2
c
      if ( qt2 .eq. r_zero ) then
          if ( q(3) .eq. r_zero ) then
             prot(1) = p(1)
             prot(2) = p(2)
             prot(3) = p(3)
          else
             psgn=dsign(r_one,q(3))
             prot(1) = p(1)*psgn
             prot(2) = p(2)*psgn
             prot(3) = p(3)*psgn
          endif
      else
          qq=sqrt(qt2+q(3)**2)
          qt=sqrt(qt2)
          p1=p(1)
          prot(1) = q(1)*q(3)/qq/qt*p1 -q(2)/qt*p(2) +q(1)/qq*p(3)
          prot(2) = q(2)*q(3)/qq/qt*p1 +q(1)/qt*p(2) +q(2)/qq*p(3)
          prot(3) =          -qt/qq*p1               +q(3)/qq*p(3)
      endif
c
      return
      end
