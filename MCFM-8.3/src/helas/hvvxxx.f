c
c ----------------------------------------------------------------------
c
      subroutine hvvxxx(v1,v2,g,smass,swidth , hvv)
c
c this subroutine computes an off-shell scalar current from the vector- 
c vector-scalar coupling.                                               
c                                                                       
c input:                                                                
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       real    g              : coupling constant                  gvvh
c       real    smass          : mass  of output scalar s               
c       real    swidth         : width of output scalar s               
c                                                                       
c output:                                                               
c       complex hvv(3)         : off-shell scalar current     j(s:v1,v2)
c
      complex*16 v1(6),v2(6),hvv(3),dg
      real*8    q(0:3),g,smass,swidth,q2
c
      real*8 r_zero
      parameter( r_zero=0.0d0 )
c
      hvv(2) = v1(5)+v2(5)
      hvv(3) = v1(6)+v2(6)
c
      q(0)=dble( hvv(2))
      q(1)=dble( hvv(3))
      q(2)=dimag(hvv(3))
      q(3)=dimag(hvv(2))
      q2=q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
c
      dg=-g/dcmplx( q2-smass**2 , max(sign( smass*swidth ,q2),r_zero) )
c
      hvv(1) = dg*(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
c
      return
      end
