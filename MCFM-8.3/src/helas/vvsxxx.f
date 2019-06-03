c
c ======================================================================
c
      subroutine vvsxxx(v1,v2,sc,g , vertex)
c
c this subroutine computes an amplitude of the vector-vector-scalar     
c coupling.                                                             
c                                                                       
c input:                                                                
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       complex sc(3)          : input  scalar                        s 
c       real    g              : coupling constant                  gvvh
c                                                                       
c output:                                                               
c       complex vertex         : amplitude                gamma(v1,v2,s)
c
      complex*16 v1(6),v2(6),sc(3),vertex
      real*8    g
c
      vertex = g*sc(1)*(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
c
      return
      end
