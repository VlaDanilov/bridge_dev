c
c ======================================================================
c
      subroutine vssxxx(vc,s1,s2,g , vertex)
c
c this subroutine computes an amplitude from the vector-scalar-scalar   
c coupling.  the coupling is absent in the minimal sm in unitary gauge. 
c                                                                       
c       complex vc(6)          : input  vector                        v 
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex g              : coupling constant (s1 charge)          
c                                                                       
c examples of the coupling constant g for susy particles are as follows:
c   -----------------------------------------------------------         
c   |    s1    | (q,i3) of s1  ||   v=a   |   v=z   |   v=w   |         
c   -----------------------------------------------------------         
c   | nu~_l    | (  0  , +1/2) ||   ---   |  gzn(1) |  gwf(1) |         
c   | e~_l     | ( -1  , -1/2) ||  gal(1) |  gzl(1) |  gwf(1) |         
c   | u~_l     | (+2/3 , +1/2) ||  gau(1) |  gzu(1) |  gwf(1) |         
c   | d~_l     | (-1/3 , -1/2) ||  gad(1) |  gzd(1) |  gwf(1) |         
c   -----------------------------------------------------------         
c   | e~_r-bar | ( +1  ,  0  ) || -gal(2) | -gzl(2) | -gwf(2) |         
c   | u~_r-bar | (-2/3 ,  0  ) || -gau(2) | -gzu(2) | -gwf(2) |         
c   | d~_r-bar | (+1/3 ,  0  ) || -gad(2) | -gzd(2) | -gwf(2) |         
c   -----------------------------------------------------------         
c where the s1 charge is defined by the flowing-out quantum number.     
c                                                                       
c output:                                                               
c       complex vertex         : amplitude                gamma(v,s1,s2)
c
      complex*16 vc(6),s1(3),s2(3),vertex,g
      real*8    p(0:3)
c
      p(0)=dble( s1(2)-s2(2))
      p(1)=dble( s1(3)-s2(3))
      p(2)=dimag(s1(3)-s2(3))
      p(3)=dimag(s1(2)-s2(2))
c
      vertex = g*s1(1)*s2(1)
     &        *(vc(1)*p(0)-vc(2)*p(1)-vc(3)*p(2)-vc(4)*p(3))
c
      return
      end
