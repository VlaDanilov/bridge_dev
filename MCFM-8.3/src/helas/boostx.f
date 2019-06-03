c
c ======================================================================
c
      subroutine boostx(p,q , pboost)
c
c this subroutine performs the lorentz boost of a four-momentum.  the   
c momentum p is assumed to be given in the rest frame of q.  pboost is  
c the momentum p boosted to the frame in which q is given.  q must be a 
c timelike momentum.                                                    
c                                                                       
c input:                                                                
c       real    p(0:3)         : four-momentum p in the q rest  frame   
c       real    q(0:3)         : four-momentum q in the boosted frame   
c                                                                       
c output:                                                               
c       real    pboost(0:3)    : four-momentum p in the boosted frame   
c
      real*8    p(0:3),q(0:3),pboost(0:3),pq,qq,m,lf
c
      real*8 r_zero
      parameter( r_zero=0.0d0 )
c
      qq=q(1)**2+q(2)**2+q(3)**2
c
      if ( qq .ne. r_zero ) then
         pq=p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
         m=sqrt(q(0)**2-qq)
         lf=((q(0)-m)*pq/qq+p(0))/m
         pboost(0) = (p(0)*q(0)+pq)/m
         pboost(1) =  p(1)+q(1)*lf
         pboost(2) =  p(2)+q(2)*lf
         pboost(3) =  p(3)+q(3)*lf
      else
         pboost(0)=p(0)
         pboost(1)=p(1)
         pboost(2)=p(2)
         pboost(3)=p(3)
      endif
c
      return
      end
