      subroutine qqbggAxbox3x12x4(p1,p2,p3,p4,p5,p6,za,zb,coeff0,coeff2)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6
      complex(dp)::coeff0(2,2),coeff2(2,2),bit0(2,2),bit2(2,2)
      
      call qqbggAxbox3x12x4sub(p1,p2,p3,p4,p5,p6,za,zb,coeff0,coeff2)
c---- flip5 operation is this exchange and a minus sign
      call qqbggAxbox3x12x4sub(p2,p1,p3,p4,p6,p5,zb,za,bit0,bit2)
      coeff0(1,1)=-bit0(2,2)
      coeff0(2,1)=-bit0(1,2)
      coeff2(1,1)=-bit2(2,2)
      coeff2(2,1)=-bit2(1,2)
      
      return
      end
      
      subroutine qqbggAxbox3x12x4sub(p1,p2,p3,p4,p5,p6,za,zb,coeff0,coeff2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer::p1,p2,p3,p4,p5,p6
      complex(dp)::coeff0(2,2),coeff2(2,2),bitpp2,zab2
      real(dp)::s3,s12,s34,s56,s123,s124

C--- begin statement function
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      bitpp2(p3,p4)=
     & 1._dp/(2*s12*s56*za(p3,p4)**3)*( za(p2,p3)**2*za(p4,p5)*(
     &  zab2(p4,p1,p2,p4)*zb(p1,p6)*zb(p2,p3)
     & -zab2(p4,p2,p3,p6)*zb(p1,p2)*zb(p3,p4)  )  
     &+1._dp/(2)*za(p2,p3)*zab2(p4,p1,p2,p3)
     &* ( za(p3,p4)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)
     &    + za(p3,p4)*za(p3,p5)*zb(p1,p4)*zb(p3,p6)  
     &     - zab2(p4,p1,p2,p4)*za(p3,p5)*zb(p1,p6)
     &    - za(p3,p4)*za(p2,p5)*zb(p1,p2)*zb(p4,p6)
     &    - 2._dp*za(p2,p3)*za(p4,p5)*zb(p1,p6)*zb(p2,p4) )  
     &+1._dp/(2)*za(p2,p3)*za(p4,p5)*zb(p1,p6)*zab2(p3,p1,p2,p4)*zab2(p4,p1,p2,p3) )
C--- end statement functions

      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s123=s3(p1,p2,p3)
      s124=s3(p1,p2,p4)

      coeff0=czip
      coeff2=czip

C     Indices are the polarizations of the gluons, 1^+  6^+ are fixed 
      coeff0(2,2) = 
     & ((s123*s124-s12*s56)*(za(p2,p3)*za(p4,p5)+za(p2,p4)*za(p3,p5))*za(p2,p5))/(4*za(p1,p2)*za(p3,p4)**3*za(p5,p6)) 
      coeff2(2,2) = bitpp2(p3,p4)-bitpp2(p4,p3)

      coeff0(1,2) =  zip
      coeff2(1,2) =  1._dp/(2*s12*s34*s56*za(p3,p4))*(  
     &        za(p2,p3)*za(p3,p4)**2*za(p3,p5)*zb(p1,p4)*zb(p3,p4)*zb(p4,p6)  
     &      + (zab2(p3,p1,p2,p4))/(2*zab2(p4,p1,p2,p3)) * (
     &         za(p2,p3)**2*za(p2,p4)*za(p4,p5)*zb(p1,p6)*zb(p2,p3)*zb(p2,p4)
     &        + za(p1,p3)*za(p2,p4)*zab2(p3,p1,p4,p6)*za(p4,p5)*zb(p1,p4)*zb(p1,p3)  
     &         - za(p2,p3)*za(p2,p4)*za(p3,p4)*zab2(p5,p2,p3,p1)*(zb(p2,p3)*zb(p4,p6)  - zb(p2,p6)*zb(p3,p4))   
     &         - za(p2,p3)*za(p2,p4)*za(p4,p5)*zab2(p3,p1,p4,p2)*zb(p1,p3)*zb(p4,p6)
     &        + za(p1,p3)*za(p2,p4)*za(p3,p5)*zab2(p4,p2,p3,p1)*zb(p1,p4)*zb(p3,p6)   
     &         + 2._dp*za(p2,p4)*za(p3,p5)*zb(p4,p6)*zab2(p3,p1,p4,p3)*zab2(p4,p2,p3,p1) 
     &        - 2._dp*za(p1,p2)*za(p3,p4)*za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p1,p6)*zb(p2,p3) )  
     &      + (zab2(p3,p1,p2,p4))/(2._dp) * (
     &         za(p2,p3)**2*za(p4,p5)*zb(p1,p4)*zb(p2,p6)
     &        - za(p2,p3)*za(p2,p4)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)
     &        - 2._dp*za(p2,p3)*za(p3,p4)*za(p3,p5)*zb(p1,p3)*zb(p4,p6)  
     &         + zab2(p2,p1,p3,p6)*za(p3,p4)*za(p3,p5)*zb(p1,p4)
     &        + za(p2,p3)*za(p4,p5)*zb(p1,p4)*zab2(p3,p1,p4,p6) ) )

      return
      end
