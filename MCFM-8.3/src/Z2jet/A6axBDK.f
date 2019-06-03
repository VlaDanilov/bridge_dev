      function A6axBDKpp(p1,p2,p3,p4,p5,p6,za,zb,xInt)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'qqbggintnames.f'
      integer p1,p2,p3,p4,p5,p6,bub123
      complex(dp)::A6axBDKpp,bitexch34,xInt(15),Ltm1,Lt0,Lt1,zab2
      real(dp)::s3,s12,s34,s56,s14,s123,s124
C--- begin statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      
      bitexch34(bub123,s123,p3,p4)=
     & -za(p2,p5)**2/(za(p1,p2)*za(p5,p6)*za(p3,p4)**2)*Ltm1(bub123,b56,s123,s56,xInt)
     & +za(p2,p4)*zb(p4,p6)*za(p2,p5)/(za(p1,p2)*za(p3,p4)**2*s56) 
     & *(s34/s56*Lt1(bub123,b56,s123,s56,xInt)+Lt0(bub123,b56,s123,s56,xInt))
     & +za(p5,p3)*zb(p3,p1)*za(p2,p5)/(za(p5,p6)*za(p3,p4)**2)*Lt0(bub123,b12,s123,s12,xInt)/s12
C--- end statement functions
      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s14=s(p1,p4)
      s123=s3(p1,p2,p3)
      s124=s3(p1,p2,p4)

      A6axBDKpp=
     & bitexch34(b123,s123,p3,p4)-bitexch34(b124,s124,p4,p3)
     & -(s14+s34)*za(p2,p5)*zb(p4,p6)/(za(p1,p3)*za(p3,p4))/(s56**2)*Lt1(b123,b56,s123,s56,xInt)
     & -za(p2,p3)*zb(p3,p1)*za(p2,p5)*zb(p3,p6)/(za(p2,p4)*za(p3,p4))/(s56**2)*Lt1(b124,b56,s124,s56,xInt)

      return
      end
      

      function A6axBDKpm(p1,p2,p3,p4,p5,p6,za,zb,xInt)
c--- the full pm amplitude is given by  A6axBDKpm + flip2[A6axBDKpm]
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'qqbggintnames.f'
      integer p1,p2,p3,p4,p5,p6
      complex(dp)::A6axBDKpm,Cax,xInt(15),Lt0,Lt1,zab2
      real(dp)::s3,s12,s34,s56,s123,s124
C--- begin statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C--- end statement functions
      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s123=s3(p1,p2,p3)
      s124=s3(p1,p2,p4)

      A6axBDKpm= -Cax(p1,p2,p3,p4,p5,p6,za,zb,xInt)
     &  +(za(p2,p4)*za(p1,p4)*zb(p4,p6)*zab2(p2,p1,p3,p6))/(za(p1,p2)*za(p1,p3)*zb(p5,p6)*zab2(p3,p1,p2,p4)) 
     & *Lt1(b56,b123,s56,s123,xInt)/(s123) 
     &   +(zab2(p2,p1,p3,p6)*zab2(p3,p1,p2,p6)*zb(p1,p3))/(zb(p5,p6)*zab2(p3,p1,p2,p4)**2) 
     & *Lt0(b123,b12,s123,s12,xInt)/(s12)  
     &   +(za(p2,p4)*zab2(p1,p2,p3,p4)*zab2(p2,p1,p3,p6)*zab2(p3,p1,p2,p6))
     &   /(za(p1,p2)*za(p1,p3)*zb(p5,p6)*zab2(p3,p1,p2,p4)**2) 
     & *Lt0(b123,b56,s123,s56,xInt)/(s56)  
     &   -(za(p2,p4)*za(p3,p5)*zab2(p4,p1,p3,p6))/(za(p1,p3)*za(p3,p4)*s56*zab2(p3,p1,p2,p4)) ! +(\rm flip)_2

      return
      end
      

      function A6axBDKmp(p1,p2,p3,p4,p5,p6,za,zb,xInt)
c--- the full mp amplitude is given by  A6axBDKmp + flip2[A6axBDKmp]
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'qqbggintnames.f'
      integer p1,p2,p3,p4,p5,p6
      complex(dp)::A6axBDKmp,Cax,xInt(15),xInt34(15),Lt0,Lt1,zab2
      real(dp)::s3,s12,s34,s56,s123,s124
C--- begin statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C--- end statement functions
      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s123=s3(p1,p2,p3)
      s124=s3(p1,p2,p4)

c--- fill bubble integral arrays with 3 and 4 exchanged
      xInt34(:)=xInt(:)
      xInt34(b123)=xInt(b124)
      xInt34(b124)=xInt(b123)

      A6axBDKmp= Cax(p1,p2,p4,p3,p5,p6,za,zb,xInt34)
     &  -(zb(p1,p4)**2*za(p4,p5)*zab2(p5,p2,p3,p1))/(zb(p1,p2)*zb(p1,p3)*za(p5,p6)*zab2(p4,p1,p2,p3))
     &  *Lt1(b56,b123,s56,s123,xInt)/(s123)   
     &   +(zab2(p5,p2,p3,p1)*zab2(p5,p1,p2,p3)*za(p2,p3))/(za(p5,p6)*zab2(p4,p1,p2,p3)**2)
     &  *Lt0(b123,b12,s123,s12,xInt)/(s12)   
     &   - (zb(p1,p4)*zab2(p4,p2,p3,p1)*zab2(p5,p2,p3,p1)*zab2(p5,p1,p2,p3))
     &    /(zb(p1,p2)*zb(p1,p3)*za(p5,p6)*zab2(p4,p1,p2,p3)**2)
     & *Lt0(b123,b56,s123,s56,xInt)/(s56) 
     &   -zb(p1,p4)**2*za(p2,p5)*zb(p3,p6)/(zb(p1,p3)*zb(p3,p4)*s56*zab2(p4,p1,p2,p3)) ! +(\rm flip)_2

      return
      end
      


      function Cax(p1,p2,p3,p4,p5,p6,za,zb,xInt)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'qqbggintnames.f'
      integer p1,p2,p3,p4,p5,p6
      complex(dp)::Cax,C1ax,xInt(15),xInt1625(15),Ltm1,zab2
      real(dp)::s3,s12,s34,s56,s123,s124,delta12,delta34,delta56,Delta3
C--- begin statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C--- end statement functions

      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      delta12=s12-s34-s56
      delta34=s34-s12-s56
      delta56=s56-s12-s34
      Delta3=s12**2+s34**2+s56**2-2._dp*(s12*s34+s34*s56+s56*s12)
      s123=s3(p1,p2,p3)
      s124=s3(p1,p2,p4)

      Cax= - (-3._dp/2._dp*(
     &          za(p5,p2)*zb(p2,p1)*za(p2,p1)*zb(p1,p6)
     &         +za(p5,p6)*zb(p6,p1)*za(p2,p5)*zb(p5,p6)
     &         -za(p5,p3)*zb(p3,p1)*za(p2,p4)*zb(p4,p6)
     &         -za(p5,p4)*zb(p4,p1)*za(p2,p3)*zb(p3,p6))
     &        * (zab2(p4,p1,p2,p3))/(zab2(p3,p1,p2,p4)*Delta3)

     &  -3*(delta34*(za(p5,p2)*zb(p2,p1)*delta12-za(p5,p6)*zb(p6,p1)*delta56)*zab2(p4,p1,p2,p3)*zab2(p2,p1,p3,p6))
     &      /(zab2(p3,p1,p2,p4)*Delta3**2)
     &      -zb(p1,p3)*za(p4,p5)*za(p2,p4)*zb(p3,p6)/(Delta3)

     &   + (zb(p1,p4)*za(p3,p5)*(s123-s124)*zab2(p4,p1,p2,p3)*zab2(p2,p1,p3,p6))/(zab2(p3,p1,p2,p4)**2*Delta3) 
     &   -1._dp/2._dp*(zb(p1,p3)*za(p4,p5)*zab2(p2,p1,p3,p6))/(s123*zab2(p3,p1,p2,p4))
     &   -1._dp/2._dp*(zab2(p2,p1,p3,p4)**2*zab2(p3,p1,p2,p6)**2-za(p2,p3)**2*zb(p4,p6)**2*s123**2)
     &     /(za(p1,p2)*zb(p5,p6)*zab2(p3,p1,p2,p4)**4)
     &  * (delta34/2._dp + s12*s56/(s123) )  )
     &    *xInt(c12_34)

     & + (zab2(p2,p1,p3,p6)**2)/(za(p1,p2)*zb(p5,p6)*zab2(p3,p1,p2,p4)**2)*Ltm1(b56,b34,s56,s34,xInt)  

     & + za(p2,p4)*zb(p3,p6)/(zab2(p3,p1,p2,p4)) 
     &  * (za(p2,p4)*zb(p4,p6)*delta34/(za(p1,p2)*zb(p5,p6)*Delta3)
     &          -za(p2,p4)*za(p3,p5)*delta56/(za(p1,p2)*za(p3,p4)*Delta3)
     & -zb(p1,p3)*zb(p4,p6)*delta12/(zb(p3,p4)*zb(p5,p6)*Delta3)
     & -2*za(p5,p3)*zb(p3,p1)/(Delta3)  
     &  +za(p2,p4)*za(p3,p5)/(za(p1,p2)*za(p3,p4)*s56))

c--- fill bubble integral array element with 1<->6 and 2<->5 exchange (only b12 and b34 used)
      xInt1625(:)=xInt(:)
      xInt1625(b12)=xInt(b56)
      Cax=Cax+C1ax(p1,p2,p3,p4,p5,p6,za,zb,xInt)+C1ax(p6,p5,p3,p4,p2,p1,za,zb,xInt1625)

c     & +C**(\ax)1+C**(\ax)1(1rightarrow 6,2rightarrow 5)  

      return
      end
      


      function C1ax(p1,p2,p3,p4,p5,p6,za,zb,xInt)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'qqbggintnames.f'
      integer p1,p2,p3,p4,p5,p6
      complex(dp)::C1ax,xInt(15),Ltm1,zab2
      real(dp)::s3,s12,s34,s56,s123,s124,delta34,Delta3
C--- begin statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C--- end statement functions

      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      delta34=s34-s12-s56
      Delta3=s12**2+s34**2+s56**2-2._dp*(s12*s34+s34*s56+s56*s12)
      s123=s3(p1,p2,p3)
      s124=s3(p1,p2,p4)

      C1ax= (
     &  -6*(zb(p1,p2)*zab2(p2,p1,p3,p6)*(za(p2,p5)*delta34-2*za(p2,p1)*zb(p1,p6)*za(p6,p5))*zab2(p4,p1,p2,p3))
     &   /(zab2(p3,p1,p2,p4)*Delta3**2)
     &  - (zb(p1,p3)*zb(p4,p6)*zab2(p2,p1,p3,p6))/(zb(p3,p4)*zb(p5,p6)*zab2(p3,p1,p2,p4)**2)
     &   +zb(p1,p4) 
     & *(zab2(p2,p1,p3,p6)*(3*zab2(p3,p1,p2,p4)*zb(p3,p6) -zb(p4,p6)*(s123-s124))*zab2(p4,p1,p2,p3))
     &   /(zb(p3,p4)*zb(p5,p6)*zab2(p3,p1,p2,p4)**2*Delta3)
     &   -zb(p1,p3)*za(p2,p4)*zb(p3,p6)**2/(zb(p3,p4)*zb(p5,p6)*Delta3))*Ltm1(b12,b34,s12,s34,xInt)

      return
      end
      


 
