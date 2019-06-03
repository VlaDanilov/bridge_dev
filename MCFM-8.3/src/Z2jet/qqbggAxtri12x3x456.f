      subroutine qqbggAxtri12x3x456(p1,p2,p3,p4,p5,p6,b3,b4,za,zb,mtsq,
     & coeff)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,b3,b4
      real(dp):: mtsq,s12,s123,s34,s1234
      complex(dp):: coeff(2,2),zab2,zba2,bitpm

C--- begin statement function
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
C--- end statement functions
      s12=s(p1,p2)
      s34=s(p3,p4)
      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s1234=s(p5,p6)

      coeff(1,1)=
     & +zb(p1,p4)**2*zb(p3,p6)*za(p4,p5)
     & +zb(p1,p3)*zb(p4,p6)*zba2(p1,p2,p3,p5)
     & +zb(p1,p4)*zb(p1,p6)*zba2(p3,p2,p1,p5)
      coeff(1,1)=coeff(1,1)
     & *za(p2,p1)*za(p3,p4)/zb(p3,p4)**2*(s123-s12)/two

!      include 'bitmp.f'
!      include 'bitpm.f'
      
!      write(6,*) 'check of bitpm'
!      write(6,*) 'include file      ',coeff(2,1)
!      write(6,*) 'statement function',bitpm(p1,p2,p3,p4,p5,p6,za,zb)
!      write(6,*)
!      write(6,*) 'check of bitmp'
!      write(6,*) 'include file      ',coeff(1,2)
!      write(6,*) 'statement function',-bitpm(p2,p1,p3,p4,p6,p5,zb,za)
!      pause
      
      coeff(2,2)=
     & -za(p2,p4)**2*za(p3,p5)*zb(p4,p6)
     & -za(p2,p3)*za(p4,p5)*zab2(p2,p1,p3,p6)
     & -za(p2,p4)*za(p2,p5)*zab2(p3,p1,p2,p6)
      coeff(2,2)=coeff(2,2)
     & *zb(p1,p2)*zb(p3,p4)/za(p3,p4)**2*(s123-s12)/two

      coeff(2,1)=+bitpm(p1,p2,p3,p4,p5,p6,za,zb)
      coeff(1,2)=-bitpm(p2,p1,p3,p4,p6,p5,zb,za)

      coeff(:,:)=coeff(:,:)/(2._dp*s12*s34*s1234)
!reinstate missing propagators and factor of 1/2 from couplings
      return
      end


      function bitpm(p1,p2,p3,p4,p5,p6,za,zb)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6
      real(dp):: s34
      complex(dp):: zab2,bitpm

C--- begin statement function
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C--- end statement functions
      s34=s(p3,p4)
      
      bitpm= (
     &  zab2(p3,p1,p2,p6)*zab2(p5,p3,p4,p1) * (  - za(p1,
     &    p2)**2*zb(p1,p2)*zb(p1,p4) + 2.D0*za(p1,p2)*za(p2,p3)*zb(p1,
     &    p2)*zb(p3,p4) )
     & + zab2(p3,p1,p2,p6) * (  - za(p2,p3)**2*
     &    za(p4,p5)*zb(p1,p2)*zb(p3,p4)**2 )
     & + zab2(p5,p1,p2,p3) * (  - 2.D0*za(p1,p2)
     &    *za(p2,p3)**2*zb(p1,p2)**2*zb(p4,p6) - 2.D0*za(p2,p3)**3*zb(
     &    p1,p2)*zb(p2,p3)*zb(p4,p6) )
     & + zab2(p5,p3,p4,p1) * ( za(p1,p2)*za(p2,
     &    p3)**2*zb(p1,p2)*zb(p2,p3)*zb(p4,p6) + za(p1,p3)*za(p2,p3)**2
     &    *zb(p1,p2)*zb(p3,p4)*zb(p3,p6) )
     & + za(p1,p2)**3*za(p3,p5)*zb(p1,p2)**2*zb(
     & p1,p4)*zb(p1,p6) + za(p1,p2)**2*za(p2,p3)*za(p2,p5)*zb(p1,p2)**3
     &    *zb(p4,p6) - 3.D0*za(p1,p2)**2*za(p2,p3)*za(p3,p5)*zb(p1,p2)
     &    **2*zb(p1,p6)*zb(p3,p4) + 3.D0*za(p1,p2)*za(p2,p3)**2*za(p3,
     &    p5)*zb(p1,p2)**2*zb(p3,p4)*zb(p3,p6) - za(p1,p2)*za(p2,p3)**2
     &    *za(p4,p5)*zb(p1,p2)*zb(p1,p4)*zb(p2,p3)*zb(p4,p6) + za(p1,p3
     &    )*za(p1,p5)*za(p2,p3)**2*zb(p1,p2)*zb(p1,p3)**2*zb(p4,p6) - 
     &    za(p1,p3)*za(p2,p3)**2*za(p4,p5)*zb(p1,p2)*zb(p1,p4)*zb(p3,p4
     &    )*zb(p3,p6) - za(p2,p3)**3*za(p2,p5)*zb(p1,p2)*zb(p2,p3)**2*
     &    zb(p4,p6) + za(p2,p3)**3*za(p3,p5)*zb(p1,p2)*zb(p2,p3)*zb(p3,
     &    p4)*zb(p3,p6)
     &  ) *s34*zab2(p3,p1,p2,p3)/zab2(p3,p1,p2,p4)**3/two
     
      return
      end
      
