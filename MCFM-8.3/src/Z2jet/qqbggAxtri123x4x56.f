      subroutine qqbggAxtri123x4x56(p1,p2,p3,p4,p5,p6,za,zb,
     & coeff0,coeff2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6
      real(dp):: s12,s123,s56
      complex(dp):: coeff0(2,2),coeff2(2,2),zab2

C--- begin statement function
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C--- end statement functions
      s12=s(p1,p2)
      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s56=s(p5,p6)

c--- Note: coeff0 will be filled later using infrared relations
      coeff0(:,:) = czip

      include 'tri123x4x56coeffs.f'

      return
      end

