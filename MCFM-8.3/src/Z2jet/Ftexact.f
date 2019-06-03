      function Ftexact(s23,mt2)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp)::Ftexact
      real(dp)::s23,mt2,musq

c--- musq is irrelevant, set it to some value
      musq=abs(s23)
      
      Ftexact=
     & -(one+six*mt2*qlI3(s23,0._dp,0._dp,mt2,mt2,mt2,musq,0)
     &   +12._dp*mt2/s23*(qlI2(s23,mt2,mt2,musq,0)-qlI2(0._dp,mt2,mt2,musq,0)))
      
      return
      end
