      function A6texact(s23,mt2)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp)::A6texact
      real(dp)::s23,mt2,musq

c--- musq is irrelevant, set it to some value
      musq=abs(s23)

c--- this is the form of the top loop vacuum polarization obtained
c--- after renormalization by subtraction at zero momentum transfer
      A6texact=-two/three*(
     &  (one+two*mt2/s23)*(qlI2(s23,mt2,mt2,musq,0)-qlI2(0._dp,mt2,mt2,musq,0))
     &   +one/three)

      return
      end
