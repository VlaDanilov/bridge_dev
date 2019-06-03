      function Bdiff(s34,s12,msq)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      complex(dp):: Bdiff
      
      include 'scale.f'
      real(dp):: s34,s12,msq

      Bdiff=qlI2(s34,msq,msq,musq,0)-qlI2(s12,msq,msq,musq,0)
      return
      end
