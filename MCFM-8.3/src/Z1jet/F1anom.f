      function F1anom(s12,s45,mt2,musq)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp)::F1anom
      complex(dp)::lnrat
      real(dp)::s12,s45,mt2,musq

      if (mt2 .eq. zero) then
      F1anom=one/two/(s45-s12)*(1.d0+s45/(s45-s12)*lnrat(-s12,-s45))
      return
      else
      F1anom=one/two/(s45-s12)
     & *(one+two*mt2*qlI3(s12,0.d0,s45,mt2,mt2,mt2,musq,0)
     & +(s45/(s45-s12))*(qlI2(s45,mt2,mt2,musq,0)-qlI2(s12,mt2,mt2,musq,0)))
      endif
            
      return
      end

