C  (C) Copr. 1986-92 Numerical Recipes Software ]2w.1,r1..

      FUNCTION pvdpythag(a,b)
      implicit none
      include 'types.f'
      real(dp):: a,b,pvdpythag
      real(dp):: absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pvdpythag=absa*sqrt(1._dp+(absb/absa)**2)
      else
        if(absb.eq.0._dp)then
          pvdpythag=0._dp
        else
          pvdpythag=absb*sqrt(1._dp+(absa/absb)**2)
        endif
      endif
      return
      END
