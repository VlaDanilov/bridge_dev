      subroutine pvsetmudim(mu)
      implicit none
      include 'types.f'
      include 'TRscale.f'
      real(dp):: mu
      scale=mu
      musq=scale**2
      return
      end
