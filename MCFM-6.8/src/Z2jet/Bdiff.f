      double complex function Bdiff(s34,s12,msq)
      implicit none
      include 'scale.f'
      double complex qlI2
      double precision s34,s12,msq
      logical first
      data first/.true./
      save first
      if (first) then
         call qlinit
         first=.false.
      endif
      Bdiff=qlI2(s34,msq,msq,musq,0)-qlI2(s12,msq,msq,musq,0)
      return
      end
