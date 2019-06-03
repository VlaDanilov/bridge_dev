      subroutine pvextclearcache
      implicit none
      include 'TRextclear.f'
      include 'pvRespectmaxcindex.f'
      include 'TRbadpoint.f'
      integer j
      do j=1,5
      clear(j)=.true.
      enddo
      pvRespectmaxcindex=.true.
      pvbadpoint=.false.
      return
      end
