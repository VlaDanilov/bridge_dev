      subroutine getptildejet(nd,pjet)
      include 'constants.f'
      include 'npart.f'
      include 'ptilde.f'
      integer nd,i,j
      double precision pjet(mxpart,4)
      
      do j=1,4
        do i=1,npart+2
        pjet(i,j)=ptildejet(nd,i,j)
        enddo
      enddo
      
      return
      end
      
