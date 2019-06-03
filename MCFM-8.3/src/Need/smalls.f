      subroutine smalls(s,npart,*)
      implicit none
      include 'types.f'
c    cut if radiated parton too close
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cutoff.f'
      include 'kprocess.f'
      integer:: npart,j,k
      real(dp):: s(mxpart,mxpart)
      
      do j=3,npart+2
      if ((-s(1,j) < cutoff_s) .or. (-s(2,j) < cutoff_s)) return 1
        do k=j+1,npart+2
        if (s(j,k) < cutoff_s) return 1
        enddo
      enddo
      
      return
      
      
      
      if (npart == 2) then
      if ( 
     &      (-s(1,4) < cutoff_s)
     & .or. (-s(2,4) < cutoff_s)
     & .or. (+s(3,4) < cutoff_s)
     & .or. (-s(1,3) < cutoff_s)
     & .or. (-s(2,3) < cutoff_s)
     & ) return 1
      
      elseif (npart == 3) then
      if ( 
     &      (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,4) < cutoff_s)
     & .or. (-s(2,4) < cutoff_s)
     & .or. (+s(4,5) < cutoff_s)
     & .or. (+s(3,4) < cutoff_s)
     & .or. (+s(3,5) < cutoff_s)
     & .or. (-s(1,3) < cutoff_s)
     & .or. (-s(2,3) < cutoff_s)
     & ) return 1

      elseif (npart == 4) then
      if ( 
     &      (-s(1,6) < cutoff_s)
     & .or. (-s(2,6) < cutoff_s)
     & .or. (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,4) < cutoff_s)
     & .or. (-s(2,4) < cutoff_s)
     & .or. (-s(1,3) < cutoff_s)
     & .or. (-s(2,3) < cutoff_s)
     & .or. (+s(3,4) < cutoff_s)
     & .or. (+s(3,5) < cutoff_s)
     & .or. (+s(3,6) < cutoff_s)
     & .or. (+s(4,5) < cutoff_s)
     & .or. (+s(4,6) < cutoff_s)
     & .or. (+s(5,6) < cutoff_s)
     & ) return 1
c      if (  (kcase==kqq_tbg) .or. (kcase==kqqtbgg)
c     & .or. (kcase==kepem3j) .or. (kcase==kW_tndk)
c     & .or. (kcase==kZ_tjet)) then
c      if ( 
c     &      (+s(3,4) < cutoff_s)
c     & .or. (+s(3,5) < cutoff_s)
c     & .or. (+s(3,6) < cutoff_s)
c     & .or. (+s(4,5) < cutoff_s)
c     & .or. (+s(4,6) < cutoff_s)
c     & ) return 1
c      endif
     
      elseif (npart == 5) then
      if ( 
     &      (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,6) < cutoff_s)
     & .or. (-s(2,6) < cutoff_s)
     & .or. (-s(1,7) < cutoff_s)
     & .or. (-s(2,7) < cutoff_s)
     & .or. (-s(1,4) < cutoff_s)
     & .or. (-s(2,4) < cutoff_s)
     & .or. (-s(1,3) < cutoff_s)
     & .or. (-s(2,3) < cutoff_s)
     & .or. (+s(3,4) < cutoff_s)
     & .or. (+s(3,5) < cutoff_s)
     & .or. (+s(3,6) < cutoff_s)
     & .or. (+s(3,7) < cutoff_s)
     & .or. (+s(4,5) < cutoff_s)
     & .or. (+s(4,6) < cutoff_s)
     & .or. (+s(4,7) < cutoff_s)
     & .or. (+s(5,6) < cutoff_s)
     & .or. (+s(5,7) < cutoff_s)
     & .or. (+s(6,7) < cutoff_s)
     & ) return 1

      elseif (npart == 6) then
      if ( 
     &      (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,6) < cutoff_s)
     & .or. (-s(2,6) < cutoff_s)
     & .or. (-s(1,7) < cutoff_s)
     & .or. (-s(2,7) < cutoff_s)
     & .or. (-s(1,8) < cutoff_s)
     & .or. (-s(2,8) < cutoff_s)
     & .or. (+s(5,6) < cutoff_s)
     & .or. (+s(5,7) < cutoff_s)
     & .or. (+s(5,8) < cutoff_s)
     & .or. (+s(6,7) < cutoff_s)
     & .or. (+s(6,8) < cutoff_s)
     & .or. (+s(7,8) < cutoff_s)
     & ) return 1
      
      elseif (npart == 7) then
      if ( 
     &      (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,6) < cutoff_s)
     & .or. (-s(2,6) < cutoff_s)
     & .or. (-s(1,7) < cutoff_s)
     & .or. (-s(2,7) < cutoff_s)
     & .or. (-s(1,8) < cutoff_s)
     & .or. (-s(2,8) < cutoff_s)
     & .or. (-s(1,9) < cutoff_s)
     & .or. (-s(2,9) < cutoff_s)
     & .or. (+s(5,6) < cutoff_s)
     & .or. (+s(5,7) < cutoff_s)
     & .or. (+s(5,8) < cutoff_s)
     & .or. (+s(5,9) < cutoff_s)
     & .or. (+s(6,7) < cutoff_s)
     & .or. (+s(6,8) < cutoff_s)
     & .or. (+s(6,9) < cutoff_s)
     & .or. (+s(7,8) < cutoff_s)
     & .or. (+s(7,9) < cutoff_s)
     & .or. (+s(8,9) < cutoff_s)
     & ) return 1
      
      endif      
      
      return
      end
