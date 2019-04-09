      subroutine zerorealhistos
c--- zero out all entries in the temporary histograms used for
c--- binning the weights in the real contribution
      implicit none
      include 'nplot.f'
      integer nplotmax,j
      common/nplotmax/nplotmax
      
      do j=1,nplotmax
      call mzero(3*maxhisto+j)
      enddo
      
      return
      end
      

      subroutine addrealhistos(wgt)
c--- add temporay histograms to the cumulative ones
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'histo.f'
      double precision wgt
      integer nplotmax
      logical added
      common/nplotmax/nplotmax

c--- loop over plots      
      do I=1,nplotmax

      added=.false.

      DO L=1,NBIN(I)

c--- add weights      
      HIST(I,L)=HIST(I,L) + HIST(3*maxhisto+I,L)
c--- add errors
      HIST(maxhisto+I,L)=HIST(maxhisto+I,L)
     . + HIST(3*maxhisto+I,L)**2*HDEL(I)
c--- the extra factor of HDEL(I) is to account for the normalization
c--- by the bin width (c.f. MFILL in mbook.f)
      
c--- count entries
      if (IHIS(3*maxhisto+I,L) .GT. 0) then
        IHIS(I,L)=IHIS(I,L)+1
        IHIS(maxhisto+I,L)=IHIS(maxhisto+I,L)+1
        added=.true.
      endif

      ENDDO

      added=.true.

c--- if any bin has been filled, increment histogram counter
      if (added) then
        IENT(I)=IENT(I)+1
        IENT(maxhisto+I)=IENT(maxhisto+I)+1
c--- otherwise, increment out of bounds counters if necessary
      else
        if     (IUSCORE(3*maxhisto+I) .GT. 0) then
          IUSCORE(I)=IUSCORE(I)+1      
          IUSCORE(maxhisto+I)=IUSCORE(maxhisto+I)+1      
        elseif (IOSCORE(3*maxhisto+I) .GT. 0) then
          IOSCORE(I)=IOSCORE(I)+1      
          IOSCORE(maxhisto+I)=IOSCORE(maxhisto+I)+1
        endif
      endif
      
      
      enddo
      
      return
      end
      
