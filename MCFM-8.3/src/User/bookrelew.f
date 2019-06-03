      subroutine bookrelew(n,tag,title,var,wt,wtnoew,xmin,xmax,xbin)
c--- routine to set up histograms that contain the relative effect
c--- of electroweak corrections; should be called from nplotter
      implicit none
      include 'types.f'
      integer n
      integer tag
      character*(*) title
      real(dp):: var,wt,wtnoew,xmin,xmax,xbin,wtwithew
      
      wtwithew=wtnoew+wt
      
      call bookplot(n,tag,trim(title)//' - no EW',var,wtnoew,wtnoew**2,
     &               xmin,xmax,xbin,'lin')
      n=n+1
      call bookplot(n,tag,trim(title)//' - with EW',var,wtwithew,wtwithew**2,
     &               xmin,xmax,xbin,'lin')
      n=n+1
      call bookplot(n,tag,trim(title)//' - +RELEW+',var,wtwithew,wtwithew**2,
     &               xmin,xmax,xbin,'lin')
      n=n+1
      
      return
      end
      
