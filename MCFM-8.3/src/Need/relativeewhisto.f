      subroutine relativeewhisto(N)
c--- Compute relative EW corrections in histogram N,
c--- assuming (NOEW+EW) is currently in N and (NOEW) is in N-2
      implicit none
      include 'types.f'
      include 'histo.f'
      integer N,j
      double precision xrel,xrelsig
      
      do j=1,NBIN(N)
      xrel=HIST(N,j)/HIST(N-2,j)-1d0
c--- at this point, histogram errors are stored in 2*maxhisto+N
      xrelsig=
     & xrel*sqrt((HIST(2*maxhisto+N,j)/HIST(N,j))**2
     &          +(HIST(2*maxhisto+N-2,j)/HIST(N-2,j))**2)
      HIST(N,j)=xrel
      HIST(2*maxhisto+N,j)=abs(xrelsig)
      enddo
      
      j=index(title(N),'+RELEW+')
      title(N)(j:j+6)='rel EW '
      
      return
      end
      
      
