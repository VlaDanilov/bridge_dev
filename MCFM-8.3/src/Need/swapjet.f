      subroutine swapjet(pjet,jetindex,i,j)
        use jettagging
      implicit none
c--- swaps jets i..j in pjet
      
      include 'constants.f'
      include 'nf.f'
      include 'cplx.h'
      include 'jetlabel.f'
      integer:: i,j,k,jetindex(mxpart)
      real(dp):: pjet(mxpart,4),tmp
      integer :: itmp
      character*2 chartmp
 
c--- escape if we're trying to swap the same jets
      if (i == j) return

      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo
 
      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp

      itmp = jetcontent(i)
      jetcontent(i) = jetcontent(j)
      jetcontent(j) = itmp
      
      return
      end
