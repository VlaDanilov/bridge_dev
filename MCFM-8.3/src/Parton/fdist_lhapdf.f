*****************
* LHAPDF version*
*****************
      subroutine fdist(ih,x,xmu,fx)
      implicit none
      double precision fx(-5:5),x,xmu,fPDF(-6:6),fphoton
      integer Iprtn,ih
      logical has_photon
c---  ih1=+1 proton 
c---  ih1=-1 pbar

C---set to zero if x out of range
      if (x .ge. 1d0) then
          do Iprtn=-5,5
             fx(Iprtn)=0d0
          enddo
          return
      endif
!$omp critical(LHApdf) 
      if (has_photon()) then
        call evolvePDFphoton(x,xmu,fPDF,fphoton)
      else
        call evolvePDF(x,xmu,fPDF)
      endif
!$omp end critical(LHApdf) 
      if (ih.eq.1) then
        do Iprtn=-5,5
          fx(+Iprtn)=fPDF(+Iprtn)/x
        enddo
      elseif(ih.eq.-1) then
        do Iprtn=-5,5
          fx(+Iprtn)=fPDF(-Iprtn)/x
        enddo
      endif
                     
      return
      end

  

