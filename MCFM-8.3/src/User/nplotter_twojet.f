      subroutine nplotter_twojet(p,wt,wt2,switch)
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c--- switch:  an integer equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'outputflags.f'
      include 'nproc.f'
      include 'ewcorr.f'
      integer i3,i4,i5,nu
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: etarap,pt,
     & etaj1,ptj1,etaj2,ptj2,etaj3,ptj3,mjj,delr,getet,
     & rar,deleta,delfi,pt3,pt4,pt5,tmp3(4),tmp4(4),tmp5(4),oldpt(3:5),
     & sumeta,wt_ew,wtplot_noew,wtplot_ew,chi,yboost,wtplot,wtplotcut
      integer switch,n,nplotmax
      integer tag
      logical, save:: first=.true.
      common/nplotmax/nplotmax
  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

c--- Initialize dummy values for all quantities that could be plotted
      ptj1=-1._dp
      ptj2=-1._dp
      ptj3=-1._dp
      etaj1=99._dp
      etaj2=99._dp
      etaj3=99._dp
      mjj=-1._dp
      delr=-1._dp
      deleta=-1._dp
      
      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        goto 99
      else
c--- Add event in histograms
        tag=tagplot
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************

c--- BEGIN: order jets according to pt
      pt3=getet(p(3,4),p(3,1),p(3,2),p(3,3))
      if (jets .gt. 1) pt4=getet(p(4,4),p(4,1),p(4,2),p(4,3))
      if (jets .gt. 2) pt5=getet(p(5,4),p(5,1),p(5,2),p(5,3))      
      i3=3
      i4=4
      i5=5
      oldpt(3)=pt3
      if (jets .gt. 1) oldpt(4)=pt4
      if (jets .gt. 2) oldpt(5)=pt5
c--- sort for 2 jets 
      if (jets .eq. 2) then          
        if (pt4 .gt. pt3) then
          i3=4
          i4=3
        endif
      endif
c--- sort for 3 jets 
      if (jets .eq. 3) then
        if ((pt3 .gt. pt4) .and. (pt3 .gt. pt5)) then
           i3=3
          if (pt4 .gt. pt5) then
            i4=4
            i5=5
          else
            i4=5
            i5=4
          endif
        endif
        if ((pt4 .gt. pt3) .and. (pt4 .gt. pt5)) then
           i3=4
          if (pt3 .gt. pt5) then
            i4=3
            i5=5
          else
            i4=5
            i5=3
          endif
        endif
        if ((pt5 .gt. pt3) .and. (pt5 .gt. pt4)) then
           i3=5
          if (pt3 .gt. pt4) then
            i4=3
            i5=4
          else
            i4=4
            i5=3
          endif
        endif
      endif
c--- perform exchange
      do nu=1,4
           tmp3(nu)=p(i3,nu)
           tmp4(nu)=p(i4,nu)
           tmp5(nu)=p(i5,nu)
      enddo
      do nu=1,4
           p(3,nu)=tmp3(nu)
           p(4,nu)=tmp4(nu)
           p(5,nu)=tmp5(nu)
      enddo
      pt3=oldpt(i3)
      if (jets .gt. 1) pt4=oldpt(i4)
      if (jets .gt. 2) pt5=oldpt(i5)
c--- END: ordering

c--- Calculate quantities to plot
      etaj1=etarap(3,p)
      ptj1=pt(3,p)
      if (jets .ge. 2) then
        etaj2=etarap(4,p)
        ptj2=pt(4,p)
        deleta=abs(etaj1-etaj2)
      endif
      if (jets .ge. 3) then
        etaj3=etarap(5,p)
        ptj3=pt(5,p)
      if     (abs(etaj1-etaj3) .gt. 
     &      max(abs(etaj1-etaj2),abs(etaj2-etaj3))) then
          deleta=abs(etaj1-etaj3)
        sumeta=etaj1+etaj3
        etaj3=etaj2
      elseif (abs(etaj2-etaj3) .gt. 
     &      max(abs(etaj1-etaj2),abs(etaj1-etaj3))) then
          deleta=abs(etaj2-etaj3)
        sumeta=etaj2+etaj3
        etaj3=etaj1
      else
        deleta=abs(etaj1-etaj2)      
        sumeta=etaj1+etaj2
      endif
      endif

      rar=(p(3,1)*p(4,1)+p(3,2)*p(4,2))
     & /sqrt((p(3,1)**2+p(3,2)**2)*(p(4,1)**2+p(4,2)**2))
      if (rar.lt.-1._dp) then 
         delfi=pi
      elseif (rar.gt.1_dp) then 
         delfi=0._dp
      else
         delfi=dacos(rar)
      endif
      delr=sqrt(deleta**2+delfi**2)

      if (jets .gt. 1) then
      mjj=sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &        -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)
      endif
            
************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

c--- Call histogram routines
   99 continue

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histograms
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/dfloat(itmx))  
      endif

c--- "n" will count the number of histograms
      n=nextnplot

c--- Syntax of "bookplot" routine is:
c
c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
c
c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale

      call bookplot(n,tag,'Jet 1 pt',ptj1,wt,wt2,
     &              0._dp,5000._dp,50._dp,'log')
      n=n+1

c--- plots for comparison with CMS EXO-15-009-PAS
      yboost=(etaj1+etaj2)/2._dp
      if (abs(yboost)<1.11_dp) then
        wtplot=wt
      else
        wtplot=zip
      endif
      chi=exp(abs(etaj1-etaj2))
      if ((Mjj > 1.9e3_dp) .and. (Mjj < 2.4e3_dp)) then
        wtplotcut=wtplot/0.5e3_dp
      else
        wtplotcut=zip
      endif
      call bookplot(n,tag,'1.9 < Mjj < 2.4: chi',chi,wtplotcut,wtplotcut**2,
     &              0._dp,16._dp,1._dp,'lin')
      n=n+1 

      if ((Mjj > 2.4e3_dp) .and. (Mjj < 3.0e3_dp)) then
        wtplotcut=wtplot/0.6e3_dp
      else
        wtplotcut=zip
      endif
      call bookplot(n,tag,'2.4 < Mjj < 3.0: chi',chi,wtplotcut,wtplotcut**2,
     &              0._dp,16._dp,1._dp,'lin') 
      n=n+1

      if ((Mjj > 3.0e3_dp) .and. (Mjj < 3.6e3_dp)) then
        wtplotcut=wtplot/0.6e3_dp
      else
        wtplotcut=zip
      endif
      call bookplot(n,tag,'3.0 < Mjj < 3.6: chi',chi,wtplotcut,wtplotcut**2,
     &              0._dp,16._dp,1._dp,'lin')
      n=n+1 
     
      if ((Mjj > 3.6e3_dp) .and. (Mjj < 4.2e3_dp)) then
        wtplotcut=wtplot/0.6e3_dp
      else
        wtplotcut=zip
      endif
      call bookplot(n,tag,'3.6 < Mjj < 4.2: chi',chi,wtplotcut,wtplotcut**2,
     &              0._dp,16._dp,1._dp,'lin')
      n=n+1 
     
      if ((Mjj > 4.2e3_dp) .and. (Mjj < 4.8e3_dp)) then
        wtplotcut=wtplot/0.6e3_dp
      else
        wtplotcut=zip
      endif
      call bookplot(n,tag,'4.2 < Mjj < 4.8: chi',chi,wtplotcut,wtplotcut**2,
     &              0._dp,16._dp,1._dp,'lin')
      n=n+1
      
      if ((Mjj > 4.8e3_dp) .and. (Mjj < 13.0e3_dp)) then
        wtplotcut=wtplot/8.2e3_dp
      else
        wtplotcut=zip
      endif
      call bookplot(n,tag,'4.8 < Mjj < 13.0: chi',chi,wtplotcut,wtplotcut**2,
     &              0._dp,16._dp,1._dp,'lin')
      n=n+1 
      
      call bookplot(n,tag,'M12 low',mjj,wt,wt2,
     &              0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'M12',mjj,wt,wt2,
     &              0._dp,6000._dp,100._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 1 pt',ptj1,wt,wt2,
     &              0._dp,3000._dp,50._dp,'log')
      n=n+1

      if (kewcorr /= knone) then

        call bookrelew(n,tag,'M12',mjj,wt,wt_noew,
     &                0._dp,6000._dp,100._dp)
        call bookrelew(n,tag,'Jet 1 pt',ptj1,wt,wt_noew,
     &                0._dp,3000._dp,50._dp)
        call bookrelew(n,tag,'Jet 2 pt',ptj2,wt,wt_noew,
     &                0._dp,3000._dp,50._dp)
        wtplot_ew=wt
        wtplot_noew=wt_noew
        if  (mjj < 2000._dp) then
          wtplot_ew=zip
          wtplot_noew=zip
        endif
        call bookrelew(n,tag,'|Del y|, Mjj > 2 TeV',
     &                deleta,wtplot_ew,wtplot_noew,0._dp,6._dp,0.25_dp)

c--- binned for comparison with Dittmaier et al., JHEP11(2012)095
        call bookrelew(n,tag,'M12 > 50',mjj,wt,wt_noew,
     &                50._dp,14050._dp,14000._dp)
        call bookrelew(n,tag,'M12 > 100',mjj,wt,wt_noew,
     &                100._dp,14100._dp,14000._dp)
        call bookrelew(n,tag,'M12 > 200',mjj,wt,wt_noew,
     &                200._dp,14200._dp,14000._dp)
        call bookrelew(n,tag,'M12 > 500',mjj,wt,wt_noew,
     &                500._dp,14500._dp,14000._dp)
        call bookrelew(n,tag,'M12 > 1000',mjj,wt,wt_noew,
     &                1000._dp,15000._dp,14000._dp)
        call bookrelew(n,tag,'M12 > 2000',mjj,wt,wt_noew,
     &                2000._dp,16000._dp,14000._dp)
        call bookrelew(n,tag,'M12 > 5000',mjj,wt,wt_noew,
     &                5000._dp,19000._dp,14000._dp)

c--- binned for comparison with Dittmaier et al., JHEP11(2012)095
        call bookrelew(n,tag,'kT1 > 25',ptj1,wt,wt_noew,
     &                25._dp,7025._dp,7000._dp)
        call bookrelew(n,tag,'kT1 > 50',ptj1,wt,wt_noew,
     &                50._dp,7050._dp,7000._dp)
        call bookrelew(n,tag,'kT1 > 100',ptj1,wt,wt_noew,
     &                100._dp,7100._dp,7000._dp)
        call bookrelew(n,tag,'kT1 > 200',ptj1,wt,wt_noew,
     &                200._dp,7200._dp,7000._dp)
        call bookrelew(n,tag,'kT1 > 500',ptj1,wt,wt_noew,
     &                500._dp,7500._dp,7000._dp)
        call bookrelew(n,tag,'kT1 > 1000',ptj1,wt,wt_noew,
     &                1000._dp,8000._dp,7000._dp)
        call bookrelew(n,tag,'kT1 > 2500',ptj1,wt,wt_noew,
     &                2500._dp,9500._dp,7000._dp)

      endif
      
      call bookplot(n,tag,'Jet 1 eta',etaj1,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 2 pt',ptj2,wt,wt2,
     &              0._dp,2000._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 2 eta',etaj2,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 3 pt log',ptj3,wt,wt2,
     &              0._dp,2000._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 3 eta',etaj3,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'(jet 1,jet 2) invariant mass',mjj,wt,wt2,
     &              0._dp,425._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'delr',delr,wt,wt2,
     &              0.35_dp,4.85_dp,0.15_dp,'lin')
      n=n+1
      call bookplot(n,tag,'deleta',deleta,wt,wt2,
     &              0._dp,4._dp,0.4_dp,'lin')
      n=n+1

************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

  777 continue

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1

c--- Ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
      return 
      end
      
