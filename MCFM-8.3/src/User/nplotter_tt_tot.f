       subroutine nplotter_tt_tot(p,wt,wt2,switch)
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
      include 'ewcorr.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: yrap,pt,yraptwo,pttwo,r
c---  Z->e+e-(31) or b bbar(33): both measured, rapidities and momenta of 3 and 4 can
c---  be calculated, also the invariant mass m34
      real(dp):: y3,y4,y5,y34,pt3,pt4,pt5,pt34,m34,r35
      real(dp):: ylep, yjet, ptlep, ptjet,costheta,
     & p3(4),p4(4),p34(4),wt_ew,wtplot_noew,wtplot_ew,dely
      real(dp):: eta,ss,etam,binedges(30)
      integer switch,n,nplotmax
      integer tag
      logical, save::first=.true.
      common/nplotmax/nplotmax

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        y3=1d3
        y4=1d3
        y5=1d3
        y34=1d3
        pt3=0._dp
        pt4=0._dp
        pt5=1d3
        pt34=0._dp
        m34=0._dp
        eta=0._dp
        r35=1d3
        jets=1
c---  (Upper) limits for the plots
        ylep=6._dp
        yjet=3.2_dp
        ptlep=80._dp
        ptjet=100._dp
        etam=1.d3
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

      y3=yrap(3,p)
      y4=yrap(4,p)
      y34=yraptwo(3,4,p)
      pt3=pt(3,p)
      pt4=pt(4,p)
      pt34=pttwo(3,4,p)
      dely=y3-y4
      m34=sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .         -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)
      ss = (p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2
     .         -(p(1,2)+p(2,2))**2-(p(1,3)+p(2,3))**2
      eta=ss/4._dp/mt**2 - 1._dp

      if(jets .gt. 0) then
         pt5=pt(5,p)
         y5=yrap(5,p)
         r35=R(p,3,5)
      else
         pt5=-1._dp
         y5=1d3
         r35=1d3
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
           
      call bookplot(n,tag,'m34',m34,wt,wt2,500._dp,10000._dp,200._dp,'lin')
      n=n+1

      call bookplot(n,tag,'ptt CHM1',pt3,wt*half,wt2*half**2,0._dp,60._dp,60._dp,'lin')
      call bookplot(n,tag,'ptt CHM1',pt4,wt*half,wt2*half**2,0._dp,60._dp,60._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CHM2',pt3,wt*half,wt2*half**2,60._dp,100._dp,40._dp,'lin')
      call bookplot(n,tag,'ptt CHM2',pt4,wt*half,wt2*half**2,60._dp,100._dp,40._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CHM3',pt3,wt*half,wt2*half**2,100._dp,150._dp,50._dp,'lin')
      call bookplot(n,tag,'ptt CHM3',pt4,wt*half,wt2*half**2,100._dp,150._dp,50._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CHM4',pt3,wt*half,wt2*half**2,150._dp,200._dp,50._dp,'lin')
      call bookplot(n,tag,'ptt CHM4',pt4,wt*half,wt2*half**2,150._dp,200._dp,50._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CHM5',pt3,wt*half,wt2*half**2,200._dp,260._dp,60._dp,'lin')
      call bookplot(n,tag,'ptt CHM5',pt4,wt*half,wt2*half**2,200._dp,260._dp,60._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CHM6',pt3,wt*half,wt2*half**2,260._dp,320._dp,60._dp,'lin')
      call bookplot(n,tag,'ptt CHM6',pt4,wt*half,wt2*half**2,260._dp,320._dp,60._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CHM7',pt3,wt*half,wt2*half**2,320._dp,400._dp,80._dp,'lin')
      call bookplot(n,tag,'ptt CHM7',pt4,wt*half,wt2*half**2,320._dp,400._dp,80._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CHM8',pt3,wt*half,wt2*half**2,400._dp,500._dp,100._dp,'lin')
      call bookplot(n,tag,'ptt CHM8',pt4,wt*half,wt2*half**2,400._dp,500._dp,100._dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'yt CHM1',y3,wt*half,wt2*half**2,-2.5_dp,-1.6_dp,0.9_dp,'lin')
      call bookplot(n,tag,'yt CHM1',y4,wt*half,wt2*half**2,-2.5_dp,-1.6_dp,0.9_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM2',y3,wt*half,wt2*half**2,-1.6_dp,-1.2_dp,0.4_dp,'lin')
      call bookplot(n,tag,'yt CHM2',y4,wt*half,wt2*half**2,-1.6_dp,-1.2_dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM3',y3,wt*half,wt2*half**2,-1.2_dp,-0.8_dp,0.4_dp,'lin')
      call bookplot(n,tag,'yt CHM3',y4,wt*half,wt2*half**2,-1.2_dp,-0.8_dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM4',y3,wt*half,wt2*half**2,-0.8_dp,-0.4_dp,0.4_dp,'lin')
      call bookplot(n,tag,'yt CHM4',y4,wt*half,wt2*half**2,-0.8_dp,-0.4_dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM5',y3,wt*half,wt2*half**2,-0.4_dp,0._dp,0.4_dp,'lin')
      call bookplot(n,tag,'yt CHM5',y4,wt*half,wt2*half**2,-0.4_dp,0._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM6',y3,wt*half,wt2*half**2,0._dp,0.4_dp,0.4_dp,'lin')
      call bookplot(n,tag,'yt CHM6',y4,wt*half,wt2*half**2,0._dp,0.4_dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM7',y3,wt*half,wt2*half**2,0.4_dp,0.8_dp,0.4_dp,'lin')
      call bookplot(n,tag,'yt CHM7',y4,wt*half,wt2*half**2,0.4_dp,0.8_dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM8',y3,wt*half,wt2*half**2,0.8_dp,1.2_dp,0.4_dp,'lin')
      call bookplot(n,tag,'yt CHM8',y4,wt*half,wt2*half**2,0.8_dp,1.2_dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM9',y3,wt*half,wt2*half**2,1.2_dp,1.6_dp,0.4_dp,'lin')
      call bookplot(n,tag,'yt CHM9',y4,wt*half,wt2*half**2,1.2_dp,1.6_dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'yt CHM10',y3,wt*half,wt2*half**2,1.6_dp,2.5_dp,0.9_dp,'lin')
      call bookplot(n,tag,'yt CHM10',y4,wt*half,wt2*half**2,1.6_dp,2.5_dp,0.9_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'mtt CFHM1',m34,wt,wt2,240._dp,412.5_dp,172.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'mtt CFHM2',m34,wt,wt2,412.5_dp,505._dp,92.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'mtt CFHM3',m34,wt,wt2,505._dp,615.5_dp,110._dp,'lin')
      n=n+1
      call bookplot(n,tag,'mtt CFHM4',m34,wt,wt2,615._dp,750._dp,135._dp,'lin')
      n=n+1
      call bookplot(n,tag,'mtt CFHM5',m34,wt,wt2,750._dp,1200._dp,450._dp,'lin')
      n=n+1
      call bookplot(n,tag,'mtt CFHM6',m34,wt,wt2,1200._dp,2000._dp,800._dp,'lin')
      n=n+1

      call bookplot(n,tag,'ptt CFHM1',pt3,wt,wt2,0._dp,45._dp,45._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CFHM2',pt3,wt,wt2,45._dp,90._dp,45._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CFHM3',pt3,wt,wt2,90._dp,140._dp,50._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CFHM4',pt3,wt,wt2,140._dp,200._dp,60._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CFHM5',pt3,wt,wt2,200._dp,300._dp,100._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CFHM6',pt3,wt,wt2,300._dp,500._dp,200._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptt CFHM7',pt3,wt,wt2,500._dp,2000._dp,1500._dp,'lin')
      n=n+1

      call bookplot(n,tag,'ayt CFHM1',abs(y3),wt,wt2,0._dp,0.25_dp,0.25_dp,'lin')
      n=n+1
      call bookplot(n,tag,'ayt CFHM2',abs(y3),wt,wt2,0.25_dp,0.5_dp,0.25_dp,'lin')
      n=n+1
      call bookplot(n,tag,'ayt CFHM3',abs(y3),wt,wt2,0.5_dp,0.75_dp,0.25_dp,'lin')
      n=n+1
      call bookplot(n,tag,'ayt CFHM4',abs(y3),wt,wt2,0.75_dp,1._dp,0.25_dp,'lin')
      n=n+1
      call bookplot(n,tag,'ayt CFHM5',abs(y3),wt,wt2,1._dp,1.25_dp,0.25_dp,'lin')
      n=n+1
      call bookplot(n,tag,'ayt CFHM6',abs(y3),wt,wt2,1.25_dp,1.5_dp,0.25_dp,'lin')
      n=n+1
      call bookplot(n,tag,'ayt CFHM7',abs(y3),wt,wt2,1.5_dp,10._dp,8.5_dp,'lin')
      n=n+1

      call bookplot(n,tag,'Dyt CFHM1',dely,wt,wt2,-8._dp,-1.5_dp,6.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dyt CFHM2',dely,wt,wt2,-1.5_dp,-1._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dyt CFHM3',dely,wt,wt2,-1._dp,-0.5_dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dyt CFHM4',dely,wt,wt2,-0.5_dp,0._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dyt CFHM5',dely,wt,wt2,0._dp,0.5_dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dyt CFHM6',dely,wt,wt2,0.5_dp,1._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dyt CFHM7',dely,wt,wt2,1._dp,1.5_dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dyt CFHM7',dely,wt,wt2,1.5_dp,8._dp,6.5_dp,'lin')
      n=n+1

      goto 777 ! no more histograms for now

      if (kewcorr /= knone) then

c--- this is the syntax for computing the relative effect of EW
c--- corrections in distributions

        call bookrelew(n,tag,'m34',m34,wt,wt_noew,346.4_dp,646.4_dp,10._dp)
        call bookrelew(n,tag,'m34',m34,wt,wt_noew,500._dp,10000._dp,200._dp)
        call bookrelew(n,tag,'m34',m34,wt,wt_noew,0._dp,5000._dp,50._dp)
        call bookrelew(n,tag,'pt3',pt3,wt,wt_noew,0._dp,5000._dp,50._dp)
        call bookrelew(n,tag,'Del y',dely,wt,wt_noew,-6._dp,6._dp,0.25_dp)
        wtplot_ew=wt
        wtplot_noew=wt_noew
        if  (m34 < 2000._dp) then
          wtplot_ew=zip
          wtplot_noew=zip
        endif
        call bookrelew(n,tag,'Del y, Mtt > 2 TeV',
     &                 dely,wtplot_ew,wtplot_noew,-6._dp,6._dp,0.25_dp)

c--- plots for both top rapidities < 2.5
        wtplot_ew=wt
        wtplot_noew=wt_noew
        if  ((abs(y3) > 2.5_dp) .or. (abs(y4) > 2.5_dp)) then
          wtplot_ew=zip
          wtplot_noew=zip
        endif
        call bookrelew(n,tag,'m34, |y|<2.5',
     &                 m34,wtplot_ew,wtplot_noew,0._dp,5000._dp,50._dp)
        call bookrelew(n,tag,'pt3, |y|<2.5',
     &                 pt3,wtplot_ew,wtplot_noew,0._dp,5000._dp,50._dp)
        if  (m34 < 2000._dp) then
          wtplot_ew=zip
          wtplot_noew=zip
        endif
        call bookrelew(n,tag,'Del y, Mtt > 2 TeV, |y|<2.5',
     &                 dely,wtplot_ew,wtplot_noew,-6._dp,6._dp,0.25_dp)

c--- plots for both top rapidities < 1.0
        wtplot_ew=wt
        wtplot_noew=wt_noew
        if  ((abs(y3) > 1._dp) .or. (abs(y4) > 1._dp)) then
          wtplot_ew=zip
          wtplot_noew=zip
        endif
        call bookrelew(n,tag,'m34, |y|<1.0',
     &                 m34,wtplot_ew,wtplot_noew,0._dp,5000._dp,50._dp)
        call bookrelew(n,tag,'pt3, |y|<1.0',
     &                 pt3,wtplot_ew,wtplot_noew,0._dp,5000._dp,50._dp)
        if  (m34 < 2000._dp) then
          wtplot_ew=zip
          wtplot_noew=zip
        endif
        call bookrelew(n,tag,'Del y, Mtt > 2 TeV, |y|<1.0',
     &                 dely,wtplot_ew,wtplot_noew,-6._dp,6._dp,0.25_dp)

c--- plots for computing asymmetry
        wtplot_noew=0._dp
        wtplot_ew=0._dp
        if (abs(y3)-abs(y4) .gt. 0._dp) then
          wtplot_noew=wt_noew
          wtplot_ew=wt
        endif
        call bookplot(n,tag,'Delta yabs > 0 - no EW',0.5_dp,
     &   wtplot_noew,wtplot_noew**2,0._dp,1._dp,1._dp,'lin')
        n=n+1
        call bookplot(n,tag,'Delta yabs > 0 - with EW',0.5_dp,
     &   wtplot_ew,wtplot_ew**2,0._dp,1._dp,1._dp,'lin')
        n=n+1
        wtplot_noew=0._dp
        wtplot_ew=0._dp
        if (abs(y3)-abs(y4) .lt. 0._dp) then
          wtplot_noew=wt_noew
          wtplot_ew=wt
        endif
        call bookplot(n,tag,'Delta yabs < 0 - no EW',0.5_dp,
     &   wtplot_noew,wtplot_noew**2,0._dp,1._dp,1._dp,'lin')
        n=n+1
        call bookplot(n,tag,'Delta yabs < 0 - with EW',0.5_dp,
     &   wtplot_ew,wtplot_ew**2,0._dp,1._dp,1._dp,'lin')
        n=n+1

      else

c--- plots for case with no EW corrections

        call bookplot(n,tag,'m34',m34,wt,wt2,0._dp,5000._dp,50._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt3',pt3,wt,wt2,0._dp,5000._dp,50._dp,'lin')
        n=n+1
        call bookplot(n,tag,'Del y',dely,wt,wt2,-6._dp,6._dp,0.25_dp,'lin')
        n=n+1
        wtplot_noew=0._dp
        if (m34 > 2000._dp) then
          wtplot_noew=wt
        endif
        call bookplot(n,tag,'Del y, Mtt > 2 TeV',dely,
     &   wtplot_noew,wtplot_noew**2,-6._dp,6._dp,0.25_dp,'lin')
        n=n+1
c--- aymmetry plots for runs without EW corrections
        wtplot_noew=0._dp
        if (dely .gt. 0._dp) then
          wtplot_noew=wt
        endif
        call bookplot(n,tag,'Delta y > 0',0.5_dp,
     &   wtplot_noew,wtplot_noew**2,0._dp,1._dp,1._dp,'lin')
        n=n+1
        wtplot_noew=0._dp
        if (dely .lt. 0._dp) then
          wtplot_noew=wt
        endif
        call bookplot(n,tag,'Delta y < 0 - no EW',0.5_dp,
     &   wtplot_noew,wtplot_noew**2,0._dp,1._dp,1._dp,'lin')
        n=n+1

      endif

      call bookplot(n,tag,'eta',eta,wt,wt2,0._dp,etam,50._dp,'lin')
      n=n+1

      call bookplot(n,tag,'y3',y3,wt,wt2,-ylep,ylep,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y4',y4,wt,wt2,-ylep,ylep,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-ylep,ylep,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt3',pt3,wt,wt2,0._dp,ptlep,2._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt4',pt4,wt,wt2,0._dp,ptlep,2._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt34',pt34,wt,wt2,0._dp,50._dp,2._dp,'lin')
      n=n+1
      call bookplot(n, tag,'m34',m34,wt,wt2,70._dp,110._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DeltaR35',r35,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y5',y5,wt,wt2,-yjet,yjet,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,0._dp,ptjet,2._dp,'lin')
      n=n+1
      
c--- compute lepton asymmetry as a function of m34  
c--- (see for example Eq.(3) of PLB718 (2013) 752)
      p3(:)=p(3,:)
      p4(:)=p(4,:)
      p34(:)=p(3,:)+p(4,:)
      costheta=p34(3)/abs(p34(3))
     & *((p3(4)+p3(3))*(p4(4)-p4(3))-(p3(4)-p3(3))*(p4(4)+p4(3)))
     & /m34/sqrt(m34**2+p34(1)**2+p34(2)**2)
c--- these histograms must be kept
      if ((costheta .gt. 0._dp) .or. (tag .eq. tagbook)) then
        call bookplot(n, tag,'m34 forward lepton',
     &   m34,wt,wt2,40._dp,200._dp,5._dp,'lin')
      endif
      n=n+1
      if ((costheta .le. 0._dp) .or. (tag .eq. tagbook)) then
        call bookplot(n, tag,'m34 backward lepton',
     &   m34,wt,wt2,40._dp,200._dp,5._dp,'lin')
      endif
      n=n+1
c--- now compute asymmetry - histograms n+1 and n+2 are only temporary
      if (tag .eq. tagbook) then
        call bookplot(n, tag,'lepton FB asymmetry',
     &   m34,wt,wt2,40._dp,200._dp,5._dp,'lin')
        call bookplot(n+1, tag,'lepton FB asymmetry',
     &   m34,wt,wt2,40._dp,200._dp,5._dp,'lin')
        call bookplot(n+2, tag,'lepton FB asymmetry',
     &   m34,wt,wt2,40._dp,200._dp,5._dp,'lin')
      endif
      call mopera(n-2,'-',n-1,n,1._dp,1._dp)
      call mopera(n-2,'+',n-1,n+1,1._dp,1._dp)
c--- this is the histogram (n) we will keep
      call mopera(n,'/',n+1,n,1._dp,1._dp)
      n=n+1

  777 continue

************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

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
      
