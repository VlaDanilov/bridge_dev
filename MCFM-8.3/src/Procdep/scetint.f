      function scetint(r,wgt)
      implicit none
      include 'types.f'
      real(dp):: scetint
      include 'constants.f'
      include 'mxpart.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'kprocess.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'npart.f'
      include 'kpart.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'taucut.f'
      include 'first.f'
      include 'scalevar.f'
      include 'scale.f'
      include 'facscale.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'nlooprun.f'
      integer itrial
      real(dp):: savescale,savefacscale,xmsqvar(6),alphas
      real(dp):: p(mxpart,4),pjet(mxpart,4),r(mxdim),W,xmsq,
     & val,val2,pswt,xjac,wgt,QB(2),z1,z2,flux,BrnRat
      integer j
      logical:: bin,includedipole
      integer, save :: iorder
      common/bin/bin
      common/BrnRat/BrnRat
!$omp threadprivate(iorder)

      scetint=0._dp

      W=sqrts**2

      p(:,:)=0._dp
      pjet(:,:)=0._dp

      call gen_lops(r,p,pswt,*999)
    
      call dotem(npart+2,p,s)

c----reject event if any s(i,j) is too small
c      call smalls(s,npart,*999)
c----reject event if any tau is too small
      call smalltau(p,npart,*999)

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
      
      if (doscalevar) then
        itrial=maxscalevar+1
        savescale=scale
        savefacscale=facscale
      endif

   66 continue

      if (dynamicscale) call scaleset(initscale,initfacscale,p)
     
! adjust scales for scale variation
      if (doscalevar) then
        if (dynamicscale .eqv. .false.) then
          scale=savescale
          facscale=savefacscale
        endif
        if      (itrial == 7) then
           facscale=facscale/two
         elseif (itrial == 6) then
           facscale=facscale*two
         elseif (itrial == 5) then
           scale=scale/two
         elseif (itrial == 4) then
           scale=scale*two
         elseif (itrial == 3) then
           scale=scale/two
           facscale=facscale/two
         elseif (itrial == 2) then
           scale=scale*two
           facscale=facscale*two
         endif
         musq=scale**2
         as=alphas(scale,amz,nlooprun)
         ason2pi=as/twopi
         ason4pi=as/fourpi
         gsq=fourpi*as
      endif
      
      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts

      z1=r(ndim-1)**2
      z2=r(ndim)**2
      xjac=four*sqrt(z1*z2)

      flux=fbGeV2/(two*xx(1)*xx(2)*W)

      QB(1)=-two*p(1,4)
      QB(2)=-two*p(2,4)
      
      if (tauboost) then
! provide beam energies in singlet c.o.m. instead if required
        QB(1)=sqrt(QB(1)*QB(2))
        QB(2)=QB(1)
      endif

c--- determine order of calculation on first call
      if (first) then
        first=.false.
        if (kpart==knnlo) then
          iorder=2
        elseif (kpart==ksnlo) then
          iorder=1
        else
          write(6,*) 'Error in scetint: kpart=',kpart
          stop
        endif
      endif

c--- Calculate the required matrix elements
      if     (kcase==kW_only) then
        call lumxmsq_w(p,xx,z1,z2,QB,iorder,xmsq)
      elseif (kcase==kZ_only) then
        call lumxmsq_z(p,xx,z1,z2,QB,iorder,xmsq)
      elseif ((kcase==kggfus0) .or. (kcase==kHigaga)) then
        call lumxmsq_h(p,xx,z1,z2,QB,iorder,xmsq)
      elseif (kcase==kHi_Zga) then
        call lumxmsq_h_Zga(p,xx,z1,z2,QB,iorder,xmsq)
      elseif ((kcase==kWHbbar) .or. (kcase==kWHgaga)
     &   .or. (kcase==kWH__WW)) then
        call lumxmsq_wh(p,xx,z1,z2,QB,iorder,xmsq)
      elseif ((kcase==kZHbbar) .or. (kcase==kZHgaga)
     &   .or. (kcase==kZH__WW)) then
        call lumxmsq_zh(p,xx,z1,z2,QB,iorder,xmsq)
      elseif (kcase==kgamgam) then
         call lumxmsq_gaga(p,xx,z1,z2,QB,iorder,xmsq)
!      elseif (kcase==kW_1jet) then
!        call lumxmsq_w1jet(p,xx,z1,z2,QB,iorder,xmsq)
!      elseif (kcase==kZ_1jet) then
!        call lumxmsq_z1jet(p,xx,z1,z2,QB,iorder,xmsq)
!      elseif (kcase==kggfus1) then
!        call lumxmsq_h1jet(p,xx,z1,z2,QB,iorder,xmsq)
      elseif (kcase==kZgamma) then
         call set_anomcoup(p)
         call lumxmsq_zgamma(p,xx,z1,z2,QB,iorder,xmsq)
      elseif (kcase==kWgamma) then
         call lumxmsq_wgamma(p,xx,z1,z2,QB,iorder,xmsq)
      else
        write(6,*) 'Process not yet available in jettiness formalism'
        stop
      endif

! compute weights for scale variation, looping if necessary
      if (doscalevar) then
        itrial=itrial-1
        if (itrial > 0) then
          xmsqvar(itrial)=xmsq
          goto 66
        endif
        if (abs(xmsq) > zip) then
          scalereweight(:)=xmsqvar(:)/xmsq
        else
          scalereweight(:)=zip
        endif
      endif

      scetint=flux*xjac*pswt*xmsq/BrnRat

      call getptildejet(0,pjet)
      val=scetint*wgt 
      val2=val**2 
      
      if (val == val) then
        continue
      else
        write(6,*) 'Discarded NaN, val=',val
        write(6,*) 'flux,xjac,pswt,xmsq,BrnRat',flux,xjac,pswt,xmsq,BrnRat
!        write(6,*) 'wgt',wgt
!        write(6,*) '      r(1)=',r(1),'_dp'
!        write(6,*) '      r(2)=',r(2),'_dp'
!        write(6,*) '      r(3)=',r(3),'_dp'
!        write(6,*) '      r(4)=',r(4),'_dp'
!        write(6,*) '      r(5)=',r(5),'_dp'
!        write(6,*) '      r(6)=',r(6),'_dp'
        call flush(6)
        goto 999
      endif

      if (bin) then
         call nplotter(pjet,val,val2,0)
      endif

      return

 999  continue
      scetint=zip
      
      return
      end


