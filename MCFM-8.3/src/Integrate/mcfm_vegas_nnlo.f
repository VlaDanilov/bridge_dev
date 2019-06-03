      subroutine mcfm_vegas_nnlo(myinit,myitmx,myncall,mybin,xinteg,xerr)
      implicit none
      include 'types.f'
************************************************************************
*                                                                      *
*  Modified version of mcfm_vegas that computes NNLO                   *
*  cross-section using dipole subtraction for NLO contribution         *
*  and SCET for the NNLO coefficient                                   *
*                                                                      *
*  This routine should perform the sweeps of vegasnr                   *
*                                                                      *
*    Input parameters:                                                 *
*       myinit  :  the vegasnr routine entry point                     *
*       myitmx  :  the number of vegasnr sweeps                        *
*      myncall  :  the number of iterations per sweep                  *
*          bin  :  whether or not the results should be histogrammed   *
*                                                                      *
*    Returned variables:                                               *
*       xinteg  :  value of integration                                *
*         xerr  :  integration error                                   *
*                                                                      *
*    coeffonly =  .true.   -->  NNLO coefficient only                  *
*    coeffonly =  .false.  -->  full NNLO result                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'gridinfo.f'
      include 'realwt.f'
      include 'scale.f'
      include 'facscale.f'
      include 'vegas_common.f'
      include 'PDFerrors.f'
      include 'frag.f'
      include 'reset.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'ipsgen.f'
      include 'kpart.f'
      include 'taucut.f'
      include 'nproc.f'
      include 'mpicommon.f'
      include 'parttypes.f'
      integer myitmx,myinit,i,j,k,mynproc,nprocabove
      integer(kind=8) myncall,myncall_save
      integer:: ierr
      logical:: mybin,bin,coeffonly_save
      real(dp):: sigr,sdr,sigdk,sddk,chidk,
     & sigfrag,sdfrag,chifrag,sigWdk,sdWdk,chiWdk,
     & xreal,xreal2,xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsigf,sumsd,sumsdr,sumsdf,
     & xcallwt,sig(5,4),sd(5,4),chi(5,4)
      integer mykpart
      character*3 getstr,psgen
      common/mykpart/mykpart
      common/bin/bin
      common/xreal/xreal,xreal2
      real(dp):: lowint,virtint,realint,fragint,scetint,qtint
      real(dp):: region(2*mxdim),lord_bypart(-1:1,-1:1)
      logical:: first,myreadin
      common/bypart/lord_bypart
      external lowint,virtint,realint,fragint,scetint,qtint
      data first/.true./
      save first
           
c--- Initialize all integration results to zero, so that the
c--- total of all contributions may be combined at the end

c--- NLO: virt
c--- NLO: real
c--- NNLO: below cut
c--- NNLO: above cut virt
c--- NNLO: above cut real

      ipsgen = 1

      sig(:,:)=zip
      sd(:,:)=zip
      chi(:,:)=zip
      
      lord_bypart(:,:)=0d0

      if (PDFerrors) then
        PDFxsec(:)=zip
      endif

c--- Controls behaviour of gen_njets: need to reset phase-space
c--- boundaries when going from virt to real (using tota)
c--- need to reset scale also, for special scalestart values
      reset=.false.
      scalereset=.false.

c--- Put the vegasnr parameters in the common block
      itmx=myitmx
      ncall=myncall
      bin=mybin
      
c--- Store value of kpart in mykpart, which will be retained;
c--- also store value of scale in myscale, which will be retained;
c--- part and scale can be changed to make sure that the tota option works.
      mykpart=kpart
      myscale=scale
      myfacscale=facscale
      mynproc=nproc
      mymb=mb
      myncall_save=myncall

c--- skip NLO calculation if coeffonly=.true.
      if (coeffonly) goto 77
      
c-------------------------------------------------------------------------------
c------------------------------ NLO calculation --------------------------------
c-------------------------------------------------------------------------------

      usescet=.false.
      abovecut=.false.
      coeffonly=.false.

c--- set up the grid info for virtual
      if (first .and. (myinit == 1)) then
c-- special input name for virtual grid
          ingridfile='dvegas_virt_below_'//ingridfile
          myreadin=readin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_virt_below.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            outgridfile='dvegas_virt_below_PS'//psgen(1:1)//'.grid'
          endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_virt_below.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            ingridfile='dvegas_virt_below_PS'//psgen(1:1)//'.grid'
          endif
        endif
      endif

c--- Virtual integration should have one extra dimensions
c--- (added and then taken away)
      kpart=kvirt
      reset=.true.
      scalereset=.true.
      ndim=ndim+1
      call boundregion(ndim,region)
      call vegasnr(region,ndim,virtint,myinit,myncall,myitmx,
     &  nprn,sig(1,ipsgen),sd(1,ipsgen),chi(1,ipsgen),snloBelow,ipsgen)
      ndim=ndim-1
            
c--- set up the grid info for real integration
      if (first .and. (myinit == 1)) then
c-- special input name for real grid
        ingridfile(8:11)='real'
        readin=myreadin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_real_below.grid'          
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            outgridfile='dvegas_real_below_PS'//psgen(1:1)//'.grid'
        endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_real_below.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            ingridfile='dvegas_real_below_PS'//psgen(1:1)//'.grid'
        endif
        endif
      endif

c--- Real integration should have three extra dimensions
      scale=myscale
      facscale=myfacscale
      kpart=kreal
      reset=.true.
      xreal=0d0
      xreal2=0d0
      adjust=(real(ndim+3,dp))/(real(ndim+1,dp))
      ncall=int(real(myncall,dp)**adjust,kind=8)/2
      if (ncall > myncall*10) ncall=myncall*10
      if (rank == 0) then
        write(6,*) 'Adjusting number of points for real to',ncall
      endif
      ndim=ndim+3
      call boundregion(ndim,region)
      call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
     &    nprn,sig(2,ipsgen),sd(2,ipsgen),chi(2,ipsgen),snloAbove,ipsgen)
      ndim=ndim-3
      if (rank == 0) write(6,*) 
      ncall=myncall

   77 continue
c-------------------------------------------------------------------------------
c-------------------------- SCET NNLO calculation ------------------------------
c-------------------------------------------------------------------------------

      coeffonly_save=coeffonly

      usescet=.true.
      coeffonly=.true.
      
c--- set up SCET variables
      call setupscet(nprocabove)

!!! DEBUG: only compute RR !!!
!      nproc=nprocabove
!      abovecut=.true.
!      call chooser
!      goto 667
!!! DEBUG: only compute RR !!!

  665 continue
      if(doipsgen) then
          psgen=getstr(ipsgen)
          write(6,*) '********* Phase space region ',ipsgen,' *********'
      endif
            
c--- SCET below-cut contribution
c--- integration should have two extra dimensions (added and then taken away)
      if (first .and. (myinit == 1)) then
c-- special input name for SCET grid
          ingridfile='dvegas_scet_below_'//ingridfile
          myreadin=readin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_scet_below.grid'
          if (doipsgen) then
            outgridfile='dvegas_scet_below_PS'//psgen(1:1)//'.grid'
          endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_scet_below.grid'
          if (doipsgen) then
            ingridfile='dvegas_scet_below_PS'//psgen(1:1)//'.grid'
          endif
        endif
      endif
      kpart=mykpart
      reset=.true.
      scalereset=.true.
      ndim=ndim+2
      abovecut=.false.
      if ((knnlopart == 0) .or. (knnlopart == knnloVV)) then
        call boundregion(ndim,region)
        call vegasnr(region,ndim,scetint,myinit,myncall,myitmx,
     &     nprn,sig(3,ipsgen),sd(3,ipsgen),chi(3,ipsgen),nnloBelow,ipsgen)
! Uncomment these lines (and comment-out the above two) to switch
! from jettiness-slicing to QT-slicing
!      call vegasnr(region,ndim,qtint,myinit,myncall,myitmx,
!     &             nprn,sig(3,ipsgen),sd(3,ipsgen),chi(3,ipsgen),nnloBelow,ipsgen)
      else
        sig(3,ipsgen)=zip
        sd(3,ipsgen)=zip
        chi(3,ipsgen)=zip
      endif
      ndim=ndim-2

      ! next ipsgen scet_below piece
      if(doipsgen .and. ipsgen < maxipsgen) then
          ipsgen = ipsgen + 1
          goto 665
      endif

! if we're only computing power corrections, nothing more to be done
        if (onlypowcorr) goto 33 

      ipsgen=1

      nproc=nprocabove
c--- scale up number of points more for W/Z/H+jet processes
      if ( (nprocabove == 22) .or. (nprocabove == 27)
     &.or. (nprocabove == 44) .or. (nprocabove == 270)
     &.or. (nprocabove == 271) .or. (nprocabove == 272) ) then
        myncall=myncall*20
      else
        myncall=myncall*10
      endif
      abovecut=.true.
      call chooser
      if (first .eqv. .true.) then
        readin=.false.
        writeout=.true.
        outgridfile='dvegas_scet_above.grid'
      else
        readin=.true.
        writeout=.false.
        ingridfile='dvegas_scet_above.grid'
      endif

  666 continue
      if(doipsgen) then
          psgen=getstr(ipsgen)
          write(6,*) '********* Phase space region ',ipsgen,' *********'
      endif
            
c--- set up the grid info for virtual
      if (first .and. (myinit == 1)) then
c-- special input name for virtual grid
          ingridfile='dvegas_virt_above_'//ingridfile
          myreadin=readin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_virt_above.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            outgridfile='dvegas_virt_above_PS'//psgen(1:1)//'.grid'
          endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_virt_above.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            ingridfile='dvegas_virt_above_PS'//psgen(1:1)//'.grid'
          endif
        endif
      endif

c--- Virtual integration should have one extra dimensions
c--- (added and then taken away)
      kpart=kvirt
      reset=.true.
      scalereset=.true.
      ndim=ndim+1
! include this part if performing a full calculation, or just RV
      if ((knnlopart == 0) .or. (knnlopart == knnloRV)) then
        call boundregion(ndim,region)
        call vegasnr(region,ndim,virtint,myinit,myncall,myitmx,
     &     nprn,sig(4,ipsgen),sd(4,ipsgen),chi(4,ipsgen),nnloVirtAbove,ipsgen)
      else
        sig(4,ipsgen)=zip
        sd(4,ipsgen)=zip
        chi(4,ipsgen)=zip
      endif
      ndim=ndim-1

      ! next ipsgen virt_above piece
      if(doipsgen .and. ipsgen < maxipsgen) then
          ipsgen = ipsgen + 1
          goto 666
      endif
      ipsgen=1

c---
c--- now handle real_above
c---

!      ! ipsgen process specific setup
!      if(nprocabove==302) then
!          doipsgen = .true.
!          ipsgen = 1
!          maxipsgen = 2
!      endif

  667 continue
      if(doipsgen) then
          psgen=getstr(ipsgen)
          write(6,*) '********* Phase space region ',ipsgen,' *********'
      endif
            
c--- set up the grid info for real integration
      if (first .and. (myinit == 1)) then
c-- special input name for real grid
        ingridfile(8:11)='real'
        readin=myreadin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_real_above.grid'          
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            outgridfile='dvegas_real_above_PS'//psgen(1:1)//'.grid'
        endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_real_above.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            ingridfile='dvegas_real_above_PS'//psgen(1:1)//'.grid'
        endif
        endif
      endif

c--- Real integration should have three extra dimensions
      scale=myscale
      facscale=myfacscale
      kpart=kreal
      reset=.true.
      xreal=0d0
      xreal2=0d0
      adjust=(real(ndim+3,dp))/(real(ndim+1,dp))
      ncall=int(real(myncall,dp)**adjust,kind=8)/2
      if (ncall > myncall*40) ncall=myncall*40
      if (rank == 0) then
        write(6,*) 'Adjusting number of points for real to',ncall
      endif
      ndim=ndim+3
! include this part if performing a full calculation, or just RR
      if ((knnlopart == 0) .or. (knnlopart == knnloRR)) then
        call boundregion(ndim,region)
        call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
     &     nprn,sig(5,ipsgen),sd(5,ipsgen),chi(5,ipsgen),nnloRealAbove,ipsgen)
      else
        sig(5,ipsgen)=zip
        sd(5,ipsgen)=zip
        chi(5,ipsgen)=zip
      endif
      ndim=ndim-3
      if (rank == 0) write(6,*) 
      ncall=myncall

      ! next ipsgen real_above piece
      if(doipsgen .and. ipsgen < maxipsgen) then
          ipsgen = ipsgen + 1
          goto 667
      endif
      ipsgen=1

   33 continue

c--- return nproc to the value from the input file
      nproc=mynproc
      if (first) call chooser

c--- calculate integration variables to be returned
      xinteg = sum(sig(:,:))
      xerr = sqrt(sum(sd(:,:)**2))
      
c--- return part, scale, myncall and coeffonly to their real values
      kpart=mykpart
      scale=myscale
      myncall=myncall_save
      coeffonly=coeffonly_save
      first=.false.

      return
      end
