      subroutine mcfm_vegas_adaptive(integ,integ_err)
          implicit none
          include 'types.f'
          include 'ipsgen.f'
          include 'reset.f'
          include 'kpart.f'
          include 'kprocess.f'
          include 'vegas_common.f'
          include 'gridinfo.f'
          include 'nproc.f'
          include 'taucut.f'
          include 'mpif.h'
          include 'mpicommon.f'

          include 'parttypes.f'
          include 'vegas_damp.f'

          real(dp), intent(out) :: integ, integ_err

          real(dp) :: region(2*mxdim)
          real(dp) :: lowint, virtint, realint, scetint
          external lowint, virtint, realint, scetint
          character*1 getstr

          ! these hold values for all maxParts integration parts and a maximum of
          ! maxIps ipsgen contributions
          real(dp) :: sig(maxParts,maxIps), sd(maxParts,maxIps),
     &                chi(maxParts,maxIps)
          ! this array holds modified sd values to give
          ! easier integrations slight preference, see below how it's used
          real(dp) :: sdMod(maxParts,maxIps)
          integer(kind=8) :: ncalls(maxParts,maxIps)
          integer :: maxipsgens(maxParts)
          logical :: computePart(maxParts)
          logical :: doFirstCall(maxParts,maxIps)
          logical :: warmupComplete(maxParts,maxIps)

          ! this sets the vegas grid dampening parameter for all parts
          ! in warmup(1) and final(2) runs
          ! setting it to 0 will disable grid adjustments
          real(dp) :: vegasDampPart(maxParts,2)

          real(dp) :: xreal,xreal2
          common/xreal/xreal,xreal2

          integer :: stage
          integer :: ierr, ierr2
          integer :: hist_readsnapshot, hist_writesnapshot
          integer :: nprocabove
          logical :: origCoeffonly

          ! number of iterations to compute in one batch
          integer :: iterBatch = 5
          ! multiplicator for number of calls if precision per batch is
          ! too low
          real(dp), parameter :: iterCallMult = 1.4_dp
          ! precision goal in percent of one batch for warmup
          real(dp), parameter :: batchPrecisionGoal = 0.25_dp
          ! precision goal of final result
          real(dp), parameter :: resultPrecisionGoal = 1e-3_dp

          ! maximum number of calls per iteration
          integer, parameter :: maxCallsPerIter = 500d6

          integer :: nextLoc(2)
          integer :: i,j,k

          logical :: bin
          common/bin/bin

          logical :: dryrun
          common/dryrun/dryrun

          character*255 runname
          common/runname/runname


          integer(kind=8), parameter, dimension(1:13)   :: s = (/1,2,3,3,4,4,5,5,5,5,5,5,6/)
          integer(kind=8), parameter, dimension(1:13)   :: a = (/0,1,1,2,1,4,2,4,7,11,13,14,1/)
          integer(kind=8), parameter, dimension(6,1:13) :: m = reshape((/1,0,0,0,0,0,
     &    					1,3,0,0,0,0,
     &    					1,3,1,0,0,0,
     &    					1,1,1,0,0,0,
     &    					1,1,3,3,0,0,
     &    					1,3,5,13,0,0,
     &    					1,1,5,5,17,0,
     &    					1,1,5,5,5,0,
     &    					1,1,7,11,19,0,
     &    					1,1,5,1,1,0,
     &    					1,1,1,3,11,0,
     &    					1,3,5,5,31,0,
     &                      1,3,3,9,7,49/), (/6,13/))
          integer(kind=8) :: toskip

          ! optional TODO:
          ! we possibly also want to look at the chisq and set a goal
          ! for the warmup

          ! no binning for warumup stage
          bin = .false.

          ! we don't need or want this here
          dryrun = .false.

          sig(:,:) = 0._dp
          sd(:,:) = 0._dp
          chi(:,:) = 0._dp

          reset = .false.
          scalereset = .false.

          ! The common block variable maxipsgen is set
          ! to the corresponding maxipsgens value for each part below.
          ! The phase space generation routines can adapt on that value.
          ipsgen = 1
          maxipsgens(:) = 1

          origKpart = kpart
          origCoeffonly = coeffonly

          doFirstCall(:,:) = .false.
          warmupComplete(:,:) = .false.
          vegasDampPart(:,:) = 1.5

c---      setup inital number of calls for warmup phase
c---      these could also be tuned to specific processes
c---
c---      note that the MC integration error is only accurate if the 
c---      number of calls is "sufficiently" large, especially for 
c---      difficult integrands
          ncalls(lord,:) = 100000
          ncalls(nloVirt,:) = 50000
          ncalls(nloReal,:) = 100000
          ncalls(nnloBelow,:) = 200000
          ncalls(nnloVirtAbove,:) = 200000
          ncalls(nnloRealAbove,:) = 1000000
          ncalls(snloBelow,:) = 200000
          ncalls(snloAbove,:) = 200000

c--- ===================================== 
c--- process specific initializations here
c--- ===================================== 

          ! doipsgen here signifies that such an ipsgen phase space
          ! splitting has been performed for possibly any part

          if (nproc == 300 .or. nproc == 3000) then
            maxipsgens(lord) = 2
            maxipsgens(nloVirt) = 2
            maxipsgens(nloReal) = 2
            maxipsgens(nnloBelow) = 2
            maxipsgens(nnloVirtAbove) = 2
            maxipsgens(nnloRealAbove) = 2
            maxipsgens(snloAbove) = 2
            maxipsgens(snloBelow) = 2
          endif

          if (nproc == 302 .or. nproc == 3002) then
            maxipsgens(nloVirt) = 2
            maxipsgens(nloReal) = 2
            maxipsgens(lord) = 2
          endif

          if (nproc == 304 .or. nproc == 3004) then
            maxipsgens(lord) = 2
          endif

          if (nproc == 164 .or. nproc==169) then
            vegasDampPart(lord,1) = 1.5
            vegasDampPart(lord,2) = 0.8
            ncalls(lord,:) = 250000

            ncalls(nloVirt,:) = 100000
            vegasDampPart(nloVirt,1) = 1.5
            vegasDampPart(nloVirt,2) = 0.8

            maxipsgens(nloReal) = 2
            ncalls(nloReal,:) = 2000000
            vegasDampPart(nloReal,1) = 1.5
            vegasDampPart(nloReal,2) = 0.8
          endif

c---
c--- determine parts to calculate
c--- 
          computePart(:) = .false.

          if (origKpart == klord) then
              computePart(lord) = .true.
          endif

          if (origKpart == kvirt) then
              origCoeffonly = .true.
              computePart(nloVirt) = .true.
          endif

          if (origKpart == kreal) then
              origCoeffonly = .true.
              computePart(nloReal) = .true.
          endif

          if (origKpart == ktota) then
              computePart(nloVirt) = .true.
              computePart(nloReal) = .true.
          endif

          if (origKpart ==  ksnlo) then
              computePart(snloBelow) = .true.
              computePart(snloAbove) = .true.
          endif

          if (origKpart == knnlo) then
            if (origCoeffonly .eqv. .false.) then
              ! alternative enable snloBelow and snloAbove
              computePart(nloVirt) = .true.
              computePart(nloReal) = .true.
            endif

            computePart(nnloBelow) = .true.
            computePart(nnloVirtAbove) = .true.
            computePart(nnloRealAbove) = .true.
          endif

          if (origKpart == ksnlo .or. origKpart == knnlo) then
            call setupscet(nprocabove)
          endif

c ================ 
c === resume integration with previous grid
c === failure to read the resumption snapshot is handled gracefully
c ================ 

          call zeroaccumhist()

          if (readin) then
            open(unit=11, file=trim(runname)//"_vegas_adaptive_snapshot.dat", status='OLD', iostat=ierr)
            if (ierr == 0) then
              partLoop: do i=1,maxParts
                do j=1,maxIps
                  read(11, '(ES25.15,ES25.15,ES25.15,I20,L4,L4)', iostat=ierr)
     &                sd(i,j), sig(i,j), chi(i,j), ncalls(i,j), doFirstCall(i,j), warmupComplete(i,j)
                  if (ierr /= 0) then
                    write (*,*) "Broken vegas_adaptive_snapshot.dat!"
                    close(unit=11)
                    call abort
                  endif
                enddo
              enddo partLoop
              close(unit=11)
            endif

            ierr2 = hist_readsnapshot()
            if (ierr2 /= 0) then
              write (*,*) "Could not find histogram_snapshot.dat"
              ! this may not be a problem, since for the warmup phase
              ! no histograms are written
            endif

            if (ierr /=  0) then
              write (*,*) "Could not find vegas_adaptive_snapshot.dat"
              readin = .false.
            endif

          endif

c ================ 
c === warmup phase
c ================ 

c-----------
c--- lord
c-----------
          if (computePart(lord)) then

          ipsLordLoop: do ipsgen=1,maxipsgens(lord)

          if(warmupComplete(lord,ipsgen) .eqv. .true.) exit ipsLordLoop

          nproc = nprocbelow
          abovecut = .false.
          call chooser

          maxipsgen = maxipsgens(lord)

          kpart = klord
          reset = .true.
          scalereset = .true.
          usescet = .false.
          coeffonly = origCoeffonly
          call boundregion(ndim,region)
          outgridfile = 
     &         "dvegas_lord_PS"//getstr(ipsgen)//".grid"
          writeout = .true.
          readin = .false.
          stage = 0

          do
            if (rank == 0)
     &      write (*,*) "lord warmup integration"

            vegas_damp = vegasDampPart(lord,1)

            call vegasnr(region,ndim,lowint,stage,
     &                   ncalls(lord,ipsgen),iterBatch,
     &                   nprn,sig(lord,ipsgen),
     &                   sd(lord,ipsgen),chi(lord,ipsgen),lord,ipsgen)

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

            if (sd(lord,ipsgen)/abs(sig(lord,ipsgen))
     &                  > batchPrecisionGoal) then
                ncalls(lord,ipsgen) = ncalls(lord,ipsgen)
     &                      * iterCallMult
                if (rank == 0)
     &          write (*,*) "Increasing calls to ",
     &                ncalls(lord,ipsgen)
            else
                if (rank == 0)
     &          write (*,*) "Reached batchPrecisionGoal with ",
     &                ncalls(lord,ipsgen), " calls per iteration"
                ! warmup complete
                exit
            endif

          enddo

          doFirstCall(lord,ipsgen) = .true.
          warmupComplete(lord,ipsgen) = .true.

          if (rank==0) call writeVegasSnapshot(sd,sig,chi,ncalls,
     &                          doFirstCall,warmupComplete)

          enddo ipsLordLoop
          endif

c-----------
c--- snlo below
c-----------
          if (computePart(snloBelow)) then

          ipsSnloBelowLoop: do ipsgen=1,maxipsgens(snloBelow)

          if(warmupComplete(snloBelow,ipsgen) .eqv. .true.) exit ipsSnloBelowLoop

          nproc = nprocbelow
          abovecut = .false.
          call chooser

          maxipsgen = maxipsgens(snloBelow)

          kpart = ksnlo
          reset = .true.
          scalereset = .true.
          usescet = .true.
          coeffonly = origCoeffonly
          abovecut = .false.
          ndim=ndim+2
          call boundregion(ndim,region)
          outgridfile = 
     &         "dvegas_snlo_below_PS"//getstr(ipsgen)//".grid"
          writeout = .true.
          readin = .false.
          stage = 0

          do
            if (rank == 0)
     &      write (*,*) "snlo below warmup integration"

            vegas_damp = vegasDampPart(snloBelow,1)

            call vegasnr(region,ndim,scetint,stage,
     &                   ncalls(snloBelow,ipsgen),iterBatch,
     &                   nprn,sig(snloBelow,ipsgen),
     &                   sd(snloBelow,ipsgen),chi(snloBelow,ipsgen),
     &                   snloBelow,ipsgen)

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

            if (sd(snloBelow,ipsgen)/abs(sig(snloBelow,ipsgen))
     &                  > batchPrecisionGoal) then
                ncalls(snloBelow,ipsgen) = ncalls(snloBelow,ipsgen)
     &                      * iterCallMult
                if (rank == 0)
     &          write (*,*) "Increasing calls to ",
     &                ncalls(snloBelow,ipsgen)
            else
                if (rank == 0)
     &          write (*,*) "Reached batchPrecisionGoal with ",
     &                ncalls(snloBelow,ipsgen), " calls per iteration"
                ! warmup complete
                exit
            endif

          enddo

          ndim=ndim-2

          doFirstCall(snloBelow,ipsgen) = .true.
          warmupComplete(snloBelow,ipsgen) = .true.

          if (rank==0) call writeVegasSnapshot(sd,sig,chi,ncalls,
     &                          doFirstCall,warmupComplete)

          enddo ipsSnloBelowLoop
          endif


c-----------
c--- snlo above
c-----------

          if (computePart(snloAbove)) then

          ipsSnloAboveLoop: do ipsgen=1,maxipsgens(snloAbove)

          if(warmupComplete(snloAbove,ipsgen) .eqv. .true.) exit ipsSnloAboveLoop

          nproc = nprocabove
          abovecut = .true.
          call chooser

          maxipsgen = maxipsgens(snloAbove)

          kpart = klord
          reset = .true.
          scalereset = .true.
          usescet = .true.
          coeffonly = origCoeffonly
          call boundregion(ndim,region)
          outgridfile = 
     &         "dvegas_snlo_above_PS"//getstr(ipsgen)//".grid"
          writeout = .true.
          readin = .false.
          stage = 0

          do
            if (rank == 0)
     &      write (*,*) "snlo above warmup integration"

            vegas_damp = vegasDampPart(snloAbove,1)

            call vegasnr(region,ndim,lowint,stage,
     &                   ncalls(snloAbove,ipsgen),iterBatch,
     &                   nprn,sig(snloAbove,ipsgen),
     &                   sd(snloAbove,ipsgen),chi(snloAbove,ipsgen),
     &                   snloAbove,ipsgen)

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

            if (sd(snloAbove,ipsgen)/abs(sig(snloAbove,ipsgen))
     &                  > batchPrecisionGoal) then
                ncalls(snloAbove,ipsgen) = ncalls(snloAbove,ipsgen)
     &                      * iterCallMult
                if (rank == 0)
     &          write (*,*) "Increasing calls to ",
     &                ncalls(snloAbove,ipsgen)
            else
                if (rank == 0)
     &          write (*,*) "Reached batchPrecisionGoal with ",
     &                ncalls(snloAbove,ipsgen), " calls per iteration"
                ! warmup complete
                exit
            endif

          enddo

          doFirstCall(snloAbove,ipsgen) = .true.
          warmupComplete(snloAbove,ipsgen) = .true.

          if (rank==0) call writeVegasSnapshot(sd,sig,chi,ncalls,
     &                          doFirstCall,warmupComplete)

          enddo ipsSnloAboveLoop
          endif


c-----------
c--- nlo virt
c-----------

          if (computePart(nloVirt)) then

          ipsNloVirtLoop: do ipsgen=1,maxipsgens(nloVirt)

          if(warmupComplete(nloVirt,ipsgen) .eqv. .true.) exit ipsNloVirtLoop

          nproc = nprocbelow
          abovecut = .false.
          call chooser

          maxipsgen = maxipsgens(nloVirt)

          kpart = kvirt
          reset = .true.
          scalereset = .true.
          usescet = .false.
          coeffonly = origCoeffonly
          ndim = ndim + 1
          call boundregion(ndim,region)

          outgridfile = "dvegas_virt_PS"//getstr(ipsgen)//".grid"
          writeout = .true.
          readin = .false.
          stage =  0

          do
            if (rank == 0)
     &      write (*,*) "virt warmup integration"

            vegas_damp = vegasDampPart(nloVirt,1)

            call vegasnr(region,ndim,virtint,stage,
     &               ncalls(nloVirt,ipsgen),iterBatch,
     &               nprn,sig(nloVirt,ipsgen),
     &               sd(nloVirt,ipsgen),chi(nloVirt,ipsgen),nloVirt,ipsgen)

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

            if (sd(nloVirt,ipsgen)/abs(sig(nloVirt,ipsgen))
     &                  > batchPrecisionGoal) then
                ncalls(nloVirt,ipsgen) = ncalls(nloVirt,ipsgen)
     &                      * iterCallMult
                if (rank == 0)
     &          write (*,*) "Increasing calls to ",
     &                ncalls(nloVirt,ipsgen)
            else
                if (rank == 0)
     &          write (*,*) "Reached batchPrecisionGoal with ",
     &                ncalls(nloVirt,ipsgen), " calls per iteration"
                ! warmup complete
                exit
            endif
          enddo

          ndim = ndim - 1

          doFirstCall(nloVirt,ipsgen) = .true.
          warmupComplete(nloVirt,ipsgen) = .true.

          if (rank==0) call writeVegasSnapshot(sd,sig,chi,ncalls,
     &                          doFirstCall,warmupComplete)

          enddo ipsNloVirtLoop
          endif

c-----------
c--- nlo real
c-----------

          if (computePart(nloReal)) then

          ipsNloRealLoop: do ipsgen=1,maxipsgens(nloReal)

          if(warmupComplete(nloReal,ipsgen) .eqv. .true.) exit ipsNloRealLoop

          nproc = nprocbelow
          abovecut = .false.
          call chooser

          maxipsgen = maxipsgens(nloReal)

          kpart = kreal
          reset = .true.
          scalereset = .true.
          usescet = .false.
          coeffonly = .false.
          xreal = 0._dp
          xreal2 = 0._dp
          ndim=ndim+3
          call boundregion(ndim,region)
          outgridfile =
     &            "dvegas_real_PS"//getstr(ipsgen)//".grid"
          writeout = .true.
          readin = .false.
          stage = 0

          do
            if (rank == 0)
     &      write (*,*) "real warmup integration"

            vegas_damp = vegasDampPart(nloReal,1)
          
            call vegasnr(region,ndim,realint,stage,
     &               ncalls(nloReal,ipsgen),iterBatch,
     &               nprn,sig(nloReal,ipsgen),
     &               sd(nloReal,ipsgen),chi(nloReal,ipsgen),nloReal,ipsgen)

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

            if (sd(nloReal,ipsgen)/abs(sig(nloReal,ipsgen))
     &                  > batchPrecisionGoal) then
                ncalls(nloReal,ipsgen) = ncalls(nloReal,ipsgen)
     &                      * iterCallMult
                if (rank == 0)
     &          write (*,*) "Increasing calls to ",
     &                ncalls(nloReal,ipsgen)
            else
                if (rank == 0)
     &          write (*,*) "Reached batchPrecisionGoal with ",
     &                ncalls(nloReal,ipsgen), " calls per iteration"
                ! warmup complete
                exit
            endif
          enddo

          ndim=ndim-3

          doFirstCall(nloReal,ipsgen) = .true.
          warmupComplete(nloReal,ipsgen) = .true.

          if (rank==0) call writeVegasSnapshot(sd,sig,chi,ncalls,
     &                          doFirstCall,warmupComplete)

          enddo ipsNloRealLoop
          endif

c--------------
c--- nnlo below
c--------------

          if (computePart(nnloBelow)) then

          ipsNnloBelowLoop: do ipsgen=1,maxipsgens(nnloBelow)

          if(warmupComplete(nnloBelow,ipsgen) .eqv. .true.) exit ipsNnloBelowLoop

          nproc = nprocbelow
          abovecut = .false.
          call chooser

          maxipsgen = maxipsgens(nnloBelow)

          kpart=knnlo
          reset = .true.
          scalereset = .true.
          usescet = .true.
          coeffonly = .true.
          abovecut = .false.
          ndim=ndim+2
          call boundregion(ndim,region)
          outgridfile = 
     &         "dvegas_nnlo_below_PS"//getstr(ipsgen)//".grid"
          writeout = .true.
          readin = .false.
          stage = 0

          do
            if (rank == 0)
     &      write (*,*) "nnlo below warmup integration"

            vegas_damp = vegasDampPart(nnloBelow,1)

            call vegasnr(region,ndim,scetint,stage,
     &                   ncalls(nnloBelow,ipsgen),iterBatch,
     &                   nprn,sig(nnloBelow,ipsgen),
     &                   sd(nnloBelow,ipsgen),chi(nnloBelow,ipsgen),
     &                   nnloBelow,ipsgen)

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

            if (sd(nnloBelow,ipsgen)/abs(sig(nnloBelow,ipsgen))
     &                  > batchPrecisionGoal) then
                ncalls(nnloBelow,ipsgen) = ncalls(nnloBelow,ipsgen)
     &                      * iterCallMult
                if (rank == 0)
     &          write (*,*) "Increasing calls to ",
     &                ncalls(nnloBelow,ipsgen)
            else
                if (rank == 0)
     &          write (*,*) "Reached batchPrecisionGoal with ",
     &                ncalls(nnloBelow,ipsgen), " calls per iteration"
                ! warmup complete
                exit
            endif

          enddo

          ndim=ndim-2

          doFirstCall(nnloBelow,ipsgen) = .true.
          warmupComplete(nnloBelow,ipsgen) = .true.

          if (rank==0) call writeVegasSnapshot(sd,sig,chi,ncalls,
     &                          doFirstCall,warmupComplete)

          enddo ipsNnloBelowLoop
          endif

c--------------------------------
c--- above cut bits, call chooser
c--------------------------------

c-------------------
c--- nnlo virt above
c-------------------
  
          if (computePart(nnloVirtAbove)) then

          ipsNnloVirtAboveLoop: do ipsgen=1,maxipsgens(nnloVirtAbove)

          if(warmupComplete(nnloVirtAbove,ipsgen) .eqv. .true.) exit ipsNnloVirtAboveLoop

          nproc = nprocabove
          abovecut = .true.
          call chooser

          maxipsgen = maxipsgens(nnloVirtAbove)

          kpart=kvirt
          reset = .true.
          scalereset = .true.
          usescet = .true.
          coeffonly = .true.
          stage = 0
          ndim=ndim+1
          call boundregion(ndim,region)
          outgridfile = 
     &           "dvegas_nnlo_virt_above_PS"//getstr(ipsgen)//".grid"
          writeout = .true.
          readin = .false.

          do
            if (rank == 0)
     &          write (*,*) "nnlo virt above warmup integration"

            vegas_damp = vegasDampPart(nnloVirtAbove,1)

            call vegasnr(region,ndim,virtint,stage,
     &                   ncalls(nnloVirtAbove,ipsgen),iterBatch,
     &                   nprn,sig(nnloVirtAbove,ipsgen),
     &                   sd(nnloVirtAbove,ipsgen),
     &                   chi(nnloVirtAbove,ipsgen), nnloVirtAbove,ipsgen)

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

            if (sd(nnloVirtAbove,ipsgen)/abs(sig(nnloVirtAbove,ipsgen))
     &                  > batchPrecisionGoal) then
             ncalls(nnloVirtAbove,ipsgen) = ncalls(nnloVirtAbove,ipsgen)
     &                      * iterCallMult
             if (rank == 0)
     &       write (*,*) "Increasing calls to ",
     &                ncalls(nnloVirtAbove,ipsgen)
            else
                if (rank == 0)
     &          write (*,*) "Reached batchPrecisionGoal with ",
     &             ncalls(nnloVirtAbove,ipsgen), " calls per iteration"
                ! warmup complete
                exit
            endif

          enddo

          ndim=ndim-1

          doFirstCall(nnloVirtAbove,ipsgen) = .true.
          warmupComplete(nnloVirtAbove,ipsgen) = .true.

          if (rank==0) call writeVegasSnapshot(sd,sig,chi,ncalls,
     &                          doFirstCall,warmupComplete)

          enddo ipsNnloVirtAboveLoop
          endif

c-------------------
c--- nnlo real above
c-------------------

          if (computePart(nnloRealAbove)) then

          ipsNnloRealAboveLoop: do ipsgen=1,maxipsgens(nnloRealAbove)

          if(warmupComplete(nnloRealAbove,ipsgen) .eqv. .true.) exit ipsNnloRealAboveLoop

          nproc = nprocabove
          abovecut = .true.
          call chooser

          maxipsgen = maxipsgens(nnloRealAbove)

          kpart=kreal
          reset = .true.
          scalereset = .true.
          usescet = .true.
          coeffonly = .true.
          xreal = 0._dp
          xreal2 = 0._dp
          stage = 0
          ndim=ndim+3
          call boundregion(ndim,region)
          outgridfile = 
     &          "dvegas_nnlo_real_above_PS"//getstr(ipsgen)//".grid"
          writeout = .true.
          readin = .false.

          do
            if (rank == 0) 
     &         write (*,*) "nnlo real above warmup integration"

            vegas_damp = vegasDampPart(nnloRealAbove,1)

            call vegasnr(region,ndim,realint,stage,
     &                   ncalls(nnloRealAbove,ipsgen),iterBatch,
     &                   nprn,sig(nnloRealAbove,ipsgen),
     &                   sd(nnloRealAbove,ipsgen),
     &                   chi(nnloRealAbove,ipsgen), nnloRealAbove,ipsgen)

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

            if (sd(nnloRealAbove,ipsgen)/abs(sig(nnloRealAbove,ipsgen))
     &                  > batchPrecisionGoal) then
             ncalls(nnloRealAbove,ipsgen) = ncalls(nnloRealAbove,ipsgen)
     &                      * iterCallMult
             if (rank == 0) write (*,*) "Increasing calls to ",
     &                ncalls(nnloRealAbove,ipsgen)
            else
                if (rank == 0) 
     &             write (*,*) "Reached batchPrecisionGoal with ",
     &             ncalls(nnloRealAbove,ipsgen), " calls per iteration"
                ! warmup complete
                exit
            endif

          enddo

          ndim=ndim-3

          doFirstCall(nnloRealAbove,ipsgen) = .true.
          warmupComplete(nnloRealAbove,ipsgen) = .true.

          if (rank==0) call writeVegasSnapshot(sd,sig,chi,ncalls,
     &                          doFirstCall,warmupComplete)

          enddo ipsNnloRealAboveLoop
          endif

          call mpi_bcast(sd,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*maxIps,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)

          if (rank == 0) then
              write (*,*) "warmup result: ", sum(sig(:,:)), "+/-",
     &                       sqrt(sum(sd(:,:)**2))
          endif

c =========================== 
c === final integration phase
c =========================== 

          ! enable calls for nplotter
          bin = .true.

          integrationLoop: do
            nextLoc(:) = 0
            ipsgen = 0

            ! first check if there has been a contribution
            ! with no final run, but only warmup integration
            loopContrib: do i=1,maxParts
                loopIpsgen: do j=1,maxIps
                    if (doFirstCall(i,j) .eqv. .true.) then
                        nextLoc(1) = i
                        nextLoc(2) = j
                        ! do not load result from saved file
                        stage = 0
                        if (rank == 0) then
                            write (*,*) "first full integration"
                        endif
                        exit loopContrib
                    endif
                enddo loopIpsgen
            enddo loopContrib

            ! otherwise pick contribution with biggest sd
            if (nextLoc(1) == 0) then

                ! we fake the "real emission" contributions uncertainty
                ! for the prioritization system here
                ! to be half as large as they really are
                ! to boost the easier to compute contributions
                sdMod(:,:) = sd(:,:)
c               sdMod(nloReal,:) = sdMod(nloReal,:) / 2._dp
c               sdMod(snloAbove,:) = sdMod(snloAbove,:) / 2._dp
c               sdMod(nnloRealAbove,:) = sdMod(nnloRealAbove,:) / 2._dp

                nextLoc = maxloc(sdMod(:,:))
                ! load previous result from saved file
                stage = 2
            endif

            ! determine corresponding ipsgen contribution
            ipsgen = nextLoc(2)
            if (rank == 0) write (*,*) "ipsgen ",ipsgen," contribution"

            ! phase space generation can depend on this
            maxipsgen = maxipsgens(nextLoc(1))


              ! calculate precision here and
              ! exit if precision goal has been reached, but only if
              ! all contributions had at least on full integration run
            if (doFirstCall(nextLoc(1),nextLoc(2)) .eqv. .true.) then
              doFirstCall(nextLoc(1),nextLoc(2)) = .false.
              ! boost number of calls w.r.t. warmup phase
              nCalls(nextLoc(1),nextLoc(2)) = ncalls(nextLoc(1),nextLoc(2))*4
            else
              integ = sum(sig(:,:))
              integ_err = sqrt(sum(sd(:,:)**2))

              ! write intermediate histogram
              if (rank == 0) then
                call finalizehist()
                call histofin(integ,integ_err,0,1)
                ! set hist(:,:) to zero, because histofin overwrites
                ! hist with finalhist
                call zerohist()
              endif

              ! integration precision abort condition here
              if (abs(integ_err/integ) < resultPrecisionGoal) then
                exit integrationLoop
              endif

              ! if all integrals are zero
              if (sig(nextLoc(1),nextLoc(2)) == 0._dp) then
                if (rank == 0) write (*,*) "all integrals zero"
                exit integrationLoop
              endif

            endif


            ! Once we have reached maxCallsPerIter we just run one iteration
            ! per batch to have frequent snapshots and updated results.
            ! We set it here again for resumed runs.
            if (ncalls(nextloc(1),nextloc(2))*1.6 > maxCallsPerIter) then
              ncalls(nextloc(1), nextloc(2)) = maxCallsPerIter
              iterBatch = 1
            endif

            reset = .true.
            scalereset = .true.
            writeout = .true.
            readin = .true.

            select case (nextLoc(1))
c--- lord
              case (lord)
                  if (rank == 0) write (*,*) "lord selected"
                  nproc = nprocbelow
                  call chooser

                  kpart = klord
                  usescet = .false.
                  coeffonly = origCoeffonly
                  call boundregion(ndim, region)
                  outgridfile =
     &                    "dvegas_lord_PS"//getstr(ipsgen)//".grid"
                  ingridfile = outgridfile

                  vegas_damp = vegasDampPart(lord,2)
                  call vegasnr(region,ndim,lowint,stage,
     &                     ncalls(lord,ipsgen),iterBatch,
     &                     nprn,sig(lord,ipsgen),
     &                     sd(lord,ipsgen),chi(lord,ipsgen),lord,ipsgen)

c--- nlo virt
              case (nloVirt)
                  if (rank == 0) write (*,*) "nlo virt selected"
                  nproc = nprocbelow
                  call chooser


                  kpart = kvirt
                  usescet = .false.
                  coeffonly = origCoeffonly
                  ndim = ndim + 1
                  call boundregion(ndim,region)
                  outgridfile =
     &                    "dvegas_virt_PS"//getstr(ipsgen)//".grid"
                  ingridfile = outgridfile

                  vegas_damp = vegasDampPart(nloVirt,2)
                  call vegasnr(region,ndim,virtint,stage,
     &                     ncalls(nloVirt,ipsgen),iterBatch,
     &                     nprn,sig(nloVirt,ipsgen),
     &                     sd(nloVirt,ipsgen),chi(nloVirt,ipsgen),
     &                     nloVirt,ipsgen)
                  ndim = ndim - 1
c--- nlo real
              case (nloReal)
                  if (rank == 0) write (*,*) "nlo real selected"
                  nproc = nprocbelow
                  call chooser

                  kpart = kreal
                  usescet = .false.
                  coeffonly = .false.
                  xreal = 0._dp
                  xreal2 = 0._dp
                  ndim=ndim+3
                  call boundregion(ndim,region)
                  outgridfile =
     &                    "dvegas_real_PS"//getstr(ipsgen)//".grid"
                  ingridfile = outgridfile
                  
                  vegas_damp = vegasDampPart(nloReal,2)
                  call vegasnr(region,ndim,realint,stage,
     &                     ncalls(nloReal,ipsgen),iterBatch,
     &                     nprn,sig(nloReal,ipsgen),
     &                     sd(nloReal,ipsgen),chi(nloReal,ipsgen),
     &                     nloReal,ipsgen)
                  ndim=ndim-3
c--- nnlo below
              case (nnloBelow)  
                  if (rank == 0) write (*,*) "nnlo below selected"
                  nproc = nprocbelow
                  call chooser

                  kpart=knnlo
                  usescet = .true.
                  coeffonly = .true.
                  abovecut = .false.
                  ndim=ndim+2
                  call boundregion(ndim,region)
                  outgridfile = 
     &                 "dvegas_nnlo_below_PS"//getstr(ipsgen)//".grid"
                  ingridfile = outgridfile

                  vegas_damp = vegasDampPart(nnloBelow,2)
                  call vegasnr(region,ndim,scetint,stage,
     &                     ncalls(nnloBelow,ipsgen),iterBatch,
     &                     nprn,sig(nnloBelow,ipsgen),
     &                     sd(nnloBelow,ipsgen),chi(nnloBelow,ipsgen),
     &                     nnloBelow,ipsgen)
                  ndim=ndim-2
c--- nnlo virt above
              case (nnloVirtAbove)
                  if (rank == 0) write (*,*) "nnlo virt above selected"
                  nproc = nprocabove
                  call chooser

                  kpart = kvirt
                  usescet = .true.
                  coeffonly = .true.
                  abovecut = .true.
                  ndim = ndim+1
                  call boundregion(ndim,region)
                  outgridfile = 
     &              "dvegas_nnlo_virt_above_PS"//getstr(ipsgen)//".grid"
                  ingridfile = outgridfile

                  vegas_damp = vegasDampPart(nnloVirtAbove,2)
                  call vegasnr(region,ndim,virtint,stage,
     &                       ncalls(nnloVirtAbove,ipsgen),iterBatch,
     &                       nprn,sig(nnloVirtAbove,ipsgen),
     &                       sd(nnloVirtAbove,ipsgen),
     &                       chi(nnloVirtAbove,ipsgen),
     &                       nnloVirtAbove,ipsgen)

                  ndim = ndim-1
c--- nnlo real above
              case (nnloRealAbove)
                  if (rank == 0) write (*,*) "nnlo real above selected"
                  nproc = nprocabove
                  call chooser

                  kpart = kreal
                  usescet = .true.
                  coeffonly = .true.
                  abovecut = .true.
                  xreal = 0._dp
                  xreal2 = 0._dp
                  ndim = ndim+3
                  call boundregion(ndim,region)
                  outgridfile = 
     &              "dvegas_nnlo_real_above_PS"//getstr(ipsgen)//".grid"
                  ingridfile = outgridfile

                  vegas_damp = vegasDampPart(nnloRealAbove,2)
                  call vegasnr(region,ndim,realint,stage,
     &                         ncalls(nnloRealAbove,ipsgen),iterBatch,
     &                         nprn,sig(nnloRealAbove,ipsgen),
     &                         sd(nnloRealAbove,ipsgen),
     &                         chi(nnloRealAbove,ipsgen),
     &                         nnloRealAbove,ipsgen)

                  ndim = ndim-3
c--- snlo below
              case (snloBelow)
                  if (rank == 0) write (*,*) "snlo below selected"
                  nproc = nprocbelow
                  call chooser

                  kpart=ksnlo
                  usescet = .true.
                  coeffonly = origCoeffonly
                  abovecut = .false.
                  ndim=ndim+2
                  call boundregion(ndim,region)
                  outgridfile = 
     &                 "dvegas_snlo_below_PS"//getstr(ipsgen)//".grid"
                  ingridfile = outgridfile

                  vegas_damp = vegasDampPart(snloBelow,2)
                  call vegasnr(region,ndim,scetint,stage,
     &                       ncalls(snloBelow,ipsgen),iterBatch,
     &                       nprn,sig(snloBelow,ipsgen),
     &                       sd(snloBelow,ipsgen),chi(snloBelow,ipsgen),
     &                       snloBelow,ipsgen)
                  ndim=ndim-2
c--- snlo above
              case (snloAbove)
                  if (rank == 0) write (*,*) "snlo above selected"
                  nproc = nprocabove
                  call chooser

                  kpart = klord
                  usescet = .true.
                  coeffonly = .true.
                  abovecut = .true.
                  call boundregion(ndim,region)
                  outgridfile = 
     &              "dvegas_snlo_above_PS"//getstr(ipsgen)//".grid"
                  ingridfile = outgridfile

                  vegas_damp = vegasDampPart(snloAbove,2)
                  call vegasnr(region,ndim,lowint,stage,
     &                       ncalls(snloAbove,ipsgen),iterBatch,
     &                       nprn,sig(snloAbove,ipsgen),
     &                       sd(snloAbove,ipsgen),chi(snloAbove,ipsgen),
     &                       snloAbove,ipsgen)
              case default
                  write (*,*) "bad next contribution"
                  stop
            end select

            ! keep increasing calls up to maxCallsPerIter
            if (ncalls(nextloc(1),nextloc(2))*1.6 > maxCallsPerIter) then
              ncalls(nextloc(1), nextloc(2)) = maxCallsPerIter
              ! once we have reached maxCallsPerIter we just run one iteration
              ! per batch to have frequent snapshots and updated results
              iterBatch = 1
            else
              ncalls(nextloc(1), nextloc(2)) = 
     &              ncalls(nextloc(1), nextloc(2)) * 1.6
            endif

            if (rank == 0)
     &      write (*,*) "chisq last: ", chi(nextloc(1),nextloc(2))

          call mpi_bcast(sd,maxParts*2,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(sig,maxParts*2,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)
          call mpi_bcast(chi,maxParts*2,mpi_double_precision,
     &                    0,mpi_comm_world,ierr)


            if (rank == 0) then
              ! only print result if all parts had at least one
              ! full integration after warmup
              if (all(doFirstCall .eqv. .false.)) then
                write (*,*) "result: ", sum(sig(:,:)), "+/-",
     &                         sqrt(sum(sd(:,:)**2))
              else
                write (*,*) "partial warmup result: ",
     &                  sum(sig(:,:)), "+/-", sqrt(sum(sd(:,:)**2))
              endif

              call writeVegasSnapshot(sd,sig,chi,ncalls,doFirstCall,warmupComplete)

              ierr2 = hist_writesnapshot()
              if (ierr /= 0) then
                write (*,*) "Problem writing histogram_snapshot.dat"
              endif

            endif

          enddo integrationLoop

      end subroutine

      subroutine writeVegasSnapshot(sd,sig,chi,ncalls,doFirstCall,warmupComplete)
        implicit none
        include 'types.f'
        include 'histo.f'
        include 'accumhist.f'

        real(dp), intent(in) :: sig(maxParts,maxIps), sd(maxParts,maxIps),
     &                chi(maxParts,maxIps)
        integer(kind=8), intent(in) :: ncalls(maxParts,maxIps)
        logical, intent(in) :: doFirstCall(maxParts,maxIps)
        logical, intent(in) :: warmupComplete(maxParts,maxIps)

        integer :: i,j

        character*255 runname
        common/runname/runname

        ! save snapshot
        open(unit=11, file=trim(runname)//"_vegas_adaptive_snapshot.dat", status='REPLACE')
        do i=1,maxParts
          do j=1,maxIps
          write(11, '(ES25.15,ES25.15,ES25.15,I20,L4,L4)')
     &        sd(i,j), sig(i,j), chi(i,j), ncalls(i,j), doFirstCall(i,j), warmupComplete(i,j)
          enddo
        enddo
        close(unit=11)

      end subroutine
