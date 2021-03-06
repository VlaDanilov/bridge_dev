      program mcfm
      
      include 'types.f'
      include "omp_lib.h"
      include 'mpicommon.f'
      
      COMMON /check/tt2,tt3,tb,te,tt1,ti1,ti2
      real(4):: t_mcfm_omp1,t_mcfm_omp2,tt2,tt3,tb,te,tt1,tsub1,tsub2,ti1,ti2

      integer:: t,threadmin,threadmax,trials,trial
      integer:: mth
      parameter (mth=240)
      real(dp):: rone,r0,r(mth),er0,er(mth)
      real(dp):: tone,t0,t1,time(mth)
      real(dp):: trialmin(mth),trialmax(mth),trialavg(mth)
      real(dp):: speedmin(mth),speedmax(mth),speedavg(mth)
      
      t_mcfm_omp1 = secnds (0.0)   

      rank=0
      size=1
     
      tone=1e0
      threadmax=omp_get_max_threads()
!set threadmin=threadmax to run at maximum number of threads only
      threadmin=1
      threadmin=threadmax
      if (threadmax>mth) then
         write(*,*) 'increase parameter mth to at least ',threadmax
         stop
      endif
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) 'maximum number of threads available: ',threadmax
      
      do t=threadmin,threadmax
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) '>>>>>>>>',t,' threads <<<<<<<<'
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         call omp_set_num_threads(t)
         t0=omp_get_wtime()
         tsub1 = secnds(0.0)
         call mcfmsub(r0,er0)
         tsub2 = secnds(tsub1)
         t1=omp_get_wtime()
         if (t==1) then 
            tone=t1-t0
            rone=r0
         endif
         time(t)=t1-t0
         r(t)=r0
         er(t)=er0
      enddo
c      write(*,*)
c      write(*,*)
c      write(*,*)
c      write(*,*)
      write(*,*)
      write(*,*) 'timing results:'
      do t=threadmin,threadmax
         if (threadmin .ne. threadmax) then
            write(*,1) t,time(t),tone/time(t),r(t)
         else
            write(*,2) t,time(t),r(t)
         endif
      enddo
      
      if (threadmin .ne. threadmax) then
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) "set 'threadmin=threadmax' in 'src/Need/mcfm_omp.f'"
      write(*,*) "for normal running with max number of threads only" 
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      endif

 1    format('threads: ',i3,' time: ',f8.2,'  ratio: ',f8.2,
     &     ' x-section: ',e18.12,' fb')
 2    format('threads: ',i3,' time: ',f8.2,
     &     '    x-section: ',e18.12,' fb')
     
      t_mcfm_omp2 = secnds (t_mcfm_omp1)   
      print *, " Time of mcfmsub.f = ",tsub2, "seconds."       
      print *, " Time of determinefilenames() = ",tt2, "seconds."
      print *, " Time of mcfmmain()  = ",tt3, "seconds."      
      print *, " Time of mcfm_init() = ",ti2, "seconds."
      print *, " Time of VEGAS THINGS = ",te, "seconds."
      print *, " Total time of mcfm_omp = ",t_mcfm_omp2, "seconds."
      
      
      end
